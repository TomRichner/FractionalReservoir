function out = fig_memory_capacity(cfg)
% FIG_MEMORY_CAPACITY Paper-ready memory-capacity strip, from a saved MC run.
%
%   out = FIG_MEMORY_CAPACITY()
%   out = FIG_MEMORY_CAPACITY('mat_file', f)
%   out = FIG_MEMORY_CAPACITY('run_dir', d)   % looks in <run_dir>/memory_capacity
%
% A 1 x 3 strip assembled from a finished run_memory_capacity ensemble:
%   (a) cumulative memory capacity against delay
%   (b) per-delay R^2
%   (c) memory horizon, paired across trials
%
% No simulation is re-run.
%
% RESOLUTION ORDER for the results file: an explicit mat_file; else the newest
% *_results.mat under <run_dir>/memory_capacity; else the newest under
% data/memory_capacity. It ERRORS rather than silently plotting nothing, and the
% message says where it looked -- this figure's source file is gitignored and is
% frequently absent on a fresh clone.
%
% See also: run_memory_capacity, plot_memory_capacity, plot_memory_capacity_combined

arguments
    cfg.mat_file    (1,:) char    = ''
    cfg.run_dir     (1,:) char    = ''
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.preset_name (1,:) char    = ''   % unused; accepted for a uniform call
end

setup_paths();
out_dir      = default_out_dir(cfg.out_dir, mfilename('fullpath'));
project_root = fileparts(which('setup_paths'));

mat_file = resolve_mc_results(cfg.mat_file, cfg.run_dir, project_root);
fprintf('[fig_memory_capacity] source: %s\n', mat_file);

%% Reference figures (shown, not saved)
% Fig1 (total-MC + horizon distributions) and Fig2 (per-delay + cumulative) are
% the working views; only the combined strip below is a paper figure. Passing ''
% displays without saving -- plot_memory_capacity guards its save on ~isempty.
replot_memory_capacity(mat_file, '');

S = load(mat_file, 'results_all');
results_all = S.results_all;

%% The paper figure
fig3 = plot_memory_capacity_combined(results_all, '');
if ~cfg.visible; set(fig3, 'Visible', 'off'); end

%% Save
fig_tag = 'Fig_Memory_Capacity';
out = struct('figs', fig3, 'files', {{}}, 'source', mat_file);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig3);
    out.files = existing_outputs(out_dir, fig_tag);
    capture_git_provenance(out_dir, project_root);

    s = results_all.settings;
    write_figure_readme(out_dir, struct( ...
        'tag',    'fig_memory_capacity', ...
        'title',  'Stability_Manuscript figure: memory capacity', ...
        'script', 'fig_memory_capacity.m', ...
        'what',   ['A 1 x 3 strip: (a) cumulative memory capacity against ' ...
                   'delay, (b) per-delay R^2, (c) memory horizon paired across ' ...
                   'trials. All four adaptation conditions, with bootstrap 95% ' ...
                   'confidence bands.'], ...
        'how',    ['Presentation replot -- no simulation is re-run. A saved ' ...
                   'run_memory_capacity ensemble is loaded and ' ...
                   'plot_memory_capacity_combined assembles the strip from the ' ...
                   'same summary statistics the working figures use. Fig1 and ' ...
                   'Fig2 are shown for reference but not written to disk.'], ...
        'source', struct('mat_file', mat_file, 'run_tag', results_all.run_tag, ...
                         'preset', getfield_or(s, 'preset_name', '(pre-refactor run)'), ...
                         'run_mode', getfield_or(s, 'run_mode', '(pre-refactor run)')), ...
        'settings', mc_settings(s), ...
        'figures', {out.files}, ...
        'notes',  ['MODEL CLASS. These figures use SRNN_ESN_reservoir, which ' ...
                   'subclasses SRNNModel2, so the memory-capacity network is ' ...
                   'NOT the SRNNCellTypePairs network every other figure in the ' ...
                   'paper shows. That is a structural constraint -- the ESN ' ...
                   'readout has not been ported to the Pairs class -- and the ' ...
                   'methods section must say so.']));
end
end

%% ------------------------------------------------------------------------
function mat_file = resolve_mc_results(explicit, run_dir, project_root)
% Newest *_results.mat, from the most specific location that has one.
if ~isempty(explicit)
    if ~isfile(explicit)
        error('fig_memory_capacity:NoSuchFile', ...
            'mat_file does not exist:\n  %s', explicit);
    end
    mat_file = explicit;
    return
end

searched = {};
if ~isempty(run_dir); searched{end+1} = fullfile(run_dir, 'memory_capacity'); end
searched{end+1} = fullfile(project_root, 'data', 'memory_capacity');
searched{end+1} = fullfile(project_root, 'data', 'memory_capacity', 'paper_ready');

for k = 1:numel(searched)
    hits = dir(fullfile(searched{k}, '*_results.mat'));
    if ~isempty(hits)
        [~, newest] = max([hits.datenum]);
        mat_file = fullfile(hits(newest).folder, hits(newest).name);
        return
    end
end

error('fig_memory_capacity:NoResults', ...
    ['No *_results.mat found. Looked in:\n%s\n' ...
     'Run run_memory_capacity first (it is gitignored, so a fresh clone has none).'], ...
    sprintf('    %s\n', searched{:}));
end

function s = mc_settings(cfg)
% The subset of a saved MC settings struct worth putting in the README.
want = {'n', 'f', 'fs', 'level_of_chaos', 'tau_d', 'activation', 'S_c', ...
        'input_type', 'T_hold', 'readout_signal', 'n_trials', ...
        'T_train_sec', 'T_test_sec', 'd_max_sec', 'ode_solver', 'std_zero_floor'};
s = struct();
for k = 1:numel(want)
    if isfield(cfg, want{k}); s.(want{k}) = cfg.(want{k}); end
end
end

function v = getfield_or(s, name, default)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end
