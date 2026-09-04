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

% INSIDE THE RUN DIRECTORY ONLY. This used to fall back to
% data/memory_capacity/ and then to a paper_ready/ subfolder that has never
% existed. On 2026-09-03 that fallback made this figure plot a .mat from
% 2026-08-22 -- a different network -- and report success, because the run's
% memory_capacity stage had failed and left nothing. To plot a standalone
% analysis, pass its location AS run_dir (fullfile(project_root, 'data')), or
% pass mat_file directly.
mat_file = resolve_data_file(cfg.mat_file, cfg.run_dir, ...
    {fullfile(cfg.run_dir, 'memory_capacity')}, ...
    '*_results.mat', ...
    'Run run_memory_capacity first');
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

end
end

%% ------------------------------------------------------------------------
% resolve_mc_results lived here. It was the one figure that resolved its data
% correctly -- run directory first, then data/ -- and so became the model for
% _common/resolve_data_file.m when fig_eig_heatmap and
% fig_memory_capacity_example were fixed to stop reading a hardcoded sibling.
% Replaced by that shared version 2026-09-01 rather than leaving one pattern in
% three near-identical copies.






