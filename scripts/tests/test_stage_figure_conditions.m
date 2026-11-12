% TEST_STAGE_FIGURE_CONDITIONS cfg.conditions on the stage-data figures (FD-006), on synthetic data.
%
% No simulation and no run directory: the test writes small synthetic stage
% files into a temporary run directory and builds three figures from them,
% each invisible and unsaved except where the table is the thing checked.
%   1. fig_dc_lle: the legend follows the REQUESTED order; empty selection
%      keeps the run's order; a name the run lacks errors
%      (select_conditions:NoSuchCondition).
%   2. fig_memory_capacity: with three of four conditions selected, the
%      table lists the three pairs among them and says six were run.
%   3. fig_single_neuron_revised: the mechanism rule of its header for
%      network names (both network selections map to all four columns), exact
%      selection for single-neuron names ({no_adaptation} draws one; a
%      reordered pair draws in the order given); the titles come from the
%      condition NAME, not the saved title.
%   4. titles_by_name: map first, then the saved title, then the name.
% Runs in a few seconds.

tmp = tempname; mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's'));
names4 = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2', 'sfa3_std2'};
set3   = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
n_fail = 0;

%% 1. fig_dc_lle
mkdir(fullfile(tmp, 'dc_lle'));
dc_lle_results = struct('conditions', {names4}, 'dc_levels', (0:4)', ...
    'level_mean', repmat(1:4, 5, 1), 'level_std', 0.1 * ones(5, 4), 'n_seeds', 2, ...
    'config', struct('preset_name', 'synthetic'));
save(fullfile(tmp, 'dc_lle', 'dc_lle_results.mat'), 'dc_lle_results');
o = fig_dc_lle('run_dir', tmp, 'save', false, 'visible', false, 'conditions', {'sfa3_std2', 'no_adaptation'});
lg = findobj(o.figs, 'Type', 'legend');
n_fail = n_fail + check('dc_lle legend in requested order', ...
    isequal(lg.String, titles_by_name({'sfa3_std2', 'no_adaptation'})));
close(o.figs);
o = fig_dc_lle('run_dir', tmp, 'save', false, 'visible', false);
lg = findobj(o.figs, 'Type', 'legend');
n_fail = n_fail + check('dc_lle no selection = run order', isequal(lg.String, titles_by_name(names4)));
close(o.figs);
try
    fig_dc_lle('run_dir', tmp, 'save', false, 'visible', false, 'conditions', {'sfa9_std9'});
    n_fail = n_fail + check('dc_lle unknown condition errors', false);
catch ME
    n_fail = n_fail + check('dc_lle unknown condition errors', strcmp(ME.identifier, 'select_conditions:NoSuchCondition'));
end

%% 2. fig_memory_capacity
mkdir(fullfile(tmp, 'memory_capacity'));
nd = 5; nt = 4;
pairs = nchoosek(1:4, 2);
stats = struct('pair', {}, 'mean_diff', {}, 'p_perm', {}, 'n_patterns', {}, 'exact', {}, 'cohens_dz', {});
for p = 1:size(pairs, 1)
    stats(p) = struct('pair', sprintf('%s vs %s', names4{pairs(p, 1)}, names4{pairs(p, 2)}), ...
        'mean_diff', p, 'p_perm', 0.5, 'n_patterns', 16, 'exact', true, 'cohens_dz', 0.1);
end
ci = struct('lo', zeros(4, nd), 'hi', ones(4, nd));
summary = struct('R2_mean', rand(4, nd), 'R2_ci95', ci, 'MC_mean', 1:4, ...
    'MC_ci95', struct('lo', 0:3, 'hi', 2:5), 'H_mean', 1:4, 'H_ci95', struct('lo', 0:3, 'hi', 2:5), ...
    'stats', stats);
settings = struct('delay_s', 0.3 * (1:nd), 'R2_threshold_for_horizon', 0.1, 'preset_name', 'synthetic', ...
    'run_mode', 'fast', 'readout_signal', 'rate', 'n_boot', 10);
results_all = struct('conditions', {names4}, 'H_trials', rand(nt, 4), 'MC_trials', rand(nt, 4), ...
    'summary', summary, 'settings', settings); %#ok<NASGU>
save(fullfile(tmp, 'memory_capacity', 'synthetic_results.mat'), 'results_all');
od = fullfile(tmp, 'mc_out');
o = fig_memory_capacity('run_dir', tmp, 'out_dir', od, 'save', true, 'visible', false, 'conditions', set3);
close(o.figs);
txt = fileread(fullfile(od, 'Fig_Memory_Capacity_table.md'));
n_rows = numel(regexp(txt, '\| \w+ vs \w+ \|', 'match'));
n_fail = n_fail + check('MC table: three pairs among the selection', n_rows == 3 && ~contains(txt, 'vs sfa3_std2'));
n_fail = n_fail + check('MC table: says six were run', contains(txt, 'Paired tests: 6 were run'));

%% 3. fig_single_neuron_revised
mkdir(fullfile(tmp, 'single_neuron'));
sn = {'no_adaptation', 'sfa1_std0', 'sfa0_std1', 'sfa0_stdc2'};
model = {'none', 'none', 'product', 'chain'};
t = linspace(0, 30, 50);
for c = 1:4
    r = struct('name', sn{c}, 'title', 'SAVED TITLE', 't', t, 'u', t, 'x', t, 'r', t, 'syn', t, ...
        'sfa', t, 'std', t, 'std_model', model{c}, 'n', 1);
    if c == 1, results = r; else, results(c) = r; end %#ok<AGROW>
end
settings = struct('display_window', [0 30]); %#ok<NASGU>
save(fullfile(tmp, 'single_neuron', 'single_neuron_data.mat'), 'results', 'settings');
cases = {set3, 4; names4, 4; {'no_adaptation'}, 1; {}, 4; {'no_adaptation', 'sfa3_std0'}, 2};
for k = 1:size(cases, 1)
    o = fig_single_neuron_revised('run_dir', tmp, 'save', false, 'conditions', cases{k, 1});
    ttl = titled_axes(o.figs);
    n_fail = n_fail + check(sprintf('single neuron: %d column(s) for {%s}', cases{k, 2}, strjoin(cases{k, 1}, ', ')), ...
        numel(ttl) == cases{k, 2} && ~any(strcmp(ttl, 'SAVED TITLE')));
    close(o.figs);
end
% Single-neuron names are taken as names, in the order given (the bundle's
% per-set lists): two columns, chain first.
o = fig_single_neuron_revised('run_dir', tmp, 'save', false, 'conditions', {'sfa0_stdc2', 'no_adaptation'});
ttl = titled_axes(o.figs);
map = srnn_condition_titles();
n_fail = n_fail + check('single neuron: exact names in the order given', ...
    isequal(reshape(ttl, 1, []), {map('sfa0_stdc2'), map('no_adaptation')}));
close(o.figs);

%% 4. titles_by_name
tt = titles_by_name({'sfa3_std2', 'not_a_name', 'also_not'}, {'old', 'Saved', ''});
map = srnn_condition_titles();
n_fail = n_fail + check('titles_by_name fallbacks', isequal(tt, {map('sfa3_std2'), 'Saved', 'also_not'}));

if n_fail > 0
    error('test_stage_figure_conditions:Failed', '%d check(s) failed.', n_fail);
end
fprintf('test_stage_figure_conditions: all checks PASS\n');

%% ------------------------------------------------------------------------
function n = check(what, ok)
if ok, fprintf('  PASS  %s\n', what); n = 0; else, fprintf('  FAIL  %s\n', what); n = 1; end
end

function s = titled_axes(fig)
% Titles of the axes that carry one, left to right.
ax = findobj(fig, 'Type', 'axes');
ax = ax(arrayfun(@(a) ~isempty(a.Title.String), ax));
[~, o] = sort(arrayfun(@(a) a.Position(1), ax));
s = arrayfun(@(a) char(a.Title.String), ax(o), 'UniformOutput', false);
end
