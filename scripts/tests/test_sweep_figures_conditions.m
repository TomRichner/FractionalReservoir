%% TEST_SWEEP_FIGURES_CONDITIONS The 'conditions' option on the sweep figures and the tau sweep's figure (FD-006, WP1b).
%
% Against the three-condition regression run (no_adaptation / sfa1_std1 /
% sfa3_std2; its tau sweep holds sfa3_std2 alone), COPIED to a temporary
% directory first, because fig_param_space_allStd writes its prep folder
% into the run directory it is given:
%
%   fig_sensitivity_medians    no selection: the preset's conditions in the
%                              legend; a subset: exactly it, requested order;
%                              a name the run lacks: NoSuchCondition
%   fig_param_space_allStd     a subset: the columns, requested order
%   fig_EI_weights_param_space a subset: the columns, requested order
%   fig_local_vs_finite_lle    a subset: the panels, requested order, the
%                              table says so; a missing name: NoSuchCondition
%   fig_sfa_EOC_allStd         the INTERSECTION rule: a selection holding
%                              conditions the tau sweep lacks draws the ones
%                              it holds (one panel, no title, as before) and
%                              the table names the rest; a selection holding
%                              none of them is NoSelectedCondition
%   write_run_parameters_md    every condition's own synapse_config, the
%                              one-pool and the product described as such
%
% fig_sensitivity_analysis_allStd (about four minutes at regression size)
% is not run here; its selection goes through replot_sensitivity, exactly as
% fig_param_space_allStd's goes through replot_param_space_analysis.
%
% The run is found as in test_select_conditions_figure (this checkout, else
% the main checkout); SKIP if neither has it. About 90 s. Figures invisible.

n_fail = 0;

src = find_regression_run();
if isempty(src)
    fprintf('SKIP: no regression run directory (data/regression_fast/reference_post_refactor)\n');
    return
end
st  = manuscript_style();
tmp = tempname; mkdir(tmp);
cleanup_tmp = onCleanup(@() rmdir(tmp, 's'));
run_dir = fullfile(tmp, 'run');
mkdir(run_dir);
for pat = {'1D_sensitivity_*', 'param_space_*', 'tau_sensitivity_*'}
    d = dir(fullfile(src, pat{1}));
    d = d([d.isdir]);
    for k = 1:numel(d); copyfile(fullfile(d(k).folder, d(k).name), fullfile(run_dir, d(k).name)); end
end
copyfile(fullfile(src, 'run_manifest.mat'), run_dir);
preset = 'regression_pairs_noStim_noise0p025_dualStd_3cond';
base = {'run_dir', run_dir, 'visible', false, 'verbose', 'near-none'};
sel  = {'sfa3_std2', 'no_adaptation'};
want_sel = cellfun(@(c) st.condition_title(c), sel, 'UniformOutput', false);

%% fig_sensitivity_medians
out = fig_sensitivity_medians(base{:}, 'preset_name', preset, 'save', false);
n_fail = n_fail + check(isequal(legend_text(out.figs(1)), ...
    cellfun(@(c) st.condition_title(c), {'no_adaptation', 'sfa1_std1', 'sfa3_std2'}, 'UniformOutput', false)), ...
    'medians, no selection: the preset''s three, in its order');
close(out.figs);
out = fig_sensitivity_medians(base{:}, 'preset_name', preset, 'save', false, 'conditions', sel);
n_fail = n_fail + check(isequal(legend_text(out.figs(1)), want_sel), ...
    'medians, subset: exactly the selection, in the requested order');
close(out.figs);
id = error_of(@() fig_sensitivity_medians(base{:}, 'preset_name', preset, 'save', false, ...
    'conditions', {'no_adaptation', 'sfa3_stdc2'}));
n_fail = n_fail + check(strcmp(id, 'select_conditions:NoSuchCondition'), ...
    'medians: a condition the run lacks -> NoSuchCondition');

%% fig_param_space_allStd
out = fig_param_space_allStd(base{:}, 'preset_name', preset, 'save', true, ...
    'out_dir', fullfile(tmp, 'ps'), 'conditions', sel);
n_fail = n_fail + check(isequal(column_titles(out.figs(1)), want_sel), ...
    'param space, subset: the columns are the selection, in order');
n_fail = n_fail + check(isempty(dir(fullfile(run_dir, 'replot_param_space_*'))), ...
    'param space: its prep folder is removed');
close(out.figs);

%% fig_EI_weights_param_space
out = fig_EI_weights_param_space(base{:}, 'preset_name', preset, 'save', false, 'conditions', sel);
n_fail = n_fail + check(isequal(column_titles(out.figs(1)), want_sel), ...
    'E:I weights, subset: the columns are the selection, in order');
close(out.figs);

%% fig_local_vs_finite_lle
out = fig_local_vs_finite_lle(base{:}, 'save', true, 'out_dir', fullfile(tmp, 'lvf'), 'conditions', sel);
n_fail = n_fail + check(isequal(column_titles(out.figs(1)), want_sel), ...
    'local vs finite, subset: the panels are the selection, in order');
tbl = fileread(fullfile(tmp, 'lvf', 'Fig_Local_vs_Finite_LLE_table.md'));
rows = regexp(tbl, '^\| (\w+) \| near-default', 'tokens', 'lineanchors');
rows = cellfun(@(c) c{1}, rows, 'UniformOutput', false);
n_fail = n_fail + check(isequal(rows, sel) && contains(tbl, 'Conditions drawn (selected, in this order): sfa3_std2, no_adaptation'), ...
    'local vs finite: table rows in the requested order, selection stated');
close(out.figs);
id = error_of(@() fig_local_vs_finite_lle(base{:}, 'save', false, 'conditions', {'sfa3_stdc2'}));
n_fail = n_fail + check(strcmp(id, 'select_conditions:NoSuchCondition'), ...
    'local vs finite: a condition the run lacks -> NoSuchCondition');

%% fig_sfa_EOC_allStd: the intersection rule
out = fig_sfa_EOC_allStd(base{:}, 'preset_name', preset, 'save', true, ...
    'out_dir', fullfile(tmp, 'eoc'), 'conditions', {'no_adaptation', 'sfa1_std1', 'sfa3_std2'});
ax = findobj(out.figs(1), 'Type', 'axes');
n_fail = n_fail + check(isscalar(ax) && isempty(ax.Title.String), ...
    'EOC: three selected, the sweep holds one -> one panel, untitled as before');
tbl = fileread(fullfile(tmp, 'eoc', 'Fig_SFA_EOC_allStd_table.md'));
n_fail = n_fail + check(contains(tbl, 'condition `sfa3_std2`') && ...
    contains(tbl, 'so not drawn: no_adaptation, sfa1_std1.'), ...
    'EOC: the table names the drawn condition and the selected ones not in the sweep');
close(out.figs(isgraphics(out.figs)));
id = error_of(@() fig_sfa_EOC_allStd(base{:}, 'preset_name', preset, 'save', false, ...
    'conditions', {'no_adaptation'}));
n_fail = n_fail + check(strcmp(id, 'fig_sfa_EOC_allStd:NoSelectedCondition'), ...
    'EOC: a selection the sweep holds none of -> NoSelectedCondition');

%% write_run_parameters_md: every condition's own synapse_config
md = fileread(write_run_parameters_md(run_dir));
n_fail = n_fail + check(~contains(md, 'identical for every'), ...
    'parameters.md: no claim that the synapse settings are shared');
n_fail = n_fail + check(contains(md, '`synapse_config` for `sfa1_std1`:') && ...
    contains(md, '`synapse_config` for `sfa3_std2`:'), ...
    'parameters.md: a synapse_config block for each depressing condition');
n_fail = n_fail + check(contains(md, 'one depression pool: `tau_rec` = ') && ...
    contains(md, 'product of 2 independent depression pools'), ...
    'parameters.md: one pool and product described as such');

if n_fail == 0
    fprintf('test_sweep_figures_conditions: all checks passed\n');
else
    error('test_sweep_figures_conditions:Failed', '%d check(s) FAILED', n_fail);
end

%% ------------------------------------------------------------------------
function d = find_regression_run()
% <project root>/data/..., else the main checkout's (a worktree's data/ is empty).
rel = fullfile('data', 'regression_fast', 'reference_post_refactor');
root = fileparts(which('setup_paths'));
d = fullfile(root, rel);
if isfolder(d); return; end
[s, common] = system(sprintf('git -C "%s" rev-parse --path-format=absolute --git-common-dir', root));
d = '';
if s == 0
    main_root = fileparts(strtrim(common));   % <main checkout>/.git -> <main checkout>
    if isfolder(fullfile(main_root, rel)); d = fullfile(main_root, rel); end
end
end

function f = check(ok, what)
if ok
    fprintf('PASS  %s\n', what);
    f = 0;
else
    fprintf('FAIL  %s\n', what);
    f = 1;
end
end

function id = error_of(fn)
id = '';
try
    fn();
catch ME
    id = ME.identifier;
end
end

function t = legend_text(fig)
lg = findobj(fig, 'Type', 'legend');
t = reshape(cellstr(lg(1).String), 1, []);
end

function t = column_titles(fig)
% The non-empty axes titles, left to right (one per column).
ax = findobj(fig, 'Type', 'axes');
ax = ax(arrayfun(@(a) ~isempty(a.Title.String), ax));
p = reshape([ax.Position], 4, [])';
[~, o] = sort(p(:, 1));
t = arrayfun(@(a) char(strjoin(cellstr(a.Title.String), ' ')), ax(o), 'UniformOutput', false)';
end
