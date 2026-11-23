%% TEST_LEADING_DIRECTION The uncoupled_reference stage, fig_leading_direction and the Figure 4 composite (Plan 2, package C).
%
% Tiny size, on a COPY of the three-condition regression run's sweeps
% (no_adaptation / sfa1_std1 / sfa3_std2, 100 neurons; found as in
% test_sweep_figures_conditions: this checkout, else the main checkout;
% SKIP if neither has it). Figures invisible.
%
%   1. uncoupled_reference_stage_plan: gains [0.01 0.05] and the repetition
%      table for every run mode of the mu7chain_AC bundle, the bundle's
%      conditions; overrides kept; an unknown field, a gain at or above the
%      sweep's lowest level, zero repetitions and a level_of_chaos in the
%      model overrides are refused (:BadConfig).
%   2. run_uncoupled_reference_stage at 100 neurons, 4 s, 2 repetitions:
%      the file, its fields, reference = the median at the first gain;
%      without adaptation the leading perturbation is all dendritic and the
%      rate is near -1/tau_d; with adaptation the rate is slower and the
%      perturbation is not dendritic.
%   3. fig_leading_direction: one panel per condition in the requested
%      order, as many points as successful runs, the zero line, NO line at
%      the uncoupled rate by default (Tom, 2026-11-21) but its numbers in the
%      stats and the table; with 'reference_line', true the dotted
%      line at the reference (symlog position), no text on the panels, the
%      table. A run without the stage: no dotted line, said in the table;
%      'require_reference' errors; 'reference_file' reads it. The column
%      layout. A run without the joint sample errors through
%      resolve_data_file.
%   4. fig_main4_connectivity (Plan 3): the six sweeps, lambda_1 by E:I
%      weight ratio and lambda_1 against rate per condition on the 9 in
%      canvas, one colour bar, letters A to C, no off-scale triangles, the
%      joint sample in one dot size, black medians and their key, files
%      and _table.md; fig_leading_direction is no longer part of it.
%
% About 2 to 3 minutes (the grouped Figure 4 rebuilds its networks).

n_fail = 0;
src = find_regression_run();
if isempty(src)
    fprintf('SKIP: no regression run directory (data/regression_fast/reference_post_refactor)\n');
    return
end
preset = 'regression_pairs_noStim_noise0p025_dualStd_3cond';
C = {'no_adaptation', 'sfa1_std1', 'sfa3_std2'};

%% 1. The plan
reps = {'essential', 2; 'fast', 4; 'fastmedium', 8; 'medium', 8; 'medium2', 8; 'production', 12};
ok = isequal(sort(reps(:, 1))', sort(run_mode_names()));
for m = 1:size(reps, 1)
    c = mu7chain_AC_config(reps{m, 1});
    pl = uncoupled_reference_stage_plan(c);
    want = bundle_run_conditions(c);
    ok = ok && isequal(pl.gains, [0.01 0.05]) && pl.n_reps == reps{m, 2} && strcmp(pl.folder, 'uncoupled_reference') && ...
        (isempty(want) || isequal(pl.conditions, want)) && strcmp(pl.preset_name, c.preset_name) && ...
        isequal(pl.ctx.model_defaults.level_of_chaos, 1);
end
n_fail = n_fail + check(ok, 'plan: gains [0.01 0.05], the repetition table, the bundle''s conditions, in every run mode');
tiny = struct('preset_name', preset, 'run_mode', 'fast', 'verbose', 'near-none');
tiny.uncoupled_reference = struct('gains', [0.02 0.1], 'n_reps', 2, ...
    'model', struct('n', 100, 'T_range', [0 4], 'lya_T_interval', [2 4], 'lya_warmup', 1));
pl = uncoupled_reference_stage_plan(tiny);
n_fail = n_fail + check(isequal(pl.gains, [0.02 0.1]) && pl.n_reps == 2 && isequal(pl.conditions, C) && ...
    isequal(pl.ctx.model_defaults.T_range, [0 4]) && pl.ctx.model_defaults.n == 100, 'plan: overrides kept');
bad = {struct('gain', 0.01), struct('gains', [0.01 0.25]), struct('gains', [0.01 0.01]), struct('n_reps', 0), ...
    struct('model', struct('level_of_chaos', 0.5))};
ok = true;
for k = 1:numel(bad)
    b = tiny; b.uncoupled_reference = bad{k};
    ok = ok && strcmp(error_of(@() uncoupled_reference_stage_plan(b)), 'uncoupled_reference_stage_plan:BadConfig');
end
b = tiny; b.run_mode = 'nonsense';
ok = ok && strcmp(error_of(@() uncoupled_reference_stage_plan(b)), 'uncoupled_reference_stage_plan:BadConfig');
n_fail = n_fail + check(ok, 'plan: unknown field, gain at the sweep''s lowest level, repeated gain, no repetitions, level_of_chaos in the model, unknown mode refused');

%% The run directory: a copy of the regression run's sweeps
tmp = tempname; mkdir(tmp);
cleanup_tmp = onCleanup(@() rmdir(tmp, 's'));
run_dir = fullfile(tmp, 'run'); mkdir(run_dir);
for pat = {'1D_sensitivity_*', 'param_space_*'}
    d = dir(fullfile(src, pat{1}));
    d = d([d.isdir]);
    for k = 1:numel(d); copyfile(fullfile(d(k).folder, d(k).name), fullfile(run_dir, d(k).name)); end
end
copyfile(fullfile(src, 'run_manifest.mat'), run_dir);
base = {'run_dir', run_dir, 'visible', false, 'verbose', 'near-none'};

%% 3a. The figure before the stage has run: no dotted line
out = fig_leading_direction(base{:}, 'out_dir', fullfile(tmp, 'f0'));
tbl = fileread(fullfile(tmp, 'f0', 'Fig_Leading_Direction_table.md'));
n_fail = n_fail + check(n_tagged(out.figs, 'uncoupled_reference') == 0 && n_tagged(out.figs, 'zero_growth_reference') == 3 && ...
    contains(tbl, 'No uncoupled reference in this run') && all(isnan([out.stats.reference])), ...
    'figure, run without the stage: no dotted line, the table says so');
close(out.figs);
n_fail = n_fail + check(strcmp(error_of(@() fig_leading_direction(base{:}, 'save', false, 'require_reference', true)), ...
    'fig_leading_direction:NoReference'), 'figure: require_reference errors without the stage');

%% 2. The stage
evalc('detail = run_uncoupled_reference_stage(tiny, run_dir, pl);');
ref_file = fullfile(run_dir, 'uncoupled_reference', 'uncoupled_reference_data.mat');
ok = isfile(ref_file) && contains(detail, ref_file);
if ok
    U = load(ref_file);
    ok = isequal({U.ref.name}, C) && isequal(U.settings.gains, [0.02 0.1]) && U.settings.n_reps == 2 && U.settings.n == 100;
    for k = 1:numel(U.ref)
        r = U.ref(k);
        ok = ok && isequal(size(r.lambda_1), [2 2]) && all(isfinite(r.lambda_1(:))) && isequal(r.n_ok, [2 2]) && ...
            abs(r.reference - median(r.lambda_1(:, 1))) < 1e-12 && r.lambda_1_min(1) <= r.reference && r.reference <= r.lambda_1_max(1);
    end
end
n_fail = n_fail + check(ok, 'stage: the file, a 2 x 2 table of lambda_1 per condition, reference = median at the first gain');
if ok
    m0 = build_from_preset(preset, 'no_adaptation', 'verbose', 'near-none'); tau_d = m0.tau_d;
    no = U.ref(1); ad = U.ref(2:3);
    fprintf('      reference lambda_1: %s (-1/tau_d = %g)\n', mat2str([U.ref.reference], 4), -1 / tau_d);
    n_fail = n_fail + check(all(no.lead_frac_x(:) == 1) && abs(no.reference + 1 / tau_d) < 0.25 / tau_d, ...
        'stage, no adaptation: all dendritic, the rate within a quarter of -1/tau_d');
    n_fail = n_fail + check(all([ad.reference] < 0) && all([ad.reference] > no.reference / 2) && ...
        all(arrayfun(@(r) median(r.lead_frac_x(:, 1)), ad) < 0.1), ...
        'stage, adapting: slower than the dendritic leak, the leading perturbation in the adaptive variables');
end

%% 3b. The figure with the reference
sel = {'sfa3_std2', 'no_adaptation'};
out = fig_leading_direction(base{:}, 'out_dir', fullfile(tmp, 'f1'), 'conditions', sel, 'reference_line', true);
ax = out.axes;
ok = numel(ax) == 2 && isequal(arrayfun(@(a) char(a.Title.String), ax, 'UniformOutput', false), ...
    grouped_condition_label(sel, 'sentence')) && isempty(findobj(out.figs, 'Type', 'text'));
for k = 1:2
    line = findall(ax(k), 'Tag', 'uncoupled_reference'); pts = findall(ax(k), 'Tag', 'leading_direction_points');
    want = U.ref(strcmp(C, sel{k})).reference;
    ok = ok && isscalar(line) && isscalar(pts) && numel(pts.XData) == out.stats(k).n && ...
        abs(line.XData(1) - symlog_lle(want, 0.05, 10)) < 1e-12 && abs(out.stats(k).reference - want) < 1e-12 && ...
        all(pts.YData >= 0 & pts.YData <= 1);
end
tbl = fileread(fullfile(tmp, 'f1', 'Fig_Leading_Direction_table.md'));
ok = ok && all(ismember({'Fig_Leading_Direction.png', 'Fig_Leading_Direction.fig'}, out.files)) && ...
    contains(tbl, 'Uncoupled reference') && contains(tbl, '0.02: ') && ~contains(tbl, 'No uncoupled reference');
n_fail = n_fail + check(ok, ['figure (reference_line on): the selection in order, a point per run, the dotted line at the reference, ' ...
    'no text on the panels, the table']);
close(out.figs);
out = fig_leading_direction('run_dir', tmp, 'data_file', fullfile(dir_of(run_dir, 'param_space_*'), 'psa_object.mat'), ...
    'reference_file', ref_file, 'save', false, 'visible', false, 'verbose', 'near-none', 'layout', 'column');
p = reshape([out.axes.Position], 4, [])';
n_fail = n_fail + check(n_tagged(out.figs, 'uncoupled_reference') == 0 && all(isfinite([out.stats.reference])) && ...
    contains(out.table_text, 'Uncoupled reference') && all(diff(p(:, 2)) < 0) && ...
    max(abs(diff(p(:, 1)))) < 1e-6, ['figure, default (no line, Tom 2026-11-21): explicit joint sample and reference file, ' ...
    'the reference in the stats and the table, no line; the column layout stacks the panels']);
close(out.figs);
empty_run = fullfile(tmp, 'empty'); mkdir(empty_run);
n_fail = n_fail + check(startsWith(error_of(@() fig_leading_direction('run_dir', empty_run, 'save', false, 'visible', false)), ...
    'resolve_data_file:'), 'figure: a run without the joint sample errors through resolve_data_file');
n_fail = n_fail + check(strcmp(error_of(@() fig_leading_direction(base{:}, 'save', false, 'conditions', {'sfa3_stdc2'})), ...
    'select_conditions:NoSuchCondition'), 'figure: a condition the run lacks -> NoSuchCondition');

%% 4. The composite (Plan 3: the sweeps, lambda_1 by E:I weight ratio, lambda_1 against rate)
st = manuscript_style();
evalc('m = fig_main4_connectivity(base{:}, ''preset_name'', preset, ''out_dir'', fullfile(tmp, ''m4''));');
m_px = getpixelposition(m.figs);
letters = sort(arrayfun(@(t) t.String(2), findall(m.figs, 'Tag', 'composite_panel_letter'))');
med = findall(m.figs, 'Tag', 'rate_bin_median');
scat = findall(m.figs, 'Type', 'scatter');
sz = arrayfun(@(h) unique(h.SizeData), scat, 'UniformOutput', false);
lgC = findall(m.figs, 'Tag', 'main4_legend_C');
ok = numel(findall(m.figs, 'Type', 'axes')) == 6 + 2 * numel(C) && isscalar(findall(m.figs, 'Type', 'colorbar')) && ...
    isequal(letters, 'ABC') && n_tagged(m.figs, 'leading_direction_points') == 0 && ...       % the old (C) has left
    isempty(findall(m.figs, 'Tag', 'offscale_marker')) && ...                                % no triangles
    numel(med) == numel(C) && all(arrayfun(@(h) isequal(h.Color, [0 0 0]), med)) && ...      % black medians
    numel(scat) == numel(C) && all(cellfun(@isscalar, sz)) && isscalar(unique([sz{:}])) && ...   % one dot size
    isscalar(lgC) && isequal(lgC.String, {'Median per rate bin'}) && ...
    abs(m_px(3) / get(groot, 'ScreenPixelsPerInch') - st.canvas_w) < 0.01 && ...
    all(ismember({'Fig_Main4_Connectivity.png', 'Fig_Main4_Connectivity.fig'}, m.files)) && ...
    isfile(fullfile(tmp, 'm4', 'Fig_Main4_Connectivity_table.md'));
if ok
    % The median's key sits above the multiple-timescale panel (the last here).
    lgC.Units = 'pixels'; own = lgC.Axes; own.Units = 'pixels';
    ok = lgC.Position(2) >= own.Position(2) + own.Position(4) && contains(own.Title.String{1}, 'Multiple') && ...
        contains(fileread(fullfile(tmp, 'm4', 'Fig_Main4_Connectivity_table.md')), '## Panel C');
end
n_fail = n_fail + check(ok, ['composite: six sweeps, lambda_1 by E:I ratio and lambda_1 against rate per condition, one color bar, ' ...
    'letters A to C, no triangles, the joint sample in one dot size, black medians with their key above the multiple-timescale panel, table']);
close(m.figs);
evalc('m = fig_main4_connectivity(base{:}, ''preset_name'', preset, ''out_dir'', fullfile(tmp, ''m4''), ''fig_tag'', ''Fig_Main4_Grouped'', ''conditions'', sel);');
n_fail = n_fail + check(any(strcmp(m.files, 'Fig_Main4_Grouped.png')) && numel(findall(m.figs, 'Type', 'scatter')) == 2, ...
    'composite: the registry''s file name; a selection reaches (B) and (C)');
close(m.figs);

if n_fail == 0
    fprintf('test_leading_direction: all checks passed\n');
else
    error('test_leading_direction:Failed', '%d check(s) FAILED', n_fail);
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

function p = dir_of(parent, pattern)
d = dir(fullfile(parent, pattern));
d = d([d.isdir]);
p = fullfile(d(1).folder, d(1).name);
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

function n = n_tagged(fig, tag)
n = numel(findall(fig, 'Tag', tag));
end
