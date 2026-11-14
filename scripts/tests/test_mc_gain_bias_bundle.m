%% test_mc_gain_bias_bundle.m
% Block C from the two entry points (FD-010, FD-024): the opt-in stage
% 'mc_gain_sweep', the bundle mu7chain_gainbias_config with its two grids,
% and the bundle of blocks A and C, mu7chain_AC_config.
%
% Checks:
%   1. The opt-in stage leaves every existing bundle's plan as it was. The
%      plans of mu7chain_config('essential'), ('fast'), mu7revisedMedium_config
%      and regression_fast_config equal the literals below, read off
%      paper_stage_plan at c9b9377 (before opt-in stages existed), field for
%      field, and none mentions 'mc_gain_sweep'.
%   2. The opt-in stage runs only when named, after the stages of the first
%      list; a typo is still refused.
%   3. The bundle: its plan is exactly {'mc_gain_sweep'} (every other stage
%      skipped, no sweep), its directories, registry (no file tag a prefix
%      of another), conditions per set, grids and variant networks, and
%      that it is mu7chain_config's otherwise.
%   4. The stage's plan in EVERY run mode: the gain x bias grid, its two
%      variants (gain only, on its gains, at their own bias) and the gain x
%      self-coupling grid, each in its own folder, with the size table's
%      levels and seeds (3 x 3 x 1, 7 x 7 x 3, 13 x 13 x 10, 17 x 17 x 10),
%      the gains log-spaced on the bias grid and evenly spaced on the
%      self-coupling grid, zero a self-coupling level, and the job count.
%   5. The form with one unnamed grid still plans as it did (the main run
%      into the run directory, variants into mc_gain_sweep_<tag>); bad
%      settings are refused before the run directory exists.
%   6. mu7chain_AC_config: block A's stages then block C's, both registries,
%      its own directories, the grouped figures pointed at its figure root.
%   7. END TO END on a tiny case: run_all_paper_analyses on the bundle
%      ('essential', 3cond) with a 50-neuron network, the shortened protocol
%      of test_mc_gain_sweep, 2 x 2 levels per grid, one seed and three
%      conditions; then make_all_paper_figures from that run directory.
%      Both grids' results and the two variants' are where the figures look,
%      every registry entry writes files under its own tag, the gain curves
%      have two columns and a row per network, and a grid or variant the
%      run directory does not hold is an error, not a search elsewhere.
%   8. A run directory laid out as before grids had tags (one unnamed grid)
%      still reads through the resolver and resumes through the stage with
%      nothing left to run.
%
% Asserts; prints one PASS line. About three minutes serial, less with a
% pool. Assumes setup_paths has run and a pool allocation is set
% (pool_set(0) will do): both entry points refuse without one.

%% 1. Existing bundles' plans are unchanged
all10 = {'sweeps', 'memory_capacity', 'mc_example', 'eig_heatmap', 'numerics', ...
    'lyapunov_spectrum', 'transient_gain', 'local_lyapunov', 'dc_lle', 'paper_illustrations'};
sw3 = {'sensitivity', 'tau', 'param_space'};
full = plan_literal(all10, cell(1, 0), sw3, cell(1, 0), false, all10);
ess_stages = {'sweeps', 'eig_heatmap', 'local_lyapunov', 'paper_illustrations'};
ess  = plan_literal(ess_stages, ...
    {'memory_capacity', 'mc_example', 'numerics', 'lyapunov_spectrum', 'transient_gain', 'dc_lle'}, ...
    {'sensitivity', 'param_space'}, {'tau'}, true, all10);
bundles = {'mu7chain_config(''essential'')', mu7chain_config('essential'), ess; ...
           'mu7chain_config(''fast'')',      mu7chain_config('fast'),      full; ...
           'mu7revisedMedium_config',        mu7revisedMedium_config(),    full; ...
           'regression_fast_config',         regression_fast_config(),     full};
for k = 1:size(bundles, 1)
    pl = paper_stage_plan(bundles{k, 2});
    assert(isequal(pl, bundles{k, 3}), '%s: the plan changed', bundles{k, 1});
    assert(~ismember('mc_gain_sweep', [pl.stages, pl.stages_skipped]) && ...
        ~isfield(pl.runs, 'mc_gain_sweep'), '%s mentions the opt-in stage', bundles{k, 1});
end

%% 2. The opt-in stage runs only when named
[st, ~, opt_in] = paper_stage_names();
assert(isequal(st, all10) && isequal(opt_in, {'mc_gain_sweep'}));
pl = paper_stage_plan(struct('illustrations', true, 'stages', {{'mc_gain_sweep', 'dc_lle'}}));
assert(isequal(pl.stages, {'dc_lle', 'mc_gain_sweep'}), 'opt-in stages run after the first list');
assert(isequal(pl.stages_skipped, setdiff(all10, {'dc_lle'}, 'stable')) && pl.runs.mc_gain_sweep);
assert(throws(@() paper_stage_plan(struct('stages', {{'mc_gain_swep'}})), 'paper_stage_plan:UnknownStage'));

%% 3. The bundle
root = fileparts(which('setup_paths'));
conds3 = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
conds4 = [conds3, {'sfa3_std2'}];
variant_tags = {'EI3to2', 'EI2to3'};
entry_names = {'fig_mc_maps_bias', 'fig_mc_maps_selfcoupling', 'fig_mc_gain_curves', ...
    'fig_mc_robustness_bias', 'fig_mc_robustness_selfcoupling', 'fig_mc_vs_lambda1', ...
    'fig_lambda1_vs_bias', 'fig_rate_vs_lambda1', 'fig_slope_share_vs_gain', ...
    'mc_statistics_bias', 'mc_statistics_selfcoupling'};
entry_grids = {'bias', 'selfcoupling', 'bias', 'bias', 'selfcoupling', 'selfcoupling', ...
    'bias', 'bias', 'bias', 'bias', 'selfcoupling'};
for md = {'essential', 'fast'}
    for fs = {'3cond', '4cond'}
        c = mu7chain_gainbias_config(md{1}, fs{1});
        base = mu7chain_config(md{1}, fs{1});
        pl = paper_stage_plan(c);
        assert(isequal(pl.stages, {'mc_gain_sweep'}) && isequal(pl.stages_skipped, all10) && ...
            isempty(pl.sweeps) && isequal(pl.sweeps_skipped, sw3) && pl.is_selection);
        assert(strcmp(c.run_dir, ['data/mu7chain_gainbias_' md{1}]));
        assert(strcmp(c.fig_root, ['figs/mu7chain_gainbias_' md{1} '_' fs{1}]));
        names = cellfun(@(f) f.name, c.figures, 'UniformOutput', false);
        assert(isequal(names, entry_names), 'the registry is: %s', strjoin(names, ', '));
        assert(isequal(cellfun(@(f) arg(f.args, 'grid'), c.figures, 'UniformOutput', false), entry_grids), ...
            'every entry names the grid it reads');
        if strcmp(fs{1}, '3cond'); want = conds3; else; want = conds4; end
        assert(isequal(c.conditions, want));
        assert(isequal(c.figure_sets, base.figure_sets));
        grids = c.mc_gain_sweep.grids;
        assert(isequal({grids.tag}, {'bias', 'selfcoupling'}) && ...
            isequal({grids.second_axis}, {'S_c', 'self_coupling'}) && ...
            isequal({grids.gain_spacing}, {'log', 'linear'}));
        assert(isequal({grids(1).variants.tag}, variant_tags) && isempty(grids(2).variants));
        assert(isequal({grids(1).variants.preset}, strcat(base.mc_preset, '_', variant_tags)));
        curves = c.figures{strcmp(names, 'fig_mc_gain_curves')};
        assert(isequal(arg(curves.args, 'variants'), variant_tags) && ~arg(curves.args, 'mc_vs_lambda1'), ...
            'the gain-curves entry names the variants the stage runs, and draws two columns');
        maps = c.figures(startsWith(names, 'fig_mc_maps_'));
        assert(all(cellfun(@(f) strcmp(arg(f.args, 'mc_scale'), 'linear') && ...
            strcmp(arg(f.args, 'lle_scale'), 'linear') && arg(f.args, 'lle_clip') == 5, maps)), ...
            'the registered maps have linear color scales, lambda_1 clipped at 5');
        scatter = c.figures{strcmp(names, 'fig_mc_vs_lambda1')};
        assert(strcmp(arg(scatter.args, 'size_by'), 'self_coupling'));
        % save_figure_stable deletes <tag>* before saving: no tag may be a
        % prefix of another, wherever the entries are saved.
        tags = cellfun(@(f) file_tag(f.args), c.figures, 'UniformOutput', false);
        for a = 1:numel(tags)
            others = tags([1:a-1, a+1:end]);
            assert(~any(startsWith(others, tags{a})), 'the tag %s is a prefix of another', tags{a});
        end
        differs = {'stages', 'sweeps', 'eig_heatmap_examples', 'mc_gain_sweep', 'run_dir', ...
            'fig_root', 'figures'};
        assert(isequal(rmfield(c, differs), rmfield(base, intersect(differs, fieldnames(base)))), ...
            'the bundle differs from mu7chain_config in more than its stage, directories and registry');
    end
end
assert(throws(@() mu7chain_gainbias_config('fast2'), 'mu7chain_gainbias_config:BadRunMode'));
assert(strcmp(mu7chain_gainbias_config().run_mode, 'essential'), 'the default is the size agents run');

%% 4. The stage's plan in every run mode: both grids, the variants, the sizes
sizes = {'essential', 3, 1; 'fast', 7, 3; 'medium', 13, 10; 'medium2', 13, 10; 'production', 17, 10};
assert(isequal(sizes(:, 1)', run_mode_names()), 'the size table below must cover run_mode_names');
for m = 1:size(sizes, 1)
    [md, n, n_seeds] = sizes{m, :};
    c = mu7chain_gainbias_config(md);
    pl = mc_gain_sweep_stage_plan(c);
    assert(isequal({pl.runs.label}, {'bias', 'bias_EI3to2', 'bias_EI2to3', 'selfcoupling'}), md);
    assert(isequal({pl.runs.folder}, {'mc_gain_sweep_bias', 'mc_gain_sweep_bias_EI3to2', ...
        'mc_gain_sweep_bias_EI2to3', 'mc_gain_sweep_selfcoupling'}), md);
    seeds = 3000 + (1:n_seeds);
    log_gains = round(logspace(-1, log10(4), n), 10);
    a = pl.runs(1).args;
    assert(strcmp(arg(a, 'preset_name'), c.mc_preset) && strcmp(arg(a, 'run_mode'), md) && ...
        strcmp(arg(a, 'ode_solver'), c.mc_ode_solver) && strcmp(arg(a, 'second_axis'), 'S_c'));
    assert(isequal(arg(a, 'gain_levels'), log_gains) && ...
        isequal(arg(a, 'bias_levels'), round(linspace(-0.6, 0.6, n), 10)) && ...
        isequal(arg(a, 'seeds'), seeds), '%s: the gain x bias grid', md);
    assert(strcmp(mc_gain_axis_scale(arg(a, 'gain_levels')), 'log'));
    for k = 2:3
        a = pl.runs(k).args;
        assert(strcmp(arg(a, 'preset_name'), [c.mc_preset '_' variant_tags{k - 1}]) && ...
            isequal(arg(a, 'gain_levels'), log_gains) && isequal(arg(a, 'bias_levels'), 0.2) && ...
            isequal(arg(a, 'seeds'), seeds), ...
            '%s, %s: gain only, on the bias grid''s gains, at its own bias', md, pl.runs(k).label);
    end
    a = pl.runs(4).args;
    s_levels = arg(a, 'second_levels');
    assert(strcmp(arg(a, 'preset_name'), c.mc_preset) && strcmp(arg(a, 'second_axis'), 'self_coupling'));
    assert(isequal(arg(a, 'gain_levels'), round(linspace(0.1, 4, n), 10)) && ...
        isequal(s_levels, round(linspace(-4, 4, n), 10)) && isequal(arg(a, 'seeds'), seeds), ...
        '%s: the gain x self-coupling grid', md);
    assert(strcmp(mc_gain_axis_scale(arg(a, 'gain_levels')), 'linear') && nnz(s_levels == 0) == 1 && ...
        s_levels(1) == -4 && s_levels(end) == 4, '%s: evenly spaced gains; s = 0 is a level', md);
    assert(isequal([pl.runs.n_jobs], [n^2, n, n, n^2] * n_seeds) && ...
        pl.n_jobs == (2 * n^2 + 2 * n) * n_seeds && abs(pl.worker_minutes - 9.5 * pl.n_jobs) < 1e-9, ...
        '%s: job counts', md);
end
% The level that reads 0.2 is 0.2 (linspace gives 0.19999999999999998).
G = mc_gain_sweep_grid('fast');
assert(any(G.second_levels == 0.2) && G.n_jobs == 147);

%% 5. One unnamed grid, as before grids had tags; and what is refused
c = mu7chain_gainbias_config('essential');
variants = c.mc_gain_sweep.grids(1).variants;
c.mc_gain_sweep = struct('gain_levels', [0.5 1], 'bias_levels', [0 0.2], 'seeds', [3001 3002], ...
    'variants', variants);
pl = mc_gain_sweep_stage_plan(c);
assert(isequal({pl.runs.label}, {'main', 'EI3to2', 'EI2to3'}));
assert(isequal({pl.runs.folder}, {'', 'mc_gain_sweep_EI3to2', 'mc_gain_sweep_EI2to3'}));
a = pl.runs(1).args;
assert(strcmp(arg(a, 'preset_name'), c.mc_preset) && strcmp(arg(a, 'run_mode'), 'essential') && ...
    strcmp(arg(a, 'ode_solver'), c.mc_ode_solver) && isequal(arg(a, 'gain_levels'), [0.5 1]) && ...
    isequal(arg(a, 'bias_levels'), [0 0.2]) && isequal(arg(a, 'seeds'), [3001 3002]));
for k = 2:3
    a = pl.runs(k).args;
    assert(strcmp(arg(a, 'preset_name'), variants(k - 1).preset) && ...
        isequal(arg(a, 'gain_levels'), [0.5 1]) && isequal(arg(a, 'bias_levels'), 0.2) && ...
        isequal(arg(a, 'seeds'), [3001 3002]), ...
        '%s: gain only, at the preset''s own bias, with the main run''s gains and seeds', pl.runs(k).label);
end
c.mc_gain_sweep.variant_seeds = 3001;
c.mc_gain_sweep.model_args = {'S_c', 0.3};
pl = mc_gain_sweep_stage_plan(c);
a = pl.runs(3).args;
assert(isequal(arg(a, 'seeds'), 3001) && isequal(arg(a, 'bias_levels'), 0.3) && ...
    isequal(arg(pl.runs(1).args, 'seeds'), [3001 3002]), ...
    'variant_seeds and a model_args override of S_c reach the variant runs, not the main one');
% A second axis with no default range needs a range or its levels.
c.mc_gain_sweep = struct('second_axis', 'self_coupling');
assert(throws(@() mc_gain_sweep_stage_plan(c), 'run_mc_gain_sweep:NoSecondLevels'));

probe = fullfile(root, 'data', '__gainbias_preflight_probe__');
assert(~isfolder(probe));
refused = { ...
    @(b) setfield(b, 'mc_gain_sweep', setfield(b.mc_gain_sweep, 'gains', [1 2])), ...                  an unknown field
        'mc_gain_sweep_stage_plan:BadConfig'; ...
    @(b) setfield(b, 'mc_gain_sweep', setfield(b.mc_gain_sweep, 'gain_levels', [1 2])), ...            a grid's field beside the list
        'mc_gain_sweep_stage_plan:BadConfig'; ...
    @(b) in_grid(b, 1, 'spacing', 'log'), 'mc_gain_sweep_stage_plan:BadConfig'; ...                    an unknown grid field
    @(b) in_grid(b, 2, 'tag', ''), 'mc_gain_sweep_stage_plan:BadConfig'; ...                           a grid with no tag
    @(b) in_grid(b, 2, 'tag', 'bias'), 'mc_gain_sweep_stage_plan:BadConfig'; ...                       two grids, one folder
    @(b) in_grid(b, 2, 'tag', 'self coupling'), 'mc_gain_sweep_folder:BadTag'; ...                     no folder name
    @(b) in_grid(b, 2, 'bias_levels', [0 1]), 'mc_gain_sweep_stage_plan:BadConfig'; ...                bias levels on another axis
    @(b) in_variant(b, 1, 'preset', 'no_such_preset'), ''; ...                                         an unknown variant preset
    @(b) in_variant(b, 1, 'tag', 'E:I 3/2'), 'mc_gain_sweep_folder:BadTag'; ...
    @(b) in_variant(b, 2, 'tag', 'EI3to2'), 'mc_gain_sweep_stage_plan:BadConfig'};                   % two variants, one folder
for k = 1:size(refused, 1)
    bad = mu7chain_gainbias_config('essential'); bad.run_dir = probe;
    bad = refused{k, 1}(bad);
    assert(throws(@() run_all_paper_analyses(bad), refused{k, 2}), 'bad setting %d was not refused', k);
end
assert(~isfolder(probe), 'a refused run must not create its run directory');

%% 6. Blocks A and C as one bundle
for md = {'essential', 'fast'}
    for fs = {'3cond', '4cond'}
        ac = mu7chain_AC_config(md{1}, fs{1});
        A = mu7chain_config(md{1}, fs{1});
        C = mu7chain_gainbias_config(md{1}, fs{1});
        pl = paper_stage_plan(ac);
        if strcmp(md{1}, 'essential'); stages_A = ess_stages; else; stages_A = all10; end
        assert(isequal(pl.stages, [stages_A, {'mc_gain_sweep'}]), ...
            '%s: block A''s stages, then block C''s', md{1});
        pl_A = paper_stage_plan(A);
        assert(isequal(pl.sweeps, pl_A.sweeps) && isequal(pl.stages_skipped, pl_A.stages_skipped));
        assert(isequal(mc_gain_sweep_stage_plan(ac), mc_gain_sweep_stage_plan(C)), ...
            'the stage plans block C''s runs');
        assert(strcmp(ac.run_dir, ['data/mu7chain_AC_' md{1}]) && ...
            strcmp(ac.fig_root, ['figs/mu7chain_AC_' md{1} '_' fs{1}]));
        names = cellfun(@(f) f.name, ac.figures, 'UniformOutput', false);
        assert(isequal(names, [cellfun(@(f) f.name, A.figures, 'UniformOutput', false), entry_names]) && ...
            numel(unique(names)) == numel(names), 'block A''s entries, then block C''s');
        grouped = ac.figures(startsWith(names, 'fig_main'));
        assert(~isempty(grouped) && all(cellfun(@(f) strcmp(arg(f.args, 'source_fig_root'), ac.fig_root), grouped)), ...
            'the grouped figures read their sources from this bundle''s figure root');
        % mc_lean: the bundle sets it so block A's two memory-capacity
        % stages take the lean path (203d4f6; the same R^2 and memory
        % capacity to the bit, test_mc_lean_path, without holding the state
        % trajectory). mu7chain_config does not set it.
        assert(isequal(ac.mc_lean, true), 'the bundle runs memory capacity on the lean path');
        differs = {'stages', 'mc_gain_sweep', 'run_dir', 'fig_root', 'figures', 'mc_lean'};
        assert(isequal(rmfield(ac, differs), rmfield(A, intersect(differs, fieldnames(A)))), ...
            'the bundle differs from mu7chain_config in more than its stages, directories, registry and the lean memory-capacity path');
    end
end
assert(strcmp(mu7chain_AC_config().run_mode, 'essential'), 'the default is the size agents run');
assert(throws(@() mu7chain_AC_config('fast2'), 'mu7chain_config:BadRunMode'));

%% 7. End to end, tiny
tmp = tempname; tmpf = tempname; tmp_old = tempname;
cleanup = onCleanup(@() cellfun(@rmdir_if, {tmp, tmpf, tmp_old}));
gains = [0.8 1.6];
tiny = struct('seeds', 3001, 'conditions', {conds3}, 'model_args', {{'n', 50, 'indegree', 25}}, ...
    'protocol', struct('T_wash_sec', 2, 'T_train_sec', 12, 'T_test_sec', 6, 'd_max_sec', 1.5), ...
    'lya_T_sec', 6, 'lya_acc_sec', 3);
cfg = mu7chain_gainbias_config('essential', '3cond');
cfg.run_dir  = tmp;
cfg.fig_root = tmpf;
cfg.verbose  = 'near-none';
grids = cfg.mc_gain_sweep.grids;
[grids.gain_levels] = deal(gains);
grids(1).second_levels = [0.1 0.3];
grids(2).second_levels = [-1 1];
cfg.mc_gain_sweep = setfield(tiny, 'grids', grids); %#ok<SFLD>
t0 = tic;
run_dir = run_all_paper_analyses(cfg);
t_run = toc(t0);
assert(strcmp(run_dir, tmp));
M = load(fullfile(tmp, 'run_manifest.mat')); M = M.run_manifest;
assert(isequal(M.stages_run, {'mc_gain_sweep'}) && isequal(M.stages_skipped, all10));
results_in = @(folder) fullfile(tmp, folder, 'mc_gain_sweep', 'mc_gain_sweep_results.mat');
L = load(results_in('mc_gain_sweep_bias'));
assert(height(L.T) == 2 * 2 * 3 && ~any(L.T.failed), 'bias grid: 2 x 2 x 1 seed x 3 conditions, none failed');
assert(strcmp(L.settings.preset_name, cfg.mc_preset) && strcmp(L.settings.ode_solver, 'sra1'));
assert(isequal(L.settings.conditions, conds3) && strcmp(L.settings.second_axis, 'S_c'));
S = load(results_in('mc_gain_sweep_selfcoupling'));
assert(height(S.T) == 2 * 2 * 3 && ~any(S.T.failed), 'self-coupling grid: 2 x 2 x 1 seed x 3 conditions');
assert(strcmp(S.settings.second_axis, 'self_coupling') && isequal(S.settings.self_coupling_levels, [-1 1]) && ...
    isequal(unique(S.T.self_coupling)', [-1 1]) && S.settings.preset_self_coupling == 0);
for e = variant_tags
    E = load(results_in(['mc_gain_sweep_bias_' e{1}]));
    assert(endsWith(E.settings.preset_name, ['_' e{1}]));
    assert(isequal(E.settings.bias_levels, 0.2) && isequal(E.settings.gain_levels, gains) && ...
        isequal(E.settings.seeds, 3001), '%s: gain only, at the preset''s bias', e{1});
    assert(height(E.T) == 2 * 3 && ~any(E.T.failed));
end
assert(~isfolder(fullfile(tmp, 'mc_gain_sweep')), 'no unnamed grid: nothing is written to <run_dir>/mc_gain_sweep');
assert(~isfolder(fullfile(tmp, 'mc_example')) && ~isfolder(fullfile(tmp, 'eig_heatmap')));

t0 = tic;
R = make_all_paper_figures(cfg);
t_fig = toc(t0);
assert(isequal({R.name}, entry_names) && all([R.ok]), 'figures failed: %s', strjoin({R(~[R.ok]).name}, ', '));
for k = 1:numel(R)
    tag = file_tag(cfg.figures{k}.args);
    [~, saved] = cellfun(@fileparts, R(k).files, 'UniformOutput', false);
    assert(all(startsWith(saved, tag)), '%s: its files are named by its tag %s', R(k).name, tag);
end
stat = R(strcmp({R.name}, 'mc_statistics_selfcoupling'));
assert(any(endsWith(stat.files, 'MC_statistics_selfcoupling_by_self_coupling.csv')) && ...
    all(cellfun(@isfile, stat.files)));

% The resolver reads the grid it is asked for, and nothing else.
quiet = {'conditions', conds3, 'save', false, 'visible', false};
D = resolve_mc_gain_sweep('run_dir', tmp, 'grid', 'selfcoupling');
assert(strcmp(D.second_axis.property, 'self_coupling'));
assert(throws(@() resolve_mc_gain_sweep('run_dir', tmp), 'resolve_data_file:NotFound'), ...
    'no unnamed grid in this run: asking for one must not find another grid''s result');
assert(throws(@() resolve_mc_gain_sweep('run_dir', tmp, 'grid', 'feedback'), 'resolve_data_file:NotFound'));
assert(throws(@() fig_mc_gain_curves('run_dir', tmp, 'grid', 'bias', 'variants', {'EI1to4'}, quiet{:}), ...
    'resolve_data_file:NotFound'), 'a variant the run directory does not hold must be an error');
assert(throws(@() fig_mc_gain_curves('run_dir', tmp, 'grid', 'selfcoupling', 'variants', {'EI3to2'}, quiet{:}), ...
    'resolve_data_file:NotFound'), 'the variants belong to the bias grid');

% Gain curves: a row per network, two columns unless the third is asked for.
o = fig_mc_gain_curves('run_dir', tmp, 'grid', 'bias', 'variants', variant_tags, quiet{:});
ax = findobj(o.figs, 'Type', 'axes');
row_titles = arrayfun(@(a) char(a.Title.String), ax, 'UniformOutput', false);
n_rows3 = numel(ax); close(o.figs);
for want = {'E : I = 1 : 1', 'E : I = 3 : 2', 'E : I = 2 : 3'}
    assert(nnz(startsWith(row_titles, want{1})) == 1, 'no row labelled %s', want{1});
end
o = fig_mc_gain_curves('run_dir', tmp, 'grid', 'bias', quiet{:});
n_row1 = numel(findobj(o.figs, 'Type', 'axes')); close(o.figs);
o = fig_mc_gain_curves('run_dir', tmp, 'grid', 'bias', 'mc_vs_lambda1', true, quiet{:});
n_three = numel(findobj(o.figs, 'Type', 'axes')); close(o.figs);
assert(n_rows3 == 6 && n_row1 == 2 && n_three == 3, ...
    'gain curves: %d axes with the variants'' rows, %d without, %d with the third column', n_rows3, n_row1, n_three);

%% 8. A run directory with one unnamed grid still reads and resumes
mkdir(tmp_old);
copyfile(fullfile(tmp, 'mc_gain_sweep_bias', 'mc_gain_sweep'), fullfile(tmp_old, 'mc_gain_sweep'));
for e = variant_tags
    copyfile(fullfile(tmp, ['mc_gain_sweep_bias_' e{1}]), fullfile(tmp_old, ['mc_gain_sweep_' e{1}]));
end
old = mu7chain_gainbias_config('essential', '3cond');
old.verbose = 'near-none';
old.mc_gain_sweep = setfield(tiny, 'variants', variants); %#ok<SFLD>
old.mc_gain_sweep.gain_levels = gains;
old.mc_gain_sweep.bias_levels = [0.1 0.3];
D_old = resolve_mc_gain_sweep('run_dir', tmp_old);
assert(isequaln(D_old.T, read_mc_gain_sweep(results_in('mc_gain_sweep_bias')).T));
D_old = resolve_mc_gain_sweep('run_dir', tmp_old, 'variant', 'EI3to2');
assert(endsWith(D_old.preset_name, '_EI3to2'));
o = fig_mc_gain_curves('run_dir', tmp_old, 'variants', variant_tags, quiet{:});
assert(numel(findobj(o.figs, 'Type', 'axes')) == 6); close(o.figs);
run_mc_gain_sweep_stage(old, tmp_old, mc_gain_sweep_stage_plan(old));
for folder = [{''}, strcat('mc_gain_sweep_', variant_tags)]
    Z = load(fullfile(tmp_old, folder{1}, 'mc_gain_sweep', 'mc_gain_sweep_results.mat'));
    assert(Z.settings.n_run == 0 && ~any(Z.T.failed), ...
        'the earlier layout must resume with nothing to run (%s ran %d)', folder{1}, Z.settings.n_run);
end
assert(isequaln(Z.T, E.T), 'the resumed table is the one the job files held');

fprintf('PASS: test_mc_gain_bias_bundle -- plans unchanged, opt-in stage, both grids in every mode, blocks A and C, end to end (run %.0f s, figures %.0f s), the earlier layout.\n', ...
    t_run, t_fig);

%% ------------------------------------------------------------------------
function p = plan_literal(stages, skipped, sweeps, sweeps_skipped, is_sel, all_stages)
p = struct();
p.stages = stages; p.stages_skipped = skipped;
p.sweeps = sweeps; p.sweeps_skipped = sweeps_skipped;
p.is_selection = is_sel;
p.runs = struct();
for k = 1:numel(all_stages); p.runs.(all_stages{k}) = ismember(all_stages{k}, stages); end
end

function v = arg(args, name)
% The value a name-value list gives name (the last one, as a call takes it).
k = find(strcmp(args(1:2:end), name), 1, 'last');
assert(~isempty(k), 'no argument %s', name);
v = args{2 * k};
end

function tag = file_tag(args)
% The file tag a registry entry saves under: a figure's fig_tag, or the
% statistics tables' tag.
names = args(1:2:end);
tag = arg(args, names{find(ismember(names, {'fig_tag', 'tag'}), 1)});
end

function b = in_grid(b, g, field, value)
% The bundle with one field of its g-th grid set.
b.mc_gain_sweep.grids(g).(field) = value;
end

function b = in_variant(b, k, field, value)
% The bundle with one field of the bias grid's k-th variant set.
b.mc_gain_sweep.grids(1).variants(k).(field) = value;
end

function tf = throws(fn, id)
% fn throws; with id non-empty, that identifier.
try
    fn();
    tf = false;
catch err
    tf = isempty(id) || strcmp(err.identifier, id);
    if ~tf; fprintf(2, 'expected %s, got %s: %s\n', id, err.identifier, err.message); end
end
end

function rmdir_if(d)
if isfolder(d); rmdir(d, 's'); end
end
