%% test_mc_gain_bias_bundle.m
% Block C from the two entry points (FD-010, FD-024, FD-027): the opt-in
% stage 'mc_gain_sweep', the bundle mu7chain_gainbias_config with its five
% grids (bias, self-coupling, and the noise grids releasenoise,
% selfcoupling_releasenoise and wienernoise) and the two noise-level time
% series, and the bundle of blocks A and C, mu7chain_AC_config.
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
%      The noise grids: a zero then the mode's log ladder on the noise axis,
%      their fixed model_args and the filtered readout, the self-coupling
%      noise grid with no gain axis; the first four runs plan exactly the
%      arguments they had. The 'zero_log' spacing of the grid table on
%      either axis, and how the maps lay it out (mc_axis_spacing).
%   5. The form with one unnamed grid still plans as it did (the main run
%      into the run directory, variants into mc_gain_sweep_<tag>); a grid's
%      model_args follow the stage's and its readout replaces the stage's;
%      bad settings are refused before the run directory exists.
%   6. mu7chain_AC_config: block A's stages then block C's, then the
%      opt-in 'noise_switch' (which neither part runs); both registries and
%      the noise-switch entry; its own directories; the grouped figures
%      pointed at its figure root. (The noise-switch stage's size table and
%      its end-to-end run are test_noise_switch_stage.)
%   7. END TO END on a tiny case: run_all_paper_analyses on the bundle
%      ('essential', 3cond) with a 50-neuron network, the shortened protocol
%      of test_mc_gain_sweep, 2 x 2 levels per grid, one seed and three
%      conditions; then make_all_paper_figures from that run directory.
%      Every grid's results and the two variants' are where the figures
%      look (the noise grids at 2 or 3 noise levels, their fixed arguments
%      and readout in the settings), the time series simulate a 50-neuron
%      network for 2 s, a noise axis is drawn as its log ladder with the
%      preset's star on the zero row, the statistics take widths along it
%      in its own units, every registry entry writes files under its own tag, the gain curves
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
assert(isequal(st, all10) && isequal(opt_in, {'mc_gain_sweep', 'noise_switch'}));
pl = paper_stage_plan(struct('illustrations', true, 'stages', {{'mc_gain_sweep', 'dc_lle'}}));
assert(isequal(pl.stages, {'dc_lle', 'mc_gain_sweep'}), 'opt-in stages run after the first list');
assert(isequal(pl.stages_skipped, setdiff(all10, {'dc_lle'}, 'stable')) && pl.runs.mc_gain_sweep);
assert(throws(@() paper_stage_plan(struct('stages', {{'mc_gain_swep'}})), 'paper_stage_plan:UnknownStage'));

%% 3. The bundle
root = fileparts(which('setup_paths'));
conds3 = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
conds4 = [conds3, {'sfa3_std2'}];
variant_tags = {'EI3to2', 'EI2to3'};
noise_grids = {'releasenoise', 'selfcoupling_releasenoise', 'wienernoise'};
entry_names = [{'fig_mc_maps_bias', 'fig_mc_maps_selfcoupling', 'fig_mc_gain_curves', ...
    'fig_mc_robustness_bias', 'fig_mc_robustness_selfcoupling', 'fig_mc_vs_lambda1', ...
    'fig_lambda1_vs_bias', 'fig_rate_vs_lambda1', 'fig_slope_share_vs_gain', ...
    'mc_statistics_bias', 'mc_statistics_selfcoupling'}, ...
    strcat('fig_mc_maps_', noise_grids), strcat('fig_mc_robustness_', noise_grids), ...
    strcat('mc_statistics_', noise_grids), {'fig_noise_levels_wiener', 'fig_noise_levels_gamma'}];
entry_grids = [{'bias', 'selfcoupling', 'bias', 'bias', 'selfcoupling', 'selfcoupling', ...
    'bias', 'bias', 'bias', 'bias', 'selfcoupling'}, noise_grids, noise_grids, noise_grids, {'', ''}];
% The noise grids (FD-027 P3): the leak, the noise ladders and what each holds fixed.
LEAK = 0.002;
gamma_ladder = [1e-4 1e-2]; wiener_ladder = [2.5e-4 0.025];
release_fixed = {'sigma_u_noise', 0, 'release_rate_all', LEAK};
wiener_fixed  = {'release_rate_all', 0, 'release_sigma_all', 0};
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
        assert(isequal(cellfun(@(f) arg_or_empty(f.args, 'grid'), c.figures, 'UniformOutput', false), entry_grids), ...
            'every entry names the grid it reads (the time series read none)');
        if strcmp(fs{1}, '3cond'); want = conds3; else; want = conds4; end
        assert(isequal(c.conditions, want));
        assert(isequal(c.figure_sets, base.figure_sets));
        grids = c.mc_gain_sweep.grids;
        assert(isequal({grids.tag}, [{'bias', 'selfcoupling'}, noise_grids]) && ...
            isequal({grids.second_axis}, {'S_c', 'self_coupling', 'release_sigma_all', ...
                'release_sigma_all', 'sigma_u_noise'}) && ...
            isequal({grids.gain_spacing}, {'log', 'linear', 'log', [], 'log'}) && ...
            isequal({grids.first_axis}, {[], [], [], 'self_coupling', []}));
        assert(isequal({grids(1).variants.tag}, variant_tags) && all(cellfun(@isempty, {grids(2:end).variants})));
        assert(all(cellfun(@isempty, {grids(1:2).model_args})) && all(cellfun(@isempty, {grids(1:2).readout_signal})), ...
            'the bias and self-coupling grids keep the stage''s arguments and readout');
        assert(isequal({grids(3:5).model_args}, {release_fixed, release_fixed, wiener_fixed}) && ...
            all(strcmp({grids(3:5).readout_signal}, 'release_filtered')) && ...
            all(strcmp({grids(3:5).second_spacing}, 'zero_log')), ...
            'the noise grids: their fixed arguments, the filtered readout, a zero then log-spaced noise axis');
        % The time-series entries: zero, the middle and the top of each noise axis.
        ts = c.figures(startsWith(names, 'fig_noise_levels_'));
        assert(all(cellfun(@(f) isequal(f.fn, @fig_noise_levels_timeseries) && ~f.in_paper && ...
            strcmp(arg(f.args, 'preset_name'), c.preset_name), ts)));
        assert(strcmp(arg(ts{1}.args, 'level_name'), 'sigma_u_noise') && ...
            isequal(arg(ts{1}.args, 'levels'), [0 0.0025 0.025]) && ...
            isequal(arg(ts{1}.args, 'model_args'), wiener_fixed), 'the Wiener time series');
        assert(strcmp(arg(ts{2}.args, 'level_name'), 'release_sigma_all') && ...
            isequal(arg(ts{2}.args, 'levels'), [0 1e-3 1e-2]) && ...
            isequal(arg(ts{2}.args, 'model_args'), release_fixed), 'the gamma time series');
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
    assert(isequal({pl.runs.label}, [{'bias', 'bias_EI3to2', 'bias_EI2to3', 'selfcoupling'}, noise_grids]), md);
    assert(isequal({pl.runs.folder}, [{'mc_gain_sweep_bias', 'mc_gain_sweep_bias_EI3to2', ...
        'mc_gain_sweep_bias_EI2to3', 'mc_gain_sweep_selfcoupling'}, strcat('mc_gain_sweep_', noise_grids)]), md);
    for k = 1:4
        a = pl.runs(k).args;
        assert(~any(strcmp(a(1:2:end), 'model_args')) && ~any(strcmp(a(1:2:end), 'readout_signal')), ...
            '%s, %s: the first four runs plan the arguments they always had', md, pl.runs(k).label);
    end
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
    % The noise grids: a zero level, then the mode's log ladder without its
    % lowest rung; fixed arguments and the readout per grid.
    zero_log = @(ladder) round([0, ladder(1) * (ladder(2) / ladder(1)).^((1:n-1) / (n - 1))], 10);
    sigma_gamma = arg(pl.runs(5).args, 'second_levels');
    assert(max(abs(sigma_gamma - zero_log(gamma_ladder))) < 1e-12 && sigma_gamma(end) == gamma_ladder(2), ...
        '%s: the release-noise levels', md);
    sigma_u = arg(pl.runs(7).args, 'second_levels');
    assert(max(abs(sigma_u - zero_log(wiener_ladder))) < 1e-12 && sigma_u(end) == wiener_ladder(2), ...
        '%s: the Wiener-noise levels', md);
    assert(strcmp(mc_axis_spacing(sigma_gamma), 'zero_log') && strcmp(mc_axis_spacing(sigma_u), 'zero_log'));
    a = pl.runs(5).args;
    assert(strcmp(arg(a, 'second_axis'), 'release_sigma_all') && isequal(arg(a, 'gain_levels'), log_gains) && ...
        isequal(arg(a, 'model_args'), release_fixed) && strcmp(arg(a, 'readout_signal'), 'release_filtered'), ...
        '%s: gain x release noise', md);
    a = pl.runs(6).args;
    assert(strcmp(arg(a, 'first_axis'), 'self_coupling') && isequal(arg(a, 'first_levels'), round(linspace(-4, 4, n), 10)) && ...
        isequal(arg(a, 'second_levels'), sigma_gamma) && isequal(arg(a, 'model_args'), release_fixed) && ...
        strcmp(arg(a, 'readout_signal'), 'release_filtered') && isequal(arg(a, 'seeds'), seeds), ...
        '%s: self-coupling x release noise, at the preset''s gain', md);
    assert(~any(strcmp(a(1:2:end), 'gain_levels')), '%s: no gain axis on the self-coupling noise grid', md);
    a = pl.runs(7).args;
    assert(strcmp(arg(a, 'second_axis'), 'sigma_u_noise') && isequal(arg(a, 'gain_levels'), log_gains) && ...
        isequal(arg(a, 'model_args'), wiener_fixed) && strcmp(arg(a, 'readout_signal'), 'release_filtered'), ...
        '%s: gain x Wiener noise', md);
    assert(isequal([pl.runs.n_jobs], [n^2, n, n, n^2, n^2, n^2, n^2] * n_seeds) && ...
        pl.n_jobs == (5 * n^2 + 2 * n) * n_seeds && abs(pl.worker_minutes - 9.5 * pl.n_jobs) < 1e-9, ...
        '%s: job counts', md);
end
% The essential stage holds 51 jobs: 24 before the noise grids, 27 of them.
assert(mc_gain_sweep_stage_plan(mu7chain_gainbias_config('essential')).n_jobs == 51);

% zero_log is a spacing of the grid table, on either axis, for any property.
G = mc_gain_sweep_grid('fast', 'second_axis', 'self_coupling', 'second_range', [0.01 1], ...
    'second_spacing', 'zero_log');
assert(isequal(G.second_levels, round([0, 10.^(-2 + 2 * (1:6) / 6)], 10)), ...
    'zero, then the 7-level ladder from 0.01 to 1 without its lowest rung, rounded to 10 decimals');
G = mc_gain_sweep_grid('essential', 'first_axis', 'sigma_u_noise', 'first_range', [1e-3 1e-1], ...
    'first_spacing', 'zero_log', 'second_levels', 0.2);
assert(max(abs(G.first_levels - [0 1e-2 1e-1])) < 1e-12, 'at three levels: zero, the middle, the top');
assert(throws(@() mc_gain_sweep_grid('fast', 'second_axis', 'sigma_u_noise', 'second_range', [0 0.1], ...
    'second_spacing', 'zero_log'), 'mc_gain_sweep_grid:BadRange'), 'a zero_log range names the ladder, lo > 0');
assert(throws(@() mc_gain_sweep_grid('fast', 'second_axis', 'sigma_u_noise', 'second_range', [1e-3 0.1], ...
    'second_spacing', 'cubic'), 'mc_gain_sweep_grid:BadSpacing'));
% How the maps lay such an axis out: evenly spaced cells, the zero one log
% step below the lowest positive level; other axes as before.
[sp, pos] = mc_axis_spacing([0 1e-3 1e-2 1e-1]);
assert(strcmp(sp, 'zero_log') && max(abs(pos([0 1e-3 1e-2 1e-1]) - (-4:-1))) < 1e-12 && pos(0) == -4);
assert(strcmp(mc_axis_spacing([0 0.5 1]), 'linear') && strcmp(mc_axis_spacing([-4 0 4]), 'linear') && ...
    strcmp(mc_axis_spacing([0.1 0.4 4]), 'log') && strcmp(mc_axis_spacing([0 1e-3 1e-2 2e-1]), 'linear'), ...
    'only a zero followed by a log ladder is zero_log');
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

% A grid's model_args come after the stage's; its readout replaces the stage's.
c = mu7chain_gainbias_config('essential');
c.mc_gain_sweep.model_args = {'n', 50, 'indegree', 25};
c.mc_gain_sweep.readout_signal = 'rate';
pl = mc_gain_sweep_stage_plan(c);
assert(isequal(arg(pl.runs(1).args, 'model_args'), {'n', 50, 'indegree', 25}) && strcmp(arg(pl.runs(1).args, 'readout_signal'), 'rate'));
assert(isequal(arg(pl.runs(2).args, 'model_args'), {'n', 50, 'indegree', 25}), 'a variant takes the stage''s model_args');
assert(isequal(arg(pl.runs(5).args, 'model_args'), [{'n', 50, 'indegree', 25}, release_fixed]) && ...
    strcmp(arg(pl.runs(5).args, 'readout_signal'), 'release_filtered'));

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
    @(b) in_grid(b, 3, 'model_args', {'n'}), 'mc_gain_sweep_stage_plan:BadConfig'; ...                 a grid's model_args not in pairs
    @(b) in_grid(b, 4, 'gain_range', [0.1 4]), 'mc_gain_sweep_stage_plan:BadConfig'; ...               gain options with self-coupling first
    @(b) in_grid(b, 5, 'second_spacing', 'cubic'), 'mc_gain_sweep_grid:BadSpacing'; ...                an unknown spacing
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
        assert(isequal(pl.stages, [stages_A, {'mc_gain_sweep', 'noise_switch'}]), ...
            '%s: block A''s stages, then block C''s, then the noise switch', md{1});
        assert(pl.runs.noise_switch && ~ismember('noise_switch', [paper_stage_plan(A).stages, ...
            paper_stage_plan(C).stages]), 'only the A + C bundle runs the noise switch');
        pl_A = paper_stage_plan(A);
        assert(isequal(pl.sweeps, pl_A.sweeps) && isequal(pl.stages_skipped, pl_A.stages_skipped));
        assert(isequal(mc_gain_sweep_stage_plan(ac), mc_gain_sweep_stage_plan(C)), ...
            'the stage plans block C''s runs');
        assert(strcmp(ac.run_dir, ['data/mu7chain_AC_' md{1}]) && ...
            strcmp(ac.fig_root, ['figs/mu7chain_AC_' md{1} '_' fs{1}]));
        names = cellfun(@(f) f.name, ac.figures, 'UniformOutput', false);
        assert(isequal(names, [cellfun(@(f) f.name, A.figures, 'UniformOutput', false), entry_names, ...
            {'fig_noise_switch_rollouts'}]) && numel(unique(names)) == numel(names), ...
            'block A''s entries, then block C''s, then the noise switch''s');
        ns = ac.figures{end};
        assert(isequal(ns.fn, @fig_noise_switch_rollouts) && isempty(ns.args) && ~ns.in_paper, ...
            'the noise-switch entry draws through fig_noise_switch_rollouts with the bundle''s arguments only');
        ns_plan = noise_switch_stage_plan(ac);
        assert(strcmp(ns_plan.folder, 'noise_switch') && strcmp(arg(ns_plan.args, 'preset_name'), ac.preset_name) && ...
            isequal(ns_plan.conditions, conds4), ...
            'the noise switch runs the bundle''s noisy preset, every condition, into noise_switch/');
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
[grids([1 2 3 5]).gain_levels] = deal(gains);
grids(4).first_levels = [-1 1];          % the self-coupling noise grid has no gain axis
grids(1).second_levels = [0.1 0.3];
grids(2).second_levels = [-1 1];
grids(3).second_levels = [0 2e-4 2e-3];   % three levels: a zero, then a log ladder
grids(4).second_levels = [0 2e-3];
grids(5).second_levels = [0 1e-3 0.01];
cfg.mc_gain_sweep = setfield(tiny, 'grids', grids); %#ok<SFLD>
% The time series simulate: on a tiny network, briefly.
for k = find(startsWith(cellfun(@(f) f.name, cfg.figures, 'UniformOutput', false), 'fig_noise_levels_'))
    a = cfg.figures{k}.args;
    cfg.figures{k}.args = [a, {'model_args', [arg(a, 'model_args'), {'n', 50, 'indegree', 25}], ...
        'T_settle', 1, 'T_show', 1, 'n_show_E', 3, 'n_show_I', 3}];
end
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
% The noise grids: each its own folder, fixed arguments and readout recorded.
N = struct();
for g = noise_grids
    N.(g{1}) = load(results_in(['mc_gain_sweep_' g{1}]));
    n_levels = numel(grids(strcmp({grids.tag}, g{1})).second_levels);
    assert(height(N.(g{1}).T) == 2 * n_levels * 3 && ~any(N.(g{1}).T.failed), '%s: 2 x %d x 3 conditions, none failed', g{1}, n_levels);
    assert(strcmp(N.(g{1}).settings.readout_signal, 'release_filtered'), '%s reads the filtered release', g{1});
end
assert(isequal(N.releasenoise.settings.model_args, [tiny.model_args, release_fixed]) && ...
    strcmp(N.releasenoise.settings.second_axis, 'release_sigma_all') && ...
    isequal(unique(N.releasenoise.T.release_sigma_all)', [0 2e-4 2e-3]) && N.releasenoise.settings.preset_release_sigma_all == 0);
assert(strcmp(N.selfcoupling_releasenoise.settings.first_axis, 'self_coupling') && ...
    isequal(unique(N.selfcoupling_releasenoise.T.self_coupling)', [-1 1]));
assert(isequal(N.wienernoise.settings.model_args, [tiny.model_args, wiener_fixed]) && ...
    isequal(unique(N.wienernoise.T.sigma_u_noise)', [0 1e-3 0.01]));
assert(~isequal(N.releasenoise.T.mc_total(N.releasenoise.T.release_sigma_all > 0), ...
    N.releasenoise.T.mc_total(N.releasenoise.T.release_sigma_all == 0)), 'the release noise changes the result');
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

% The noise axis of a map is drawn as its log ladder: evenly spaced cells,
% the zero one step below the lowest positive level, ticks labelled 0 and
% the levels, and the caption says the axis is broken.
quiet = {'conditions', conds3, 'save', false, 'visible', false};
o = fig_mc_gain_bias_maps('run_dir', tmp, 'grid', 'wienernoise', 'mc_scale', 'linear', 'lle_scale', 'linear', quiet{:});
ax = findobj(o.figs, 'Type', 'axes');
assert(all(arrayfun(@(a) max(abs(a.YTick - log10([1e-4 1e-3 1e-2]))) < 1e-12, ax)), ...
    'the Wiener axis {0, 0.001, 0.01}: evenly spaced, 0 one log step below 0.001');
labelled = ax(arrayfun(@(a) ~isempty(a.YTickLabel), ax));
assert(~isempty(labelled) && all(arrayfun(@(a) isequal(cellstr(a.YTickLabel)', {'0', '0.001', '0.01'}), labelled)), ...
    'the noise ticks are labelled 0 and the levels');
caption = ax(1).Parent.Parent.Title.String;
assert(any(contains(caption, 'broken')), 'the caption states the broken noise axis');
close(o.figs);
D_rn = resolve_mc_gain_sweep('run_dir', tmp, 'grid', 'releasenoise');
o = fig_mc_gain_bias_maps('data', D_rn, quiet{:});
star = findobj(o.figs, 'Marker', 'pentagram');
assert(~isempty(star) && all(abs([star.YData] - log10(2e-5)) < 1e-12), ...
    'the preset (no release noise) sits on the zero row, one log step below 2e-4');
close(o.figs);
[Sn, ~, SBn] = mc_gain_bias_statistics('data', D_rn, 'conditions', conds3, 'verbose', 'near-none');
assert(all(Sn.fit_release_sigma_all == 0) && ismember('release_sigma_all_range_best', Sn.Properties.VariableNames) && ...
    ismember('release_sigma_all', SBn.Properties.VariableNames), ...
    'statistics along the noise axis: in sigma''s units, the fit at the preset''s sigma, 0');

% The resolver reads the grid it is asked for, and nothing else.
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

fprintf('PASS: test_mc_gain_bias_bundle -- plans unchanged, opt-in stage, five grids in every mode, blocks A and C, end to end (run %.0f s, figures %.0f s), the earlier layout.\n', ...
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

function v = arg_or_empty(args, name)
% The value a name-value list gives name, or '' when it gives none.
v = '';
if any(strcmp(args(1:2:end), name)); v = arg(args, name); end
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
