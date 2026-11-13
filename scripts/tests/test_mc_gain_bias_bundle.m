%% test_mc_gain_bias_bundle.m
% Block C from the two entry points (FD-010): the opt-in stage
% 'mc_gain_sweep' and the bundle mu7chain_gainbias_config.
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
%      skipped, no sweep), its directories, registry, conditions per set,
%      variant networks, and that it is mu7chain_config's otherwise; the
%      stage's plan (the main run, then one gain-only run per variant at
%      its own bias, in its own folder); a bad mode, an unknown
%      cfg.mc_gain_sweep field, an unknown variant preset and a tag that
%      is no folder name are refused before the run directory exists.
%   4. END TO END on a tiny case: run_all_paper_analyses on the bundle
%      ('essential', 3cond) with a 50-neuron network, the shortened protocol
%      of test_mc_gain_sweep, a 2 x 2 grid, one seed and three conditions,
%      serial; then make_all_paper_figures from that run directory. The main
%      result and the two variants' gain-only results (one bias level each,
%      the preset's 0.2) are where the figures look; all four registry
%      entries write files; the gain-curves figure has one row per network
%      named, labelled by its E : I ratio, and a variant the run directory
%      does not hold is an error.
%
% Asserts; prints one PASS line. About two minutes, serial, no pool needed.
% Assumes setup_paths has run and a pool allocation is set (pool_set(0)
% will do): both entry points refuse without one.

%% 1. Existing bundles' plans are unchanged
all10 = {'sweeps', 'memory_capacity', 'mc_example', 'eig_heatmap', 'numerics', ...
    'lyapunov_spectrum', 'transient_gain', 'local_lyapunov', 'dc_lle', 'paper_illustrations'};
sw3 = {'sensitivity', 'tau', 'param_space'};
full = plan_literal(all10, cell(1, 0), sw3, cell(1, 0), false, all10);
ess  = plan_literal({'sweeps', 'eig_heatmap', 'local_lyapunov', 'paper_illustrations'}, ...
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
        assert(isequal(names, {'fig_mc_gain_bias_maps', 'fig_mc_gain_curves', ...
            'fig_mc_robustness', 'mc_gain_bias_statistics'}));
        if strcmp(fs{1}, '3cond'); want = conds3; else; want = conds4; end
        assert(isequal(c.conditions, want));
        assert(isequal(c.figure_sets, base.figure_sets));
        assert(isequal({c.mc_gain_sweep.variants.tag}, {'EI3to2', 'EI2to3'}));
        assert(isequal({c.mc_gain_sweep.variants.preset}, ...
            {[base.mc_preset '_EI3to2'], [base.mc_preset '_EI2to3']}));
        curves = c.figures{strcmp(names, 'fig_mc_gain_curves')};
        assert(isequal(curves.args, {'variants', {'EI3to2', 'EI2to3'}}), ...
            'the gain-curves entry must name the variants the stage runs');
        differs = {'stages', 'sweeps', 'eig_heatmap_examples', 'mc_gain_sweep', 'run_dir', ...
            'fig_root', 'figures'};
        assert(isequal(rmfield(c, differs), rmfield(base, intersect(differs, fieldnames(base)))), ...
            'the bundle differs from mu7chain_config in more than its stage, directories and registry');
    end
end
assert(throws(@() mu7chain_gainbias_config('fast2'), 'mu7chain_gainbias_config:BadRunMode'));
assert(strcmp(mu7chain_gainbias_config().run_mode, 'essential'), 'the default is the size agents run');

% The stage's plan: the main run into the run directory, then each variant
% along gain only, at its own bias, into its own folder, with the main
% run's options (and seeds, unless variant_seeds says otherwise).
c = mu7chain_gainbias_config('essential');
c.mc_gain_sweep.gain_levels = [0.5 1];
c.mc_gain_sweep.bias_levels = [0 0.2];
c.mc_gain_sweep.seeds = [3001 3002];
pl = mc_gain_sweep_stage_plan(c);
assert(isequal({pl.runs.label}, {'main', 'EI3to2', 'EI2to3'}));
assert(isequal({pl.runs.folder}, {'', 'mc_gain_sweep_EI3to2', 'mc_gain_sweep_EI2to3'}));
a = pl.runs(1).args;
assert(strcmp(arg(a, 'preset_name'), c.mc_preset) && strcmp(arg(a, 'run_mode'), 'essential') && ...
    strcmp(arg(a, 'ode_solver'), c.mc_ode_solver) && isequal(arg(a, 'gain_levels'), [0.5 1]) && ...
    isequal(arg(a, 'bias_levels'), [0 0.2]) && isequal(arg(a, 'seeds'), [3001 3002]));
for k = 2:3
    a = pl.runs(k).args;
    assert(strcmp(arg(a, 'preset_name'), c.mc_gain_sweep.variants(k - 1).preset) && ...
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

probe = fullfile(root, 'data', '__gainbias_preflight_probe__');
assert(~isfolder(probe));
bad = mu7chain_gainbias_config('essential'); bad.run_dir = probe;
bad.mc_gain_sweep.gains = [1 2];
assert(throws(@() run_all_paper_analyses(bad), 'mc_gain_sweep_stage_plan:BadConfig'));
bad = mu7chain_gainbias_config('essential'); bad.run_dir = probe;
bad.mc_gain_sweep.variants(1).preset = 'no_such_preset';
assert(throws(@() run_all_paper_analyses(bad), ''));
bad = mu7chain_gainbias_config('essential'); bad.run_dir = probe;
bad.mc_gain_sweep.variants(1).tag = 'E:I 3/2';
assert(throws(@() run_all_paper_analyses(bad), 'mc_gain_sweep_variant_folder:BadTag'));
bad = mu7chain_gainbias_config('essential'); bad.run_dir = probe;
bad.mc_gain_sweep.variants(2).tag = bad.mc_gain_sweep.variants(1).tag;
assert(throws(@() run_all_paper_analyses(bad), 'mc_gain_sweep_stage_plan:BadConfig'));
assert(~isfolder(probe), 'a refused run must not create its run directory');

%% 4. End to end, tiny
tmp = tempname; tmpf = tempname;
cleanup = onCleanup(@() cellfun(@rmdir_if, {tmp, tmpf}));
cfg = mu7chain_gainbias_config('essential', '3cond');
cfg.run_dir  = tmp;
cfg.fig_root = tmpf;
cfg.verbose  = 'near-none';
cfg.mc_gain_sweep.gain_levels = [0.8 1.6];
cfg.mc_gain_sweep.bias_levels = [0.1 0.3];
cfg.mc_gain_sweep.seeds       = 3001;
cfg.mc_gain_sweep.conditions  = conds3;
cfg.mc_gain_sweep.model_args  = {'n', 50, 'indegree', 25};
cfg.mc_gain_sweep.protocol    = struct('T_wash_sec', 2, 'T_train_sec', 12, 'T_test_sec', 6, 'd_max_sec', 1.5);
cfg.mc_gain_sweep.lya_T_sec   = 6;
cfg.mc_gain_sweep.lya_acc_sec = 3;
t0 = tic;
run_dir = run_all_paper_analyses(cfg);
t_run = toc(t0);
assert(strcmp(run_dir, tmp));
M = load(fullfile(tmp, 'run_manifest.mat')); M = M.run_manifest;
assert(isequal(M.stages_run, {'mc_gain_sweep'}) && isequal(M.stages_skipped, all10));
L = load(fullfile(tmp, 'mc_gain_sweep', 'mc_gain_sweep_results.mat'));
assert(height(L.T) == 2 * 2 * 3 && ~any(L.T.failed), 'main: 2 x 2 x 1 seed x 3 conditions, none failed');
assert(strcmp(L.settings.preset_name, cfg.mc_preset) && strcmp(L.settings.ode_solver, 'sra1'));
assert(isequal(L.settings.conditions, conds3));
for e = {'EI3to2', 'EI2to3'}
    f = fullfile(tmp, ['mc_gain_sweep_' e{1}], 'mc_gain_sweep', 'mc_gain_sweep_results.mat');
    E = load(f);
    assert(endsWith(E.settings.preset_name, ['_' e{1}]));
    assert(isequal(E.settings.bias_levels, 0.2) && isequal(E.settings.gain_levels, [0.8 1.6]) && ...
        isequal(E.settings.seeds, 3001), '%s: gain only, at the preset''s bias', e{1});
    assert(height(E.T) == 2 * 3 && ~any(E.T.failed));
end
assert(~isfolder(fullfile(tmp, 'mc_example')) && ~isfolder(fullfile(tmp, 'eig_heatmap')));

t0 = tic;
R = make_all_paper_figures(cfg);
t_fig = toc(t0);
assert(all([R.ok]) && numel(R) == 4, 'figures failed: %s', strjoin({R(~[R.ok]).name}, ', '));
stat = R(strcmp({R.name}, 'mc_gain_bias_statistics'));
assert(any(endsWith(stat.files, 'MC_gain_bias_statistics.csv')) && all(cellfun(@isfile, stat.files)));
quiet = {'conditions', conds3, 'save', false, 'visible', false};
o = fig_mc_gain_curves('run_dir', tmp, 'variants', {'EI3to2', 'EI2to3'}, quiet{:});
ax = findobj(o.figs, 'Type', 'axes');
row_titles = arrayfun(@(a) char(a.Title.String), ax, 'UniformOutput', false);
n3 = numel(ax); close(o.figs);
for want = {'E : I = 1 : 1', 'E : I = 3 : 2', 'E : I = 2 : 3'}
    assert(nnz(startsWith(row_titles, want{1})) == 1, 'no row labelled %s', want{1});
end
o = fig_mc_gain_curves('run_dir', tmp, quiet{:});
n1 = numel(findobj(o.figs, 'Type', 'axes')); close(o.figs);
assert(n3 == 9 && n1 == 3, 'gain curves: %d axes with the variants'' rows, %d without', n3, n1);
assert(throws(@() fig_mc_gain_curves('run_dir', tmp, 'variants', {'EI1to4'}, quiet{:}), ...
    'resolve_data_file:NotFound'), 'a variant the run directory does not hold must be an error');

fprintf('PASS: test_mc_gain_bias_bundle -- plans unchanged, opt-in stage, bundle, end to end (run %.0f s, figures %.0f s).\n', ...
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
