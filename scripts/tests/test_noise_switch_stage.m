%% test_noise_switch_stage.m
% The noise-switch rollouts as an opt-in stage of run_all_paper_analyses
% ('noise_switch'; noise_switch_stage_plan, run_noise_switch_stage) and as
% a figure entry of the A + C bundle (Tom, 2026-11-14).
%
% Checks:
%   1. The plan in EVERY run mode, from mu7chain_AC_config: the size table
%      (copies and rollout per mode, 10 s shared), the bundle's noisy preset
%      and integrator mode, every condition the preset states, the folder.
%   2. cfg.noise_switch overrides the size and the conditions; a bad
%      setting is refused by the plan, and through run_all_paper_analyses
%      before the run directory exists (an unknown field, one copy, a
%      condition the preset lacks, a noise-free network).
%   3. END TO END on a tiny case: a bundle naming only the stage, the
%      mu7chain network at n = 60 (test_noise_switch_rollouts' overrides),
%      3 copies, 2 s shared, 1 s rollout, the 3cond set's conditions;
%      run_all_paper_analyses writes <run_dir>/noise_switch/, the manifest
%      records the stage, and make_all_paper_figures draws the registered
%      entry: one time-series figure per condition and the summary, none a
%      prefix of another.
%
% Asserts; prints one PASS line. Under a minute with a small pool. Assumes
% setup_paths has run and a pool allocation is set (pool_set(0) will do):
% both entry points refuse without one.
%
% See also: noise_switch_stage_plan, run_noise_switch_stage,
%           test_noise_switch_rollouts, test_mc_gain_bias_bundle

%% 1. The plan in every run mode
root = fileparts(which('setup_paths'));
sizes = {'essential', 4, 5; 'fast', 8, 10; 'medium', 16, 20; 'medium2', 16, 20; 'production', 32, 20};
assert(isequal(sizes(:, 1)', run_mode_names()), 'the size table below must cover run_mode_names');
for m = 1:size(sizes, 1)
    [md, n_copies, t_rollout] = sizes{m, :};
    ac = mu7chain_AC_config(md);
    [~, ~, conditions] = srnn_param_preset(ac.preset_name);
    all_conds = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
    pl = noise_switch_stage_plan(ac);
    assert(pl.n_copies == n_copies && pl.t_rollout == t_rollout && pl.t_shared == 10, ...
        '%s: %d copies x %g s rollout after 10 s shared', md, pl.n_copies, pl.t_rollout);
    assert(strcmp(pl.folder, 'noise_switch') && isequal(pl.conditions, all_conds));
    assert(strcmp(arg(pl.args, 'preset_name'), ac.preset_name) && strcmp(arg(pl.args, 'run_mode'), md) && ...
        arg(pl.args, 'n_copies') == n_copies && arg(pl.args, 't_rollout') == t_rollout && ...
        isequal(arg(pl.args, 'conditions'), all_conds) && isempty(arg(pl.args, 'model_args')));
end

%% 2. Overrides, and what is refused
c = mu7chain_AC_config('essential', '3cond');
c.noise_switch = struct('n_copies', 3, 't_rollout', 2.5, 't_shared', 4, 'conditions', {c.conditions});
pl = noise_switch_stage_plan(c);
assert(pl.n_copies == 3 && pl.t_rollout == 2.5 && pl.t_shared == 4 && isequal(pl.conditions, c.conditions));

probe = fullfile(root, 'data', '__noise_switch_preflight_probe__');
assert(~isfolder(probe));
refused = { ...
    struct('copies', 3), ...                                   an unknown field
    struct('n_copies', 1), ...                                 one copy has nothing to compare
    struct('conditions', {{'no_adaptation', 'sfa9_std9'}}), ... a condition the preset lacks
    struct('model_args', {{'sigma_u_noise', 0}})};            % a noise-free network
for k = 1:numel(refused)
    bad = mu7chain_AC_config('essential');
    bad.run_dir = probe;
    bad.noise_switch = refused{k};
    assert(throws(@() noise_switch_stage_plan(bad), 'noise_switch_stage_plan:BadConfig'), 'setting %d', k);
    assert(throws(@() run_all_paper_analyses(bad), 'noise_switch_stage_plan:BadConfig'), ...
        'setting %d was not refused by the entry point', k);
end
bad = mu7chain_AC_config('essential');
bad.preset_name = bad.mc_preset;                              % the noise-free memory-capacity twin
assert(throws(@() noise_switch_stage_plan(bad), 'noise_switch_stage_plan:BadConfig'));
assert(~isfolder(probe), 'a refused run must not create its run directory');

%% 3. End to end, tiny
tmp = tempname; tmpf = tempname;
cleanup = onCleanup(@() cellfun(@rmdir_if, {tmp, tmpf}));
cfg = mu7chain_AC_config('essential', '3cond');
cfg.stages = {'noise_switch'};
cfg.sweeps = {};
cfg.eig_heatmap_examples = [];
cfg.run_dir  = tmp;
cfg.fig_root = tmpf;
cfg.verbose  = 'near-none';
cfg.noise_switch = struct('n_copies', 3, 't_shared', 2, 't_rollout', 1, ...
    'conditions', {cfg.conditions}, 'model_args', {{'n', 60, 'indegree', 20}});
names = cellfun(@(f) f.name, cfg.figures, 'UniformOutput', false);
cfg.figures = cfg.figures(strcmp(names, 'fig_noise_switch_rollouts'));
assert(isscalar(cfg.figures));

t0 = tic;
run_dir = run_all_paper_analyses(cfg);
t_run = toc(t0);
assert(strcmp(run_dir, tmp));
M = load(fullfile(tmp, 'run_manifest.mat')); M = M.run_manifest;
assert(isequal(M.stages_run, {'noise_switch'}));
D = load(fullfile(tmp, 'noise_switch', 'noise_switch_rollouts_data.mat'));
assert(isequal({D.results.name}, cfg.conditions) && D.settings.n_copies == 3 && ...
    D.settings.t_rollout == 1 && D.settings.n == 60 && strcmp(D.settings.preset_name, cfg.preset_name));

t0 = tic;
R = make_all_paper_figures(cfg);
t_fig = toc(t0);
assert(isscalar(R) && R.ok, 'the figure entry failed: %s', R.err);
tags = [strcat('Fig_noise_switch_timeseries_', cfg.conditions, '_copies'), {'Fig_noise_switch_divergence'}];
out_dir = fullfile(tmpf, 'fig_noise_switch_rollouts');
assert(all(cellfun(@(tg) isfile(fullfile(out_dir, [tg '.png'])), tags)) && ...
    isfile(fullfile(out_dir, 'Fig_noise_switch_divergence_table.md')), ...
    'one time-series figure per condition of the set, the summary and its table');
for a = 1:numel(tags)
    assert(~any(startsWith(tags([1:a-1, a+1:end]), tags{a})), 'the tag %s is a prefix of another', tags{a});
end

fprintf('PASS: test_noise_switch_stage -- the size table in every mode, overrides and refusals, end to end (run %.0f s, figures %.0f s).\n', ...
    t_run, t_fig);

%% ------------------------------------------------------------------------
function v = arg(args, name)
% The value a name-value list gives name (the last one, as a call takes it).
k = find(strcmp(args(1:2:end), name), 1, 'last');
assert(~isempty(k), 'no argument %s', name);
v = args{2 * k};
end

function tf = throws(fn, id)
% fn throws that identifier.
try
    fn();
    tf = false;
catch err
    tf = strcmp(err.identifier, id);
    if ~tf; fprintf(2, 'expected %s, got %s: %s\n', id, err.identifier, err.message); end
end
end

function rmdir_if(d)
if isfolder(d); rmdir(d, 's'); end
end
