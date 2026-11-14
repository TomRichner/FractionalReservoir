% test_noise_switch_rollouts.m - run_noise_switch_rollouts and
% fig_noise_switch_rollouts, end to end on a small network.
%
% The mu7chain preset at n = 60, two conditions (no adaptation and the
% chain), 3 copies, 2 s shared and 1 s of rollout, into a temp run
% directory. Checks:
%   1. The stage returns its directory; the .mat and git_provenance.txt exist.
%   2. Settings record the protocol; the rate and integrator are the sweeps'.
%   3. Blocks: no adaptation has x, r, theta only; the chain has a and b too,
%      b labeled as the releasable pool.
%   4. Before the switch every statistic says "one trajectory": RMS
%      difference and variance exactly 0, correlation 1. After it the RMS
%      difference of x is positive at the first stored sample.
%   5. rms_diff equals the root mean square over copy pairs, computed here
%      from an independent rerun of the copies' x at the end of the rollout
%      (which also checks the stage against the model directly).
%   6. The example series: two copies, 5 + 5 neurons, identical before the
%      switch, different after.
%   7. The fit fields are finite for x; lambda_1 is finite.
%   8. The figure writes one time-series figure per condition, the summary
%      and its table; no tag is a prefix of another; `conditions` selects;
%      a run directory without the stage errors.
%
% About 20 s, serial. Prints PASS/FAIL per check and a final banner.
% Assumes setup_paths.
%
% See also: run_noise_switch_rollouts, fig_noise_switch_rollouts,
%           test_noise_seed_schedule

fprintf('=== Testing run_noise_switch_rollouts + figure ===\n\n');
all_passed = true;
P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain';
conds = {'no_adaptation', 'sfa3_stdc2'};
small = {'n', 60, 'indegree', 20};
t_shared = 2; t_rollout = 1; M = 3;

run_dir = fullfile(tempdir, 'noise_switch_rollouts_test');
if exist(run_dir, 'dir'), rmdir(run_dir, 's'); end
mkdir(run_dir);
stage_dir = fullfile(run_dir, 'noise_switch_rollouts');

evalc(['d = run_noise_switch_rollouts(''preset_name'', P, ''conditions'', conds, ' ...
    '''t_shared'', t_shared, ''t_rollout'', t_rollout, ''n_copies'', M, ''n_show_E'', 5, ' ...
    '''n_show_I'', 5, ''model_args'', small, ''out_dir'', stage_dir);']);
mat_file = fullfile(stage_dir, 'noise_switch_rollouts_data.mat');
all_passed = check('stage returns its directory, with the .mat and the provenance', ...
    strcmp(d, stage_dir) && isfile(mat_file) && isfile(fullfile(stage_dir, 'git_provenance.txt'))) && all_passed;
D = load(mat_file);
S = D.settings;

ref = analysis_run_config('param_space', 'medium', srnn_param_preset(P));
all_passed = check('settings record the protocol', S.t_shared == t_shared && S.t_rollout == t_rollout && ...
    S.n_copies == M && isequal(S.rollout_seeds, S.noise_seed + (1:M)) && S.n == 60) && all_passed;
all_passed = check('rate and integrator are the sweeps'' for this preset', ...
    S.fs == ref.model.fs && strcmp(S.ode_solver, ref.model.ode_solver)) && all_passed;

R0 = D.results(1); R1 = D.results(2);
all_passed = check('conditions in the order asked', isequal({D.results.name}, conds)) && all_passed;
all_passed = check('no adaptation: x, r, theta only', isequal(R0.blocks, {'x', 'r', 'theta'})) && all_passed;
all_passed = check('chain: a and b too, b the releasable pool', ...
    isequal(R1.blocks, {'x', 'r', 'theta', 'a', 'b'}) && contains(R1.block_labels.b, 'releasable')) && all_passed;

before = R1.t <= t_shared; first_after = find(~before, 1);
ok = true;
for k = 1:numel(R1.blocks)
    st = R1.stats.(R1.blocks{k});
    ok = ok && all(st.rms_diff(before) == 0) && all(st.variance(before) == 0) && ...
        all(abs(st.corr(before) - 1) < 1e-12) && st.scale > 0;
end
all_passed = check('before the switch: RMS difference and variance 0, correlation 1', ok) && all_passed;
all_passed = check('x differs at the first stored sample after the switch', ...
    R1.stats.x.rms_diff(first_after) > 0) && all_passed;

% An independent rerun of the copies, straight from the model.
X_end = zeros(60, M);
for m = 1:M
    mdl = build_from_preset(P, 'sfa3_stdc2', 'rng_seeds', S.rng_seeds, 'T_range', S.T_range, ...
        'fs', S.fs, 'ode_solver', S.ode_solver, 'noise_seed', S.noise_seed, ...
        'noise_seed_schedule', [t_shared, S.rollout_seeds(m)], 'lya_method', 'none', ...
        'store_full_state', true, small{:});
    evalc('mdl.run();');
    X_end(:, m) = mdl.S_out(end, mdl.cached_params.state_layout.x)';
end
sq = [];
for i = 1:M - 1
    for j = i + 1:M
        sq(end + 1, :) = (X_end(:, i) - X_end(:, j))' .^ 2; %#ok<SAGROW>
    end
end
rms_pairs = sqrt(mean(sq, 'all'));
all_passed = check('rms_diff is the RMS over copy pairs and neurons of an independent rerun', ...
    abs(R1.stats.x.rms_diff(end) - rms_pairs) < 1e-12 * rms_pairs) && all_passed;

E = R1.example;
eb = E.t <= t_shared;
all_passed = check('example: two copies of 5 + 5 neurons per block', ...
    isequal(size(E.series.x, [1 3]), [10 2]) && isequal(E.type.x, [ones(1, 5), 2 * ones(1, 5)])) && all_passed;
all_passed = check('example copies identical before the switch, different after', ...
    isequal(E.series.x(:, eb, 1), E.series.x(:, eb, 2)) && ...
    ~isequal(E.series.x(:, end, 1), E.series.x(:, end, 2))) && all_passed;

F = R1.stats.x.fit;
all_passed = check('fit of x and lambda_1 are finite', all(isfinite([F.saturation, F.t_half, ...
    F.growth_rate, F.rise_exponent, R0.lambda_1, R1.lambda_1]))) && all_passed;

%% Figure
fig_dir = fullfile(run_dir, 'figs');
evalc('out = fig_noise_switch_rollouts(''run_dir'', run_dir, ''out_dir'', fig_dir, ''visible'', false);');
close(out.figs);
tags = out.tags;
want = [cellfun(@(c) sprintf('Fig_noise_switch_timeseries_%s_copies', c), conds, 'UniformOutput', false), ...
    {'Fig_noise_switch_divergence'}];
all_passed = check('one time-series figure per condition and the summary', isequal(tags, want)) && all_passed;
all_passed = check('every figure wrote its png; the table exists', ...
    all(cellfun(@(tg) isfile(fullfile(fig_dir, [tg '.png'])), tags)) && ...
    isfile(fullfile(fig_dir, 'Fig_noise_switch_divergence_table.md'))) && all_passed;
is_prefix = false;
for i = 1:numel(tags)
    for j = 1:numel(tags)
        is_prefix = is_prefix || (i ~= j && startsWith(tags{j}, tags{i}));
    end
end
all_passed = check('no tag is a prefix of another', ~is_prefix) && all_passed;

evalc('one = fig_noise_switch_rollouts(''run_dir'', run_dir, ''conditions'', {''sfa3_stdc2''}, ''save'', false, ''visible'', false);');
close(one.figs);
all_passed = check('conditions selects', numel(one.figs) == 2) && all_passed;

empty_dir = fullfile(run_dir, 'empty'); mkdir(empty_dir);
all_passed = check('a run directory without the stage errors', ...
    throws_id(@() fig_noise_switch_rollouts('run_dir', empty_dir, 'save', false, 'visible', false), ...
    'resolve_data_file:NotFound')) && all_passed;

rmdir(run_dir, 's');

%% Result
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');

%% Helpers
function tf = throws_id(fcn, id)
tf = false;
try
    fcn();
catch err
    tf = strcmp(err.identifier, id);
end
end

function passed = check(name, condition)
if condition
    fprintf('  %s: PASS\n', name);
    passed = true;
else
    fprintf('  %s: FAIL\n', name);
    passed = false;
end
end
