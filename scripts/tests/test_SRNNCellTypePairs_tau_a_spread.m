% test_SRNNCellTypePairs_tau_a_spread.m - Per-neuron SFA ladders.
%
% tau_a_spread (1 x C, log-normal spread of the ladder's two ends, geometric
% interpolation between) is drawn at build() into the read-only tau_a_matrix
% (1 x C cell, n_q x n_a(q) per type), and get_params() hands the physics
% params.tau_a_matrix at that size whether or not a spread was requested.
%
% Checks:
%   1. Defaults are inert: tau_a_matrix empty, params.tau_a_matrix is the
%      nominal row replicated, and a reference run (sfa3_std2, n = 40, seed
%      [1 2], captured before the feature existed) is BIT-IDENTICAL: final
%      state and top-10 spectrum, max abs diff exactly 0. (Skipped when the
%      reference file from that session is not on this machine.)
%   2. Shape contract: n_q x n_a(q); n_q x 0 without SFA; a type with zero
%      spread keeps the nominal row exactly while another type is spread.
%   3. Draw statistics at n = 900: at each end, log(tau_ik / tau_k) has mean
%      0 and std sigma within 3 s.e.; the two ends are uncorrelated.
%   4. Interpolation: interior rungs of 3- and 5-rung log ladders equal the
%      log-linear interpolation of the drawn endpoints to 1e-12; a
%      NON-logspaced nominal is reproduced exactly at sigma = 0 and keeps its
%      own log positions at sigma > 0.
%   5. K = 1 is the slow draw alone; K = 2 is the two endpoints.
%   6. Validation: too wide a spread errors TauSpreadTooWide; just under the
%      bound is accepted; negative / wrong-length spreads error.
%   7. Seed control: reproducible; tau_a_seed and rng_seeds both move it;
%      W, S0, u_ex and S_c_vec are untouched by the draw.
%   8. Jacobian at sigma = 0.2 on the all-blocks network: analytic vs central
%      finite differences < 2e-6, and jacobian_times == J * Y to 1e-12.
%   9. SRNN_ESN_reservoir builds with a spread.
%  10. The paper preset physics at n = 40 builds all three regimes with a
%      spread: no_adaptation has empty ladders, sfa1_std1 an n_q x 1 draw,
%      sfa3_std2 n_q x 3; and the stable regime's lambda_1 moves toward
%      -1/tau_slowest.
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths has run.
%
% See also: SRNNCellTypePairs, test_SRNNCellTypePairs_S_c_heterogeneity,
%           test_jacobian_times

fprintf('=== Testing SRNNCellTypePairs per-neuron SFA ladders (tau_a_spread) ===\n\n');
all_passed = true;
scratch = fullfile(tempdir, 'claude', 'C--Users-m218089-Desktop-github-repos-FractionalReservoir', ...
    '90c5825b-19c1-4546-90a4-703d562404ba', 'scratchpad', 'tau_spread_baseline.mat');
ladder = @(a, b, k) log_ladder(a, b, k);

%% 1. Defaults are inert
m0 = tiny_model();
m0.build();
p0 = m0.get_params();
all_passed = check('tau_a_matrix is empty with the properties at their defaults', ...
    isempty(m0.tau_a_matrix)) && all_passed;
all_passed = check('tau_a_spread completes to zeros(1, C)', ...
    isequal(m0.tau_a_spread, zeros(1, m0.n_cellTypes))) && all_passed;
ok = true;
for q = 1:m0.n_cellTypes
    want = repmat(m0.tau_a{q}, m0.n_per_type(q), 1);
    ok = ok && isequal(p0.tau_a_matrix{q}, want);
end
all_passed = check('params.tau_a_matrix is the nominal row replicated per neuron', ok) && all_passed;

P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25';
base = {'n', 40, 'indegree', 8, 'F_tracks_network', true, 'sigma_u_noise', 0, ...
    'ode_solver', 'rk4', 'rng_seeds', [1 2], 'T_range', [0 12], 'lya_T_interval', [4 12], ...
    'lya_warmup', 4, 'lya_dt', 0.1, 'store_full_state', true, 'lya_method', 'topk', 'lya_K', 10};
if exist(scratch, 'file')
    ref = load(scratch);
    mB = run_preset(P, 'sfa3_std2', base);
    d_S  = max(abs(mB.S_out(end, :) - ref.S_end));
    d_LE = max(abs(mB.lya_results.LE_spectrum - ref.LE));
    all_passed = check(sprintf('BIT-IDENTICAL to the pre-feature reference: final state (diff %.1e)', d_S), d_S == 0) && all_passed;
    all_passed = check(sprintf('BIT-IDENTICAL to the pre-feature reference: top-10 spectrum (diff %.1e)', d_LE), d_LE == 0) && all_passed;
else
    fprintf('  (no pre-feature reference file at %s; bit-identity check skipped)\n', scratch);
end

%% 2. Shape contract
m = tiny_model('tau_a_spread', [0.1 0 0.05]);
m.build();
M = m.tau_a_matrix;
all_passed = check('tau_a_matrix is a 1 x C cell', iscell(M) && isequal(size(M), [1 m.n_cellTypes])) && all_passed;
ok = true;
for q = 1:m.n_cellTypes
    ok = ok && isequal(size(M{q}), [m.n_per_type(q), m.n_a(q)]);
end
all_passed = check('each cell is n_q x n_a(q) (n_q x 0 for the type without SFA)', ok) && all_passed;
all_passed = check('a type with zero spread keeps the nominal row exactly', ...
    isequal(M{2}, repmat(m.tau_a{2}, m.n_per_type(2), 1))) && all_passed;
all_passed = check('a type with a spread differs from the nominal in every row', ...
    all(any(M{1} ~= m.tau_a{1}, 2)) && all(M{1}(:) > 0)) && all_passed;
p = m.get_params();
all_passed = check('params.tau_a_matrix is the realised draw when a spread is on', ...
    isequal(p.tau_a_matrix, M)) && all_passed;

%% 3. Draw statistics
sig = [0.10 0.05 0.20];
m = tiny_model('n', 900, 'indegree', 100, 'tau_a', {ladder(0.25, 10, 3), ladder(0.5, 5, 2), 2}, ...
    'tau_a_spread', sig);
m.build();
stat_ok = true; corr_ok = true;
for q = 1:m.n_cellTypes
    Mq = m.tau_a_matrix{q}; tau = m.tau_a{q}; nq = size(Mq, 1);
    lf = log(Mq(:, 1) / tau(1)); ls = log(Mq(:, end) / tau(end));
    for v = {lf, ls}
        stat_ok = stat_ok && abs(mean(v{1})) < 3 * sig(q) / sqrt(nq) && ...
            abs(std(v{1}) - sig(q)) < 3 * sig(q) / sqrt(2 * nq);
    end
    if numel(tau) >= 2
        corr_ok = corr_ok && abs(corr(lf, ls)) < 3 / sqrt(nq);
    end
end
all_passed = check('log(tau_ik / tau_k) at each end has mean 0 and std = spread (3 s.e.)', stat_ok) && all_passed;
all_passed = check('the fast and slow draws are uncorrelated', corr_ok) && all_passed;

%% 4. Interpolation
interp_ok = true;
for K = [3 5]
    m = tiny_model('tau_a', {ladder(0.25, 10, K), ladder(0.25, 10, K), []}, 'tau_a_spread', 0.1);
    m.build();
    for q = 1:2
        Mq = m.tau_a_matrix{q}; tau = m.tau_a{q};
        w = (log(tau) - log(tau(1))) / (log(tau(end)) - log(tau(1)));
        want = exp(log(Mq(:, 1)) .* (1 - w) + log(Mq(:, end)) .* w);
        interp_ok = interp_ok && max(abs(Mq(:) - want(:))) < 1e-12 * max(tau);
    end
end
all_passed = check('interior rungs are the log-linear interpolation of the drawn ends (3- and 5-rung)', interp_ok) && all_passed;

odd = [0.25 0.5 10];
m = tiny_model('tau_a', {odd, [], []}, 'tau_a_spread', 0);
m.build(); p = m.get_params();
all_passed = check('a non-logspaced nominal is reproduced exactly at zero spread', ...
    isequal(p.tau_a_matrix{1}, repmat(odd, m.n_per_type(1), 1))) && all_passed;
m = tiny_model('tau_a', {odd, [], []}, 'tau_a_spread', 0.1);
m.build(); Mq = m.tau_a_matrix{1};
w2 = (log(0.5) - log(0.25)) / (log(10) - log(0.25));
want = exp(log(Mq(:, 1)) * (1 - w2) + log(Mq(:, 3)) * w2);
all_passed = check('a non-logspaced nominal keeps its own log position under a spread', ...
    max(abs(Mq(:, 2) - want)) < 1e-12 * 10) && all_passed;

%% 5. K = 1 and K = 2
m = tiny_model('tau_a', {10, ladder(0.25, 10, 2), []}, 'tau_a_spread', 0.1, 'tau_a_seed', 3);
m.build();
seed_state = rng; rng(3, 'twister'); z = randn(m.n, 2); rng(seed_state);
i1 = m.type_indices{1}; i2 = m.type_indices{2};
all_passed = check('K = 1 is the slow draw alone', ...
    max(abs(m.tau_a_matrix{1} - 10 * exp(0.1 * z(i1, 2)))) < 1e-12) && all_passed;
want2 = [0.25 * exp(0.1 * z(i2, 1)), 10 * exp(0.1 * z(i2, 2))];
all_passed = check('K = 2 is the two endpoint draws', ...
    max(abs(m.tau_a_matrix{2}(:) - want2(:))) < 1e-12) && all_passed;

%% 6. Validation
bound = log(40) / (4 * sqrt(2));                       % ladder end ratio 40
all_passed = check('a spread past the ordering margin errors TauSpreadTooWide', ...
    throws_id(@() build_it(tiny_model('tau_a_spread', bound * 1.05)), ...
    'SRNNCellTypePairs:TauSpreadTooWide')) && all_passed;
all_passed = check('a spread just under the margin is accepted', ...
    ~throws_id(@() build_it(tiny_model('tau_a_spread', bound * 0.95)), ...
    'SRNNCellTypePairs:TauSpreadTooWide')) && all_passed;
all_passed = check('a negative spread errors', ...
    throws_id(@() build_it(tiny_model('tau_a_spread', -0.1)), 'SRNNCellTypePairs:InvalidParams')) && all_passed;
all_passed = check('a wrong-length spread errors', ...
    throws_id(@() build_it(tiny_model('tau_a_spread', [0.1 0.1])), 'SRNNCellTypePairs:InvalidParams')) && all_passed;

%% 7. Seed control and independence
args = {'tau_a_spread', 0.05, 'sigma_S_c', 0.05};
a = tiny_model(args{:}); a.build();
b = tiny_model(args{:}); b.build();
all_passed = check('the same configuration reproduces the same draw', ...
    isequal(a.tau_a_matrix, b.tau_a_matrix)) && all_passed;
c1 = tiny_model(args{:}, 'tau_a_seed', 11); c1.build();
c2 = tiny_model(args{:}, 'tau_a_seed', 12); c2.build();
all_passed = check('a different tau_a_seed gives a different draw', ...
    ~isequal(c1.tau_a_matrix, c2.tau_a_matrix)) && all_passed;
d1 = tiny_model(args{:}, 'rng_seeds', [5 6]); d1.build();
d2 = tiny_model(args{:}, 'rng_seeds', [7 8]); d2.build();
all_passed = check('a different network seed also moves the ladders', ...
    ~isequal(d1.tau_a_matrix, d2.tau_a_matrix)) && all_passed;
plain = tiny_model('rng_seeds', [3 4], 'sigma_S_c', 0.05);       plain.build();
het   = tiny_model('rng_seeds', [3 4], args{:});                 het.build();
all_passed = check('W is bit-identical with and without the draw', ...
    max(abs(full(plain.W) - full(het.W)), [], 'all') == 0) && all_passed;
all_passed = check('S0 is bit-identical with and without the draw', ...
    max(abs(plain.S0 - het.S0)) == 0) && all_passed;
all_passed = check('the external input is bit-identical too', ...
    max(abs(plain.u_ex - het.u_ex), [], 'all') == 0) && all_passed;
all_passed = check('the setpoint draw S_c_vec is bit-identical too', ...
    isequal(plain.S_c_vec, het.S_c_vec)) && all_passed;
e1 = tiny_model('tau_a_spread', [0.05 0 0], 'tau_a_seed', 11); e1.build();
e2 = tiny_model('tau_a_spread', [0.05 0 0], 'tau_a_seed', 12); e2.build();
all_passed = check('a type with zero spread is unaffected by the seed', ...
    isequal(e1.tau_a_matrix{2}, e2.tau_a_matrix{2}) && ~isequal(e1.tau_a_matrix{1}, e2.tau_a_matrix{1})) && all_passed;
n0 = tiny_model('tau_a', {[], [], []}, 'c', [0 0 0], 'tau_a_spread', 0.1); n0.build();
all_passed = check('a spread with no adapting type leaves tau_a_matrix empty', ...
    isempty(n0.tau_a_matrix)) && all_passed;

%% 8. Jacobian with a spread
rng(0, 'twister');
model = make_pair_model(0.2);
model.build();
params = model.get_params();
params.u_interpolant = model.u_interpolant;
S = random_state(model.S0, params);
for floor_on = [false true]
    params.std_zero_floor = floor_on;
    J = SRNNCellTypePairs.compute_Jacobian_fast(S, params);
    J_fd = finite_difference_jacobian(S, params, 1e-6);
    e_fd = max(abs(full(J(:)) - J_fd(:)));
    all_passed = check(sprintf('analytic Jacobian vs finite differences, spread 0.2, std_zero_floor = %d (%.1e)', floor_on, e_fd), e_fd < 2e-6) && all_passed;
    Y = randn(model.N_sys_eqs, 7);
    e_jt = max(abs(SRNNCellTypePairs.jacobian_times(S, Y, params) - J * Y), [], 'all') / max(1, max(abs(J * Y), [], 'all'));
    all_passed = check(sprintf('jacobian_times == J * Y with a spread (rel %.1e)', e_jt), e_jt < 1e-12) && all_passed;
end
all_passed = check('every neuron of a spread type has its own ladder in the Jacobian', ...
    numel(unique(model.tau_a_matrix{1}(:, end))) == model.n_per_type(1)) && all_passed;

%% 9. ESN subclass
sc = struct();
sc.E.E.std = struct('tau_rec', 1, 'tau_rel', 0.25); sc.E.I.std = sc.E.E.std;
esn = SRNN_ESN_reservoir('n', 20, 'indegree', 5, 'n_cellTypes', 2, 'cell_type_names', {'E', 'I'}, ...
    'f', [0.5 0.5], 'mu_tilde_relative', [0.1 -0.1], 'sigma_tilde_relative', [0.01 0.01], ...
    'tau_a', {ladder(0.25, 10, 3), []}, 'c', [0.1 0], 'synapse_config', sc, ...
    'tau_a_spread', 0.05, 'T_range', [0 0.5], 'fs', 200);
evalc('esn.build();');
all_passed = check('SRNN_ESN_reservoir builds with a spread and carries the ladders', ...
    esn.is_built && isequal(size(esn.tau_a_matrix{1}), [10 3])) && all_passed;

%% 10. The paper preset physics
[~, ~, conditions] = srnn_param_preset(P);
names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
short = {'n', 40, 'indegree', 8, 'F_tracks_network', true, 'sigma_u_noise', 0, 'ode_solver', 'rk4', ...
    'rng_seeds', [1 2], 'T_range', [0 3], 'lya_method', 'none', 'store_full_state', true, 'tau_a_spread', 0.05};
shape_ok = true;
for i = 1:numel(names)
    m = run_preset(P, names{i}, short);
    na = m.n_a;
    if all(na == 0)
        shape_ok = shape_ok && isempty(m.tau_a_matrix);
    else
        for q = 1:m.n_cellTypes
            shape_ok = shape_ok && isequal(size(m.tau_a_matrix{q}), [m.n_per_type(q), na(q)]);
        end
    end
    if isempty(m.tau_a_matrix)
        shapes = 'empty';
    else
        shapes = strjoin(cellfun(@(x) mat2str(size(x)), m.tau_a_matrix, 'UniformOutput', false), ' ');
    end
    fprintf('  %-14s n_a = %s, ladders %s\n', names{i}, mat2str(na), shapes);
end
all_passed = check('all three regimes build with a spread and the ladders have the right shapes', shape_ok) && all_passed;

mS = run_preset(P, 'sfa3_std2', [base, {'tau_a_spread', 0.1}]);
tau_slowest = max(cellfun(@(x) max([x(:); 0]), mS.tau_a_matrix));
fprintf('  sfa3_std2, spread 0.1: lambda_1 %+.4f (nominal-ladder run %+.4f), -1/tau_slowest = %+.4f, band lambda_1 - lambda_10 = %.4f\n', ...
    mS.lya_results.LLE, mB.lya_results.LLE, -1 / tau_slowest, mS.lya_results.LE_spectrum(1) - mS.lya_results.LE_spectrum(10));
all_passed = check('with a spread the stable regime''s lambda_1 moves toward -1/tau_slowest', ...
    mS.lya_results.LLE > mB.lya_results.LLE && mS.lya_results.LLE < 0) && all_passed;

fprintf('\n');
if all_passed
    fprintf('=== ALL tau_a_spread TESTS PASSED ===\n');
else
    fprintf('=== SOME tau_a_spread TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function m = tiny_model(varargin)
% Three cell types, SFA on two of them, so nothing here depends on E/I.
m = SRNNCellTypePairs('n', 60, 'indegree', 10, ...
    'n_cellTypes', 3, 'cell_type_names', {'A', 'B', 'C'}, ...
    'f', [0.4 0.3 0.3], ...
    'mu_tilde_relative', [0.1 -0.2 0.05], 'sigma_tilde_relative', [0.01 0.02 0.01], ...
    'tau_a', {log_ladder(0.25, 10, 3), log_ladder(0.25, 10, 2), []}, ...
    'c', [0.1 0.1 0], 'T_range', [0 0.2], 'fs', 200, 'lya_method', 'none', ...
    varargin{:});
end

function build_it(m) %#ok<INUSD>  used inside evalc
evalc('m.build();');
end

function m = run_preset(P, cond, args) %#ok<INUSD>  used inside evalc
m = [];
evalc('m = build_from_preset(P, cond, args{:});');
evalc('m.run();');
end

function model = make_pair_model(spread)
% The all-blocks network of test_SRNNCellTypePairs, plus a ladder spread.
synapse_config = struct();
synapse_config.E.E.std = struct('tau_rec', [0.3 1], 'tau_rel', 0.2);
synapse_config.E.PV.stf = struct('tau_dec', 0.5, 'tau_fac', 0.25, 'G', 2);
synapse_config.E.SST.std = struct('tau_rec', 0.8, 'tau_rel', 0.25);
synapse_config.E.SST.stf = struct( ...
    'tau_dec', [0.4 1.2], 'tau_fac', [0.2 0.3], 'G', [1.5 2]);
synapse_config.PV.PV.std = struct('tau_rec', 0.6, 'tau_rel', 0.3);
synapse_config.SST.PV.stf = struct('tau_dec', 0.9, 'tau_fac', 0.4, 'G', 1.8);
model = SRNNCellTypePairs( ...
    'n', 12, 'indegree', 4, ...
    'n_cellTypes', 3, 'cell_type_names', {'E', 'PV', 'SST'}, ...
    'f', [0.5 0.3 0.2], ...
    'mu_tilde_relative', [0.1 -0.2 -0.08], ...
    'sigma_tilde_relative', [0.01 0.02 0.015], ...
    'tau_a', {[0.25 1], [], 2}, ...
    'c', [0.05 0 0.03], 'synapse_config', synapse_config, ...
    'tau_a_spread', spread, ...
    'T_range', [0 0.1], 'fs', 200, 'lya_method', 'none', ...
    'store_full_state', true);
end

function S = random_state(S0, params)
S = S0;
for q = 1:params.n_cellTypes
    S(params.state_layout.a{q}) = 0.3 * rand(numel(params.state_layout.a{q}), 1);
end
for pre = 1:params.n_cellTypes
    npre = params.n_per_type(pre);
    for post = 1:params.n_cellTypes
        row_b = params.state_layout.b{pre, post};
        S(row_b) = 0.4 + 0.6 * rand(numel(row_b), 1);
        ng = params.n_g_pairs(pre, post);
        if ng > 0
            G = repmat(params.G{pre, post}, npre, 1);
            row_g = params.state_layout.g{pre, post};
            S(row_g) = 1 + (G(:) - 1) .* rand(numel(row_g), 1);
        end
    end
end
S(params.state_layout.x) = 0.5 * randn(params.n, 1);
end

function J = finite_difference_jacobian(S, params, h)
N = numel(S);
J = zeros(N, N);
for k = 1:N
    plus = S; minus = S;
    plus(k) = plus(k) + h;
    minus(k) = minus(k) - h;
    J(:, k) = (SRNNCellTypePairs.dynamics_fast(0, plus, params) - ...
        SRNNCellTypePairs.dynamics_fast(0, minus, params)) / (2 * h);
end
end

function ok = throws_id(fn, id)
ok = false;
try
    fn();
catch ME
    ok = strcmp(ME.identifier, id);
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
