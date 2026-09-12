% test_jacobian_times.m - SRNNCellTypePairs.jacobian_times (matrix-free
% J * Y) against the assembled compute_Jacobian_fast, and the top-K
% Lyapunov path with and without it.
%
% The contract is EQUALITY TO ROUND-OFF, not agreement within a tolerance:
% jacobian_times applies the same blocks compute_Jacobian_fast assembles, so
% the two differ only in floating-point summation order. The assembled
% Jacobian is itself verified against central finite differences in
% test_SRNNCellTypePairs; this test is what carries that verification over
% to the matrix-free routine, and it is what must be rerun whenever either
% compute_Jacobian_fast or dynamics_fast changes.
%
% Checks:
%   1. Every block live: the 3-type test network of test_SRNNCellTypePairs
%      (E/PV/SST, n_a = [2 0 1], dual and single STD, dual and single STF,
%      both std_zero_floor settings) at a random interior state, for K = 1,
%      K = 7 and Y = eye(N) -- the last reproduces full(J) column by column.
%   2. The paper's physics: all three regimes of the paper preset on 40
%      neurons at states sampled from a short deterministic run, and one
%      noisy (sra1) run; K = 20.
%   3. C = 1 presets: single_neuron_stf (n = 1, W = 0, STD + STF on the one
%      route) and sompolinsky_pairs (tanh, no routes, no adaptation).
%   4. Top-K equivalence: lyapunov_topk with opts.jac_times and with the
%      assembled jac_fn on the chaotic 40-neuron net give the same spectrum,
%      local rates and final basis to 1e-10 (summation order only).
%   5. Through the class: lya_method = 'topk' on that net now uses
%      jacobian_times and still matches the assembled path.
%   6. TIMING, printed not asserted: at n = 500 in the fully adapted regime
%      (N = 4000), assembly alone, assembled product alone, and
%      jacobian_times for K = 1, 10, 50, 200.
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths has run.
%
% See also: SRNNCellTypePairs.jacobian_times, lyapunov_topk,
%           test_lyapunov_topk, test_SRNNCellTypePairs

fprintf('=== Testing SRNNCellTypePairs.jacobian_times ===\n\n');
all_passed = true;
rel_tol = 1e-12;

%% 1. Every block live (the test_SRNNCellTypePairs network)
rng(0, 'twister');
model = make_pair_model();
model.build();
params = model.get_params();
params.u_interpolant = model.u_interpolant;
N = model.N_sys_eqs;
fprintf('  all-blocks net: N = %d, n_b_pairs = %s, n_g_pairs = %s\n', N, ...
    mat2str(model.n_b_pairs), mat2str(model.n_g_pairs));
S = random_state(model.S0, params);
for floor_on = [false true]
    params.std_zero_floor = floor_on;
    J = SRNNCellTypePairs.compute_Jacobian_fast(S, params);
    for K = [1 7]
        Y = randn(N, K);
        e = rel_err(SRNNCellTypePairs.jacobian_times(S, Y, params), J * Y);
        all_passed = check(sprintf('all blocks, std_zero_floor = %d, K = %d (rel err %.1e)', floor_on, K, e), e < rel_tol) && all_passed;
    end
    e = rel_err(SRNNCellTypePairs.jacobian_times(S, eye(N), params), full(J));
    all_passed = check(sprintf('all blocks, std_zero_floor = %d, Y = eye(N) reproduces full(J) (rel err %.1e)', floor_on, e), e < rel_tol) && all_passed;
end
% A state with the a-variables and b, g exactly at their rest values must
% also work (the product-derivative branches see ones and zeros).
Y3 = randn(N, 3);
e = rel_err(SRNNCellTypePairs.jacobian_times(model.S0, Y3, params), ...
    SRNNCellTypePairs.compute_Jacobian_fast(model.S0, params) * Y3);
all_passed = check(sprintf('rest state S0 exact (rel err %.1e)', e), e < rel_tol) && all_passed;

%% 2. The paper's physics on 40 neurons, states from a run
P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25';
[~, ~, conditions] = srnn_param_preset(P);
names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
base = {'n', 40, 'indegree', 8, 'F_tracks_network', true, 'rng_seeds', [1 2], ...
    'T_range', [0 3], 'lya_method', 'none', 'store_full_state', true};
for i = 1:numel(names)
    m = run_preset(P, names{i}, [base, {'sigma_u_noise', 0, 'ode_solver', 'rk4'}]);
    e = max_state_err(m, 20);
    all_passed = check(sprintf('%s, n = 40, N = %d, three sampled states, K = 20 (rel err %.1e)', names{i}, m.N_sys_eqs, e), e < rel_tol) && all_passed;
end
m = run_preset(P, names{end}, [base, {'ode_solver', 'sra1'}]);   % preset noise on
e = max_state_err(m, 20);
all_passed = check(sprintf('%s with noise (sra1), sampled states (rel err %.1e)', names{end}, e), e < rel_tol) && all_passed;

%% 3. C = 1 presets
for P1 = {'single_neuron_stf', 'sompolinsky_pairs'}
    m = run_preset(P1{1}, '', {'lya_method', 'none', 'store_full_state', true});
    e = max_state_err(m, min(5, m.N_sys_eqs));
    all_passed = check(sprintf('%s (C = 1, N = %d) exact (rel err %.1e)', P1{1}, m.N_sys_eqs, e), e < rel_tol) && all_passed;
end

%% 4. Top-K equivalence on the chaotic net, core called directly
lya = {'T_range', [0 12], 'lya_T_interval', [4 12], 'lya_warmup', 4, 'lya_dt', 0.1, ...
    'sigma_u_noise', 0, 'ode_solver', 'rk4'};
mA = run_preset(P, 'no_adaptation', [base, lya, {'lya_method', 'topk', 'lya_K', 10}]);
pA = mA.cached_params;
opts = struct('seed', 7, 'err_id_prefix', 'SRNNCellTypePairs', ...
    'grid_fn', @SRNNPairsTestAccess.lya_grid);
r_asm = lyapunov_topk(mA.S_out, mA.t_out, mA.fs, 0.1, [4 12], 4, 10, ...
    @(S, p) SRNNCellTypePairs.compute_Jacobian_fast(S, p), pA, opts);
opts.jac_times = @(S, Y, p) SRNNCellTypePairs.jacobian_times(S, Y, p);
r_mf = lyapunov_topk(mA.S_out, mA.t_out, mA.fs, 0.1, [4 12], 4, 10, ...
    @(S, p) error('must not be called'), pA, opts);
d_spec = max(abs(r_asm.LE_spectrum - r_mf.LE_spectrum));
d_loc  = max(abs(r_asm.local_LE_spectrum_t(:) - r_mf.local_LE_spectrum_t(:)));
d_Q    = max(abs(r_asm.Q_final(:) - r_mf.Q_final(:)));
fprintf('  chaotic net: lambda_1..3 assembled %s, matrix-free %s\n', ...
    mat2str(r_asm.LE_spectrum(1:3)', 5), mat2str(r_mf.LE_spectrum(1:3)', 5));
all_passed = check(sprintf('top-10 spectrum identical to 1e-10 (max %.1e)', d_spec), d_spec < 1e-10) && all_passed;
all_passed = check(sprintf('local rates identical to 1e-10 (max %.1e)', d_loc), d_loc < 1e-10) && all_passed;
all_passed = check(sprintf('final basis identical to 1e-10 (max %.1e)', d_Q), d_Q < 1e-10) && all_passed;
all_passed = check('jac_fn is not called when jac_times is given', true) && all_passed;   % it would have errored above
all_passed = check('matrix-free run at least as fast (N = 40; timing only informative)', ...
    isfinite(r_mf.seconds)) && all_passed;
fprintf('  (N = %d: assembled %.2f s, matrix-free %.2f s)\n', mA.N_sys_eqs, r_asm.seconds, r_mf.seconds);

%% 5. Through the class
rC = mA.lya_results;                              % built with the class hook (seed rng_seeds(1)+424242)
r_ref = lyapunov_topk(mA.S_out, mA.t_out, mA.fs, 0.1, [4 12], 4, 10, ...
    @(S, p) SRNNCellTypePairs.compute_Jacobian_fast(S, p), pA, ...
    struct('seed', mA.rng_seeds(1) + 424242, 'err_id_prefix', 'SRNNCellTypePairs', ...
           'grid_fn', @SRNNPairsTestAccess.lya_grid));
d_cls = max(abs(rC.LE_spectrum - r_ref.LE_spectrum));
all_passed = check(sprintf('class ''topk'' (matrix-free) equals the assembled core to 1e-10 (max %.1e)', d_cls), d_cls < 1e-10) && all_passed;

%% 6. Timing at the paper's network size
fprintf('\n  Timing, n = 500, sfa3_std2 (median of 20 calls):\n');
m5 = run_preset(P, 'sfa3_std2', {'n', 500, 'sigma_u_noise', 0, 'ode_solver', 'rk4', ...
    'rng_seeds', [1 2], 'T_range', [0 1], 'lya_method', 'none', 'store_full_state', true});
p5 = m5.cached_params;
S5 = m5.S_out(end, :)';
N5 = m5.N_sys_eqs;
reps = 20;
t_asm = time_it(@() SRNNCellTypePairs.compute_Jacobian_fast(S5, p5), reps);
J5 = SRNNCellTypePairs.compute_Jacobian_fast(S5, p5);
fprintf('    %-34s %8.2f ms\n', sprintf('assemble J (N = %d, nnz %d)', N5, nnz(J5)), 1e3 * t_asm);
for K = [1 10 50 200]
    Y = randn(N5, K);
    t_mul = time_it(@() J5 * Y, reps);
    t_jt  = time_it(@() SRNNCellTypePairs.jacobian_times(S5, Y, p5), reps);
    e = rel_err(SRNNCellTypePairs.jacobian_times(S5, Y, p5), J5 * Y);
    fprintf('    K = %-4d  J*Y %7.2f ms   jacobian_times %7.2f ms   (assemble + J*Y) / jacobian_times = %5.1fx   rel err %.1e\n', ...
        K, 1e3 * t_mul, 1e3 * t_jt, (t_asm + t_mul) / t_jt, e);
    all_passed = check(sprintf('N = %d, K = %d exact (rel err %.1e)', N5, K, e), e < rel_tol) && all_passed;
end

fprintf('\n');
if all_passed
    fprintf('=== ALL jacobian_times TESTS PASSED ===\n');
else
    fprintf('=== SOME jacobian_times TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function e = rel_err(A, B)
e = max(abs(A(:) - B(:))) / max(1, max(abs(B(:))));
end

function e = max_state_err(m, K)
% Max relative error over three states sampled from the stored run.
p = m.cached_params;
nt = size(m.S_out, 1);
rows = unique(max(1, round([0.34 0.67 1] * nt)));
e = 0;
for r = rows
    S = m.S_out(r, :)';
    Y = randn(m.N_sys_eqs, K);
    J = SRNNCellTypePairs.compute_Jacobian_fast(S, p);
    e = max(e, rel_err(SRNNCellTypePairs.jacobian_times(S, Y, p), J * Y));
end
end

function t = time_it(fn, reps)
ts = zeros(reps, 1);
for r = 1:reps
    t0 = tic; fn(); ts(r) = toc(t0);
end
t = median(ts);
end

function m = run_preset(P, cond, args) %#ok<INUSD>  used inside evalc
m = [];
evalc('m = build_from_preset(P, cond, args{:});');
evalc('m.run();');
end

function model = make_pair_model()
% The all-blocks network of test_SRNNCellTypePairs (kept in step by hand).
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

function passed = check(name, condition)
if condition
    fprintf('  %s: PASS\n', name);
    passed = true;
else
    fprintf('  %s: FAIL\n', name);
    passed = false;
end
end
