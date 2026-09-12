% test_lyapunov_topk.m - The top-K Lyapunov method (lya_method = 'topk',
% shared core lyapunov_topk) against its three references: itself at a larger
% K, the verified full-spectrum 'qr' path, and Benettin's K = 1 estimate.
%
% Two networks, both built from the paper's preset physics on 40 neurons:
%   A  no_adaptation, seed [1 2]: CHAOTIC (two positive exponents, N = 40).
%      Alignment is fast here, so finite-time exponents from different
%      initial bases agree and the ode45 QR is a fair reference. Found by
%      scanning seeds: most 40-neuron no-adaptation nets are dead (-10) or
%      marginal, which is itself a finding (the numerics report, sec. 4.4).
%   B  sfa3_std2, same seed: STABLE, N = 320, every Jacobian block live.
%      Used for the exact checks that do not care about regime.
%
% Claims checked:
%   1. NESTING (exact, on the UNSORTED columns). The first 5 raw columns of
%      a K = N run equal the K = 5 run to round-off at every segment: the
%      initial basis is a column-sequential QR of a column-major Gaussian
%      draw, propagation is linear and column-wise, thin QR is
%      column-sequential. (The SORTED outputs need not nest: in a degenerate
%      band the top-5 by final value are not the first 5 raw columns.)
%   2. FULL SPECTRUM vs 'qr' on A: the top 5 and their partial sums within
%      0.1 of the ode45 variational integration (same window, same lya_dt,
%      different initial basis and integrator). On a STABLE net this
%      comparison is ill-posed: the slow band is nearly degenerate, alignment
%      takes ~1/(lambda_1 - lambda_2) ~ 100 s, and 8 s finite-time exponents
%      differ by 0.2 between any two bases -- which is why A is chaotic.
%   3. LIOUVILLE, on A and B and on SRNNModel2: the sum of all N exponents
%      equals the time-average of trace(J) over the accumulation window. An
%      identity for the exact flow (volume contraction = divergence), so it
%      checks the propagator and the accumulation end to end, in any regime,
%      independently of 'qr'. This is the decisive test.
%   4. lambda_1 vs Benettin on A, same run, within 0.1 (finite perturbation
%      vs tangent vector over 8 s of chaos).
%   5. Orthonormality after re-orthonormalisation, finite local rates (i.e.
%      positive diag(R)), conditioning recorded and >= 1.
%   6. WINDOW CONTRACT, as test_lya_window pins for the other two methods.
%   7. Derived quantities: LLE = lambda_1; h_KS = sum of positives (> 0 on A);
%      D_KY on A agrees with the 'qr' branch's; K = 1 on A is unresolved
%      (lambda_1 > 0 so the sum never crosses); B resolves to 0.
%   8. Runs on a NOISY trajectory (sra1), all finite.
%   9. SRNNModel2: Liouville and lambda_1 vs its Benettin.
%  10. Argument errors: K > N, non-integer, negative -> :InvalidLyapunovK;
%      lya_K = 0 means N.
%  11. TIMING TABLE, printed not asserted.
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths has run.
%
% See also: lyapunov_topk, test_benettin_vs_qr, test_lya_window,
%           docs/notes/Lyapunov_estimation_methods.md

fprintf('=== Testing lya_method = ''topk'' ===\n\n');
all_passed = true;

P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25';
base = {'n', 40, 'indegree', 8, 'F_tracks_network', true, 'sigma_u_noise', 0, ...
    'ode_solver', 'rk4', 'rng_seeds', [1 2], 'T_range', [0 12], ...
    'lya_T_interval', [4 12], 'lya_warmup', 4, 'lya_dt', 0.1, 'store_full_state', true};
A = @(varargin) run_preset(P, 'no_adaptation', base, varargin{:});
B = @(varargin) run_preset(P, 'sfa3_std2',     base, varargin{:});

mA = A('lya_method', 'topk', 'lya_K', 0);  rA = mA.lya_results;  NA = mA.N_sys_eqs;
mB = B('lya_method', 'topk', 'lya_K', 0);  rB = mB.lya_results;  NB = mB.N_sys_eqs;
fprintf('  A: N = %d, top 4 %s, %d positive, D_KY %.2f\n', NA, mat2str(rA.LE_spectrum(1:4)', 3), rA.n_positive, rA.D_KY);
fprintf('  B: N = %d, lambda_1 %+.4f, %d positive, D_KY %.2f\n', NB, rB.LLE, rB.n_positive, rB.D_KY);
all_passed = check('network A is chaotic (>= 2 positive exponents)', rA.n_positive >= 2) && all_passed;

%% 1. Nesting on the raw columns
m5 = A('lya_method', 'topk', 'lya_K', 5);  r5 = m5.lya_results;
rawN = unsort(rA);  raw5 = unsort(r5);
d_nest = max(abs(raw5(:) - reshape(rawN(:, 1:5), [], 1)));
all_passed = check(sprintf('raw columns of K = 5 equal the first 5 of K = N at every segment (max %.1e)', d_nest), ...
    d_nest < 1e-10) && all_passed;

%% 2. Full spectrum vs the verified 'qr' path (chaotic net)
mQ = A('lya_method', 'qr');  rQ = mQ.lya_results;
d_top = max(abs(rA.LE_spectrum(1:5) - rQ.LE_spectrum(1:5)));
d_sum = max(abs(cumsum(rA.LE_spectrum(1:5)) - cumsum(rQ.LE_spectrum(1:5))));
fprintf('  (topk vs qr, top 5: %s vs %s)\n', mat2str(rA.LE_spectrum(1:5)', 3), mat2str(rQ.LE_spectrum(1:5)', 3));
all_passed = check(sprintf('top 5 match the ode45 QR within 0.1 (max %.3f)', d_top), d_top < 0.1) && all_passed;
all_passed = check(sprintf('partial sums of the top 5 within 0.1 (max %.3f)', d_sum), d_sum < 0.1) && all_passed;
all_passed = check('same segment grid as qr', isequal(rA.t_lya, rQ.t_lya)) && all_passed;

%% 3. Liouville
[dl, sA, trA] = liouville(mA, @SRNNCellTypePairs.compute_Jacobian_fast, 4, 12);
fprintf('  (A: sum %.4f vs mean trace(J) %.4f, rel %.1e)\n', sA, trA, dl);
all_passed = check('A: sum of all exponents equals the mean divergence (rel < 1e-3)', dl < 1e-3) && all_passed;
[dl, sB, trB] = liouville(mB, @SRNNCellTypePairs.compute_Jacobian_fast, 4, 12);
fprintf('  (B: sum %.4f vs mean trace(J) %.4f, rel %.1e)\n', sB, trB, dl);
all_passed = check('B: sum of all exponents equals the mean divergence (rel < 1e-3)', dl < 1e-3) && all_passed;
[dl, sQ] = deal(abs(sum(rQ.LE_spectrum) - trA) / abs(trA), sum(rQ.LE_spectrum));
fprintf('  (for reference, the qr branch on A: sum %.4f, rel %.1e)\n', sQ, dl);

%% 4. lambda_1 vs Benettin
mb = A('lya_method', 'benettin');
fprintf('  (lambda_1: topk %+.4f, benettin %+.4f, qr %+.4f)\n', rA.LLE, mb.lya_results.LLE, rQ.LE_spectrum(1));
all_passed = check('lambda_1 agrees with Benettin within 0.1', abs(rA.LLE - mb.lya_results.LLE) < 0.1) && all_passed;

%% 5. Orthonormality, signs, conditioning
all_passed = check('orthonormality defect after QR < 1e-12 (A and B)', rA.orth_defect_max < 1e-12 && rB.orth_defect_max < 1e-12) && all_passed;
all_passed = check('all local rates finite (positive diag(R))', all(isfinite(rA.local_LE_spectrum_t(:))) && all(isfinite(rB.local_LE_spectrum_t(:)))) && all_passed;
all_passed = check('conditioning recorded, finite and >= 1', ...
    all(isfinite(rB.cond_t)) && all(rB.cond_t >= 1) && rB.cond_max == max(rB.cond_t)) && all_passed;

%% 6. Window contract
first_acc = find(~isnan(rB.finite_LE_spectrum_t(:, 1)), 1);
last_acc  = find(~isnan(rB.finite_LE_spectrum_t(:, 1)), 1, 'last');
all_passed = check('iteration starts lya_warmup before the window', abs(rB.t_lya(1) - 0) < 1e-9) && all_passed;
all_passed = check('accumulation starts at lya_T_interval(1)', abs(rB.t_lya(first_acc) - 4) < 1e-9) && all_passed;
all_passed = check('nothing accumulated past lya_T_interval(2)', rB.t_lya(last_acc) + 0.1 <= 12 + 1e-9) && all_passed;
all_passed = check('lya_dt resolved to 0.1 and recorded', rB.lya_dt == 0.1 && rB.lya_fs == 10) && all_passed;

%% 7. Derived quantities
all_passed = check('LLE is the largest exponent', rA.LLE == rA.LE_spectrum(1)) && all_passed;
all_passed = check('h_KS is the sum of the positive exponents, > 0 on A, bits = nats / ln 2', ...
    rA.h_KS > 0 && abs(rA.h_KS - sum(rA.LE_spectrum(rA.LE_spectrum > 0))) < 1e-12 && ...
    abs(rA.h_KS_bits - rA.h_KS / log(2)) < 1e-12) && all_passed;
D_qr = SRNNPairsTestAccess.kaplan_yorke(rQ.LE_spectrum);
fprintf('  (D_KY on A: topk %.3f, qr %.3f)\n', rA.D_KY, D_qr);
all_passed = check('D_KY on A agrees with the qr branch (within 0.2)', rA.D_KY_resolved && abs(rA.D_KY - D_qr) < 0.2) && all_passed;
m1 = A('lya_method', 'topk', 'lya_K', 1);
all_passed = check('K = 1 on A: D_KY unresolved (NaN, flag false)', ~m1.lya_results.D_KY_resolved && isnan(m1.lya_results.D_KY)) && all_passed;
all_passed = check('B: D_KY = 0 and resolved (stable)', rB.D_KY_resolved && rB.D_KY == 0 && rB.h_KS == 0) && all_passed;

%% 8. Noisy trajectory
mn = B('lya_method', 'topk', 'lya_K', 5, 'sigma_u_noise', 0.05, 'ode_solver', 'sra1');
all_passed = check('runs on a noisy (sra1) trajectory with finite results', ...
    all(isfinite(mn.lya_results.LE_spectrum)) && numel(mn.lya_results.LE_spectrum) == 5) && all_passed;

%% 9. SRNNModel2
m2c = {'n', 40, 'indegree', 20, 'n_a_E', 3, 'n_b_E', 1, 'level_of_chaos', 3, 'fs', 200, ...
    'ode_solver', 'rk4', 'T_range', [0 12], 'lya_T_interval', [4 12], 'lya_warmup', 4, ...
    'lya_dt', 0.1, 'store_full_state', true};
s_top = SRNNModel2(m2c{:}, 'lya_method', 'topk', 'lya_K', 0); s_top.build(); evalc('s_top.run();');
[dl, s2, tr2] = liouville(s_top, @SRNNModel2.compute_Jacobian_fast, 4, 12);
fprintf('  (SRNNModel2: N = %d, sum %.4f vs mean trace(J) %.4f, rel %.1e)\n', s_top.N_sys_eqs, s2, tr2, dl);
all_passed = check('SRNNModel2: Liouville (rel < 1e-3)', dl < 1e-3) && all_passed;
s_b = SRNNModel2(m2c{:}, 'lya_method', 'benettin'); s_b.build(); evalc('s_b.run();');
fprintf('  (SRNNModel2 lambda_1: topk %+.4f, benettin %+.4f)\n', s_top.lya_results.LLE, s_b.lya_results.LLE);
all_passed = check('SRNNModel2: lambda_1 within 0.1 of Benettin', abs(s_top.lya_results.LLE - s_b.lya_results.LLE) < 0.1) && all_passed;
all_passed = check('SRNNModel2 error prefix on a bad K', ...
    throws_id(@() run_srnn2(m2c, 'lya_method', 'topk', 'lya_K', s_top.N_sys_eqs + 1), 'SRNNModel:InvalidLyapunovK')) && all_passed;

%% 10. Argument errors
all_passed = check('K > N errors', throws_id(@() A('lya_method', 'topk', 'lya_K', NA + 1), 'SRNNCellTypePairs:InvalidLyapunovK')) && all_passed;
all_passed = check('non-integer K errors', throws_id(@() A('lya_method', 'topk', 'lya_K', 2.5), 'SRNNCellTypePairs:InvalidLyapunovK')) && all_passed;
all_passed = check('negative K errors', throws_id(@() A('lya_method', 'topk', 'lya_K', -1), 'SRNNCellTypePairs:InvalidLyapunovK')) && all_passed;
all_passed = check('lya_K = 0 means all N', rA.K == NA && numel(rA.LE_spectrum) == NA) && all_passed;

%% 11. Timing (printed, not asserted)
fprintf('\n-- timing, whole run() including the 12 s trajectory (identical for every method) --\n');
fprintf('  %-28s %8s\n', 'method', 'seconds');
fprintf('  A (N = %d):\n', NA);
for spec = {{'benettin', []}, {'topk', 1}, {'topk', 5}, {'topk', 0}, {'qr', []}}
    [name, t] = time_method(A, spec{1}{1}, spec{1}{2}, NA);
    fprintf('  %-28s %8.2f\n', name, t);
end
fprintf('  B (N = %d):\n', NB);
for spec = {{'benettin', []}, {'topk', 1}, {'topk', 10}, {'topk', 50}, {'topk', 0}}
    [name, t] = time_method(B, spec{1}{1}, spec{1}{2}, NB);
    fprintf('  %-28s %8.2f\n', name, t);
end
fprintf('  %-28s %8s\n', 'B qr (ode45)', 'skipped (minutes at N = 320)');
C = @(varargin) run_preset(P, 'sfa3_std2', [base, {'n', 100, 'indegree', 20}], varargin{:});
mC = C('lya_method', 'topk', 'lya_K', 1);
fprintf('  C = sfa3_std2 at n = 100 (N = %d):\n', mC.N_sys_eqs);
for spec = {{'benettin', []}, {'topk', 10}, {'topk', 50}}
    [name, t] = time_method(C, spec{1}{1}, spec{1}{2}, mC.N_sys_eqs);
    fprintf('  %-28s %8.2f\n', name, t);
end

%% Summary
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');

%% ------------------------------------------------------------------------
function m = run_preset(P, cond, base, varargin) %#ok<INUSD>  used inside evalc
m = [];
evalc('m = build_from_preset(P, cond, base{:}, varargin{:});');
evalc('m.run();');
end

function m = run_srnn2(common, varargin)
m = SRNNModel2(common{:}, varargin{:});
m.build();
evalc('m.run();');
end

function raw = unsort(r)
% Undo the descending sort so column j is the j-th basis vector's rates.
raw = zeros(size(r.local_LE_spectrum_t));
raw(:, r.sort_idx) = r.local_LE_spectrum_t;
end

function [rel, s, tr_mean] = liouville(m, jac_fn, t0, t1)
% Sum of the exponents vs the mean divergence over the window [t0, t1].
r = m.lya_results;
params = m.cached_params;
idx = find(m.t_out >= t0 & m.t_out <= t1);
tr = zeros(numel(idx), 1);
for k = 1:numel(idx)
    J = jac_fn(m.S_out(idx(k), :)', params);
    tr(k) = full(sum(diag(J)));
end
tr_mean = mean(tr);
s = sum(r.LE_spectrum);
rel = abs(s - tr_mean) / abs(tr_mean);
end

function [name, s] = time_method(runner, method, K, N)
args = {'lya_method', method};
name = method;
if ~isempty(K)
    args = [args, {'lya_K', K}];
    if K == 0; name = sprintf('topk K = N (%d)', N); else; name = sprintf('topk K = %d', K); end
end
if strcmp(method, 'qr'); name = 'qr (ode45, K = N)'; end
t0 = tic;
runner(args{:});
s = toc(t0);
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
