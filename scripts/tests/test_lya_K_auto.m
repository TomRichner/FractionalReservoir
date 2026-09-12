% test_lya_K_auto.m - The top-K retry on an unresolved dimension, and the
% lya_summary bundle both estimators fill.
%
% Network: the chaotic 40-neuron no-adaptation net of test_lyapunov_topk
% (seed [1 2], two positive exponents, D_KY ~ 2.3, N = 40), and the stable
% sfa3_std2 net (N = 320) for the block fractions.
%
% Checks:
%   1. lya_K = 15 resolves D_KY at once: K_used 15, retries 0.
%   2. lya_K = 2, lya_K_max = 8: unresolved at 2 (two positive exponents),
%      retried to 4, resolved there: K_used 4, retries 1; the first two
%      exponents equal the K = 2 run's (nesting).
%   3. lya_K = 2, lya_K_max = 2: stays unresolved, flags it, no retry.
%   4. lya_K_auto = false never retries.
%   5. lya_K = 0 (all N) never retries.
%   6. lya_summary: every field of lya_summary_fields present and scalar;
%      LLE equals lya_results.LLE; block fractions sum to 1; transient
%      scalars finite and in range; the local series covers the
%      accumulation window only; D_KY_resolved and n_positive_at_K are 0/1.
%   7. Benettin path: lya_summary gives LLE, lambda_1_drift and the
%      transient scalars, NaN for the spectrum-derived ones and the blocks.
%   8. transient_divergence on a synthetic series with known runs.
%   9. Stable net: leading vector lives in the SFA block (lead_frac_sfa is
%      the largest fraction), the sign of the tau_max claim.
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths has run.
%
% See also: lyapunov_topk, test_lyapunov_topk, SRNNCellTypePairs.lya_summary

fprintf('=== Testing lya_K_auto and lya_summary ===\n\n');
all_passed = true;

P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25';
base = {'n', 40, 'indegree', 8, 'F_tracks_network', true, 'sigma_u_noise', 0, ...
    'ode_solver', 'rk4', 'rng_seeds', [1 2], 'T_range', [0 12], ...
    'lya_T_interval', [4 12], 'lya_warmup', 4, 'lya_dt', 0.05, 'store_full_state', true};
A = @(varargin) run_preset(P, 'no_adaptation', base, varargin{:});
B = @(varargin) run_preset(P, 'sfa3_std2',     base, varargin{:});

%% 1-5. The retry
m15 = A('lya_method', 'topk', 'lya_K', 15);
r = m15.lya_results;
fprintf('  A: N = %d, K = 15: D_KY %.3f, n_pos %d\n', m15.N_sys_eqs, r.D_KY, r.n_positive);
all_passed = check('K = 15 resolves at once (K_used 15, retries 0)', ...
    r.D_KY_resolved && r.K == 15 && r.retries == 0 && r.K_requested == 15) && all_passed;

m2 = A('lya_method', 'topk', 'lya_K', 2, 'lya_K_auto', false);
m2a = A('lya_method', 'topk', 'lya_K', 2, 'lya_K_max', 8);
r2 = m2.lya_results; r2a = m2a.lya_results;
all_passed = check('K = 2 without auto is unresolved (two positive exponents)', ...
    ~r2.D_KY_resolved && isnan(r2.D_KY) && r2.retries == 0) && all_passed;
all_passed = check(sprintf('K = 2 with auto retries to 4 and resolves (K_used %d, retries %d)', r2a.K, r2a.retries), ...
    r2a.D_KY_resolved && r2a.K == 4 && r2a.retries == 1 && r2a.K_requested == 2) && all_passed;
d_nest = max(abs(r2a.LE_spectrum(1:2) - r2.LE_spectrum(1:2)));
all_passed = check(sprintf('the retried run nests the K = 2 exponents (max %.1e)', d_nest), d_nest < 1e-10) && all_passed;

m2c = A('lya_method', 'topk', 'lya_K', 2, 'lya_K_max', 2);
all_passed = check('a cap equal to K stays unresolved and flags it', ...
    ~m2c.lya_results.D_KY_resolved && m2c.lya_results.retries == 0 && m2c.lya_results.K == 2) && all_passed;

mN = A('lya_method', 'topk', 'lya_K', 0);
all_passed = check('lya_K = 0 (all N) never retries', mN.lya_results.retries == 0 && mN.lya_results.K == mN.N_sys_eqs) && all_passed;

%% 6. lya_summary on the top-K run
S = m15.lya_summary();
names = SRNNCellTypePairs.lya_summary_fields();
have = all(cellfun(@(f) isfield(S, f) && isscalar(S.(f)), names));
all_passed = check(sprintf('lya_summary has all %d scalar fields', numel(names)), have) && all_passed;
all_passed = check('LLE equals lya_results.LLE; gap = lambda_1 - lambda_2', ...
    S.LLE == r.LLE && abs(S.lambda_gap - (r.LE_spectrum(1) - r.LE_spectrum(2))) < 1e-12) && all_passed;
all_passed = check('spectrum-derived fields copied (n_positive, h_KS, D_KY, K_used, cond_max)', ...
    S.n_positive == r.n_positive && S.h_KS_bits == r.h_KS_bits && S.D_KY == r.D_KY && ...
    S.K_used == 15 && S.cond_max == r.cond_max && S.D_KY_resolved == 1 && S.n_positive_at_K == 0) && all_passed;
fsum = S.lead_frac_x + S.lead_frac_sfa + S.lead_frac_std + S.lead_frac_stf;
all_passed = check(sprintf('block fractions sum to 1 (%.12f); no adaptation -> all in x', fsum), ...
    abs(fsum - 1) < 1e-12 && S.lead_frac_x == 1) && all_passed;
all_passed = check('transient scalars finite and in range', ...
    S.frac_local_positive >= 0 && S.frac_local_positive <= 1 && isfinite(S.p95_finite_0p2s) && ...
    (isnan(S.mean_positive_excursion_s) || S.mean_positive_excursion_s > 0)) && all_passed;
all_passed = check('local series covers the accumulation window only', ...
    numel(S.local_rate_lead) == numel(S.t_lya_lead) && min(S.t_lya_lead) >= 4 - 1e-9 && ...
    max(S.t_lya_lead) < 12 && numel(S.local_rate_lead) == nnz(~isnan(r.finite_LE_spectrum_t(:, 1)))) && all_passed;
fw = r.finite_LE_spectrum_t(~isnan(r.finite_LE_spectrum_t(:, 1)), 1);
all_passed = check('lambda_1_drift is end minus three-quarter finite value', ...
    abs(S.lambda_1_drift - (fw(end) - fw(round(0.75 * numel(fw))))) < 1e-12) && all_passed;
fprintf('  (chaotic: frac positive %.3f, p95 over 0.2 s %.3f, mean excursion %.3f s, drift %+.4f)\n', ...
    S.frac_local_positive, S.p95_finite_0p2s, S.mean_positive_excursion_s, S.lambda_1_drift);

%% 7. Benettin path
mb = A('lya_method', 'benettin');
Sb = mb.lya_summary();
all_passed = check('benettin: LLE, drift and transient scalars filled', ...
    Sb.LLE == mb.lya_results.LLE && isfinite(Sb.lambda_1_drift) && ...
    isfinite(Sb.frac_local_positive) && isfinite(Sb.p95_finite_0p2s)) && all_passed;
all_passed = check('benettin: spectrum-derived fields and blocks are NaN', ...
    isnan(Sb.n_positive) && isnan(Sb.h_KS_bits) && isnan(Sb.D_KY) && isnan(Sb.K_used) && ...
    isnan(Sb.lead_frac_x) && isnan(Sb.lambda_gap)) && all_passed;
all_passed = check('benettin: local series present', ...
    ~isempty(Sb.local_rate_lead) && numel(Sb.local_rate_lead) == numel(Sb.t_lya_lead)) && all_passed;

%% 8. transient_divergence on a synthetic series
dt = 0.05;
loc = [-ones(10, 1); ones(4, 1); -ones(6, 1); ones(2, 1); -ones(8, 1)];   % 30 segments, runs of 4 and 2
T = SRNNCellTypePairs.transient_divergence(loc, dt, 0.2);
all_passed = check('synthetic: fraction positive 6/30', abs(T.frac_positive - 6 / 30) < 1e-12) && all_passed;
all_passed = check('synthetic: mean excursion (4 + 2)/2 segments = 0.15 s', abs(T.mean_excursion_s - 0.15) < 1e-12) && all_passed;
p95_want = prctile(movmean(loc, 4, 'Endpoints', 'discard'), 95);   % one full +1 window of 27
all_passed = check(sprintf('synthetic: window is 4 segments and p95 is the sliding-window percentile (%.3f)', T.p95_finite), ...
    T.window_segments == 4 && abs(T.p95_finite - p95_want) < 1e-12 && T.p95_finite > 0.5 && T.p95_finite <= 1) && all_passed;
T0 = SRNNCellTypePairs.transient_divergence(-ones(30, 1), dt, 0.2);
all_passed = check('synthetic: no positives -> fraction 0, excursion NaN, p95 = -1', ...
    T0.frac_positive == 0 && isnan(T0.mean_excursion_s) && abs(T0.p95_finite + 1) < 1e-12) && all_passed;
Ts = SRNNCellTypePairs.transient_divergence(ones(2, 1), dt, 0.2);
all_passed = check('synthetic: shorter than a window -> p95 NaN', isnan(Ts.p95_finite) && Ts.frac_positive == 1) && all_passed;

%% 9. Stable net: the leading direction is in the SFA block
mS = B('lya_method', 'topk', 'lya_K', 15);
SS = mS.lya_summary();
fprintf('  B (stable, N = %d): lambda_1 %+.4f, lead fractions x %.3f sfa %.3f std %.3f\n', ...
    mS.N_sys_eqs, SS.LLE, SS.lead_frac_x, SS.lead_frac_sfa, SS.lead_frac_std);
all_passed = check('stable net: leading vector lives mostly in the SFA block', ...
    SS.lead_frac_sfa > SS.lead_frac_x && SS.lead_frac_sfa > SS.lead_frac_std && SS.lead_frac_sfa > 0.5) && all_passed;
all_passed = check('stable net: D_KY = 0 resolved, no positives, K_used 15', ...
    SS.D_KY == 0 && SS.D_KY_resolved == 1 && SS.n_positive == 0 && SS.K_used == 15) && all_passed;

fprintf('\n');
if all_passed
    fprintf('=== ALL lya_K_auto TESTS PASSED ===\n');
else
    fprintf('=== SOME lya_K_auto TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function m = run_preset(P, cond, base, varargin) %#ok<INUSD>  used inside evalc
m = [];
evalc('m = build_from_preset(P, cond, base{:}, varargin{:});');
evalc('m.run();');
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
