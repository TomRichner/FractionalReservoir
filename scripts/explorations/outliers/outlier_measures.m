function m = outlier_measures(R, opts)
%OUTLIER_MEASURES From one outlier run to the measures: how far, how many, how long, one mode or many, against theory.
%
%   m = outlier_measures(R)                 R as saved by outlier_run
%   m = outlier_measures(R, 'overlap_thr', 0.7)
%
% DEFINITIONS (lambda is an eigenvalue of the frozen-gain effective
% connectivity M, in units of W; its rate is (lambda - 1)/tau_d in 1/s)
%
%   bulk edge     edge(t) = kappa R_eff(t). R_eff(t) is the radius
%                 predicted from the gains at that state (variance profile
%                 of M); kappa is the run's median, over the calibration
%                 states, of the 95th percentile of the rightmost real part
%                 of the null draws divided by R_eff. So the edge is the
%                 predicted disk, corrected once per run for finite size.
%                 kappa_abs is the same for the largest modulus, and
%                 agree_abs = median null radius / R_eff says how well the
%                 prediction and the null agree before any percentile.
%   outlier       an eigenvalue of M with Re > edge (right outliers, the
%                 ones that can reach the abscissa). n_out counts them per
%                 state; n_out_abs counts moduli beyond kappa_abs R_eff
%                 (any direction). Silent and saturated neurons put
%                 eigenvalues at exactly 0 (rate -1/tau_d), the center of
%                 the disk, and the slow variables are not in M at all, so
%                 neither can be counted.
%   near the abscissa   the rightmost eigenvalue lambda_1(t), measured two
%                 ways: excess(t) = (Re lambda_1 - edge)/tau_d beyond the
%                 bulk, rate_top(t) = (Re lambda_1 - 1)/tau_d beyond zero.
%   exposure      x(t) = max(0, Re lambda_1 - max(edge, 1))/tau_d: the part
%                 of the rightmost rate that is BOTH outside the bulk and
%                 unstable. Its time average X (1/s) is the proposed
%                 summary measure: zero when the outlier is inside the bulk
%                 or the network is locally stable, and it grows with how
%                 far out, and for how long, the outlier is.
%   dwell         lengths of the runs of consecutive sampled states with
%                 excess > 0, with rate_top > 0, and with x > 0. Resolution
%                 is R.dt: an excursion shorter than that can be missed and
%                 two separated by less are merged. Runs touching either
%                 end of the window are flagged censored.
%   identity      o11(k), the overlap |v_1(k)' v_1(k+1)| of the rightmost
%                 unit eigenvector at consecutive states, and oany(k), its
%                 largest overlap with any of the n_top rightmost at the
%                 next state. Above overlap_thr it is the same mode.
%                 base_overlap is the median overlap of states more than
%                 1 s apart (what "unrelated" looks like here; eigenvectors
%                 of a non-normal matrix are not orthogonal). N_eff is the
%                 participation ratio of the time-averaged projector
%                 sum_k v_1(k) v_1(k)': 1 if one direction leads
%                 throughout, large if many take turns.
%   theory        K1(t), the largest real part of the eigenvalues of
%                 K_eff(t): where the note's section 2 puts the outlier.
%                 Compared with Re lambda_1 where K1 > edge (an outlier is
%                 predicted), and the fraction of states where one is
%                 predicted is given beside the fraction where one is seen.
%   detachment    rho(t) = Re lambda_1 / edge (the table's `ratio`): 1 at
%                 the edge of the bulk. rho_uniform is the same at the
%                 calibration states with every gain set to the rms gain
%                 (the structure's own detachment), and
%                 reined = 1 - (rho - 1)/(rho_uniform - 1) the fraction of
%                 it the gains removed. frac_on_E, frac_on_I (fraction of
%                 neurons off the floor and ceiling) say whether that was
%                 done with the population still responsive or by pinning
%                 it. THIS, read beside X_full, replaced exposure as the
%                 summary after the first runs (see fig_outliers).
%   detachment against the THEORETICAL radius (Tom, 2026-11-20: the bulk is
%                 the unscaled predicted radius; enlarging it by kappa turns
%                 close outliers into non-outliers)
%                 ratio_theory(t) = Re lambda_1(t) / R_eff(t): 1 at the
%                 predicted radius, with no finite-size correction. Because
%                 kappa is one number per run, ratio_theory = kappa ratio
%                 state by state, so its median and percentiles are kappa
%                 times those of `ratio` (to round-off); it is computed
%                 from the per-state values all the same.
%                 rho_uniform_theory is the uniform-gain rightmost real part
%                 over R_eff at each calibration state (median over them).
%                 That one is NOT kappa rho_uniform: rho_uniform divides by
%                 each state's own null percentile, and the median of those
%                 percentiles over R_eff is what kappa is.
%                 reined_theory is `reined` from these two. It is given for
%                 the same runs as `reined` (whether the structure has an
%                 outlier is decided against the null, rho_uniform > 1.05,
%                 since a bulk of finite size reaches past R_eff with no
%                 structure at all: agree_re is how far).
%                 The fields with `ratio`, `rho` and `reined` and no
%                 `_theory` keep their meaning: against the calibrated edge.
%   continuation  at the dense states, the eigenvalue of the FULL Jacobian
%                 whose dendritic pattern matches the rightmost eigenvector
%                 of M (outlier_dense_spectrum): cont_re, cont_im (1/s).
%   X_full        time average of max(0, rightmost real part of the full
%                 Jacobian): the instability that is left when the slow
%                 variables are allowed to respond (1/s).
%   full Jacobian the rightmost eigenvalues with the slow variables free.
%                 An eigenvalue to the right of the kinetic edge is a
%                 network eigenvalue (outlier_full_spectrum). abscissa_full
%                 is the rightmost real part; its dwell beyond zero is at
%                 the coarser R.full.dt.
%
% Output: struct m with the series (1 x n_s) and the scalars; m.row is the
% one-line summary outlier_table prints.
%
% See also: outlier_run, outlier_table, fig_outliers

arguments
    R (1,1) struct
    opts.overlap_thr (1,1) double = 0.7
end

tau_d = R.tau_d;
dt = R.dt;
n_s = numel(R.t);

%% The bulk edge
calR = [R.cal.R_eff];
ok = isfinite(calR) & calR > 1e-6;
if any(ok)
    kappa     = median([R.cal(ok).q95_re]  ./ calR(ok));
    kappa_abs = median([R.cal(ok).q95_abs] ./ calR(ok));
    agree_abs = median([R.cal(ok).med_abs] ./ calR(ok));
    agree_re  = median([R.cal(ok).med_re]  ./ calR(ok));
else
    kappa = 1; kappa_abs = 1; agree_abs = NaN; agree_re = NaN;
end
edge = kappa .* R.R_eff;

%% How far
lam1     = real(R.ev_top(1, :));
rate_top = (lam1 - 1) ./ tau_d;
rate_edge = (edge - 1) ./ tau_d;
excess   = (lam1 - edge) ./ tau_d;
x        = max(0, lam1 - max(edge, 1)) ./ tau_d;
ratio    = nan(1, n_s);
okr = R.R_eff > 1e-6;
ratio(okr) = lam1(okr) ./ edge(okr);
ratio_theory = nan(1, n_s);                            % against the theoretical radius, state by state
ratio_theory(okr) = lam1(okr) ./ R.R_eff(okr);

%% How many
n_out     = sum(R.re_keep  > edge, 1);
n_out_abs = sum(R.abs_keep > kappa_abs .* R.R_eff, 1);

%% How long
[dw_bulk, cz_bulk] = run_lengths(excess > 0, dt);
[dw_zero, cz_zero] = run_lengths(rate_top > 0, dt);
[dw_x,    cz_x]    = run_lengths(x > 0, dt);

%% One mode or many
A = double(squeeze(R.V_top(:, 1, :)));                 % n x n_s
o11 = abs(sum(conj(A(:, 1:end - 1)) .* A(:, 2:end), 1));
oany = o11;
for q = 2:size(R.V_top, 2)
    B = double(squeeze(R.V_top(:, q, :)));
    oany = max(oany, abs(sum(conj(A(:, 1:end - 1)) .* B(:, 2:end), 1)));
end
lag = max(1, round(1 / dt));
if n_s > 2 * lag
    i1 = 1:(n_s - lag);
    base_overlap = median(abs(sum(conj(A(:, i1)) .* A(:, i1 + lag), 1)), 'omitnan');
else
    base_overlap = NaN;
end
same = o11 > opts.overlap_thr;
[life, ~] = run_lengths(same, dt);
turnover = sum(~same) / (n_s * dt);                    % changes of leader per second
det = excess > 0;
if nnz(det) >= 20; use = find(det); else; use = 1:n_s; end
use = use(unique(round(linspace(1, numel(use), min(400, numel(use))))));
G = abs(A(:, use)' * A(:, use)).^2;
N_eff = numel(use)^2 / sum(G(:));
N_eff_on = ternary(nnz(det) >= 20, 'detached states', 'all states');

%% Theory: K_eff
K1 = max(real(R.K_eig), [], 1);
pred_out = K1 > edge;
seen_out = excess > 0;
if nnz(pred_out) >= 5
    K_err_med = median(lam1(pred_out) - K1(pred_out));
    K_err_abs = median(abs(lam1(pred_out) - K1(pred_out)));
else
    K_err_med = NaN; K_err_abs = NaN;
end
K_corr = corr_(K1, lam1);

%% Full Jacobian
F = R.full;
absc_full = max(real(F.lam), [], 1, 'omitnan');
kin_edge = max(R.kin.edge);
if isfinite(kin_edge)
    n_net_full = sum(real(F.lam) > kin_edge, 1);       % capped at k_full
else
    n_net_full = nan(1, numel(F.t));                   % no slow variables: every eigenvalue is a network one
end
[dw_full, cz_full] = run_lengths(absc_full > 0, F.dt);
rate_top_at_full = interp1(R.t, rate_top, F.t, 'nearest', 'extrap');
full_minus_fast = absc_full - rate_top_at_full;

%% Detachment against what the structure alone would give, and collapse
% rho(t) = Re lambda_1 / edge is the detachment (1 at the edge of the
% bulk). rho_uniform is the same ratio at the calibration states with every
% gain replaced by the rms gain: the bulk the same size, the E-I difference
% of the gains removed, so the mean structure is back at its structural
% ratio to the bulk. reined is the fraction of that structural detachment
% the actual gains have removed (1: back at the edge; more than 1: pushed
% inside the bulk; NaN where the structure has no outlier to rein in).
q95 = [R.cal.q95_re];
rho_cal     = [R.cal.obs_re] ./ q95;
rho_uniform = [R.cal.uniform_re] ./ q95;
rho_uniform_med = median(rho_uniform(ok));
if rho_uniform_med > 1.05
    reined = 1 - (median(ratio, 'omitnan') - 1) / (rho_uniform_med - 1);
else
    reined = NaN;
end
% The same against the theoretical radius: each calibration state's own
% R_eff in place of its null percentile. The gate is the calibrated one, so
% reined and reined_theory are given for the same runs.
rho_cal_theory     = nan(size(calR));
rho_uniform_theory = nan(size(calR));
rho_cal_theory(ok)     = [R.cal(ok).obs_re]     ./ calR(ok);
rho_uniform_theory(ok) = [R.cal(ok).uniform_re] ./ calR(ok);
rho_uniform_theory_med = median(rho_uniform_theory(ok));
if rho_uniform_med > 1.05
    reined_theory = 1 - (median(ratio_theory, 'omitnan') - 1) / (rho_uniform_theory_med - 1);
else
    reined_theory = NaN;
end
x_full = mean(max(0, absc_full), 'omitnan');           % instability exposure of the FULL Jacobian (1/s)

%% The dense pass: what became of the rightmost frozen-gain mode
dn = struct('cont', NaN, 'cont_ov', NaN, 'n_pos', NaN, 'abscissa', NaN, 'rate1_M', NaN, 't', []);
if isfield(R, 'dense') && ~isempty(R.dense)
    dn = struct('cont', [R.dense.cont], 'cont_ov', [R.dense.cont_ov], 'n_pos', [R.dense.n_pos], ...
        'abscissa', [R.dense.abscissa], 'rate1_M', [R.dense.rate1_M], 't', R.cal_t);
end

%% Collect
m = struct();
m.network = R.network; m.condition = R.condition; m.seed = R.seed;
m.t = R.t; m.dt = dt; m.tau_d = tau_d;
m.kappa = kappa; m.kappa_abs = kappa_abs; m.agree_abs = agree_abs; m.agree_re = agree_re;
m.edge = edge; m.lam1 = lam1; m.rate_top = rate_top; m.rate_edge = rate_edge;
m.excess = excess; m.x = x; m.ratio = ratio; m.ratio_theory = ratio_theory;
m.n_out = n_out; m.n_out_abs = n_out_abs;
m.dwell = struct('bulk', dw_bulk, 'bulk_censored', cz_bulk, 'zero', dw_zero, ...
    'zero_censored', cz_zero, 'x', dw_x, 'x_censored', cz_x, ...
    'full_zero', dw_full, 'full_zero_censored', cz_full);
m.o11 = o11; m.oany = oany; m.base_overlap = base_overlap;
m.identity_life = life; m.turnover = turnover; m.N_eff = N_eff; m.N_eff_on = N_eff_on;
m.coh1 = R.coh(1, :);
m.K1 = K1; m.pred_out = pred_out; m.K_err_med = K_err_med; m.K_err_abs = K_err_abs; m.K_corr = K_corr;
m.full = struct('t', F.t, 'abscissa', absc_full, 'n_net', n_net_full, 'kin_edge', kin_edge, ...
    'minus_fast', full_minus_fast, 'xfrac1', F.xfrac(1, :), 'pr1', F.pr(1, :), ...
    'coh1', F.coh(1, :), 'n_converged', F.n_converged);
m.shuffle_pct = [R.cal.obs_pct_in_shuffle];
m.cal = struct('obs_re', [R.cal.obs_re], 'typemean_re', [R.cal.typemean_re], ...
    'uniform_re', [R.cal.uniform_re], 'shuffle_med', arrayfun(@(c) median(c.shuf_re), R.cal), ...
    'R_eff', calR, 'q95_re', q95, 'rho', rho_cal, 'rho_uniform', rho_uniform, ...
    'rho_theory', rho_cal_theory, 'rho_uniform_theory', rho_uniform_theory);
m.dense = dn;
m.frac_on = R.frac_on;
m.structural_ratio = max(real(R.structural.K_eig_realized)) / R.structural.R_realized;

lle = NaN;
if isfield(R.lya, 'LLE'); lle = R.lya.LLE; end
E = 1; I = min(2, size(R.dbar, 1));
m.row = struct( ...
    'network', R.network, 'condition', R.condition, 'seed', R.seed(1), ...
    'lambda_1', lle, ...
    'kappa', kappa, 'agree_abs', agree_abs, ...
    'R_eff_med', median(R.R_eff), ...
    'rate_edge_med', median(rate_edge), ...
    'rate_top_med', median(rate_top), 'rate_top_p95', pct(rate_top, 0.95), 'rate_top_max', max(rate_top), ...
    'excess_med', median(excess), 'excess_p95', pct(excess, 0.95), ...
    'ratio_med', median(ratio, 'omitnan'), ...
    'frac_detached', mean(excess > 0), 'frac_unstable', mean(rate_top > 0), 'frac_both', mean(x > 0), ...
    'X', mean(x), ...
    'n_out_mean', mean(n_out), 'n_out_max', max(n_out), 'n_out_abs_mean', mean(n_out_abs), ...
    'n_unstable_mean', mean(R.n_unstable), ...
    'dwell_bulk_med', median(dw_bulk), 'dwell_bulk_max', max([dw_bulk 0]), 'n_exc_bulk', numel(dw_bulk), ...
    'dwell_zero_med', median(dw_zero), 'dwell_zero_max', max([dw_zero 0]), 'n_exc_zero', numel(dw_zero), ...
    'dwell_x_med', median(dw_x), 'dwell_x_max', max([dw_x 0]), 'n_exc_x', numel(dw_x), ...
    'o11_med', median(o11, 'omitnan'), 'frac_same', mean(same), 'base_overlap', base_overlap, ...
    'identity_life_med', median(life), 'turnover_per_s', turnover, 'N_eff', N_eff, ...
    'coh1_med', median(R.coh(1, :), 'omitnan'), ...
    'K1_med', median(K1), 'frac_pred_out', mean(pred_out), 'frac_seen_out', mean(seen_out), ...
    'K_err_med', K_err_med, 'K_corr', K_corr, ...
    'dbar_E', median(squeeze(R.dbar(E, E, :))), 'dbar_I', median(squeeze(R.dbar(E, I, :))), ...
    'rate_E', median(R.rate(E, :)), 'rate_I', median(R.rate(I, :)), ...
    'full_abscissa_med', median(absc_full, 'omitnan'), 'full_abscissa_max', max(absc_full), ...
    'full_frac_unstable', mean(absc_full > 0), 'full_dwell_zero_med', median(dw_full), ...
    'full_n_net_mean', mean(n_net_full, 'omitnan'), 'kin_edge', kin_edge, ...
    'full_minus_fast_med', median(full_minus_fast, 'omitnan'), ...
    'full_states', numel(F.t), ...
    'shuffle_pct_med', median(m.shuffle_pct), ...
    'uniform_minus_obs', median(m.cal.uniform_re - m.cal.obs_re), ...
    'ratio_p05', pct(ratio, 0.05), 'ratio_p95', pct(ratio, 0.95), ...
    'rho_uniform', rho_uniform_med, 'reined', reined, ...
    'structural_ratio', m.structural_ratio, ...
    'agree_re', agree_re, ...
    'frac_on_E', median(R.frac_on(E, :)), 'frac_on_I', median(R.frac_on(I, :)), ...
    'gain_ratio_EI', median(squeeze(R.dbar(E, E, :)) ./ squeeze(R.dbar(E, I, :)), 'omitnan'), ...
    'X_full', x_full, 'full_abscissa_p05', pct(absc_full, 0.05), 'full_abscissa_p95', pct(absc_full, 0.95), ...
    'full_dwell_zero_max', max([dw_full 0]), 'full_n_exc_zero', numel(dw_full), ...
    'cont_re_med', median(real(dn.cont), 'omitnan'), 'cont_im_med', median(imag(dn.cont), 'omitnan'), ...
    'cont_ov_med', median(dn.cont_ov, 'omitnan'), 'dense_n_pos_med', median(dn.n_pos, 'omitnan'), ...
    'dense_frac_unstable', mean(dn.n_pos > 0), ...
    'ratio_theory_med', median(ratio_theory, 'omitnan'), ...
    'ratio_theory_p05', pct(ratio_theory, 0.05), 'ratio_theory_p95', pct(ratio_theory, 0.95), ...
    'rho_uniform_theory', rho_uniform_theory_med, 'reined_theory', reined_theory);
end

%% ------------------------------------------------------------------------
function [len, censored] = run_lengths(mask, dt)
% Durations (s) of the runs of true in a logical row, and whether each
% touches an end of the series.
mask = logical(mask(:)');
d = diff([false mask false]);
on = find(d == 1); off = find(d == -1) - 1;
len = (off - on + 1) .* dt;
censored = (on == 1) | (off == numel(mask));
end

function q = pct(x, p)
x = sort(x(isfinite(x)));
if isempty(x); q = NaN; return; end
q = x(min(numel(x), max(1, round(p * numel(x)))));
end

function r = corr_(a, b)
ok = isfinite(a) & isfinite(b);
a = a(ok) - mean(a(ok)); b = b(ok) - mean(b(ok));
den = sqrt(sum(a.^2) * sum(b.^2));
if den == 0; r = NaN; else; r = sum(a .* b) / den; end
end

function v = ternary(c, a, b)
if c; v = a; else; v = b; end
end
