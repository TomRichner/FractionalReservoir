function [row, mm] = ostruct_mode_measures(R, m)
%OSTRUCT_MODE_MEASURES The measures of the INTENDED mode of one run of the outlier_structures stage, against the theoretical bulk radius.
%
%   row = ostruct_mode_measures(R)            R from ospec_load_run (a run of the stage)
%   [row, mm] = ostruct_mode_measures(R, m)   m = ospec_measures(R), when already computed
%
% Ported from scripts/explorations/outliers/outlier_mode_measures.m
% (FD-043). THE ONE CHANGE OF DEFINITION: the bulk is the UNSCALED
% THEORETICAL radius R_eff(t) of each state (Tom, 2026-11-20), not the
% null-calibrated edge kappa R_eff(t) the exploration divided by. So
% every detachment here is the exploration's times kappa (about 1.1 at
% 500 neurons), and a mode exactly at the theoretical radius has
% detachment 1. The rest is the exploration's: the eigenvalue of the
% frozen-gain effective connectivity M that carries the structural mode's
% pattern (ospec_snapshot), followed every R.opts.dt_mode seconds, and its
% continuation into the full Jacobian at the calibration states
% (ospec_dense_spectrum). Eigenvalues of M are in units of W; a rate is
% (lambda - 1)/tau_d in 1/s; a rotation rate is Im lambda / tau_d in rad/s.
%
% THE THREE PICTURES of the mode
%   structural  lam0 (the eigenvalue of W) times a uniform gain, the rms
%               d_rms of the observed gains at each state (the bulk the
%               same size as observed): what the weights alone give. This
%               is the uniform-gain diamond of the figures.
%   frozen      the eigenvalue of M that carries the mode, slow variables held
%   full        its continuation in the full Jacobian, slow variables free
%               (medians over the calibration states)
%
% FIELDS OF row (medians over the followed states unless said)
%   network, condition, seed1, seed2, target
%   lam0_re_R, lam0_im_R, lam0_abs_R   the structural eigenvalue over the
%               structural bulk radius; coh0, pr0, cond0, p_E (share of the
%               participation on E neurons), struct_is_rightmost
%   lambda_1    leading Lyapunov exponent over the window (1/s)
%   mode_rho    DETACHMENT: Re lambda_mode / R_eff(t); mode_rho_p05,
%               mode_rho_p95 its 5th and 95th percentiles over the states
%   mode_rho_uniform   the same in the structural picture,
%               Re(lam0) d_rms(t) / R_eff(t): the uniform-gain value
%   mode_reined the fraction of the structural detachment the gains
%               removed, 1 - (mode_rho - 1)/(mode_rho_uniform - 1): 1 back
%               at the radius, more inside the bulk; NaN when
%               mode_rho_uniform <= 1.05 (no detachment to remove)
%   mode_rho_abs, mode_rho_abs_uniform, mode_reined_abs   the same by modulus
%   mode_ov     overlap of the followed eigenvector with the structural
%               one; mode_frac_rightmost: fraction of followed states where
%               the mode is the rightmost eigenvalue of M
%   rate_struct, rate_frozen, rate_full   real part of the mode as a rate
%               in the three pictures (1/s)
%   rot_struct, rot_frozen, rot_full   its ROTATION RATE (rad/s) in the
%               three pictures
%   slow_shift  rate_full minus the frozen rate at the same states
%   cont_ov     overlap of the continuation with the frozen mode
%   gain_ratio  the PARTICIPATION-WEIGHTED GAIN over the mean gain (below
%               1: the carriers lost more gain than the average neuron)
%   mean_over_rms   mean gain over rms gain (what the spread of the gains
%               alone takes from a mode spread evenly over the neurons)
%   shift_obs_re, shift_obs_im, shift_pred_re, shift_pred_im   H1: the
%               relative shift lam / (d_rms lam0) - 1, observed and
%               PREDICTED TO FIRST ORDER from the carriers' gains,
%               l0' M r0 / (l0' r0) (ospec_snapshot)
%   on_top, on_rest   RESPONSIVE FRACTION (nonzero gain) of the mode's top
%               10 % of carriers and of the rest; gain_top_rest the ratio
%               of their mean gains
%   frac_on_E, frac_on_I   responsive fraction of the E and of the I neurons
%   d_rms       rms gain
%   right_rho, right_rho_uniform, right_reined   the same three for the
%               RIGHTMOST eigenvalue of M (ospec_measures' _theory fields)
%   gain_ratio_EI, rate_E, rate_I   from ospec_measures
%   full_abscissa, full_frac_unstable   the rightmost real part of the
%               full Jacobian over the calibration states (median), and the
%               fraction of those states with an eigenvalue right of zero
%   n_out103_mean   eigenvalues of M with real part beyond 1.03 R_eff, mean
%               over the states
%   kappa       the null-calibrated factor of the run (to compare with the
%               exploration's numbers)
%
% mm: the series behind the medians, at the followed states: t, lam,
% R_eff, rho, rho_uniform, rho_abs, rho_abs_uniform, ov, is_rightmost,
% lam_pred, shift_obs, shift_pred, gain_ratio, on_top, on_rest, d_rms.
%
% See also: ostruct_load, ospec_measures, ospec_snapshot, ospec_load_run

arguments
    R (1,1) struct
    m = []
end
if isempty(m); m = ospec_measures(R); end
if ~isfield(R, 'mode') || ~isfield(R.mode, 'lam') || isempty(R.mode.lam)
    error('ostruct_mode_measures:NoMode', ...
        'Run %s | %s did not follow the intended mode (dt_mode was 0).', R.network, R.condition);
end
tau = R.tau_d;
ref = R.ref;
md = R.mode;

[tf, at] = ismember(md.t, R.t);
if ~all(tf)
    error('ostruct_mode_measures:Times', 'The followed states are not a subset of the sampled states.');
end
Reff = R.R_eff(at);
ok = Reff > 1e-6;
lam = md.lam;
d_rms = md.d_rms;
rho = nan(size(lam)); rho_u = rho; rho_abs = rho; rho_abs_u = rho;
rho(ok)       = real(lam(ok)) ./ Reff(ok);
rho_u(ok)     = real(ref.lam0) .* d_rms(ok) ./ Reff(ok);
rho_abs(ok)   = abs(lam(ok)) ./ Reff(ok);
rho_abs_u(ok) = abs(ref.lam0) .* d_rms(ok) ./ Reff(ok);
reined     = reined_fraction(median(rho, 'omitnan'), median(rho_u, 'omitnan'));
reined_abs = reined_fraction(median(rho_abs, 'omitnan'), median(rho_abs_u, 'omitnan'));

% The calibration states: the continuation in the full Jacobian.
cont = complex(NaN, NaN); cont_ov = NaN; slow_shift = NaN; full_abscissa = NaN; full_frac_unstable = NaN;
if isfield(R, 'dense') && ~isempty(R.dense)
    D = R.dense;
    mc = [D.mode_cont];
    cont = complex(median(real(mc)), median(imag(mc)));
    cont_ov = median([D.mode_cont_ov]);
    slow_shift = median(real(mc) - (real([D.mode_lam_M]) - 1) ./ tau);
    full_abscissa = median([D.abscissa]);
    full_frac_unstable = mean([D.n_pos] > 0);
end

r = m.row;
row = struct( ...
    'network', R.network, 'condition', R.condition, 'seed1', R.seed(1), 'seed2', R.seed(2), ...
    'target', ref.target, ...
    'lam0_re_R', real(ref.lam0) / ref.R_struct, 'lam0_im_R', imag(ref.lam0) / ref.R_struct, ...
    'lam0_abs_R', ref.bulk_ratio, 'coh0', ref.coh, 'pr0', ref.pr, 'cond0', ref.cond, ...
    'p_E', ref.p_type(1), 'struct_is_rightmost', ref.is_rightmost, ...
    'lambda_1', r.lambda_1, ...
    'mode_rho', median(rho, 'omitnan'), 'mode_rho_p05', pct(rho, 0.05), 'mode_rho_p95', pct(rho, 0.95), ...
    'mode_rho_uniform', median(rho_u, 'omitnan'), 'mode_reined', reined, ...
    'mode_rho_abs', median(rho_abs, 'omitnan'), 'mode_rho_abs_uniform', median(rho_abs_u, 'omitnan'), ...
    'mode_reined_abs', reined_abs, ...
    'mode_ov', median(md.ov), 'mode_frac_rightmost', mean(md.is_rightmost), ...
    'rate_struct', (real(ref.lam0) * median(d_rms) - 1) / tau, ...
    'rate_frozen', median(real(lam) - 1) / tau, 'rate_full', real(cont), ...
    'rot_struct', imag(ref.lam0) * median(d_rms) / tau, ...
    'rot_frozen', median(imag(lam)) / tau, 'rot_full', imag(cont), ...
    'slow_shift', slow_shift, 'cont_ov', cont_ov, ...
    'gain_ratio', median(md.gain_ratio), 'mean_over_rms', median(md.d_mean ./ d_rms), ...
    'shift_obs_re', median(real(md.shift_obs)), 'shift_obs_im', median(imag(md.shift_obs)), ...
    'shift_pred_re', median(real(md.shift_pred)), 'shift_pred_im', median(imag(md.shift_pred)), ...
    'on_top', median(md.on_top), 'on_rest', median(md.on_rest), ...
    'gain_top_rest', median(md.gain_top ./ max(md.gain_rest, realmin)), ...
    'frac_on_E', r.frac_on_E, 'frac_on_I', r.frac_on_I, 'd_rms', median(d_rms), ...
    'right_rho', r.ratio_theory_med, 'right_rho_uniform', r.rho_uniform_theory, 'right_reined', r.reined_theory, ...
    'gain_ratio_EI', r.gain_ratio_EI, 'rate_E', r.rate_E, 'rate_I', r.rate_I, ...
    'full_abscissa', full_abscissa, 'full_frac_unstable', full_frac_unstable, ...
    'n_out103_mean', r.n_out_theory103_mean, 'kappa', m.kappa);

mm = struct('t', md.t, 'lam', lam, 'R_eff', Reff, 'rho', rho, 'rho_uniform', rho_u, ...
    'rho_abs', rho_abs, 'rho_abs_uniform', rho_abs_u, 'ov', md.ov, 'is_rightmost', md.is_rightmost, ...
    'lam_pred', md.lam_pred, 'shift_obs', md.shift_obs, 'shift_pred', md.shift_pred, ...
    'gain_ratio', md.gain_ratio, 'on_top', md.on_top, 'on_rest', md.on_rest, 'd_rms', d_rms);
end

%% ------------------------------------------------------------------------
function f = reined_fraction(rho, rho_uniform)
% The fraction of the structural detachment removed; NaN when the structure
% has none worth the name.
if rho_uniform > 1.05
    f = 1 - (rho - 1) / (rho_uniform - 1);
else
    f = NaN;
end
end

function q = pct(x, p)
x = sort(x(isfinite(x)));
if isempty(x); q = NaN; return; end
q = x(min(numel(x), max(1, round(p * numel(x)))));
end
