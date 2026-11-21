function [row, mm] = outlier_mode_measures(R, m)
%OUTLIER_MODE_MEASURES The measures of the INTENDED mode of one outlier run, beside those of the rightmost eigenvalue.
%
%   row = outlier_mode_measures(R)             R as saved by outlier_run (with R.mode, R.ref)
%   [row, mm] = outlier_mode_measures(R, m)    m = outlier_measures(R), when already computed
%
% outlier_measures follows the RIGHTMOST eigenvalue of the frozen-gain
% effective connectivity M. This follows the mode the network was built to
% have (outlier_mode_reference), by its eigenvector, wherever it goes, and
% says when the two differ. Eigenvalues of M are in units of W; a rate is
% (lambda - 1)/tau_d in 1/s.
%
% THE THREE PICTURES of the mode
%   structural  the eigenvalue lam0 of W times a uniform gain, the rms of
%               the observed gains (the bulk the same size as observed):
%               what the weights alone give
%   frozen      the eigenvalue of M that carries the mode's pattern, slow
%               variables held (outlier_mode_snapshot)
%   full        the eigenvalue of the full Jacobian that carries the frozen
%               mode's dendritic pattern, slow variables free
%               (outlier_dense_spectrum, at the calibration states)
%
% Fields of row (medians over the sampled states unless said)
%   lam0_re_R, lam0_im_R, lam0_abs_R   the structural eigenvalue over the
%               structural bulk radius; coh0, pr0, cond0, p_E as in
%               outlier_mode_reference (p_E: share of the participation on
%               the excitatory neurons); struct_is_rightmost
%   mode_rho    detachment of the mode by real part: Re lambda / edge, edge
%               the null-calibrated bulk edge kappa R_eff(t) of
%               outlier_measures. mode_rho_uniform: the same in the
%               structural picture. mode_reined: the fraction of the
%               structural detachment the gains removed,
%               1 - (mode_rho - 1)/(mode_rho_uniform - 1) (1: at the edge,
%               more: inside the bulk; NaN when the structure has no
%               detachment to remove).
%   mode_rho_abs, mode_rho_abs_uniform, mode_reined_abs   the same by
%               modulus, against kappa_abs R_eff (the disk, any direction)
%   mode_ov     overlap of the followed eigenvector with the structural
%               one; mode_frac_rightmost: fraction of states where the
%               followed mode is the rightmost eigenvalue of M
%   rate_struct, rate_frozen, rate_full   real part of the mode in the
%               three pictures (1/s); rot_struct, rot_frozen, rot_full its
%               rotation rate, the imaginary part (rad/s)
%   slow_shift  rate_full minus the frozen rate at the same states: what
%               letting the slow variables respond does to the real part
%   cont_ov     overlap of the continuation with the frozen mode
%   gain_ratio  participation-weighted gain over the mean gain (below 1:
%               the carriers have lost more gain than the average neuron)
%   mean_over_rms   mean gain over rms gain. The bulk radius follows the
%               rms of the gains and a mode spread evenly over the neurons
%               follows their mean, so this is what the SPREAD of the
%               gains alone takes from such a mode relative to the bulk
%               (1 with equal gains; sqrt(f) when a fraction f of the
%               neurons share one gain and the rest have none). To first
%               order 1 + shift_pred = mean_over_rms x gain_ratio for a
%               mode with real, positive participations.
%   shift_obs_re, shift_obs_im, shift_pred_re, shift_pred_im   H1: the
%               relative shift of the eigenvalue from the structural
%               picture, lam / (d_rms lam0) - 1, observed and predicted to
%               first order from the carriers' gains
%   on_top, on_rest   responsive fraction (nonzero gain) among the top 10%
%               of carriers and among the rest; gain_top_rest: mean gain
%               of the first over that of the second
%   and, from outlier_measures: lambda_1, the rightmost eigenvalue's
%   detachment (right_rho, right_rho_uniform, right_reined), frac_on_E,
%   frac_on_I, gain_ratio_EI, rate_E, rate_I, right_rate_frozen,
%   full_abscissa (median rightmost real part of the full Jacobian over the
%   dense states), full_frac_unstable (fraction of dense states with an
%   eigenvalue right of zero), right_cont_re, right_cont_im.
%
% mm holds the series behind the medians (at R.mode_t).
%
% See also: outlier_run, outlier_measures, outlier_mode_snapshot, fig_outlier_structures

if nargin < 2 || isempty(m); m = outlier_measures(R); end
if ~isfield(R, 'mode') || isempty(R.mode)
    error('outlier_mode_measures:NoMode', ...
        'Run %s | %s has no mode pass (made before 2026-11-20): rerun it with outlier_run.', R.network, R.condition);
end
tau = R.tau_d;
ref = R.ref;
md = R.mode;

% The bulk edge at the mode-pass states (a subset of the sampled states).
[tf, at] = ismember(R.mode_t, R.t);
if ~all(tf); error('outlier_mode_measures:Times', 'The mode-pass states are not a subset of the sampled states.'); end
edge = m.edge(at);
edge_abs = m.kappa_abs .* R.R_eff(at);
ok = R.R_eff(at) > 1e-6;

lam = [md.lam];
d_rms = [md.d_rms];
rho = nan(size(lam)); rho_u = rho; rho_abs = rho; rho_abs_u = rho;
rho(ok)       = real(lam(ok)) ./ edge(ok);
rho_u(ok)     = real(ref.lam0) .* d_rms(ok) ./ edge(ok);
rho_abs(ok)   = abs(lam(ok)) ./ edge_abs(ok);
rho_abs_u(ok) = abs(ref.lam0) .* d_rms(ok) ./ edge_abs(ok);
reined     = reined_fraction(median(rho, 'omitnan'), median(rho_u, 'omitnan'));
reined_abs = reined_fraction(median(rho_abs, 'omitnan'), median(rho_abs_u, 'omitnan'));

% The dense states: the continuation in the full Jacobian.
cont = NaN; cont_ov = NaN; slow_shift = NaN; full_abscissa = NaN; full_frac_unstable = NaN;
if isfield(R, 'dense') && ~isempty(R.dense)
    D = R.dense;
    cont = median(real([D.mode_cont])) + 1i * median(imag([D.mode_cont]));
    cont_ov = median([D.mode_cont_ov]);
    slow_shift = median(real([D.mode_cont]) - (real([D.mode_lam_M]) - 1) ./ tau);
    full_abscissa = median([D.abscissa]);
    full_frac_unstable = mean([D.n_pos] > 0);
end

r = m.row;
row = struct( ...
    'network', R.network, 'condition', R.condition, 'seed', R.seed(1), ...
    'target', ref.target, ...
    'lam0_re_R', real(ref.lam0) / ref.R_struct, 'lam0_im_R', imag(ref.lam0) / ref.R_struct, ...
    'lam0_abs_R', ref.bulk_ratio, 'coh0', ref.coh, 'pr0', ref.pr, 'cond0', ref.cond, ...
    'p_E', ref.p_type(1), 'struct_is_rightmost', ref.is_rightmost, ...
    'lambda_1', r.lambda_1, ...
    'mode_rho', median(rho, 'omitnan'), 'mode_rho_uniform', median(rho_u, 'omitnan'), 'mode_reined', reined, ...
    'mode_rho_abs', median(rho_abs, 'omitnan'), 'mode_rho_abs_uniform', median(rho_abs_u, 'omitnan'), ...
    'mode_reined_abs', reined_abs, ...
    'mode_ov', median([md.ov]), 'mode_frac_rightmost', mean([md.is_rightmost]), ...
    'rate_struct', (real(ref.lam0) * median(d_rms) - 1) / tau, ...
    'rate_frozen', median(real(lam) - 1) / tau, 'rate_full', real(cont), ...
    'rot_struct', imag(ref.lam0) * median(d_rms) / tau, ...
    'rot_frozen', median(imag(lam)) / tau, 'rot_full', imag(cont), ...
    'slow_shift', slow_shift, 'cont_ov', cont_ov, ...
    'gain_ratio', median([md.gain_ratio]), ...
    'shift_obs_re', median(real([md.shift_obs])), 'shift_obs_im', median(imag([md.shift_obs])), ...
    'shift_pred_re', median(real([md.shift_pred])), 'shift_pred_im', median(imag([md.shift_pred])), ...
    'on_top', median([md.on_top]), 'on_rest', median([md.on_rest]), ...
    'gain_top_rest', median([md.gain_top] ./ max([md.gain_rest], realmin)), ...
    'd_rms', median(d_rms), 'mean_over_rms', median([md.d_mean] ./ d_rms), ...
    'right_rho', r.ratio_med, 'right_rho_uniform', r.rho_uniform, 'right_reined', r.reined, ...
    'frac_on_E', r.frac_on_E, 'frac_on_I', r.frac_on_I, 'gain_ratio_EI', r.gain_ratio_EI, ...
    'rate_E', r.rate_E, 'rate_I', r.rate_I, 'right_rate_frozen', r.rate_top_med, ...
    'full_abscissa', full_abscissa, 'full_frac_unstable', full_frac_unstable, ...
    'right_cont_re', r.cont_re_med, 'right_cont_im', r.cont_im_med, ...
    'kappa', m.kappa, 'rate_edge', r.rate_edge_med);

mm = struct('t', R.mode_t, 'lam', lam, 'edge', edge, 'edge_abs', edge_abs, ...
    'rho', rho, 'rho_uniform', rho_u, 'rho_abs', rho_abs, 'rho_abs_uniform', rho_abs_u, ...
    'ov', [md.ov], 'is_rightmost', [md.is_rightmost], 'lam_pred', [md.lam_pred], ...
    'shift_obs', [md.shift_obs], 'shift_pred', [md.shift_pred], 'gain_ratio', [md.gain_ratio], ...
    'on_top', [md.on_top], 'on_rest', [md.on_rest], 'd_rms', d_rms);
end

function f = reined_fraction(rho, rho_uniform)
% The fraction of the structural detachment removed; NaN when the structure
% has none worth the name (outlier_measures' threshold).
if rho_uniform > 1.05
    f = 1 - (rho - 1) / (rho_uniform - 1);
else
    f = NaN;
end
end
