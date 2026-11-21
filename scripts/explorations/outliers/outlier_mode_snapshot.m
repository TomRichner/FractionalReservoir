function s = outlier_mode_snapshot(S, params, ws, ref)
%OUTLIER_MODE_SNAPSHOT The intended mode at one state: where its eigenvalue is in the frozen-gain matrix, and what the gains of its carriers predict.
%
%   s = outlier_mode_snapshot(S, params, ws, ref)     ref = outlier_mode_reference(...)
%
% THE MODE IS FOLLOWED BY ITS PATTERN, not by being the rightmost: among
% all eigenvectors v_k of the frozen-gain effective connectivity M at this
% state (outlier_effective_matrix), the one with the largest overlap
% |r0' v_k| with the structural right eigenvector. With uniform gains that
% is the structural mode exactly. The overlap is returned: when it is low
% the mode has not survived as one eigenvector (its carriers are silent or
% saturated, or it has mixed into the bulk) and the eigenvalue reported is
% only the best match.
%
% H1, THE FIRST-ORDER PREDICTION. With uniform gain d_ref the eigenvalue
% is d_ref lam0. Actual gains change the columns of M, and to first order
%   lam_pred = l0' M r0 / (l0' r0)  =  lam0 sum_j p_j d_j     (one gain per neuron)
% (the first form is used: it is right also when a neuron's gain differs
% by postsynaptic type). d_ref is the rms of all gains, the uniform
% control of outlier_null: the bulk the same size. shift_obs and
% shift_pred are lam / (d_ref lam0) - 1 for the followed and the predicted
% eigenvalue: 0 if the mode moved as the bulk did, negative real part if
% it lost more than the bulk (reined in), and the imaginary part is the
% change of its angle.
%
% Fields (eigenvalues in units of W)
%   lam          the followed eigenvalue (Im >= 0), ov its overlap with r0
%   is_rightmost whether it is the rightmost eigenvalue of M
%   lam_right    the rightmost eigenvalue of M (Im >= 0)
%   lam_pred     the first-order prediction
%   d_rms, d_mean   rms and mean of all gains (C x n)
%   gain_w       gain averaged over neurons with weights |p_j|
%   gain_ratio   gain_w / mean per-neuron gain
%   shift_obs, shift_pred   complex, see above
%   on_top, on_rest   fraction of responsive neurons (nonzero gain) among
%                the mode's top carriers (ref.top) and among the rest
%   gain_top, gain_rest   mean gain of the two groups
%
% See also: outlier_mode_reference, outlier_effective_matrix, outlier_run

[Mf, d] = outlier_effective_matrix(S, params, ws);
[V, D] = eig(Mf);
ev = diag(D);
ov_all = abs(ref.r0' * V);                    % eig returns unit columns
[ov, k] = max(ov_all);
[~, kr] = max(real(ev));
lam = complex(real(ev(k)), abs(imag(ev(k))));
lam_right = complex(real(ev(kr)), abs(imag(ev(kr))));

% One gain per neuron: the average over its routes, weighted by the number
% of postsynaptic neurons of each type.
n_per = cellfun(@numel, ws.ti);
dj = (n_per(:)' * d) ./ sum(n_per);            % 1 x n
ap = abs(ref.p(:))';
d_rms = sqrt(mean(d(:).^2));
lam_pred = (ref.l0' * Mf * ref.r0) / (ref.l0' * ref.r0);
on = any(d > 0, 1);
rest = true(1, numel(dj)); rest(ref.top) = false;

s = struct('lam', lam, 'ov', ov, 'is_rightmost', abs(lam - lam_right) < 1e-9 * max(1, abs(lam_right)), ...
    'lam_right', lam_right, 'lam_pred', lam_pred, ...
    'd_rms', d_rms, 'd_mean', mean(d(:)), ...
    'gain_w', sum(ap .* dj) / sum(ap), 'gain_ratio', (sum(ap .* dj) / sum(ap)) / mean(dj), ...
    'shift_obs', lam / (d_rms * ref.lam0) - 1, 'shift_pred', lam_pred / (d_rms * ref.lam0) - 1, ...
    'on_top', mean(on(ref.top)), 'on_rest', mean(on(rest)), ...
    'gain_top', mean(dj(ref.top)), 'gain_rest', mean(dj(rest)));
end
