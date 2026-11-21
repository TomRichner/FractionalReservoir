function s = ospec_snapshot(S, params, ws, n_top, n_keep, ref)
%OSPEC_SNAPSHOT The spectrum of the frozen-gain effective connectivity at one state, reduced to what the outlier measures use.
%
%   s = ospec_snapshot(S, params, ws, n_top, n_keep)
%   s = ospec_snapshot(S, params, ws, n_top, n_keep, ref)   also follows the INTENDED mode
%
% One dense eig of the n x n matrix M (ospec_effective_matrix), with its
% eigenvectors. A conjugate pair counts as one mode: the member with
% Im >= 0 represents it.
%
% Ported from scripts/explorations/outliers/outlier_snapshot.m and
% outlier_mode_snapshot.m (FD-040). TWO CHANGES. (1) The eigenvectors of
% the rightmost modes come from the same eig call as the eigenvalues; the
% exploration took the eigenvalues alone and then found each vector by
% inverse iteration from a random start. At n = 500 that cost 0.20 s of
% the 0.26 s a state took, against 0.013 s for asking eig for the vectors
% (measured 2026-11-20, MacBook Air M3, 8 threads). The eigenvalues agree
% with the exploration's to round-off; a vector is the same direction up
% to a phase, which nothing downstream reads (every use is an absolute
% overlap), and no random number is drawn. (2) With ref, the intended mode is followed from the
% same decomposition (the exploration's separate mode pass did a second
% eig of the same matrix).
%
% Nothing here depends on where the bulk edge is drawn: the predicted
% radius R_eff is returned, and the sorted real parts and moduli are kept
% (the n_keep largest of each), so outliers can be counted afterwards
% against any edge.
%
% Fields (units of W unless said; a rate in 1/s is (value - 1)/tau_d)
%   ev_top   n_top x 1 complex: the rightmost modes (Im >= 0), by real part
%   V_top    n x n_top single complex: their unit right eigenvectors
%   res      n_top x 1: residual norm(M v - lambda v) of each
%   coh      n_top x 1: fraction of each eigenvector's norm^2 in the span of
%            the population-uniform patterns (a random vector has C/n)
%   re_keep  n_keep x 1: the largest real parts of ALL eigenvalues, sorted
%   abs_keep n_keep x 1: the largest moduli, sorted
%   n_unstable   eigenvalues with Re > 1 (frozen-gain rate > 0)
%   R_eff, K_eig, K_eff   the predictions (ospec_effective_matrix)
%   dbar, drms, rate, frac_on   the population descriptors (same)
%   mode     [] without ref; with ref (ospec_mode_reference) a struct, the
%            intended mode at this state, FOLLOWED BY ITS PATTERN, not by
%            being the rightmost: among all eigenvectors v_k of M the one
%            with the largest overlap |r0' v_k| with the structural right
%            eigenvector. When the overlap is low the mode has not survived
%            as one eigenvector and the eigenvalue is only the best match.
%              lam          the followed eigenvalue (Im >= 0), ov its overlap
%              is_rightmost whether it is the rightmost eigenvalue of M
%              lam_right    the rightmost eigenvalue of M (Im >= 0)
%              lam_pred     the first-order prediction l0' M r0 / (l0' r0)
%              d_rms, d_mean   rms and mean of all gains (C x n)
%              gain_w       gain averaged over neurons with weights |p_j|
%              gain_ratio   gain_w / mean per-neuron gain
%              shift_obs, shift_pred   lam / (d_rms lam0) - 1 for the
%                           followed and the predicted eigenvalue: 0 if the
%                           mode moved as the bulk did
%              on_top, on_rest   fraction of responsive neurons among the
%                           mode's top carriers (ref.top) and the rest
%              gain_top, gain_rest   mean gain of the two groups
%
% See also: ospec_effective_matrix, ospec_null, ospec_mode_reference, ospec_run

if nargin < 4 || isempty(n_top);  n_top = 3;   end
if nargin < 5 || isempty(n_keep); n_keep = 30; end
if nargin < 6; ref = []; end

[Mf, d, ~, pred, pop] = ospec_effective_matrix(S, params, ws);
n = size(Mf, 1);
n_keep = min(n_keep, n);

[V, D] = eig(Mf);
ev = diag(D);
re_sorted  = sort(real(ev), 'descend');
abs_sorted = sort(abs(ev), 'descend');

upper = find(imag(ev) >= 0);
[~, order] = sort(real(ev(upper)), 'descend');
upper = upper(order);
n_top = min(n_top, numel(upper));
pick = upper(1:n_top);
ev_top = ev(pick);
V_top = V(:, pick);                          % eig returns unit columns
res = zeros(n_top, 1); coh = zeros(n_top, 1);
for q = 1:n_top
    v = V_top(:, q);
    res(q) = norm(Mf * v - ev_top(q) * v);
    coh(q) = sum(abs(ws.U' * v).^2);
end

mode = [];
if ~isempty(ref)
    [ov, k] = max(abs(ref.r0' * V));
    [~, kr] = max(real(ev));
    lam = complex(real(ev(k)), abs(imag(ev(k))));
    lam_right = complex(real(ev(kr)), abs(imag(ev(kr))));
    % One gain per neuron: the average over its routes, weighted by the
    % number of postsynaptic neurons of each type.
    n_per = cellfun(@numel, ws.ti);
    dj = (n_per(:)' * d) ./ sum(n_per);        % 1 x n
    ap = abs(ref.p(:))';
    d_rms = sqrt(mean(d(:).^2));
    lam_pred = (ref.l0' * Mf * ref.r0) / (ref.l0' * ref.r0);
    on = any(d > 0, 1);
    rest = true(1, numel(dj)); rest(ref.top) = false;
    mode = struct('lam', lam, 'ov', ov, ...
        'is_rightmost', abs(lam - lam_right) < 1e-9 * max(1, abs(lam_right)), ...
        'lam_right', lam_right, 'lam_pred', lam_pred, ...
        'd_rms', d_rms, 'd_mean', mean(d(:)), ...
        'gain_w', sum(ap .* dj) / sum(ap), 'gain_ratio', (sum(ap .* dj) / sum(ap)) / mean(dj), ...
        'shift_obs', lam / (d_rms * ref.lam0) - 1, 'shift_pred', lam_pred / (d_rms * ref.lam0) - 1, ...
        'on_top', mean(on(ref.top)), 'on_rest', mean(on(rest)), ...
        'gain_top', mean(dj(ref.top)), 'gain_rest', mean(dj(rest)));
end

s = struct('ev_top', ev_top, 'V_top', single(V_top), 'res', res, 'coh', coh, ...
    're_keep', re_sorted(1:n_keep), 'abs_keep', abs_sorted(1:n_keep), ...
    'n_unstable', sum(real(ev) > 1), ...
    'R_eff', pred.R_eff, 'K_eig', pred.K_eig, 'K_eff', pred.K_eff, ...
    'dbar', pop.dbar, 'drms', pop.drms, 'rate', pop.rate, 'frac_on', pop.frac_on, ...
    'mode', mode);
end
