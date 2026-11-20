function f = outlier_full_spectrum(S, params, ws, k, want_kinetic)
%OUTLIER_FULL_SPECTRUM The rightmost eigenvalues of the FULL Jacobian at one state, with what each eigenvector is made of.
%
%   f = outlier_full_spectrum(S, params, ws, k, want_kinetic)
%
% WHY. The frozen-gain matrix M (outlier_effective_matrix) holds the
% adaptation and depression variables still. The full Jacobian lets them
% respond, which matters exactly near the abscissa: a network mode at a
% rate of order 1/s is as slow as the fastest adaptation variables. This
% returns the k rightmost eigenvalues of the full Jacobian (eigs,
% 'largestreal', on the sparse matrix; a dense eig when the state is small,
% as in the condition without adaptation, where the full Jacobian IS the
% dendritic block).
%
% WHICH OF THEM ARE NETWORK MODES. With the recurrent weights removed the
% Jacobian is block triangular and its spectrum is -1/tau_d (once per
% neuron) and the eigenvalues of the slow block J(slow, slow): the bands of
% adaptation and depression rates, all negative. Call the rightmost of
% those the kinetic edge. An eigenvalue of the full Jacobian to the RIGHT
% of the kinetic edge cannot be made by the slow variables alone, so it is
% a product of the recurrent coupling: that is the definition of a network
% eigenvalue used here. The eigenvector content is returned beside it to
% describe the mode, not to define it: its norm^2 fraction in x, and the
% number of neurons it spreads over (participation ratio of the per-neuron
% norm^2, all of a neuron's variables together; a single-neuron slow or
% dead-end mode has about 1).
%
% The kinetic edge costs a dense eig of the slow block (about 3500 x 3500
% in the multiple-timescale condition), so it is computed only when
% want_kinetic is true, at a few states per run.
%
% Fields
%   lam     k x 1 complex, by real part (NaN where eigs did not converge)
%   xfrac   k x 1, norm^2 fraction of the eigenvector in the x rows
%   pr      k x 1, participation ratio over neurons (1 to n)
%   coh     k x 1, fraction of the x part's norm^2 in the population-uniform span
%   n_converged, seconds
%   kin_edge   rightmost real part of the slow block (-Inf with no slow
%              variables; NaN when not asked for)
%   kin_band   [5 50 95] percentiles of the slow block's real parts
%
% See also: outlier_effective_matrix, outlier_run

if nargin < 4 || isempty(k); k = 8; end
if nargin < 5; want_kinetic = false; end

t0 = tic;
S = S(:);
J = SRNNCellTypePairs.compute_Jacobian_fast(S, params);
N = size(J, 1);
n = params.n;
k = min(k, N - 2);

if N <= 1000
    [V, D] = eig(full(J));
    lam = diag(D);
    [~, order] = sort(real(lam), 'descend');
    lam = lam(order(1:k));
    V = V(:, order(1:k));
else
    w = warning('off', 'MATLAB:eigs:NotAllEigsConverged');
    restore = onCleanup(@() warning(w));   %#ok<NASGU>  restores on exit
    [V, D] = eigs(J, k, 'largestreal', 'Tolerance', 1e-6, ...
        'MaxIterations', 300, 'SubspaceDimension', min(N, max(60, 6 * k)));
    lam = diag(D);
    [~, order] = sort(real(lam), 'descend', 'MissingPlacement', 'last');
    lam = lam(order);
    V = V(:, order);
end

kk = numel(lam);
xfrac = nan(kk, 1); pr = nan(kk, 1); coh = nan(kk, 1);
for q = 1:kk
    v = V(:, q);
    if ~all(isfinite(v)) || ~isfinite(lam(q)); continue; end
    e = abs(v).^2;
    xfrac(q) = sum(e(ws.is_x)) / sum(e);
    per_neuron = accumarray(ws.neuron_of_state, e, [n 1]);
    pr(q) = sum(per_neuron)^2 / sum(per_neuron.^2);
    vx = v(ws.is_x);
    coh(q) = sum(abs(ws.U' * vx).^2) / max(sum(abs(vx).^2), realmin);
end

kin_edge = NaN; kin_band = nan(1, 3);
if want_kinetic
    slow = find(~ws.is_x);
    if isempty(slow)
        kin_edge = -Inf;
    else
        es = sort(real(eig(full(J(slow, slow)))));
        kin_edge = es(end);
        m = numel(es);
        kin_band = es(max(1, round([0.05 0.5 0.95] * m)))';
    end
end

f = struct('lam', lam, 'xfrac', xfrac, 'pr', pr, 'coh', coh, ...
    'n_converged', sum(isfinite(lam)), 'seconds', toc(t0), ...
    'kin_edge', kin_edge, 'kin_band', kin_band);
end
