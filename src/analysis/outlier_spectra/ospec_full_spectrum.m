function f = ospec_full_spectrum(S, params, ws, k, x_thr)
%OSPEC_FULL_SPECTRUM The rightmost eigenvalues of the FULL Jacobian at one state, which of them are network modes, and the kinetic edge.
%
%   f = ospec_full_spectrum(S, params, ws, k, x_thr)
%
% WHY. The frozen-gain matrix M (ospec_effective_matrix) holds the
% adaptation and depression variables still. The full Jacobian lets them
% respond, which matters exactly near the abscissa: a network mode at a
% rate of order 1/s is as slow as the fastest adaptation variables. This
% returns the k rightmost eigenvalues of the full Jacobian (eigs,
% 'largestreal', on the sparse matrix; a dense eig when the state is small,
% as in the condition without adaptation, where the full Jacobian IS the
% dendritic block).
%
% WHICH OF THEM ARE NETWORK MODES (the rule of this stage, FD-040). The
% full Jacobian of an adapting network has, beside its network modes, one
% slow eigenvalue per adaptation and depression variable: a single
% neuron's variable relaxing at its own rate, with an eigenvector that
% sits in that neuron's slow rows. In a stable network those are the
% rightmost eigenvalues (the slowest relax at about -0.05 1/s), so "the
% rightmost eigenvalue" is then a single-neuron rate that says nothing
% about the network. The two kinds are told apart by the EIGENVECTOR:
%   network mode    xfrac >= x_thr (default 0.1): at least that fraction
%                   of the eigenvector's norm^2 is in the dendritic rows.
%                   It is the rule ospec_dense_spectrum uses for the whole
%                   spectrum.
% In the exploration's runs the split is wide among the rightmost
% eigenvalues (four runs, 5600 eigenvalues, 2026-11-19): slow modes had
% xfrac below 0.03 and spread over fewer than 3 neurons (participation
% ratio pr), network modes had xfrac above 0.47 and pr above 100, nothing
% between. Deeper in the band of slow rates the two kinds mix and the rule
% is a convention (ospec_dense_spectrum has the numbers).
%
% The eigenvalue alone does not separate them. With the recurrent weights
% removed the spectrum is -1/tau_d and the eigenvalues of the slow block
% J(slow, slow), all negative; the rightmost of those is the KINETIC EDGE
% (kin_edge). The exploration called an eigenvalue to the right of it a
% network eigenvalue. But the coupling moves the slow eigenvalues a
% little, to either side: its runs hold single-neuron modes at -0.0446 1/s
% with the kinetic edge at -0.0489. kin_edge is returned for reference.
% What does follow from it: kin_edge < 0, so AN EIGENVALUE WITH POSITIVE
% REAL PART IS A NETWORK MODE under either rule.
%
% WHEN NONE OF THE k IS A NETWORK MODE (a stable network), the rightmost
% network mode is to the LEFT of all k: net is NaN and bound, the smallest
% real part returned, is an upper bound on its real part. Finding it then
% takes the whole spectrum (ospec_dense_spectrum, at the calibration
% states): in the exploration's E-to-E x 1.5 run with multiple timescales
% it was eigenvalue number 679 from the right.
%
% Fields
%   lam     k x 1 complex, by real part (NaN where eigs did not converge)
%   xfrac   k x 1, norm^2 fraction of the eigenvector in the x rows
%   pr      k x 1, participation ratio over neurons (1 to n; all of a
%           neuron's variables together)
%   coh     k x 1, fraction of the x part's norm^2 in the population-uniform span
%   n_converged, seconds
%   abscissa   the rightmost real part, of any kind (1/s)
%   net        the rightmost NETWORK mode among the k (complex, Im >= 0;
%              NaN when none of them is one)
%   net_found  whether one was found
%   bound      smallest real part among the converged eigenvalues: when
%              net_found is false and n_converged == k, every network mode
%              has real part below it
%   kin_edge   rightmost real part of the slow block (-Inf with no slow
%              variables), by ospec_kinetic_edge
%
% Ported from scripts/explorations/outliers/outlier_full_spectrum.m
% (FD-040). Changed: the network-mode fields are new; the kinetic edge is
% computed at every call (the slow block is block diagonal by neuron, so
% it costs a few ms; the exploration took a dense eig of the whole block
% and asked for it at two states per run).
%
% See also: ospec_kinetic_edge, ospec_dense_spectrum, ospec_effective_matrix, ospec_run

if nargin < 4 || isempty(k); k = 8; end
if nargin < 5 || isempty(x_thr); x_thr = 0.1; end

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
    restore = onCleanup(@() warning(w));   % restores on exit
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

ok = isfinite(lam);
is_net = ok & xfrac >= x_thr;
net = complex(NaN, NaN);
q = find(is_net, 1, 'first');                  % lam is sorted by real part
if ~isempty(q)
    net = complex(real(lam(q)), abs(imag(lam(q))));
end
bound = NaN;
if any(ok); bound = min(real(lam(ok))); end

f = struct('lam', lam, 'xfrac', xfrac, 'pr', pr, 'coh', coh, ...
    'n_converged', sum(ok), 'abscissa', max(real(lam(ok))), ...
    'net', net, 'net_found', ~isempty(q), 'bound', bound, ...
    'kin_edge', ospec_kinetic_edge(J, ws), 'seconds', toc(t0));
if isempty(f.abscissa); f.abscissa = NaN; end
end
