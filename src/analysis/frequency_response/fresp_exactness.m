function E = fresp_exactness(S, params, ws, opts)
%FRESP_EXACTNESS Is the elimination of the slow variables exact and local? The reduced problem checked against the dense full Jacobian at one state.
%
%   E = fresp_exactness(S, params, ws)
%   E = fresp_exactness(S, params, ws, opts)
%
% EXACT. All eigenvalues z_k of the dense full Jacobian are computed. z_k
% is an eigenvalue of the full Jacobian exactly when 1 + tau_d z_k is an
% eigenvalue of M(z_k) = W .* d(z_k), d the closed-form gains (fresp_gain).
% For the n_lead rightmost eigenvalues, the one nearest opts.cont (the
% leftover mode the outlier_spectra stage followed) and n_sample drawn at
% random:
%   sigma  smallest singular value of M(z_k) - (1 + tau_d z_k) I, over the
%          largest singular value of M(z_k)  (zero when 1 + tau_d z_k is an
%          eigenvalue of M(z_k))
%   eigd   the distance from 1 + tau_d z_k to the nearest eigenvalue of
%          M(z_k), over |1 + tau_d z_k|
%   dslow  the distance from z_k to the nearest eigenvalue of the slow
%          block J_ss. The elimination divides by z - J_ss: an eigenvalue
%          of J that is ALSO one of J_ss (a slow mode of a neuron whose
%          slope is zero) is not a root of the reduced problem and does
%          not need to be. Judge sigma where dslow > 1e-6.
%
% LOCAL. M(z) from the Jacobian's own blocks (fresp_schur; nothing there
% knows the equations) against the closed form, at the frequencies
% opts.exact_f:
%   closed_diff  largest |closed - Schur| gain over the largest |Schur| gain
%   offdiag      how far the Schur M(z) is from "the weights times one gain
%                per presynaptic neuron and route" (Frobenius, relative)
% Both are zero to round-off: a neuron's slow variables are driven by its
% own rate and by nothing else.
%
% A state without slow variables (N = n) has nothing to eliminate: the
% eigenvalues are returned and the checks are empty.
%
%   S        N x 1 state
%   params, ws   model.get_params(), fresp_weight_stats(params)
%   opts     .n_lead (20), .n_sample (200), .exact_f ([0.001 0.03 0.5 3 30]
%            Hz), .cont ([] or the leftover mode, 1/s)
%
% E fields
%   lam        N x 1 complex: all eigenvalues of the full Jacobian, by real
%              part, descending (double)
%   abscissa   the rightmost real part
%   z, kind    the eigenvalues tested and why (1 leading, 2 nearest the
%              leftover mode, 3 random)
%   sigma, eigd, dslow   as above, one per tested eigenvalue
%   exact_f, closed_diff, offdiag   as above, one per frequency
%   seconds
%
% See also: fresp_gain, fresp_schur, fresp_matrix, fresp_run_one

if nargin < 4; opts = struct(); end
n_lead   = get_opt(opts, 'n_lead', 20);
n_sample = get_opt(opts, 'n_sample', 200);
exact_f  = get_opt(opts, 'exact_f', [0.001 0.03 0.5 3 30]);
cont     = get_opt(opts, 'cont', []);

t0 = tic;
S = S(:);
tau_d = params.tau_d;
n = params.n;
J = SRNNCellTypePairs.compute_Jacobian_fast(S, params);
zk = eig(full(J));
[~, order] = sort(real(zk), 'descend');
zk = zk(order);

E = struct('lam', zk, 'abscissa', real(zk(1)), 'z', zeros(0, 1), 'kind', zeros(0, 1), ...
    'sigma', zeros(0, 1), 'eigd', zeros(0, 1), 'dslow', zeros(0, 1), ...
    'exact_f', exact_f, 'closed_diff', nan(1, numel(exact_f)), 'offdiag', nan(1, numel(exact_f)), ...
    'seconds', 0);
slow = find(~ws.is_x);
if isempty(slow)
    E.seconds = toc(t0);
    return;
end
z_slow = eig(full(J(slow, slow)));

n_lead = min(n_lead, numel(zk));
pick = (1:n_lead)';
kind = ones(n_lead, 1);
if ~isempty(cont) && isfinite(cont)
    [~, ic] = min(abs(zk - cont));
    pick(end + 1) = ic; kind(end + 1) = 2;
end
stream = RandStream('twister', 'Seed', 1);
rest = setdiff((1:numel(zk))', pick);
extra = rest(randperm(stream, numel(rest), min(n_sample, numel(rest))));
pick = [pick; extra(:)]; kind = [kind; 3 * ones(numel(extra), 1)];

m = numel(pick);
sigma = zeros(m, 1); eigd = zeros(m, 1); dslow = zeros(m, 1);
I = eye(n);
for q = 1:m
    zq = zk(pick(q));
    Mz = fresp_matrix(fresp_gain(S, params, zq), ws);
    target = 1 + tau_d * zq;
    sv = svd(Mz - target * I);
    sigma(q) = sv(end) / norm(Mz);
    eigd(q) = min(abs(eig(Mz) - target)) / abs(target);
    dslow(q) = min(abs(zq - z_slow));
end

zc = 2i * pi * exact_f;
ok = ws.col_sq > 0;
for q = 1:numel(zc)
    [~, dn, E.offdiag(q)] = fresp_schur(J, params, ws, zc(q));
    dc = fresp_gain(S, params, zc(q));
    E.closed_diff(q) = max(abs(dc(ok) - dn(ok))) / max(abs(dn(ok)));
end

E.z = zk(pick); E.kind = kind; E.sigma = sigma; E.eigd = eigd; E.dslow = dslow;
E.seconds = toc(t0);
end

function v = get_opt(s, name, default)
v = default;
if isfield(s, name); v = s.(name); end
end
