function s = outlier_snapshot(S, params, ws, n_top, n_keep)
%OUTLIER_SNAPSHOT The spectrum of the frozen-gain effective connectivity at one state, reduced to what the outlier measures use.
%
%   s = outlier_snapshot(S, params, ws, n_top, n_keep)
%
% One dense eig of the n x n matrix M (outlier_effective_matrix), then the
% right eigenvectors of the n_top rightmost modes by inverse iteration
% (cheaper than asking eig for all n vectors). A conjugate pair counts as
% one mode: the member with Im >= 0 represents it.
%
% Nothing here depends on where the bulk edge is drawn: the edge factor
% comes afterwards from the null draws (outlier_null, outlier_measures),
% so the sorted real parts and moduli are kept (the n_keep largest of
% each) and the outliers are counted then.
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
%   R_eff, K_eig, K_eff   the predictions (outlier_effective_matrix)
%   dbar, drms, rate, frac_on   the population descriptors (same)
%
% See also: outlier_effective_matrix, outlier_null, outlier_run

if nargin < 4 || isempty(n_top);  n_top = 3;   end
if nargin < 5 || isempty(n_keep); n_keep = 30; end

[Mf, ~, ~, pred, pop] = outlier_effective_matrix(S, params, ws);
n = size(Mf, 1);
n_keep = min(n_keep, n);

ev = eig(Mf);
re_sorted  = sort(real(ev), 'descend');
abs_sorted = sort(abs(ev), 'descend');

modes = ev(imag(ev) >= 0);
[~, order] = sort(real(modes), 'descend');
modes = modes(order);
n_top = min(n_top, numel(modes));
ev_top = modes(1:n_top);

% Inverse iteration at a slightly displaced shift. The matrix is close to
% singular on purpose, so the two conditioning warnings are silenced here.
w1 = warning('off', 'MATLAB:nearlySingularMatrix');
w2 = warning('off', 'MATLAB:singularMatrix');
restore = onCleanup(@() warning([w1; w2]));   %#ok<NASGU>  restores on exit
V_top = complex(zeros(n, n_top));
res = zeros(n_top, 1); coh = zeros(n_top, 1);
I = eye(n);
for q = 1:n_top
    lam = ev_top(q);
    A = Mf - (lam + 1e-7 * (1 + abs(lam))) * I;
    v = complex(randn(n, 1), randn(n, 1));
    for it = 1:3
        v = A \ v;
        v = v / norm(v);
    end
    if ~all(isfinite(v))
        v = complex(nan(n, 1));
    end
    V_top(:, q) = v;
    res(q) = norm(Mf * v - lam * v);
    coh(q) = sum(abs(ws.U' * v).^2);
end

s = struct('ev_top', ev_top, 'V_top', single(V_top), 'res', res, 'coh', coh, ...
    're_keep', re_sorted(1:n_keep), 'abs_keep', abs_sorted(1:n_keep), ...
    'n_unstable', sum(real(ev) > 1), ...
    'R_eff', pred.R_eff, 'K_eig', pred.K_eig, 'K_eff', pred.K_eff, ...
    'dbar', pop.dbar, 'drms', pop.drms, 'rate', pop.rate, 'frac_on', pop.frac_on);
end
