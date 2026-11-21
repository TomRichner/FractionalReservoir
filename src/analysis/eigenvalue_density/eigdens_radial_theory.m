function T = eigdens_radial_theory(s2, r)
%EIGDENS_RADIAL_THEORY Radial eigenvalue distribution of a random matrix whose columns carry different variances.
%
%   T = eigdens_radial_theory(s2)
%   T = eigdens_radial_theory(s2, r)
%
% THE OBJECT. M = X diag(d): X has independent zero-mean entries of variance
% v / n, and column j is scaled by the gain d_j. Write s2(j) = v d_j^2, the
% variance of column j times n. A low-rank mean added to M moves at most a
% few eigenvalues (the outliers) and leaves this distribution as it is.
%
% THE RESULT (Ahmadian, Fumarola and Miller 2015, Eqs. 2.13 to 2.18, the
% case M = 0 of A = M + L J R with R L = diag(sqrt(s2)); the two-population
% case is Rajan and Abbott 2006, and Harris et al. 2023, Eqs. 19 to 23, is
% that case with the variances of a sparse matrix). The spectrum is
% isotropic about the origin and fills a disk of radius
%     R = sqrt(mean(s2))                                    (2.13)
% and the proportion of eigenvalues FURTHER than r from the origin,
% n_>(r) = g^2, solves
%     mean( s2 ./ (r^2 + g^2 s2) ) = 1,   0 < r < R         (2.15, 2.16)
% The density is  rho(r) = -(1 / (2 pi r)) d(g^2)/dr       (2.14)
% with the two ends in closed form:
%     rho(0) = mean(1 ./ s2) / pi                           (2.17)
%     rho(R) = R^2 / (pi mean(s2.^2))                       (2.18)
% Relative to a uniform disk of the same radius (density 1 / (pi R^2)):
%     center   mean(s2) mean(1 ./ s2)  >= 1
%     edge     mean(s2)^2 / mean(s2.^2) <= 1
% The edge value is the participation ratio of the squared gains: the
% fraction of neurons that effectively carry the bulk.
%
% ZERO GAINS. A column with s2 = 0 is a zero column of M, so M has an exact
% zero eigenvalue for each: an atom of mass p0 = mean(s2 == 0) at the
% origin. Equation 2.15 gives it without special treatment: as r -> 0,
% g^2 -> 1 - p0. (2.17 is then infinite; rho0_rel below is for the
% columns with nonzero gain.)
%
% Inputs
%   s2   vector, n x column variances (any scale; R is in the same units)
%   r    radii at which to evaluate (default: 400 points from 0 to R)
%
% Output: struct with fields
%   R         disk radius, sqrt(mean(s2))
%   r         the radii (row)
%   n_less    proportion of eigenvalues with modulus <= r, the atom at
%             zero included (so n_less(r = 0+) = p0)
%   rho_rel   density relative to the uniform disk of radius R, rho(r)
%             pi R^2, the atom excluded (numerical derivative of n_less
%             with respect to (r / R)^2); NaN beyond R
%   p0        mass of the atom at the origin
%   edge_rel  rho(R) pi R^2 = mean(s2)^2 / mean(s2.^2)
%   rho0_rel  density at the center relative to uniform, from the columns
%             with nonzero gain: (1 - p0) mean(s2) mean(1 ./ s2(s2 > 0))
%
% See also: eigdens_harris_density, eigdens_check_theory

s2 = double(s2(:));
if any(s2 < 0) || ~any(s2 > 0)
    error('eigdens_radial_theory:BadVariances', 's2 must be nonnegative with at least one positive entry.');
end
R = sqrt(mean(s2));
if nargin < 2 || isempty(r)
    r = linspace(0, R, 400);
end
r = double(r(:)');
p0 = mean(s2 == 0);
sp = s2(s2 > 0);
w = numel(sp) / numel(s2);          % weight of the nonzero columns in the mean

% Solve mean(s2 ./ (r^2 + q s2)) = 1 for q = g^2 by bisection; the left side
% falls monotonically in q, is above 1 at q -> 0 when r < R and equals
% (1 - p0) / q < 1 at q = 1 when r > 0.
q = zeros(size(r));
inside = r < R;
lo = zeros(1, nnz(inside)); hi = ones(1, nnz(inside));
rr2 = r(inside).^2;
for it = 1:60
    mid = (lo + hi) / 2;
    f = w * mean(sp ./ (rr2 + mid .* sp), 1) - 1;
    up = f > 0;
    lo(up) = mid(up); hi(~up) = mid(~up);
end
q(inside) = (lo + hi) / 2;
q(r == 0) = 1 - p0;
n_less = 1 - q;

% Density relative to uniform: d n_less / d (r/R)^2. For a uniform disk
% n_less = (r/R)^2, so this is 1.
u = (r / R).^2;
rho_rel = nan(size(r));
if numel(r) > 2
    rho_rel = gradient(n_less, u);
    rho_rel(~inside) = NaN;
    k0 = find(r == 0);
    if ~isempty(k0) && numel(r) > k0
        rho_rel(k0) = (n_less(k0 + 1) - p0) / u(k0 + 1);   % one-sided, the atom excluded
    end
end

T = struct('R', R, 'r', r, 'n_less', n_less, 'rho_rel', rho_rel, 'p0', p0, ...
    'edge_rel', mean(s2)^2 / mean(s2.^2), ...
    'rho0_rel', (1 - p0) * mean(s2) * mean(1 ./ sp));
end
