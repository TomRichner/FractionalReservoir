function T = edens_radial_theory_blocks(s2, type_of, d, r)
%EDENS_RADIAL_THEORY_BLOCKS Radial eigenvalue distribution for block variances by (post, pre) type and a gain on every column.
%
%   T = edens_radial_theory_blocks(s2, type_of, d)
%   T = edens_radial_theory_blocks(s2, type_of, d, r)
%
% THE OBJECT. M(i, j) has zero mean and variance S(i, j) = s2(a, b) d_j^2,
% a the type of the postsynaptic neuron i and b that of the presynaptic
% neuron j. When the rows of s2 are equal this is a pure column scaling and
% edens_radial_theory applies (and the two agree: edens_check_theory).
% When they are not (the E-to-E block strengthened, say), the variance
% profile does not separate into a row factor times a column factor and
% the general equations are needed.
%
% THE EQUATIONS. For independent entries with variance profile S the
% spectrum is isotropic, and the proportion of eigenvalues beyond radius r
% comes from the positive solution v1, v2 (n x 1 each) of
%     v1 = (S' v1) ./ ((S v2) .* (S' v1) + r^2)
%     v2 = (S  v2) ./ ((S v2) .* (S' v1) + r^2),     sum(v1) = sum(v2)
% as  n_>(r) = (1 / n) v1' (S v2),  inside the disk of radius
% R = sqrt(spectral radius of S). These are the zero-regularization limit
% of the Dyson equation of the Hermitized matrix (the inhomogeneous
% circular law; Aljadeff, Stern and Sharpee 2015 treat the block case).
% For S(i, j) = s2_j / n they reduce to Ahmadian et al. (2015) Eq. 2.15,
% with n_> = g^2. This form of the equations was derived for this
% exploration and is CHECKED NUMERICALLY against sampled matrices in
% edens_check_theory; it is not copied from a reference.
%
% Because S depends on i only through its type, S v2 takes one value per
% type (A_a) and S' v1 is d_j^2 times one value per type (B_b), so the
% unknowns are 2 C numbers. They are found by damped fixed-point iteration
% in the logarithm, warm-started from the neighboring radius.
%
% Inputs
%   s2       C x C (post, pre) entry variances of W (zeros included)
%   type_of  n x 1 type of each neuron
%   d        n x 1 gains
%   r        radii (default 400 points from 0 to R)
%
% Output: struct with R, r, n_less (the atom of zero-gain columns
% included), rho_rel (density relative to the uniform disk of radius R,
% atom excluded), p0, n_iter (iterations used at each radius) and resid
% (the largest relative change in the last iteration).
%
% See also: edens_radial_theory, edens_check_theory

d = double(d(:)); type_of = type_of(:);
n = numel(d); C = size(s2, 1);
d2 = d.^2;
onehot = double(type_of == 1:C);                 % n x C
V = s2 .* (onehot' * d2)';                       % V(a, b) = s2(a, b) sum_{j in b} d_j^2
R = sqrt(max(real(eig(V))));
if nargin < 4 || isempty(r)
    r = linspace(0, R, 400);
end
r = double(r(:)');
p0 = mean(d == 0);

n_more = zeros(size(r)); n_iter = zeros(size(r)); resid = zeros(size(r));
A = ones(C, 1); B = ones(C, 1);
order = find(r > 0 & r < R);
[~, o] = sort(r(order));                          % from the center outwards
order = order(o);
max_iter = 20000; tol = 1e-10; damp = 0.5;
for k = order
    r2 = r(k)^2;
    for it = 1:max_iter
        At = A(type_of); Bt = B(type_of);
        den = At .* d2 .* Bt + r2;
        v1 = d2 .* Bt ./ den;
        v2 = At ./ den;
        c = sqrt(sum(v2) / sum(v1));               % fix the scale: sum(v1) = sum(v2)
        v1 = v1 * c; v2 = v2 / c;
        A_new = s2 * (onehot' * (d2 .* v2));      % A_a = sum_b s2(a, b) sum_{j in b} d_j^2 v2_j
        B_new = s2' * (onehot' * v1);             % B_b = sum_a s2(a, b) sum_{i in a} v1_i
        change = max(abs([log(A_new ./ A); log(B_new ./ B)]));
        A = exp((1 - damp) * log(A) + damp * log(A_new));
        B = exp((1 - damp) * log(B) + damp * log(B_new));
        if change < tol; break; end
    end
    At = A(type_of); Bt = B(type_of);
    den = At .* d2 .* Bt + r2;
    n_more(k) = mean((d2 .* Bt ./ den) .* At);
    n_iter(k) = it; resid(k) = change;
end
n_more(r == 0) = 1 - p0;
n_less = 1 - n_more;

u = (r / R).^2;
rho_rel = nan(size(r));
if numel(r) > 2
    rho_rel = gradient(n_less, u);
    rho_rel(r >= R) = NaN;
    k0 = find(r == 0);
    if ~isempty(k0) && numel(r) > k0
        rho_rel(k0) = (n_less(k0 + 1) - p0) / u(k0 + 1);
    end
end
T = struct('R', R, 'r', r, 'n_less', n_less, 'rho_rel', rho_rel, 'p0', p0, ...
    'n_iter', n_iter, 'resid', resid);
end
