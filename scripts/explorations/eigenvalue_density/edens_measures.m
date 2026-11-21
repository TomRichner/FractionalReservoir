function m = edens_measures(kind, varargin)
%EDENS_MEASURES The expressivity measures, for a spectrum, for a radial distribution, and for a uniform disk.
%
%   m = edens_measures('spectrum', ev, R)            ev: eigenvalues of M (units of W)
%   m = edens_measures('radial', r, n_less, n)       a radial law: proportion within r
%   m = edens_measures('uniform', R, n)              the uniform disk of radius R
%   m = edens_measures('rates', lam)                 lam: eigenvalues of a Jacobian (1/s)
%   s = edens_measures('settings')                   the constants below
%
% WHAT IS MEASURED. M = W diag(d) is the frozen-gain effective
% connectivity; an eigenvalue lambda of M is the rate (lambda - 1) / tau_d
% of the dendritic block of the Jacobian, so the stability boundary is the
% line Re(lambda) = 1.
%
%   shell     proportion of the eigenvalues with modulus above shell_r
%             (0.8) times the bulk radius R: the outer 36% of the disk's
%             area. SHAPE ONLY: it does not move when every gain is scaled
%             by one number. A uniform disk gives 0.36; shell / 0.36 is the
%             occupancy of the rim relative to a uniform disk. Blind to
%             where the disk sits relative to the boundary: a flat disk
%             that lies wholly inside the stable half-plane scores 1.
%   r50       the median modulus divided by R (uniform disk: 0.707)
%   n_band    number of eigenvalues with 1 - band < Re(lambda) <= 1
%             (band = 0.2: decay rates slower than 2 1/s at tau_d = 0.1 s):
%             the modes one small change of gain away from expanding.
%   n_unstable  number with Re(lambda) > 1: the modes that expand now,
%             while the gains stay where they are.
%             Both counts mix SIZE and SHAPE (a larger disk puts more of
%             itself near the boundary), which is what the dynamics see,
%             and both are blind to what the eigenvectors are (non-normal
%             growth, which neurons a mode involves) and to how long the
%             gains stay put.
%
% For a radial law the counts are expectations: an isotropic distribution
% puts the fraction acos(c / r) / pi of the eigenvalues at modulus r to the
% right of Re = c, so P(Re > c) is the integral of that against the radial
% distribution. For 'rates' the same two counts are taken on a Jacobian's
% own eigenvalues: n_band is -band / tau_d < Re <= 0 (1/s; tau_d = 0.1 s)
% and n_unstable is Re > 0.
%
% See also: explore_eigenvalue_density, edens_radial_theory

S = struct('shell_r', 0.8, 'band', 0.2, 'tau_d', 0.1);
switch kind
    case 'settings'
        m = S;
    case 'spectrum'
        [ev, R] = varargin{:};
        a = abs(ev(:));
        m = struct('shell', mean(a > S.shell_r * R), 'r50', median(a) / R, ...
            'n_band', nnz(real(ev) > 1 - S.band & real(ev) <= 1), ...
            'n_unstable', nnz(real(ev) > 1));
    case 'radial'
        [r, n_less, n] = varargin{:};
        r = r(:)'; n_less = n_less(:)';
        R = max(r);
        p_right = @(c) right_of(r, n_less, c);
        m = struct('shell', 1 - interp1(r, n_less, S.shell_r * R), ...
            'r50', interp1(n_less + 1e-12 * (1:numel(r)), r, 0.5, 'linear', 0) / R, ...
            'n_band', n * (p_right(1 - S.band) - p_right(1)), ...
            'n_unstable', n * p_right(1));
    case 'uniform'
        [R, n] = varargin{:};
        p_right = @(c) uniform_right_of(c / R);
        m = struct('shell', 1 - S.shell_r^2, 'r50', sqrt(0.5), ...
            'n_band', n * (p_right(1 - S.band) - p_right(1)), ...
            'n_unstable', n * p_right(1));
    case 'rates'
        lam = varargin{1};
        m = struct('n_band', nnz(real(lam) > -S.band / S.tau_d & real(lam) <= 0), ...
            'n_unstable', nnz(real(lam) > 0));
    otherwise
        error('edens_measures:Kind', 'Unknown kind ''%s''.', kind);
end
end

function p = right_of(r, n_less, c)
% P(Re > c) for an isotropic law with cumulative radial distribution
% n_less(r), c > 0: sum over the radial increments of acos(c / r) / pi.
dF = diff(n_less);
rm = (r(1:end - 1) + r(2:end)) / 2;
frac = zeros(size(rm));
ok = rm > c;
frac(ok) = acos(c ./ rm(ok)) / pi;
p = sum(dF .* frac);
end

function p = uniform_right_of(u)
% Fraction of a uniform unit disk to the right of Re = u.
u = min(max(u, -1), 1);
p = (acos(u) - u .* sqrt(1 - u.^2)) / pi;
end
