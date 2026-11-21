function C = eigdens_check_theory(opts)
%EIGDENS_CHECK_THEORY The radial-density formulas against sampled random matrices, and against each other.
%
%   C = eigdens_check_theory()
%   C = eigdens_check_theory('n', 500, 'n_draws', 20, 'seed', 1)
%
% Five cases, each a matrix X diag(d) with X i.i.d. Gaussian of variance
% 1 / n (cases 1 to 4) or a block variance profile (case 5):
%   1 uniform       every gain 1: the circular law
%   2 broad         gains log-normal, coefficient of variation 1
%   3 zeros         half the gains 0, the rest 1: an atom of mass 1/2 at
%                   the origin and a uniform disk of radius sqrt(1/2)
%   4 two types     Harris et al. (2023), Fig. 2(c): f = 0.8, mu_e =
%                   sigma_e = 1, mu_i = -4, sigma_i = 4, alpha = 0.3
%                   (variances only; the mean is left out)
%   5 blocks        variances that depend on (post, pre) type, with broad
%                   gains and a fifth of them zero: the non-separable case
%                   that eigdens_radial_theory_blocks is for
% For each: the sampled cumulative radial distribution (n_draws matrices
% of size n) and the theory's, the largest difference between them
% (max_diff, to be read against the sampling error 1 / sqrt(n n_draws),
% returned as C(c).sampling_error), and for cases 1 to 4 the agreement of
% eigdens_radial_theory with eigdens_radial_theory_blocks
% (max_diff_solvers) and, for case 4, with Harris's closed form
% (max_diff_harris, in density relative to uniform).
%
% Ported from the exploration's edens_check_theory (FD-042) with the
% figure taken out: fig_eigdens_theory_check draws it. The draws use a
% private stream (mt19937ar, opts.seed), so the global generator is not
% touched and the result is the same in a session and on a pool worker.
% n must be even (the block case splits the neurons in two halves).
%
% Output: struct array, one element per case, with fields name, u (radius
% over the bulk radius R), sampled, theory, R, max_diff, sampling_error,
% max_diff_solvers, max_diff_harris (NaN where not applicable) and
% edge_mean_max (mean largest modulus over R).
%
% See also: eigdens_radial_theory, eigdens_radial_theory_blocks,
%           eigdens_harris_density, fig_eigdens_theory_check

arguments
    opts.n       (1,1) double {mustBePositive, mustBeInteger} = 500
    opts.n_draws (1,1) double {mustBePositive, mustBeInteger} = 20
    opts.seed    (1,1) double = 1
end
n = opts.n;
if mod(n, 2) ~= 0
    error('eigdens_check_theory:OddSize', 'n must be even (the block case splits the neurons in two halves); got %d.', n);
end
stream = RandStream('mt19937ar', 'Seed', opts.seed);
u = linspace(0, 1, 201);
type_of = [ones(n / 2, 1); 2 * ones(n / 2, 1)];

H = eigdens_harris_density(0.3, 0.8, 1, -4, 1, 4);
n_e = round(0.8 * n);
d_broad = exp(randn(stream, n, 1) * sqrt(log(2)));       % log-normal, cv 1
d_zero = [zeros(n / 2, 1); ones(n / 2, 1)];
d_mixed = rand(stream, n, 1).^2; d_mixed(rand(stream, n, 1) < 0.2) = 0;
cases = { ...
    'uniform',   ones(n, 1),                                         ones(2) / n; ...
    'broad',     d_broad,                                            ones(2) / n; ...
    'zeros',     d_zero(randperm(stream, n)),                        ones(2) / n; ...
    'two types', [sqrt(H.v_e) * ones(n_e, 1); sqrt(H.v_i) * ones(n - n_e, 1)], ones(2) / n; ...
    'blocks',    d_mixed,                                            [0.10 0.046; 0.046 0.03]};

C = struct('name', cases(:, 1)', 'u', u, 'sampled', [], 'theory', [], 'R', [], ...
    'max_diff', [], 'sampling_error', 1 / sqrt(n * opts.n_draws), ...
    'max_diff_solvers', NaN, 'max_diff_harris', NaN, 'edge_mean_max', []);
for c = 1:size(cases, 1)
    d = cases{c, 2}; s2 = cases{c, 3};
    T = eigdens_radial_theory_blocks(s2, type_of, d, []);
    T = eigdens_radial_theory_blocks(s2, type_of, d, u * T.R);
    sd = sqrt(s2(type_of, type_of)) .* d';
    acc = zeros(size(u)); mx = 0;
    for k = 1:opts.n_draws
        a = sort(abs(eig(randn(stream, n) .* sd)));
        acc = acc + sum(a <= u * T.R, 1) / n / opts.n_draws;
        mx = mx + a(end) / opts.n_draws;
    end
    C(c).sampled = acc; C(c).theory = T.n_less; C(c).R = T.R;
    C(c).max_diff = max(abs(acc - T.n_less));
    C(c).edge_mean_max = mx / T.R;
    if c <= 4
        T1 = eigdens_radial_theory(d.^2, u * T.R);
        C(c).max_diff_solvers = max(abs(T1.n_less - T.n_less));
    end
    if c == 4
        T1 = eigdens_radial_theory(d.^2, H.r);
        C(c).max_diff_harris = max(abs(T1.rho_rel(2:end - 1) - H.rho_rel(2:end - 1)));
    end
end
end
