function res = eigdens_run_measures(X, grids)
%EIGDENS_RUN_MEASURES The eigenvalue density of one run at its calibration states: measured, null, theory, and the measures.
%
%   res = eigdens_run_measures(X)
%   res = eigdens_run_measures(X, grids)
%
% X is one run as eigdens_load_runs returns it. At every calibration state
% k, with d the saved per-neuron gains:
%   1. the theory: eigdens_radial_theory_blocks(X.s2, X.type_of, d), the
%      radial law of M = W diag(d) from the gains and the block variances of
%      W alone, and its bulk radius R_eff (the unit of every radius below);
%   2. the cumulative radial distribution (cdf_*) on grids.u and the density
%      relative to a uniform disk of radius R_eff (dens_*) in equal-area bins
%      grids.u2_edges of (r / R_eff)^2, for the measured spectrum of M
%      (ev_M), the null draws (averaged) and the theory; the cumulative
%      distribution also for W under one uniform gain (ev_uniform);
%   3. the measures of eigdens_measures ('shell', 'r50', 'n_band',
%      'n_unstable') of the same five: measured, null (averaged over the
%      draws), theory, a uniform disk of radius R_eff, and the uniform gain;
%   4. on the full Jacobian (X.dense, when present): eigdens_measures
%      'rates' on all eigenvalues and on those whose eigenvector has at least
%      half its norm^2 in the dendritic rows (the fast modes: n_band_x,
%      n_unstable_x);
%   5. descriptors of the gains (p0 the zero fraction, cv_nz, pr =
%      <d^2>^2 / <d^4>, mean, rms) per state, and by type (p0, mean and cv
%      of the nonzero gains, all states pooled).
% This is the loop body of the exploration's explore_eigenvalue_density
% (FD-042), with the gains READ rather than recovered from eigenvectors
% (the exploration's edens_run_gains is not ported) and every per-state
% curve kept, so the figures can take medians across seeds and states.
%
% grids, optional (the defaults are the exploration's):
%   u          0:0.01:1.4    radius over R_eff, for the cumulative curves
%   u2_edges   0:0.04:1.44   bins of (r / R_eff)^2 (equal areas)
%   u_theory   sqrt(linspace(0, 1, 101))   where the theory is solved
%
% Output res (K states): network, condition, seed, n, N, tau_d, lambda_1,
% t, R_eff, R_eff_saved, cdf_M, cdf_null, cdf_unigain, cdf_theory (K x
% numel(u)), dens_M, dens_null, dens_theory (K x bins), theory_gap (1 x K:
% max over r <= R_eff of |cdf_M - cdf_theory|), measure_names, m_M, m_null,
% m_theory, m_uniform, m_unigain (4 x K), full (struct of 1 x K counts, NaN
% without a dense pass), p0, cv_nz, pr, mean_gain, rms_gain (1 x K),
% gain_by_type (C x 3), gains (n x K, single), type_of, ev_example and
% R_eff_example (the spectrum of M at the middle state, for the plane
% figure), seconds.
%
% See also: eigdens_load_runs, run_eigenvalue_density_stage,
%           eigdens_radial_theory_blocks, eigdens_measures

if nargin < 2 || isempty(grids)
    grids = struct('u', 0:0.01:1.4, 'u2_edges', 0:0.04:1.44, 'u_theory', sqrt(linspace(0, 1, 101)));
end
t0 = tic;
S = eigdens_measures('settings');
if abs(X.tau_d - S.tau_d) > 1e-12
    % eigdens_measures states its rate band at tau_d = 0.1 s (the paper's).
    error('eigdens_run_measures:TauD', 'Run %s has tau_d = %g s; the measures are written for %g s.', ...
        X.file, X.tau_d, S.tau_d);
end
u = grids.u; u2_edges = grids.u2_edges; u_theory = grids.u_theory;
n = X.n; K = size(X.gain, 2); nb = numel(u2_edges) - 1;
names = {'shell', 'r50', 'n_band', 'n_unstable'};

cdf_M = zeros(K, numel(u)); cdf_null = cdf_M; cdf_unigain = cdf_M; cdf_theory = cdf_M;
dens_M = zeros(K, nb); dens_null = dens_M; dens_theory = dens_M;
R_eff = zeros(1, K); p0 = R_eff; cv_nz = R_eff; pr = R_eff; mean_gain = R_eff; rms_gain = R_eff;
theory_gap = R_eff;
mM = zeros(4, K); mNull = mM; mTheory = mM; mUniform = mM; mUnigain = mM;
inside = u <= 1;
for k = 1:K
    d = X.gain(:, k);
    Tk = eigdens_radial_theory_blocks(X.s2, X.type_of, d, []);
    Tk = eigdens_radial_theory_blocks(X.s2, X.type_of, d, u_theory * Tk.R);
    Re = Tk.R;
    R_eff(k) = Re;
    g = gain_stats(d);
    p0(k) = g.p0; cv_nz(k) = g.cv_nz; pr(k) = g.pr; mean_gain(k) = g.mean; rms_gain(k) = g.rms;
    if ~(Re > 0)
        % Every gain zero: M = 0, no bulk, nothing to measure in units of
        % its radius. The state is kept, as NaN, so the counts stay honest.
        [cdf_M(k, :), cdf_null(k, :), cdf_unigain(k, :), cdf_theory(k, :)] = deal(NaN);
        [dens_M(k, :), dens_null(k, :), dens_theory(k, :)] = deal(NaN);
        theory_gap(k) = NaN;
        [mM(:, k), mNull(:, k), mTheory(:, k), mUniform(:, k), mUnigain(:, k)] = deal(NaN);
        continue;
    end

    ev_M = X.ev_M(:, k); ev_u = X.ev_uniform(:, k); ev_null = X.ev_null(:, :, k);
    cdf_M(k, :) = radial_cdf(ev_M, Re, u);
    cdf_unigain(k, :) = radial_cdf(ev_u, Re, u);
    cdf_theory(k, :) = interp1([u_theory, 1.5], [Tk.n_less, 1], u);
    theory_gap(k) = max(abs(cdf_M(k, inside) - cdf_theory(k, inside)));
    dens_M(k, :) = radial_density(ev_M, Re, u2_edges);
    th = interp1([u_theory.^2, 2], [Tk.n_less, 1], u2_edges);
    th(1) = Tk.p0;                                    % leave the atom out of the first bin
    dens_theory(k, :) = diff(th) ./ diff(u2_edges);
    nd = size(ev_null, 2);
    m_null = zeros(4, nd);
    for q = 1:nd
        cdf_null(k, :) = cdf_null(k, :) + radial_cdf(ev_null(:, q), Re, u) / nd;
        dens_null(k, :) = dens_null(k, :) + radial_density(ev_null(:, q), Re, u2_edges) / nd;
        m_null(:, q) = as_column(eigdens_measures('spectrum', ev_null(:, q), Re), names);
    end
    mM(:, k) = as_column(eigdens_measures('spectrum', ev_M, Re), names);
    mNull(:, k) = mean(m_null, 2);
    mTheory(:, k) = as_column(eigdens_measures('radial', Tk.r, Tk.n_less, n), names);
    mUniform(:, k) = as_column(eigdens_measures('uniform', Re, n), names);
    mUnigain(:, k) = as_column(eigdens_measures('spectrum', ev_u, Re), names);
end

% The full Jacobian at the same states: all its eigenvalues, and those
% carried mostly (at least half the norm^2) by the dendritic variables.
full = struct('n_band', nan(1, K), 'n_unstable', nan(1, K), ...
    'n_band_x', nan(1, K), 'n_unstable_x', nan(1, K));
if ~isempty(X.dense)
    for k = 1:K
        lam = X.dense(k).lam; xf = X.dense(k).xfrac;
        a = eigdens_measures('rates', lam);
        b = eigdens_measures('rates', lam(xf >= 0.5));
        full.n_band(k) = a.n_band; full.n_unstable(k) = a.n_unstable;
        full.n_band_x(k) = b.n_band; full.n_unstable_x(k) = b.n_unstable;
    end
end

by_type = zeros(X.C, 3);                  % p0, mean of the nonzero gains, cv of the nonzero gains
for b = 1:X.C
    gb = gain_stats(reshape(X.gain(X.type_of == b, :), [], 1));
    by_type(b, :) = [gb.p0, gb.mean_nz, gb.cv_nz];
end

k_ex = ceil(K / 2);
res = struct('network', X.network, 'condition', X.condition, 'seed', X.seed, 'file', X.file, ...
    'n', n, 'N', X.N, 'tau_d', X.tau_d, 'lambda_1', X.lambda_1, 't', X.t, ...
    'R_eff', R_eff, 'R_eff_saved', X.R_eff_saved, ...
    'cdf_M', cdf_M, 'cdf_null', cdf_null, 'cdf_unigain', cdf_unigain, 'cdf_theory', cdf_theory, ...
    'dens_M', dens_M, 'dens_null', dens_null, 'dens_theory', dens_theory, 'theory_gap', theory_gap, ...
    'measure_names', {names}, 'm_M', mM, 'm_null', mNull, 'm_theory', mTheory, ...
    'm_uniform', mUniform, 'm_unigain', mUnigain, 'full', full, ...
    'p0', p0, 'cv_nz', cv_nz, 'pr', pr, 'mean_gain', mean_gain, 'rms_gain', rms_gain, ...
    'gain_by_type', by_type, 'gains', single(X.gain), 'type_of', X.type_of, ...
    'ev_example', X.ev_M(:, k_ex), 'R_eff_example', R_eff(k_ex), 'k_example', k_ex, ...
    'seconds', toc(t0));
end

%% ------------------------------------------------------------------------
function g = gain_stats(d)
% Descriptors of one gain vector. pr is the participation ratio of the
% squared gains, mean(d^2)^2 / mean(d^4): the density at the rim of the
% bulk relative to a uniform disk (Ahmadian et al. 2015, Eq. 2.18) when
% the weights have one variance. NaN where no gain is nonzero.
d = d(:); nz = d(d ~= 0);
g = struct('p0', mean(d == 0), 'mean', mean(d), 'rms', sqrt(mean(d.^2)), ...
    'mean_nz', NaN, 'cv_nz', NaN, 'pr', NaN);
if ~isempty(nz)
    g.mean_nz = mean(nz); g.cv_nz = std(nz, 1) / mean(nz);
    g.pr = mean(d.^2)^2 / mean(d.^4);
end
end

function c = radial_cdf(ev, R, u)
a = sort(abs(ev(:))) / R;
c = sum(a <= u, 1) / numel(a);
end

function h = radial_density(ev, R, u2_edges)
% Density relative to a uniform disk of radius R, in equal-area bins, the
% exact zeros (the zero-gain columns) left out.
a2 = (abs(ev(:)) / R).^2;
a2 = a2(a2 > 1e-16);
h = histcounts(a2, u2_edges) ./ numel(ev) ./ diff(u2_edges);
end

function v = as_column(m, names)
v = cellfun(@(f) m.(f), names(:));
end
