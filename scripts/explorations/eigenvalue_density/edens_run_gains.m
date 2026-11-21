function G = edens_run_gains(R, idx)
%EDENS_RUN_GAINS The per-neuron frozen gains of a saved outlier run, recovered from its stored eigenvectors (no simulation).
%
%   G = edens_run_gains(R)          at the run's calibration states (R.cal_t)
%   G = edens_run_gains(R, idx)     at the sampled states idx (columns of R.ev_top)
%
% WHY. outlier_run saves the spectrum of the frozen-gain effective
% connectivity M = W diag(d) at its calibration states, but of the gains
% only their mean and rms by type. The density of M depends on the whole
% distribution of d, so it is recovered here.
%
% HOW. W is redrawn from the run's preset, condition and seed (build only;
% nothing is integrated), and the run stored the n_top rightmost
% eigenvalues of M with their right eigenvectors at every sampled state.
% Each pair satisfies W diag(v) d = lambda v, which is linear in d: 2 n
% real equations per eigenvector for n unknowns. The stacked system is
% solved by least squares. It is exact up to the single precision of the
% stored eigenvectors: at the calibration states the eigenvalues of
% W diag(d) reproduce the stored spectrum of M (G.check_ev, of order 1e-6),
% and the mean and rms of d by type reproduce R.dbar and R.drms
% (G.check_moments). Gains below 1e-5 in magnitude are set to exactly 0 (a
% neuron on the floor or the ceiling of the nonlinearity); the count is
% checked against the number of zero eigenvalues of the stored M
% (G.check_zeros, the difference in counts).
%
% ASSUMES one gain per neuron (a neuron's routes to the two postsynaptic
% types carry the same factor), which holds in these runs: R.dbar and
% R.drms have equal rows. The check on the spectrum would fail otherwise.
% Networks with added low-rank structure (OutlierNet) are not supported;
% see network_args below for how the run's network is rebuilt.
%
% Output: struct with
%   d        n x numel(idx) gains
%   idx, t   the sampled states and their times
%   W        n x n full, the redrawn connectivity
%   ti       1 x C cell, neurons by type;  type_of  n x 1
%   mu, s2   C x C (post, pre) realized block mean and variance of W,
%            zeros included (outlier_weight_stats's definitions)
%   v_col    n x 1: n times the entry variance of column j, averaged over
%            the postsynaptic neurons, so that v_col .* d.^2 is what
%            edens_radial_theory takes WHEN the rows of s2 are equal (the
%            reference network); otherwise use edens_radial_theory_blocks
%   R_struct the structural bulk radius of the redrawn W (all gains 1)
%   check_ev, check_moments, check_zeros   as above (NaN where the run has
%            no stored spectrum at that state)
%   n_unseen 1 x numel(idx): neurons absent from every stored eigenvector
%            at that state, whose gain could not be recovered and is left
%            at 0 (nonzero only when nearly every neuron is pinned; a
%            check_moments above about 1e-6 at that state says one of them
%            had a nonzero gain)
%
% See also: outlier_run, outlier_effective_matrix, edens_radial_theory

if isfield(R, 'net') && ~isempty(R.net.added_structure)
    error('edens_run_gains:AddedStructure', 'Runs with added low-rank structure are not supported.');
end
[~, cal_idx] = min(abs(R.t(:) - R.cal_t(:)'), [], 1);
if nargin < 2 || isempty(idx)
    idx = cal_idx;
end

net_args = network_args(R);
model = build_from_preset(R.preset, R.condition, net_args{:}, ...
    'rng_seeds', R.seed, 'T_range', [0 1], 'fs', R.opts.fs, 'verbose', 'near-none');
params = model.get_params();
if abs(model.R - R.structural.R) > 1e-9 * max(1, R.structural.R)
    error('edens_run_gains:WrongNetwork', 'The redrawn network has R = %.6g, the run had %.6g.', model.R, R.structural.R);
end
W = full(params.W);
n = params.n; C = params.n_cellTypes; ti = params.type_indices;
type_of = zeros(n, 1);
mu = zeros(C); s2 = zeros(C);
for a = 1:C
    type_of(ti{a}) = a;
    for b = 1:C
        blk = W(ti{a}, ti{b});
        mu(a, b) = mean(blk(:));
        s2(a, b) = var(blk(:), 1);
    end
end
n_per = cellfun(@numel, ti);
v_col = zeros(n, 1);
for b = 1:C
    v_col(ti{b}) = sum(n_per(:) .* s2(:, b));
end
% The structural bulk radius is the root of the largest eigenvalue of the
% block variances times the type sizes (outlier_weight_stats's R_struct);
% sqrt(mean(v_col)) equals it only when the rows of s2 are equal.
R_struct = sqrt(max(real(eig(s2 .* n_per(:)'))));
if abs(R_struct - R.structural.R_realized) > 1e-6 * R_struct
    error('edens_run_gains:WrongWeights', 'The redrawn W has realized radius %.6g, the run had %.6g.', R_struct, R.structural.R_realized);
end

n_s = numel(idx);
d = zeros(n, n_s);
check_ev = nan(1, n_s); check_zeros = nan(1, n_s); check_moments = zeros(1, n_s);
n_unseen = zeros(1, n_s);
for k = 1:n_s
    lam = double(R.ev_top(:, idx(k)));
    V = double(R.V_top(:, :, idx(k)));
    A = zeros(2 * n * numel(lam), n); rhs = zeros(2 * n * numel(lam), 1);
    for q = 1:numel(lam)
        Aq = W .* V(:, q).';
        bq = lam(q) * V(:, q);
        rows = (q - 1) * 2 * n + (1:2 * n);
        A(rows, :) = [real(Aq); imag(Aq)];
        rhs(rows) = [real(bq); imag(bq)];
    end
    % A neuron whose entry is zero in every stored eigenvector cannot be
    % seen: its column of A vanishes. That happens when almost every gain
    % is zero, so that M has only a few nonzero columns and a neuron that
    % receives from none of them is absent from every eigenvector with a
    % nonzero eigenvalue. Its gain is left at 0 and it is counted in
    % n_unseen; whether that was right shows in check_moments.
    col = sqrt(sum(A.^2, 1));
    seen = col > 1e-6 * max(col);
    n_unseen(k) = nnz(~seen);
    dk = zeros(n, 1);
    dk(seen) = A(:, seen) \ rhs;
    dk(abs(dk) < 1e-5) = 0;
    d(:, k) = dk;
    for b = 1:C
        check_moments(k) = max([check_moments(k), ...
            abs(mean(dk(ti{b})) - R.dbar(1, b, idx(k))), ...
            abs(sqrt(mean(dk(ti{b}).^2)) - R.drms(1, b, idx(k)))]);
    end
    c = find(cal_idx == idx(k), 1);
    if ~isempty(c)
        ev_M = R.cal(c).ev_M;
        ev = eig(W .* dk');
        check_ev(k) = max(abs(sort(abs(ev)) - sort(abs(ev_M))));
        check_zeros(k) = nnz(dk == 0) - nnz(abs(ev_M) < 1e-9);
    end
end

G = struct('d', d, 'idx', idx, 't', R.t(idx), 'W', W, 'ti', {ti}, 'type_of', type_of, ...
    'mu', mu, 's2', s2, 'v_col', v_col, 'R_struct', R_struct, ...
    'check_ev', check_ev, 'check_moments', check_moments, 'check_zeros', check_zeros, ...
    'n_unseen', n_unseen);
end

function net_args = network_args(R)
% The constructor arguments that make the run's network. A run saved since
% 2026-11-20 carries them (R.net, from outlier_network). An earlier run
% does not, and what it stored as structural.mu_tilde_relative is the
% preset's value, not the network's, so its name is parsed as outlier_run
% did then: 'reference' changes nothing, and '<XY>x<a>p<b>' sets the
% property mu_<XY>_relative to a.b times the preset's 7. Whichever is
% used, the caller checks the redrawn network against the run's radius
% and, at the calibration states, against its stored spectrum.
if isfield(R, 'net')
    net_args = R.net.net_args;
    return;
end
net_args = {};
tokens = strsplit(R.network, '-');
for k = 1:numel(tokens)
    tok = tokens{k};
    if strcmp(tok, 'reference'); continue; end
    parts = regexp(tok, '^(EE|EI|IE|II)x(\d+)p(\d+)$', 'tokens', 'once');
    if isempty(parts)
        error('edens_run_gains:Network', 'Network token ''%s'' is not supported.', tok);
    end
    net_args = [net_args, {['mu_' parts{1} '_relative'], 7 * str2double([parts{2} '.' parts{3}])}]; %#ok<AGROW>
end
end
