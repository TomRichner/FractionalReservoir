function ws = fresp_weight_stats(params)
%FRESP_WEIGHT_STATS What the frequency response needs from W and the state layout, computed once per network.
%
%   ws = fresp_weight_stats(params)      params = model.get_params()
%
% Fields
%   Wf        n x n, full(params.W): the matrix the dynamics use (gain and
%             any route scale included)
%   ti        1 x C cell, the neurons of each type
%   mu        C x C (post, pre): REALIZED mean of all entries of the block,
%             zeros included
%   s2        C x C (post, pre): REALIZED variance of all entries of the
%             block, zeros included: the per-entry variance the bulk sees
%   col_sq    C x n: sum over the post-type's rows of W(i, j)^2, the
%             denominator of the gain fit in fresp_schur
%   neuron_of_state   N x 1: the neuron every state variable belongs to
%             (x_j, a_jk, and the route variables of presynaptic neuron j)
%   is_x      N x 1 logical: the dendritic rows
%   K_struct_eig   the eigenvalues of the structural mean matrix
%             K(a, b) = n_b mu(a, b) (all gains 1), from the realized W
%   R_struct  the structural bulk radius from s2 (all gains 1), realized
%   K_theory_eig, R_theory   the same two for the ENSEMBLE the network was
%             drawn from: block means alpha mu_tilde and block variances
%             alpha (1 - alpha) mu_tilde^2 + alpha sigma_tilde^2 (Harris
%             Eq. 16 blockwise), times the gain level_of_chaos
%   has_structural_outlier   true when an eigenvalue of the ensemble's mean
%             matrix has real part beyond R_theory: the network is built
%             with a population mode outside the bulk before any gain acts
%             (E-to-E x 1.5 has one; the balanced reference network, whose
%             mean matrix is nilpotent, does not). Judged on the ensemble,
%             not on the realized W: in a 60-neuron network the realized
%             mean matrix of the reference network can reach past the
%             radius by chance (5.09 against 4.84 on one seed).
%
% The subset of the outlier exploration's outlier_weight_stats that the
% frequency response uses, under the stage's own name (FD-041).
%
% See also: fresp_gain, fresp_radius, fresp_schur

C  = params.n_cellTypes;
n  = params.n;
ti = params.type_indices;
L  = params.state_layout;

Wf = full(params.W);
mu = zeros(C); s2 = zeros(C); col_sq = zeros(C, n);
for a = 1:C
    for b = 1:C
        blk = Wf(ti{a}, ti{b});
        mu(a, b) = mean(blk(:));
        s2(a, b) = var(blk(:), 1);
    end
    col_sq(a, :) = sum(Wf(ti{a}, :).^2, 1);
end

% The neuron behind every state variable. a{q} is timescale-major and a
% route's rows are component-major, so each is the type's neuron list
% repeated.
neuron_of_state = zeros(L.N, 1);
neuron_of_state(L.x) = (1:n)';
for q = 1:C
    rows = L.a{q};
    if ~isempty(rows)
        neuron_of_state(rows) = repmat(ti{q}(:), numel(rows) / numel(ti{q}), 1);
    end
end
for pre = 1:C
    for post = 1:C
        for kind = {'b', 'g'}
            rows = L.(kind{1}){pre, post};
            if ~isempty(rows)
                neuron_of_state(rows) = repmat(ti{pre}(:), numel(rows) / numel(ti{pre}), 1);
            end
        end
    end
end
if any(neuron_of_state == 0)
    error('fresp_weight_stats:Layout', 'Some state variables were not assigned to a neuron.');
end
is_x = false(L.N, 1);
is_x(L.x) = true;

n_per = cellfun(@numel, ti);
K_struct_eig = eig(mu .* n_per(:)');                 % column b scaled by n_b
R_struct = sqrt(max(0, max(real(eig(s2 .* n_per(:)')))));

alpha = params.alpha;
gain = params.level_of_chaos;
K_theory_eig = eig((alpha .* params.mu_tilde .* gain) .* n_per(:)');
s2_theory = (alpha * (1 - alpha) * params.mu_tilde.^2 + alpha * params.sigma_tilde.^2) .* gain.^2;
R_theory = sqrt(max(0, max(real(eig(s2_theory .* n_per(:)')))));

ws = struct('Wf', Wf, 'ti', {ti}, 'mu', mu, 's2', s2, 'col_sq', col_sq, ...
    'neuron_of_state', neuron_of_state, 'is_x', is_x, ...
    'K_struct_eig', K_struct_eig, 'R_struct', R_struct, ...
    'K_theory_eig', K_theory_eig, 'R_theory', R_theory, ...
    'has_structural_outlier', max(real(K_theory_eig)) > R_theory);
end
