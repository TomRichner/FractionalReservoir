function ws = ospec_weight_stats(params)
%OSPEC_WEIGHT_STATS What the outlier measures need from W and the state layout, computed once per network.
%
%   ws = ospec_weight_stats(params)      params = model.get_params()
%
% Fields
%   Wf        n x n, full(params.W): the matrix the dynamics use (gain and
%             any self-coupling included)
%   ti        1 x C cell, the neurons of each type
%   type_of   n x 1, the type of each neuron
%   mu        C x C (post, pre): REALIZED mean of all entries of the block,
%             zeros included. This is mu_s of the note (after sparsity).
%   s2        C x C (post, pre): REALIZED variance of all entries of the
%             block, zeros included: the per-entry variance the bulk sees
%   s2_theory C x C (post, pre): the ensemble's, Harris Eq. 16 blockwise,
%             alpha (1 - alpha) mu^2 + alpha sigma^2, times level_of_chaos^2
%   col_sq    C x n: sum over the post-type's rows of W(i, j)^2, the
%             denominator of the gain fit in ospec_effective_matrix
%   U         n x C, orthonormal type indicators (column q is 1/sqrt(n_q)
%             on type q): the span of the population-uniform patterns
%   neuron_of_state   N x 1: the neuron every state variable belongs to
%             (x_j, a_jk, and the route variables of presynaptic neuron j)
%   is_x      N x 1 logical: the dendritic rows
%   K_struct, R_struct   the structural 2 x 2 (C x C) mean matrix
%             K(a, b) = n_b mu(a, b), its eigenvalues, and the structural
%             bulk radius from s2 (all gains 1), from the realized W
%
% Ported from scripts/explorations/outliers/outlier_weight_stats.m (FD-040), unchanged.
%
% See also: ospec_effective_matrix, ospec_snapshot

C  = params.n_cellTypes;
n  = params.n;
ti = params.type_indices;
L  = params.state_layout;

Wf = full(params.W);
type_of = zeros(n, 1);
U = zeros(n, C);
for q = 1:C
    type_of(ti{q}) = q;
    U(ti{q}, q) = 1 / sqrt(numel(ti{q}));
end

mu = zeros(C); s2 = zeros(C); col_sq = zeros(C, n);
for a = 1:C
    for b = 1:C
        blk = Wf(ti{a}, ti{b});
        mu(a, b) = mean(blk(:));
        s2(a, b) = var(blk(:), 1);
    end
    col_sq(a, :) = sum(Wf(ti{a}, :).^2, 1);
end

alpha = params.alpha;
s2_theory = (alpha * (1 - alpha) * params.mu_tilde.^2 + alpha * params.sigma_tilde.^2) ...
    .* params.level_of_chaos.^2;

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
    error('ospec_weight_stats:Layout', 'Some state variables were not assigned to a neuron.');
end
is_x = false(L.N, 1);
is_x(L.x) = true;

n_per = cellfun(@numel, ti);
K_struct = mu .* n_per(:)';                 % column b scaled by n_b
V_struct = s2 .* n_per(:)';

ws = struct('Wf', Wf, 'ti', {ti}, 'type_of', type_of, 'mu', mu, 's2', s2, ...
    's2_theory', s2_theory, 'col_sq', col_sq, 'U', U, ...
    'neuron_of_state', neuron_of_state, 'is_x', is_x, ...
    'K_struct', K_struct, 'K_struct_eig', eig(K_struct), ...
    'R_struct', sqrt(max(0, max(real(eig(V_struct))))), ...
    'R_struct_theory', sqrt(max(0, max(real(eig(s2_theory .* n_per(:)'))))));
end
