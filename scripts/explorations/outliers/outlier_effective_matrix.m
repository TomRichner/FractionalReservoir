function [Mf, d, J, pred, pop] = outlier_effective_matrix(S, params, ws)
%OUTLIER_EFFECTIVE_MATRIX The frozen-gain effective connectivity at one state, its gains, and what theory predicts from them.
%
%   [Mf, d, J, pred, pop] = outlier_effective_matrix(S, params, ws)
%
% THE OBJECT. The dendritic block of the Jacobian is J_xx = (-I + M)/tau_d
% with M = W .* (gain of the presynaptic neuron on that route): the weights
% each neuron sends, scaled by the slope of its nonlinearity times its
% route's output factor (the available vesicle fraction), with the
% adaptation and depression variables held where they are. M is read off
% the class's own Jacobian, M = I + tau_d J_xx, so it is right for every
% route model.
%
% Outputs
%   Mf    n x n full, the effective connectivity (units of W: an eigenvalue
%         lambda of M is a rate (lambda - 1)/tau_d in 1/s)
%   d     C x n: d(a, j) is the gain of neuron j on its route to
%         postsynaptic type a, recovered exactly from M and W
%         (M(i, j) = W(i, j) d(a, j) for i of type a; least squares over
%         the column, which is exact). 0 where the neuron has no synapse
%         onto that type.
%   J     the sparse full Jacobian, for callers that want it
%   pred  struct, the predictions of the note's section 2 from these gains:
%           K_eff   C x C (post, pre): mu(a, b) sum_{j in b} d(a, j), the
%                   mean part of M
%           K_eig   its eigenvalues: where an outlier sits if it is outside
%                   the bulk
%           R_eff   bulk radius: sqrt of the Perron root of
%                   V(a, b) = s2(a, b) sum_{j in b} d(a, j)^2 (realized s2)
%   pop   struct of population descriptors at this state:
%           dbar    C x C (post, pre) mean gain, drms likewise (rms)
%           rate    1 x C mean rate by type
%           frac_on 1 x C fraction of each type with nonzero slope
%
% See also: outlier_weight_stats, outlier_snapshot, outlier_null

S  = S(:);
C  = params.n_cellTypes;
n  = params.n;
ti = ws.ti;
L  = params.state_layout;

J  = SRNNCellTypePairs.compute_Jacobian_fast(S, params);
Mf = full(J(L.x, L.x)) .* params.tau_d;
Mf(1:n + 1:end) = Mf(1:n + 1:end) + 1;

d = zeros(C, n);
for a = 1:C
    num = sum(Mf(ti{a}, :) .* ws.Wf(ti{a}, :), 1);
    den = ws.col_sq(a, :);
    ok = den > 0;
    d(a, ok) = num(ok) ./ den(ok);
end

K_eff = zeros(C); V = zeros(C); dbar = zeros(C); drms = zeros(C);
for a = 1:C
    for b = 1:C
        g = d(a, ti{b});
        K_eff(a, b) = ws.mu(a, b) * sum(g);
        V(a, b)     = ws.s2(a, b) * sum(g.^2);
        dbar(a, b)  = mean(g);
        drms(a, b)  = sqrt(mean(g.^2));
    end
end
pred = struct('K_eff', K_eff, 'K_eig', eig(K_eff), ...
    'R_eff', sqrt(max(0, max(real(eig(V))))));

if nargout >= 5
    % The rate and the slope, as compute_Jacobian_fast forms them.
    x_eff = S(L.x);
    for q = 1:C
        if params.n_a(q) > 0
            a_q = reshape(S(L.a{q}), params.n_per_type(q), params.n_a(q));
            x_eff(ti{q}) = x_eff(ti{q}) - params.c_eff(q) .* sum(a_q, 2);
        end
    end
    r  = params.activation_function(x_eff);
    rp = params.activation_function_derivative(x_eff);
    rate = zeros(1, C); frac_on = zeros(1, C);
    for q = 1:C
        rate(q) = mean(r(ti{q}));
        frac_on(q) = mean(rp(ti{q}) > 0);
    end
    pop = struct('dbar', dbar, 'drms', drms, 'rate', rate, 'frac_on', frac_on);
end
end
