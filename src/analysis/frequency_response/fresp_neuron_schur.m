function d = fresp_neuron_schur(J, params, ws, j, z, free)
%FRESP_NEURON_SCHUR Numerical: one neuron's gain on its route to the excitatory type, from the Jacobian's blocks for that neuron alone.
%
%   d = fresp_neuron_schur(J, params, ws, j, z, free)
%
% With s the slow variables of neuron j (its adaptation variables and the
% variables of its routes) and i a postsynaptic neuron of type 1 that j has
% a synapse onto (the one with the largest weight),
%
%   d(z) = (tau_d / W(i, j)) [ J(x_i, x_j) + J(x_i, s) (z I - J(s, s))^(-1) J(s, x_j) ]
%
% Nothing here knows the model's equations: it is the check on the closed
% form of fresp_gain for a single neuron, at every frequency of z.
%
%   J        N x N sparse, the full Jacobian at the state
%   params, ws   model.get_params() and fresp_weight_stats(params)
%   j        the neuron
%   z        complex frequencies (1 x nz)
%   free     struct with .sfa and .std: which slow variables respond; the
%            rows of one held are left out
%   d        1 x nz complex
%
% See also: fresp_gain, fresp_schur, fresp_run_one

L = params.state_layout;
mine = find(ws.neuron_of_state == j & ~ws.is_x)';
keep = false(size(mine));
if free.sfa
    keep = keep | ismember(mine, [L.a{:}]);
end
if free.std
    route_rows = [];
    for pre = 1:params.n_cellTypes
        for post = 1:params.n_cellTypes
            route_rows = [route_rows, params.route{pre, post}.idx(:)']; %#ok<AGROW>
        end
    end
    keep = keep | ismember(mine, route_rows);
end
s = mine(keep);
post = params.type_indices{1};
post = post(post ~= j);
[w_best, best] = max(abs(ws.Wf(post, j)));
if isempty(w_best) || w_best == 0
    d = complex(nan(1, numel(z)));          % no synapse onto the excitatory type
    return;
end
i = post(best);
xi = L.x(i); xj = L.x(j);
Jss = full(J(s, s)); Jxs = full(J(xi, s)); Jsx = full(J(s, xj));
d = complex(zeros(1, numel(z)));
for q = 1:numel(z)
    v = full(J(xi, xj));
    if ~isempty(s)
        v = v + Jxs * ((z(q) * eye(numel(s)) - Jss) \ Jsx);
    end
    d(q) = params.tau_d * v / ws.Wf(i, j);
end
end
