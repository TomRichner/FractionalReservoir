function [d, d_inf, parts] = freqresp_gain(S, params, z, free)
%FREQRESP_GAIN Closed form: each neuron's gain at complex frequency z, with its slow variables eliminated.
%
%   [d, d_inf, parts] = freqresp_gain(S, params, z)
%   [...] = freqresp_gain(S, params, z, struct('sfa', true, 'std', false))
%
% THE OBJECT. Linearize the model about the state S and let every variable
% go as exp(z t). A neuron's adaptation and depression variables are
% driven by that neuron's own rate and by nothing else, so they can be
% solved for neuron by neuron, and what is left is the dendritic equation
%
%   (1 + tau_d z) dx_i = sum_j W(i, j) d(a, j; z) dx_j ,   a the type of i ,
%
% with one complex gain per presynaptic neuron j and postsynaptic type a
% (per route, because depression is a route's). z is an eigenvalue of the
% full Jacobian exactly when 1 + tau_d z is an eigenvalue of
% M(z) = W .* d(., .; z)  (freqresp_matrix). The gain is a product,
%
%   d(a, j; z) = phi'_j  *  G_sfa(j; z)  *  F(a, j)  *  G_std(a, j; z)
%
%   frozen gain (z -> infinity):  d_inf = phi' F, slope times the route's
%   output factor (the available vesicle fraction b_1, or the product).
%
%   adaptation, K timescales, c_eff = c / K:
%     G_sfa = 1 / (1 + phi' c_eff sum_k 1 / (1 + z tau_k))
%
%   depression, one pool (or a product of M pools), rate r:
%     G_std = 1 - sum_m (r / tau_rel_m) / (z + 1/tau_rec_m + r/tau_rel_m)
%   which for one pool is (z + 1/tau_rec) / (z + 1/tau_rec + r/tau_rel)
%
%   depression, chained pools (releasable b1, reserve b2):
%     G_std = 1 - (r/tau_rel) (z - A22) / ((z - A11)(z - A22) - A12 A21)
%     A11 = -b2/tau_rec1 - r/tau_rel     A12 = (1 - b1)/tau_rec1
%     A21 = b2/(n2 tau_rec1)             A22 = -1/tau_rec2 - (1 - b1)/(n2 tau_rec1)
%
% These are derived from the equations (Equations_stability_paper.md), not
% read off the code's Jacobian: freqresp_schur is the numerical route
% through the Jacobian, and the two are compared in
% explore_frequency_response.
%
% Inputs
%   S       N x 1 state
%   params  model.get_params()
%   z       complex frequencies, any shape (taken as a row of nz); 2 pi i f
%           for a real frequency f in Hz
%   free    which slow variables respond: .sfa and .std (both true by
%           default). One held is frozen where it is in S: its factor is 1.
%
% Outputs
%   d       C x n x nz complex: d(a, j, q) the gain of neuron j on its
%           route to postsynaptic type a at z(q)
%   d_inf   C x n real: the frozen gain
%   parts   struct: slope (n x 1), rate (n x 1), F (C x n), G_sfa (n x nz),
%           G_std (C x n x nz)
%
% Not covered (errors): facilitation, the geometric-mean combine, the zero
% floor, spontaneous release.
%
% See also: freqresp_schur, freqresp_matrix, freqresp_radius

if nargin < 4 || isempty(free); free = struct(); end
if ~isfield(free, 'sfa'); free.sfa = true; end
if ~isfield(free, 'std'); free.std = true; end

S = S(:);
z = reshape(z, 1, []);
nz = numel(z);
C = params.n_cellTypes;
n = params.n;
L = params.state_layout;
ti = params.type_indices;
if any(params.release_rate ~= 0) || params.std_zero_floor
    error('freqresp_gain:NotCovered', 'Spontaneous release and the zero floor are not covered.');
end

% The rate and the slope at this state.
x_eff = S(L.x);
for q = 1:C
    if params.n_a(q) > 0
        a_q = reshape(S(L.a{q}), params.n_per_type(q), params.n_a(q));
        x_eff(ti{q}) = x_eff(ti{q}) - params.c_eff(q) .* sum(a_q, 2);
    end
end
rate  = params.activation_function(x_eff);
slope = params.activation_function_derivative(x_eff);

% Adaptation.
G_sfa = ones(n, nz);
if free.sfa
    for q = 1:C
        if params.n_a(q) == 0; continue; end
        tau = params.tau_a_matrix{q};                     % n_q x K
        H = zeros(numel(ti{q}), nz);
        for k = 1:size(tau, 2)
            H = H + 1 ./ (1 + tau(:, k) * z);
        end
        G_sfa(ti{q}, :) = 1 ./ (1 + params.c_eff(q) .* slope(ti{q}) .* H);
    end
end

% Depression, route by route.
F = ones(C, n);
G_std = ones(C, n, nz);
for pre = 1:C
    r = rate(ti{pre});
    npre = numel(ti{pre});
    for post = 1:C
        P = params.route{pre, post}.P;
        if P.ns == 0; continue; end
        if P.n_g > 0 || (isfield(P, 'geomean') && any(P.geomean))
            error('freqresp_gain:NotCovered', 'Facilitation and the geometric-mean combine are not covered.');
        end
        s = reshape(S(params.route{pre, post}.idx), npre, P.ns);
        switch P.model
            case 'product'
                F(post, ti{pre}) = prod(s, 2);
                g = ones(npre, nz);
                for m = 1:P.n_b
                    use = r ./ P.tau_drive(:, m);
                    g = g - use ./ (z + 1 ./ P.tau_rest(:, m) + use);
                end
            case 'pool'
                b1 = s(:, 1);
                F(post, ti{pre}) = b1;
                use = r ./ P.tau_rel;
                if P.n_b == 1
                    g = 1 - use ./ (z + 1 ./ P.tau_rec1 + use);
                else
                    b2 = s(:, 2);
                    A11 = -b2 ./ P.tau_rec1 - use;
                    A12 = (1 - b1) ./ P.tau_rec1;
                    A21 = b2 ./ (P.n2 .* P.tau_rec1);
                    A22 = -1 ./ P.tau_rec2 - (1 - b1) ./ (P.n2 .* P.tau_rec1);
                    g = 1 - use .* (z - A22) ./ ((z - A11) .* (z - A22) - A12 .* A21);
                end
            otherwise
                error('freqresp_gain:NotCovered', 'Route model ''%s'' is not covered.', P.model);
        end
        if free.std
            G_std(post, ti{pre}, :) = reshape(g, 1, npre, nz);
        end
    end
end

d_inf = F .* slope';
d = d_inf .* reshape(G_sfa, 1, n, nz) .* G_std;
if nargout >= 3
    parts = struct('slope', slope, 'rate', rate, 'F', F, 'G_sfa', G_sfa, 'G_std', G_std);
end
end
