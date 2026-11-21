function [Mz, d, offdiag] = freqresp_schur(J, params, ws, z, free)
%FREQRESP_SCHUR Numerical: the slow variables eliminated from the model's own Jacobian at one complex frequency.
%
%   [Mz, d, offdiag] = freqresp_schur(J, params, ws, z, free)
%
% Partition the full Jacobian J (SRNNCellTypePairs.compute_Jacobian_fast)
% into the dendritic rows x and the slow rows s (adaptation and route
% variables). For a mode going as exp(z t),
%
%   M(z) = I + tau_d [ J_xx + J_xs (z I - J_ss)^(-1) J_sx ]
%
% and z is an eigenvalue of J exactly when 1 + tau_d z is an eigenvalue of
% M(z) (for z not an eigenvalue of J_ss). Nothing here knows the model's
% equations: it is the check on the closed form of freqresp_gain.
%
% Inputs
%   J       N x N sparse, the full Jacobian at the state
%   params, ws   model.get_params() and outlier_weight_stats(params)
%   z       one complex frequency
%   free    which slow variables respond: .sfa and .std (both true by
%           default). The rows of one held are left out, which is what
%           holding it means.
%
% Outputs
%   Mz      n x n full complex
%   d       C x n complex: the gain of each presynaptic neuron on its
%           route to each postsynaptic type, fitted to Mz by least squares
%           over the column (as outlier_effective_matrix does for the
%           frozen gain): Mz(i, j) = W(i, j) d(a, j) for i of type a
%   offdiag norm(Mz - W .* d, 'fro') / norm(Mz, 'fro'): how far M(z) is
%           from "the weights times one gain per neuron and route". Zero
%           to round-off when the elimination is local.
%
% See also: freqresp_gain, freqresp_matrix

if nargin < 5 || isempty(free); free = struct(); end
if ~isfield(free, 'sfa'); free.sfa = true; end
if ~isfield(free, 'std'); free.std = true; end

L = params.state_layout;
C = params.n_cellTypes;
n = params.n;
slow = [];
if free.sfa
    slow = [slow, L.a{:}];
end
if free.std
    for pre = 1:C
        for post = 1:C
            slow = [slow, params.route{pre, post}.idx]; %#ok<AGROW>
        end
    end
end
slow = slow(:)';
x = L.x;

Mz = full(J(x, x));
if ~isempty(slow)
    ns = numel(slow);
    Mz = Mz + full(J(x, slow) * ((z .* speye(ns) - J(slow, slow)) \ J(slow, x)));
end
Mz = params.tau_d .* Mz;
Mz(1:n + 1:end) = Mz(1:n + 1:end) + 1;

d = complex(zeros(C, n));
fit = complex(zeros(n));
for a = 1:C
    rows = ws.ti{a};
    num = sum(Mz(rows, :) .* ws.Wf(rows, :), 1);
    den = ws.col_sq(a, :);
    ok = den > 0;
    d(a, ok) = num(ok) ./ den(ok);
    fit(rows, :) = ws.Wf(rows, :) .* d(a, :);
end
offdiag = norm(Mz - fit, 'fro') / norm(Mz, 'fro');
end
