function S = fresp_matched_state(params, r0)
%FRESP_MATCHED_STATE A synthetic state: every neuron at rate r0 with its slow variables at their steady state for that rate.
%
%   S = fresp_matched_state(params, r0)
%
% For comparing conditions at one operating point. The rate r0 is put on
% every neuron (0.1 < r0 < 0.9 is the linear part of the nonlinearity at
% S_a = 0.8, slope 1), each adaptation variable is at r0, and each route's
% variables are at the route model's own steady state for r0. The
% dendritic variable is whatever gives that rate. On one network seed a
% neuron has the same fast adaptation constant and the same
% releasable-pool times in the single- and multiple-timescale conditions,
% so the comparison is neuron by neuron.
%
% It is not a state the network visits (there the rates differ by neuron
% and the slow variables lag): it isolates what the timescale structure
% does at a given slope and rate.
%
% See also: fresp_gain, fresp_matched

C = params.n_cellTypes;
n = params.n;
L = params.state_layout;
ti = params.type_indices;

% x_eff with phi(x_eff) = r0, by bisection (phi is monotone; its center
% may be per neuron and is inside the handle).
lo = -5 * ones(n, 1); hi = 5 * ones(n, 1);
for it = 1:80
    mid = (lo + hi) / 2;
    above = params.activation_function(mid) > r0;
    hi(above) = mid(above);
    lo(~above) = mid(~above);
end
x_eff = (lo + hi) / 2;

S = zeros(L.N, 1);
x = x_eff;
for q = 1:C
    if params.n_a(q) > 0
        S(L.a{q}) = r0;
        x(ti{q}) = x(ti{q}) + params.c_eff(q) * params.n_a(q) * r0;
    end
end
S(L.x) = x;
for pre = 1:C
    for post = 1:C
        R = params.route{pre, post};
        if R.P.ns == 0; continue; end
        ss = R.model.steady_state(r0 * ones(numel(ti{pre}), 1), R.P);
        S(R.idx) = ss.s(:);
    end
end
end
