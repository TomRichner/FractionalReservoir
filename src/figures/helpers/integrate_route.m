function S = integrate_route(model, P, r_t, dt, r_half, r_end)
% INTEGRATE_ROUTE Classical RK4 on a route model's own equations (model.rhs).
%
%   S = INTEGRATE_ROUTE(model, P, r_t, dt)
%   S = INTEGRATE_ROUTE(model, P, r_t, dt, r_half, r_end)
%
%   model, P  a compiled route (SRNNCellTypePairs.compile_route: R.model, R.P),
%             started from model.initial(P) (every pool full).
%   r_t       npre x nt, the rate of each neuron at the START of each step
%             (a 1 x nt row for one neuron).
%   dt        the step (s).
%   r_half, r_end  (optional) npre x nt, the rate at the middle and at the END
%             of each step, for a rate that varies within a step (a sinusoid).
%             Omitted, the rate is held at r_t over the step: a piecewise-
%             constant input that changes only at step boundaries, which is
%             exact for a staircase.
%
%   S         nt x ns x npre, the state at the END of each step, so that
%             S(i, :, j) is neuron j at t = i*dt. With one neuron that is the
%             nt x ns matrix fig_STD_steady_state reads.
%
% Every number is the route model's: the right-hand side is model.rhs, the
% one place a route's equations are coded (route_model, route_pool,
% route_product). A chain has no closed-form time course (its equations are
% bilinear), so every model is integrated the same way.
%
% Moved out of fig_STD_steady_state (where it was a local function holding
% the rate constant within each step) when fig_std_1TS_vs_MTS needed it for a
% sinusoidal rate and for many neurons at once. Called with four arguments
% and one neuron it performs exactly the arithmetic it did there.
%
% See also: fig_STD_steady_state, fig_std_1TS_vs_MTS, route_model

if nargin < 5
    r_half = r_t;
    r_end  = r_t;
end
[npre, nt] = size(r_t);
S = zeros(nt, P.ns, npre);
s = repmat(model.initial(P), npre, 1);
for i = 1:nt
    ri = r_t(:, i);
    rh = r_half(:, i);
    re = r_end(:, i);
    k1 = model.rhs(s, ri, P);
    k2 = model.rhs(s + 0.5 * dt * k1, rh, P);
    k3 = model.rhs(s + 0.5 * dt * k2, rh, P);
    k4 = model.rhs(s + dt * k3, re, P);
    s = s + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
    S(i, :, :) = reshape(s.', 1, P.ns, npre);
end
end
