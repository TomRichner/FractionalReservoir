function [t, Y] = ode_rk4(odefun, tspan, y0, ~)
%ODE_RK4 Fixed-step classic RK4 with an @ode45-compatible call signature.
%   [t, Y] = ode_rk4(odefun, tspan, y0, opts)
%
%   Drop-in replacement for ode45 as used by SRNNModel2:
%       solver(rhs, t_ex, S0, opts)
%   Steps RK4 at the NATIVE spacing of the supplied time vector tspan
%   (assumed to be a full uniform fs grid) and returns the solution at
%   exactly those times, so SRNNModel2.run()'s output-time check passes.
%   The 4th argument (odeset options) is accepted for signature
%   compatibility and IGNORED (RK4 is fixed-step; no tolerances/Jacobian).
%
%   REQUIREMENT: tspan must be the full time grid, not a 2-point [t0 tf]
%   span. A 2-point span would be integrated as a single giant step and
%   silently produce garbage, so it is rejected below. (The QR Lyapunov
%   variational path passes such 2-point spans; use ode45 there.)
%
%   Much faster than adaptive ode45 with noisy forcing (no step-size
%   control thrashing). odefun is evaluated 4x per step.

    if numel(tspan) < 3
        error('ode_rk4:SpanTooShort', ...
            ['ode_rk4 requires a full time grid (numel(tspan) >= 3); got %d ', ...
             'point(s). A 2-point [t0 tf] span (e.g. the QR variational-eq ', ...
             'path) is not supported by this fixed-step integrator -- use ode45.'], ...
            numel(tspan));
    end

    t  = tspan(:);
    nt = numel(t);
    y  = y0(:);
    Y  = zeros(nt, numel(y));
    Y(1, :) = y.';
    for k = 1:nt-1
        h  = t(k+1) - t(k);
        tk = t(k);
        k1 = odefun(tk,       y);
        k2 = odefun(tk + h/2, y + (h/2)*k1);
        k3 = odefun(tk + h/2, y + (h/2)*k2);
        k4 = odefun(tk + h,   y + h*k3);
        y  = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        Y(k+1, :) = y.';
    end
end
