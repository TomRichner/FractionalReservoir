function [t, Y] = sde_fixed_step(odefun, tspan, y0, ~, noise, scheme)
%SDE_FIXED_STEP Fixed-step SDE integrators for ADDITIVE noise.
%   [t, Y] = sde_fixed_step(odefun, tspan, y0, opts, noise, scheme)
%
%   Integrates   dy = f(t,y) dt + sigma dW   where sigma is a constant scalar
%   applied to the state components named by noise.idx. Same call shape as
%   ode_rk4 plus the two trailing arguments; the odeset options are accepted
%   for signature compatibility and IGNORED (fixed step, no tolerances).
%
%   SCHEME is 'euler', 'heun' or 'sra1'. For additive noise the diffusion is
%   state-independent, so Ito and Stratonovich coincide, the Milstein
%   correction vanishes, and NOISE NEVER ENTERS f -- every f evaluation below
%   is the plain deterministic drift, and the stochastic term appears only in
%   the update formula. That is the whole reason a model's dynamics function
%   needs no changes to be simulated as an SDE.
%
%   Strong orders, for additive noise:
%     'euler'  Euler-Maruyama.   1 drift eval/step.  Strong order 1.0.
%     'heun'   stochastic Heun.  2 drift evals/step. Strong order 1.0, but a
%              ~2x smaller noise error constant than Euler and second-order
%              drift, because the trapezoidal corrector gets the mean of the
%              I(1,0) term right.
%     'sra1'   Roessler SRA1.    2 drift evals/step. Strong order 1.5, weak 2.
%
%   SRA1 is the same cost as Heun and half an order more accurate. It manages
%   that by representing the second stochastic integral I(1,0) = int_0^h W ds
%   EXACTLY rather than approximating it: (dW, I(1,0)) is jointly Gaussian, so
%   it is generated from two independent normals via the Kloeden-Platen
%   identity  I(1,0) = (h/2)*(dW + dZ/sqrt(3)).  Knowing where the Brownian
%   path ENDS does not determine the area UNDER it, which is why dZ exists;
%   it is not a second noise source, just a second number describing one.
%
%   The tableau is Roessler's constructSRA1
%     c0 = [0, 3/4]; A0 = [0 0; 3/4 0]; B0 = [0 0; 3/2 0];
%     alpha = [1/3, 2/3]; beta1 = [1, 0]; beta2 = [-1, 1]
%   specialised to time-independent additive noise, where g() is the same at
%   every node so the beta2 contributions cancel (-1 + 1 = 0). Two checks on
%   that specialisation, both pinned in test_sde_integrators: the drift part is
%   a valid second-order RK (b1+b2 = 1, b2*c2 = 1/2), so at sigma = 0 SRA1 must
%   equal RK2 with node 3/4 exactly; and the dW coefficient is exactly sigma.
%
%   NOISE is a struct, or [] for a purely deterministic run:
%     .xi1, .xi2  unit-variance standard normals, numel(idx) x n_steps.
%                 Increments are formed as dW = sqrt(h)*xi1 at use time, so
%                 the stored numbers are step-size independent -- which is what
%                 lets a convergence study build coarse-step noise by summing
%                 adjacent fine-step columns (exact for Brownian increments).
%     .t0, .fs    absolute time of column 1, and the sampling rate the columns
%                 are laid out on.
%     .sigma      scalar diffusion coefficient (already in raw dy/dt units).
%     .idx        state indices that receive noise.
%
%   ABSOLUTE-TIME INDEXING is what makes Benettin work. Columns are addressed
%   by where tspan sits on the global grid, not from the start of this call:
%
%       i0 = round((tspan(1) - noise.t0) * noise.fs) + 1
%
%   so a re-integrated trajectory segment consumes exactly the increments the
%   original run consumed over the same interval. Both trajectories then see
%   one common noise path, the additive noise cancels in their difference to
%   first order, and the Benettin estimate stays valid at any noise level.
%
%   See also ODE_RK4, SRNNMODEL2.

    if numel(tspan) < 3
        error('sde_fixed_step:SpanTooShort', ...
            ['sde_fixed_step requires a full time grid (numel(tspan) >= 3); ' ...
             'got %d point(s). A 2-point [t0 tf] span is not supported by a ' ...
             'fixed-step integrator.'], numel(tspan));
    end

    t  = tspan(:);
    nt = numel(t);
    y  = y0(:);
    Y  = zeros(nt, numel(y));
    Y(1, :) = y.';

    switch lower(scheme)
        case 'euler', scheme_code = 1;
        case 'heun',  scheme_code = 2;
        case 'sra1',  scheme_code = 3;
        otherwise
            error('sde_fixed_step:UnknownScheme', ...
                'Unknown scheme ''%s''. Valid: euler, heun, sra1.', ...
                char(string(scheme)));
    end

    has_noise = ~isempty(noise) && noise.sigma ~= 0;
    i0 = 1; idx = []; s = 0;
    if has_noise
        sde_fixed_step_check_noise(noise, t);
        idx = noise.idx;
        s   = noise.sigma;
        i0  = round((t(1) - noise.t0) * noise.fs) + 1;
        last = i0 + nt - 2;
        if i0 < 1 || last > size(noise.xi1, 2)
            error('sde_fixed_step:NoiseOutOfRange', ...
                ['tspan [%g, %g] needs noise columns %d..%d, but the tensor ' ...
                 'starts at t0 = %g and holds %d columns. The requested span ' ...
                 'is not covered by the pre-generated noise.'], ...
                t(1), t(end), i0, last, noise.t0, size(noise.xi1, 2));
        end
    end

    dW = 0; dZ = 0;
    for k = 1:nt-1
        h  = t(k+1) - t(k);
        tk = t(k);
        if has_noise
            j  = i0 + k - 1;
            rh = sqrt(h);
            dW = (s * rh) * noise.xi1(:, j);
            dZ = (s * rh) * noise.xi2(:, j);
        end

        switch scheme_code
            case 1      % Euler-Maruyama
                k1 = odefun(tk, y);
                y  = y + h * k1;
                if has_noise, y(idx) = y(idx) + dW; end

            case 2      % stochastic Heun (trapezoidal predictor-corrector)
                k1 = odefun(tk, y);
                yp = y + h * k1;
                if has_noise, yp(idx) = yp(idx) + dW; end
                k2 = odefun(tk + h, yp);
                y  = y + (h / 2) * (k1 + k2);
                if has_noise, y(idx) = y(idx) + dW; end

            case 3      % Roessler SRA1
                k1 = odefun(tk, y);
                ys = y + (3 * h / 4) * k1;
                if has_noise
                    % (3/2)*sigma*I(1,0)/h with I(1,0)/h = (dW + dZ/sqrt(3))/2
                    ys(idx) = ys(idx) + 0.75 * (dW + dZ / sqrt(3));
                end
                k2 = odefun(tk + 3 * h / 4, ys);
                y  = y + (h / 3) * (k1 + 2 * k2);
                if has_noise, y(idx) = y(idx) + dW; end
        end

        Y(k+1, :) = y.';
    end
end

function sde_fixed_step_check_noise(noise, t)
% Validate the noise struct and that the grid matches the one it was laid out
% on. A mismatched grid would silently pair increments with the wrong steps,
% which is the failure mode most likely to look like "the model is just noisy".
    required = {'xi1', 'xi2', 't0', 'fs', 'sigma', 'idx'};
    missing = required(~isfield(noise, required));
    if ~isempty(missing)
        error('sde_fixed_step:BadNoiseStruct', ...
            'noise is missing field(s): %s.', strjoin(missing, ', '));
    end
    if ~isequal(size(noise.xi1), size(noise.xi2))
        error('sde_fixed_step:BadNoiseStruct', ...
            'noise.xi1 and noise.xi2 must be the same size ([%s] vs [%s]).', ...
            num2str(size(noise.xi1)), num2str(size(noise.xi2)));
    end
    if size(noise.xi1, 1) ~= numel(noise.idx)
        error('sde_fixed_step:BadNoiseStruct', ...
            ['noise.xi1 has %d rows but noise.idx names %d states; the tensor ' ...
             'must carry one row per noise-driven state.'], ...
            size(noise.xi1, 1), numel(noise.idx));
    end

    dt_nominal = 1 / noise.fs;
    dt_actual = diff(t);
    if max(abs(dt_actual - dt_nominal)) > 1e-9 * dt_nominal
        error('sde_fixed_step:GridMismatch', ...
            ['tspan is not on the 1/fs = %g s grid the noise was generated ' ...
             'on (max deviation %.3g s). Increments would be paired with the ' ...
             'wrong steps.'], dt_nominal, max(abs(dt_actual - dt_nominal)));
    end
end
