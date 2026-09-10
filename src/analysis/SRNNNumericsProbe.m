classdef SRNNNumericsProbe < SRNNCellTypePairs
    % SRNNNumericsProbe SRNNCellTypePairs with three doors opened for
    % run_numerics_verification. It adds NO physics and NO analysis logic; the
    % reshooting loop lives in the analysis function. The class exists because
    % three things that stage needs are protected on the parent:
    %
    %   * noise_increments is SetAccess = protected, and run() calls
    %     build_noise() unconditionally, so a caller cannot make a run consume a
    %     Brownian path it supplies -- which is what a strong-convergence study
    %     needs (the SAME path, rebuilt on a coarser grid by coarsen_noise).
    %     set_noise() stores a struct and the build_noise override uses it.
    %   * run() clears noise_increments once the trajectory and Lyapunov steps
    %     are done (it is the largest thing on the object). The override keeps
    %     a handle in consumed_noise so the stage can read back the path the
    %     REFERENCE run used and coarsen it.
    %   * The trajectory integrator is reachable via the public integrate(), but
    %     the right-hand side it needs is assembled inside run() from
    %     cached_params and u_interpolant. integrate_segment() assembles it the
    %     same way (see compute_lyapunov) so a segment can be re-integrated from
    %     a supplied state at an absolute time, exactly as Benettin does.
    %
    % A model class is otherwise untouched: ParamSpaceAnalysis2 validates
    % model_defaults against the parent's property list and freezes it into
    % resolved_defaults, so a public knob on SRNNCellTypePairs would show up in
    % every run's provenance for something no sweep sets. Precedent for the
    % pattern: SRNNPairsTestAccess (test-only) and SRNN_ESN_reservoir.
    %
    % See also: run_numerics_verification, coarsen_noise, sde_fixed_step,
    %           SRNNCellTypePairs

    properties (SetAccess = protected)
        % Noise struct (xi1/xi2/t0/fs/sigma/idx, see sde_fixed_step) to use in
        % place of build_noise's own draw. Empty means draw as the parent does.
        injected_noise = []
        % The noise struct the last run() consumed, kept after run() clears
        % noise_increments. Empty for a deterministic run.
        consumed_noise = []
    end

    methods
        function set_noise(obj, nz)
            % SET_NOISE Make the next run() consume this Brownian path.
            %   The struct must lie on THIS model's grid: fs equal to obj.fs and
            %   t0 equal to the start of the time vector, else sde_fixed_step
            %   would pair increments with the wrong steps (it checks the grid,
            %   but only after build() has laid out t_ex, so check here too).
            required = {'xi1', 'xi2', 't0', 'fs', 'sigma', 'idx'};
            missing = required(~isfield(nz, required));
            if ~isempty(missing)
                error('SRNNNumericsProbe:BadNoiseStruct', ...
                    'noise is missing field(s): %s.', strjoin(missing, ', '));
            end
            if abs(nz.fs - obj.fs) > 1e-9 * obj.fs
                error('SRNNNumericsProbe:NoiseGridMismatch', ...
                    'noise.fs = %g but the model runs at fs = %g.', nz.fs, obj.fs);
            end
            if abs(nz.t0 - obj.T_range(1)) > 1e-9
                error('SRNNNumericsProbe:NoiseGridMismatch', ...
                    'noise.t0 = %g but the model starts at T_range(1) = %g.', ...
                    nz.t0, obj.T_range(1));
            end
            obj.injected_noise = nz;
        end

        function clear_noise(obj)
            % CLEAR_NOISE Go back to drawing the path from noise_seed.
            obj.injected_noise = [];
        end

        function build_noise(obj)
            % BUILD_NOISE Use the injected path if there is one, else the
            % parent's draw. Either way keep what run() is about to consume.
            if isempty(obj.injected_noise)
                build_noise@SRNNCellTypePairs(obj);
            else
                if obj.sigma_u_noise == 0
                    error('SRNNNumericsProbe:NoiseWithoutSigma', ...
                        ['A noise path was injected but sigma_u_noise = 0. ' ...
                         'Set sigma_u_noise so the integrator is stochastic.']);
                end
                nz = obj.injected_noise;
                n_steps = numel(obj.t_ex) - 1;
                if size(nz.xi1, 2) < n_steps
                    error('SRNNNumericsProbe:NoiseTooShort', ...
                        ['The injected path holds %d steps but the run needs ' ...
                         '%d (fs = %g over T_range = [%g, %g]).'], ...
                        size(nz.xi1, 2), n_steps, obj.fs, obj.T_range);
                end
                nz.sigma = obj.sigma_x_raw;
                nz.idx   = obj.cached_params.state_layout.x;
                obj.noise_increments = nz;
            end
            obj.consumed_noise = obj.noise_increments;
        end

        function [t_out, S_out] = integrate_segment(obj, tspan, S0)
            % INTEGRATE_SEGMENT Integrate from state S0 over the grid tspan with
            % the model's own integrator, driven by the model's own stimulus.
            %   tspan must be a full grid at 1/fs spacing inside the built time
            %   vector (the stimulus interpolant does not extrapolate). With a
            %   stochastic integrator, noise_increments must be populated: pass
            %   the path via set_noise or call arm_noise() after build().
            if ~obj.is_built
                error('SRNNNumericsProbe:NotBuilt', ...
                    'Call build() before integrate_segment().');
            end
            params = obj.cached_params;
            params.u_interpolant = obj.u_interpolant;
            rhs = @(t, S) SRNNCellTypePairs.dynamics_fast(t, S, params);
            [t_out, S_out] = obj.integrate(rhs, tspan, S0);
        end

        function arm_noise(obj)
            % ARM_NOISE Populate noise_increments without running, so
            % integrate_segment can be used on a model that never called run().
            % Uses the injected path when there is one.
            if ~obj.is_built
                error('SRNNNumericsProbe:NotBuilt', ...
                    'Call build() before arm_noise().');
            end
            obj.build_noise();
        end

        function disarm_noise(obj)
            % DISARM_NOISE Drop the (large) noise tensor once the segments are done.
            obj.noise_increments = [];
        end
    end
end
