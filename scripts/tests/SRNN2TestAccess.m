classdef SRNN2TestAccess < SRNNModel2
    % SRNN2TestAccess Test-only subclass that exposes protected statics so
    % verification scripts can exercise them directly:
    %   rhs      - dynamics_fast, to finite-difference against the public
    %              analytic Jacobian (compute_Jacobian_fast)
    %   lya_grid - lyapunov_sample_grid, the segment layout shared by the
    %              Benettin and QR paths (see test_lya_window)
    methods
        function build_noise_for_test(obj)
            % Expose the protected build_noise so a test can inspect the drawn
            % increments without running a whole simulation.
            obj.build_noise();
        end
    end

    methods (Static)
        function dS = rhs(t, S, params)
            dS = SRNNModel2.dynamics_fast(t, S, params);
        end

        function [idx, t_lya, acc_start] = lya_grid(t, dt, deci, tau, T, warmup)
            [idx, t_lya, acc_start] = SRNNModel2.lyapunov_sample_grid( ...
                t, dt, deci, tau, T, warmup, 'SRNNModel');
        end
    end
end
