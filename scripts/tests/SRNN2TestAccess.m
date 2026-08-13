classdef SRNN2TestAccess < SRNNModel2
    % SRNN2TestAccess Test-only subclass that exposes protected statics so
    % verification scripts can exercise them directly:
    %   rhs      - dynamics_fast, to finite-difference against the public
    %              analytic Jacobian (compute_Jacobian_fast)
    %   lya_grid - lyapunov_sample_grid, the segment layout shared by the
    %              Benettin and QR paths (see test_lya_window)
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
