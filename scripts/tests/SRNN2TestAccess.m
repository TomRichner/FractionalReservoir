classdef SRNN2TestAccess < SRNNModel2
    % SRNN2TestAccess Test-only subclass that exposes the protected static
    % RHS (dynamics_fast) so verification scripts can finite-difference it
    % against the public analytic Jacobian (compute_Jacobian_fast).
    methods (Static)
        function dS = rhs(t, S, params)
            dS = SRNNModel2.dynamics_fast(t, S, params);
        end
    end
end
