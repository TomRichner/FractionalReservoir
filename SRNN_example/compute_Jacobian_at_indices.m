function J_array = compute_Jacobian_at_indices(S_out, J_times, params)
% COMPUTE_JACOBIAN_AT_INDICES computes Jacobian matrices at multiple time indices
%
% J_array = compute_Jacobian_at_indices(S_out, J_times, params)
%
% Evaluates the Jacobian matrix of the SRNN system at specified time indices
% from a trajectory.
%
% Inputs:
%   S_out   - State trajectory matrix (nt × N_sys_eqs) from ODE solver
%   J_times - Vector of time indices to evaluate Jacobian at (e.g., [100, 500, 1000])
%   params  - Parameter structure (same as used in compute_Jacobian and SRNN_reservoir)
%
% Output:
%   J_array - 3D array (N_sys_eqs × N_sys_eqs × length(J_times))
%             J_array(:,:,i) is the Jacobian at time index J_times(i)
%
% Example:
%   % After running ODE integration
%   J_times = [100, 500, 1000];
%   J_array = compute_Jacobian_at_indices(S_out, J_times, params);
%   % Access Jacobian at first time index:
%   J_at_t1 = J_array(:,:,1);

    % Get dimensions
    N_sys_eqs = size(S_out, 2);
    n_times = length(J_times);
    
    % Validate time indices
    nt = size(S_out, 1);
    if any(J_times < 1) || any(J_times > nt)
        error('J_times contains invalid indices. Must be between 1 and %d', nt);
    end
    
    % Initialize output array
    J_array = zeros(N_sys_eqs, N_sys_eqs, n_times);
    
    % Loop over time indices and compute Jacobian at each
    for i = 1:n_times
        % Extract state vector at this time index (convert to column vector)
        S = S_out(J_times(i), :)';
        
        % Compute Jacobian at this state (sparse fast version -> dense)
        J_array(:,:,i) = full(compute_Jacobian_fast(S, params));
    end
    
end

