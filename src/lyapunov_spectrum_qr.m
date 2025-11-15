function [LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya_vec] = lyapunov_spectrum_qr(X_fid_traj, t_fid_traj, lya_dt_interval, params, ode_solver, ode_options_main, jacobian_func_handle, T_full_interval, N_states_sys, fs_fid)
    % Calculates the Lyapunov spectrum using the QR decomposition method.
    % Integrates variational equations along a pre-computed fiducial trajectory.
    %
    % Inputs:
    %   X_fid_traj        : Fiducial trajectory (N_samples x N_states_sys)
    %   t_fid_traj        : Time vector for fiducial trajectory
    %   lya_dt_interval   : Rescaling time interval (tau_lya)
    %   params            : System parameters
    %   ode_options_main  : ODE solver options for variational equations
    %   jacobian_func_handle : Handle to the system's Jacobian function
    %   T_full_interval   : Original time interval [T_start, T_end] for global LE calculation
    %   N_states_sys      : Dimension of the system
    %   fs_fid            : Sampling frequency of the fiducial trajectory
    %
    % Outputs:
    %   LE_spectrum           : Vector of global Lyapunov exponents (N_states_sys x 1)
    %   local_LE_spectrum_t   : Matrix of local LEs over time (N_timesteps_lya x N_states_sys)
    %   finite_LE_spectrum_t  : Matrix of finite-time LEs over time (N_timesteps_lya x N_states_sys)
    %   t_lya_vec             : Time vector for Lyapunov exponent estimates

    % Create an interpolant for the fiducial trajectory for use in ODE solver
    fiducial_interpolants = cell(N_states_sys, 1);
    for i = 1:N_states_sys
        fiducial_interpolants{i} = griddedInterpolant(t_fid_traj, X_fid_traj(:,i), 'pchip');
    end

    dt_fid = 1/fs_fid;
    deci_lya = round(lya_dt_interval / dt_fid);
    if deci_lya == 0
        error('lya_dt_interval is too small compared to the fiducial trajectory sampling time, leading to zero samples per interval.');
    end
    tau_lya = dt_fid * deci_lya; 

    t_lya_indices = 1:deci_lya:length(t_fid_traj);
    t_lya_vec = t_fid_traj(t_lya_indices);

    if ~isempty(t_lya_vec) && (t_lya_vec(end) + tau_lya > t_fid_traj(end) + eps(t_fid_traj(end)))
        t_lya_vec(end) = [];
        t_lya_indices(end) = []; 
    end
    
    nt_lya = numel(t_lya_vec);
    if nt_lya == 0
        warning('No Lyapunov intervals could be formed. Check T, lya_dt_interval, and fs_fid.');
        LE_spectrum = nan(N_states_sys,1); local_LE_spectrum_t = []; finite_LE_spectrum_t = [];
        return;
    end

    Q_current = eye(N_states_sys); 
    sum_log_R_diag = zeros(N_states_sys, 1); 
    
    local_LE_spectrum_t  = zeros(nt_lya, N_states_sys);
    finite_LE_spectrum_t = nan(nt_lya, N_states_sys);
    
    total_positive_time_accumulated = 0; 
    ode_options_var = ode_options_main;

    for k = 1:nt_lya
        t_start_segment = t_lya_vec(k);
        t_end_segment = t_start_segment + tau_lya;
        t_end_segment = min(t_end_segment, t_fid_traj(end));
        
        current_segment_duration = t_end_segment - t_start_segment;
        if current_segment_duration <= eps 
            if k > 1
                local_LE_spectrum_t(k,:) = local_LE_spectrum_t(k-1,:);
                finite_LE_spectrum_t(k,:) = finite_LE_spectrum_t(k-1,:);
            else
                local_LE_spectrum_t(k,:) = NaN;
                finite_LE_spectrum_t(k,:) = NaN;
            end
            continue;
        end
        
        t_span_ode = [t_start_segment, t_end_segment];
        Psi0_vec = reshape(Q_current, [], 1);
        
        % [~, Psi_solution_vec] = ode15s(@variational_eqs_ode, t_span_ode, Psi0_vec, ode_options_var);
        [~, Psi_solution_vec] = ode_solver(@variational_eqs_ode, t_span_ode, Psi0_vec, ode_options_var);
        
        Psi_evolved_matrix = reshape(Psi_solution_vec(end,:)', [N_states_sys, N_states_sys]);
        
        % Check for divergence in the integrated variational equations
        if any(~isfinite(Psi_evolved_matrix(:)))
            warning('lyapunov_spectrum_qr:Divergence', 'System diverged. Variational equation integration produced non-finite values at t=%f. Truncating results.', t_start_segment);
            
            % Calculate final spectrum based on data accumulated so far
            if total_positive_time_accumulated > eps
                LE_spectrum = sum_log_R_diag / total_positive_time_accumulated;
            else
                LE_spectrum = nan(N_states_sys,1);
            end

            % Truncate output arrays to exclude invalid data
            t_lya_vec(k:end) = [];
            local_LE_spectrum_t(k:end,:) = [];
            finite_LE_spectrum_t(k:end,:) = [];

            return; % Stop computation and return partial results
        end
        
        [Q_new, R_segment] = qr(Psi_evolved_matrix);
        diag_R = diag(R_segment);
        
        log_abs_diag_R = log(abs(diag_R)); 
        valid_diag_R = abs(diag_R) > eps; 
        
        current_local_LEs = zeros(N_states_sys,1);
        current_local_LEs(valid_diag_R) = log_abs_diag_R(valid_diag_R) / current_segment_duration;
        current_local_LEs(~valid_diag_R) = -Inf; 
        local_LE_spectrum_t(k,:) = current_local_LEs';
        
        if t_start_segment >= -eps(0) 
            sum_log_R_diag(valid_diag_R) = sum_log_R_diag(valid_diag_R) + log_abs_diag_R(valid_diag_R);
            total_positive_time_accumulated = total_positive_time_accumulated + current_segment_duration;
        end
        
        if total_positive_time_accumulated > eps
            finite_LE_spectrum_t(k,:) = (sum_log_R_diag / total_positive_time_accumulated)';
        elseif k > 1 && t_start_segment >= -eps(0)
             finite_LE_spectrum_t(k,:) = finite_LE_spectrum_t(k-1,:); 
        else
            finite_LE_spectrum_t(k,:) = NaN; 
        end

        Q_current = Q_new;
    end

    if total_positive_time_accumulated > eps
        LE_spectrum = sum_log_R_diag / total_positive_time_accumulated;
    else
        warning('QR:NoPositiveTime', 'No accumulation over positive time for global LEs.');
        LE_spectrum = nan(N_states_sys,1);
    end
    
    % Nested function definition moved to the end of the parent function
    function dPsi_vec_dt = variational_eqs_ode(tt, current_Psi_vec)
        % This nested function defines the variational ODE system:
        % d(Psi)/dt = J(X_fid(t)) * Psi
        % It has access to variables from the parent function's workspace,
        % such as fiducial_interpolants, N_states_sys, jacobian_func_handle, and params.

        % Interpolate fiducial state X_fid at current time tt
        X_fid_at_tt = zeros(N_states_sys, 1);
        for state_idx_loop = 1:N_states_sys % Renamed loop variable to avoid conflict if N_states was used
            X_fid_at_tt(state_idx_loop) = fiducial_interpolants{state_idx_loop}(tt);
        end
        
        % Calculate Jacobian at X_fid(tt) using the provided function handle
        J_matrix = jacobian_func_handle(tt, X_fid_at_tt, params);
        
        % Reshape Psi_vec (input from ODE solver) to matrix form
        Psi_matrix = reshape(current_Psi_vec, [N_states_sys, N_states_sys]);
        
        % Calculate d(Psi_matrix)/dt = J * Psi
        dPsi_matrix_dt = J_matrix * Psi_matrix;
        
        % Reshape back to vector for ODE solver output
        dPsi_vec_dt = reshape(dPsi_matrix_dt, [], 1);
    end % End of nested function variational_eqs_ode

end % End of lyapunov_spectrum_qr function