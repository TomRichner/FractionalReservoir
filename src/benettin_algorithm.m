function [LLE, local_lya, finite_lya, t_lya] = benettin_algorithm(X, t, dt, fs, d0, T, lya_dt, params, ode_options, dynamics_func, t_ex, u_ex, ode_solver)
    % Benettin's algorithm to compute the largest Lyapunov exponent
    % reshoots small segments to compute the divergence rate along the system trajectory in X

    % Input Validations
    if ~isscalar(lya_dt) || ~isnumeric(lya_dt) || lya_dt <= 0
        error('benettin_algorithm:InvalidLyaDt', 'lya_dt must be a positive scalar.');
    end

    deci_lya = round(lya_dt*fs);     % samples per Lyapunov interval
    if deci_lya < 1
        error('benettin_algorithm:InvalidDeciLya', 'lya_dt * fs must result in at least 1 sample per Lyapunov interval (deci_lya >= 1). Check lya_dt and fs values.');
    end

    tau_lya = dt*deci_lya;    % Lyapunov rescaling time interval (integration time between rescalings)
    t_lya    = t(1:deci_lya:end);     % direct decimation of `t`

    if t_lya(end) + tau_lya > T(2)    % keep segments fully inside [T(1),T(2)]
        t_lya(end) = [];
    end
    nt_lya  = numel(t_lya);           % number of Lyapunov intervals

    local_lya  = zeros(nt_lya,1);
    finite_lya = nan(nt_lya,1);
    sum_log_stretching_factors = 0;

    % Initial perturbation
    n_state = size(X,2);
    rnd_IC = randn(n_state,1);
    pert = (rnd_IC./norm(rnd_IC)).*d0;

    % Get bounds for the full state vector X once, as they are static.
    % This now expects params to contain all necessary info for get_minMaxRange
    min_max_range = get_minMaxRange(params); 
    % Ensure current_min_max_range matches the structure of X_k_pert (n_state x 2)
    
    min_bnds = min_max_range(:, 1); % n_state x 1 column vector
    max_bnds = min_max_range(:, 2); % n_state x 1 column vector

    for k = 1:nt_lya
        idx_start = (k-1)*deci_lya + 1;      % index of t_lya(k) in `t`
        idx_end   = idx_start + deci_lya;    % => t_lya(k) + tau_lya

        X_start = X(idx_start,:).';      % fiducial state at t_lya(k), n_state x 1

        X_k_pert = X_start + pert; % Rescale perturbation from previous normalized delta, n_state x 1

        % Vectorized clipping of X_k_pert to ensure it's within bounds
        % Apply lower bounds where they are defined and violated
        idx_violates_min = ~isnan(min_bnds) & (X_k_pert < min_bnds);
        X_k_pert(idx_violates_min) = min_bnds(idx_violates_min);

        % Apply upper bounds where they are defined and violated
        idx_violates_max = ~isnan(max_bnds) & (X_k_pert > max_bnds);
        X_k_pert(idx_violates_max) = max_bnds(idx_violates_max);
        
        % Integrate ONLY the perturbed trajectory over [t_lya(k), t_lya(k)+tau_lya]
        t_seg_detailed = t_lya(k) + (0:dt:tau_lya); % tau lya is always a multiple of dt, so this should be fine
        
        % Pass the now-bounded X_k_pert and the detailed t_seg_detailed to the ODE solver
        % ode_options are passed along; for ode45 they dictate MaxStep=dt, 
        % for ode4/ode_RKn they are mostly ignored for step control but passed to odefun if needed.
        [~, X_pert_output_all_steps] = ode_solver(@(tt,XX) dynamics_func(tt,XX,t_ex,u_ex,params), t_seg_detailed, X_k_pert, ode_options);
        
        X_pert_end = X_pert_output_all_steps(end,:).';
        
        X_end   = X(idx_end,:).';        % fiducial state tau_lya later

        % Local exponent
        delta   = X_pert_end - X_end;
        d_k     = norm(delta);
        local_lya(k) = log(d_k/d0)/tau_lya;

        % Check for divergence (NaN or Inf in local Lyapunov exponent)
        if ~isfinite(local_lya(k))
            warning('benettin_algorithm:Divergence', 'System diverged or produced non-finite Lyapunov exponent at t=%f. Truncating results.', t_lya(k));
            
            if k > 1
                % Use the last valid finite Lyapunov exponent as the final LLE
                last_valid_finite_lya = finite_lya(1:k-1);
                last_valid_finite_lya = last_valid_finite_lya(~isnan(last_valid_finite_lya));
                if ~isempty(last_valid_finite_lya)
                    LLE = last_valid_finite_lya(end);
                else
                    LLE = 0; % No valid finite LLE computed before divergence
                end
            else
                LLE = 0; % Diverged on the first step
            end
            
            % Truncate output arrays to remove invalid data
            local_lya(k:end) = [];
            finite_lya(k:end) = [];
            t_lya(k:end) = [];
            
            return; % Stop computation and return partial results
        end

        pert = (delta./d_k).*d0; % rescaled perturbation for next step

        % Accumulate stretching only from t >= 0
        if t_lya(k) >= 0
            sum_log_stretching_factors = sum_log_stretching_factors + log(d_k/d0);
            finite_lya(k,1) = sum_log_stretching_factors / max(t_lya(k)+tau_lya, eps);
        end
    end

    % Final LLE is the last computed finite Lyapunov exponent.
    last_valid_finite_lya = finite_lya(~isnan(finite_lya));
    if ~isempty(last_valid_finite_lya)
        LLE = last_valid_finite_lya(end);
    else
        LLE = 0; % No valid finite LLE was computed
    end
end