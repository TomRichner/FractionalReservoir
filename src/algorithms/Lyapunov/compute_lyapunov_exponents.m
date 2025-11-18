function lya_results = compute_lyapunov_exponents(Lya_method, S_out, t_out, dt, fs, T_interval, params, opts, ode_solver, rhs_func, t_ex, u_ex)
% compute_lyapunov_exponents - Compute Lyapunov exponents using various methods
%
% Syntax:
%   lya_results = compute_lyapunov_exponents(Lya_method, S_out, t_out, dt, fs, T_interval, params, opts, ode_solver, rhs_func, t_ex, u_ex)
%
% Description:
%   Computes Lyapunov exponents using either Benettin's algorithm (single
%   largest exponent) or QR decomposition method (full spectrum). Returns
%   an empty struct if method is 'none'.
%
% Inputs:
%   Lya_method  - 'benettin', 'qr', or 'none'
%   S_out       - State trajectory (nt x N_sys_eqs)
%   t_out       - Time vector (nt x 1)
%   dt          - Integration time step (seconds)
%   fs          - Sampling frequency (Hz)
%   T_interval  - [T_start, T_end] time interval for analysis
%   params      - SRNN parameters struct
%   opts        - ODE solver options
%   ode_solver  - ODE solver function handle (e.g., @ode45)
%   rhs_func    - Right-hand side function handle for integration
%   t_ex        - External input time vector
%   u_ex        - External input matrix
%
% Outputs:
%   lya_results - Struct containing Lyapunov analysis results:
%                 For 'benettin': LLE, local_lya, finite_lya, t_lya
%                 For 'qr': LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t,
%                           t_lya, sort_idx, params.N_sys_eqs
%                 For 'none': empty struct
%
% Example:
%   lya_results = compute_lyapunov_exponents('benettin', S_out, t_out, ...
%       dt, fs, [3, 150], params, opts, @ode45, @SRNN_reservoir, t_ex, u_ex);

    lya_results = struct();
    
    if strcmpi(Lya_method, 'none')
        return;
    end
    
    % Adjust lya_dt based on method: QR needs longer interval
    if strcmpi(Lya_method, 'qr')
        lya_dt = 0.1;  % Longer interval for QR method
    elseif strcmpi(Lya_method, 'benettin')
        lya_dt = 0.02;  % Standard interval for Benettin
    else
        lya_dt = 0.1;
    end
    
    % Check lya_dt is a nice multiple of dt
    if abs(round(lya_dt/dt) - lya_dt/dt) > 1e-11
        error('lya_dt must be a multiple of dt');
    end
    
    lya_fs = 1 / lya_dt;
    
    switch lower(Lya_method)
        case 'benettin'
            fprintf('Computing largest Lyapunov exponent using Benettin''s algorithm...\n');
            d0 = 1e-3;
            tic
            [LLE, local_lya, finite_lya, t_lya] = benettin_algorithm(S_out, t_out, dt, fs, d0, T_interval, lya_dt, params, opts, rhs_func, t_ex, u_ex, ode_solver);
            toc
            fprintf('Largest Lyapunov Exponent: %.4f\n', LLE);
            lya_results.LLE = LLE;
            lya_results.local_lya = local_lya;
            lya_results.finite_lya = finite_lya;
            lya_results.t_lya = t_lya;
            lya_results.lya_dt = lya_dt;
            lya_results.lya_fs = lya_fs;
            
        case 'qr'
            fprintf('Computing full Lyapunov spectrum using QR decomposition method...\n');
            tic
            [LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya] = lyapunov_spectrum_qr(S_out, t_out, lya_dt, params, ode_solver, opts, @SRNN_Jacobian_wrapper, T_interval, params.N_sys_eqs, fs);
            toc
            fprintf('Lyapunov Dimension: %.2f\n', compute_kaplan_yorke_dimension(LE_spectrum));
            lya_results.LE_spectrum = LE_spectrum;
            lya_results.local_LE_spectrum_t = local_LE_spectrum_t;
            lya_results.finite_LE_spectrum_t = finite_LE_spectrum_t;
            lya_results.t_lya = t_lya;
            lya_results.params.N_sys_eqs = params.N_sys_eqs;
            
            % Sort LE spectra by real component (descending) and keep index map
            [sorted_LE, sort_idx] = sort(real(lya_results.LE_spectrum), 'descend');
            lya_results.LE_spectrum = sorted_LE;
            lya_results.local_LE_spectrum_t = lya_results.local_LE_spectrum_t(:, sort_idx);
            lya_results.finite_LE_spectrum_t = lya_results.finite_LE_spectrum_t(:, sort_idx);
            lya_results.sort_idx = sort_idx;
            lya_results.lya_dt = lya_dt;
            lya_results.lya_fs = lya_fs;
            fprintf('Largest Lyapunov Exponent (sorted): %.4f\n', lya_results.LE_spectrum(1));
            
        otherwise
            error('Unknown Lyapunov method: %s. Use ''benettin'', ''qr'', or ''none''.', Lya_method);
    end
end

%% Helper function to compute Kaplan-Yorke dimension
function D_KY = compute_kaplan_yorke_dimension(lambda)
    % Compute Kaplan-Yorke (Lyapunov) dimension from spectrum
    % lambda: sorted Lyapunov exponents (descending order)
    
    lambda = sort(lambda, 'descend');
    cumsum_lambda = cumsum(lambda);
    
    % Find largest j such that sum of first j exponents is non-negative
    j = find(cumsum_lambda >= 0, 1, 'last');
    
    if isempty(j)
        D_KY = 0;
    elseif j == length(lambda)
        D_KY = length(lambda);
    else
        D_KY = j + cumsum_lambda(j) / abs(lambda(j+1));
    end
end

%% Jacobian wrapper for lyapunov_spectrum_qr
function J_jac = SRNN_Jacobian_wrapper(tt, S, params)
    J_jac = compute_Jacobian_fast(S, params);
end

