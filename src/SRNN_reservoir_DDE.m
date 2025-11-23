function dS_dt = SRNN_reservoir_DDE(t, S, S_delay, t_ex, u_ex, params)
% SRNN_RESERVOIR_DDE Implements the DDE version of the SRNN
%
% Syntax:
%   dS_dt = SRNN_reservoir_DDE(t, S, S_delay, t_ex, u_ex, params)
%
% Inputs:
%   t               - Current time
%   S               - Current state vector
%   S_delay         - Delayed state vectors (column k corresponds to params.lags(k))
%   t_ex    - External input time vector
%   u_ex    - External input matrix
%   params  - Parameter struct
%
% Params must contain:
%   params.W_components - Cell array where W_components{1} is instantaneous connectivity,
%                         and W_components{k+1} corresponds to params.lags(k).

    persistent u_interpolant t_ex_last;

    % Initialize interpolant if needed
    if isempty(u_interpolant) || isempty(t_ex_last) || ...
       numel(t_ex_last) ~= numel(t_ex) || t_ex_last(1) ~= t_ex(1) || t_ex_last(end) ~= t_ex(end)
        u_interpolant = griddedInterpolant(t_ex, u_ex', 'linear', 'none');
        t_ex_last = t_ex;
    end

    %% 1. Interpolate external input
    u = u_interpolant(t)'; 

    %% 2. Compute current firing rates and unpack state
    [r_current, x_current, a_current, b_current] = compute_rates_and_unpack(S, params);

    %% 3. Compute recurrent input (Instantaneous + Delayed)
    
    % Instantaneous contribution (Lag 0)
    W_inst = params.W_components{1};
    input_recurrent = W_inst * r_current;
    
    % Delayed contributions
    % params.lags matches columns of Z
    if isfield(params, 'lags') && ~isempty(params.lags)
        for k = 1:length(params.lags)
            S_delayed_k = S_delay(:, k);
            % We only need r from the delayed state
            [r_delayed, ~, ~, ~] = compute_rates_and_unpack(S_delayed_k, params);
            
            W_delayed = params.W_components{k+1};
            input_recurrent = input_recurrent + W_delayed * r_delayed;
        end
    end

    %% 4. Compute Derivatives
    
    % Parameters
    tau_d = params.tau_d;
    n_E = params.n_E; n_I = params.n_I;
    n_a_E = params.n_a_E; n_a_I = params.n_a_I;
    n_b_E = params.n_b_E; n_b_I = params.n_b_I;
    E_indices = params.E_indices;
    I_indices = params.I_indices;
    
    tau_a_E = params.tau_a_E;
    tau_a_I = params.tau_a_I;
    tau_b_E_rec = params.tau_b_E_rec; tau_b_E_rel = params.tau_b_E_rel;
    tau_b_I_rec = params.tau_b_I_rec; tau_b_I_rel = params.tau_b_I_rel;

    % dx/dt = (-x + Input + u) / tau_d
    dx_dt = (-x_current + input_recurrent + u) / tau_d;

    % Adaptation derivatives (depend on CURRENT firing rate r_current)
    
    % da_E/dt
    da_E_dt = [];
    if n_E > 0 && n_a_E > 0 && ~isempty(a_current{1})
        % a_current is returned as cell {a_E, a_I} to handle unpacking logic cleanly
        % but compute_rates_and_unpack might just return the raw matrices. 
        % Let's stick to the raw matrix unpacking inside the helper or here.
        % Re-using logic from SRNN_reservoir
        
        % Helper returns a_current as {a_E, a_I}
        val_a_E = a_current{1};
        da_E_dt = (r_current(E_indices) - val_a_E) ./ tau_a_E;
    end

    % da_I/dt
    da_I_dt = [];
    if n_I > 0 && n_a_I > 0 && ~isempty(a_current{2})
        val_a_I = a_current{2};
        da_I_dt = (r_current(I_indices) - val_a_I) ./ tau_a_I;
    end

    % STD derivatives (depend on CURRENT firing rate)
    
    % db_E/dt
    db_E_dt = [];
    if n_E > 0 && n_b_E > 0 && ~isempty(b_current{1})
        val_b_E = b_current{1};
        db_E_dt = (1 - val_b_E) / tau_b_E_rec - (val_b_E .* r_current(E_indices)) / tau_b_E_rel;
    end

    % db_I/dt
    db_I_dt = [];
    if n_I > 0 && n_b_I > 0 && ~isempty(b_current{2})
        val_b_I = b_current{2};
        db_I_dt = (1 - val_b_I) / tau_b_I_rec - (val_b_I .* r_current(I_indices)) / tau_b_I_rel;
    end

    %% 5. Pack Derivatives
    dS_dt = [da_E_dt(:); da_I_dt(:); db_E_dt(:); db_I_dt(:); dx_dt];

end

function [r, x, a_cell, b_cell] = compute_rates_and_unpack(S, params)
    % Unpacks state and computes firing rates
    
    n = params.n;
    n_E = params.n_E; n_I = params.n_I;
    E_indices = params.E_indices; I_indices = params.I_indices;
    
    n_a_E = params.n_a_E; n_a_I = params.n_a_I;
    n_b_E = params.n_b_E; n_b_I = params.n_b_I;
    
    c_E = params.c_E; c_I = params.c_I;
    activation_function = params.activation_function;
    
    %% Unpack
    current_idx = 0;
    
    % a_E
    len_a_E = n_E * n_a_E;
    if len_a_E > 0
        a_E = reshape(S(current_idx + (1:len_a_E)), n_E, n_a_E);
    else
        a_E = [];
    end
    current_idx = current_idx + len_a_E;
    
    % a_I
    len_a_I = n_I * n_a_I;
    if len_a_I > 0
        a_I = reshape(S(current_idx + (1:len_a_I)), n_I, n_a_I);
    else
        a_I = [];
    end
    current_idx = current_idx + len_a_I;
    
    % b_E
    len_b_E = n_E * n_b_E;
    if len_b_E > 0
        b_E = S(current_idx + (1:len_b_E));
    else
        b_E = [];
    end
    current_idx = current_idx + len_b_E;
    
    % b_I
    len_b_I = n_I * n_b_I;
    if len_b_I > 0
        b_I = S(current_idx + (1:len_b_I));
    else
        b_I = [];
    end
    current_idx = current_idx + len_b_I;
    
    % x
    x = S(current_idx + (1:n));
    
    %% Compute Rates
    x_eff = x;
    
    % Adaptation
    if ~isempty(a_E)
        x_eff(E_indices) = x_eff(E_indices) - c_E * sum(a_E, 2);
    end
    if ~isempty(a_I)
        x_eff(I_indices) = x_eff(I_indices) - c_I * sum(a_I, 2);
    end
    
    % STD
    b = ones(n, 1);
    if ~isempty(b_E)
        b(E_indices) = b_E;
    end
    if ~isempty(b_I)
        b(I_indices) = b_I;
    end
    
    r = b .* activation_function(x_eff);
    
    % Pack aux outputs for derivative calc
    a_cell = {a_E, a_I};
    b_cell = {b_E, b_I};
end

