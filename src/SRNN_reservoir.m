function [dS_dt] = SRNN_reservoir(t, S, t_ex, u_ex, params)
% SRNN_reservoir implements a rate network with spike-frequency adaptation
% and short-term synaptic depression
%
% Implements the following equations:
%   dx_i/dt = (-x_i + sum_j(w_ij * r_j) + u_i) / tau_d
%   r_i = b_i * phi(x_i - c * sum_k(a_i,k))
%   da_i,k/dt = (-a_i,k + r_i) / tau_k
%   db_i/dt = (1 - b_i) / tau_rec - (b_i * r_i) / tau_rel
%
% where c = c_E for excitatory neurons and c = c_I for inhibitory neurons
%
% State organization: S = [a_E(:); a_I(:); b_E(:); b_I(:); x(:)]

    persistent u_interpolant t_ex_last;

    % To improve performance, create a griddedInterpolant for the external
    % input u_ex and store it in a persistent variable. This avoids
    % repeatedly setting up the interpolation on every function call.
    % The interpolant is rebuilt only if the time vector t_ex appears
    % to have changed between simulations.
    if isempty(u_interpolant) || isempty(t_ex_last) || ...
       numel(t_ex_last) ~= numel(t_ex) || t_ex_last(1) ~= t_ex(1) || t_ex_last(end) ~= t_ex(end)
        
        % We use 'none' for extrapolation to match the behavior of the
        % previous interp1qr implementation, which returns NaN for
        % out-of-bounds queries. This can help catch errors if the
        % ODE solver attempts to step outside the defined time range of u_ex.
        u_interpolant = griddedInterpolant(t_ex, u_ex', 'linear', 'none');
        t_ex_last = t_ex;
    end

    %% interpolate u vector
    u = u_interpolant(t)'; % u_interpolant(t) is 1-by-n, so we transpose to n x 1

    %% load parameters
    n = params.n; % total number of neurons
    n_E = params.n_E; % number of excitatory neurons
    n_I = params.n_I; % number of inhibitory neurons
    E_indices = params.E_indices; % indices of E neurons
    I_indices = params.I_indices; % indices of I neurons
    
    n_a_E = params.n_a_E; % number of adaptation time constants for E neurons
    n_a_I = params.n_a_I; % number of adaptation time constants for I neurons
    n_b_E = params.n_b_E; % number of STD timescales for E neurons (0 or 1)
    n_b_I = params.n_b_I; % number of STD timescales for I neurons (0 or 1)
    
    W = params.W; % connection matrix (n x n)
    tau_d = params.tau_d; % dendritic time constant (scalar)
    tau_a_E = params.tau_a_E; % adaptation time constants for E neurons (1 x n_a_E)
    tau_a_I = params.tau_a_I; % adaptation time constants for I neurons (1 x n_a_I)
    tau_b_E_rec = params.tau_b_E_rec; % STD recovery time constant for E neurons (scalar)
    tau_b_E_rel = params.tau_b_E_rel; % STD release time constant for E neurons (scalar)
    tau_b_I_rec = params.tau_b_I_rec; % STD recovery time constant for I neurons (scalar)
    tau_b_I_rel = params.tau_b_I_rel; % STD release time constant for I neurons (scalar)
    
    % Adaptation scaling parameters
    if isfield(params, 'c_E')
        c_E = params.c_E; % adaptation scaling for E neurons (scalar)
    else
        c_E = 1.0; % default to 1.0
    end
    
    if isfield(params, 'c_I')
        c_I = params.c_I; % adaptation scaling for I neurons (scalar)
    else
        c_I = 1.0; % default to 1.0
    end
    
    % Activation function (nonlinearity)
    if isfield(params, 'activation_function') && isa(params.activation_function, 'function_handle')
        activation_function = params.activation_function;
    else
        error('SRNN_reservoir:MissingActivationFunction', ...
              'params.activation_function must be provided as a function handle');
    end

    %% unpack state variables
    % State organization: S = [a_E(:); a_I(:); b_E(:); b_I(:); x(:)]
    % S is N_sys_eqs x 1 here.
    current_idx = 0;

    % --- Adaptation states for E neurons (a_E) ---
    len_a_E = n_E * n_a_E;
    if len_a_E > 0
        a_E = reshape(S(current_idx + (1:len_a_E)), n_E, n_a_E);
    else
        a_E = [];
    end
    current_idx = current_idx + len_a_E;

    % --- Adaptation states for I neurons (a_I) ---
    len_a_I = n_I * n_a_I;
    if len_a_I > 0
        a_I = reshape(S(current_idx + (1:len_a_I)), n_I, n_a_I);
    else
        a_I = [];
    end
    current_idx = current_idx + len_a_I;

    % --- STD states for E neurons (b_E) ---
    len_b_E = n_E * n_b_E;
    if len_b_E > 0
        b_E = S(current_idx + (1:len_b_E));
    else
        b_E = [];
    end
    current_idx = current_idx + len_b_E;

    % --- STD states for I neurons (b_I) ---
    len_b_I = n_I * n_b_I;
    if len_b_I > 0
        b_I = S(current_idx + (1:len_b_I));
    else
        b_I = [];
    end
    current_idx = current_idx + len_b_I;

    % --- Dendritic states (x) ---
    x = S(current_idx + (1:n));

    %% compute firing rates
    x_eff = x; % n x 1, effective dendritic potential before activation function

    % Apply adaptation effect to E neurons (scaled by c_E)
    if n_E > 0 && n_a_E > 0 && ~isempty(a_E)
        % sum(a_E, 2) is n_E x 1, summing across all adaptation variables
        x_eff(E_indices) = x_eff(E_indices) - c_E * sum(a_E, 2);
    end
    
    % Apply adaptation effect to I neurons (scaled by c_I)
    if n_I > 0 && n_a_I > 0 && ~isempty(a_I)
        % sum(a_I, 2) is n_I x 1, summing across all adaptation variables
        x_eff(I_indices) = x_eff(I_indices) - c_I * sum(a_I, 2);
    end
    
    % Apply STD effect (b multiplicative factor)
    b = ones(n, 1);  % Initialize b = 1 for all neurons (no depression)
    if n_b_E > 0 && ~isempty(b_E)
        b(E_indices) = b_E;
    end
    if n_b_I > 0 && ~isempty(b_I)
        b(I_indices) = b_I;
    end
    
    r = b .* activation_function(x_eff); % n x 1, firing rate

    %% compute derivatives
    % dx/dt = -x/tau_d + W*r + u
    dx_dt = (-x + W * r + u) / tau_d;

    % da_E/dt = (r_E - a_E) / tau_a_E
    da_E_dt = [];
    if n_E > 0 && n_a_E > 0 && ~isempty(a_E)
        % r(E_indices) is n_E x 1. tau_a_E is 1 x n_a_E.
        % Broadcasting makes (r_E - a_E) ./ tau_a_E valid (n_E x n_a_E).
        da_E_dt = (r(E_indices) - a_E) ./ tau_a_E;
    end

    % da_I/dt = (r_I - a_I) / tau_a_I
    da_I_dt = [];
    if n_I > 0 && n_a_I > 0 && ~isempty(a_I)
        % r(I_indices) is n_I x 1. tau_a_I is 1 x n_a_I.
        % Broadcasting makes (r_I - a_I) ./ tau_a_I valid (n_I x n_a_I).
        da_I_dt = (r(I_indices) - a_I) ./ tau_a_I;
    end

    % db_E/dt = (1 - b_E) / tau_b_E_rec - (b_E .* r_E) / tau_b_E_rel
    db_E_dt = [];
    if n_E > 0 && n_b_E > 0 && ~isempty(b_E)
        db_E_dt = (1 - b_E) / tau_b_E_rec - (b_E .* r(E_indices)) / tau_b_E_rel;
    end

    % db_I/dt = (1 - b_I) / tau_b_I_rec - (b_I .* r_I) / tau_b_I_rel
    db_I_dt = [];
    if n_I > 0 && n_b_I > 0 && ~isempty(b_I)
        db_I_dt = (1 - b_I) / tau_b_I_rec - (b_I .* r(I_indices)) / tau_b_I_rel;
    end

    %% pack derivatives into output vector
    % State organization: S = [a_E(:); a_I(:); b_E(:); b_I(:); x(:)]
    dS_dt = [da_E_dt(:); da_I_dt(:); db_E_dt(:); db_I_dt(:); dx_dt];

end

