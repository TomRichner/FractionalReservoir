function [dX_dt] = SRNN_NL(t, X, t_ex, u_ex, params)
% SRNN_NL supports nonlinear activation functions other than just ReLU

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
    % u = interp1(t_ex, u_ex', t', 'linear')'; % u will be n x 1
    % u = interp1qr(t_ex, u_ex',t')'; % faster version from the file exchange
    u = u_interpolant(t)'; % u_interpolant(t) is 1-by-n, so we transpose.

    %% load parameters
    n = params.n; % total number of neurons
    M = params.M; % connection matrix
    tau_d = params.tau_d; % s, scalar, dendritic time constant
    
    if isfield(params, 'activation_function') && isa(params.activation_function, 'function_handle')
        activation_function = params.activation_function;
    else
        activation_function = @(x) max(0, x);
    end
    
    % E/I specific parameters
    n_E = params.n_E;
    n_I = params.n_I;
    E_indices = params.E_indices;
    I_indices = params.I_indices;

    n_a_E = params.n_a_E; 
    n_a_I = params.n_a_I;
    n_b_E = params.n_b_E;
    n_b_I = params.n_b_I;

    tau_a_E = params.tau_a_E; % 1 x n_a_E (or empty)
    tau_a_I = params.tau_a_I; % 1 x n_a_I (or empty)
    tau_b_E = params.tau_b_E; % 1 x n_b_E (or empty)
    tau_b_I = params.tau_b_I; % 1 x n_b_I (or empty)
    
    c_SFA = params.c_SFA; % n x 1, SFA strength coefficients
    F_STD = params.F_STD; % n x 1, STD strength coefficients
    tau_STD = params.tau_STD; % scalar, STD recovery time const for b's derivative

    %% unpack state variables by inlining unpack_SRNN_state.m for performance

    % unpack_SRNN_state will return:
    % a_E: n_E x n_a_E (or empty)
    % a_I: n_I x n_a_I (or empty)
    % b_E: n_E x n_b_E (or empty)
    % b_I: n_I x n_b_I (or empty)
    % u_d: n x 1
    % [a_E, a_I, b_E, b_I, u_d] = unpack_SRNN_state(X, params); % previously unpacked with a function, but this is slower than inline below

    % X is N_sys_eqs x 1 here.
    current_idx = 0;

    % --- SFA states for E neurons (a_E) ---
    len_a_E = n_E * n_a_E;
    if len_a_E > 0
        a_E = reshape(X(current_idx + (1:len_a_E)), n_E, n_a_E);
    else
        a_E = [];
    end
    current_idx = current_idx + len_a_E;

    % --- SFA states for I neurons (a_I) ---
    len_a_I = n_I * n_a_I;
    if len_a_I > 0
        a_I = reshape(X(current_idx + (1:len_a_I)), n_I, n_a_I);
    else
        a_I = [];
    end
    current_idx = current_idx + len_a_I;

    % --- STD states for E neurons (b_E) ---
    len_b_E = n_E * n_b_E;
    if len_b_E > 0
        b_E = reshape(X(current_idx + (1:len_b_E)), n_E, n_b_E);
    else
        b_E = [];
    end
    current_idx = current_idx + len_b_E;

    % --- STD states for I neurons (b_I) ---
    len_b_I = n_I * n_b_I;
    if len_b_I > 0
        b_I = reshape(X(current_idx + (1:len_b_I)), n_I, n_b_I);
    else
        b_I = [];
    end
    current_idx = current_idx + len_b_I;

    % --- Dendrite states (u_d) ---
    u_d = X(current_idx + (1:n));

    %% make dependent variables from state variables and parameters
    u_eff = u_d; % n x 1, effective dendritic potential before relu

    % Apply SFA effect to E neurons
    if n_E > 0 && n_a_E > 0 && ~isempty(a_E)
        % c_SFA(E_indices) is n_E x 1. sum(a_E, 2) is n_E x 1.
        u_eff(E_indices) = u_eff(E_indices) - c_SFA(E_indices) .* sum(a_E, 2);
    end
    % Apply SFA effect to I neurons
    if n_I > 0 && n_a_I > 0 && ~isempty(a_I)
        % c_SFA(I_indices) is n_I x 1. sum(a_I, 2) is n_I x 1.
        u_eff(I_indices) = u_eff(I_indices) - c_SFA(I_indices) .* sum(a_I, 2);
    end
    
    r = activation_function(u_eff); % Hz, spike rate, n x 1

    p = r; % n x 1, axonal output, initially same as r

    % Apply STD effect to E neurons
    if n_E > 0 && n_b_E > 0 && ~isempty(b_E)
        % prod(b_E, 2) is n_E x 1
        p(E_indices) = p(E_indices) .* prod(b_E, 2);
    end
    % Apply STD effect to I neurons
    if n_I > 0 && n_b_I > 0 && ~isempty(b_I)
        % prod(b_I, 2) is n_I x 1
        p(I_indices) = p(I_indices) .* prod(b_I, 2);
    end

    %% derivatives
    % da_E_dt should be n_E x n_a_E (or empty)
    % da_I_dt should be n_I x n_a_I (or empty)
    % db_E_dt should be n_E x n_b_E (or empty)
    % db_I_dt should be n_I x n_b_I (or empty)
    % u_d_dt should be n x 1

    da_E_dt = [];
    if n_E > 0 && n_a_E > 0 && ~isempty(a_E)
        % r(E_indices) is n_E x 1. tau_a_E is 1 x n_a_E. Broadcasting makes (r_E - a_E) ./ tau_a_E valid.
        da_E_dt = (r(E_indices) - a_E) ./ tau_a_E;
        % Enforce no adaptation for E neurons where c_SFA is 0
        da_E_dt(c_SFA(E_indices)==0, :) = 0; 
    end

    da_I_dt = [];
    if n_I > 0 && n_a_I > 0 && ~isempty(a_I)
        da_I_dt = (r(I_indices) - a_I) ./ tau_a_I;
        da_I_dt(c_SFA(I_indices)==0, :) = 0;
    end

    db_E_dt = [];
    if n_E > 0 && n_b_E > 0 && ~isempty(b_E)
        % (1-b_E) is n_E x n_b_E. tau_b_E is 1 x n_b_E.
        % F_STD(E_indices) is n_E x 1. p(E_indices) is n_E x 1. (F_STD_E .* p_E) is n_E x 1.
        % Broadcasting makes (F_STD_E .* p_E) ./ tau_STD valid against (1-b_E)./tau_b_E if tau_STD is scalar.
        db_E_dt = (1 - b_E) ./ tau_b_E - (F_STD(E_indices) .* p(E_indices)) ./ tau_STD;
        db_E_dt(F_STD(E_indices)==0, :) = 0;
    end

    db_I_dt = [];
    if n_I > 0 && n_b_I > 0 && ~isempty(b_I)
        db_I_dt = (1 - b_I) ./ tau_b_I - (F_STD(I_indices) .* p(I_indices)) ./ tau_STD;
        db_I_dt(F_STD(I_indices)==0, :) = 0;
    end
    
    u_d_dt = (-u_d + u + M * p) ./ tau_d;

    %% load derivatives into dXdt in the correct order
    dX_dt = [da_E_dt(:); da_I_dt(:); db_E_dt(:); db_I_dt(:); u_d_dt];

end