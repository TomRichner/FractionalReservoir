function J = compute_Jacobian(S, params)
% COMPUTE_JACOBIAN computes the Jacobian matrix of the SRNN system
%
% J = compute_Jacobian(S, params)
%
% Computes the Jacobian matrix ∂(dS/dt)/∂S for the SRNN system at a given state.
%
% Inputs:
%   S      - State vector (N_sys_eqs × 1) with organization [a_E(:); a_I(:); x(:)]
%   params - Parameter structure containing:
%            .n, .n_E, .n_I, .E_indices, .I_indices
%            .n_a_E, .n_a_I, .tau_d, .tau_a_E, .tau_a_I
%            .W (connectivity matrix)
%            .c_E, .c_I (adaptation scaling parameters, optional, default = 1.0)
%            .activation_function_derivative (function handle for φ'(x))
%
% Output:
%   J - Jacobian matrix (N_sys_eqs × N_sys_eqs)
%
% The Jacobian has the block structure:
%   J = [ ∂(da_E/dt)/∂a_E    ∂(da_E/dt)/∂a_I    ∂(da_E/dt)/∂x   ]
%       [ ∂(da_I/dt)/∂a_E    ∂(da_I/dt)/∂a_I    ∂(da_I/dt)/∂x   ]
%       [ ∂(dx/dt)/∂a_E      ∂(dx/dt)/∂a_I      ∂(dx/dt)/∂x     ]

    %% Load parameters
    n = params.n;
    n_E = params.n_E;
    n_I = params.n_I;
    E_indices = params.E_indices;
    I_indices = params.I_indices;
    
    n_a_E = params.n_a_E;
    n_a_I = params.n_a_I;
    
    W = params.W;
    tau_d = params.tau_d;
    tau_a_E = params.tau_a_E;
    tau_a_I = params.tau_a_I;
    
    % Adaptation scaling parameters
    if isfield(params, 'c_E')
        c_E = params.c_E;
    else
        c_E = 1.0;
    end
    
    if isfield(params, 'c_I')
        c_I = params.c_I;
    else
        c_I = 1.0;
    end
    
    % Activation function derivative
    if isfield(params, 'activation_function_derivative') && isa(params.activation_function_derivative, 'function_handle')
        phi_prime = params.activation_function_derivative;
    else
        % Default to tanh derivative
        phi_prime = @(x) 1 - tanh(x).^2;
    end
    
    %% Unpack state variables
    % State organization: S = [a_E(:); a_I(:); x(:)]
    current_idx = 0;
    
    % Adaptation states for E neurons (a_E)
    len_a_E = n_E * n_a_E;
    if len_a_E > 0
        a_E = reshape(S(current_idx + (1:len_a_E)), n_E, n_a_E);
    else
        a_E = [];
    end
    current_idx = current_idx + len_a_E;
    
    % Adaptation states for I neurons (a_I)
    len_a_I = n_I * n_a_I;
    if len_a_I > 0
        a_I = reshape(S(current_idx + (1:len_a_I)), n_I, n_a_I);
    else
        a_I = [];
    end
    current_idx = current_idx + len_a_I;
    
    % Dendritic states (x)
    x = S(current_idx + (1:n));
    
    %% Compute effective dendritic potential and its derivative
    x_eff = x; % n x 1
    
    % Apply adaptation effect to E neurons (scaled by c_E)
    if n_E > 0 && n_a_E > 0 && ~isempty(a_E)
        x_eff(E_indices) = x_eff(E_indices) - c_E * sum(a_E, 2);
    end
    
    % Apply adaptation effect to I neurons (scaled by c_I)
    if n_I > 0 && n_a_I > 0 && ~isempty(a_I)
        x_eff(I_indices) = x_eff(I_indices) - c_I * sum(a_I, 2);
    end
    
    % Compute derivative of activation function at x_eff
    phi_prime_x_eff = phi_prime(x_eff); % n x 1
    
    %% Initialize Jacobian matrix
    N_sys_eqs = len_a_E + len_a_I + n;
    J = zeros(N_sys_eqs, N_sys_eqs);
    
    %% Compute Jacobian blocks
    
    % Block row indices
    row_a_E_start = 1;
    row_a_E_end = len_a_E;
    row_a_I_start = len_a_E + 1;
    row_a_I_end = len_a_E + len_a_I;
    row_x_start = len_a_E + len_a_I + 1;
    row_x_end = N_sys_eqs;
    
    % Block column indices
    col_a_E_start = 1;
    col_a_E_end = len_a_E;
    col_a_I_start = len_a_E + 1;
    col_a_I_end = len_a_E + len_a_I;
    col_x_start = len_a_E + len_a_I + 1;
    col_x_end = N_sys_eqs;
    
    %% Block 1: ∂(da_E/dt)/∂a_E
    % da_{i,k}/dt = (-a_{i,k} + r_i) / τ_k
    % ∂(da_{i,k}/dt)/∂a_{i,j} = -δ_{kj}/τ_k - c_E*φ'(x_eff_i)/τ_k
    if len_a_E > 0
        for i = 1:n_E
            neuron_idx = E_indices(i);
            for k = 1:n_a_E
                row_idx = (i-1)*n_a_E + k;
                for j = 1:n_a_E
                    col_idx = (i-1)*n_a_E + j;
                    if k == j
                        J(row_idx, col_idx) = -1/tau_a_E(k) - c_E*phi_prime_x_eff(neuron_idx)/tau_a_E(k);
                    else
                        J(row_idx, col_idx) = -c_E*phi_prime_x_eff(neuron_idx)/tau_a_E(k);
                    end
                end
            end
        end
    end
    
    %% Block 2: ∂(da_E/dt)/∂a_I
    % No coupling between E and I adaptation variables
    % This block remains zero
    
    %% Block 3: ∂(da_E/dt)/∂x
    % ∂(da_{i,k}/dt)/∂x_j = (φ'(x_eff_i) * δ_{ij}) / τ_k
    if len_a_E > 0
        for i = 1:n_E
            neuron_idx = E_indices(i);
            for k = 1:n_a_E
                row_idx = (i-1)*n_a_E + k;
                J(row_idx, col_x_start + neuron_idx - 1) = phi_prime_x_eff(neuron_idx) / tau_a_E(k);
            end
        end
    end
    
    %% Block 4: ∂(da_I/dt)/∂a_E
    % No coupling between I and E adaptation variables
    % This block remains zero
    
    %% Block 5: ∂(da_I/dt)/∂a_I
    % da_{i,k}/dt = (-a_{i,k} + r_i) / τ_k
    % ∂(da_{i,k}/dt)/∂a_{i,j} = -δ_{kj}/τ_k - c_I*φ'(x_eff_i)/τ_k
    if len_a_I > 0
        for i = 1:n_I
            neuron_idx = I_indices(i);
            for k = 1:n_a_I
                row_idx = len_a_E + (i-1)*n_a_I + k;
                for j = 1:n_a_I
                    col_idx = len_a_E + (i-1)*n_a_I + j;
                    if k == j
                        J(row_idx, col_idx) = -1/tau_a_I(k) - c_I*phi_prime_x_eff(neuron_idx)/tau_a_I(k);
                    else
                        J(row_idx, col_idx) = -c_I*phi_prime_x_eff(neuron_idx)/tau_a_I(k);
                    end
                end
            end
        end
    end
    
    %% Block 6: ∂(da_I/dt)/∂x
    % ∂(da_{i,k}/dt)/∂x_j = (φ'(x_eff_i) * δ_{ij}) / τ_k
    if len_a_I > 0
        for i = 1:n_I
            neuron_idx = I_indices(i);
            for k = 1:n_a_I
                row_idx = len_a_E + (i-1)*n_a_I + k;
                J(row_idx, col_x_start + neuron_idx - 1) = phi_prime_x_eff(neuron_idx) / tau_a_I(k);
            end
        end
    end
    
    %% Block 7: ∂(dx/dt)/∂a_E
    % dx_i/dt = -x_i/τ_d + Σ_j w_ij r_j + u_i
    % ∂(dx_i/dt)/∂a_{j,k} = w_ij * φ'(x_eff_j) * (-c_E) for j in E_indices
    if len_a_E > 0
        for i = 1:n
            row_idx = row_x_start + i - 1;
            for j = 1:n_E
                neuron_j = E_indices(j);
                for k = 1:n_a_E
                    col_idx = (j-1)*n_a_E + k;
                    J(row_idx, col_idx) = -c_E * W(i, neuron_j) * phi_prime_x_eff(neuron_j);
                end
            end
        end
    end
    
    %% Block 8: ∂(dx/dt)/∂a_I
    % ∂(dx_i/dt)/∂a_{j,k} = w_ij * φ'(x_eff_j) * (-c_I) for j in I_indices
    if len_a_I > 0
        for i = 1:n
            row_idx = row_x_start + i - 1;
            for j = 1:n_I
                neuron_j = I_indices(j);
                for k = 1:n_a_I
                    col_idx = len_a_E + (j-1)*n_a_I + k;
                    J(row_idx, col_idx) = -c_I * W(i, neuron_j) * phi_prime_x_eff(neuron_j);
                end
            end
        end
    end
    
    %% Block 9: ∂(dx/dt)/∂x
    % ∂(dx_i/dt)/∂x_j = -δ_{ij}/τ_d + w_ij * φ'(x_eff_j)
    for i = 1:n
        row_idx = row_x_start + i - 1;
        for j = 1:n
            col_idx = col_x_start + j - 1;
            if i == j
                J(row_idx, col_idx) = -1/tau_d + W(i, j) * phi_prime_x_eff(j);
            else
                J(row_idx, col_idx) = W(i, j) * phi_prime_x_eff(j);
            end
        end
    end
    
end

