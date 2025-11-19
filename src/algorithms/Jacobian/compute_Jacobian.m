function J = compute_Jacobian(S, params)
% COMPUTE_JACOBIAN computes the Jacobian matrix of the SRNN system
%
% J = compute_Jacobian(S, params)
%
% Computes the Jacobian matrix ∂(dS/dt)/∂S for the SRNN system at a given state.
%
% Inputs:
%   S      - State vector (N_sys_eqs × 1) with organization [a_E(:); a_I(:); b_E(:); b_I(:); x(:)]
%   params - Parameter structure containing:
%            .n, .n_E, .n_I, .E_indices, .I_indices
%            .n_a_E, .n_a_I, .n_b_E, .n_b_I
%            .tau_d, .tau_a_E, .tau_a_I
%            .tau_b_E_rec, .tau_b_E_rel, .tau_b_I_rec, .tau_b_I_rel
%            .W (connectivity matrix)
%            .c_E, .c_I (adaptation scaling parameters, optional, default = 1.0)
%            .activation_function_derivative (function handle for φ'(x))
%
% Output:
%   J - Jacobian matrix (N_sys_eqs × N_sys_eqs)
%
% The Jacobian has the block structure:
%   J = [ ∂(da_E/dt)/∂a_E  ∂(da_E/dt)/∂a_I  ∂(da_E/dt)/∂b_E  ∂(da_E/dt)/∂b_I  ∂(da_E/dt)/∂x ]
%       [ ∂(da_I/dt)/∂a_E  ∂(da_I/dt)/∂a_I  ∂(da_I/dt)/∂b_E  ∂(da_I/dt)/∂b_I  ∂(da_I/dt)/∂x ]
%       [ ∂(db_E/dt)/∂a_E  ∂(db_E/dt)/∂a_I  ∂(db_E/dt)/∂b_E  ∂(db_E/dt)/∂b_I  ∂(db_E/dt)/∂x ]
%       [ ∂(db_I/dt)/∂a_E  ∂(db_I/dt)/∂a_I  ∂(db_I/dt)/∂b_E  ∂(db_I/dt)/∂b_I  ∂(db_I/dt)/∂x ]
%       [ ∂(dx/dt)/∂a_E    ∂(dx/dt)/∂a_I    ∂(dx/dt)/∂b_E    ∂(dx/dt)/∂b_I    ∂(dx/dt)/∂x   ]

    %% Load parameters
    n = params.n;
    n_E = params.n_E;
    n_I = params.n_I;
    E_indices = params.E_indices;
    I_indices = params.I_indices;
    
    n_a_E = params.n_a_E;
    n_a_I = params.n_a_I;
    n_b_E = params.n_b_E;
    n_b_I = params.n_b_I;
    
    W = params.W;
    tau_d = params.tau_d;
    tau_a_E = params.tau_a_E;
    tau_a_I = params.tau_a_I;
    tau_b_E_rec = params.tau_b_E_rec;
    tau_b_E_rel = params.tau_b_E_rel;
    tau_b_I_rec = params.tau_b_I_rec;
    tau_b_I_rel = params.tau_b_I_rel;
    
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
        error('compute_Jacobian:MissingActivationFunctionDerivative', ...
              'params.activation_function_derivative must be provided as a function handle');
    end
    
    %% Unpack state variables
    % State organization: S = [a_E(:); a_I(:); b_E(:); b_I(:); x(:)]
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
    
    % STD states for E neurons (b_E)
    len_b_E = n_E * n_b_E;
    if len_b_E > 0
        b_E = S(current_idx + (1:len_b_E));
    else
        b_E = [];
    end
    current_idx = current_idx + len_b_E;
    
    % STD states for I neurons (b_I)
    len_b_I = n_I * n_b_I;
    if len_b_I > 0
        b_I = S(current_idx + (1:len_b_I));
    else
        b_I = [];
    end
    current_idx = current_idx + len_b_I;
    
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
    
    % Construct b vector for STD effect
    b = ones(n, 1);  % Initialize b = 1 for all neurons (no depression)
    if n_b_E > 0 && ~isempty(b_E)
        b(E_indices) = b_E;
    end
    if n_b_I > 0 && ~isempty(b_I)
        b(I_indices) = b_I;
    end
    
    % Compute activation function and its derivative at x_eff
    phi_prime_x_eff = phi_prime(x_eff); % n x 1
    
    % Compute activation function φ(x_eff) - needed for STD blocks
    if isfield(params, 'activation_function') && isa(params.activation_function, 'function_handle')
        phi_x_eff = params.activation_function(x_eff);
    else
        error('compute_Jacobian:MissingActivationFunction', ...
              'params.activation_function must be provided as a function handle');
    end
    
    % Compute firing rates r = b .* φ(x_eff) - needed for STD Jacobian blocks
    r_vec = b .* phi_x_eff;  % n x 1
    
    %% Initialize Jacobian matrix
    N_sys_eqs = len_a_E + len_a_I + len_b_E + len_b_I + n;
    J = zeros(N_sys_eqs, N_sys_eqs);
    
    %% Compute Jacobian blocks
    
    % Block row indices
    row_a_E_start = 1;
    row_a_E_end = len_a_E;
    row_a_I_start = len_a_E + 1;
    row_a_I_end = len_a_E + len_a_I;
    row_b_E_start = len_a_E + len_a_I + 1;
    row_b_E_end = len_a_E + len_a_I + len_b_E;
    row_b_I_start = len_a_E + len_a_I + len_b_E + 1;
    row_b_I_end = len_a_E + len_a_I + len_b_E + len_b_I;
    row_x_start = len_a_E + len_a_I + len_b_E + len_b_I + 1;
    row_x_end = N_sys_eqs;
    
    % Block column indices
    col_a_E_start = 1;
    col_a_E_end = len_a_E;
    col_a_I_start = len_a_E + 1;
    col_a_I_end = len_a_E + len_a_I;
    col_b_E_start = len_a_E + len_a_I + 1;
    col_b_E_end = len_a_E + len_a_I + len_b_E;
    col_b_I_start = len_a_E + len_a_I + len_b_E + 1;
    col_b_I_end = len_a_E + len_a_I + len_b_E + len_b_I;
    col_x_start = len_a_E + len_a_I + len_b_E + len_b_I + 1;
    col_x_end = N_sys_eqs;
    
    %% Block 1: ∂(da_E/dt)/∂a_E
    % da_{i,k}/dt = (-a_{i,k} + r_i) / τ_k, where r_i = b_i * φ(x_eff_i)
    % ∂(da_{i,k}/dt)/∂a_{i,j} = -δ_{kj}/τ_k - b_i*c_E*φ'(x_eff_i)/τ_k
    if len_a_E > 0
        for i = 1:n_E
            neuron_idx = E_indices(i);
            b_i = b(neuron_idx);
            for k = 1:n_a_E
                row_idx = (i-1)*n_a_E + k;
                for j = 1:n_a_E
                    col_idx = (i-1)*n_a_E + j;
                    if k == j
                        J(row_idx, col_idx) = -1/tau_a_E(k) - b_i*c_E*phi_prime_x_eff(neuron_idx)/tau_a_E(k);
                    else
                        J(row_idx, col_idx) = -b_i*c_E*phi_prime_x_eff(neuron_idx)/tau_a_E(k);
                    end
                end
            end
        end
    end
    
    %% Block 2: ∂(da_E/dt)/∂a_I
    % No coupling between E and I adaptation variables
    % This block remains zero
    
    %% Block 2b: ∂(da_E/dt)/∂b_E
    % da_{i,k}/dt = (-a_{i,k} + b_i*φ(x_eff_i)) / τ_k
    % ∂(da_{i,k}/dt)/∂b_i = φ(x_eff_i) / τ_k
    if len_a_E > 0 && len_b_E > 0
        for i = 1:n_E
            neuron_idx = E_indices(i);
            for k = 1:n_a_E
                row_idx = (i-1)*n_a_E + k;
                col_idx = col_b_E_start + (i-1)*n_b_E;  % b_E is only 1 element per E neuron when n_b_E=1
                J(row_idx, col_idx) = phi_x_eff(neuron_idx) / tau_a_E(k);
            end
        end
    end
    
    %% Block 2c: ∂(da_E/dt)/∂b_I
    % No coupling between E adaptation and I STD
    % This block remains zero
    
    %% Block 3: ∂(da_E/dt)/∂x
    % ∂(da_{i,k}/dt)/∂x_j = (b_i * φ'(x_eff_i) * δ_{ij}) / τ_k
    if len_a_E > 0
        for i = 1:n_E
            neuron_idx = E_indices(i);
            b_i = b(neuron_idx);
            for k = 1:n_a_E
                row_idx = (i-1)*n_a_E + k;
                J(row_idx, col_x_start + neuron_idx - 1) = b_i * phi_prime_x_eff(neuron_idx) / tau_a_E(k);
            end
        end
    end
    
    %% Block 4: ∂(da_I/dt)/∂a_E
    % No coupling between I and E adaptation variables
    % This block remains zero
    
    %% Block 5: ∂(da_I/dt)/∂a_I
    % da_{i,k}/dt = (-a_{i,k} + r_i) / τ_k, where r_i = b_i * φ(x_eff_i)
    % ∂(da_{i,k}/dt)/∂a_{i,j} = -δ_{kj}/τ_k - b_i*c_I*φ'(x_eff_i)/τ_k
    if len_a_I > 0
        for i = 1:n_I
            neuron_idx = I_indices(i);
            b_i = b(neuron_idx);
            for k = 1:n_a_I
                row_idx = len_a_E + (i-1)*n_a_I + k;
                for j = 1:n_a_I
                    col_idx = len_a_E + (i-1)*n_a_I + j;
                    if k == j
                        J(row_idx, col_idx) = -1/tau_a_I(k) - b_i*c_I*phi_prime_x_eff(neuron_idx)/tau_a_I(k);
                    else
                        J(row_idx, col_idx) = -b_i*c_I*phi_prime_x_eff(neuron_idx)/tau_a_I(k);
                    end
                end
            end
        end
    end
    
    %% Block 5b: ∂(da_I/dt)/∂b_E
    % No coupling between I adaptation and E STD
    % This block remains zero
    
    %% Block 5c: ∂(da_I/dt)/∂b_I
    % da_{i,k}/dt = (-a_{i,k} + b_i*φ(x_eff_i)) / τ_k
    % ∂(da_{i,k}/dt)/∂b_i = φ(x_eff_i) / τ_k
    if len_a_I > 0 && len_b_I > 0
        for i = 1:n_I
            neuron_idx = I_indices(i);
            for k = 1:n_a_I
                row_idx = len_a_E + (i-1)*n_a_I + k;
                col_idx = col_b_I_start + (i-1)*n_b_I;  % b_I is only 1 element per I neuron when n_b_I=1
                J(row_idx, col_idx) = phi_x_eff(neuron_idx) / tau_a_I(k);
            end
        end
    end
    
    %% Block 6: ∂(da_I/dt)/∂x
    % ∂(da_{i,k}/dt)/∂x_j = (b_i * φ'(x_eff_i) * δ_{ij}) / τ_k
    if len_a_I > 0
        for i = 1:n_I
            neuron_idx = I_indices(i);
            b_i = b(neuron_idx);
            for k = 1:n_a_I
                row_idx = len_a_E + (i-1)*n_a_I + k;
                J(row_idx, col_x_start + neuron_idx - 1) = b_i * phi_prime_x_eff(neuron_idx) / tau_a_I(k);
            end
        end
    end
    
    %% NEW BLOCKS FOR b DYNAMICS
    
    %% Block 7new: ∂(db_E/dt)/∂a_E
    % db_i/dt = (1-b_i)/τ_rec - (b_i*r_i)/τ_rel, where r_i = b_i*φ(x_eff_i)
    % ∂(db_i/dt)/∂a_{i,k} = -b_i * ∂r_i/∂a_{i,k} / τ_rel = -b_i * b_i * (-c_E) * φ'(x_eff_i) / τ_rel
    %                      = b_i^2 * c_E * φ'(x_eff_i) / τ_rel
    if len_b_E > 0 && len_a_E > 0
        for i = 1:n_E
            neuron_idx = E_indices(i);
            b_i = b(neuron_idx);
            row_idx = row_b_E_start + (i-1)*n_b_E;
            for k = 1:n_a_E
                col_idx = (i-1)*n_a_E + k;
                J(row_idx, col_idx) = b_i^2 * c_E * phi_prime_x_eff(neuron_idx) / tau_b_E_rel;
            end
        end
    end
    
    %% Block 8new: ∂(db_E/dt)/∂a_I
    % No coupling between E STD and I adaptation
    % This block remains zero
    
    %% Block 9new: ∂(db_E/dt)/∂b_E
    % db_i/dt = (1-b_i)/τ_rec - (b_i*r_i)/τ_rel
    % ∂(db_i/dt)/∂b_i = -1/τ_rec - r_i/τ_rel - b_i * ∂r_i/∂b_i / τ_rel
    %                 = -1/τ_rec - r_i/τ_rel - b_i * φ(x_eff_i) / τ_rel
    %                 = -1/τ_rec - 2*r_i/τ_rel  (since r_i = b_i*φ(x_eff_i))
    if len_b_E > 0
        for i = 1:n_E
            neuron_idx = E_indices(i);
            row_idx = row_b_E_start + (i-1)*n_b_E;
            col_idx = col_b_E_start + (i-1)*n_b_E;
            J(row_idx, col_idx) = -1/tau_b_E_rec - 2*r_vec(neuron_idx)/tau_b_E_rel;
        end
    end
    
    %% Block 10new: ∂(db_E/dt)/∂b_I
    % No coupling between E STD and I STD
    % This block remains zero
    
    %% Block 11new: ∂(db_E/dt)/∂x
    % db_i/dt = (1-b_i)/τ_rec - (b_i*r_i)/τ_rel
    % ∂(db_i/dt)/∂x_j = -b_i * ∂r_i/∂x_j / τ_rel = -b_i * b_i * φ'(x_eff_i) * δ_{ij} / τ_rel
    %                  = -b_i^2 * φ'(x_eff_i) / τ_rel
    if len_b_E > 0
        for i = 1:n_E
            neuron_idx = E_indices(i);
            b_i = b(neuron_idx);
            row_idx = row_b_E_start + (i-1)*n_b_E;
            J(row_idx, col_x_start + neuron_idx - 1) = -b_i^2 * phi_prime_x_eff(neuron_idx) / tau_b_E_rel;
        end
    end
    
    %% Block 12new: ∂(db_I/dt)/∂a_E
    % No coupling between I STD and E adaptation
    % This block remains zero
    
    %% Block 13new: ∂(db_I/dt)/∂a_I
    % db_i/dt = (1-b_i)/τ_rec - (b_i*r_i)/τ_rel
    % Similar to Block 7new but for I neurons
    if len_b_I > 0 && len_a_I > 0
        for i = 1:n_I
            neuron_idx = I_indices(i);
            b_i = b(neuron_idx);
            row_idx = row_b_I_start + (i-1)*n_b_I;
            for k = 1:n_a_I
                col_idx = col_a_I_start + (i-1)*n_a_I + k - 1;
                J(row_idx, col_idx) = b_i^2 * c_I * phi_prime_x_eff(neuron_idx) / tau_b_I_rel;
            end
        end
    end
    
    %% Block 14new: ∂(db_I/dt)/∂b_E
    % No coupling between I STD and E STD
    % This block remains zero
    
    %% Block 15new: ∂(db_I/dt)/∂b_I
    % Similar to Block 9new but for I neurons
    if len_b_I > 0
        for i = 1:n_I
            neuron_idx = I_indices(i);
            row_idx = row_b_I_start + (i-1)*n_b_I;
            col_idx = col_b_I_start + (i-1)*n_b_I;
            J(row_idx, col_idx) = -1/tau_b_I_rec - 2*r_vec(neuron_idx)/tau_b_I_rel;
        end
    end
    
    %% Block 16new: ∂(db_I/dt)/∂x
    % Similar to Block 11new but for I neurons
    if len_b_I > 0
        for i = 1:n_I
            neuron_idx = I_indices(i);
            b_i = b(neuron_idx);
            row_idx = row_b_I_start + (i-1)*n_b_I;
            J(row_idx, col_x_start + neuron_idx - 1) = -b_i^2 * phi_prime_x_eff(neuron_idx) / tau_b_I_rel;
        end
    end
    
    %% ORIGINAL BLOCKS FOR dx/dt (updated to account for b in r)
    
    %% Block 7: ∂(dx/dt)/∂a_E
    % dx_i/dt = (-x_i + Σ_j w_ij r_j + u_i) / τ_d, where r_j = b_j*φ(x_eff_j)
    % ∂(dx_i/dt)/∂a_{j,k} = (w_ij * b_j * φ'(x_eff_j) * (-c_E)) / τ_d for j in E_indices
    if len_a_E > 0
        for i = 1:n
            row_idx = row_x_start + i - 1;
            for j = 1:n_E
                neuron_j = E_indices(j);
                b_j = b(neuron_j);
                for k = 1:n_a_E
                    col_idx = (j-1)*n_a_E + k;
                    J(row_idx, col_idx) = (-c_E * W(i, neuron_j) * b_j * phi_prime_x_eff(neuron_j)) / tau_d;
                end
            end
        end
    end
    
    %% Block 8: ∂(dx/dt)/∂a_I
    % ∂(dx_i/dt)/∂a_{j,k} = (w_ij * b_j * φ'(x_eff_j) * (-c_I)) / τ_d for j in I_indices
    if len_a_I > 0
        for i = 1:n
            row_idx = row_x_start + i - 1;
            for j = 1:n_I
                neuron_j = I_indices(j);
                b_j = b(neuron_j);
                for k = 1:n_a_I
                    col_idx = col_a_I_start + (j-1)*n_a_I + k - 1;
                    J(row_idx, col_idx) = (-c_I * W(i, neuron_j) * b_j * phi_prime_x_eff(neuron_j)) / tau_d;
                end
            end
        end
    end
    
    %% Block 8b: ∂(dx/dt)/∂b_E
    % dx_i/dt = (-x_i + Σ_j w_ij r_j + u_i) / τ_d, where r_j = b_j*φ(x_eff_j)
    % ∂(dx_i/dt)/∂b_j = (w_ij * φ(x_eff_j)) / τ_d for j in E_indices
    if len_b_E > 0
        for i = 1:n
            row_idx = row_x_start + i - 1;
            for j = 1:n_E
                neuron_j = E_indices(j);
                col_idx = col_b_E_start + (j-1)*n_b_E;
                J(row_idx, col_idx) = (W(i, neuron_j) * phi_x_eff(neuron_j)) / tau_d;
            end
        end
    end
    
    %% Block 8c: ∂(dx/dt)/∂b_I
    % ∂(dx_i/dt)/∂b_j = (w_ij * φ(x_eff_j)) / τ_d for j in I_indices
    if len_b_I > 0
        for i = 1:n
            row_idx = row_x_start + i - 1;
            for j = 1:n_I
                neuron_j = I_indices(j);
                col_idx = col_b_I_start + (j-1)*n_b_I;
                J(row_idx, col_idx) = (W(i, neuron_j) * phi_x_eff(neuron_j)) / tau_d;
            end
        end
    end
    
    %% Block 9: ∂(dx/dt)/∂x
    % ∂(dx_i/dt)/∂x_j = (-δ_{ij} + w_ij * b_j * φ'(x_eff_j)) / τ_d
    for i = 1:n
        row_idx = row_x_start + i - 1;
        for j = 1:n
            col_idx = col_x_start + j - 1;
            b_j = b(j);
            if i == j
                J(row_idx, col_idx) = (-1 + W(i, j) * b_j * phi_prime_x_eff(j)) / tau_d;
            else
                J(row_idx, col_idx) = (W(i, j) * b_j * phi_prime_x_eff(j)) / tau_d;
            end
        end
    end
    
end

