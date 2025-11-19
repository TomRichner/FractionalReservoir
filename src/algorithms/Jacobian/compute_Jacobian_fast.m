function J = compute_Jacobian_fast(S, params)
% COMPUTE_JACOBIAN_FAST Sparse/vectorized Jacobian assembly for the SRNN system.
%
%   J = compute_Jacobian_fast(S, params)
%
% This version mirrors compute_Jacobian but assembles the matrix using sparse
% block operations (kron, spdiags) for improved scalability in Lyapunov
% calculations and other workflows that require frequent Jacobian evaluations.
%
% Assumptions:
%   - At most one short-term depression state per neuron (n_b_E, n_b_I ∈ {0,1}),
%     matching the current SRNN_reservoir dynamics.

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
    
    if n_b_E > 1 || n_b_I > 1
        error('compute_Jacobian_fast:UnsupportedSTDStates', ...
              'Fast Jacobian currently supports at most one STD state per neuron.');
    end
    
    W = params.W;
    tau_d = params.tau_d;
    tau_a_E = params.tau_a_E;
    tau_a_I = params.tau_a_I;
    tau_b_E_rec = params.tau_b_E_rec;
    tau_b_E_rel = params.tau_b_E_rel;
    tau_b_I_rec = params.tau_b_I_rec;
    tau_b_I_rel = params.tau_b_I_rel;
    
    c_E = safe_get(params, 'c_E', 1.0);
    c_I = safe_get(params, 'c_I', 1.0);
    
    if ~isfield(params, 'activation_function_derivative') || ...
       ~isa(params.activation_function_derivative, 'function_handle')
        error('compute_Jacobian_fast:MissingActivationFunctionDerivative', ...
              'params.activation_function_derivative must be provided as a function handle');
    end
    phi_prime = params.activation_function_derivative;
    
    if ~isfield(params, 'activation_function') || ...
       ~isa(params.activation_function, 'function_handle')
        error('compute_Jacobian_fast:MissingActivationFunction', ...
              'params.activation_function must be provided as a function handle');
    end
    phi_fun = params.activation_function;
    
    %% Unpack state variables
    current_idx = 0;
    len_a_E = n_E * n_a_E;
    len_a_I = n_I * n_a_I;
    len_b_E = n_E * n_b_E;
    len_b_I = n_I * n_b_I;
    
    if len_a_E > 0
        a_E = reshape(S(current_idx + (1:len_a_E)), n_E, n_a_E);
    else
        a_E = [];
    end
    current_idx = current_idx + len_a_E;
    
    if len_a_I > 0
        a_I = reshape(S(current_idx + (1:len_a_I)), n_I, n_a_I);
    else
        a_I = [];
    end
    current_idx = current_idx + len_a_I;
    
    if len_b_E > 0
        b_E = S(current_idx + (1:len_b_E));
    else
        b_E = [];
    end
    current_idx = current_idx + len_b_E;
    
    if len_b_I > 0
        b_I = S(current_idx + (1:len_b_I));
    else
        b_I = [];
    end
    current_idx = current_idx + len_b_I;
    
    x = S(current_idx + (1:n));
    
    %% Effective potentials and rates
    x_eff = x;
    if len_a_E > 0
        x_eff(E_indices) = x_eff(E_indices) - c_E * sum(a_E, 2);
    end
    if len_a_I > 0
        x_eff(I_indices) = x_eff(I_indices) - c_I * sum(a_I, 2);
    end
    
    b = ones(n, 1);
    if len_b_E > 0
        b(E_indices) = b_E;
    end
    if len_b_I > 0
        b(I_indices) = b_I;
    end
    
    phi_x_eff = phi_fun(x_eff);
    phi_prime_x_eff = phi_prime(x_eff);
    r_vec = b .* phi_x_eff;
    
    %% Dimensions and indexing
    N_sys_eqs = len_a_E + len_a_I + len_b_E + len_b_I + n;
    
    row_a_E = 1:len_a_E;
    row_a_I = len_a_E + (1:len_a_I);
    row_b_E = len_a_E + len_a_I + (1:len_b_E);
    row_b_I = len_a_E + len_a_I + len_b_E + (1:len_b_I);
    row_x   = len_a_E + len_a_I + len_b_E + len_b_I + (1:n);
    
    col_a_E = row_a_E;
    col_a_I = row_a_I;
    col_b_E = row_b_E;
    col_b_I = row_b_I;
    col_x   = row_x;
    
    W_sparse = sparse(W);
    J = sparse(N_sys_eqs, N_sys_eqs);
    
    %% Helper Kronecker scaffolds
    if len_a_E > 0
        tau_inv_E = 1 ./ tau_a_E(:);
        diag_block_E = kron(speye(n_E), spdiags(-tau_inv_E, 0, n_a_E, n_a_E));
        gamma_E = c_E * (b(E_indices) .* phi_prime_x_eff(E_indices));
        row_template_E = sparse(tau_inv_E * ones(1, n_a_E));
        coupling_block_E = kron(spdiags(-gamma_E, 0, n_E, n_E), row_template_E);
        J(row_a_E, col_a_E) = diag_block_E + coupling_block_E;
        
        beta_E = b(E_indices) .* phi_prime_x_eff(E_indices);
        vals = kron(beta_E, tau_inv_E);
        rows = (1:len_a_E)';
        cols = repelem(E_indices(:), n_a_E);
        J(row_a_E, col_x) = sparse(rows, cols, vals, len_a_E, n);
        
        if len_b_E > 0
            phi_E = phi_x_eff(E_indices);
            J(row_a_E, col_b_E) = kron(spdiags(phi_E, 0, n_E, n_E), sparse(tau_inv_E));
        end
    end
    
    if len_a_I > 0
        tau_inv_I = 1 ./ tau_a_I(:);
        diag_block_I = kron(speye(n_I), spdiags(-tau_inv_I, 0, n_a_I, n_a_I));
        gamma_I = c_I * (b(I_indices) .* phi_prime_x_eff(I_indices));
        row_template_I = sparse(tau_inv_I * ones(1, n_a_I));
        coupling_block_I = kron(spdiags(-gamma_I, 0, n_I, n_I), row_template_I);
        J(row_a_I, col_a_I) = diag_block_I + coupling_block_I;
        
        beta_I = b(I_indices) .* phi_prime_x_eff(I_indices);
        vals = kron(beta_I, tau_inv_I);
        rows = (1:len_a_I)';
        cols = repelem(I_indices(:), n_a_I);
        J(row_a_I, col_x) = sparse(rows, cols, vals, len_a_I, n);
        
        if len_b_I > 0
            phi_I = phi_x_eff(I_indices);
            J(row_a_I, col_b_I) = kron(spdiags(phi_I, 0, n_I, n_I), sparse(tau_inv_I));
        end
    end
    
    %% STD blocks (E)
    if len_b_E > 0
        phi_prime_E = phi_prime_x_eff(E_indices);
        coeff_a_E = (b(E_indices).^2) * c_E .* phi_prime_E / tau_b_E_rel;
        if len_a_E > 0
            J(row_b_E, col_a_E) = kron(spdiags(coeff_a_E, 0, n_E, n_E), sparse(ones(1, n_a_E)));
        end
        diag_vals_b_E = -1/tau_b_E_rec - 2 * r_vec(E_indices) / tau_b_E_rel;
        J(row_b_E, col_b_E) = spdiags(diag_vals_b_E, 0, len_b_E, len_b_E);
        J(row_b_E, col_x) = sparse(1:n_E, E_indices, - (b(E_indices).^2) .* phi_prime_E / tau_b_E_rel, n_E, n);
    end
    
    %% STD blocks (I)
    if len_b_I > 0
        phi_prime_I = phi_prime_x_eff(I_indices);
        coeff_a_I = (b(I_indices).^2) * c_I .* phi_prime_I / tau_b_I_rel;
        if len_a_I > 0
            J(row_b_I, col_a_I) = kron(spdiags(coeff_a_I, 0, n_I, n_I), sparse(ones(1, n_a_I)));
        end
        diag_vals_b_I = -1/tau_b_I_rec - 2 * r_vec(I_indices) / tau_b_I_rel;
        J(row_b_I, col_b_I) = spdiags(diag_vals_b_I, 0, len_b_I, len_b_I);
        J(row_b_I, col_x) = sparse(1:n_I, I_indices, - (b(I_indices).^2) .* phi_prime_I / tau_b_I_rel, n_I, n);
    end
    
    %% dx/dt blocks
    if len_a_E > 0
        replicate_a_E = kron(speye(n_E), ones(1, n_a_E));
        block = -c_E * W_sparse(:, E_indices) * spdiags(b(E_indices) .* phi_prime_x_eff(E_indices), 0, n_E, n_E);
        J(row_x, col_a_E) = (block * replicate_a_E) / tau_d;
    end
    
    if len_a_I > 0
        replicate_a_I = kron(speye(n_I), ones(1, n_a_I));
        block = -c_I * W_sparse(:, I_indices) * spdiags(b(I_indices) .* phi_prime_x_eff(I_indices), 0, n_I, n_I);
        J(row_x, col_a_I) = (block * replicate_a_I) / tau_d;
    end
    
    if len_b_E > 0
        replicate_b_E = kron(speye(n_E), ones(1, max(1, n_b_E)));
        block = W_sparse(:, E_indices) * spdiags(phi_x_eff(E_indices), 0, n_E, n_E);
        J(row_x, col_b_E) = (block * replicate_b_E) / tau_d;
    end
    
    if len_b_I > 0
        replicate_b_I = kron(speye(n_I), ones(1, max(1, n_b_I)));
        block = W_sparse(:, I_indices) * spdiags(phi_x_eff(I_indices), 0, n_I, n_I);
        J(row_x, col_b_I) = (block * replicate_b_I) / tau_d;
    end
    
    diag_term = spdiags(-ones(n,1)/tau_d, 0, n, n);
    gain_diag = spdiags(b .* phi_prime_x_eff, 0, n, n);
    J(row_x, col_x) = diag_term + (W_sparse * gain_diag) / tau_d;
end

function value = safe_get(params, field, default_value)
    if isfield(params, field)
        value = params.(field);
    else
        value = default_value;
    end
end

