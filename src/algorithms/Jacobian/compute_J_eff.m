function J_eff = compute_J_eff(S, params)
% COMPUTE_J_EFF Compute the effective Jacobian of x-dynamics w.r.t. x
%
% Syntax:
%   J_eff = compute_J_eff(S, params)
%
% Description:
%   Computes the Jacobian of dx/dt with respect to x, treating adaptation
%   variables a and synaptic depression variables b as constants (frozen
%   at the current state). This is the effective connectivity matrix that
%   describes the linearized dynamics of x around the current state.
%
%   Based on the derivation in J_eff_notes.md:
%     J_eff(x,a,b) = (1/tau_d) * (-I + W * G)
%   where
%     G = diag(b_i * phi'(x_i - c * sum_k(a_{ik})))
%
% Inputs:
%   S      - State vector (N_sys_eqs x 1) at a single time point
%            State organization: [a_E(:); a_I(:); b_E(:); b_I(:); x(:)]
%   params - Struct containing network parameters, including:
%            .n - total number of neurons
%            .W - connectivity matrix (n x n)
%            .tau_d - dendritic time constant
%            .activation_function_derivative - function handle for phi'(x)
%            .c_E, .c_I - adaptation scaling parameters
%            .E_indices, .I_indices - neuron type indices
%
% Outputs:
%   J_eff  - Effective Jacobian matrix (n x n)
%
% Example:
%   J_eff = compute_J_eff(S_out(100,:)', params);
%   eigenvalues = eig(J_eff);
%
% See also: compute_Jacobian_fast, SRNN_reservoir

    %% Load parameters
    n = params.n;
    W = params.W;
    tau_d = params.tau_d;
    
    % Get adaptation scaling parameters
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
    
    % Get activation function derivative
    if ~isfield(params, 'activation_function_derivative') || ...
       ~isa(params.activation_function_derivative, 'function_handle')
        error('compute_J_eff:MissingActivationFunctionDerivative', ...
              'params.activation_function_derivative must be provided as a function handle');
    end
    phi_prime = params.activation_function_derivative;
    
    %% Unpack state variables to get actual a, b, and x values
    current_idx = 0;
    
    % Unpack adaptation states for E neurons (a_E)
    len_a_E = params.n_E * params.n_a_E;
    if len_a_E > 0
        a_E = reshape(S(current_idx + (1:len_a_E)), params.n_E, params.n_a_E);
    else
        a_E = [];
    end
    current_idx = current_idx + len_a_E;
    
    % Unpack adaptation states for I neurons (a_I)
    len_a_I = params.n_I * params.n_a_I;
    if len_a_I > 0
        a_I = reshape(S(current_idx + (1:len_a_I)), params.n_I, params.n_a_I);
    else
        a_I = [];
    end
    current_idx = current_idx + len_a_I;
    
    % Unpack STD states for E neurons (b_E)
    len_b_E = params.n_E * params.n_b_E;
    if len_b_E > 0
        b_E = S(current_idx + (1:len_b_E));
    else
        b_E = [];
    end
    current_idx = current_idx + len_b_E;
    
    % Unpack STD states for I neurons (b_I)
    len_b_I = params.n_I * params.n_b_I;
    if len_b_I > 0
        b_I = S(current_idx + (1:len_b_I));
    else
        b_I = [];
    end
    current_idx = current_idx + len_b_I;
    
    % Unpack dendritic states (x)
    x = S(current_idx + (1:n));
    
    %% Compute effective dendritic state x_eff = x - c * sum(a)
    x_eff = x;  % n x 1
    
    % Apply adaptation effect to E neurons (scaled by c_E)
    if params.n_E > 0 && params.n_a_E > 0 && ~isempty(a_E)
        % sum(a_E, 2) sums across all adaptation time constants
        x_eff(params.E_indices) = x_eff(params.E_indices) - c_E * sum(a_E, 2);
    end
    
    % Apply adaptation effect to I neurons (scaled by c_I)
    if params.n_I > 0 && params.n_a_I > 0 && ~isempty(a_I)
        % sum(a_I, 2) sums across all adaptation time constants
        x_eff(params.I_indices) = x_eff(params.I_indices) - c_I * sum(a_I, 2);
    end
    
    %% Construct b vector (defaults to ones, but replace if STD is present)
    b = ones(n, 1);
    if params.n_b_E > 0 && ~isempty(b_E)
        b(params.E_indices) = b_E;
    end
    if params.n_b_I > 0 && ~isempty(b_I)
        b(params.I_indices) = b_I;
    end
    
    %% Compute gain vector g = b .* phi'(x_eff)
    g = b .* phi_prime(x_eff);  % n x 1
    
    %% Construct gain matrix G = diag(g)
    G = spdiags(g, 0, n, n);  % n x n sparse diagonal matrix
    
    %% Compute J_eff = (1/tau_d) * (-I + W * G)
    J_eff = (-speye(n, n) + W * G) / tau_d;
    
end

