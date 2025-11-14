function [r_ts, p_ts] = compute_dependent_variables(a_E_ts, a_I_ts, b_E_ts, b_I_ts, u_d_ts, params)
    % Computes firing rates (r_ts) and axonal outputs (p_ts) for all neurons over time,
    % vectorized across time, considering separate E/I adaptation and depression states.
    %
    % Inputs:
    %   a_E_ts   - SFA variables for E neurons (n_E x n_a_E x nt or empty)
    %   a_I_ts   - SFA variables for I neurons (n_I x n_a_I x nt or empty)
    %   b_E_ts   - STD variables for E neurons (n_E x n_b_E x nt or empty)
    %   b_I_ts   - STD variables for I neurons (n_I x n_b_I x nt or empty)
    %   u_d_ts   - Dendritic potential (n x nt)
    %   params   - Struct containing parameters including:
    %              n, n_E, n_I, n_a_E, n_a_I, n_b_E, n_b_I, 
    %              c_SFA (n x 1), F_STD (n x 1), E_indices, I_indices,
    %              activation_function (function handle).
    %
    % Outputs:
    %   r_ts     - Firing rates (n x nt)
    %   p_ts     - Axonal outputs (n x nt)

    if isempty(u_d_ts)
        r_ts = [];
        p_ts = [];
        return;
    end

    n = params.n;
    nt = size(u_d_ts, 2);

    if isfield(params, 'activation_function') && isa(params.activation_function, 'function_handle')
        activation_function = params.activation_function;
    else
        activation_function = @(x) max(0, x);
    end

    u_eff_ts = u_d_ts; % n x nt, effective dendritic potential

    E_idx = params.E_indices;
    I_idx = params.I_indices;
    c_SFA_full = params.c_SFA; % n x 1

    % Apply SFA effect to E neurons (vectorized)
    if params.n_E > 0 && params.n_a_E > 0 && ~isempty(a_E_ts)
        % a_E_ts is n_E x n_a_E x nt
        % sum(a_E_ts, 2) is n_E x 1 x nt
        % squeeze_dim(sum(a_E_ts, 2), 2) is n_E x nt
        squeezed_sum_a_E = squeeze_dim(sum(a_E_ts, 2), 2); % Result: n_E x nt
        % c_SFA_full(E_idx) is n_E x 1. MATLAB broadcasts for element-wise multiplication.
        u_eff_ts(E_idx, :) = u_eff_ts(E_idx, :) - c_SFA_full(E_idx) .* squeezed_sum_a_E;
    end

    % Apply SFA effect to I neurons (vectorized)
    if params.n_I > 0 && params.n_a_I > 0 && ~isempty(a_I_ts)
        % a_I_ts is n_I x n_a_I x nt
        % sum(a_I_ts, 2) is n_I x 1 x nt
        % squeeze_dim(sum(a_I_ts, 2), 2) is n_I x nt
        squeezed_sum_a_I = squeeze_dim(sum(a_I_ts, 2), 2); % Result: n_I x nt
        u_eff_ts(I_idx, :) = u_eff_ts(I_idx, :) - c_SFA_full(I_idx) .* squeezed_sum_a_I;
    end
    
    r_ts = activation_function(u_eff_ts); % n x nt, (firing rates)
    
    p_ts = r_ts; % n x nt, axonal output, initially same as r_ts

    % Apply STD effect to E neurons (vectorized)
    if params.n_E > 0 && params.n_b_E > 0 && ~isempty(b_E_ts)
        % b_E_ts is n_E x n_b_E x nt
        % prod(b_E_ts, 2) is n_E x 1 x nt
        % squeeze_dim(prod(b_E_ts, 2), 2) is n_E x nt
        squeezed_prod_b_E = squeeze_dim(prod(b_E_ts, 2), 2); % Result: n_E x nt
        p_ts(E_idx, :) = p_ts(E_idx, :) .* squeezed_prod_b_E;
    end

    % Apply STD effect to I neurons (vectorized)
    if params.n_I > 0 && params.n_b_I > 0 && ~isempty(b_I_ts)
        % b_I_ts is n_I x n_b_I x nt
        % prod(b_I_ts, 2) is n_I x 1 x nt
        % squeeze_dim(prod(b_I_ts, 2), 2) is n_I x nt
        squeezed_prod_b_I = squeeze_dim(prod(b_I_ts, 2), 2); % Result: n_I x nt
        p_ts(I_idx, :) = p_ts(I_idx, :) .* squeezed_prod_b_I;
    end
end