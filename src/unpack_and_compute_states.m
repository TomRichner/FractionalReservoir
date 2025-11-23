function [x, a, b, r, br] = unpack_and_compute_states(S_out, params, a_zeros_b_ones)
% unpack_and_compute_states - Unpack state vector and compute dependent variables
%
% Syntax:
%   [x, a, b, r, br] = unpack_and_compute_states(S_out, params)
%   [x, a, b, r, br] = unpack_and_compute_states(S_out, params, a_zeros_b_ones)
%
% Description:
%   Unpacks the state trajectory S_out into individual state variables,
%   splits them into excitatory and inhibitory components, and computes
%   the firing rate r and synaptic output br.
%
% Inputs:
%   S_out          - State trajectory (nt x N_sys_eqs) or column vector (N_sys_eqs x 1)
%                    State organization: [a_E(:); a_I(:); b_E(:); b_I(:); x(:)]
%   params         - Struct containing network parameters
%   a_zeros_b_ones - (Optional) If true, returns a as zeros(n, nt) and b as ones(n, nt)
%                    Default: false
%
% Outputs:
%   x  - Struct with fields .E and .I - dendritic states
%   a  - Struct with fields .E and .I - adaptation variables
%   b  - Struct with fields .E and .I - STD variables
%   r  - Struct with fields .E and .I - firing rates (raw, phi(x_eff))
%   br - Struct with fields .E and .I - synaptic output (b .* r)
%
%   If a_zeros_b_ones is true, all outputs are returned as simple arrays/zeros/ones.
%
% Example:
%   [x, a, b, r] = unpack_and_compute_states(S_out, params);
%   plot(t, x.E');  % Plot excitatory dendritic states
%
%   [x, a, b, r] = unpack_and_compute_states(S_out, params, true);
%   % a and b are now simple arrays for Jacobian computation

    % Handle optional parameter
    if nargin < 3
        a_zeros_b_ones = false;
    end
    
    nt = size(S_out, 1);
    current_idx = 0;
    
    %% Unpack adaptation states
    
    % Unpack adaptation states for E neurons (a_E)
    len_a_E = params.n_E * params.n_a_E;
    if len_a_E > 0
        a_E_ts = reshape(S_out(:, current_idx + (1:len_a_E))', params.n_E, params.n_a_E, nt);
    else
        a_E_ts = [];
    end
    current_idx = current_idx + len_a_E;
    
    % Unpack adaptation states for I neurons (a_I)
    len_a_I = params.n_I * params.n_a_I;
    if len_a_I > 0
        a_I_ts = reshape(S_out(:, current_idx + (1:len_a_I))', params.n_I, params.n_a_I, nt);
    else
        a_I_ts = [];
    end
    current_idx = current_idx + len_a_I;
    
    %% Unpack STD variables (b)
    
    % Unpack b states for E neurons (b_E)
    len_b_E = params.n_E * params.n_b_E;
    if len_b_E > 0
        b_E_ts = S_out(:, current_idx + (1:len_b_E))';  % n_E x nt
    else
        b_E_ts = [];
    end
    current_idx = current_idx + len_b_E;
    
    % Unpack b states for I neurons (b_I)
    len_b_I = params.n_I * params.n_b_I;
    if len_b_I > 0
        b_I_ts = S_out(:, current_idx + (1:len_b_I))';  % n_I x nt
    else
        b_I_ts = [];
    end
    current_idx = current_idx + len_b_I;
    
    %% Unpack dendritic states (x)
    x_ts = S_out(:, current_idx + (1:params.n))';  % n x nt
    
    %% Compute firing rates with adaptation and STD
    
    % Compute effective dendritic state (x_eff = x - c * sum(a))
    x_eff_ts = x_ts;  % n x nt
    
    % Apply adaptation effect to E neurons (scaled by c_E)
    if params.n_E > 0 && params.n_a_E > 0 && ~isempty(a_E_ts)
        % sum(a_E_ts, 2) is n_E x 1 x nt, need to squeeze to n_E x nt
        sum_a_E = squeeze(sum(a_E_ts, 2));  % n_E x nt
        if size(sum_a_E, 1) ~= params.n_E  % Handle case where nt=1
            sum_a_E = sum_a_E';
        end
        x_eff_ts(params.E_indices, :) = x_eff_ts(params.E_indices, :) - params.c_E * sum_a_E;
    end
    
    % Apply adaptation effect to I neurons (scaled by c_I)
    if params.n_I > 0 && params.n_a_I > 0 && ~isempty(a_I_ts)
        sum_a_I = squeeze(sum(a_I_ts, 2));  % n_I x nt
        if size(sum_a_I, 1) ~= params.n_I
            sum_a_I = sum_a_I';
        end
        x_eff_ts(params.I_indices, :) = x_eff_ts(params.I_indices, :) - params.c_I * sum_a_I;
    end
    
    % Apply STD effect (b multiplicative factor)
    b_ts = ones(params.n, nt);  % Initialize b = 1 for all neurons (no depression)
    if params.n_b_E > 0 && ~isempty(b_E_ts)
        b_ts(params.E_indices, :) = b_E_ts;
    end
    if params.n_b_I > 0 && ~isempty(b_I_ts)
        b_ts(params.I_indices, :) = b_I_ts;
    end
    
    % Compute firing rates: r = phi(x_eff) (raw rate)
    r_ts = params.activation_function(x_eff_ts);  % n x nt
    
    % Compute synaptic output: br = b .* r (presynaptically depressed)
    br_ts = b_ts .* r_ts; % n x nt
    
    %% Split into E and I components and package into structs
    
    % x: dendritic states
    x.E = x_ts(params.E_indices, :);  % n_E x nt
    x.I = x_ts(params.I_indices, :);  % n_I x nt
    
    % a: adaptation variables
    a.E = a_E_ts;  % n_E x n_a_E x nt (or empty)
    a.I = a_I_ts;  % n_I x n_a_I x nt (or empty)
    
    % b: STD variables
    if isempty(b_E_ts)
        b.E = ones(params.n_E, nt);  % Default to no depression
    else
        b.E = b_E_ts;  % n_E x nt
    end
    
    if isempty(b_I_ts)
        b.I = ones(params.n_I, nt);  % Default to no depression
    else
        b.I = b_I_ts;  % n_I x nt
    end
    
    % r: firing rates
    r.E = r_ts(params.E_indices, :);  % n_E x nt
    r.I = r_ts(params.I_indices, :);  % n_I x nt
    
    % br: synaptic output
    br.E = br_ts(params.E_indices, :);
    br.I = br_ts(params.I_indices, :);
    
    %% Override with zeros/ones if requested (for Jacobian computation)
    if a_zeros_b_ones
        % Return x as simple array instead of struct (n x nt)
        x = x_ts;
        
        % Replace a with zeros (n x nt)
        a = zeros(params.n, nt);
        
        % Replace b with ones (n x nt)
        b = ones(params.n, nt);
        
        % Return r as simple array instead of struct (n x nt)
        r = r_ts;
        
        % Return br as simple array (equal to r since b=1)
        br = r_ts;
    end
end

