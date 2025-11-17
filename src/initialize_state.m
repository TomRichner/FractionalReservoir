function S0 = initialize_state(params)
% initialize_state - Initialize state vector for SRNN
%
% Syntax:
%   S0 = initialize_state(params)
%
% Description:
%   Initializes the complete state vector for an SRNN with adaptation and
%   short-term depression. The state is organized as:
%   S0 = [a_E(:); a_I(:); b_E(:); b_I(:); x(:)]
%   
%   - a_E, a_I: Adaptation variables (initialized to 0)
%   - b_E, b_I: STD variables (initialized to 1, no depression)
%   - x: Dendritic states (initialized to small random values)
%
% Inputs:
%   params - Struct containing:
%            .n       - Total number of neurons
%            .n_E     - Number of excitatory neurons
%            .n_I     - Number of inhibitory neurons
%            .n_a_E   - Number of adaptation timescales for E neurons
%            .n_a_I   - Number of adaptation timescales for I neurons
%            .n_b_E   - Number of STD timescales for E neurons (0 or 1)
%            .n_b_I   - Number of STD timescales for I neurons (0 or 1)
%
% Outputs:
%   S0 - Initial state vector (column vector)
%
% Example:
%   params.n = 100;
%   params.n_E = 50;
%   params.n_I = 50;
%   params.n_a_E = 1;
%   params.n_a_I = 0;
%   params.n_b_E = 1;
%   params.n_b_I = 0;
%   S0 = initialize_state(params);

    % Initialize adaptation states for E neurons
    a0_E = [];
    if params.n_a_E > 0
        a0_E = zeros(params.n_E * params.n_a_E, 1);
    end
    
    % Initialize adaptation states for I neurons
    a0_I = [];
    if params.n_a_I > 0
        a0_I = zeros(params.n_I * params.n_a_I, 1);
    end
    
    % Initialize STD states for E neurons (b variables start at 1.0, no depression)
    b0_E = [];
    if params.n_b_E > 0
        b0_E = ones(params.n_E * params.n_b_E, 1);
    end
    
    % Initialize STD states for I neurons
    b0_I = [];
    if params.n_b_I > 0
        b0_I = ones(params.n_I * params.n_b_I, 1);
    end
    
    % Initialize dendritic states (small random values)
    x0 = 0.1 .* randn(params.n, 1);
    
    % Pack state vector: [a_E; a_I; b_E; b_I; x]
    S0 = [a0_E; a0_I; b0_E; b0_I; x0];
end

