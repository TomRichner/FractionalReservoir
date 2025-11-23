function [W, M, G, Z] = create_paired_W_matrix(params)
% CREATE_PAIRED_W_MATRIX Creates a connectivity matrix with paired E-I structure
%
% The structure enforces that for every Excitatory neuron i, there is an
% Inhibitory neuron i (index i + n_E) that it drives.
% The Inhibitory neuron then projects to the same targets as the Excitatory neuron,
% providing feed-forward inhibition.
%
% Usage:
%   [W, M, G, Z] = create_paired_W_matrix(params)
%
% Inputs:
%   params.n_E - Number of Excitatory neurons
%   params.n_I - Number of Inhibitory neurons (must equal n_E)
%   params.n   - Total number of neurons (must be n_E + n_I)

    n_E = params.n_E;
    n_I = params.n_I;
    n = params.n;

    % Validate inputs
    if n_E ~= n_I
        error('create_paired_W_matrix:UnequalPopulations', ...
              'Paired connectivity requires equal E and I populations (n_E=%d, n_I=%d)', n_E, n_I);
    end
    if n ~= (n_E + n_I)
        error('create_paired_W_matrix:SizeMismatch', ...
              'params.n (%d) must equal params.n_E + params.n_I (%d)', n, n_E + n_I);
    end

    % Initialize matrices
    W = zeros(n);
    
    %% 1. E -> E connections (Sparse Random)
    % Generate random E-to-E connectivity
    p_connect = params.indegree / n_E; % Maintain similar indegree from E population
    
    % Mean structure (if mu_E is used)
    M_EE = params.mu_E * ones(n_E);
    
    % Gaussian perturbation
    G_EE = params.G_stdev * randn(n_E);
    
    % Base weights
    W_EE = M_EE + G_EE;
    
    % Sparsity mask (1 means removed)
    Z_EE = rand(n_E) > p_connect;
    W_EE(Z_EE) = 0;
    
    %% 2. E -> I (Driving Pairs)
    % Each E_i drives I_i strongly.
    % This is a 1-to-1 connection.
    w_drive = 5.0; % Strong driving force
    W_IE = eye(n_E) * w_drive;
    
    %% 3. I -> E (Feed-forward Inhibition)
    % I_i projects to the same targets as E_i.
    % Ratio of inhibition strength relative to excitation
    beta = 1.2; 
    W_EI = -abs(W_EE) * beta;
    
    %% 4. I -> I (Recurrent Inhibition)
    M_II = params.mu_I * ones(n_I);
    G_II = params.G_stdev * randn(n_I);
    W_II = M_II + G_II;
    Z_II = rand(n_I) > p_connect;
    W_II(Z_II) = 0;
    
    % Ensure negative
    W_II = -abs(W_II);

    %% Assemble W
    W(1:n_E, 1:n_E)       = W_EE;
    W(1:n_E, (n_E+1):n)   = W_EI;
    W((n_E+1):n, 1:n_E)   = W_IE;
    W((n_E+1):n, (n_E+1):n) = W_II;
    
    %% Assemble Z (1 where connection is MISSING)
    % Note: For 1-to-1 blocks, the off-diagonals are "missing".
    
    Z = ones(n); % Start assuming everything is missing
    
    % E->E block
    Z(1:n_E, 1:n_E) = Z_EE; 
    
    % E->I block (Diagonal exists, so off-diagonals are missing)
    Z((n_E+1):n, 1:n_E) = ~eye(n_E);
    
    % I->E block (Shadow of E->E)
    % If W_EI matches W_EE sparsity, then Z_EI = Z_EE
    Z(1:n_E, (n_E+1):n) = Z_EE; 
    
    % I->I block
    Z((n_E+1):n, (n_E+1):n) = Z_II;

    %% Row-mean centering (optional)
    % Center E->E
    nonzero_mask = W_EE ~= 0;
    row_means = sum(W_EE, 2) ./ max(1, sum(nonzero_mask, 2));
    W_EE = W_EE - bsxfun(@times, row_means, nonzero_mask);
    W(1:n_E, 1:n_E) = W_EE;

    % Center I->I
    nonzero_mask = W_II ~= 0;
    row_means = sum(W_II, 2) ./ max(1, sum(nonzero_mask, 2));
    W_II = W_II - bsxfun(@times, row_means, nonzero_mask);
    W((n_E+1):n, (n_E+1):n) = W_II;
    
    % Update W with centered blocks
    W(1:n_E, 1:n_E)       = W_EE;
    W((n_E+1):n, (n_E+1):n) = W_II;
    
    % Return auxiliary matrices (M, G are placeholders here)
    M = zeros(n); 
    G = zeros(n);

end
