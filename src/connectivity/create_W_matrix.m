function [W, M, G, Z] = create_W_matrix(params)
% create_W_matrix - Generate connectivity matrix for SRNN
%
% Syntax:
%   [W, M, G, Z] = create_W_matrix(params)
%
% Description:
%   Creates a sparse connectivity matrix W with structured excitatory/inhibitory
%   balance and row-mean centering. The matrix is composed of:
%   - M: Mean connectivity structure (E-to-all and I-to-all)
%   - G: Gaussian random perturbations
%   - Z: Binary sparsification mask
%   - Row-mean centering applied to non-zero elements
%
% Inputs:
%   params - Struct containing:
%            .n         - Total number of neurons
%            .n_E       - Number of excitatory neurons
%            .n_I       - Number of inhibitory neurons
%            .mu_E      - Mean excitatory connection strength
%            .mu_I      - Mean inhibitory connection strength
%            .G_stdev   - Standard deviation of Gaussian perturbations
%            .indegree  - Expected in-degree (number of inputs per neuron)
%
% Outputs:
%   W - n x n connectivity matrix (sparse, row-mean centered)
%   M - n x n mean structure matrix
%   G - n x n Gaussian perturbation matrix
%   Z - n x n binary sparsification mask (1 = connection removed)
%
% Example:
%   params.n = 100;
%   params.n_E = 50;
%   params.n_I = 50;
%   params.mu_E = 1;
%   params.mu_I = -1;
%   params.G_stdev = 2;
%   params.indegree = 15;
%   [W, M, G, Z] = create_W_matrix(params);

    % Create mean structure matrix M
    % First n_E columns get mu_E (excitatory), remaining get mu_I (inhibitory)
    M = [params.mu_E .* ones(params.n, params.n_E), ...
         params.mu_I .* ones(params.n, params.n_I)];
    
    % Add Gaussian perturbations
    G = params.G_stdev * randn(params.n, params.n);
    W = M + G;
    
    % Apply sparsification
    d = params.indegree / params.n;  % density
    Z = rand(params.n, params.n) > d;  % 1 where connections are removed
    W(Z) = 0;
    
    % Zero out row mean of non-zero elements
    nonzero_mask = ~Z;
    row_counts = sum(nonzero_mask, 2);
    row_sums = sum(W, 2);
    row_means = zeros(size(row_sums));
    valid_rows = row_counts > 0;
    row_means(valid_rows) = row_sums(valid_rows) ./ row_counts(valid_rows);
    W = W - bsxfun(@times, row_means, nonzero_mask);
end

