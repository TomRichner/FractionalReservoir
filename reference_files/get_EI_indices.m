function [E_indices, I_indices, n_E, n_I] = get_EI_indices(EI_vec)
% get_EI_indices - Get indices and counts of excitatory and inhibitory neurons
%
% Inputs:
%   EI_vec   - Vector (n x 1) with 1 for excitatory, -1 for inhibitory neurons.
%
% Outputs:
%   E_indices - Indices of excitatory neurons.
%   I_indices - Indices of inhibitory neurons.
%   n_E       - Number of excitatory neurons.
%   n_I       - Number of inhibitory neurons.

    E_indices = find(EI_vec == 1);
    I_indices = find(EI_vec == -1);
    n_E = numel(E_indices);
    n_I = numel(I_indices);
end