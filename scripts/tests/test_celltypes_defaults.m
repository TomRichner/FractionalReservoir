% test_celltypes_defaults.m - Run SRNNModelCellTypes with Campagnola data, n=50
%
% Cell-typed analogue of test_SRNN2_defaults.m. Constructs SRNNModelCellTypes at
% defaults (SFA + STD + STF on), which loads the 4 cell types (Pyr/Pvalb/Sst/Vip)
% from the committed Campagnola matrices: block connectivity, per-type STD/STF
% release-model parameters, the per-type fitted SFA tau_a (sfa_tau_per_type.csv),
% and per-type n_a (non-adapting Pvalb carries no SFA state). n=50 neurons.

close all; clear; clc;

% Add paths
setup_paths();

%% Create model: 50 neurons, all 4 cell types from Campagnola data
% Piecewise (hard) sigmoid nonlinearity, S_a=0.9 (shape), S_c=0.35 (center).
model = SRNNModelCellTypes('n', 300, 'lya_method', 'benettin', ...
    'S_a', 0.9, 'S_c', 0.35, ...
    'activation_function', @(x) SRNNModelBase.piecewiseSigmoid(x, 0.9, 0.35), ...
    'activation_function_derivative', @(x) SRNNModelBase.piecewiseSigmoidDerivative(x, 0.9, 0.35));

%% Build, run, and plot
model.build();

% Report the data-derived per-type SFA (fitted tau_a + which types adapt)
fprintf('\nPer-type SFA loaded from Campagnola data:\n');
for T = 1:model.n_types
    j = find(model.type_of == T, 1);
    fprintf('  %-6s adaptation_index=%.4f  adapting=%d  tau_a=%.3f s  tau_d=%.4f s\n', ...
        model.type_names{T}, model.adapt_index(T), ismember(j, model.ad_idx), ...
        model.tau_a_type(T), model.tau_d_type(T));
end

model.run();
model.plot();
