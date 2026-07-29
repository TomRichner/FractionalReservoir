% test_SRNNCellTypes_defaults.m - Two-type counterpart to test_SRNN2_defaults.
%
% Reproduces the SRNNModel2 default parameterization using the generalized
% cell-type API: n=300, indegree=100, equal E/I fractions, three SFA
% timescales on E, one STD timescale on E, T_range=[0,50], and Benettin LLE.

close all; clear; clc;

%% Match SRNNModel2 connectivity defaults
n = 400;
indegree = 100;
alpha = indegree / n;
F = 1 / sqrt(n * alpha * (2 - alpha));

%% Use the rounded piecewise sigmoid
S_a = 0.9;
S_c = 0.35;
phi = @(x) SRNNCellTypes.piecewiseSigmoid(x, S_a, S_c);
phi_derivative = @(x) SRNNCellTypes.piecewiseSigmoidDerivative(x, S_a, S_c);

%% Create the equivalent two-cell-type model
model = SRNNCellTypes( ...
    'n', n, ...
    'indegree', indegree, ...
    'n_cellTypes', 2, ...
    'cell_type_names', {'E', 'I'}, ...
    'f', [0.5 0.5], ...
    'mu_tilde', F .* [3 -4], ...
    'sigma_tilde', F .* [1 1], ...
    'n_a', [3 0], ...
    'n_b', [1 0], ...
    'S_a', S_a, ...
    'S_c', S_c, ...
    'activation_function', phi, ...
    'activation_function_derivative', phi_derivative);

% Match SRNNModel2's default type-selective sparse stimulus.
model.input_config.step_density = struct('E', 0.15, 'I', 0);

%% Build, run, and plot
model.build();
model.run();
model.plot();
model.plot_celltypes();