% run_sensitivity_analysis.m
% Example script demonstrating how to use the SensitivityAnalysis class
%
% This script performs a sensitivity analysis on the SRNN model by varying
% specified parameters across multiple levels and repetitions, comparing
% different adaptation conditions.
%
% See also: SensitivityAnalysis, SRNNModel

clear;
clc;
close all;

%% Setup paths
setup_paths();

%% Create SensitivityAnalysis object
% Configure the number of levels and repetitions
% (reduce these for faster testing)
sa = SensitivityAnalysis(...
    'n_levels', 11, ...       % Number of parameter values to test
    'n_reps', 10, ...         % Number of repetitions per value (for statistics)
    'note', 'example_run', ...% Optional note for output folder naming
    'verbose', true ...       % Print progress during execution
);

%% Add parameters to sweep
% Specify which parameters to vary and their ranges
% You can add any property of SRNNModel

% Network structure parameters
sa.add_parameter('level_of_chaos', [0.5, 2.5]);  % Abscissa scaling
sa.add_parameter('n', [50, 150]);                 % Number of neurons
% sa.add_parameter('indegree', [10, 40]);         % Expected in-degree

% Dynamics parameters
% sa.add_parameter('tau_d', [0.05, 0.2]);         % Dendritic time constant
% sa.add_parameter('u_ex_scale', [1.0, 3.0]);     % External input scaling

% Adaptation parameters
% sa.add_parameter('c_E', [0.01, 0.1]);           % SFA strength for E neurons

%% Configure model defaults (optional)
% Set any SRNNModel properties that should be constant across all runs
sa.model_defaults.T_range = [0, 30];       % Shorter simulation for faster runs
sa.model_defaults.fs = 200;                 % Sampling frequency
sa.model_defaults.lya_method = 'benettin';  % Lyapunov computation method

%% Configure conditions (optional)
% By default, all four adaptation conditions are tested:
%   - no_adaptation: n_a_E=0, n_b_E=0
%   - sfa_only:      n_a_E=3, n_b_E=0
%   - std_only:      n_a_E=0, n_b_E=1
%   - sfa_and_std:   n_a_E=3, n_b_E=1
%
% To use only a subset of conditions, uncomment and modify below:
%
% sa.set_conditions({
%     struct('name', 'no_adaptation', 'n_a_E', 0, 'n_b_E', 0),
%     struct('name', 'sfa_and_std',   'n_a_E', 3, 'n_b_E', 1)
% });

%% Run the sensitivity analysis
% This will loop over all conditions and parameters using parfor
% Results are automatically saved to disk during execution
sa.run();

%% Plot results
% Generate heatmaps showing LLE distribution across parameter values
sa.plot('metric', 'LLE');

% Optionally plot mean firing rate
sa.plot('metric', 'mean_rate');

%% Display summary
fprintf('\n=== Analysis Summary ===\n');
fprintf('Output directory: %s\n', sa.output_dir);
fprintf('Parameters analyzed: %s\n', strjoin(fieldnames(sa.param_ranges), ', '));
fprintf('Conditions: %s\n', strjoin(cellfun(@(c) c.name, sa.conditions, 'UniformOutput', false), ', '));

%% Optional: Access raw results programmatically
% The results are stored in sa.results structure:
%
%   sa.results.<condition_name>.<param_name>.results_reshaped
%   sa.results.<condition_name>.<param_name>.metadata
%
% Example: Get LLE values for 'level_of_chaos' under 'sfa_only' condition:
%
%   res = sa.results.sfa_only.level_of_chaos.results_reshaped;
%   LLE_values = cellfun(@(r) r.LLE, res);
%   param_levels = sa.results.sfa_only.level_of_chaos.metadata.param_levels;
%
%   figure;
%   plot(param_levels, mean(LLE_values, 2), 'o-');
%   xlabel('level\_of\_chaos');
%   ylabel('Mean LLE');

fprintf('\nDone! Results saved to: %s\n', sa.output_dir);

