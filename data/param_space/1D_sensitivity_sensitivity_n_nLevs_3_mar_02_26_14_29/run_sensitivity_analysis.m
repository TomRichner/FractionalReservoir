% run_sensitivity_analysis.m
% Sensitivity analysis using ParamSpaceAnalysis2 in 1D mode
%
% Sweeps individual parameters across multiple levels and repetitions,
% comparing different adaptation conditions. Each parameter is swept
% independently (one PSA run per parameter).
%
% See also: ParamSpaceAnalysis2, SRNNModel2

clear;
clc;
close all;

%% Setup paths
setup_paths();

%% Analysis Configuration
n_levels = 3;      % Number of parameter values to test
n_reps = 3;        % Number of repetitions per level (for statistics)
note = 'sensitivity';

% Parameters to sweep: {param_name, [min, max]}
params_to_sweep = {
    'n',              [50, 250];
    'f',              [0.4, 0.6];
};

%% Run sensitivity analysis for each parameter
all_output_dirs = {};

for p_idx = 1:size(params_to_sweep, 1)
    param_name = params_to_sweep{p_idx, 1};
    param_range = params_to_sweep{p_idx, 2};

    fprintf('\n========================================\n');
    fprintf('=== Sensitivity: %s [%.3g, %.3g] (%d/%d) ===\n', ...
        param_name, param_range(1), param_range(2), p_idx, size(params_to_sweep, 1));
    fprintf('========================================\n');

    % Create PSA for this parameter
    psa = ParamSpaceAnalysis2(...
        'n_levels', n_levels, ...
        'batch_size', 25, ...
        'note', sprintf('%s_%s', note, param_name), ...
        'randomize_order', false, ...
        'verbose', true ...
        );
    psa.folder_prefix = '1D_sensitivity';
    if exist('master_output_dir', 'var')
        psa.output_dir = master_output_dir;
    end

    % Override defaults to ensure indegree <= min(n)
    psa.model_defaults.indegree = 50;

    % Add the swept parameter and reps
    psa.add_grid_parameter(param_name, param_range);
    psa.add_grid_parameter('reps', 1:n_reps);

    % Run
    psa.run();

    % Copy this script for reproducibility
    copyfile([mfilename('fullpath') '.m'], psa.output_dir);

    % Plot sensitivity heatmaps
    psa.plot_sensitivity('metric', 'LLE');
    psa.plot_sensitivity('metric', 'mean_rate');

    % Save PSA object
    save(fullfile(psa.output_dir, 'psa_object.mat'), 'psa');

    all_output_dirs{end+1} = psa.output_dir; %#ok<SAGROW>
end

%% Summary
fprintf('\n=== Sensitivity Analysis Complete ===\n');
fprintf('Parameters analyzed:\n');
for p_idx = 1:size(params_to_sweep, 1)
    fprintf('  %s: %s\n', params_to_sweep{p_idx, 1}, all_output_dirs{p_idx});
end
fprintf('Done!\n');
