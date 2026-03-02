% run_tau_sensitivity_analysis.m
% Sensitivity analysis for tau_a_E(end) and tau_b_E_rec parameters
%
% Analyzes the effect of the maximum SFA time constant (tau_a_E_max)
% and STD recovery time constant (tau_b_E_rec) on SRNN dynamics.
%
% Uses ParamSpaceAnalysis2 with add_vector_parameter for tau_a_E
% and add_grid_parameter for tau_b_E_rec.
%
% See also: ParamSpaceAnalysis2, SRNNModel2

clear;
clc;
close all;

%% Setup paths
setup_paths();

%% Analysis Configuration
n_levels = 25;
n_reps = 100;
note = 'tau_timescales';

% Condition: SFA + STD (n_a_E=3, n_b_E=1)
condition = {struct('name', 'sfa_and_std', 'n_a_E', 3, 'n_b_E', 1)};

%% 1. tau_a_E(end) sweep — vector parameter
fprintf('\n========================================\n');
fprintf('=== Tau Sensitivity: tau_a_E(end) [5, 60] ===\n');
fprintf('========================================\n');

psa_tau_a = ParamSpaceAnalysis2(...
    'n_levels', n_levels, ...
    'batch_size', 25, ...
    'note', sprintf('%s_tau_a_E_max', note), ...
    'randomize_order', false, ...
    'verbose', true ...
    );
psa_tau_a.folder_prefix = 'tau_sensitivity';
if exist('master_output_dir', 'var')
    psa_tau_a.output_dir = master_output_dir;
end

psa_tau_a.set_conditions(condition);

% tau_a_E is a vector of length n_a_E=3, logspaced from 0.25 to max
% We sweep the max (last element) from 5 to 60
psa_tau_a.add_vector_parameter('tau_a_E', ...
    'vary_element', 'last', ...
    'fixed_value', 0.25, ...
    'vary_range', [5, 60], ...
    'n_elements', 3, ...
    'spacing', 'log', ...
    'level_spacing', 'linear');

psa_tau_a.add_grid_parameter('reps', 1:n_reps);

psa_tau_a.run();

% Copy script for reproducibility
copyfile([mfilename('fullpath') '.m'], psa_tau_a.output_dir);

% Plot
psa_tau_a.plot_sensitivity('metric', 'LLE', 'hist_range', [-0.3, 0.1]);
psa_tau_a.plot_sensitivity('metric', 'mean_rate');

save(fullfile(psa_tau_a.output_dir, 'psa_object.mat'), 'psa_tau_a');

%% 2. tau_b_E_rec sweep — scalar parameter
fprintf('\n========================================\n');
fprintf('=== Tau Sensitivity: tau_b_E_rec [5, 60] ===\n');
fprintf('========================================\n');

psa_tau_b = ParamSpaceAnalysis2(...
    'n_levels', n_levels, ...
    'batch_size', 25, ...
    'note', sprintf('%s_tau_b_E_rec', note), ...
    'randomize_order', false, ...
    'verbose', true ...
    );
psa_tau_b.folder_prefix = 'tau_sensitivity';
if exist('master_output_dir', 'var')
    psa_tau_b.output_dir = master_output_dir;
end

psa_tau_b.set_conditions(condition);

psa_tau_b.add_grid_parameter('tau_b_E_rec', [5, 60]);
psa_tau_b.add_grid_parameter('reps', 1:n_reps);

psa_tau_b.run();

copyfile([mfilename('fullpath') '.m'], psa_tau_b.output_dir);

psa_tau_b.plot_sensitivity('metric', 'LLE', 'hist_range', [-0.3, 0.1]);
psa_tau_b.plot_sensitivity('metric', 'mean_rate');

save(fullfile(psa_tau_b.output_dir, 'psa_object.mat'), 'psa_tau_b');

%% Summary
fprintf('\n========================================\n');
fprintf('=== Tau Sensitivity Analysis Complete ===\n');
fprintf('tau_a_E results: %s\n', psa_tau_a.output_dir);
fprintf('tau_b_E_rec results: %s\n', psa_tau_b.output_dir);
fprintf('========================================\n');
beep; pause(0.5); beep; pause(0.2); beep
