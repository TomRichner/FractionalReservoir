%% run_all_analyses.m
% Master script to run all analysis pipelines in sequence
%
% This script executes the following analyses in order:
%   1. Tau sensitivity analysis (tau_a and tau_b parameter sweeps)
%   2. General sensitivity analysis (SFA/STD parameter exploration)
%   3. Parameter space analysis
%
% Each analysis saves its results to disk automatically.

%% Setup
fprintf('=== Starting All Analyses ===\n');
fprintf('Start time: %s\n\n', datetime('now'));
tic;

%% 1. Tau Sensitivity Analysis
fprintf('========================================\n');
fprintf('[1/3] Running Tau Sensitivity Analysis...\n');
fprintf('========================================\n');
run_tau_sensitivity_analysis;

%% 2. General Sensitivity Analysis
fprintf('========================================\n');
fprintf('[2/3] Running Sensitivity Analysis...\n');
fprintf('========================================\n');
run_sensitivity_analysis;

%% 3. Parameter Space Analysis
fprintf('========================================\n');
fprintf('[3/3] Running Parameter Space Analysis...\n');
fprintf('========================================\n');
run_param_space_analysis2;

%% Summary
fprintf('========================================\n');
fprintf('=== All Analyses Complete ===\n');
fprintf('Total runtime: %.2f minutes\n', toc/60);
fprintf('End time: %s\n', datetime('now'));
fprintf('========================================\n');
