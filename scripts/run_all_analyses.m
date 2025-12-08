%% run_all_analyses.m
% Master script to run all analysis pipelines in sequence
%
% This script executes the following analyses in order:
%   1. Tau sensitivity analysis (tau_a and tau_b parameter sweeps)
%   2. General sensitivity analysis (SFA/STD parameter exploration)
%   3. Paired-pulse mutual information analysis
%   4. Parameter space analysis
%
% Each analysis saves its results to disk automatically.

%% Setup
clear; clc;
fprintf('=== Starting All Analyses ===\n');
fprintf('Start time: %s\n\n', datetime('now'));
total_start = tic;

%% 1. Tau Sensitivity Analysis
fprintf('========================================\n');
fprintf('[1/4] Running Tau Sensitivity Analysis...\n');
fprintf('========================================\n');
step_start = tic;
run_tau_sensitivity_analysis;
fprintf('Tau Sensitivity Analysis completed in %.2f minutes.\n\n', toc(step_start)/60);

%% 2. General Sensitivity Analysis
fprintf('========================================\n');
fprintf('[2/4] Running Sensitivity Analysis...\n');
fprintf('========================================\n');
step_start = tic;
run_sensitivity_analysis;
fprintf('Sensitivity Analysis completed in %.2f minutes.\n\n', toc(step_start)/60);

%% 3. Paired-Pulse Mutual Information Analysis
fprintf('========================================\n');
fprintf('[3/4] Running Paired-Pulse MI Analysis...\n');
fprintf('========================================\n');
step_start = tic;
run_paired_pulse_mi;
fprintf('Paired-Pulse MI Analysis completed in %.2f minutes.\n\n', toc(step_start)/60);

%% 4. Parameter Space Analysis
fprintf('========================================\n');
fprintf('[4/4] Running Parameter Space Analysis...\n');
fprintf('========================================\n');
step_start = tic;
run_param_space_analysis;
fprintf('Parameter Space Analysis completed in %.2f minutes.\n\n', toc(step_start)/60);

%% Summary
fprintf('========================================\n');
fprintf('=== All Analyses Complete ===\n');
fprintf('Total runtime: %.2f minutes\n', toc(total_start)/60);
fprintf('End time: %s\n', datetime('now'));
fprintf('========================================\n');
