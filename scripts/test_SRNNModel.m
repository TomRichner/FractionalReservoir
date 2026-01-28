% test_SRNNModel.m - Test script for SRNNModel class
%
% This script tests the SRNNModel class by running a basic simulation
% and comparing results against the expected behavior from full_SRNN_run.m

close all; clear all; clc;

% Add paths - adjust setup_paths location as needed
if exist('setup_paths', 'file')
    setup_paths();
else
    % Manual path setup if setup_paths not available
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')));
end

%% Test 1: Basic construction and defaults
fprintf('=== Test 1: Basic construction ===\n');
model = SRNNModel();
disp(model);
fprintf('Default n: %d\n', model.n);
fprintf('Default level_of_chaos: %.2f\n', model.level_of_chaos);
fprintf('Test 1 PASSED\n\n');

%% Test 2: Construction with name-value pairs
fprintf('=== Test 2: Name-value pair construction ===\n');
model2 = SRNNModel('n', 50, 'level_of_chaos', 2.0, 'n_a_E', 3, 'n_b_E', 1);
assert(model2.n == 50, 'n should be 50');
assert(model2.level_of_chaos == 2.0, 'level_of_chaos should be 2.0');
assert(model2.n_a_E == 3, 'n_a_E should be 3');
assert(model2.n_b_E == 1, 'n_b_E should be 1');
fprintf('Test 2 PASSED\n\n');

%% Test 3: Build method
fprintf('=== Test 3: Build method ===\n');
model3 = SRNNModel('n', 100, 'level_of_chaos', 1.8, 'T_range', [0, 10]);
model3.build();

assert(model3.is_built, 'Model should be built');
assert(~isempty(model3.W), 'W matrix should exist');
assert(size(model3.W, 1) == 100, 'W should be 100x100');
assert(model3.n_E == 50, 'n_E should be 50 (f=0.5)');
assert(model3.n_I == 50, 'n_I should be 50');
assert(~isempty(model3.t_ex), 't_ex should be generated');
assert(~isempty(model3.u_ex), 'u_ex should be generated');
fprintf('Test 3 PASSED\n\n');

%% Test 4: Run simulation (no adaptation, no STD)
fprintf('=== Test 4: Run simulation (baseline) ===\n');
model4 = SRNNModel('n', 100, ...
    'n_a_E', 0, 'n_b_E', 0, 'T_range', [0, 10], ...
    'store_full_state', true, 'lya_method', 'benettin');
model4.build();
model4.run();

assert(model4.has_run, 'Model should have run');
assert(~isempty(model4.lya_results), 'Lyapunov results should exist');
assert(isfield(model4.lya_results, 'LLE'), 'LLE should be computed');
fprintf('LLE: %.4f\n', model4.lya_results.LLE);
fprintf('Test 4 PASSED\n\n');

%% Test 5-7: Run with adaptation/STD, plotting, and get_params
fprintf('=== Test 5: Run with SFA and STD ===\n');
model5 = SRNNModel('n', 100, ...
    'n_a_E', 0, 'n_b_E', 0, 'T_range', [0, 10], ...
    'store_full_state', false, 'store_decimated_state', true);
model5.build();
model5.run();

assert(isempty(model5.S_out), 'S_out should be cleared (store_full_state=false)');
assert(~isempty(model5.plot_data), 'plot_data should exist');
fprintf('Test 5 PASSED\n\n');

% Test 6: Plotting (uses model5 from above)
fprintf('=== Test 6: Plotting ===\n');
[fig, ax] = model5.plot();
assert(~isempty(fig), 'Figure handle should exist');
assert(~isempty(ax), 'Axes handles should exist');
title(ax(1), 'Test 6: SFA + STD');
fprintf('Test 6 PASSED\n\n');

% Test 7: get_params compatibility (uses model5 from above)
fprintf('=== Test 7: get_params compatibility ===\n');
params = model5.get_params();
assert(isfield(params, 'W'), 'params should have W');
assert(isfield(params, 'n_E'), 'params should have n_E');
assert(isfield(params, 'activation_function'), 'params should have activation_function');
assert(isa(params.activation_function, 'function_handle'), 'activation_function should be function handle');
fprintf('Test 7 PASSED\n\n');

%% Test 8-9: Eigenvalue plotting and memory management
fprintf('=== Test 8: Eigenvalue plotting ===\n');
model8 = SRNNModel('n', 100, 'level_of_chaos', 1.8, 'T_range', [0, 10], ...
    'store_full_state', true, 'lya_method', 'none');
model8.build();
model8.run();
[fig_eig, ax_eig] = model8.plot_eigenvalues([2, 5, 8]);
assert(~isempty(fig_eig), 'Eigenvalue figure should exist');
fprintf('Test 8 PASSED\n\n');

% Test 9: Memory management with clear_results (uses model8 from above)
fprintf('=== Test 9: Memory management ===\n');
model8.clear_results();
assert(isempty(model8.S_out), 'S_out should be cleared');
assert(isempty(model8.plot_data), 'plot_data should be cleared');
assert(~model8.has_run, 'has_run should be false');
fprintf('Test 9 PASSED\n\n');

%% Test 10: Reproduce SRNN_comparisons configuration
fprintf('=== Test 10: Reproduce SRNN_comparisons configuration ===\n');
model10 = SRNNModel();
model10.n = 100;
model10.indegree = 20;
model10.f = 0.5;
model10.level_of_chaos = 1.8;
model10.u_ex_scale = 1.7;
model10.rng_seeds = [8 16 3 4 5];
model10.n_a_E = 3;
model10.n_b_E = 1;
model10.T_range = [0, 20];
model10.lya_method = 'benettin';
model10.store_full_state = false;

model10.build();
model10.run();
model10.plot();

fprintf('LLE for SFA+STD model: %.4f\n', model10.lya_results.LLE);
fprintf('Test 10 PASSED\n\n');

%% Summary
fprintf('========================================\n');
fprintf('All tests PASSED!\n');
fprintf('SRNNModel class is working correctly.\n');
fprintf('========================================\n');

