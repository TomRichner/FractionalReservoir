%% replot_all_analyses.m
% Master replot script -- regenerates figures for all three analyses in a
% prior run_all_analyses.m run, without re-running any simulations.
%
% Calls (in order):
%   1. replot_tau_sensitivity
%   2. replot_sensitivity      (+ assemble_sensitivity_figure)
%   3. replot_param_space_analysis
%
% Each function writes into its own replot_<...>_<dt>/ subfolder under
% data_root, so originals are never overwritten.

%% Setup
setup_paths();
close all;

fprintf('=== Starting Replot of All Analyses ===\n');
fprintf('Start time: %s\n\n', datetime('now'));
tic;

%% Configuration -- point this at the run_all_<dt> folder you want to replot
% Derive project_root from setup_paths.m (at the repo root) so this tolerates living
% in a subdirectory such as scripts/run_all_analyses/replot/.
project_root = fileparts(which('setup_paths'));
data_root    = fullfile(project_root, 'data', 'param_space', ...
    'run_all_mar_02_26_17_12');

% LLE histogram y-range applied to the (1D and tau) sensitivity replots.
% param_space replot does not use this (psa.plot has no hist_range arg).
lle_hist_range = [-1, 1];

fprintf('Replotting from: %s\n', data_root);
fprintf('LLE hist_range:  [%g, %g]\n\n', lle_hist_range(1), lle_hist_range(2));

%% 1. Tau Sensitivity
fprintf('========================================\n');
fprintf('[1/3] Replotting Tau Sensitivity...\n');
fprintf('========================================\n');
replot_tau_sensitivity(data_root, lle_hist_range);

%% 2. 1D Sensitivity (n, f, S_c)
fprintf('========================================\n');
fprintf('[2/3] Replotting Sensitivity...\n');
fprintf('========================================\n');
sens_replot_dir = replot_sensitivity(data_root, lle_hist_range);

% Assemble per-param LLE sensitivity figs into a single stacked figure
assemble_sensitivity_figure(sens_replot_dir, 'LLE');

%% 3. Parameter Space Analysis
fprintf('========================================\n');
fprintf('[3/3] Replotting Parameter Space Analysis...\n');
fprintf('========================================\n');
replot_param_space_analysis(data_root);

%% Summary
fprintf('========================================\n');
fprintf('=== All Replots Complete ===\n');
fprintf('Total runtime: %.2f minutes\n', toc/60);
fprintf('End time: %s\n', datetime('now'));
fprintf('Source data root: %s\n', data_root);
fprintf('========================================\n');
