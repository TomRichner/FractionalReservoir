%% run_all_analyses.m
% Master script to run all analysis pipelines in sequence
%
% This script executes the following analyses in order:
%   1. Tau sensitivity analysis (tau_a and tau_b parameter sweeps)
%   2. General sensitivity analysis (SFA/STD parameter exploration)
%   3. Parameter space analysis
%   4. Fig 2: Fraction excitatory analysis
%
% Each analysis saves its results to disk automatically.

%% Setup
setup_paths();

fprintf('=== Starting All Analyses ===\n');
fprintf('Start time: %s\n\n', datetime('now'));
tic;

% Create shared output directory for this run.
% Derive project_root from setup_paths.m (which lives in scripts/) rather than
% from this file's own location, so run_all_analyses.m can live in any
% subdirectory of scripts/ without breaking the output paths.
project_root = fileparts(fileparts(which('setup_paths')));
dt_str = lower(datestr(now, 'mmm_dd_yy_HH_MM')); %#ok<TNOW1,DATST>
master_output_dir = fullfile(project_root, 'data', 'param_space', ...
    sprintf('run_all_%s', dt_str));
if ~exist(master_output_dir, 'dir')
    mkdir(master_output_dir);
end
fprintf('Master output directory: %s\n\n', master_output_dir);

% Capture git provenance so the run can be retraced later
capture_git_provenance(master_output_dir, project_root);

% Control figure saving across all sub-scripts:
%   'save_all_figs'            - Override all scripts to save figures
%   'save_no_figs'             - Override all scripts to NOT save figures
%   'follow_scripts_save_figs' - Let each script use its own save_figs setting
master_save_figs = 'save_all_figs';

% Run mode for all three sub-analyses (controls n_levels / n_reps):
%   'fast'       - few levels/reps; finishes quickly (for testing)
%   'production' - full-size sweeps (for real results)
% Defaults to 'production'. To do a quick run WITHOUT editing any file, set the
% variable first in the console:   run_mode = 'fast'; run_all_analyses
if ~exist('run_mode', 'var'); run_mode = 'production'; end
fprintf('Run mode: %s\n\n', run_mode);

%% 1. Tau Sensitivity Analysis
fprintf('========================================\n');
fprintf('[1/4] Running Tau Sensitivity Analysis...\n');
fprintf('========================================\n');
run_tau_sensitivity_analysis;

%% 2. General Sensitivity Analysis
fprintf('========================================\n');
fprintf('[2/4] Running Sensitivity Analysis...\n');
fprintf('========================================\n');
run_sensitivity_analysis;

%% 3. Parameter Space Analysis
fprintf('========================================\n');
fprintf('[3/4] Running Parameter Space Analysis...\n');
fprintf('========================================\n');
run_param_space_analysis2;

% %% 4. Fig 2: Fraction Excitatory Analysis
% fprintf('========================================\n');
% fprintf('[4/4] Running Fig 2 Fraction Excitatory Analysis...\n');
% fprintf('========================================\n');
% Fig_2_fraction_excitatory_analysis;

%% Summary
fprintf('========================================\n');
fprintf('=== All Analyses Complete ===\n');
fprintf('Total runtime: %.2f minutes\n', toc/60);
fprintf('End time: %s\n', datetime('now'));
fprintf('========================================\n');
