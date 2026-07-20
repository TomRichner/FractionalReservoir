%% run_all_analyses.m
% Master script to run all analysis pipelines in sequence
%
% This script executes the following analyses in order:
%   1. General sensitivity analysis (SFA/STD parameter exploration)
%   2. Tau sensitivity analysis (tau_a and tau_b parameter sweeps)
%   3. Parameter space analysis
%   4. DC Lyapunov analysis (local LLE vs DC stim level, across seeds)
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
prov = capture_git_provenance(master_output_dir, project_root);

% Control figure saving across all sub-scripts:
%   'save_all_figs'            - Override all scripts to save figures
%   'save_no_figs'             - Override all scripts to NOT save figures
%   'follow_scripts_save_figs' - Let each script use its own save_figs setting
master_save_figs = 'save_all_figs';

% Run mode for the sub-analyses (controls n_levels / n_reps, fs, T_range, and
% the LLE window):
%   'fast'       - few levels/reps, fs=200, T_range=[0,20]; finishes quickly (testing)
%   'medium'     - ode45, fs=400, 15 levels x 50 reps, T_range=[0,30] (tau: [0,50]); ~halfway
%   'production' - full-size sweeps, fs=400, T_range=[0,50] (for real results)
% Defaults to 'medium'. To pick another WITHOUT editing this file, set the
% variable first in the console:   run_mode = 'fast'; run_all_analyses
if ~exist('run_mode', 'var'); run_mode = 'medium'; end
fprintf('Run mode: %s\n\n', run_mode);

% Machine-readable run manifest for reproducibility + self-documentation, and
% so combine_runs can verify pooled runs used the SAME nonlinearity. We
% fingerprint the CLASS-DEFAULT activation via a throwaway SRNNModel2
% (constructor sets it in set_defaults; no build needed) -- the run_all
% sub-scripts all use that default and never override model parameters in
% model_defaults, so it is not otherwise recorded.
m0 = SRNNModel2();
run_manifest = struct( ...
    'run_mode', run_mode, ...
    'master_save_figs', master_save_figs, ...
    'activation', func2str(m0.activation_function), ...
    'activation_derivative', func2str(m0.activation_function_derivative), ...
    'S_a', m0.S_a, ...
    'S_c', m0.S_c, ...
    'git_commit', prov.commit, ...
    'git_commit_short', prov.commit_short, ...
    'git_branch', prov.branch, ...
    'git_dirty', prov.is_dirty, ...
    'matlab_version', version, ...
    'timestamp', dt_str);
save(fullfile(master_output_dir, 'run_manifest.mat'), 'run_manifest');
clear m0;

%% 1. General Sensitivity Analysis
fprintf('========================================\n');
fprintf('[1/4] Running Sensitivity Analysis...\n');
fprintf('========================================\n');
run_sensitivity_analysis;

%% 1b. Assemble 1D sensitivity figures into one stacked figure
% run_sensitivity_analysis saves each swept parameter's LLE figure into its own
% 1D_sensitivity_* subfolder. replot_sensitivity gathers them into a single
% replot_sensitivity_<dt>/figures/ folder, then assemble_sensitivity_figure
% stacks the per-parameter LLE figures into sensitivity_LLE_combined.{fig,png}.
if ~strcmp(master_save_figs, 'save_no_figs')
    fprintf('========================================\n');
    fprintf('Assembling 1D sensitivity figures...\n');
    fprintf('========================================\n');
    sens_replot_dir = replot_sensitivity(master_output_dir);
    assemble_sensitivity_figure(sens_replot_dir, 'LLE');
end

%% 2. Tau Sensitivity Analysis
fprintf('========================================\n');
fprintf('[2/4] Running Tau Sensitivity Analysis...\n');
fprintf('========================================\n');
run_tau_sensitivity_analysis;

%% 3. Parameter Space Analysis
fprintf('========================================\n');
fprintf('[3/4] Running Parameter Space Analysis...\n');
fprintf('========================================\n');
run_param_space_analysis2;

%% 4. DC Lyapunov Analysis
fprintf('========================================\n');
fprintf('[4/4] Running DC Lyapunov Analysis...\n');
fprintf('========================================\n');
run_dc_lle_analysis;

%% 4b. Replot the DC Lyapunov figure from saved data
% run_dc_lle_analysis saves dc_lle_results.mat into a dc_lle_nSeeds_* subfolder.
% replot_dc_lle rebuilds the DC-vs-LLE confplot figure from that .mat alone.
if ~strcmp(master_save_figs, 'save_no_figs')
    fprintf('========================================\n');
    fprintf('Replotting DC Lyapunov figure...\n');
    fprintf('========================================\n');
    replot_dc_lle(master_output_dir);
end

% %% Fig 2: Fraction Excitatory Analysis
% fprintf('========================================\n');
% fprintf('Running Fig 2 Fraction Excitatory Analysis...\n');
% fprintf('========================================\n');
% Fig_2_fraction_excitatory_analysis;

%% Summary
fprintf('========================================\n');
fprintf('=== All Analyses Complete ===\n');
fprintf('Total runtime: %.2f minutes\n', toc/60);
fprintf('End time: %s\n', datetime('now'));
fprintf('========================================\n');
