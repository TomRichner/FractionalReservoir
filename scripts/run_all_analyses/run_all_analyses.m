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
% Derive project_root from setup_paths.m (which lives at the repo root) rather than
% from this file's own location, so run_all_analyses.m can live in any
% subdirectory of scripts/ without breaking the output paths.
project_root = fileparts(which('setup_paths'));
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
%   'fast'       - few levels/reps, fs=200, short T_range; finishes quickly (testing)
%   'fast2'      - 'fast' but with twice the reps on the 1-D sensitivity sweeps,
%                  so their histograms have enough samples per level to read
%   'medium'     - fs=400, 11 levels, T_range=[0,20] (tau: ode45 over [0,50])
%   'production' - full-size sweeps, fs=400, T_range=[0,50] (for real results)
% See analysis_run_config.m for the actual per-analysis numbers.
% Defaults to 'medium'. To pick another WITHOUT editing this file, set the
% variable first in the console:   run_mode = 'fast'; run_all_analyses
if ~exist('run_mode', 'var'); run_mode = 'medium'; end
fprintf('Run mode: %s\n\n', run_mode);

% Which network to simulate, as opposed to how much compute to spend on it
% (that is run_mode). Edit this ONE line to run the whole pipeline on a
% different parameter set; see src/srnn_param_preset.m for the available names
% and for what may and may not go in a preset.
%
% The preset also names the MODEL CLASS its overrides are written for -- the two
% classes have disjoint parameter vocabularies -- and the ADAPTATION CONDITIONS
% to sweep them under, since a preset with its own depression timescales can only
% express them through a condition. The sub-scripts pick their sweep axes and
% plot settings off master_model_class.
if ~exist('preset_name', 'var'); preset_name = 'celltype_pairs_S_c_by_type'; end
[master_model_overrides, master_model_class, master_conditions] = ...
    srnn_param_preset(preset_name);
fprintf('Parameter preset: %s (%d override(s)) on %s\n', ...
    preset_name, numel(fieldnames(master_model_overrides)), master_model_class);
fprintf('Conditions: %s\n\n', ...
    strjoin(cellfun(@(c) c.name, master_conditions, 'UniformOutput', false), ', '));

% Machine-readable run manifest for reproducibility + self-documentation. The
% activation fields record the nonlinearity in force: the preset's, where it
% overrides one, otherwise the class default.
%
% This used to default-construct an SRNNModel2 and func2str its handles.
% Neither survives the move to multiple model classes: SRNNCellTypePairs has
% five required constructor arguments and so cannot be default-constructed at
% all, and the nonlinearity is named data now, so the NAME is the honest record.
%
% NOTE: these fields are largely redundant. Each PSA freezes its FULL effective
% parameter set into psa.resolved_defaults (saved in psa_object.mat and
% param_space_summary.mat), and since the nonlinearity is named data
% (`activation` + S_a/S_c) rather than a function handle, same_config compares
% it EXACTLY. The manifest fields are kept only for continuity with runs made
% before resolved_defaults existed.
activation_fields = struct('activation', [], 'S_a', [], 'S_c', []);
for fn = fieldnames(activation_fields)'
    if isfield(master_model_overrides, fn{1})
        activation_fields.(fn{1}) = master_model_overrides.(fn{1});
    else
        % class_default errors for a class with required constructor arguments
        % (SRNNCellTypePairs). Nothing here is worth failing a run over, and in
        % practice such a preset has to name the nonlinearity itself anyway.
        try
            activation_fields.(fn{1}) = ...
                ParamSpaceAnalysis2.class_default(fn{1}, master_model_class);
        catch
            activation_fields.(fn{1}) = 'unknown';
        end
    end
end
run_manifest = struct( ...
    'run_mode', run_mode, ...
    'preset_name', preset_name, ...
    'model_class', master_model_class, ...
    'master_save_figs', master_save_figs, ...
    'activation', activation_fields.activation, ...
    'S_a', activation_fields.S_a, ...
    'S_c', activation_fields.S_c, ...
    'git_commit', prov.commit, ...
    'git_commit_short', prov.commit_short, ...
    'git_branch', prov.branch, ...
    'git_dirty', prov.is_dirty, ...
    'matlab_version', version, ...
    'timestamp', dt_str);
save(fullfile(master_output_dir, 'run_manifest.mat'), 'run_manifest');
clear activation_fields fn;

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

    % A second sheet with just the four connectivity blocks, in (post, pre)
    % reading order rather than the alphabetical order the call above uses.
    % They are only interpretable against each other -- what matters is whether
    % E->E and I->I behave alike, and whether the two cross terms do -- and on
    % the all-parameters sheet they are interleaved with n, f_E and
    % level_of_chaos. Skipped with a warning on a model class that has no mu
    % blocks, so this costs nothing on an SRNNModel2 run.
    assemble_sensitivity_figure(sens_replot_dir, 'LLE', ...
        {'mu_EE_relative', 'mu_EI_relative', 'mu_IE_relative', 'mu_II_relative'}, ...
        'mu_blocks');
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

% %% 4. DC Lyapunov Analysis
% fprintf('========================================\n');
% fprintf('[4/4] Running DC Lyapunov Analysis...\n');
% fprintf('========================================\n');
% run_dc_lle_analysis;

% %% 4b. Replot the DC Lyapunov figure from saved data
% % run_dc_lle_analysis saves dc_lle_results.mat into a dc_lle_nSeeds_* subfolder.
% % replot_dc_lle rebuilds the DC-vs-LLE confplot figure from that .mat alone.
% if ~strcmp(master_save_figs, 'save_no_figs')
%     fprintf('========================================\n');
%     fprintf('Replotting DC Lyapunov figure...\n');
%     fprintf('========================================\n');
%     replot_dc_lle(master_output_dir);
% end

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
