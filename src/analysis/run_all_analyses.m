function run_dir = run_all_analyses(preset_name, run_mode, opts)
% RUN_ALL_ANALYSES The sweep pipeline: sensitivity, tau sensitivity, param space.
%
%   run_dir = RUN_ALL_ANALYSES()
%   run_dir = RUN_ALL_ANALYSES(preset_name, run_mode)
%   run_dir = RUN_ALL_ANALYSES(preset_name, run_mode, 'save_figs', false)
%
% Runs the three ParamSpaceAnalysis2 sweeps into a single dated directory under
% data/param_space/, with a git-provenance record, a machine-readable
% run_manifest.mat and a human-readable parameters.md, so the directory
% describes itself.
%
% This is the SWEEP pipeline only. Memory capacity and the eig-heatmap sampling
% are figure-specific compute and live one layer up, in
% scripts/paper/run_all_paper_analyses.m, which calls this.
%
% INPUTS
%   preset_name  which network to simulate. See srnn_param_preset.
%                Carries the model class and the adaptation conditions too.
%   run_mode     how much compute to spend: 'fast' | 'fast2' | 'medium' |
%                'medium2' | 'production'. See analysis_run_config.m for the
%                per-analysis numbers. 'fast' is the smoke test.
%   save_figs    write extra figure formats alongside each sweep (default true).
%   output_dir   override the dated directory (default: create one).
%   assemble     stack the 1-D sensitivity sheets after the sweeps (default true).
%
% WAS A SCRIPT that ran three other scripts, passing settings through
% base-workspace `master_*` variables the callees read with exist(). Everything
% is arguments and return values now; see resolve_run_context for why.
%
% TWO DEFAULTS CHANGED with that conversion, both deliberately:
%   * preset_name defaults to the paper's operating point rather than
%     'celltype_pairs_S_c_by_type', which is three generations upstream of it
%     and was left behind when the target moved.
%   * save_figs is a logical. The old master_save_figs had a third value,
%     'follow_scripts_save_figs', which meant "let each sub-script use its own
%     local save_figs flag". There are no local flags any more -- that was the
%     base-workspace protocol -- so the option had nothing left to mean.
%
% See also: resolve_run_context, run_sensitivity_analysis,
%           run_tau_sensitivity_analysis, run_param_space_analysis,
%           srnn_param_preset, analysis_run_config, write_run_parameters_md

arguments
    preset_name       (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    run_mode          (1,:) char    = 'medium'
    opts.save_figs    (1,1) logical = true
    opts.output_dir   (1,:) char    = ''
    opts.assemble     (1,1) logical = true
    opts.verbose      (1,1) logical = true
end

setup_paths();

fprintf('=== Starting All Analyses ===\n');
fprintf('Start time: %s\n\n', datetime('now'));
t_start = tic;

%% Output directory
% project_root comes from setup_paths.m's own location (the repo root), not from
% this file's, so this function tolerates living in any subdirectory.
project_root = fileparts(which('setup_paths'));
if isempty(opts.output_dir)
    dt_str  = lower(datestr(now, 'mmm_dd_yy_HH_MM')); %#ok<TNOW1,DATST>
    run_dir = fullfile(project_root, 'data', 'param_space', ...
        sprintf('run_all_%s', dt_str));
else
    run_dir = opts.output_dir;
    [~, folder] = fileparts(run_dir);
    dt_str = regexprep(folder, '^run_all_', '');
end
if ~exist(run_dir, 'dir')
    mkdir(run_dir);
end
fprintf('Master output directory: %s\n\n', run_dir);

prov = capture_git_provenance(run_dir, project_root);

%% Resolve the preset once, for the manifest
% Each sub-analysis resolves its own ctx from the same preset NAME, which is
% deterministic -- srnn_param_preset is a pure lookup. Resolving here as well is
% only so the manifest can record the class and the conditions before any
% compute starts.
[preset_defaults, model_class, conditions] = srnn_param_preset(preset_name);
fprintf('Parameter preset: %s (%d override(s)) on %s\n', ...
    preset_name, numel(fieldnames(preset_defaults)), model_class);
fprintf('Conditions: %s\n', ...
    strjoin(cellfun(@(c) c.name, conditions, 'UniformOutput', false), ', '));
fprintf('Run mode: %s\n\n', run_mode);

save_manifest(run_dir, preset_name, preset_defaults, model_class, run_mode, ...
    opts.save_figs, prov, dt_str);

%% Shared context arguments for every sub-analysis
ctx_args = {'preset_name', preset_name, 'run_mode', run_mode, ...
            'output_dir', run_dir, 'save_figs', opts.save_figs, ...
            'verbose', opts.verbose};

%% 1. Sensitivity
fprintf('========================================\n');
fprintf('[1/3] Running Sensitivity Analysis...\n');
fprintf('========================================\n');
% A fresh pool before each stage. The three analyses run back to back for hours
% against one pool otherwise; see restart_parpool for why that is worth
% avoiding and for what is and is not established about the aug_13 failure.
restart_parpool();
run_sensitivity_analysis(resolve_run_context('sensitivity', ctx_args{:}));

%% 1b. Assemble the 1-D sensitivity sheets
% run_sensitivity_analysis saves each swept parameter's figure into its own
% 1D_sensitivity_* subfolder. replot_sensitivity gathers them into a single
% replot_sensitivity_<dt>/figures/ folder, then assemble_sensitivity_figure
% stacks them.
%
% TWO sheets, split by what the parameters MEAN rather than one sheet of all
% seven. A row is only worth reading against the rows it is comparable to: the
% network-scale parameters answer "how does this network respond to being made
% bigger, more excitatory, more chaotic", while the connectivity blocks answer
% "which route does the stability actually hinge on". Interleaving them
% alphabetically -- f_E, level_of_chaos, mu_EE, mu_EI, mu_IE, mu_II, n -- buries
% both questions.
if opts.assemble && opts.save_figs
    fprintf('========================================\n');
    fprintf('Assembling 1D sensitivity figures...\n');
    fprintf('========================================\n');
    sens_replot_dir = replot_sensitivity(run_dir);

    if strcmp(model_class, 'SRNNCellTypePairs')
        f_param = 'f_E';
    else
        f_param = 'f';
    end
    assemble_sensitivity_figure(sens_replot_dir, 'LLE', ...
        {f_param, 'n', 'level_of_chaos'}, 'network');

    % The four connectivity blocks, in (post, pre) reading order rather than
    % alphabetical. What matters is whether E->E and I->I behave alike and
    % whether the two cross terms do, which only reads side by side. Skipped
    % with a warning on a class that has no mu blocks, so this costs nothing on
    % an SRNNModel2 run.
    assemble_sensitivity_figure(sens_replot_dir, 'LLE', ...
        {'mu_EE_relative', 'mu_EI_relative', 'mu_IE_relative', 'mu_II_relative'}, ...
        'mu_blocks');
end

%% 2. Tau sensitivity
fprintf('========================================\n');
fprintf('[2/3] Running Tau Sensitivity Analysis...\n');
fprintf('========================================\n');
restart_parpool();
run_tau_sensitivity_analysis(resolve_run_context('tau_sensitivity', ctx_args{:}));

%% 3. Parameter space
fprintf('========================================\n');
fprintf('[3/3] Running Parameter Space Analysis...\n');
fprintf('========================================\n');
restart_parpool();
run_param_space_analysis(resolve_run_context('param_space', ctx_args{:}));

%% Human-readable parameter record
% parameters.md next to run_manifest.mat: the preset and run mode, what the
% preset set, every parameter in effect (as run under sfa_and_std), the four
% conditions, and each analysis's sweep axes and timings. Built from the saved
% artifacts alone, so it can be regenerated for any past run directory with
% write_run_parameters_md(<dir>).
%
% Wrapped because a failure to write a description must not be the last thing a
% completed overnight run says.
try
    params_md = write_run_parameters_md(run_dir);
    fprintf('Parameter record: %s\n', params_md);
catch ME
    warning('run_all_analyses:ParametersMdFailed', ...
        'Could not write parameters.md: %s', ME.message);
end

%% Summary
fprintf('========================================\n');
fprintf('=== All Analyses Complete ===\n');
fprintf('Total runtime: %.2f minutes\n', toc(t_start)/60);
fprintf('End time: %s\n', datetime('now'));
fprintf('Output: %s\n', run_dir);
fprintf('========================================\n');
end

%% ------------------------------------------------------------------------
function save_manifest(run_dir, preset_name, preset_defaults, model_class, ...
    run_mode, save_figs, prov, dt_str)
% Machine-readable run manifest, for reproducibility and self-documentation.
%
% The activation fields record the nonlinearity in force: the preset's, where it
% overrides one, otherwise the class default.
%
% NOTE: these fields are largely redundant. Each PSA freezes its FULL effective
% parameter set into psa.resolved_defaults (saved in psa_object.mat and
% param_space_summary.mat), and since the nonlinearity is named data
% (`activation` + S_a/S_c) rather than a function handle, same_config compares it
% EXACTLY. They are kept for continuity with runs made before resolved_defaults
% existed -- and because preset_default_values and write_run_parameters_md read
% preset_name and model_class from here.
activation_fields = struct('activation', [], 'S_a', [], 'S_c', []);
for fn = fieldnames(activation_fields)'
    if isfield(preset_defaults, fn{1})
        activation_fields.(fn{1}) = preset_defaults.(fn{1});
    else
        % class_default errors for a class with required constructor arguments
        % (SRNNCellTypePairs). Nothing here is worth failing a run over, and in
        % practice such a preset has to name the nonlinearity itself anyway.
        try
            activation_fields.(fn{1}) = ...
                ParamSpaceAnalysis2.class_default(fn{1}, model_class);
        catch
            activation_fields.(fn{1}) = 'unknown';
        end
    end
end

if save_figs
    save_figs_str = 'save_all_figs';
else
    save_figs_str = 'save_no_figs';
end

run_manifest = struct( ...
    'run_mode', run_mode, ...
    'preset_name', preset_name, ...
    'model_class', model_class, ...
    'master_save_figs', save_figs_str, ...
    'activation', activation_fields.activation, ...
    'S_a', activation_fields.S_a, ...
    'S_c', activation_fields.S_c, ...
    'git_commit', prov.commit, ...
    'git_commit_short', prov.commit_short, ...
    'git_branch', prov.branch, ...
    'git_dirty', prov.is_dirty, ...
    'matlab_version', version, ...
    'timestamp', dt_str);
save(fullfile(run_dir, 'run_manifest.mat'), 'run_manifest');
end
