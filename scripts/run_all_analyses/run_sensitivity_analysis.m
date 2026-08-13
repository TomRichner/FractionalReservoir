% run_sensitivity_analysis.m
% Sensitivity analysis using ParamSpaceAnalysis2 in 1D mode
%
% Sweeps individual parameters across multiple levels and repetitions,
% comparing different adaptation conditions. Each parameter is swept
% independently (one PSA run per parameter).
%
% See also: ParamSpaceAnalysis2, SRNNModel2

% No clear/clc/close all on standalone runs, so base-workspace settings such as
% run_mode and save_figs (set in the console before running) survive into this
% script. run_all_analyses relies on the same: it never clears the sub-scripts.

%% Figure saving configuration
if exist('master_save_figs', 'var')
    if strcmp(master_save_figs, 'save_all_figs')
        save_figs = true;
    elseif strcmp(master_save_figs, 'save_no_figs')
        save_figs = false;
    end
end
if ~exist('save_figs', 'var')
    save_figs = false;
end

%% Setup paths
setup_paths();

%% Analysis Configuration
% Run mode: 'fast' for quick checks, 'production' for full-size runs.
% Set run_mode in the base workspace (or via run_all_analyses); defaults to
% 'production' when this script is run standalone.
if ~exist('run_mode', 'var'); run_mode = 'production'; end
cfg = analysis_run_config('sensitivity', run_mode);
n_levels = cfg.n_levels;
n_reps = cfg.n_reps;
fprintf('[run_sensitivity_analysis] run_mode=%s, n_levels=%d, n_reps=%d, ode_solver=%s, fs=%d, T_range=[%g %g]\n', ...
    run_mode, n_levels, n_reps, func2str(cfg.model.ode_solver), cfg.model.fs, ...
    cfg.model.T_range(1), cfg.model.T_range(2));

% Which network to simulate, and with which model class. run_all_analyses sets
% both from srnn_param_preset; standalone runs get SRNNModel2 class defaults.
if exist('master_model_overrides', 'var')
    preset_defaults = master_model_overrides;
else
    preset_defaults = srnn_param_preset('default');
end
if exist('master_model_class', 'var')
    model_class = master_model_class;
else
    model_class = 'SRNNModel2';
end
% Conditions come from the preset when there is one, because a preset may carry
% its own depression routes, which reach the model only through a condition.
if exist('master_conditions', 'var')
    conditions = master_conditions;
else
    conditions = srnn_adaptation_conditions(model_class);
end

note = 'sensitivity';

% LLE histogram y-axis range for the sensitivity heatmaps (plot_sensitivity).
lle_hist_range = [-2, 2];

% Parameters to sweep: {param_name, [min, max]}
% The fraction-excitatory axis is named differently per class: SRNNModel2's f is
% the scalar fraction, whereas SRNNCellTypePairs' f is a 1 x C row, so the sweep
% goes through its scalar alias f_E (which sets f = [f_E, 1-f_E]).
params_to_sweep = {
    'n',              [100, 1000];
    'f',              [0.2, 0.8];
    'level_of_chaos', [0.25, 4];
    };
if strcmp(model_class, 'SRNNCellTypePairs')
    params_to_sweep{2, 1} = 'f_E';
end

% params_to_sweep = {
%     'level_of_chaos', [0.5, 2];
% };

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
    psa.model_class = model_class;
    % The class default carries SRNNModel2's adaptation counts; only n and
    % indegree are integer-valued axes common to both classes.
    psa.integer_params = {'n', 'indegree'};
    if exist('master_output_dir', 'var')
        psa.output_dir = master_output_dir;
    end
    % Parameter preset first, run_mode timings second, so run_mode keeps final
    % say over ode_solver / fs / T_range / lya_T_interval.
    psa.model_defaults = merge_struct(preset_defaults, cfg.model);

    % --- Adaptation conditions --------------------------------------------
    psa.set_conditions(conditions);

    if strcmp(model_class, 'SRNNModel2')
        % STD with two recovery timescales. Give the STD conditions n_b_E = 2
        % (two E depression timescales, product of the two b columns) with a
        % 1x2 recovery-time-constant vector; the release constant tau_b_E_rel
        % stays scalar/shared. n_b_E must come from the condition
        % (ParamSpaceAnalysis2 ignores it in model_defaults), while tau_b_E_rec
        % flows through model_defaults. STD is on E neurons only.
        %
        % SRNNCellTypePairs says the same thing per route inside
        % synapse_config, so there is nothing to override here for that class --
        % and tau_b_E_rec is not one of its properties, so setting it would be a
        % hard validate_model_defaults error rather than a no-op.
        for ci = 1:numel(psa.conditions)
            if psa.conditions{ci}.n_b_E > 0
                psa.conditions{ci}.n_b_E = 2;
            end
        end
        psa.model_defaults.tau_b_E_rec = [0.5, 5];   % two E recovery timescales (s)
    end
    % ----------------------------------------------------------------------
    
    % Add the swept parameter and reps
    psa.add_grid_parameter(param_name, param_range);
    psa.add_grid_parameter('reps', 1:n_reps);
    
    % Run
    psa.run();
    
    % Copy this script for reproducibility
    copyfile([mfilename('fullpath') '.m'], psa.output_dir);
    
    % Plot sensitivity heatmaps
    psa.plot_sensitivity('metric', 'LLE', 'hist_range', lle_hist_range);
    psa.plot_sensitivity('metric', 'mean_rate');
    
    % Save PSA object
    save(fullfile(psa.output_dir, 'psa_object.mat'), 'psa');
    
    all_output_dirs{end+1} = psa.output_dir; %#ok<SAGROW>
    
    % Save figures in additional formats
    if save_figs
        fig_dir = fullfile(psa.output_dir, 'figures');
        save_some_figs_to_folder_2(fig_dir, ...
            sprintf('sensitivity_%s', param_name), [], {'fig', 'png'});
        fprintf('Figures saved to %s\n', fig_dir);
    end
    close all;
end

%% Summary
fprintf('\n=== Sensitivity Analysis Complete ===\n');
fprintf('Parameters analyzed:\n');
for p_idx = 1:size(params_to_sweep, 1)
    fprintf('  %s: %s\n', params_to_sweep{p_idx, 1}, all_output_dirs{p_idx});
end
fprintf('Done!\n');
