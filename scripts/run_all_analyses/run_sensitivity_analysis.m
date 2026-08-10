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
switch run_mode
    case 'fast'
        % Fast iteration: fewer levels/reps, half the sample rate, short time
        % range. ode_rk4 at fs=200, T_range=[0,10] with an explicit 5 s LLE
        % window [5,10].
        n_levels = 4;  n_reps = 3;  ode_solver_mode = @ode_rk4;
        fs_mode = 200;  T_range_mode = [0, 10];  lya_T_interval_mode = [5, 10];
    case 'medium'
        % Medium: roughly halfway between fast and production. ode_rk4 at
        % fs=200, T_range=[0,20] with a 10 s LLE window [10,20], 11 levels x
        % 15 reps. fs=200 gives dt=0.005, so benettin's lya_dt=0.1 clears the
        % lya_dt/dt >= 3 guard with room to spare (ratio 20).
        n_levels = 11; n_reps = 15; ode_solver_mode = @ode_rk4;
        fs_mode = 200;  T_range_mode = [0, 20];  lya_T_interval_mode = [10, 20];
    case 'production'
        n_levels = 25; n_reps = 50; ode_solver_mode = @ode45;
        fs_mode = 400;  T_range_mode = [0, 50];  lya_T_interval_mode = [20 50];
    otherwise, error('run_sensitivity_analysis:badMode', ...
            'Unknown run_mode ''%s'' (expected ''fast'', ''medium'', or ''production'').', run_mode);
end
fprintf('[run_sensitivity_analysis] run_mode=%s, n_levels=%d, n_reps=%d, ode_solver=%s, fs=%d, T_range=[%g %g]\n', ...
    run_mode, n_levels, n_reps, func2str(ode_solver_mode), fs_mode, T_range_mode(1), T_range_mode(2));
note = 'sensitivity';

% LLE histogram y-axis range for the sensitivity heatmaps (plot_sensitivity).
lle_hist_range = [-2, 2];

% Parameters to sweep: {param_name, [min, max]}
params_to_sweep = {
    'n',              [100, 1000];
    'f',              [0.2, 0.8];
    'level_of_chaos', [0.25, 4];
    };

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
    if exist('master_output_dir', 'var')
        psa.output_dir = master_output_dir;
    end
    psa.model_defaults.ode_solver = ode_solver_mode;  % fast=ode_rk4, production=ode45
    psa.model_defaults.fs = fs_mode;                   % fast=200 (default 400)
    psa.model_defaults.T_range = T_range_mode;         % fast=[0,20] (default [0,50])
    if ~isempty(lya_T_interval_mode)
        psa.model_defaults.lya_T_interval = lya_T_interval_mode;  % fast: LLE window [10,20]
    end
    
    % --- STD with one recovery timescale ----------------------------------
    % n_b_E must come from the condition (ParamSpaceAnalysis2 ignores it in
    % model_defaults). tau_b_E_rec / tau_b_E_rel are left at the class defaults
    % (1 s / 0.25 s), matching the pre-multi-timescale runs. STD is on E only.
    psa.set_conditions({ ...
        struct('name', 'no_adaptation', 'n_a_E', 0, 'n_b_E', 0), ...
        struct('name', 'sfa_only',      'n_a_E', 3, 'n_b_E', 0), ...
        struct('name', 'std_only',      'n_a_E', 0, 'n_b_E', 1), ...
        struct('name', 'sfa_and_std',   'n_a_E', 3, 'n_b_E', 1) ...
        });
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
