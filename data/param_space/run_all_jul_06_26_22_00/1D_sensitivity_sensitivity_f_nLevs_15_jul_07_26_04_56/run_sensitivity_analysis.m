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

%% STD zero-floor configuration
% Honor the master override from run_all_analyses; default off when standalone.
if exist('master_std_zero_floor', 'var')
    std_zero_floor = master_std_zero_floor;
end
if ~exist('std_zero_floor', 'var')
    std_zero_floor = false;
end

%% SFA strength (c_E) configuration
% Honor the master override from run_all_analyses; default to the class default.
if exist('master_c_E', 'var')
    c_E = master_c_E;
end
if ~exist('c_E', 'var')
    c_E = 0.15/3;   % SRNNModel2 default
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
        % range. fs=200 keeps Benettin's lya_dt/dt guard satisfied (4>=3);
        % T_range=[0,20] with an explicit 10 s LLE window [10,20].
        n_levels = 5;  n_reps = 4;  ode_solver_mode = @ode_rk4;
        fs_mode = 200;  T_range_mode = [0, 20];  lya_T_interval_mode = [10, 20];
    case 'medium'
        % Medium: roughly halfway between fast and production. ode45 at fs=400,
        % T_range=[0,30] with a 15 s LLE window [15,30], 15 levels x 50 reps.
        n_levels = 15; n_reps = 50; ode_solver_mode = @ode45;
        fs_mode = 400;  T_range_mode = [0, 30];  lya_T_interval_mode = [15, 30];
    case 'production'
        n_levels = 25; n_reps = 50; ode_solver_mode = @ode45;
        fs_mode = 400;  T_range_mode = [0, 50];  lya_T_interval_mode = [];
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
    'f',              [0.25, 0.75];
    'level_of_chaos', [0.5, 3];
};

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
    if exist('master_seed_offset', 'var')
        psa.network_seed_offset = master_seed_offset;  % fresh networks per run (pooling)
    end
    psa.model_defaults.ode_solver = ode_solver_mode;  % fast=ode_rk4, production=ode45
    psa.model_defaults.std_zero_floor = std_zero_floor;
    psa.model_defaults.c_E = c_E;                      % SFA strength (run_all: 0.5/3)
    psa.model_defaults.fs = fs_mode;                   % fast=200 (default 400)
    psa.model_defaults.T_range = T_range_mode;         % fast=[0,20] (default [0,50])
    if ~isempty(lya_T_interval_mode)
        psa.model_defaults.lya_T_interval = lya_T_interval_mode;  % fast: LLE window [10,20]
    end

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
