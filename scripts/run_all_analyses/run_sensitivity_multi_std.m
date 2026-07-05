% run_sensitivity_multi_std.m
% Fast 1D sensitivity analysis with MULTI-TIMESCALE short-term depression.
%
% Dedicated wrapper around the same ParamSpaceAnalysis2 sweep used by
% run_sensitivity_analysis.m, but with two STD timescales:
%   tau_b_E_rec = [1, 3] s  (=> n_b_E = 2 for the STD conditions)
%   tau_b_E_rel =  0.25 s   (default, shared release across timescales)
% The synaptic depression factor is the product prod_m b_{i,m}.
%
% Runs in 'fast' mode (few levels/reps, ode_rk4). Leaves the shared
% orchestrator and run_sensitivity_analysis.m untouched.
%
% See also: run_sensitivity_analysis, ParamSpaceAnalysis2, SRNNModel2

close all; clear; clc;

%% Setup paths
setup_paths();

%% Configuration
run_mode = 'fast';                 % few levels/reps for a quick check
save_figs = true;                 % set true to also write fig/png
figs_stay_open = true;            % true: leave figures open after the run (rapid iteration)
use_stimulus = false;             % false: autonomous network (no external input) => intrinsic LLE
std_zero_floor = true;            % match standalone sensitivity default
tau_b_E_rec = [0.5, 2];              % two STD recovery timescales (s)
n_b_E_multi = numel(tau_b_E_rec);  % => 2

switch run_mode
    case 'fast'
        % Fast iteration: fewer levels/reps, half the sample rate, and a
        % short time range. fs=200 keeps Benettin's lya_dt/dt guard
        % satisfied (0.02/0.005 = 4 >= 3). T_range=[0,20] with an explicit
        % 10 s transient skip => a 10 s LLE window [10,20].
        n_levels = 7;  n_reps = 3;  ode_solver_mode = @ode_rk4;
        fs_mode = 200;  T_range_mode = [0, 20];  n_mode = 200;
        lya_T_interval_mode = [10, 20];   % skip first 10 s (overrides default 15 s)
    case 'production'
        n_levels = 25; n_reps = 50; ode_solver_mode = @ode45;
        fs_mode = 400;  T_range_mode = [0, 50];  n_mode = 300;
        lya_T_interval_mode = [];         % [] => class default (skip first 15 s)
    otherwise, error('run_sensitivity_multi_std:badMode', ...
        'Unknown run_mode ''%s'' (expected ''fast'' or ''production'').', run_mode);
end
fprintf('[run_sensitivity_multi_std] run_mode=%s, n_levels=%d, n_reps=%d, ode_solver=%s, fs=%d, T_range=[%g %g], n=%d, use_stimulus=%d\n', ...
    run_mode, n_levels, n_reps, func2str(ode_solver_mode), fs_mode, T_range_mode(1), T_range_mode(2), n_mode, use_stimulus);
fprintf('[run_sensitivity_multi_std] tau_b_E_rec = [%s], n_b_E = %d\n', ...
    strjoin(compose('%.3g', tau_b_E_rec), ', '), n_b_E_multi);
note = 'sensitivity_multi_std';

% LLE histogram y-axis range for the sensitivity heatmaps (plot_sensitivity).
lle_hist_range = [-2, 2];

% Adaptation conditions: same four regimes as the default sensitivity run,
% but the STD conditions use n_b_E = 2 (two timescales) instead of 1.
% multi_std_conditions = { ...
%     struct('name', 'no_adaptation', 'n_a_E', 0, 'n_b_E', 0), ...
%     struct('name', 'sfa_only',      'n_a_E', 3, 'n_b_E', 0), ...
%     struct('name', 'std_only',      'n_a_E', 0, 'n_b_E', n_b_E_multi), ...
%     struct('name', 'sfa_and_std',   'n_a_E', 3, 'n_b_E', n_b_E_multi) ...
%     };
multi_std_conditions = { ...
    struct('name', 'sfa_and_std',   'n_a_E', 3, 'n_b_E', n_b_E_multi) ...
    };

% Parameters to sweep: {param_name, [min, max]}
params_to_sweep = {
    'level_of_chaos', [0.5, 2];
};

%% Run sensitivity analysis for each parameter
all_output_dirs = {};

for p_idx = 1:size(params_to_sweep, 1)
    param_name = params_to_sweep{p_idx, 1};
    param_range = params_to_sweep{p_idx, 2};

    fprintf('\n========================================\n');
    fprintf('=== Multi-STD sensitivity: %s [%.3g, %.3g] (%d/%d) ===\n', ...
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
    psa.folder_prefix = 'sensitivity_new_defaults';
    psa.model_defaults.ode_solver = ode_solver_mode;   % fast=ode_rk4
    psa.model_defaults.fs = fs_mode;                    % fast=200 (default 400)
    psa.model_defaults.T_range = T_range_mode;          % fast=[0,30] (default [0,50])
    psa.model_defaults.n = n_mode;                      % fast=200 (default 300)
    if ~isempty(lya_T_interval_mode)
        psa.model_defaults.lya_T_interval = lya_T_interval_mode;  % fast: LLE window [10,20]
    end
    if ~use_stimulus
        psa.model_defaults.u_ex_scale = 0;              % autonomous: zero external input (no OFF/ON/OFF)
    end
    psa.model_defaults.std_zero_floor = std_zero_floor;
    psa.model_defaults.tau_b_E_rec = tau_b_E_rec;       % two STD timescales
    psa.model_defaults.c_E = 1/3;
    psa.set_conditions(multi_std_conditions);           % STD conditions use n_b_E=2

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
            sprintf('sensitivity_multi_std_%s', param_name), [], {'fig', 'png'});
        fprintf('Figures saved to %s\n', fig_dir);
    end
    if ~figs_stay_open
        close all;   % skipped when figs_stay_open=true so results stay on screen
    end
end

%% Summary
fprintf('\n=== Multi-STD Sensitivity Analysis Complete ===\n');
fprintf('Parameters analyzed:\n');
for p_idx = 1:size(params_to_sweep, 1)
    fprintf('  %s: %s\n', params_to_sweep{p_idx, 1}, all_output_dirs{p_idx});
end
fprintf('Done!\n');
