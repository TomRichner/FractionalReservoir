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

close all; clearvars -except rng_seeds; clc;

%% Setup paths
setup_paths();

%% seed
if exist('rng_seeds', 'var')
    rng_seeds = rng_seeds + 1
else
    rng_seeds = [1 2]
end

%% Configuration
run_mode = 'fast';                 % few levels/reps for a quick check
save_figs = true;                 % set true to also write fig/png
figs_stay_open = true;            % true: leave figures open after the run (rapid iteration)
use_stimulus = false;             % false: autonomous network (no external input) => intrinsic LLE
show_example_timeseries = true;   % after the sweep, build+plot() one example network
example_level_of_chaos = 10.0;     % level_of_chaos for that example time series
activation_type = 'logistic';     % 'piecewise' (r_max-scaled) or 'logistic' (slope 1 at x=0)
% Firing-rate ceiling for the activation. The nonlinearity is the standard
% piecewise sigmoid stretched on both axes: r(x) = r_max*phi((x-c)/r_max,a,0),
% so the linear-region slope stays 1 but saturation moves to r_max instead of
% 1 (=> higher rates give STD/SFA more dynamic range). r_max=1 reproduces the
% default nonlinearity exactly. S_a/S_c match the SRNNModel2 defaults.
r_max = 1;
S_a_val = 0.9;                    % piecewise linear fraction (SRNNModel2 default)
% Center scaled by r_max so the resting operating point stays at the same
% fraction of the ceiling as the default (0.35/1). Without this scaling a
% larger r_max would push the baseline rate up toward half-max.
S_c_val = 0.35 * r_max;           % activation center (0.35 at r_max=1)
std_zero_floor = true;            % match standalone sensitivity default
tau_b_E_rec = [1];              % two STD recovery timescales (s)
n_b_E_multi = numel(tau_b_E_rec);  % => 2

% Build the activation handles once (empty => use the SRNNModel2 default).
% Reused by both the parameter sweep and the example time series below.
switch activation_type
    case 'piecewise'
        if r_max ~= 1
            % Slope-1 sigmoid saturating at r_max (see r_max comment above).
            act_fn = @(x) r_max * SRNNModel2.piecewiseSigmoid((x - S_c_val) ./ r_max, S_a_val, 0);
            act_fn_deriv = @(x) SRNNModel2.piecewiseSigmoidDerivative((x - S_c_val) ./ r_max, S_a_val, 0);
        else
            act_fn = [];  act_fn_deriv = [];   % r_max=1 => SRNNModel2 default
        end
    case 'logistic'
        % Logistic sigmoid with center S_c_val and unit slope at its center:
        %   sigma(x)  = 1/(1+exp(-4*(x - S_c_val)))
        %   => sigma(S_c_val)=0.5, sigma'(S_c_val)=1, sigma'(x)=4*sigma*(1-sigma)
        % Range (0,1), smooth (never fully saturates), r>=0 (valid rate).
        % S_c_val = 0 => slope 1 at x=0. S_c_val > 0 shifts the inflection to
        % positive x, so at the resting point (x~0) the network sits on the
        % LOWER part of the curve (baseline rate < 0.5, gain < 1).
        act_fn = @(x) 1 ./ (1 + exp(-4 * (x - S_c_val)));
        act_fn_deriv = @(x) 4 .* (1 ./ (1 + exp(-4 * (x - S_c_val)))) .* (1 - 1 ./ (1 + exp(-4 * (x - S_c_val))));
    otherwise
        error('run_sensitivity_multi_std:badActivation', ...
            'Unknown activation_type ''%s'' (expected ''piecewise'' or ''logistic'').', activation_type);
end

switch run_mode
    case 'fast'
        % Fast iteration: fewer levels/reps, half the sample rate, and a
        % short time range. fs=200 keeps Benettin's lya_dt/dt guard
        % satisfied (0.02/0.005 = 4 >= 3). T_range=[0,20] with an explicit
        % 10 s transient skip => a 10 s LLE window [10,20].
        n_levels = 7;  n_reps = 5;  ode_solver_mode = @ode_rk4;
        fs_mode = 200;  T_range_mode = [0, 20];  n_mode = 300;
        lya_T_interval_mode = [10, 20];   % skip first 10 s (overrides default 15 s)
    case 'medium'
        % Medium: roughly halfway between fast and production. Fixed-step
        % ode_rk4 at fs=200, T_range=[0,30] with a 15 s LLE window [15,30].
        n_levels = 15; n_reps = 20; ode_solver_mode = @ode_rk4;
        fs_mode = 200;  T_range_mode = [0, 30];  n_mode = 300;
        lya_T_interval_mode = [15, 30];
    case 'production'
        n_levels = 25; n_reps = 50; ode_solver_mode = @ode45;
        fs_mode = 400;  T_range_mode = [0, 50];  n_mode = 300;
        lya_T_interval_mode = [];         % [] => class default (skip first 15 s)
    otherwise, error('run_sensitivity_multi_std:badMode', ...
        'Unknown run_mode ''%s'' (expected ''fast'', ''medium'', or ''production'').', run_mode);
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
multi_std_conditions = { ...
    struct('name', 'no_adaptation', 'n_a_E', 0, 'n_b_E', 0), ...
    struct('name', 'sfa_only',      'n_a_E', 3, 'n_b_E', 0), ...
    struct('name', 'std_only',      'n_a_E', 0, 'n_b_E', n_b_E_multi), ...
    struct('name', 'sfa_and_std',   'n_a_E', 3, 'n_b_E', n_b_E_multi) ...
    };
% multi_std_conditions = { ...
%     struct('name', 'sfa_and_std',   'n_a_E', 3, 'n_b_E', n_b_E_multi) ...
%     };

% Parameters to sweep: {param_name, [min, max]}
params_to_sweep = {
    'level_of_chaos', [0.5, 3];
};
% params_to_sweep = {
%     'n',              [100, 1000];
%     'f',              [0.25, 0.75];
%     'level_of_chaos', [0.5, 4];
% };

%% Example time series (same model settings, fixed level_of_chaos)
% Build one network with the sweep's settings and call model.plot() so the
% dynamics can be eyeballed before the sweep runs. Nothing is written to disk.
if show_example_timeseries
    fprintf('\n=== Example time series (level_of_chaos = %.2f) ===\n', example_level_of_chaos);
    example_args = { ...
        'n', n_mode, 'indegree', 100, ...
        'rng_seeds', rng_seeds, ...            % fresh network each run (see seed block)
        'level_of_chaos', example_level_of_chaos, ...
        'n_a_E', 3, 'n_b_E', n_b_E_multi, 'c_E', 1/3, ...
        'tau_b_E_rec', tau_b_E_rec, 'std_zero_floor', std_zero_floor, ...
        'fs', fs_mode, 'T_range', T_range_mode, 'ode_solver', ode_solver_mode, ...
        'store_full_state', true};
    if ~isempty(lya_T_interval_mode)
        example_args = [example_args, {'lya_T_interval', lya_T_interval_mode}];
    end
    if ~use_stimulus
        example_args = [example_args, {'u_ex_scale', 0}];
    end
    if ~isempty(act_fn)
        example_args = [example_args, ...
            {'activation_function', act_fn, 'activation_function_derivative', act_fn_deriv}];
    end
    example_model = SRNNModel2(example_args{:});
    example_model.build();
    example_model.run();
    fprintf('[example] level_of_chaos=%.2f, LLE=%+.4f\n', ...
        example_level_of_chaos, example_model.lya_results.LLE);
    example_model.plot();
end

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
    if ~isempty(act_fn)
        psa.model_defaults.activation_function = act_fn;
        psa.model_defaults.activation_function_derivative = act_fn_deriv;
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
