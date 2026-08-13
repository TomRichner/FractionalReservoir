% run_param_space_analysis2.m
% Example script demonstrating how to use the ParamSpaceAnalysis2 class
%
% This script performs a multi-dimensional parameter space exploration
% comparing all four adaptation conditions on networks with the same
% connectivity matrix W.
%
% Key differences from SensitivityAnalysis:
%   - Multi-dimensional grid (all parameter combinations)
%   - Same network W used across all 4 conditions for fair comparison
%   - Randomized execution order for representative early-stopping
%   - Batched execution with resume capability
%
% See also: ParamSpaceAnalysis2, SRNNModel2, SensitivityAnalysis

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

%% Create ParamSpaceAnalysis2 object
% Configure the analysis parameters
% NOTE: Total simulations = n_levels^(n_params) * n_conditions
%       e.g., 4 levels x 3 params x 4 conditions = 256 simulations

% Run mode: 'fast' for quick checks, 'production' for full-size runs.
% Set run_mode in the base workspace (or via run_all_analyses); defaults to
% 'production' when this script is run standalone.
if ~exist('run_mode', 'var'); run_mode = 'production'; end
cfg = analysis_run_config('param_space', run_mode);
n_levels = cfg.n_levels;
fprintf('[run_param_space_analysis2] run_mode=%s, n_levels=%d, ode_solver=%s, fs=%d, T_range=[%g %g]\n', ...
    run_mode, n_levels, func2str(cfg.model.ode_solver), cfg.model.fs, ...
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

psa = ParamSpaceAnalysis2(...
    'n_levels', n_levels, ...   % set by run_mode (fast=3, production=5)
    'batch_size', 25, ...       % configs per batch (for checkpointing)
    'note', 'test_refactor', ...         % Optional note for folder naming
    'verbose', true ...         % Print progress during execution
    );
psa.model_class = model_class;
psa.integer_params = {'n', 'indegree'};
if exist('master_output_dir', 'var')
    psa.output_dir = master_output_dir;
end
% Parameter preset first, run_mode timings second, so run_mode keeps final say
% over ode_solver / fs / T_range / lya_T_interval.
psa.model_defaults = merge_struct(preset_defaults, cfg.model);

%% Add parameters to the grid
% All combinations of these parameters will be tested
% The order in which parameters are added doesn't matter

% Network structure parameters
% Same swept variables and ranges as run_sensitivity_analysis.m (n, f,
% level_of_chaos), but here as a joint multi-dimensional grid.
% The fraction-excitatory axis is f on SRNNModel2 (a scalar property) and f_E on
% SRNNCellTypePairs (a scalar alias onto the 1 x C row f).
if strcmp(model_class, 'SRNNCellTypePairs')
    f_param = 'f_E';
else
    f_param = 'f';
end
psa.add_grid_parameter('n',              [100, 1000]);   % network size
psa.add_grid_parameter(f_param,          [0.25, 0.75]);  % fraction excitatory
psa.add_grid_parameter('level_of_chaos', [0.5, 3]);      % W scaling (edge of chaos)

% Dynamics parameters (uncomment to include)
% psa.add_grid_parameter('tau_d', [0.05, 0.2]);           % Dendritic time constant
% psa.add_grid_parameter('u_ex_scale', [1.0, 3.0]);       % External input scaling

% Adaptation parameters (uncomment to include)
% psa.add_grid_parameter('c_E', [0.01, 0.5]);             % SFA strength

%% Configure model defaults (optional)
% Set any SRNNModel2 properties that should be constant across all runs
% psa.model_defaults.n = 300;                   % Number of neurons
% psa.model_defaults.T_range = [-10, 20];       % With settling time, similar to example
% psa.model_defaults.fs = 200;                  % Sampling frequency
% psa.model_defaults.c_E = 0.15/3;              % SFA strength (≈0.05), matches example
% psa.model_defaults.tau_b_E_rec = 1;           % STD recovery time for E neurons
% psa.model_defaults.tau_b_I_rec = 1;           % STD recovery time for I neurons
% psa.model_defaults.S_c = 0.4;                 % Activation function center
% psa.model_defaults.u_ex_scale = 1.5;          % External input scaling
% psa.model_defaults.lya_method = 'benettin';   % Lyapunov computation method
% psa.model_defaults.level_of_chaos = 1.0;
% psa.model_defaults.activation = 'logistic';   % or 'piecewise' / 'tanh'
% psa.store_local_lya = true;                   % Store decimated local LLE time series
% psa.store_local_lya_dt = 0.01;                 % Time resolution for local_lya (seconds)

%% Configure conditions
% All four adaptation conditions are tested. Each configuration (set of grid
% parameter values + network seed) is simulated under ALL of them, allowing
% direct comparison on the same W.
%
% These come from the shared helper rather than from PSA's built-in defaults,
% which are spelled in SRNNModel2's n_a_E / n_b_E vocabulary and would be passed
% verbatim to whatever model_class is set -- silently wrong for any other class.
% See srnn_adaptation_conditions for how each regime is spelled per class. They
% come from the preset when there is one, because a preset may carry its own
% depression routes, which reach the model only through a condition.
if exist('master_conditions', 'var')
    psa.set_conditions(master_conditions);
else
    psa.set_conditions(srnn_adaptation_conditions(model_class));
end

%% Run the parameter space analysis
% This will:
%   1. Generate all parameter combinations
%   2. Randomize execution order
%   3. Run batched parfor with checkpoint files (resume-capable)
%   4. Consolidate results into per-condition MAT files
%
% Results are automatically saved during execution

psa.run();

% Copy this script to the output directory for reproducibility
copyfile([mfilename('fullpath') '.m'], psa.output_dir);

%% Save the PSA object
save_file = fullfile(psa.output_dir, 'psa_object.mat');
save(save_file, 'psa');
fprintf('PSA object saved to: %s\n', save_file);

%% Plot results
% Generate histograms showing metric distributions across the parameter space

% Colour by the fraction-excitatory axis. This has to be named explicitly for
% SRNNCellTypePairs: plot's default color_by is 'f', which there is a 1 x C row
% and blows up the scalar assignment the histogram colouring does.
psa.plot('metric', 'LLE',       'color_by', f_param);
psa.plot('metric', 'mean_rate', 'color_by', f_param);

if save_figs
    fig_dir = fullfile(psa.output_dir, 'figures');
    save_some_figs_to_folder_2(fig_dir, 'param_space', [], {'fig', 'png'});
    fprintf('Figures saved to %s\n', fig_dir);
end

%% Display summary
fprintf('\n=== Parameter Space Analysis Summary ===\n');
fprintf('Output directory: %s\n', psa.output_dir);
fprintf('Grid parameters: %s\n', strjoin(psa.grid_params, ', '));
fprintf('Levels per parameter: %d\n', psa.n_levels);
fprintf('Total combinations: %d^%d = %d\n', ...
    psa.n_levels, length(psa.grid_params), psa.n_levels^length(psa.grid_params));
fprintf('Conditions: %s\n', strjoin(cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false), ', '));

%% Optional: Access results programmatically
% Results are stored in psa.results structure:
%
%   psa.results.<condition_name>{config_idx}
%
% Each result contains:
%   - success: boolean
%   - config: struct with parameter values
%   - config_idx: index in the grid
%   - network_seed: RNG seed for W (same across conditions)
%   - LLE: largest Lyapunov exponent
%   - mean_rate: mean firing rate
%
% Example: Extract LLE values for 'sfa_only' condition:
%
%   results = psa.results.sfa_only;
%   LLEs = cellfun(@(r) r.LLE, results(~cellfun(@isempty, results)));
%   histogram(LLEs);

%% Optional: Load results from a previous run
% psa_loaded = ParamSpaceAnalysis2();
% psa_loaded.load_results('/path/to/param_space_...');
% psa_loaded.plot('metric', 'LLE');

fprintf('\nDone! Results saved to: %s\n', psa.output_dir);
