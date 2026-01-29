% run_param_space_analysis.m
% Example script demonstrating how to use the ParamSpaceAnalysis class
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
% See also: ParamSpaceAnalysis, SRNNModel, SensitivityAnalysis

clear all;
clc;
close all;

%% Setup paths
setup_paths();

%% Create ParamSpaceAnalysis object
% Configure the analysis parameters
% NOTE: Total simulations = n_levels^(n_params) * n_conditions
%       e.g., 4 levels x 3 params x 4 conditions = 256 simulations

psa = ParamSpaceAnalysis(...
    'n_levels', 5, ...          % Number of levels per parameter
    'batch_size', 25, ...       % Configs per batch (for checkpointing)
    'note', 'test_refactor', ...         % Optional note for folder naming
    'verbose', true ...         % Print progress during execution
    );

%% Add parameters to the grid
% All combinations of these parameters will be tested
% The order in which parameters are added doesn't matter

% Network structure parameters
N = 300;
indegree = 100;
alpha = indegree / N;

default_tilde_val = 1 / sqrt(N * alpha * (2 - alpha));

% psa.add_grid_parameter('E_W', [-2 2] .* default_tilde_val);  % Mean offset (scaled by 1/sqrt(n))
psa.add_grid_parameter('f', [0.4 0.6]);     % fraction of neurons that are E

% Repetition index (creates unique network seeds per parameter combo)

psa.add_grid_parameter('rep_idx', [1, 2]);

% Dynamics parameters (uncomment to include)
% psa.add_grid_parameter('tau_d', [0.05, 0.2]);           % Dendritic time constant
% psa.add_grid_parameter('u_ex_scale', [1.0, 3.0]);       % External input scaling

% Adaptation parameters (uncomment to include)
% psa.add_grid_parameter('c_E', [0.01, 0.5]);             % SFA strength

%% Configure model defaults (optional)
% Set any SRNNModel properties that should be constant across all runs
% Parameters tuned to edge-of-chaos from full_SRNN_run_SRNNModel.m

% Network architecture
psa.model_defaults.n = N;                     % Number of neurons
psa.model_defaults.indegree = indegree;            % Sparse connectivity
% psa.model_defaults.f = 0.5;                   % E/I fraction

% Timing
psa.model_defaults.T_range = [-15, 45];       % Match full_SRNN_run_SRNNModel timing
psa.model_defaults.fs = 100;                  % Sampling frequency
psa.model_defaults.tau_d = 0.1;               % Dendritic time constant

% RMT tilde-parameters (Harris 2023)
psa.model_defaults.mu_E_tilde = 3.5 * default_tilde_val;
psa.model_defaults.mu_I_tilde = -3.5 * default_tilde_val;
psa.model_defaults.sigma_E_tilde = default_tilde_val;
psa.model_defaults.sigma_I_tilde = default_tilde_val;
psa.model_defaults.E_W = -0.5 / sqrt(N * alpha * (2 - alpha));
% psa.model_defaults.zrs_mode = 'Partial_SZRS';
psa.model_defaults.zrs_mode = 'none';
psa.model_defaults.level_of_chaos = 1.0;      % Edge of chaos
psa.model_defaults.rescale_by_abscissa = false;

% Adaptation parameters
psa.model_defaults.c_E = 0.25/3;              % SFA strength
psa.model_defaults.tau_a_E = logspace(log10(0.1), log10(10), 3);  % Match full_SRNN_run range
psa.model_defaults.tau_b_E_rec = 1;           % STD recovery time for E neurons
psa.model_defaults.tau_b_E_rel = 0.5;         % STD release time for E neurons

% Activation function
S_a = 0.9;
S_c = 0.40;
psa.model_defaults.S_a = S_a;
psa.model_defaults.S_c = S_c;
psa.model_defaults.activation_function = @(x) piecewiseSigmoid(x, S_a, S_c);
psa.model_defaults.activation_function_derivative = @(x) piecewiseSigmoidDerivative(x, S_a, S_c);

% Input configuration
psa.model_defaults.u_ex_scale = 1.0;          % Match full_SRNN_run_SRNNModel
psa.model_defaults.input_config.n_steps = 3;
psa.model_defaults.input_config.positive_only = true;
psa.model_defaults.input_config.step_density = 0.15;    % Fallback density (required)
psa.model_defaults.input_config.step_density_E = 0.15;
psa.model_defaults.input_config.step_density_I = 0;
psa.model_defaults.input_config.amp = 0.5;
psa.model_defaults.input_config.no_stim_pattern = [true false true];  % Stim on 2nd period only
psa.model_defaults.input_config.intrinsic_drive = zeros(N, 1);  % Required field

% Lyapunov settings
psa.model_defaults.lya_method = 'benettin';
psa.store_local_lya = true;                   % Store decimated local LLE time series
psa.store_local_lya_dt = 0.1;                 % Time resolution for local_lya (seconds)

%% Configure conditions (optional)
% By default, all four adaptation conditions are tested:
%   - no_adaptation: n_a_E=0, n_b_E=0
%   - sfa_only:      n_a_E=3, n_b_E=0
%   - std_only:      n_a_E=0, n_b_E=1
%   - sfa_and_std:   n_a_E=3, n_b_E=1
%
% Each configuration (set of grid parameter values + network seed) is
% simulated under ALL conditions, allowing direct comparison.
%
% To use only a subset of conditions:
%
% psa.set_conditions({
%     struct('name', 'no_adaptation', 'n_a_E', 0, 'n_b_E', 0),
%     struct('name', 'sfa_and_std',   'n_a_E', 3, 'n_b_E', 1)
% });

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

psa.plot('metric', 'LLE');
psa.plot('metric', 'mean_rate');
load_and_plot_lle_by_stim_period(psa.output_dir, 'transient_skip', 5, 'periods_to_plot', [0 1 1]);

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
% psa_loaded = ParamSpaceAnalysis();
% psa_loaded.load_results('/path/to/param_space_...');
% psa_loaded.plot('metric', 'LLE');

fprintf('\nDone! Results saved to: %s\n', psa.output_dir);
