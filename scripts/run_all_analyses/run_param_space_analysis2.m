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
switch run_mode
    case 'fast'
        % Fast iteration: fewer levels, half the sample rate, short time range.
        % fs=200 keeps Benettin's lya_dt/dt guard satisfied (4>=3);
        % T_range=[0,20] with an explicit 10 s LLE window [10,20].
        n_levels = 3; ode_solver_mode = @ode_rk4;
        fs_mode = 200;  T_range_mode = [0, 20];  lya_T_interval_mode = [10, 20];
    case 'medium'
        % Medium: roughly halfway between fast and production. ode45 at fs=400,
        % T_range=[0,30] with a 15 s LLE window [15,30].
        % n_levels=4 (this is a multi-dim grid with no reps axis, so 4^n_params
        % configs) sits between fast (3) and production (5).
        n_levels = 4; ode_solver_mode = @ode45;
        fs_mode = 400;  T_range_mode = [0, 30];  lya_T_interval_mode = [15, 30];
    case 'production'
        n_levels = 5; ode_solver_mode = @ode45;
        fs_mode = 400;  T_range_mode = [0, 50];  lya_T_interval_mode = [];
    otherwise, error('run_param_space_analysis2:badMode', ...
        'Unknown run_mode ''%s'' (expected ''fast'', ''medium'', or ''production'').', run_mode);
end
fprintf('[run_param_space_analysis2] run_mode=%s, n_levels=%d, ode_solver=%s, fs=%d, T_range=[%g %g]\n', ...
    run_mode, n_levels, func2str(ode_solver_mode), fs_mode, T_range_mode(1), T_range_mode(2));

psa = ParamSpaceAnalysis2(...
    'n_levels', n_levels, ...   % set by run_mode (fast=3, production=5)
    'batch_size', 25, ...       % configs per batch (for checkpointing)
    'note', 'test_refactor', ...         % Optional note for folder naming
    'verbose', true ...         % Print progress during execution
    );
if exist('master_output_dir', 'var')
    psa.output_dir = master_output_dir;
end
psa.model_defaults.ode_solver = ode_solver_mode;  % fast=ode_rk4, production=ode45
psa.model_defaults.fs = fs_mode;                   % fast=200 (default 400)
psa.model_defaults.T_range = T_range_mode;         % fast=[0,20] (default [0,50])
if ~isempty(lya_T_interval_mode)
    psa.model_defaults.lya_T_interval = lya_T_interval_mode;  % fast: LLE window [10,20]
end

%% Add parameters to the grid
% All combinations of these parameters will be tested
% The order in which parameters are added doesn't matter

% Network structure parameters
% Same swept variables and ranges as run_sensitivity_analysis.m (n, f,
% level_of_chaos), but here as a joint multi-dimensional grid.
psa.add_grid_parameter('n',              [100, 1000]);   % network size
psa.add_grid_parameter('f',              [0.25, 0.75]);  % fraction excitatory
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
% psa.model_defaults.activation_function = @logisticSigmoid;
% psa.model_defaults.activation_function_derivative = @logisticSigmoidDerivative;
% psa.store_local_lya = true;                   % Store decimated local LLE time series
% psa.store_local_lya_dt = 0.01;                 % Time resolution for local_lya (seconds)

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
