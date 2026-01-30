%% run_parameter_exploration.m
% Parameter Exploration Script for Reservoir Dynamics
%
% This script performs comprehensive parameter space exploration using
% Latin Hypercube Sampling (LHS) to efficiently cover the widest parameter
% space. It detects unstable and dead simulations, supports GPU acceleration,
% and provides optional result saving.
%
% Usage:
%   run_parameter_exploration
%   run_parameter_exploration(config)
%
% See also: ParameterExploration, SRNN_ESN

clear;
clc;
close all;

%% Setup paths
setup_paths();

%% Configuration
fprintf('=================================================================\n');
fprintf('  Parameter Exploration for Reservoir Dynamics\n');
fprintf('=================================================================\n\n');

% Create configuration struct
config = struct();

% === Exploration Settings ===
config.n_samples = 250;              % Number of LHS samples
config.use_gpu = true;                 % Use GPU if available (hybrid mode)
config.save_results = false;            % Save results to disk
config.checkpoint_interval = 100;      % Save checkpoint every N samples
config.batch_size = 8;                % Batch size for parfor
config.verbose = true;                 % Print progress
config.random_seed = 42;               % Random seed for reproducibility
config.note = 'comprehensive';         % Optional note for output folder

% === Learning Evaluation Settings ===
config.evaluate_learning = true;       % Enable learning capability evaluation
config.learning_task = 'prediction';   % Task: 'prediction', 'narma10', 'memory'
config.learning_n_train = 2000;        % Training samples for learning task
config.learning_n_test = 500;          % Test samples for learning task
config.learning_washout = 100;         % Washout steps for learning

% === Parameter Ranges ===
% Define ranges for all parameters to explore
config.param_ranges = getDefaultParamRanges();

% === Fixed Parameters ===
% Parameters that remain constant across all simulations
config.fixed_params = struct();
config.fixed_params.dt = 0.1;                    % Integration time step
config.fixed_params.n_steps = 300;              % Number of time steps to stimulate
config.fixed_params.N = 10;                      % Silence before stimulus (seconds)
config.fixed_params.M = 200;                    % Silence after stimulus (seconds)
config.fixed_params.input_amplitude = 1;         % Input amplitude
config.fixed_params.S_a = 0.85;                 % Activation function parameter
config.fixed_params.S_c = 0.4;                  % Activation function parameter

% === Override from command line or user modifications ===
% Users can modify config here or pass it as input
% To customize, modify config struct above or create a custom config and pass it:
%   custom_config = config;
%   custom_config.n_samples = 500;
%   run_parameter_exploration(custom_config);

%% Display Configuration
fprintf('Configuration:\n');
fprintf('  Number of samples: %d\n', config.n_samples);
fprintf('  Batch size: %d\n', config.batch_size);
fprintf('  Checkpoint interval: %d\n', config.checkpoint_interval);
fprintf('  GPU computing: %s\n', mat2str(config.use_gpu));
fprintf('  Save results: %s\n', mat2str(config.save_results));
fprintf('  Random seed: %d\n', config.random_seed);
fprintf('\n');

fprintf('Parameters to explore:\n');
param_names = fieldnames(config.param_ranges);
for i = 1:length(param_names)
    param_name = param_names{i};
    param_range = config.param_ranges.(param_name);
    if iscell(param_range)
        fprintf('  %s: %s\n', param_name, strjoin(param_range, ', '));
    else
        fprintf('  %s: [%.3g, %.3g]\n', param_name, param_range(1), param_range(2));
    end
end
fprintf('\n');

%% Run Exploration
fprintf('Starting parameter exploration...\n\n');

explorer = ParameterExploration(config);
results = explorer.runExploration();

%% Display Summary
fprintf('\n=================================================================\n');
fprintf('  EXPLORATION SUMMARY\n');
fprintf('=================================================================\n');
fprintf('Total samples: %d\n', config.n_samples);
fprintf('Success rate: %.1f%%\n', results.success_rate * 100);
fprintf('Unstable rate: %.1f%%\n', results.unstable_rate * 100);
fprintf('Dead rate: %.1f%%\n', results.dead_rate * 100);
fprintf('Stable rate: %.1f%%\n', results.stable_rate * 100);
fprintf('\n');

if isfield(results, 'metrics_stats')
    fprintf('Key Dynamics Metrics (mean ± std):\n');
    if isfield(results.metrics_stats, 'mean_firing_rate')
        mfr = results.metrics_stats.mean_firing_rate;
        fprintf('  Mean firing rate: %.4f ± %.4f Hz\n', mfr.mean, mfr.std);
    end
    if isfield(results.metrics_stats, 'participation_ratio')
        pr = results.metrics_stats.participation_ratio;
        fprintf('  Participation ratio: %.2f ± %.2f\n', pr.mean, pr.std);
    end
    if isfield(results.metrics_stats, 'EI_balance')
        eib = results.metrics_stats.EI_balance;
        if isfinite(eib.mean)
            fprintf('  E/I balance: %.2f ± %.2f\n', eib.mean, eib.std);
        end
    end
    
    % Learning metrics
    fprintf('\nLearning Capability Metrics (mean ± std):\n');
    if isfield(results.metrics_stats, 'learning_test_nrmse')
        nrmse = results.metrics_stats.learning_test_nrmse;
        fprintf('  Test NRMSE: %.4f ± %.4f\n', nrmse.mean, nrmse.std);
    end
    if isfield(results.metrics_stats, 'learning_R2')
        r2 = results.metrics_stats.learning_R2;
        fprintf('  Test R²: %.4f ± %.4f\n', r2.mean, r2.std);
    end
    if isfield(results.metrics_stats, 'memory_capacity')
        mc = results.metrics_stats.memory_capacity;
        fprintf('  Memory Capacity: %.2f ± %.2f\n', mc.mean, mc.std);
    end
end

% Display learning evaluation summary
if isfield(results, 'learning_evaluated_count')
    fprintf('\nLearning Evaluation:\n');
    fprintf('  Evaluated: %d / %d stable simulations\n', ...
        results.learning_evaluated_count, results.stable_count);
    if isfield(results, 'learning_success_count')
        fprintf('  Successful: %d (%.1f%%)\n', ...
            results.learning_success_count, ...
            100 * results.learning_success_count / max(1, results.learning_evaluated_count));
    end
end

fprintf('\n');

if config.save_results
    fprintf('Results saved to: %s\n', explorer.output_dir);
end

fprintf('=================================================================\n\n');

%% Generate Visualization Figures
fprintf('Generating visualization figures...\n\n');

fig_handles = explorer.plot(results, 'save_figs', config.save_results);

% Display figures
fprintf('\nFigures generated:\n');
fig_names = fieldnames(fig_handles);
for i = 1:length(fig_names)
    fprintf('  - %s\n', fig_names{i});
end

if config.save_results
    fprintf('\nFigures saved to: %s/figures/\n', explorer.output_dir);
end

fprintf('\n=================================================================\n');
fprintf('  Parameter exploration complete!\n');
fprintf('=================================================================\n\n');

%% Helper Function: Get Default Parameter Ranges
function param_ranges = getDefaultParamRanges()
    %GETDEFAULTPARAMRANGES Get default parameter ranges for exploration
    %
    % Returns a struct with parameter ranges covering the widest space possible
    
    param_ranges = struct();
    
    % === Architecture Parameters ===
    param_ranges.n = [50, 250];                    % Total number of neurons
    param_ranges.fraction_E = [0.3, 0.7];         % Fraction of excitatory neurons
    
    % === Dynamics Parameters ===
    param_ranges.level_of_chaos = [0.5, 4.0];     % Chaos level
    param_ranges.tau_d = [0.1, 1.0];              % Dendritic time constant
    param_ranges.lags = [0.01, 0.05];             % Synaptic delay
    
    % === Adaptation Parameters ===
    param_ranges.n_a_E = [0, 3];                  % Number of E adaptation timescales
    param_ranges.n_a_I = [0, 3];                  % Number of I adaptation timescales
    param_ranges.c_E = [0, 0.5];                  % E adaptation strength
    param_ranges.c_I = [0, 0.5];                  % I adaptation strength
    
    % === Short-term Depression Parameters ===
    param_ranges.n_b_E = [0, 1];                  % Enable STD for E (0 or 1)
    param_ranges.n_b_I = [0, 1];                  % Enable STD for I (0 or 1)
    param_ranges.tau_b_E_rec = [0.1, 1.0];       % E STD recovery time
    param_ranges.tau_b_E_rel = [0.05, 0.5];      % E STD release time
    param_ranges.tau_b_I_rec = [0.1, 1.0];       % I STD recovery time
    param_ranges.tau_b_I_rel = [0.05, 0.5];     % I STD release time
    
    % === Activation Function Parameters ===
    param_ranges.S_a = [0.5, 1.0];                % Activation function parameter
    param_ranges.S_c = [0.2, 0.6];                % Activation function parameter
    
    % === Input Configuration ===
    param_ranges.input_scaling = [0.1, 2.0];      % Input weight scaling
    param_ranges.input_type = {'pulse', 'sine', 'step', 'noise', 'mackey'};  % Input type (categorical)
end

