%% Compute half of the Stability_Manuscript example memory-capacity figure.
% Runs the (expensive) memory-capacity simulation for the 4 adaptation
% conditions and saves the results to mc_example_data.mat in this folder. The
% plotting is done separately by Fig_memory_capacity_example.m, so the figure
% look can be iterated without re-simulating.
%
% Model settings match looped_memory_capacity.m. See also: SRNN_ESN_reservoir.

clear; clc; close all;

% Assumes setup_paths.m has already been run (src/ + scripts/ on the MATLAB path).
% this_dir is only used to locate the output .mat next to this script.
this_dir = fileparts(mfilename('fullpath'));

%% Common parameters
n = 300;                    % Number of neurons
f = 0.6;                    % Fraction excitatory (off perfect E/I balance so no-adaptation isn't favored)
level_of_chaos = 2.5;         % Higher for the logistic: its mean slope < 1 pushes the edge of chaos to larger level_of_chaos
rng_seed_net = 42+15;          % Fixed seed for network reproducibility
rng_seed_stim = 43;         % Fixed seed for stimulus reproducibility

% Sampling frequency
fs = 200;                     % Sampling frequency (Hz)

% MC protocol parameters (defined in seconds, converted to samples).
% These are shortened for fast dial-in; increase for final results (in
% particular T_wash should exceed a few x the slowest tau_a_E for the SFA
% conditions -- here tau_a_E max = 10 s).
T_wash_sec = 10;              % Washout duration (s); discarded transient before training
T_train_sec = 300;           % Training duration (seconds); matches looped_memory_capacity
T_test_sec = 100;            % Test duration (seconds); matches looped_memory_capacity

T_wash = T_wash_sec * fs;     % Washout samples
T_train = T_train_sec * fs;   % Training samples
T_test = T_test_sec * fs;     % Test samples
d_max = 15*fs;                 % Maximum delay

% Input type: 'white' (standard ESN), 'bandlimited' (fair for systems with tau_d),
%             'one_over_f' (1/f^alpha noise), or 'sample_hold' (i.i.d. values held
%             for ~a few tau_d; MC is then measured in hold units -- fair for a
%             low-pass reservoir and free of autocorrelation inflation).
input_type = 'sample_hold'; % options: 'white', 'bandlimited', 'one_over_f', 'sample_hold'
u_f_cutoff = 5;               % Cutoff frequency for bandlimited input (Hz)
u_alpha = 1;                  % Spectral exponent for 1/f^alpha noise (1=pink, 2=red/Brownian)
tau_d = 0.1;                  % Dendritic time constant (s)
T_hold = 3 * tau_d;           % sample_hold: hold each i.i.d. input value this long (= 3*tau_d)

%% Base ESN configuration (shared across all conditions)
% All parameters here are identical for every condition. Condition-specific
% overrides (n_a_E, n_b_E) are applied below via the condition_args cell array.
base_args = { ...
    'n', n, ...
    'f', f, ...                    % fraction excitatory (0.6; off perfect balance)
    'fs', fs, ...
    'level_of_chaos', level_of_chaos, ...
    'ode_solver', @ode_rk4, ...    % fixed-step solver (fast); default is @ode45
    'rng_seeds', [rng_seed_net, rng_seed_stim], ...
    'tau_d', tau_d, ...            % Dendritic time constant (s)
    'S_c', 0.35, ...               % Nonlinearity bias (center); matches SRNNModel2 default
    'S_a', 0.9, ...                % Fraction of nonlinearity with slope 1 (unused by the logistic default)
    'n_a_I', 0, ...                % No SFA for I neurons (all conditions)
    'n_b_I', 0, ...                % No STD for I neurons (all conditions)
    'c_E', 0.5/3, ...              % Adaptation strength for E neurons (matches looped_memory_capacity)
    'tau_a_E', [0.1, 1.0, 10], ... % Adaptation time constants (s)
    'tau_b_E_rec', 1.0, ...        % STD recovery time constant (s)
    'tau_b_E_rel', 0.25, ...       % STD release time constant (s)
    'std_zero_floor', true, ...    % rescale STD so (prod b)*r reaches 0 at full depression (matches run_all)
    'input_type', input_type, ...
    'u_f_cutoff', u_f_cutoff, ...
    'u_alpha', u_alpha, ...
    'T_hold', T_hold, ...          % hold duration for input_type='sample_hold'
    'T_wash', T_wash, ...
    'T_train', T_train, ...
    'T_test', T_test, ...
    'd_max', d_max};

%% Condition-specific overrides (only what differs)
condition_names = {'Baseline (no adaptation)', 'SFA only', 'STD only', 'SFA + STD'};
condition_args = { ...
    {'n_a_E', 0, 'n_b_E', 0}, ...   % Baseline: no adaptation
    {'n_a_E', 3, 'n_b_E', 0}, ...   % SFA only: 3 adaptation timescales
    {'n_a_E', 0, 'n_b_E', 1}, ...   % STD only: short-term depression, no SFA
    {'n_a_E', 3, 'n_b_E', 1}, ...   % SFA + STD: adaptation + depression
    };
n_cond = numel(condition_names);

%% Build all conditions
esn = cell(1, n_cond);
for i = 1:n_cond
    fprintf('\n==============================\n');
    fprintf('Building %s...\n', condition_names{i});
    fprintf('==============================\n');
    esn{i} = SRNN_ESN_reservoir(base_args{:}, condition_args{i}{:});
    esn{i}.build();
end

%% Verify shared structure (Option E)
% Asserts that W, W_in, u_scalar, u_ex, t_ex are identical across conditions,
% all public config properties match, and n_a_E/n_b_E/tau_a_E actually differ.
SRNN_ESN_reservoir.verify_shared_build(esn, ...
    {'n_a_E', 'n_b_E'}, ...
    {'W', 'W_in', 'u_scalar', 'u_ex', 't_ex'});

%% Run all conditions
MC = zeros(1, n_cond);
R2 = cell(1, n_cond);
results = cell(1, n_cond);   % per-condition mc_results structs (public getter)
for i = 1:n_cond
    fprintf('\n==============================\n');
    fprintf('CONDITION %d: %s\n', i, condition_names{i});
    fprintf('==============================\n');
    [MC(i), R2{i}, results{i}] = esn{i}.run_memory_capacity();
end

%% Summary
fprintf('\n==============================\n');
fprintf('SUMMARY\n');
fprintf('==============================\n');
fprintf('Memory Capacity Results:\n');
for i = 1:n_cond
    fprintf('  %s: MC = %.2f\n', condition_names{i}, MC(i));
end

%% Save results for the plotting script
% Saves the per-condition mc_results structs (which carry input, predictions,
% targets, and full state time series) plus the summary + config, so
% Fig_memory_capacity_example.m can re-render without re-simulating. The ESN
% objects themselves are NOT saved (S_out alone is ~GB; mc_results has all the
% plotting data). File is gitignored.
delay_s = results{1}.delay_s;   % delay axis in seconds (shared across conditions)
save_file = fullfile(this_dir, 'mc_example_data.mat');
save(save_file, 'results', 'MC', 'R2', 'delay_s', 'condition_names', ...
    'base_args', 'condition_args', '-v7.3');
fprintf('\nSaved: %s\n', save_file);
