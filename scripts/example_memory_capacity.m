%% Example: Memory Capacity Measurement with SRNN_ESN_reservoir
% This script demonstrates how to use the SRNN_ESN_reservoir class to measure
% memory capacity under different adaptation conditions (baseline, SFA, SFA+STD).
%
% Uses Option A (base config + condition-specific overrides) to avoid
% duplicating configuration, and Option E (verify_shared_build) to assert
% that all conditions share identical network structure.
%
% The memory capacity protocol:
%   1. Drive reservoir with scalar random input u(t) ~ U(0,1)
%   2. Train linear readouts for each delay d to reconstruct u(t-d)
%   3. Compute R^2_d and sum to get total memory capacity
%
% See also: SRNN_ESN_reservoir, Memory_capacity_protocol.md

clear; clc; close all;

%% Add paths
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'src')));

%% Common parameters
n = 300;                    % Number of neurons
level_of_chaos = 1.0;       % Moderate chaos level
rng_seed_net = 42;          % Fixed seed for network reproducibility
rng_seed_stim = 43;         % Fixed seed for stimulus reproducibility

% Sampling frequency
fs = 200;                     % Sampling frequency (Hz)

% MC protocol parameters (defined in seconds, converted to samples)
T_wash_sec = 20;              % Washout duration (seconds)
T_train_sec = 50;            % Training duration (seconds)
T_test_sec = 50;             % Test duration (seconds)

T_wash = T_wash_sec * fs;     % Washout samples
T_train = T_train_sec * fs;   % Training samples
T_test = T_test_sec * fs;     % Test samples
d_max = 3*fs;                 % Maximum delay

% Input type: 'white' (standard ESN), 'bandlimited' (fair for systems with tau_d),
%             or 'one_over_f' (1/f^alpha noise, mimics SEEG/EEG power spectrum)
% Bandlimited uses low-pass filtered noise matching the system bandwidth
input_type = 'bandlimited';
u_f_cutoff = 5;               % Cutoff frequency for bandlimited input (Hz)
u_alpha = 1;                  % Spectral exponent for 1/f^alpha noise (1=pink, 2=red/Brownian)

%% Base ESN configuration (shared across all conditions)
% All parameters here are identical for every condition. Condition-specific
% overrides (n_a_E, n_b_E) are applied below via the condition_args cell array.
% Setting tau_a_E, c_E, tau_b_E_rec, tau_b_E_rel here is harmless when
% n_a_E=0 or n_b_E=0 — those parameters are simply ignored by the dynamics.
base_args = { ...
    'n', n, ...
    'fs', fs, ...
    'level_of_chaos', level_of_chaos, ...
    'rng_seeds', [rng_seed_net, rng_seed_stim], ...
    'tau_d', 0.1, ...              % Dendritic time constant (s)
    'S_c', 0.4, ...                % Nonlinearity bias (center)
    'S_a', 0.9, ...                % Fraction of nonlinearity with slope 1
    'n_a_I', 0, ...                % No SFA for I neurons (all conditions)
    'n_b_I', 0, ...                % No STD for I neurons (all conditions)
    'c_E', 0.15/3, ...             % Adaptation strength for E neurons
    'tau_a_E', [0.1, 1.0, 10], ... % Adaptation time constants (s)
    'tau_b_E_rec', 1.0, ...        % STD recovery time constant (s)
    'tau_b_E_rel', 0.25, ...       % STD release time constant (s)
    'input_type', input_type, ...
    'u_f_cutoff', u_f_cutoff, ...
    'u_alpha', u_alpha, ...
    'T_wash', T_wash, ...
    'T_train', T_train, ...
    'T_test', T_test, ...
    'd_max', d_max};

%% Condition-specific overrides (only what differs)
condition_names = {'Baseline (no adaptation)', 'SFA only', 'SFA + STD'};
condition_args = { ...
    {'n_a_E', 0, 'n_b_E', 0}, ...   % Baseline: no adaptation
    {'n_a_E', 3, 'n_b_E', 0}, ...   % SFA only: 3 adaptation timescales
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
verify_shared_build(esn, ...
    {'n_a_E', 'n_b_E'}, ...
    {'W', 'W_in', 'u_scalar', 'u_ex', 't_ex'});

%% Run all conditions
MC = zeros(1, n_cond);
R2 = cell(1, n_cond);
results = cell(1, n_cond);
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

%% Comparison Plot
colors = [0.7, 0.7, 0.7;   % Gray for baseline
    0.3, 0.6, 0.9;   % Blue for SFA
    0.9, 0.4, 0.3];  % Red for SFA+STD

figure();

% Plot 1: R^2 vs delay for all conditions
subplot(1, 3, 1);
hold on;
bar_data = zeros(d_max, n_cond);
for i = 1:n_cond
    bar_data(:, i) = R2{i}(:);
end
b = bar(1:d_max, bar_data);
for i = 1:n_cond
    b(i).FaceColor = colors(i, :);
end
xlabel('Delay d (samples)');
ylabel('R^2_d');
title('Memory Capacity by Delay');
legend(condition_names, 'Location', 'northeast');
grid on;
hold off;

% Plot 2: Cumulative MC
subplot(1, 3, 2);
hold on;
line_styles = {'k-', 'b-', 'r-'};
for i = 1:n_cond
    plot(1:d_max, cumsum(R2{i}), line_styles{i}, 'LineWidth', 2, 'DisplayName', condition_names{i});
end
xlabel('Delay d (samples)');
ylabel('Cumulative MC');
title('Cumulative Memory Capacity');
legend('Location', 'southeast');
grid on;
hold off;

% Plot 3: Bar chart of total MC
subplot(1, 3, 3);
b = bar(1:n_cond, MC);
b.FaceColor = 'flat';
for i = 1:n_cond
    b.CData(i, :) = colors(i, :);
end
set(gca, 'XTickLabel', condition_names);
ylabel('Total Memory Capacity');
title('Total MC Comparison');
grid on;

% Add value labels on bars
for i = 1:n_cond
    text(i, MC(i) + 0.5, sprintf('%.1f', MC(i)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

sgtitle('Memory Capacity Analysis: Effect of Spike-Frequency Adaptation and Short-Term Depression');

%% Time Series Plots for Each Condition
delays_to_plot = [1, 50, 100, 200, 400, 800, 1600, 3200];

fprintf('\nGenerating time series plots...\n');

for i = 1:n_cond
    esn{i}.plot_esn_timeseries(delays_to_plot, 'title', condition_names{i});
end

fprintf('Time series plots generated.\n');

%% Save results
results_all = struct();
results_all.conditions = {condition_names};
results_all.MC = MC;
results_all.R2 = {R2};
for i = 1:n_cond
    results_all.(sprintf('cond%d', i)) = results{i};
end

% Save results
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
output_dir = fullfile(project_root, 'data', 'memory_capacity');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

save_file = fullfile(output_dir, 'memory_capacity_results.mat');
save(save_file, 'results_all');
fprintf('\nResults saved to memory_capacity_results.mat\n');
