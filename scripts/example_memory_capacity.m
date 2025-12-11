%% Example: Memory Capacity Measurement with SRNN_ESN_reservoir
% This script demonstrates how to use the SRNN_ESN_reservoir class to measure
% memory capacity under different adaptation conditions (baseline, SFA, SFA+STD).
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
n = 200;                    % Number of neurons
level_of_chaos = 1.8;       % Moderate chaos level
rng_seed_net = 42;          % Fixed seed for reproducibility

% MC protocol parameters (reduced for faster demo)
T_wash = 200*10;               % Washout samples
T_train = 200*60;             % Training samples  
T_test = 200*60;              % Test samples
d_max = 200*5;                % Maximum delay

% Input type: 'white' (standard ESN) or 'bandlimited' (fair for systems with tau_d)
% Bandlimited uses low-pass filtered noise matching the system bandwidth
input_type = 'bandlimited';
u_f_cutoff = 5;               % Cutoff frequency for bandlimited input (Hz)

%% Condition 1: Baseline (no SFA, no STD)
fprintf('\n==============================\n');
fprintf('CONDITION 1: Baseline (no adaptation)\n');
fprintf('==============================\n');

esn_baseline = SRNN_ESN_reservoir(...
    'n', n, ...
    'level_of_chaos', level_of_chaos, ...
    'rng_seeds', [rng_seed_net, 123], ...
    'tau_d', 0.025, ...        % Dendritic time constant (s)
    'n_a_E', 0, ...            % No SFA for E neurons
    'n_a_I', 0, ...            % No SFA for I neurons
    'n_b_E', 0, ...            % No STD for E neurons
    'n_b_I', 0, ...            % No STD for I neurons
    'input_type', input_type, ...
    'u_f_cutoff', u_f_cutoff, ...
    'T_wash', T_wash, ...
    'T_train', T_train, ...
    'T_test', T_test, ...
    'd_max', d_max);

esn_baseline.build();
[MC_baseline, R2_baseline, results_baseline] = esn_baseline.run_memory_capacity();

%% Condition 2: SFA only
fprintf('\n==============================\n');
fprintf('CONDITION 2: SFA only\n');
fprintf('==============================\n');

esn_sfa = SRNN_ESN_reservoir(...
    'n', n, ...
    'level_of_chaos', level_of_chaos, ...
    'rng_seeds', [rng_seed_net, 123], ...  % Same seeds as baseline
    'tau_d', 0.025, ...        % Dendritic time constant (s)
    'n_a_E', 3, ...            % 3 adaptation timescales for E neurons
    'tau_a_E', [0.1, 1.0, 10], ... % Adaptation time constants (s)
    'n_a_I', 0, ...            % No SFA for I neurons
    'n_b_E', 0, ...            % No STD
    'n_b_I', 0, ...
    'c_E', 0.5, ...            % Adaptation strength
    'input_type', input_type, ...
    'u_f_cutoff', u_f_cutoff, ...
    'T_wash', T_wash, ...
    'T_train', T_train, ...
    'T_test', T_test, ...
    'd_max', d_max);

esn_sfa.build();
[MC_sfa, R2_sfa, results_sfa] = esn_sfa.run_memory_capacity();

%% Condition 3: SFA + STD
fprintf('\n==============================\n');
fprintf('CONDITION 3: SFA + STD\n');
fprintf('==============================\n');

esn_full = SRNN_ESN_reservoir(...
    'n', n, ...
    'level_of_chaos', level_of_chaos, ...
    'rng_seeds', [rng_seed_net, 123], ...  % Same seeds
    'tau_d', 0.025, ...        % Dendritic time constant (s)
    'n_a_E', 3, ...            % SFA for E neurons
    'tau_a_E', [0.1, 1.0, 10], ... % Adaptation time constants (s)
    'n_a_I', 0, ...
    'n_b_E', 1, ...            % STD for E neurons
    'n_b_I', 0, ...
    'c_E', 0.5, ...
    'tau_b_E_rec', 1.0, ...    % STD recovery time
    'tau_b_E_rel', 0.25, ...   % STD release time
    'input_type', input_type, ...
    'u_f_cutoff', u_f_cutoff, ...
    'T_wash', T_wash, ...
    'T_train', T_train, ...
    'T_test', T_test, ...
    'd_max', d_max);

esn_full.build();
[MC_full, R2_full, results_full] = esn_full.run_memory_capacity();

%% Summary and Comparison Plot
fprintf('\n==============================\n');
fprintf('SUMMARY\n');
fprintf('==============================\n');
fprintf('Memory Capacity Results:\n');
fprintf('  Baseline (no adaptation): MC = %.2f\n', MC_baseline);
fprintf('  SFA only:                 MC = %.2f\n', MC_sfa);
fprintf('  SFA + STD:                MC = %.2f\n', MC_full);

% Create comparison figure
figure('Position', [100, 100, 1200, 500]);

% Plot 1: R^2 vs delay for all conditions
subplot(1, 3, 1);
hold on;
bar_data = [R2_baseline; R2_sfa; R2_full]';
b = bar(1:d_max, bar_data);
b(1).FaceColor = [0.7, 0.7, 0.7];  % Gray for baseline
b(2).FaceColor = [0.3, 0.6, 0.9];  % Blue for SFA
b(3).FaceColor = [0.9, 0.4, 0.3];  % Red for SFA+STD
xlabel('Delay d (samples)');
ylabel('R^2_d');
title('Memory Capacity by Delay');
legend('Baseline', 'SFA only', 'SFA+STD', 'Location', 'northeast');
grid on;
hold off;

% Plot 2: Cumulative MC
subplot(1, 3, 2);
hold on;
plot(1:d_max, cumsum(R2_baseline), 'k-', 'LineWidth', 2, 'DisplayName', 'Baseline');
plot(1:d_max, cumsum(R2_sfa), 'b-', 'LineWidth', 2, 'DisplayName', 'SFA only');
plot(1:d_max, cumsum(R2_full), 'r-', 'LineWidth', 2, 'DisplayName', 'SFA+STD');
xlabel('Delay d (samples)');
ylabel('Cumulative MC');
title('Cumulative Memory Capacity');
legend('Location', 'southeast');
grid on;
hold off;

% Plot 3: Bar chart of total MC
subplot(1, 3, 3);
bar_heights = [MC_baseline, MC_sfa, MC_full];
b = bar(1:3, bar_heights);
b.FaceColor = 'flat';
b.CData(1,:) = [0.7, 0.7, 0.7];
b.CData(2,:) = [0.3, 0.6, 0.9];
b.CData(3,:) = [0.9, 0.4, 0.3];
set(gca, 'XTickLabel', {'Baseline', 'SFA', 'SFA+STD'});
ylabel('Total Memory Capacity');
title('Total MC Comparison');
grid on;

% Add value labels on bars
for i = 1:3
    text(i, bar_heights(i) + 0.5, sprintf('%.1f', bar_heights(i)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

sgtitle('Memory Capacity Analysis: Effect of Spike-Frequency Adaptation and Short-Term Depression');

%% Time Series Plots for Each Condition
delays_to_plot = [1, 50, 100, 200, 400, 800, 1600, 3200];

fprintf('\nGenerating time series plots...\n');

% Baseline time series
esn_baseline.plot_esn_timeseries(delays_to_plot, 'title', 'Baseline (No Adaptation)');

% SFA only time series
esn_sfa.plot_esn_timeseries(delays_to_plot, 'title', 'SFA Only');

% SFA + STD time series  
esn_full.plot_esn_timeseries(delays_to_plot, 'title', 'SFA + STD');

fprintf('Time series plots generated.\n');

%% Save results
results = struct();
results.baseline = results_baseline;
results.sfa = results_sfa;
results.full = results_full;
results.MC = [MC_baseline, MC_sfa, MC_full];
results.conditions = {'Baseline', 'SFA only', 'SFA+STD'};

save('memory_capacity_results.mat', 'results');
fprintf('\nResults saved to memory_capacity_results.mat\n');
