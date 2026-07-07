%% Plot half of the Stability_Manuscript example memory-capacity figure.
% Loads mc_example_data.mat (written by compute_memory_capacity_example.m) and
% renders the figures -- no simulation is re-run, so the look can be iterated
% quickly. Run compute_memory_capacity_example.m first.
%
% Produces:
%   Fig 1 -- 1x3 comparison: R^2 by delay, cumulative MC, total-MC bars.
%   Fig 2 -- per-condition reconstruction: the test input and, for a few delays,
%            the target u(t-d) vs the readout y_d(t). Built from the saved
%            mc_results structs (predictions/t_pred/R2_d), a clean struct-based
%            stand-in for the diagnostic plot_esn_timeseries method.

clear; clc; close all;

% Assumes setup_paths.m has already been run (src/ + scripts/ on the MATLAB path).
% this_dir is only used to locate mc_example_data.mat next to this script.
this_dir = fileparts(mfilename('fullpath'));

%% Load precomputed results
data_file = fullfile(this_dir, 'mc_example_data.mat');
assert(isfile(data_file), ...
    'Missing %s -- run compute_memory_capacity_example.m first.', data_file);
S = load(data_file);
results         = S.results;
MC              = S.MC;
R2              = S.R2;
delay_s         = S.delay_s;
condition_names = S.condition_names;
n_cond          = numel(condition_names);

%% Plotting parameters (edit here to restyle)
colors = [0.7, 0.7, 0.7;   % Gray:  baseline
          0.3, 0.6, 0.9;   % Blue:  SFA only
          0.3, 0.7, 0.4;   % Green: STD only
          0.9, 0.4, 0.3];  % Red:   SFA + STD
delays_to_plot = [1, 2, 3, 5, 8];   % hold-delays to show in the reconstruction

%% ==================== Fig 1: MC comparison (1x3) ====================
figure('Color', 'w', 'Position', [100 100 1200 380]);
nd = numel(delay_s);

% Plot 1: R^2 vs delay for all conditions
subplot(1, 3, 1); hold on;
bar_data = zeros(nd, n_cond);
for i = 1:n_cond
    bar_data(:, i) = R2{i}(:);
end
b = bar(delay_s, bar_data);
for i = 1:n_cond
    b(i).FaceColor = colors(i, :);
end
xlabel('Delay (s)'); ylabel('R^2_d');
title('Memory Capacity by Delay');
legend(condition_names, 'Location', 'northeast');
grid on; hold off;

% Plot 2: Cumulative MC
subplot(1, 3, 2); hold on;
line_styles = {'k-', 'b-', 'g-', 'r-'};
for i = 1:n_cond
    plot(delay_s, cumsum(R2{i}), line_styles{i}, 'LineWidth', 2, 'DisplayName', condition_names{i});
end
xlabel('Delay (s)'); ylabel('Cumulative MC');
title('Cumulative Memory Capacity');
legend('Location', 'southeast');
grid on; hold off;

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
for i = 1:n_cond
    text(i, MC(i) + 0.5, sprintf('%.1f', MC(i)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end
sgtitle('Memory Capacity Analysis: Effect of Spike-Frequency Adaptation and Short-Term Depression');

%% ============ Fig 2: per-condition reconstruction (from struct) ============
% For each condition: the test input on top, then target u(t-d) vs readout
% y_d(t) for each delay -- rebuilt from the saved mc_results (no object needed).
for i = 1:n_cond
    mcr = results{i};
    d_avail = delays_to_plot(delays_to_plot <= numel(mcr.predictions));
    nD = numel(d_avail);

    figure('Color', 'w', 'Position', [120 120 720 130*(nD+1)]);
    tiledlayout(nD + 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Input being remembered
    nexttile;
    plot(mcr.t_pred, mcr.u_test, 'k-', 'LineWidth', 1);
    ylabel('u(t)'); title('Test input'); grid on;
    set(gca, 'XTickLabel', []);

    % Delay reconstructions
    for k = 1:nD
        d = d_avail(k);
        p = mcr.predictions(d);
        t_delay = mcr.t_pred(p.t_indices);
        nexttile; hold on;
        plot(t_delay, p.y_true, 'k-', 'LineWidth', 0.9);
        plot(t_delay, p.y_pred, '-', 'Color', colors(i, :), 'LineWidth', 1.4);
        hold off; grid on;
        ylabel(sprintf('u(t-%d)', d));
        title(sprintf('Delay d = %d:  R^2 = %.3f', d, mcr.R2_d(d)), 'FontWeight', 'normal');
        if k < nD; set(gca, 'XTickLabel', []); end
        if k == 1
            legend({'target', 'readout'}, 'Location', 'northeast', 'FontSize', 8);
        end
    end
    xlabel('Time (s)');
    sgtitle(sprintf('Memory reconstruction -- %s', condition_names{i}));
end

fprintf('Rendered MC comparison + %d reconstruction figures.\n', n_cond);
