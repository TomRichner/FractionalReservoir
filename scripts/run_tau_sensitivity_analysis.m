% run_tau_sensitivity_analysis.m
% Sensitivity analysis for tau_a_E(end) and tau_b_E_rec parameters
%
% This script analyzes the effect of the maximum SFA time constant (tau_a_E_max)
% and STD recovery time constant (tau_b_E_rec) on SRNN dynamics, producing
% heatmap plots of LLE distribution vs. parameter values.
%
% Similar to: SRNN/fig_scripts/SRNN_paper/tau_a_E_edge_of_chaos/tau_a_sensitivity_analysis.m
%
% See also: SensitivityAnalysis, SRNNModel

clear;
clc;
close all;

%% Setup paths
setup_paths();

%% Analysis Configuration
n_levels = 25;      % Number of parameter values to test
n_reps = 50;        % Number of repetitions per value
note = 'tau_timescales';

% Parameter ranges (matching SRNN analysis)
param_config = struct();
param_config.tau_a_E_max = struct('range', [5, 60], 'type', 'vector_max');
param_config.tau_b_E_rec = struct('range', [5, 60], 'type', 'scalar');

param_names = fieldnames(param_config);

% Condition: SFA + STD (n_a_E=3, n_b_E=1)
condition = struct('name', 'sfa_and_std', 'n_a_E', 3, 'n_b_E', 1);

% Model defaults
model_defaults = struct();
model_defaults.T_range = [0, 20];       % Simulation time
model_defaults.fs = 200;                % Sampling frequency
model_defaults.lya_method = 'benettin'; % Lyapunov method

%% Create output directory
dt_str = lower(strrep(datestr(now, 'mmm_dd_yy_HH_MM_AM'), ':', '_'));
output_dir_base = fullfile(pwd, sprintf('sensitivity_%s_nLevs_%d_nReps_%d_%s', ...
    note, n_levels, n_reps, dt_str));

if ~exist(output_dir_base, 'dir')
    mkdir(output_dir_base);
end

% Copy this script to output directory for reproducibility
copyfile([mfilename('fullpath') '.m'], output_dir_base);

% Create condition subdirectory
condition_dir = fullfile(output_dir_base, condition.name);
if ~exist(condition_dir, 'dir')
    mkdir(condition_dir);
end

%% Print analysis summary
fprintf('\n========================================\n');
fprintf('=== Tau Sensitivity Analysis ===\n');
fprintf('========================================\n');
fprintf('Parameters: %s\n', strjoin(param_names, ', '));
fprintf('Condition: %s (n_a_E=%d, n_b_E=%d)\n', condition.name, condition.n_a_E, condition.n_b_E);
fprintf('Levels per parameter: %d\n', n_levels);
fprintf('Repetitions per level: %d\n', n_reps);
fprintf('Total runs per parameter: %d\n', n_levels * n_reps);
fprintf('Output directory: %s\n', output_dir_base);
fprintf('========================================\n\n');

overall_start = tic;
all_results_summary = struct();

%% Run sensitivity analysis for each parameter
for p_idx = 1:length(param_names)
    param_name = param_names{p_idx};
    config = param_config.(param_name);
    param_range = config.range;
    param_type = config.type;
    
    fprintf('\n--- Processing: %s [%.1f, %.1f] (%d/%d) ---\n', ...
        param_name, param_range(1), param_range(2), p_idx, length(param_names));
    
    param_start = tic;
    
    % Create parameter levels
    param_levels = linspace(param_range(1), param_range(2), n_levels);
    
    total_runs = n_levels * n_reps;
    results = cell(total_runs, 1);
    
    % Extract condition values for parfor
    n_a_E_val = condition.n_a_E;
    n_b_E_val = condition.n_b_E;
    
    parfor k = 1:total_runs
        level_idx = floor((k-1) / n_reps) + 1;
        rep_idx = mod(k-1, n_reps) + 1;
        
        param_value = param_levels(level_idx);
        sim_seed = k + (p_idx-1)*total_runs*100;
        
        % Progress reporting (sparse)
        if mod(k, max(1, floor(total_runs/10))) == 1 || k == total_runs
            fprintf('  Run %d/%d: %s=%.1f, rep=%d, seed=%d\n', ...
                k, total_runs, param_name, param_value, rep_idx, sim_seed);
        end
        
        run_start = tic;
        try
            % Build model arguments
            model_args = {'n_a_E', n_a_E_val, 'n_b_E', n_b_E_val, ...
                          'rng_seeds', [sim_seed, sim_seed+1]};
            
            % Add model defaults
            default_fields = fieldnames(model_defaults);
            for d_idx = 1:length(default_fields)
                field_name = default_fields{d_idx};
                model_args = [model_args, {field_name, model_defaults.(field_name)}];
            end
            
            % Handle parameter type
            if strcmp(param_type, 'vector_max')
                % For tau_a_E_max: generate logspaced vector from 0.25 to max value
                tau_a_E_vec = logspace(log10(0.25), log10(param_value), n_a_E_val);
                model_args = [model_args, {'tau_a_E', tau_a_E_vec}];
            else
                % For scalar parameters like tau_b_E_rec
                model_args = [model_args, {param_name, param_value}];
            end
            
            % Create and run model
            model = SRNNModel(model_args{:});
            model.build();
            model.run();
            
            % Extract results
            result = struct();
            result.success = true;
            result.param_value = param_value;
            result.seed = sim_seed;
            result.run_duration = toc(run_start);
            
            % Extract LLE
            if ~isempty(model.lya_results) && isfield(model.lya_results, 'LLE')
                result.LLE = model.lya_results.LLE;
            else
                result.LLE = NaN;
            end
            
            % Extract mean firing rate
            if ~isempty(model.plot_data) && isfield(model.plot_data, 'r')
                r_E = model.plot_data.r.E;
                r_I = model.plot_data.r.I;
                all_rates = [r_E(:); r_I(:)];
                result.mean_rate = mean(all_rates(~isnan(all_rates)));
            else
                result.mean_rate = NaN;
            end
            
            results{k} = result;
            
        catch ME
            run_duration = toc(run_start);
            
            results{k} = struct(...
                'success', false, ...
                'error_message', ME.message, ...
                'param_value', param_value, ...
                'seed', sim_seed, ...
                'run_duration', run_duration, ...
                'LLE', NaN, ...
                'mean_rate', NaN ...
            );
            
            fprintf('  ERROR in run %d: %s\n', k, ME.message);
        end
    end
    
    % Reshape results to n_levels x n_reps
    results_reshaped = reshape(results, [n_reps, n_levels])';
    
    % Calculate statistics
    success_mask = cellfun(@(x) isfield(x, 'success') && x.success, results);
    n_success = sum(success_mask);
    n_failed = total_runs - n_success;
    
    param_elapsed = toc(param_start);
    fprintf('Parameter %s completed in %.2f minutes.\n', param_name, param_elapsed/60);
    fprintf('Success rate: %d/%d (%.1f%%)\n', n_success, total_runs, 100*n_success/total_runs);
    
    % Create metadata
    metadata = struct();
    metadata.param_name = param_name;
    metadata.param_levels = param_levels;
    metadata.param_range = param_range;
    metadata.param_type = param_type;
    metadata.n_levels = n_levels;
    metadata.n_reps = n_reps;
    metadata.n_success = n_success;
    metadata.n_failed = n_failed;
    metadata.success_rate = n_success / total_runs;
    metadata.elapsed_time = param_elapsed;
    metadata.condition = condition;
    metadata.model_defaults = model_defaults;
    metadata.analysis_date = datestr(now);
    
    % Save results
    save_filename = fullfile(condition_dir, sprintf('sensitivity_%s.mat', param_name));
    save(save_filename, 'results_reshaped', 'metadata', '-v7.3');
    fprintf('Saved to: %s\n', save_filename);
    
    all_results_summary.(param_name) = metadata;
end

%% Save summary
overall_elapsed = toc(overall_start);
fprintf('\n========================================\n');
fprintf('=== Analysis Complete ===\n');
fprintf('Total time: %.2f hours\n', overall_elapsed/3600);
fprintf('========================================\n');

summary_data = struct();
summary_data.all_conditions_summary = struct();
summary_data.all_conditions_summary.(condition.name) = all_results_summary;
summary_data.conditions = {condition};
summary_data.param_config = param_config;
summary_data.n_levels = n_levels;
summary_data.n_reps = n_reps;
summary_data.analysis_parameters = struct('n_levels', n_levels, 'n_reps', n_reps);
summary_data.model_defaults = model_defaults;
summary_data.total_elapsed_time = overall_elapsed;
summary_data.analysis_completed = datestr(now);

summary_file = fullfile(output_dir_base, 'sensitivity_analysis_summary_all_conditions.mat');
save(summary_file, 'summary_data', '-v7.3');
fprintf('Summary saved to: %s\n', summary_file);

%% Generate heatmap plots
fprintf('\nGenerating plots...\n');

% LLE histogram bins
lle_range = [-0.3, 0.1];
n_bins = 35;
lle_bins = [-inf, linspace(lle_range(1), lle_range(2), n_bins), inf];

% Create figure
n_params = length(param_names);
fig_lle = figure('Name', 'Tau Sensitivity - LLE', ...
    'Position', [100, 100, 400 * n_params, 400]);

for p_idx = 1:n_params
    param_name = param_names{p_idx};
    
    % Load results
    load_file = fullfile(condition_dir, sprintf('sensitivity_%s.mat', param_name));
    data = load(load_file);
    results_reshaped = data.results_reshaped;
    metadata = data.metadata;
    
    ax = subplot(1, n_params, p_idx);
    plot_sensitivity_heatmap(ax, results_reshaped, metadata, lle_bins, 'LLE', ...
        '$\lambda_1$', param_name, p_idx == 1);
end

% Save figure
fig_dir = fullfile(output_dir_base, 'figs');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end
saveas(fig_lle, fullfile(fig_dir, 'sensitivity_LLE_tau_params.png'));
saveas(fig_lle, fullfile(fig_dir, 'sensitivity_LLE_tau_params.fig'));
fprintf('Figures saved to: %s\n', fig_dir);

fprintf('\nDone! Results saved to: %s\n', output_dir_base);
beep; pause(0.5); beep; pause(0.2); beep  % Notification sound

%% Helper function for heatmap plotting
function plot_sensitivity_heatmap(ax, results_reshaped, metadata, hist_bins, metric, y_label, param_name, show_ylabel)
    n_levels = metadata.n_levels;
    n_reps = metadata.n_reps;
    param_levels = metadata.param_levels;
    
    num_bins = length(hist_bins) - 1;
    histogram_matrix = zeros(num_bins, n_levels);
    median_values = NaN(n_levels, 1);
    
    % Build histogram for each level
    for level_idx = 1:n_levels
        values_level = [];
        
        for rep_idx = 1:n_reps
            result = results_reshaped{level_idx, rep_idx};
            
            if isfield(result, 'success') && result.success && isfield(result, metric)
                val = result.(metric);
                if ~isnan(val)
                    values_level(end+1) = val;
                end
            end
        end
        
        if ~isempty(values_level)
            [counts, ~] = histcounts(values_level, hist_bins);
            histogram_matrix(:, level_idx) = counts;
            median_values(level_idx) = median(values_level);
        end
    end
    
    % Compute y-coordinates
    finite_edges = hist_bins(~isinf(hist_bins));
    step_size = finite_edges(2) - finite_edges(1);
    y_coords = zeros(num_bins, 1);
    y_coords(1) = finite_edges(1) - step_size/2;
    for k = 2:length(finite_edges)
        y_coords(k) = (finite_edges(k-1) + finite_edges(k)) / 2;
    end
    y_coords(end) = finite_edges(end) + step_size/2;
    
    % Plot
    imagesc(ax, param_levels, y_coords, histogram_matrix);
    hold(ax, 'on');
    yline(ax, 0, '--', 'Color', [0 0.7 0], 'LineWidth', 4, 'Alpha', 0.5);
    plot(ax, param_levels, median_values, 'b-', 'LineWidth', 4, 'Color', [0 0 1 0.55]);
    hold(ax, 'off');
    
    colormap(ax, flipud(gray));
    caxis(ax, [0, n_reps]);
    axis(ax, 'xy');
    box(ax, 'on');
    
    % Labels
    xlabel(ax, strrep(param_name, '_', '\_'), 'FontSize', 14);
    if show_ylabel
        ylabel(ax, y_label, 'Interpreter', 'latex', 'FontSize', 18);
    end
    
    title(ax, sprintf('%s', strrep(param_name, '_', ' ')), 'FontSize', 12);
    
    % Set y-ticks
    yticks(ax, [-0.3, -0.15, 0, 0.1]);
end
