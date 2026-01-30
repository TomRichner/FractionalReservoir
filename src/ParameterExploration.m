classdef ParameterExploration < handle
    % PARAMETEREXPLORATION Latin Hypercube Sampling-based parameter exploration
    %
    % This class performs comprehensive parameter space exploration using
    % Latin Hypercube Sampling (LHS) to efficiently cover the widest parameter
    % space. It detects unstable and dead simulations, supports GPU acceleration,
    % and provides optional result saving with full scientific rigor.
    %
    % Key features:
    %   - LHS sampling across all parameters
    %   - Automatic detection of unstable/dead simulations
    %   - Hybrid GPU/CPU computing support
    %   - Checkpoint/resume capability
    %   - Optional result saving
    %   - Comprehensive metrics extraction
    %
    % Usage:
    %   config = getDefaultExplorationConfig();
    %   explorer = ParameterExploration(config);
    %   results = explorer.runExploration();
    %
    % See also: SRNN_ESN, checkReservoirStability, extractReservoirMetrics
    
    %% Configuration Properties
    properties
        n_samples = 1000              % Number of LHS samples
        param_ranges = struct()      % Parameter ranges [min, max] for each param
        use_gpu = true                % Use GPU if available (hybrid mode)
        save_results = true           % Save results to disk
        output_dir = ''               % Output directory (empty = auto-generate)
        checkpoint_interval = 100     % Save checkpoint every N samples
        batch_size = 50               % Batch size for parfor
        verbose = true                % Print progress
        random_seed = 42              % Random seed for reproducibility
        note = ''                     % Optional note for output folder
        
        % Learning evaluation settings
        evaluate_learning = true      % Enable learning capability evaluation
        learning_task = 'prediction'  % Task: 'prediction', 'narma10', 'memory'
        learning_n_train = 2000       % Number of training samples for learning
        learning_n_test = 500         % Number of test samples for learning
        learning_washout = 100        % Washout steps for learning evaluation
    end
    
    %% Fixed Parameters (not explored)
    properties
        fixed_params = struct()       % Parameters that remain constant
    end
    
    %% Results (SetAccess = private)
    properties (SetAccess = private)
        results = struct()            % Stored results
        has_run = false               % Flag indicating if exploration has run
        analysis_start_time          % Timestamp when exploration started
        lhs_samples                  % Generated LHS samples
        param_names                  % List of parameter names being explored
        integer_params = {'n', 'n_a_E', 'n_a_I', 'n_b_E', 'n_b_I'}  % Integer parameters
        categorical_params = {'input_type'}  % Categorical parameters
    end
    
    %% Constructor
    methods
        function obj = ParameterExploration(config)
            % PARAMETEREXPLORATION Constructor
            %
            % Input:
            %   config - struct with configuration options (see properties)
            
            if nargin < 1
                config = struct();
            end
            
            % Set properties from config
            props = properties(obj);
            for i = 1:length(props)
                prop_name = props{i};
                if isfield(config, prop_name)
                    obj.(prop_name) = config.(prop_name);
                end
            end
            
            % Set default output directory if not specified
            if isempty(obj.output_dir) && obj.save_results
                obj.output_dir = pwd;
            end
        end
    end
    
    %% Public Methods
    methods
        function results = runExploration(obj)
            % RUNEXPLORATION Execute the parameter exploration
            %
            % Output:
            %   results - struct with all results and metrics
            
            % Validate configuration
            if isempty(fieldnames(obj.param_ranges))
                error('ParameterExploration:NoParameters', ...
                    'No parameter ranges defined. Set param_ranges first.');
            end
            
            % Initialize
            obj.analysis_start_time = datetime('now');
            obj.param_names = fieldnames(obj.param_ranges);
            
            % Generate LHS samples
            if obj.verbose
                fprintf('\n========================================\n');
                fprintf('=== Parameter Exploration (LHS) ===\n');
                fprintf('========================================\n');
                fprintf('Number of samples: %d\n', obj.n_samples);
                fprintf('Parameters: %s\n', strjoin(obj.param_names, ', '));
                fprintf('Batch size: %d\n', obj.batch_size);
                if obj.evaluate_learning
                    fprintf('Learning evaluation: ENABLED (%s task)\n', obj.learning_task);
                else
                    fprintf('Learning evaluation: DISABLED\n');
                end
            end
            
            obj.generateLHSSamples();
            
            % Create output directory if saving
            if obj.save_results
                obj.createOutputDirectory();
            end
            
            % Run simulations
            obj.runSimulations();
            
            % Aggregate results
            results = obj.aggregateResults();
            
            % Print final error summary if there were failures
            if obj.verbose
                obj.printErrorSummary();
            end
            
            % Save results if requested
            if obj.save_results
                obj.saveResults(results);
            end
            
            obj.has_run = true;
            
            if obj.verbose
                fprintf('\n========================================\n');
                fprintf('=== Exploration Complete ===\n');
                fprintf('========================================\n');
            end
        end
        
        function setParamRanges(obj, param_ranges)
            % SETPARAMRANGES Set parameter ranges for exploration
            %
            % Input:
            %   param_ranges - struct with fields for each parameter
            %                   Each field should be [min, max] or cell array
            %                   for categorical parameters
            
            obj.param_ranges = param_ranges;
            obj.param_names = fieldnames(param_ranges);
        end
        
        function setFixedParams(obj, fixed_params)
            % SETFIXEDPARAMS Set parameters that remain constant
            %
            % Input:
            %   fixed_params - struct with constant parameter values
            
            obj.fixed_params = fixed_params;
        end
        
        function fig_handles = plot(obj, results, varargin)
            % PLOT Generate visualization figures from exploration results
            %
            % Usage:
            %   explorer.plot(results)
            %   explorer.plot(results, 'save_figs', true)
            %   fig_handles = explorer.plot(results, 'metric', 'mean_firing_rate')
            %
            % Optional arguments:
            %   'metric' - Metric to focus on (default: 'mean_firing_rate')
            %   'save_figs' - Save figures to output directory (default: true)
            
            % Parse optional arguments
            save_figs = true;
            metric = 'mean_firing_rate';
            for i = 1:2:length(varargin)
                if strcmpi(varargin{i}, 'save_figs')
                    save_figs = varargin{i+1};
                elseif strcmpi(varargin{i}, 'metric')
                    metric = varargin{i+1};
                end
            end
            
            if obj.verbose
                fprintf('\nGenerating visualization figures...\n');
            end
            
            fig_handles = struct();
            
            % Figure 1: Summary Overview
            fig_handles.summary = obj.plotSummaryOverview(results);
            
            % Figure 2: Parameter Distributions
            fig_handles.distributions = obj.plotParameterDistributions(results);
            
            % Figure 3: Parameter vs Metric Scatter Plots
            fig_handles.scatter = obj.plotParameterVsMetric(results, metric);
            
            % Figure 4: Stability Analysis
            fig_handles.stability = obj.plotStabilityAnalysis(results);
            
            % Figure 5: Correlation Matrix
            fig_handles.correlation = obj.plotCorrelationMatrix(results);
            
            % Save figures if requested
            if save_figs && ~isempty(obj.output_dir)
                obj.saveFigures(fig_handles);
            end
            
            if obj.verbose
                fprintf('Visualization complete!\n');
            end
        end
    end
    
    %% Plotting Methods
    methods (Access = public)
        function fig = plotSummaryOverview(obj, results)
            % PLOTSUMMARYOVERVIEW Create summary overview figure
            
            fig = figure('Name', 'Parameter Exploration Summary', ...
                'Color', 'w');
            
            % Subplot 1: Success/Failure Pie Chart
            subplot(2, 3, 1);
            labels = {'Stable', 'Unstable', 'Dead', 'Failed'};
            n_samples = length(results.samples);
            n_stable = sum(cellfun(@(s) isstruct(s) && s.success && ~s.is_unstable && ~s.is_dead, results.samples));
            n_unstable = sum(cellfun(@(s) isstruct(s) && s.success && s.is_unstable, results.samples));
            n_dead = sum(cellfun(@(s) isstruct(s) && s.success && s.is_dead, results.samples));
            n_failed = n_samples - n_stable - n_unstable - n_dead;
            
            values = [n_stable, n_unstable, n_dead, n_failed];
            valid_idx = values > 0;
            if any(valid_idx)
                pie(values(valid_idx), labels(valid_idx));
                title('Simulation Outcomes', 'FontSize', 12, 'FontWeight', 'bold');
            end
            
            % Subplot 2: Summary Statistics Text
            subplot(2, 3, 2);
            axis off;
            text_str = {
                sprintf('Total Samples: %d', n_samples);
                sprintf('Stable: %d (%.1f%%)', n_stable, 100*n_stable/n_samples);
                sprintf('Unstable: %d (%.1f%%)', n_unstable, 100*n_unstable/n_samples);
                sprintf('Dead: %d (%.1f%%)', n_dead, 100*n_dead/n_samples);
                sprintf('Failed: %d (%.1f%%)', n_failed, 100*n_failed/n_samples);
                '';
                'Key Metrics:';
            };
            
            if isfield(results, 'metrics_stats') && isfield(results.metrics_stats, 'mean_firing_rate')
                mfr = results.metrics_stats.mean_firing_rate;
                text_str{end+1} = sprintf('Mean FR: %.3f +/- %.3f Hz', mfr.mean, mfr.std);
            end
            if isfield(results, 'metrics_stats') && isfield(results.metrics_stats, 'participation_ratio')
                pr = results.metrics_stats.participation_ratio;
                text_str{end+1} = sprintf('Part. Ratio: %.2f +/- %.2f', pr.mean, pr.std);
            end
            
            text(0.1, 0.9, text_str, 'Units', 'normalized', 'FontSize', 11, ...
                'VerticalAlignment', 'top', 'FontName', 'FixedWidth');
            title('Summary Statistics', 'FontSize', 12, 'FontWeight', 'bold');
            
            % Subplot 3: Mean Firing Rate Distribution
            subplot(2, 3, 3);
            firing_rates = obj.extractMetricFromResults(results, 'mean_firing_rate');
            if ~isempty(firing_rates)
                histogram(firing_rates, 30, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'none');
                xlabel('Mean Firing Rate (Hz)');
                ylabel('Count');
                title('Firing Rate Distribution', 'FontSize', 12, 'FontWeight', 'bold');
                grid on;
            end
            
            % Subplot 4: Participation Ratio Distribution
            subplot(2, 3, 4);
            part_ratios = obj.extractMetricFromResults(results, 'participation_ratio');
            if ~isempty(part_ratios)
                histogram(part_ratios, 30, 'FaceColor', [0.9 0.5 0.3], 'EdgeColor', 'none');
                xlabel('Participation Ratio');
                ylabel('Count');
                title('Dimensionality Distribution', 'FontSize', 12, 'FontWeight', 'bold');
                grid on;
            end
            
            % Subplot 5: E/I Balance Distribution
            subplot(2, 3, 5);
            ei_balance = obj.extractMetricFromResults(results, 'EI_balance');
            if ~isempty(ei_balance)
                valid_ei = ei_balance(isfinite(ei_balance) & ei_balance < 10);
                if ~isempty(valid_ei)
                    histogram(valid_ei, 30, 'FaceColor', [0.5 0.8 0.5], 'EdgeColor', 'none');
                    xlabel('E/I Balance Ratio');
                    ylabel('Count');
                    title('E/I Balance Distribution', 'FontSize', 12, 'FontWeight', 'bold');
                    grid on;
                end
            end
            
            % Subplot 6: Stability Reason Breakdown
            subplot(2, 3, 6);
            reasons = obj.extractReasonsFromResults(results);
            if ~isempty(reasons)
                [unique_reasons, ~, idx] = unique(reasons);
                counts = accumarray(idx, 1);
                barh(counts);
                set(gca, 'YTick', 1:length(unique_reasons), 'YTickLabel', strrep(unique_reasons, '_', ' '));
                xlabel('Count');
                title('Stability Reasons', 'FontSize', 12, 'FontWeight', 'bold');
                grid on;
            end
            
            sgtitle('Parameter Exploration Summary', 'FontSize', 14, 'FontWeight', 'bold');
        end
        
        function fig = plotParameterDistributions(obj, results)
            % PLOTPARAMETERDISTRIBUTIONS Plot distributions of sampled parameters
            
            n_params = length(obj.param_names);
            n_cols = ceil(sqrt(n_params));
            n_rows = ceil(n_params / n_cols);
            
            fig = figure('Name', 'Parameter Distributions', ...
                'Color', 'w');
            
            for i = 1:n_params
                subplot(n_rows, n_cols, i);
                
                param_name = obj.param_names{i};
                param_values = results.params(:, i);
                
                % Check if categorical
                if ismember(param_name, obj.categorical_params)
                    param_range = obj.param_ranges.(param_name);
                    if iscell(param_range)
                        histogram(param_values, length(param_range), 'FaceColor', [0.4 0.6 0.8]);
                        set(gca, 'XTick', 1:length(param_range), 'XTickLabel', param_range);
                        xtickangle(45);
                    end
                else
                    histogram(param_values, 20, 'FaceColor', [0.4 0.6 0.8], 'EdgeColor', 'none');
                end
                
                xlabel(strrep(param_name, '_', '\_'));
                ylabel('Count');
                title(strrep(param_name, '_', '\_'), 'FontSize', 10);
                grid on;
            end
            
            sgtitle('LHS Parameter Distributions', 'FontSize', 14, 'FontWeight', 'bold');
        end
        
        function fig = plotParameterVsMetric(obj, results, metric)
            % PLOTPARAMETERVMETRIC Plot parameter vs metric scatter plots
            
            metric_values = obj.extractMetricFromResults(results, metric);
            if isempty(metric_values)
                fig = figure('Visible', 'off');
                return;
            end
            
            % Get stability flags for coloring
            is_stable = cellfun(@(s) isstruct(s) && s.success && ~s.is_unstable && ~s.is_dead, results.samples);
            
            n_params = min(length(obj.param_names), 12);  % Limit to 12 params
            n_cols = ceil(sqrt(n_params));
            n_rows = ceil(n_params / n_cols);
            
            fig = figure('Name', sprintf('Parameters vs %s', metric), ...
                'Color', 'w');
            
            for i = 1:n_params
                subplot(n_rows, n_cols, i);
                
                param_name = obj.param_names{i};
                param_values = results.params(:, i);
                
                % Get valid indices (where we have both param and metric)
                valid_idx = ~isnan(metric_values) & isfinite(metric_values);
                
                if sum(valid_idx) > 0
                    scatter(param_values(valid_idx & is_stable), metric_values(valid_idx & is_stable), ...
                        20, [0.2 0.6 0.2], 'filled', 'MarkerFaceAlpha', 0.5);
                    hold on;
                    scatter(param_values(valid_idx & ~is_stable), metric_values(valid_idx & ~is_stable), ...
                        20, [0.8 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.5);
                    hold off;
                    
                    % Add trend line for stable samples
                    if sum(valid_idx & is_stable) > 5
                        hold on;
                        p = polyfit(param_values(valid_idx & is_stable), metric_values(valid_idx & is_stable), 1);
                        x_fit = linspace(min(param_values), max(param_values), 100);
                        y_fit = polyval(p, x_fit);
                        plot(x_fit, y_fit, 'k--', 'LineWidth', 1.5);
                        hold off;
                    end
                end
                
                xlabel(strrep(param_name, '_', '\_'));
                ylabel(strrep(metric, '_', '\_'));
                title(strrep(param_name, '_', '\_'), 'FontSize', 10);
                grid on;
            end
            
            % Add legend to first subplot
            subplot(n_rows, n_cols, 1);
            legend('Stable', 'Unstable/Dead', 'Location', 'best');
            
            sgtitle(sprintf('Parameters vs %s', strrep(metric, '_', ' ')), ...
                'FontSize', 14, 'FontWeight', 'bold');
        end
        
        function fig = plotStabilityAnalysis(obj, results)
            % PLOTSTABILITYANALYSIS Detailed stability analysis figure
            
            fig = figure('Name', 'Stability Analysis', ...
                'Color', 'w');
            
            % Extract stability data
            is_stable = cellfun(@(s) isstruct(s) && s.success && ~s.is_unstable && ~s.is_dead, results.samples);
            is_unstable = cellfun(@(s) isstruct(s) && s.success && s.is_unstable, results.samples);
            is_dead = cellfun(@(s) isstruct(s) && s.success && s.is_dead, results.samples);
            
            % Key parameters for stability
            key_params = {'level_of_chaos', 'fraction_E', 'c_E', 'input_scaling'};
            available_params = intersect(key_params, obj.param_names);
            
            n_plots = length(available_params);
            if n_plots == 0
                available_params = obj.param_names(1:min(4, length(obj.param_names)));
                n_plots = length(available_params);
            end
            
            for i = 1:n_plots
                subplot(2, n_plots, i);
                
                param_name = available_params{i};
                param_idx = find(strcmp(obj.param_names, param_name));
                param_values = results.params(:, param_idx);
                
                % Stability histogram
                edges = linspace(min(param_values), max(param_values), 15);
                
                [n_stable_hist, ~] = histcounts(param_values(is_stable), edges);
                [n_unstable_hist, ~] = histcounts(param_values(is_unstable), edges);
                [n_dead_hist, ~] = histcounts(param_values(is_dead), edges);
                
                bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
                
                bar(bin_centers, [n_stable_hist; n_unstable_hist; n_dead_hist]', 'stacked');
                colormap(gca, [0.3 0.7 0.3; 0.9 0.3 0.3; 0.5 0.5 0.5]);
                xlabel(strrep(param_name, '_', '\_'));
                ylabel('Count');
                title(sprintf('Stability vs %s', strrep(param_name, '_', '\_')), 'FontSize', 10);
                if i == 1
                    legend('Stable', 'Unstable', 'Dead', 'Location', 'best');
                end
                grid on;
            end
            
            % Bottom row: Stability rate vs parameter
            for i = 1:n_plots
                subplot(2, n_plots, n_plots + i);
                
                param_name = available_params{i};
                param_idx = find(strcmp(obj.param_names, param_name));
                param_values = results.params(:, param_idx);
                
                % Compute stability rate in bins
                edges = linspace(min(param_values), max(param_values), 10);
                bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
                
                stability_rate = zeros(length(bin_centers), 1);
                for b = 1:length(bin_centers)
                    in_bin = param_values >= edges(b) & param_values < edges(b+1);
                    if sum(in_bin) > 0
                        stability_rate(b) = sum(is_stable & in_bin) / sum(in_bin);
                    end
                end
                
                plot(bin_centers, stability_rate * 100, 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
                    'MarkerFaceColor', [0.3 0.6 0.9]);
                xlabel(strrep(param_name, '_', '\_'));
                ylabel('Stable Rate (%)');
                ylim([0 100]);
                title(sprintf('Stability Rate vs %s', strrep(param_name, '_', '\_')), 'FontSize', 10);
                grid on;
            end
            
            sgtitle('Stability Analysis', 'FontSize', 14, 'FontWeight', 'bold');
        end
        
        function fig = plotCorrelationMatrix(obj, results)
            % PLOTCORRELATIONMATRIX Plot correlation matrix between parameters and metrics
            
            fig = figure('Name', 'Correlation Matrix', ...
                'Color', 'w');
            
            % Extract metrics for stable samples
            is_stable = cellfun(@(s) isstruct(s) && s.success && ~s.is_unstable && ~s.is_dead, results.samples);
            
            metric_names = {'mean_firing_rate', 'participation_ratio', 'EI_balance', ...
                'silence_mean_rate', 'coefficient_of_variation'};
            
            % Build data matrix
            n_params = length(obj.param_names);
            n_metrics = length(metric_names);
            
            param_data = results.params(is_stable, :);
            metric_data = zeros(sum(is_stable), n_metrics);
            
            for i = 1:n_metrics
                metric_vals = obj.extractMetricFromResults(results, metric_names{i});
                if ~isempty(metric_vals)
                    metric_data(:, i) = metric_vals(is_stable);
                else
                    metric_data(:, i) = NaN;
                end
            end
            
            % Combine and compute correlation
            all_data = [param_data, metric_data];
            all_names = [obj.param_names; metric_names'];
            
            % Remove columns with NaN or Inf
            valid_cols = all(isfinite(all_data), 1);
            all_data = all_data(:, valid_cols);
            all_names = all_names(valid_cols);
            
            if size(all_data, 2) > 2
                R = corrcoef(all_data, 'Rows', 'complete');
                
                imagesc(R);
                colormap(bluewhitered_colormap(256));
                colorbar;
                caxis([-1 1]);
                
                set(gca, 'XTick', 1:length(all_names), 'XTickLabel', strrep(all_names, '_', ' '));
                set(gca, 'YTick', 1:length(all_names), 'YTickLabel', strrep(all_names, '_', ' '));
                xtickangle(45);
                
                title('Parameter-Metric Correlation Matrix', 'FontSize', 14, 'FontWeight', 'bold');
                axis square;
            else
                text(0.5, 0.5, 'Insufficient data for correlation', ...
                    'HorizontalAlignment', 'center', 'Units', 'normalized');
            end
        end
    end
    
    %% Helper methods for plotting
    methods (Access = private)
        function values = extractMetricFromResults(obj, results, metric_name)
            % Extract a specific metric from all results
            n_samples = length(results.samples);
            values = NaN(n_samples, 1);
            
            for i = 1:n_samples
                sample = results.samples{i};
                if isstruct(sample) && isfield(sample, 'metrics') && ...
                        isfield(sample.metrics, metric_name)
                    val = sample.metrics.(metric_name);
                    if isscalar(val) && isfinite(val)
                        values(i) = val;
                    end
                end
            end
        end
        
        function reasons = extractReasonsFromResults(obj, results)
            % Extract stability reasons from all results
            n_samples = length(results.samples);
            reasons = cell(n_samples, 1);
            
            for i = 1:n_samples
                sample = results.samples{i};
                if isstruct(sample) && isfield(sample, 'reason')
                    reasons{i} = sample.reason;
                else
                    reasons{i} = 'unknown';
                end
            end
        end
        
        function saveFigures(obj, fig_handles)
            % Save all figures to output directory
            
            fig_dir = fullfile(obj.output_dir, 'figures');
            if ~exist(fig_dir, 'dir')
                mkdir(fig_dir);
            end
            
            fig_names = fieldnames(fig_handles);
            for i = 1:length(fig_names)
                fig_name = fig_names{i};
                fig = fig_handles.(fig_name);
                
                if isvalid(fig)
                    saveas(fig, fullfile(fig_dir, sprintf('%s.png', fig_name)));
                    saveas(fig, fullfile(fig_dir, sprintf('%s.fig', fig_name)));
                end
            end
            
            if obj.verbose
                fprintf('Figures saved to: %s\n', fig_dir);
            end
        end
    end
    
    %% Private Methods
    methods (Access = private)
        function generateLHSSamples(obj)
            % GENERATELHSSAMPLES Generate Latin Hypercube samples
            
            n_params = length(obj.param_names);
            
            % Set random seed for reproducibility
            rng(obj.random_seed);
            
            % Generate LHS design (uniform in [0, 1])
            lhs_design = lhsdesign(obj.n_samples, n_params);
            
            % Scale to parameter ranges
            obj.lhs_samples = zeros(obj.n_samples, n_params);
            
            for i = 1:n_params
                param_name = obj.param_names{i};
                param_range = obj.param_ranges.(param_name);
                
                if iscell(param_range)
                    % Categorical parameter - map to indices
                    n_categories = length(param_range);
                    indices = round(lhs_design(:, i) * (n_categories - 1)) + 1;
                    indices = max(1, min(n_categories, indices));  % Clamp
                    obj.lhs_samples(:, i) = indices;
                else
                    % Numeric parameter
                    if ismember(param_name, obj.integer_params)
                        % Integer parameter - round after scaling
                        scaled = lhs_design(:, i) * (param_range(2) - param_range(1)) + param_range(1);
                        obj.lhs_samples(:, i) = round(scaled);
                    else
                        % Continuous parameter
                        obj.lhs_samples(:, i) = lhs_design(:, i) * (param_range(2) - param_range(1)) + param_range(1);
                    end
                end
            end
            
            if obj.verbose
                fprintf('Generated %d LHS samples\n', obj.n_samples);
            end
        end
        
        function createOutputDirectory(obj)
            % CREATEOUTPUTDIRECTORY Create timestamped output directory
            
            dt_str = lower(strrep(datestr(obj.analysis_start_time, 'mmm_dd_yy_HH_MM_AM'), ':', '_'));
            
            if ~isempty(obj.note)
                folder_name = sprintf('param_exploration_%s_nSamples_%d_%s', ...
                    obj.note, obj.n_samples, dt_str);
            else
                folder_name = sprintf('param_exploration_nSamples_%d_%s', ...
                    obj.n_samples, dt_str);
            end
            
            obj.output_dir = fullfile(obj.output_dir, folder_name);
            if ~exist(obj.output_dir, 'dir')
                mkdir(obj.output_dir);
            end
            
            if obj.verbose
                fprintf('Output directory: %s\n', obj.output_dir);
            end
        end
        
        function runSimulations(obj)
            % RUNSIMULATIONS Execute simulations in batches
            
            n_batches = ceil(obj.n_samples / obj.batch_size);
            
            % Initialize results storage
            obj.results.samples = cell(obj.n_samples, 1);
            obj.results.params = obj.lhs_samples;
            obj.results.param_names = obj.param_names;
            
            % Check for existing checkpoint
            checkpoint_file = fullfile(obj.output_dir, 'checkpoint.mat');
            start_batch = 1;
            
            if exist(checkpoint_file, 'file') && obj.save_results
                if obj.verbose
                    fprintf('Loading checkpoint...\n');
                end
                loaded = load(checkpoint_file);
                if isfield(loaded, 'results')
                    obj.results = loaded.results;
                    start_batch = loaded.last_batch + 1;
                    if obj.verbose
                        fprintf('Resuming from batch %d\n', start_batch);
                    end
                end
            end
            
            % Check GPU availability
            use_gpu_flag = obj.use_gpu && gpuReservoirOps('check');
            if obj.verbose
                if use_gpu_flag
                    fprintf('GPU computing enabled (hybrid mode)\n');
                else
                    fprintf('CPU computing only\n');
                end
            end
            
            % Extract data needed for parfor (avoid accessing obj inside parfor)
            all_lhs_samples = obj.lhs_samples;
            param_names_local = obj.param_names;
            param_ranges_local = obj.param_ranges;
            fixed_params_local = obj.fixed_params;
            categorical_params_local = obj.categorical_params;
            integer_params_local = obj.integer_params;
            verbose_local = obj.verbose;
            
            % Learning evaluation settings
            learning_settings = struct();
            learning_settings.evaluate = obj.evaluate_learning;
            learning_settings.task = obj.learning_task;
            learning_settings.n_train = obj.learning_n_train;
            learning_settings.n_test = obj.learning_n_test;
            learning_settings.washout = obj.learning_washout;
            
            % Run batches
            for batch_idx = start_batch:n_batches
                batch_start = (batch_idx - 1) * obj.batch_size + 1;
                batch_end = min(batch_idx * obj.batch_size, obj.n_samples);
                batch_indices = batch_start:batch_end;
                current_batch_size = length(batch_indices);
                
                if verbose_local
                    fprintf('\n--- Batch %d/%d (samples %d-%d) ---\n', ...
                        batch_idx, n_batches, batch_start, batch_end);
                end
                
                % Extract batch data for parfor
                batch_param_values = all_lhs_samples(batch_indices, :);
                
                % Run batch in parallel
                batch_results = cell(current_batch_size, 1);
                
                parfor i = 1:current_batch_size
                    param_values = batch_param_values(i, :);
                    
                    try
                        result = runSingleSimulationStatic(param_values, use_gpu_flag, ...
                            param_names_local, param_ranges_local, fixed_params_local, ...
                            categorical_params_local, integer_params_local, learning_settings);
                        batch_results{i} = result;
                    catch ME
                        % Create detailed error result
                        result = struct();
                        result.success = false;
                        result.error_message = ME.message;
                        result.error_identifier = ME.identifier;
                        result.error_stack = ME.stack;
                        result.param_values = param_values;
                        result.is_unstable = true;
                        result.is_dead = false;
                        result.reason = 'simulation_error';
                        result.metrics = struct();  % Empty metrics
                        result.stability = struct('is_unstable', true, 'is_dead', false, 'reason', 'simulation_error', 'metrics', struct());
                        batch_results{i} = result;
                    end
                end
                
                % Store batch results
                for i = 1:current_batch_size
                    obj.results.samples{batch_indices(i)} = batch_results{i};
                end
                
                % Save checkpoint
                if obj.save_results && mod(batch_idx, ceil(obj.checkpoint_interval / obj.batch_size)) == 0
                    obj.saveCheckpoint(batch_idx);
                end
                
                % Report results and errors
                if verbose_local
                    n_success = sum(cellfun(@(r) isstruct(r) && isfield(r, 'success') && r.success, batch_results));
                    n_failed = current_batch_size - n_success;
                    fprintf('Batch %d complete: %d/%d successful, %d failed\n', batch_idx, n_success, current_batch_size, n_failed);
                    
                    % Print error details for failed simulations
                    if n_failed > 0
                        fprintf('\n  === ERROR DETAILS ===\n');
                        for i = 1:current_batch_size
                            r = batch_results{i};
                            if isstruct(r) && isfield(r, 'success') && ~r.success
                                sample_idx = batch_indices(i);
                                fprintf('  Sample %d: %s\n', sample_idx, r.error_message);
                                if isfield(r, 'error_identifier') && ~isempty(r.error_identifier)
                                    fprintf('    Identifier: %s\n', r.error_identifier);
                                end
                                if isfield(r, 'error_stack') && ~isempty(r.error_stack)
                                    fprintf('    Stack trace:\n');
                                    for j = 1:min(5, length(r.error_stack))  % Show first 5 stack frames
                                        fprintf('      %s (line %d)\n', r.error_stack(j).name, r.error_stack(j).line);
                                    end
                                end
                                if isfield(r, 'param_values') && ~isempty(r.param_values)
                                    fprintf('    Parameters: ');
                                    for k = 1:min(5, length(param_names_local))
                                        fprintf('%s=%.3g ', param_names_local{k}, r.param_values(k));
                                    end
                                    fprintf('...\n');
                                end
                                fprintf('\n');
                            end
                        end
                        fprintf('  === END ERROR DETAILS ===\n\n');
                    end
                end
            end
            
            % Final checkpoint
            if obj.save_results
                obj.saveCheckpoint(n_batches);
            end
        end
        
        function result = runSingleSimulation(obj, param_values, use_gpu)
            % RUNSINGLESIMULATION Run a single simulation with given parameters
            %
            % Input:
            %   param_values - vector of parameter values (same order as param_names)
            %   use_gpu - whether to use GPU
            %
            % Output:
            %   result - struct with simulation results
            
            % Build parameter struct
            params = obj.buildParameterStruct(param_values);
            
            % Run simulation
            [states_r, S_history, metrics, stability] = obj.executeSimulation(params, use_gpu);
            
            % Package result
            result = struct();
            result.success = true;
            result.params = params;
            result.states_r = states_r;
            result.S_history = S_history;
            result.metrics = metrics;
            result.stability = stability;
            result.is_unstable = stability.is_unstable;
            result.is_dead = stability.is_dead;
            result.reason = stability.reason;
        end
        
        function params = buildParameterStruct(obj, param_values)
            % BUILDPARAMETERSTRUCT Build parameter struct from LHS sample values
            
            params = obj.fixed_params;  % Start with fixed parameters
            
            % Add sampled parameters
            for i = 1:length(obj.param_names)
                param_name = obj.param_names{i};
                param_value = param_values(i);
                
                % Handle categorical parameters
                if ismember(param_name, obj.categorical_params)
                    param_range = obj.param_ranges.(param_name);
                    if iscell(param_range)
                        idx = round(param_value);
                        idx = max(1, min(length(param_range), idx));
                        params.(param_name) = param_range{idx};
                    else
                        params.(param_name) = param_value;
                    end
                else
                    params.(param_name) = param_value;
                end
            end
            
            % Set defaults for required parameters if not specified
            params = obj.setDefaultParameters(params);
        end
        
        function params = setDefaultParameters(obj, params)
            % SETDEFAULTPARAMETERS Set default values for required parameters
            
            defaults = struct();
            defaults.dt = 0.1;
            defaults.n_steps = 300;
            defaults.N = 10;  % Silence before stimulus (seconds)
            defaults.M = 200;  % Silence after stimulus (seconds)
            defaults.input_amplitude = 1;
            defaults.S_a = 0.85;
            defaults.S_c = 0.4;
            defaults.pulse_period = 20;
            defaults.pulse_width = 5;
            defaults.input_frequency = 0.001;
            
            default_fields = fieldnames(defaults);
            for i = 1:length(default_fields)
                field_name = default_fields{i};
                if ~isfield(params, field_name)
                    params.(field_name) = defaults.(field_name);
                end
            end
            
            % Ensure required fields exist
            if ~isfield(params, 'n') || ~isfield(params, 'fraction_E')
                error('ParameterExploration:MissingRequired', ...
                    'n and fraction_E are required parameters');
            end
            
            % Compute derived parameters
            params.n_E = round(params.n * params.fraction_E);
            params.n_I = params.n - params.n_E;
            params.E_indices = 1:params.n_E;
            params.I_indices = (params.n_E+1):params.n;
            
            % Set activation function
            if ~isfield(params, 'activation_function')
                params.activation_function = @(x) piecewiseSigmoid(x, params.S_a, params.S_c);
                params.activation_function_derivative = @(x) piecewiseSigmoidDerivative(x, params.S_a, params.S_c);
            end
            
            % Set adaptation timescales if not specified
            if isfield(params, 'n_a_E') && params.n_a_E > 0 && ~isfield(params, 'tau_a_E')
                params.tau_a_E = logspace(log10(0.25), log10(25), params.n_a_E);
            end
            if isfield(params, 'n_a_I') && params.n_a_I > 0 && ~isfield(params, 'tau_a_I')
                params.tau_a_I = logspace(log10(0.25), log10(25), params.n_a_I);
            end
        end
        
        function [states_r, S_history, metrics, stability] = executeSimulation(obj, params, use_gpu)
            % EXECUTESIMULATION Execute a single reservoir simulation
            %
            % This method extracts the simulation logic from test_reservoir_dynamics.m
            
            % Generate input signal
            [U, t, silence_start_idx] = obj.generateInputSignal(params);
            
            % Build reservoir
            [W, W_in] = obj.buildReservoir(params);
            
            % Create ESN parameters
            esn_params = obj.buildESNParams(params, W, W_in);
            
            % Create and run ESN
            esn = SRNN_ESN(esn_params);
            esn.resetState();
            
            % Run reservoir
            [~, S_history] = esn.runReservoir(U);
            
            % Extract firing rates
            states_r = obj.extractFiringRates(S_history, params);
            
            % Check stability
            stability = checkReservoirStability(states_r, silence_start_idx);
            
            % Extract metrics
            metrics = extractReservoirMetrics(states_r, S_history, params, silence_start_idx);
        end
        
        function [U, t, silence_start_idx] = generateInputSignal(obj, params)
            % GENERATEINPUTSIGNAL Generate input signal based on parameters
            
            dt = params.dt;
            n_steps = params.n_steps;
            N = params.N;  % Silence before
            M = params.M;  % Silence after
            
            n_silence_steps_start = round(N / dt);
            n_silence_steps_end = round(M / dt);
            total_steps = n_steps + n_silence_steps_start + n_silence_steps_end;
            
            t = (0:n_steps-1)' * dt;
            U = zeros(total_steps, 1);
            
            % Generate input based on type
            input_type = getFieldOrDefault(params, 'input_type', 'pulse');
            input_amplitude = getFieldOrDefault(params, 'input_amplitude', 1);
            
            switch lower(input_type)
                case 'sine'
                    input_frequency = getFieldOrDefault(params, 'input_frequency', 0.001);
                    U(n_silence_steps_start+1:n_silence_steps_start+n_steps) = ...
                        input_amplitude * sin(2*pi*input_frequency*t);
                    
                case 'pulse'
                    pulse_period = getFieldOrDefault(params, 'pulse_period', 20);
                    pulse_width = getFieldOrDefault(params, 'pulse_width', 5);
                    pulse_signal = zeros(n_steps, 1);
                    pulse_signal(mod(0:n_steps-1, pulse_period) < pulse_width) = input_amplitude;
                    U(n_silence_steps_start+1:n_silence_steps_start+n_steps) = pulse_signal;
                    
                case 'step'
                    step_signal = zeros(n_steps, 1);
                    step_signal(round(n_steps/4):end) = input_amplitude;
                    U(n_silence_steps_start+1:n_silence_steps_start+n_steps) = step_signal;
                    
                case 'noise'
                    rng('shuffle');  % Different seed for each simulation
                    U(n_silence_steps_start+1:n_silence_steps_start+n_steps) = ...
                        input_amplitude * randn(n_steps, 1);
                    
                case 'mackey'
                    [~, mg] = generate_mackey_glass('tau', 17, 'n_samples', n_steps+1000, ...
                        'dt', dt, 'discard', 1000);
                    U(n_silence_steps_start+1:n_silence_steps_start+n_steps) = ...
                        input_amplitude * mg(1:n_steps);
                    
                otherwise
                    error('Unknown input_type: %s', input_type);
            end
            
            U = U(:);
            t = (0:length(U)-1)' * dt;
            silence_start_idx = n_silence_steps_start + n_steps + 1;
        end
        
        function [W, W_in] = buildReservoir(obj, params)
            % BUILDRESERVOIR Build reservoir weight matrices
            
            n = params.n;
            n_E = params.n_E;
            n_I = params.n_I;
            level_of_chaos = params.level_of_chaos;
            input_scaling = getFieldOrDefault(params, 'input_scaling', 0.75);
            
            % Create recurrent weight matrix
            rng('shuffle');  % Different seed for each simulation
            W = randn(n, n);
            
            % Apply Dale's law
            W(:, 1:n_E) = abs(W(:, 1:n_E));           % E columns positive
            W(:, n_E+1:end) = -abs(W(:, n_E+1:end));  % I columns negative
            
            % Center rows
            W = W - mean(W, 2);
            
            % Apply chaos scaling
            W_eigs = eig(W);
            abscissa_0 = max(real(W_eigs));
            if abs(abscissa_0) > 1e-10
                gamma = 1 / abscissa_0;
                W = level_of_chaos * gamma * W;
            end
            
            % Create input weight matrix
            n_inputs = 1;  % Single input dimension
            W_in = (2 * rand(n, n_inputs) - 1) * input_scaling;
            W_in(rand(n, n_inputs) > 0.2) = 0;  % Sparse input
        end
        
        function esn_params = buildESNParams(obj, params, W, W_in)
            % BUILDESNPARAMS Build parameters struct for SRNN_ESN
            
            esn_params = struct();
            esn_params.n = params.n;
            esn_params.n_E = params.n_E;
            esn_params.n_I = params.n_I;
            esn_params.W = W;
            esn_params.W_in = W_in;
            esn_params.tau_d = params.tau_d;
            esn_params.activation_function = params.activation_function;
            esn_params.activation_function_derivative = params.activation_function_derivative;
            esn_params.E_indices = params.E_indices;
            esn_params.I_indices = params.I_indices;
            esn_params.which_states = 'x';
            esn_params.include_input = false;
            esn_params.lambda = 1e-6;
            esn_params.dt = params.dt;
            
            % Optional parameters
            if isfield(params, 'n_a_E')
                esn_params.n_a_E = params.n_a_E;
            end
            if isfield(params, 'n_a_I')
                esn_params.n_a_I = params.n_a_I;
            end
            if isfield(params, 'tau_a_E')
                esn_params.tau_a_E = params.tau_a_E;
            end
            if isfield(params, 'tau_a_I')
                esn_params.tau_a_I = params.tau_a_I;
            end
            if isfield(params, 'n_b_E')
                esn_params.n_b_E = params.n_b_E;
            end
            if isfield(params, 'n_b_I')
                esn_params.n_b_I = params.n_b_I;
            end
            if isfield(params, 'tau_b_E_rec')
                esn_params.tau_b_E_rec = params.tau_b_E_rec;
            end
            if isfield(params, 'tau_b_E_rel')
                esn_params.tau_b_E_rel = params.tau_b_E_rel;
            end
            if isfield(params, 'tau_b_I_rec')
                esn_params.tau_b_I_rec = params.tau_b_I_rec;
            end
            if isfield(params, 'tau_b_I_rel')
                esn_params.tau_b_I_rel = params.tau_b_I_rel;
            end
            if isfield(params, 'c_E')
                esn_params.c_E = params.c_E;
            end
            if isfield(params, 'c_I')
                esn_params.c_I = params.c_I;
            end
            if isfield(params, 'lags')
                esn_params.lags = params.lags;
            end
        end
        
        function states_r = extractFiringRates(obj, S_history, params)
            % EXTRACTFIRINGRATES Extract firing rates from state history
            
            n = params.n;
            n_E = params.n_E;
            n_I = params.n_I;
            n_a_E = getFieldOrDefault(params, 'n_a_E', 0);
            n_a_I = getFieldOrDefault(params, 'n_a_I', 0);
            n_b_E = getFieldOrDefault(params, 'n_b_E', 0);
            n_b_I = getFieldOrDefault(params, 'n_b_I', 0);
            c_E = getFieldOrDefault(params, 'c_E', 0);
            c_I = getFieldOrDefault(params, 'c_I', 0);
            activation_function = params.activation_function;
            
            total_steps = size(S_history, 1);
            states_r = zeros(total_steps, n);
            
            len_a_E = n_E * n_a_E;
            len_a_I = n_I * n_a_I;
            len_b_E = n_E * n_b_E;
            len_b_I = n_I * n_b_I;
            
            for i = 1:total_steps
                current_idx = 1;
                
                % Extract adaptation variables a_E
                if n_a_E > 0
                    a_E = reshape(S_history(i, current_idx:current_idx+len_a_E-1), n_a_E, n_E)';
                    current_idx = current_idx + len_a_E;
                else
                    a_E = zeros(n_E, 0);
                end
                
                % Extract adaptation variables a_I
                if n_a_I > 0
                    a_I = reshape(S_history(i, current_idx:current_idx+len_a_I-1), n_a_I, n_I)';
                    current_idx = current_idx + len_a_I;
                else
                    a_I = zeros(n_I, 0);
                end
                
                % Skip depression variables
                current_idx = current_idx + len_b_E + len_b_I;
                
                % Extract dendritic states x
                x = S_history(i, current_idx:current_idx+n-1)';
                
                % Compute effective input with adaptation
                x_eff = x;
                if n_a_E > 0 && n_E > 0
                    x_eff(1:n_E) = x_eff(1:n_E) - c_E * sum(a_E, 2);
                end
                if n_a_I > 0 && n_I > 0
                    x_eff(n_E+1:end) = x_eff(n_E+1:end) - c_I * sum(a_I, 2);
                end
                
                % Apply activation function
                states_r(i, :) = activation_function(x_eff);
            end
        end
        
        function saveCheckpoint(obj, last_batch)
            % SAVECHECKPOINT Save checkpoint for resume capability
            
            if ~obj.save_results
                return;
            end
            
            checkpoint_file = fullfile(obj.output_dir, 'checkpoint.mat');
            save(checkpoint_file, 'last_batch', '-v7.3');
            
            % Also save results structure
            results = obj.results;
            save(fullfile(obj.output_dir, 'results_checkpoint.mat'), 'results', '-v7.3');
        end
        
        function results = aggregateResults(obj)
            % AGGREGATERESULTS Aggregate and analyze all results
            
            results = struct();
            results.params = obj.results.params;
            results.param_names = obj.results.param_names;
            results.samples = obj.results.samples;
            
            % Extract metrics and stability flags
            n_samples = length(obj.results.samples);
            success_flags = false(n_samples, 1);
            unstable_flags = false(n_samples, 1);
            dead_flags = false(n_samples, 1);
            learning_evaluated = false(n_samples, 1);
            learning_success = false(n_samples, 1);
            
            all_metrics = cell(n_samples, 1);
            
            for i = 1:n_samples
                sample = obj.results.samples{i};
                if isstruct(sample) && isfield(sample, 'success') && sample.success
                    success_flags(i) = true;
                    unstable_flags(i) = sample.is_unstable;
                    dead_flags(i) = sample.is_dead;
                    all_metrics{i} = sample.metrics;
                    
                    % Track learning evaluation
                    if isfield(sample, 'learning') && isstruct(sample.learning)
                        if isfield(sample.learning, 'evaluated') && sample.learning.evaluated
                            learning_evaluated(i) = true;
                        end
                        if isfield(sample.learning, 'success') && sample.learning.success
                            learning_success(i) = true;
                        end
                    end
                end
            end
            
            results.success_rate = mean(success_flags);
            results.unstable_rate = mean(unstable_flags);
            results.dead_rate = mean(dead_flags);
            results.stable_rate = mean(success_flags & ~unstable_flags & ~dead_flags);
            results.stable_count = sum(success_flags & ~unstable_flags & ~dead_flags);
            
            % Learning evaluation counts
            results.learning_evaluated_count = sum(learning_evaluated);
            results.learning_success_count = sum(learning_success);
            
            % Aggregate metrics statistics (includes learning metrics if available)
            if any(success_flags)
                results.metrics_stats = obj.computeMetricStatistics(all_metrics(success_flags));
            end
            
            results.config = struct();
            results.config.n_samples = obj.n_samples;
            results.config.random_seed = obj.random_seed;
            results.config.analysis_start_time = obj.analysis_start_time;
            results.config.analysis_end_time = datetime('now');
            results.config.evaluate_learning = obj.evaluate_learning;
            results.config.learning_task = obj.learning_task;
        end
        
        function stats = computeMetricStatistics(obj, metrics_cell)
            % COMPUTEMETRICSTATISTICS Compute statistics across all metrics
            %
            % Handles cases where metrics may have different fields (e.g.,
            % learning metrics only present for stable simulations)
            
            % Get valid (non-empty) metrics
            valid_metrics = metrics_cell(~cellfun(@isempty, metrics_cell));
            if isempty(valid_metrics)
                stats = struct();
                return;
            end
            
            % Collect all unique field names from all metrics
            all_field_names = {};
            for i = 1:length(valid_metrics)
                if isstruct(valid_metrics{i})
                    all_field_names = union(all_field_names, fieldnames(valid_metrics{i}));
                end
            end
            
            stats = struct();
            
            for i = 1:length(all_field_names)
                field_name = all_field_names{i};
                
                % Extract values only from metrics that have this field
                values = {};
                for j = 1:length(valid_metrics)
                    m = valid_metrics{j};
                    if isstruct(m) && isfield(m, field_name)
                        values{end+1} = m.(field_name); %#ok<AGROW>
                    end
                end
                
                % Check which values are numeric, scalar, and finite
                numeric_values = cellfun(@(v) isnumeric(v) && isscalar(v) && isfinite(v), values);
                
                if any(numeric_values)
                    numeric_array = cell2mat(values(numeric_values));
                    stats.(field_name) = struct();
                    stats.(field_name).mean = mean(numeric_array);
                    stats.(field_name).std = std(numeric_array);
                    stats.(field_name).min = min(numeric_array);
                    stats.(field_name).max = max(numeric_array);
                    stats.(field_name).median = median(numeric_array);
                    stats.(field_name).prctile_25 = prctile(numeric_array, 25);
                    stats.(field_name).prctile_75 = prctile(numeric_array, 75);
                    stats.(field_name).count = sum(numeric_values);  % How many samples had this metric
                end
            end
        end
        
        function saveResults(obj, results)
            % SAVERESULTS Save results to disk
            
            if ~obj.save_results
                return;
            end
            
            % Save full results
            save_file = fullfile(obj.output_dir, 'results_full.mat');
            save(save_file, 'results', '-v7.3');
            
            % Save summary
            summary = struct();
            summary.n_samples = obj.n_samples;
            summary.success_rate = results.success_rate;
            summary.unstable_rate = results.unstable_rate;
            summary.dead_rate = results.dead_rate;
            summary.stable_rate = results.stable_rate;
            summary.param_names = obj.param_names;
            summary.param_ranges = obj.param_ranges;
            summary.fixed_params = obj.fixed_params;
            summary.config = results.config;
            
            if isfield(results, 'metrics_stats')
                summary.metrics_stats = results.metrics_stats;
            end
            
            summary_file = fullfile(obj.output_dir, 'summary.mat');
            save(summary_file, 'summary', '-v7.3');
            
            if obj.verbose
                fprintf('Results saved to: %s\n', obj.output_dir);
            end
        end
        
        function printErrorSummary(obj)
            % PRINTERRORSUMMARY Print a summary of all errors encountered
            
            if isempty(obj.results) || ~isfield(obj.results, 'samples')
                return;
            end
            
            fprintf('\n========================================\n');
            fprintf('=== Error Summary ===\n');
            fprintf('========================================\n');
            
            n_total = length(obj.results.samples);
            n_success = 0;
            n_failed = 0;
            error_types = containers.Map('KeyType', 'char', 'ValueType', 'double');
            failed_indices = [];
            
            for i = 1:n_total
                r = obj.results.samples{i};
                if isempty(r)
                    n_failed = n_failed + 1;
                    failed_indices(end+1) = i; %#ok<AGROW>
                    continue;
                end
                
                if isfield(r, 'success') && r.success
                    n_success = n_success + 1;
                else
                    n_failed = n_failed + 1;
                    failed_indices(end+1) = i; %#ok<AGROW>
                    
                    % Categorize error
                    if isfield(r, 'error_message') && ~isempty(r.error_message)
                        % Extract first line of error message as category
                        msg = r.error_message;
                        newline_idx = find(msg == newline, 1, 'first');
                        if ~isempty(newline_idx)
                            msg = msg(1:newline_idx-1);
                        end
                        if length(msg) > 100
                            msg = [msg(1:100) '...'];
                        end
                        
                        if isKey(error_types, msg)
                            error_types(msg) = error_types(msg) + 1;
                        else
                            error_types(msg) = 1;
                        end
                    else
                        if isKey(error_types, 'Unknown error')
                            error_types('Unknown error') = error_types('Unknown error') + 1;
                        else
                            error_types('Unknown error') = 1;
                        end
                    end
                end
            end
            
            fprintf('Total simulations: %d\n', n_total);
            fprintf('Successful: %d (%.1f%%)\n', n_success, 100*n_success/n_total);
            fprintf('Failed: %d (%.1f%%)\n', n_failed, 100*n_failed/n_total);
            
            if n_failed > 0
                fprintf('\n--- Error Types ---\n');
                error_keys = keys(error_types);
                for i = 1:length(error_keys)
                    err_key = error_keys{i};
                    err_count = error_types(err_key);
                    fprintf('  [%d occurrences] %s\n', err_count, err_key);
                end
                
                % Print detailed info for first few failures
                max_details = 5;
                fprintf('\n--- First %d Failed Simulations (Details) ---\n', min(max_details, n_failed));
                detail_count = 0;
                for i = 1:length(failed_indices)
                    if detail_count >= max_details
                        break;
                    end
                    
                    sample_idx = failed_indices(i);
                    r = obj.results.samples{sample_idx};
                    
                    if isempty(r)
                        fprintf('\nSample %d: Empty result (no simulation attempted)\n', sample_idx);
                        detail_count = detail_count + 1;
                        continue;
                    end
                    
                    fprintf('\n--- Sample %d ---\n', sample_idx);
                    if isfield(r, 'error_message')
                        fprintf('Error: %s\n', r.error_message);
                    end
                    if isfield(r, 'error_identifier') && ~isempty(r.error_identifier)
                        fprintf('Identifier: %s\n', r.error_identifier);
                    end
                    if isfield(r, 'error_stack') && ~isempty(r.error_stack)
                        fprintf('Stack trace:\n');
                        for j = 1:min(10, length(r.error_stack))
                            fprintf('  %s (line %d)\n', r.error_stack(j).name, r.error_stack(j).line);
                        end
                    end
                    if isfield(r, 'param_values') && ~isempty(r.param_values)
                        fprintf('Parameters:\n');
                        for k = 1:length(obj.param_names)
                            fprintf('  %s = %.6g\n', obj.param_names{k}, r.param_values(k));
                        end
                    end
                    
                    detail_count = detail_count + 1;
                end
            end
            
            fprintf('\n========================================\n');
        end
    end
end

function val = getFieldOrDefault(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end

function result = runSingleSimulationStatic(param_values, use_gpu, param_names, param_ranges, fixed_params, categorical_params, integer_params, learning_settings)
    %RUNSINGLESIMULATIONSTATIC Static function to run a single simulation (for parfor)
    %
    % This is a standalone function that can be called from within parfor loops
    % since parfor cannot call instance methods directly.
    
    % Default integer_params if not provided
    if nargin < 7 || isempty(integer_params)
        integer_params = {'n', 'n_a_E', 'n_a_I', 'n_b_E', 'n_b_I'};
    end
    
    % Default learning_settings if not provided
    if nargin < 8 || isempty(learning_settings)
        learning_settings = struct('evaluate', false);
    end
    
    result = struct();
    result.param_values = param_values;  % Store raw parameter values for debugging
    
    try
        % Step 1: Build parameter struct
        params = buildParamStructStatic(param_values, param_names, param_ranges, fixed_params, categorical_params, integer_params);
        result.params = params;  % Store params even if later steps fail
        
        % Step 2: Run dynamics simulation
        [states_r, S_history, metrics, stability, esn] = executeSimulationStatic(params, use_gpu);
        
        % Package dynamics results
        result.success = true;
        result.states_r = [];  % Don't store large arrays to save memory
        result.S_history = [];  % Don't store large arrays to save memory
        result.metrics = metrics;
        result.stability = stability;
        result.is_unstable = stability.is_unstable;
        result.is_dead = stability.is_dead;
        result.reason = stability.reason;
        result.error_message = '';
        result.error_identifier = '';
        result.error_stack = [];
        
        % Step 3: Learning evaluation (only if stable and not dead)
        result.learning = struct();
        result.learning.evaluated = false;
        result.learning.success = false;
        
        if learning_settings.evaluate && ~stability.is_unstable && ~stability.is_dead
            try
                learning_options = struct();
                learning_options.task = learning_settings.task;
                learning_options.n_train = learning_settings.n_train;
                learning_options.n_test = learning_settings.n_test;
                learning_options.washout = learning_settings.washout;
                learning_options.verbose = false;
                
                learning_metrics = evaluateLearningCapability(esn, params, learning_options);
                
                result.learning = learning_metrics;
                result.learning.evaluated = true;
                
                % Add key learning metrics to main metrics struct
                if learning_metrics.success
                    result.metrics.learning_test_mse = learning_metrics.test_mse;
                    result.metrics.learning_test_nrmse = learning_metrics.test_nrmse;
                    result.metrics.learning_R2 = learning_metrics.R2;
                    if isfield(learning_metrics, 'memory_capacity') && ~isnan(learning_metrics.memory_capacity)
                        result.metrics.memory_capacity = learning_metrics.memory_capacity;
                    end
                end
            catch ME_learn
                result.learning.evaluated = true;
                result.learning.success = false;
                result.learning.error_message = ME_learn.message;
            end
        elseif stability.is_unstable || stability.is_dead
            % Mark as not learnable
            result.learning.evaluated = false;
            result.learning.reason = 'reservoir_unstable_or_dead';
        end
        
    catch ME
        % Package error result with full details
        result.success = false;
        result.error_message = ME.message;
        result.error_identifier = ME.identifier;
        result.error_stack = ME.stack;
        result.is_unstable = true;
        result.is_dead = false;
        result.reason = 'simulation_error';
        result.metrics = struct();
        result.stability = struct('is_unstable', true, 'is_dead', false, 'reason', 'simulation_error', 'metrics', struct());
        result.states_r = [];
        result.S_history = [];
        result.learning = struct('evaluated', false, 'success', false, 'reason', 'simulation_error');
        
        % Try to add more context if we have params
        if ~isfield(result, 'params')
            result.params = struct();
        end
    end
end

function params = buildParamStructStatic(param_values, param_names, param_ranges, fixed_params, categorical_params, integer_params)
    %BUILDPARAMSTRUCTSTATIC Build parameter struct from values (static version)
    
    % Default integer_params if not provided
    if nargin < 6 || isempty(integer_params)
        integer_params = {'n', 'n_a_E', 'n_a_I', 'n_b_E', 'n_b_I'};
    end
    
    params = fixed_params;  % Start with fixed parameters
    
    % Add sampled parameters
    for i = 1:length(param_names)
        param_name = param_names{i};
        param_value = param_values(i);
        
        % Handle categorical parameters
        if ismember(param_name, categorical_params)
            param_range = param_ranges.(param_name);
            if iscell(param_range)
                idx = round(param_value);
                idx = max(1, min(length(param_range), idx));
                params.(param_name) = param_range{idx};
            else
                params.(param_name) = param_value;
            end
        % Handle integer parameters - ensure they are integers
        elseif ismember(param_name, integer_params)
            params.(param_name) = round(param_value);
        else
            params.(param_name) = param_value;
        end
    end
    
    % Set defaults for required parameters
    params = setDefaultParamsStatic(params);
end

function params = setDefaultParamsStatic(params)
    %SETDEFAULTPARAMSSTATIC Set default values for required parameters (static version)
    
    defaults = struct();
    defaults.dt = 0.1;
    defaults.n_steps = 300;
    defaults.N = 10;  % Silence before stimulus (seconds)
    defaults.M = 200;  % Silence after stimulus (seconds)
    defaults.input_amplitude = 1;
    defaults.S_a = 0.85;
    defaults.S_c = 0.4;
    defaults.pulse_period = 20;
    defaults.pulse_width = 5;
    defaults.input_frequency = 0.001;
    
    default_fields = fieldnames(defaults);
    for i = 1:length(default_fields)
        field_name = default_fields{i};
        if ~isfield(params, field_name)
            params.(field_name) = defaults.(field_name);
        end
    end
    
    % Ensure required fields exist
    if ~isfield(params, 'n') || ~isfield(params, 'fraction_E')
        error('ParameterExploration:MissingRequired', ...
            'n and fraction_E are required parameters');
    end
    
    % Ensure integer parameters are actually integers
    params.n = round(params.n);
    if isfield(params, 'n_a_E')
        params.n_a_E = round(params.n_a_E);
    else
        params.n_a_E = 0;
    end
    if isfield(params, 'n_a_I')
        params.n_a_I = round(params.n_a_I);
    else
        params.n_a_I = 0;
    end
    if isfield(params, 'n_b_E')
        params.n_b_E = round(params.n_b_E);
    else
        params.n_b_E = 0;
    end
    if isfield(params, 'n_b_I')
        params.n_b_I = round(params.n_b_I);
    else
        params.n_b_I = 0;
    end
    
    % Compute derived parameters
    params.n_E = round(params.n * params.fraction_E);
    params.n_I = params.n - params.n_E;
    params.E_indices = 1:params.n_E;
    params.I_indices = (params.n_E+1):params.n;
    
    % Set activation function
    if ~isfield(params, 'activation_function')
        params.activation_function = @(x) piecewiseSigmoid(x, params.S_a, params.S_c);
        params.activation_function_derivative = @(x) piecewiseSigmoidDerivative(x, params.S_a, params.S_c);
    end
    
    % Set adaptation timescales if not specified
    if params.n_a_E > 0 && ~isfield(params, 'tau_a_E')
        params.tau_a_E = logspace(log10(0.25), log10(25), params.n_a_E);
    end
    if params.n_a_I > 0 && ~isfield(params, 'tau_a_I')
        params.tau_a_I = logspace(log10(0.25), log10(25), params.n_a_I);
    end
end

function [states_r, S_history, metrics, stability, esn] = executeSimulationStatic(params, use_gpu)
    %EXECUTESIMULATIONSTATIC Execute a single reservoir simulation (static version)
    %
    % Each step is wrapped with error context for better debugging.
    % Returns the ESN object for subsequent learning evaluation.
    
    current_step = 'initialization';
    debug_info = '';
    esn = [];  % Initialize to empty in case of early error
    
    try
        % Step 1: Generate input signal
        current_step = 'generateInputSignal';
        [U, ~, silence_start_idx] = generateInputSignalStatic(params);
        debug_info = sprintf('U size: [%d, %d]', size(U, 1), size(U, 2));
        
        % Step 2: Build reservoir
        current_step = 'buildReservoir';
        [W, W_in] = buildReservoirStatic(params);
        debug_info = sprintf('%s, W size: [%d, %d], W_in size: [%d, %d]', ...
            debug_info, size(W, 1), size(W, 2), size(W_in, 1), size(W_in, 2));
        
        % Step 3: Create ESN parameters
        current_step = 'buildESNParams';
        esn_params = buildESNParamsStatic(params, W, W_in);
        debug_info = sprintf('%s, n=%d, n_E=%d, n_I=%d, n_a_E=%d, n_a_I=%d, n_b_E=%d, n_b_I=%d', ...
            debug_info, esn_params.n, esn_params.n_E, esn_params.n_I, ...
            getFieldOrDefault(esn_params, 'n_a_E', 0), getFieldOrDefault(esn_params, 'n_a_I', 0), ...
            getFieldOrDefault(esn_params, 'n_b_E', 0), getFieldOrDefault(esn_params, 'n_b_I', 0));
        
        % Step 4: Create and initialize ESN
        current_step = 'createESN';
        esn = SRNN_ESN(esn_params);
        
        current_step = 'resetState';
        esn.resetState();
        debug_info = sprintf('%s, S0 size: [%d, %d]', debug_info, size(esn.S0, 1), size(esn.S0, 2));
        
        % Step 5: Run reservoir
        current_step = 'runReservoir';
        % Clear persistent variables in reservoir functions to avoid dimension mismatch
        % when different simulations have the same time steps but different n
        clear SRNN_reservoir SRNN_reservoir_DDE;
        [~, S_history] = esn.runReservoir(U);
        
        % Step 6: Extract firing rates
        current_step = 'extractFiringRates';
        states_r = extractFiringRatesStatic(S_history, params);
        
        % Step 7: Check stability
        current_step = 'checkStability';
        stability = checkReservoirStability(states_r, silence_start_idx);
        
        % Step 8: Extract metrics
        current_step = 'extractMetrics';
        metrics = extractReservoirMetrics(states_r, S_history, params, silence_start_idx);
        
    catch ME
        % Re-throw with additional context about which step failed
        % Format debug info compactly
        param_info = sprintf('n=%d n_E=%d n_I=%d n_a_E=%d n_a_I=%d n_b_E=%d n_b_I=%d tau_d=%.3f lags=%.4f', ...
            params.n, params.n_E, params.n_I, ...
            getFieldOrDefault(params, 'n_a_E', 0), getFieldOrDefault(params, 'n_a_I', 0), ...
            getFieldOrDefault(params, 'n_b_E', 0), getFieldOrDefault(params, 'n_b_I', 0), ...
            params.tau_d, getFieldOrDefault(params, 'lags', 0));
        error('ParameterExploration:SimulationFailed', ...
            'Step "%s" failed: %s | %s | Sizes: %s | Params: %s', ...
            current_step, ME.message, ME.identifier, debug_info, param_info);
    end
end

function [U, t, silence_start_idx] = generateInputSignalStatic(params)
    %GENERATEINPUTSIGNALSTATIC Generate input signal based on parameters (static version)
    
    dt = params.dt;
    n_steps = params.n_steps;
    N = params.N;  % Silence before
    M = params.M;  % Silence after
    
    n_silence_steps_start = round(N / dt);
    n_silence_steps_end = round(M / dt);
    total_steps = n_steps + n_silence_steps_start + n_silence_steps_end;
    
    t = (0:n_steps-1)' * dt;
    U = zeros(total_steps, 1);
    
    % Generate input based on type
    input_type = getFieldOrDefault(params, 'input_type', 'pulse');
    input_amplitude = getFieldOrDefault(params, 'input_amplitude', 1);
    
    switch lower(input_type)
        case 'sine'
            input_frequency = getFieldOrDefault(params, 'input_frequency', 0.001);
            U(n_silence_steps_start+1:n_silence_steps_start+n_steps) = ...
                input_amplitude * sin(2*pi*input_frequency*t);
            
        case 'pulse'
            pulse_period = getFieldOrDefault(params, 'pulse_period', 20);
            pulse_width = getFieldOrDefault(params, 'pulse_width', 5);
            pulse_signal = zeros(n_steps, 1);
            pulse_signal(mod(0:n_steps-1, pulse_period) < pulse_width) = input_amplitude;
            U(n_silence_steps_start+1:n_silence_steps_start+n_steps) = pulse_signal;
            
        case 'step'
            step_signal = zeros(n_steps, 1);
            step_signal(round(n_steps/4):end) = input_amplitude;
            U(n_silence_steps_start+1:n_silence_steps_start+n_steps) = step_signal;
            
        case 'noise'
            U(n_silence_steps_start+1:n_silence_steps_start+n_steps) = ...
                input_amplitude * randn(n_steps, 1);
            
        case 'mackey'
            try
                [~, mg] = generate_mackey_glass('tau', 17, 'n_samples', n_steps+1000, ...
                    'dt', dt, 'discard', 1000);
                U(n_silence_steps_start+1:n_silence_steps_start+n_steps) = ...
                    input_amplitude * mg(1:n_steps);
            catch
                % Fallback to pulse if mackey-glass fails
                pulse_signal = zeros(n_steps, 1);
                pulse_signal(mod(0:n_steps-1, 20) < 5) = input_amplitude;
                U(n_silence_steps_start+1:n_silence_steps_start+n_steps) = pulse_signal;
            end
            
        otherwise
            % Default to pulse
            pulse_signal = zeros(n_steps, 1);
            pulse_signal(mod(0:n_steps-1, 20) < 5) = input_amplitude;
            U(n_silence_steps_start+1:n_silence_steps_start+n_steps) = pulse_signal;
    end
    
    U = U(:);
    t = (0:length(U)-1)' * dt;
    silence_start_idx = n_silence_steps_start + n_steps + 1;
end

function [W, W_in] = buildReservoirStatic(params)
    %BUILDRESERVOIRSTATIC Build reservoir weight matrices (static version)
    
    n = params.n;
    n_E = params.n_E;
    level_of_chaos = params.level_of_chaos;
    input_scaling = getFieldOrDefault(params, 'input_scaling', 0.75);
    
    % Create recurrent weight matrix
    W = randn(n, n);
    
    % Apply Dale's law
    W(:, 1:n_E) = abs(W(:, 1:n_E));           % E columns positive
    W(:, n_E+1:end) = -abs(W(:, n_E+1:end));  % I columns negative
    
    % Center rows
    W = W - mean(W, 2);
    
    % Apply chaos scaling
    W_eigs = eig(W);
    abscissa_0 = max(real(W_eigs));
    if abs(abscissa_0) > 1e-10
        gamma = 1 / abscissa_0;
        W = level_of_chaos * gamma * W;
    end
    
    % Create input weight matrix
    n_inputs = 1;  % Single input dimension
    W_in = (2 * rand(n, n_inputs) - 1) * input_scaling;
    W_in(rand(n, n_inputs) > 0.2) = 0;  % Sparse input
end

function esn_params = buildESNParamsStatic(params, W, W_in)
    %BUILDESNPARAMSSTATIC Build parameters struct for SRNN_ESN (static version)
    
    esn_params = struct();
    esn_params.n = params.n;
    esn_params.n_E = params.n_E;
    esn_params.n_I = params.n_I;
    esn_params.W = W;
    esn_params.W_in = W_in;
    esn_params.tau_d = params.tau_d;
    esn_params.activation_function = params.activation_function;
    esn_params.activation_function_derivative = params.activation_function_derivative;
    esn_params.E_indices = params.E_indices;
    esn_params.I_indices = params.I_indices;
    esn_params.which_states = 'x';
    esn_params.include_input = false;
    esn_params.lambda = 1e-6;
    esn_params.dt = params.dt;
    
    % Optional parameters
    if isfield(params, 'n_a_E')
        esn_params.n_a_E = params.n_a_E;
    end
    if isfield(params, 'n_a_I')
        esn_params.n_a_I = params.n_a_I;
    end
    if isfield(params, 'tau_a_E')
        esn_params.tau_a_E = params.tau_a_E;
    end
    if isfield(params, 'tau_a_I')
        esn_params.tau_a_I = params.tau_a_I;
    end
    if isfield(params, 'n_b_E')
        esn_params.n_b_E = params.n_b_E;
    end
    if isfield(params, 'n_b_I')
        esn_params.n_b_I = params.n_b_I;
    end
    if isfield(params, 'tau_b_E_rec')
        esn_params.tau_b_E_rec = params.tau_b_E_rec;
    end
    if isfield(params, 'tau_b_E_rel')
        esn_params.tau_b_E_rel = params.tau_b_E_rel;
    end
    if isfield(params, 'tau_b_I_rec')
        esn_params.tau_b_I_rec = params.tau_b_I_rec;
    end
    if isfield(params, 'tau_b_I_rel')
        esn_params.tau_b_I_rel = params.tau_b_I_rel;
    end
    if isfield(params, 'c_E')
        esn_params.c_E = params.c_E;
    end
    if isfield(params, 'c_I')
        esn_params.c_I = params.c_I;
    end
    if isfield(params, 'lags')
        esn_params.lags = params.lags;
    end
end

function states_r = extractFiringRatesStatic(S_history, params)
    %EXTRACTFIRINGRATESSTATIC Extract firing rates from state history (static version)
    
    n = params.n;
    n_E = params.n_E;
    n_I = params.n_I;
    n_a_E = getFieldOrDefault(params, 'n_a_E', 0);
    n_a_I = getFieldOrDefault(params, 'n_a_I', 0);
    n_b_E = getFieldOrDefault(params, 'n_b_E', 0);
    n_b_I = getFieldOrDefault(params, 'n_b_I', 0);
    c_E = getFieldOrDefault(params, 'c_E', 0);
    c_I = getFieldOrDefault(params, 'c_I', 0);
    activation_function = params.activation_function;
    
    total_steps = size(S_history, 1);
    states_r = zeros(total_steps, n);
    
    len_a_E = n_E * n_a_E;
    len_a_I = n_I * n_a_I;
    len_b_E = n_E * n_b_E;
    len_b_I = n_I * n_b_I;
    
    for i = 1:total_steps
        current_idx = 1;
        
        % Extract adaptation variables a_E
        if n_a_E > 0
            a_E = reshape(S_history(i, current_idx:current_idx+len_a_E-1), n_a_E, n_E)';
            current_idx = current_idx + len_a_E;
        else
            a_E = zeros(n_E, 0);
        end
        
        % Extract adaptation variables a_I
        if n_a_I > 0
            a_I = reshape(S_history(i, current_idx:current_idx+len_a_I-1), n_a_I, n_I)';
            current_idx = current_idx + len_a_I;
        else
            a_I = zeros(n_I, 0);
        end
        
        % Skip depression variables
        current_idx = current_idx + len_b_E + len_b_I;
        
        % Extract dendritic states x
        x = S_history(i, current_idx:current_idx+n-1)';
        
        % Compute effective input with adaptation
        x_eff = x;
        if n_a_E > 0 && n_E > 0
            x_eff(1:n_E) = x_eff(1:n_E) - c_E * sum(a_E, 2);
        end
        if n_a_I > 0 && n_I > 0
            x_eff(n_E+1:end) = x_eff(n_E+1:end) - c_I * sum(a_I, 2);
        end
        
        % Apply activation function
        states_r(i, :) = activation_function(x_eff);
    end
end

