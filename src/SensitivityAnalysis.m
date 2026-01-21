classdef SensitivityAnalysis < handle
    % SENSITIVITYANALYSIS Performs parameter sweeps using SRNNModel
    %
    % This class runs sensitivity analysis by varying specified parameters
    % across multiple levels and repetitions, optionally comparing different
    % adaptation conditions (no_adaptation, sfa_only, std_only, sfa_and_std).
    %
    % Usage:
    %   sa = SensitivityAnalysis('n_levels', 11, 'n_reps', 20);
    %   sa.add_parameter('level_of_chaos', [0.5, 3.0]);
    %   sa.add_parameter('n', [50, 200]);
    %   sa.run();
    %   sa.save_results();
    %   sa.plot();
    %
    % See also: SRNNModel

    %% Configuration Properties
    properties
        param_ranges = struct()     % struct: param_name -> [min, max]
        n_levels = 21               % Number of levels per parameter
        n_reps = 50                 % Number of repetitions per level
        conditions                  % Cell array of condition structs
        integer_params = {'n', 'indegree', 'n_a_E', 'n_a_I', 'n_b_E', 'n_b_I'}
    end

    %% Model Default Properties
    properties
        model_defaults = struct()   % Default SRNNModel properties to apply
        verbose = true              % Whether to print progress
    end

    %% Output Properties
    properties
        output_dir                  % Base directory for saving results
        note = ''                   % Optional note for output folder naming
    end

    %% Results (SetAccess = private)
    properties (SetAccess = private)
        results = struct()          % Stored results after run
        has_run = false             % Flag indicating if analysis has run
        analysis_start_time         % Timestamp when analysis started
    end

    %% Constructor
    methods
        function obj = SensitivityAnalysis(varargin)
            % SENSITIVITYANALYSIS Constructor with name-value pairs
            %
            % Usage:
            %   sa = SensitivityAnalysis()  % All defaults
            %   sa = SensitivityAnalysis('n_levels', 11, 'n_reps', 20)

            % Set default conditions
            obj.set_default_conditions();

            % Parse name-value pairs
            for i = 1:2:length(varargin)
                if isprop(obj, varargin{i})
                    obj.(varargin{i}) = varargin{i+1};
                else
                    warning('SensitivityAnalysis:UnknownProperty', ...
                        'Unknown property: %s', varargin{i});
                end
            end

            % Set default output directory
            if isempty(obj.output_dir)
                % Default to 'data/sensitivity' in the project root
                % Assuming this class is in src/, go up one level to project root
                src_path = fileparts(mfilename('fullpath'));
                project_root = fileparts(src_path);
                obj.output_dir = fullfile(project_root, 'data', 'sensitivity');
            end
        end
    end

    %% Public Methods
    methods
        function add_parameter(obj, param_name, param_range)
            % ADD_PARAMETER Add a parameter to sweep
            %
            % Usage:
            %   sa.add_parameter('level_of_chaos', [0.5, 3.0])
            %   sa.add_parameter('n', [50, 200])

            if ~ischar(param_name) && ~isstring(param_name)
                error('SensitivityAnalysis:InvalidInput', ...
                    'param_name must be a string or char array');
            end

            if ~isnumeric(param_range) || length(param_range) ~= 2
                error('SensitivityAnalysis:InvalidInput', ...
                    'param_range must be a 2-element numeric array [min, max]');
            end

            if param_range(2) <= param_range(1)
                error('SensitivityAnalysis:InvalidInput', ...
                    'param_range(2) must be > param_range(1)');
            end

            obj.param_ranges.(param_name) = param_range;

            if obj.verbose
                fprintf('Added parameter: %s, range: [%.3g, %.3g]\n', ...
                    param_name, param_range(1), param_range(2));
            end
        end

        function remove_parameter(obj, param_name)
            % REMOVE_PARAMETER Remove a parameter from the sweep
            %
            % Usage:
            %   sa.remove_parameter('level_of_chaos')

            if isfield(obj.param_ranges, param_name)
                obj.param_ranges = rmfield(obj.param_ranges, param_name);
                if obj.verbose
                    fprintf('Removed parameter: %s\n', param_name);
                end
            else
                warning('SensitivityAnalysis:ParamNotFound', ...
                    'Parameter %s not found in param_ranges', param_name);
            end
        end

        function set_conditions(obj, conditions_cell)
            % SET_CONDITIONS Set custom conditions for the analysis
            %
            % Usage:
            %   sa.set_conditions({
            %       struct('name', 'custom1', 'n_a_E', 2, 'n_b_E', 0),
            %       struct('name', 'custom2', 'n_a_E', 0, 'n_b_E', 2)
            %   })

            obj.conditions = conditions_cell;
            if obj.verbose
                fprintf('Set %d custom conditions\n', length(conditions_cell));
            end
        end

        function run(obj)
            % RUN Execute the full sensitivity analysis
            %
            % This method:
            %   1. Creates output directory with timestamp
            %   2. Loops over conditions and parameters
            %   3. Uses parfor for parallel execution
            %   4. Stores results in obj.results

            % Validate
            param_names = fieldnames(obj.param_ranges);
            if isempty(param_names)
                error('SensitivityAnalysis:NoParameters', ...
                    'No parameters to sweep. Use add_parameter() first.');
            end

            if isempty(obj.conditions)
                error('SensitivityAnalysis:NoConditions', ...
                    'No conditions defined. Set conditions or use defaults.');
            end

            % Create output directory with timestamp
            obj.analysis_start_time = datetime('now');
            dt_str = lower(strrep(datestr(obj.analysis_start_time, 'mmm_dd_yy_HH_MM_AM'), ':', '_'));

            if ~isempty(obj.note)
                folder_name = sprintf('sensitivity_%s_nLevs_%d_nReps_%d_%s', ...
                    obj.note, obj.n_levels, obj.n_reps, dt_str);
            else
                folder_name = sprintf('sensitivity_nLevs_%d_nReps_%d_%s', ...
                    obj.n_levels, obj.n_reps, dt_str);
            end

            obj.output_dir = fullfile(obj.output_dir, folder_name);
            if ~exist(obj.output_dir, 'dir')
                mkdir(obj.output_dir);
            end

            % Print analysis summary
            fprintf('\n========================================\n');
            fprintf('=== SRNN Sensitivity Analysis ===\n');
            fprintf('========================================\n');
            fprintf('Parameters to sweep: %s\n', strjoin(param_names, ', '));
            fprintf('Conditions: %s\n', strjoin(cellfun(@(c) c.name, obj.conditions, 'UniformOutput', false), ', '));
            fprintf('Levels per parameter: %d\n', obj.n_levels);
            fprintf('Repetitions per level: %d\n', obj.n_reps);
            fprintf('Total runs per parameter per condition: %d\n', obj.n_levels * obj.n_reps);
            fprintf('Output directory: %s\n', obj.output_dir);
            fprintf('========================================\n\n');

            overall_start = tic;

            % Initialize results structure
            obj.results = struct();

            % Loop over conditions
            for c_idx = 1:length(obj.conditions)
                condition = obj.conditions{c_idx};
                condition_name = condition.name;

                fprintf('\n--- Condition: %s (n_a_E=%d, n_b_E=%d) ---\n', ...
                    upper(condition_name), condition.n_a_E, condition.n_b_E);

                % Create condition output directory
                condition_dir = fullfile(obj.output_dir, condition_name);
                if ~exist(condition_dir, 'dir')
                    mkdir(condition_dir);
                end

                % Initialize results for this condition
                obj.results.(condition_name) = struct();

                % Loop over parameters
                for p_idx = 1:length(param_names)
                    param_name = param_names{p_idx};
                    param_range = obj.param_ranges.(param_name);

                    fprintf('\nParameter: %s [%.3g, %.3g] (%d/%d)\n', ...
                        param_name, param_range(1), param_range(2), p_idx, length(param_names));

                    param_start = tic;

                    % Create parameter levels
                    if ismember(param_name, obj.integer_params)
                        param_levels = round(linspace(param_range(1), param_range(2), obj.n_levels));
                    else
                        param_levels = linspace(param_range(1), param_range(2), obj.n_levels);
                    end

                    % Run the parameter sweep
                    [results_reshaped, metadata] = obj.run_parameter_sweep(...
                        param_name, param_levels, condition, c_idx, p_idx);

                    % Store results
                    obj.results.(condition_name).(param_name) = struct(...
                        'results_reshaped', {results_reshaped}, ...
                        'metadata', metadata);

                    % Save immediately
                    save_filename = fullfile(condition_dir, sprintf('sensitivity_%s.mat', param_name));
                    save(save_filename, 'results_reshaped', 'metadata', '-v7.3');

                    param_elapsed = toc(param_start);
                    fprintf('Parameter %s completed in %.2f minutes. Saved to %s\n', ...
                        param_name, param_elapsed/60, save_filename);
                end
            end

            overall_elapsed = toc(overall_start);
            fprintf('\n========================================\n');
            fprintf('=== Analysis Complete ===\n');
            fprintf('Total time: %.2f hours\n', overall_elapsed/3600);
            fprintf('========================================\n');

            obj.has_run = true;

            % Save summary
            obj.save_summary();
        end

        function save_results(obj)
            % SAVE_RESULTS Save all results to disk (called automatically by run())
            %
            % This method is provided for manual re-saving if needed.

            if ~obj.has_run
                error('SensitivityAnalysis:NotRun', ...
                    'Analysis has not been run yet. Call run() first.');
            end

            obj.save_summary();
            fprintf('Results saved to: %s\n', obj.output_dir);
        end

        function plot(obj, varargin)
            % PLOT Generate sensitivity plots
            %
            % Usage:
            %   sa.plot()
            %   sa.plot('metric', 'LLE')
            %   sa.plot('metric', 'mean_rate')

            if ~obj.has_run
                error('SensitivityAnalysis:NotRun', ...
                    'Analysis has not been run yet. Call run() first.');
            end

            % Parse optional arguments
            metric = 'LLE';  % Default metric
            for i = 1:2:length(varargin)
                if strcmpi(varargin{i}, 'metric')
                    metric = varargin{i+1};
                end
            end

            % Define histogram bins based on metric
            if strcmpi(metric, 'LLE')
                hist_bins = [-inf, linspace(-1.5, 1.5, 25), inf];
                y_label = '$\lambda_1$';
            elseif strcmpi(metric, 'mean_rate')
                hist_bins = [linspace(0, 50, 25), inf];
                y_label = 'Mean Rate (Hz)';
            else
                hist_bins = [-inf, linspace(-10, 10, 25), inf];
                y_label = metric;
            end

            % Get parameters and conditions
            condition_names = fieldnames(obj.results);
            if isempty(condition_names)
                error('SensitivityAnalysis:NoResults', 'No results to plot.');
            end

            param_names = fieldnames(obj.results.(condition_names{1}));

            n_conditions = length(condition_names);
            n_params = length(param_names);

            % Create figure
            fig = figure('Name', sprintf('%s Sensitivity', metric), ...
                'Position', [100, 100, 300 * n_params, 250 * n_conditions]);

            for c_idx = 1:n_conditions
                cond_name = condition_names{c_idx};

                for p_idx = 1:n_params
                    param_name = param_names{p_idx};

                    subplot_idx = (c_idx - 1) * n_params + p_idx;
                    ax = subplot(n_conditions, n_params, subplot_idx);

                    % Get results
                    res = obj.results.(cond_name).(param_name);
                    results_reshaped = res.results_reshaped;
                    metadata = res.metadata;

                    % Plot heatmap
                    obj.plot_sensitivity_heatmap(ax, results_reshaped, metadata, ...
                        hist_bins, metric, y_label, param_name, cond_name, ...
                        p_idx == 1, c_idx == n_conditions);
                end
            end

            % Save figure
            if ~isempty(obj.output_dir)
                fig_dir = fullfile(obj.output_dir, sprintf('%s_figs', metric));
                if ~exist(fig_dir, 'dir')
                    mkdir(fig_dir);
                end
                saveas(fig, fullfile(fig_dir, sprintf('sensitivity_%s.png', metric)));
                saveas(fig, fullfile(fig_dir, sprintf('sensitivity_%s.fig', metric)));
                fprintf('Figures saved to: %s\n', fig_dir);
            end
        end
    end

    %% Private Methods
    methods (Access = private)
        function set_default_conditions(obj)
            % SET_DEFAULT_CONDITIONS Initialize default adaptation conditions

            obj.conditions = { ...
                struct('name', 'no_adaptation', 'n_a_E', 0, 'n_b_E', 0), ...
                struct('name', 'sfa_only',      'n_a_E', 3, 'n_b_E', 0), ...
                struct('name', 'std_only',      'n_a_E', 0, 'n_b_E', 1), ...
                struct('name', 'sfa_and_std',   'n_a_E', 3, 'n_b_E', 1) ...
                };
        end

        function [results_reshaped, metadata] = run_parameter_sweep(obj, param_name, param_levels, condition, c_idx, p_idx)
            % RUN_PARAMETER_SWEEP Execute parfor loop for one parameter

            n_levels_local = obj.n_levels;
            n_reps_local = obj.n_reps;
            total_runs = n_levels_local * n_reps_local;

            % Pre-allocate results cell array
            results = cell(total_runs, 1);

            % Extract condition values
            n_a_E_val = condition.n_a_E;
            n_b_E_val = condition.n_b_E;

            % Extract model defaults
            model_defaults_local = obj.model_defaults;
            verbose_local = obj.verbose;

            % Progress tracking
            n_success = 0;
            n_failed = 0;

            parfor k = 1:total_runs
                level_idx = floor((k-1) / n_reps_local) + 1;
                rep_idx = mod(k-1, n_reps_local) + 1;

                param_value = param_levels(level_idx);

                % Unique seed for each run
                sim_seed = k + (c_idx-1)*total_runs*100 + (p_idx-1)*total_runs*10;

                % Progress reporting (sparse)
                if verbose_local && (mod(k, max(1, floor(total_runs/10))) == 1 || k == total_runs)
                    fprintf('  Run %d/%d: %s=%.3g, rep=%d, seed=%d\n', ...
                        k, total_runs, param_name, param_value, rep_idx, sim_seed);
                end

                run_start = tic;
                try
                    % Create model with condition settings
                    model_args = {'n_a_E', n_a_E_val, 'n_b_E', n_b_E_val, ...
                        'rng_seeds', [sim_seed, sim_seed+1]};

                    % Add the swept parameter
                    model_args = [model_args, {param_name, param_value}];

                    % Add model defaults
                    default_fields = fieldnames(model_defaults_local);
                    for d_idx = 1:length(default_fields)
                        field_name = default_fields{d_idx};
                        % Don't override the swept parameter or condition params
                        if ~strcmp(field_name, param_name) && ...
                                ~strcmp(field_name, 'n_a_E') && ~strcmp(field_name, 'n_b_E')
                            model_args = [model_args, {field_name, model_defaults_local.(field_name)}];
                        end
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

                    % Extract mean firing rate from plot_data
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

                    if verbose_local
                        fprintf('  ERROR in run %d: %s\n', k, ME.message);
                    end
                end
            end

            % Reshape results to n_levels x n_reps
            results_reshaped = reshape(results, [n_reps_local, n_levels_local])';

            % Count successes
            success_mask = cellfun(@(x) isfield(x, 'success') && x.success, results);
            n_success = sum(success_mask);
            n_failed = total_runs - n_success;

            % Create metadata
            metadata = struct();
            metadata.param_name = param_name;
            metadata.param_levels = param_levels;
            metadata.param_range = [param_levels(1), param_levels(end)];
            metadata.n_levels = n_levels_local;
            metadata.n_reps = n_reps_local;
            metadata.n_success = n_success;
            metadata.n_failed = n_failed;
            metadata.success_rate = n_success / total_runs;
            metadata.condition = condition;
            metadata.analysis_date = datestr(now);

            fprintf('  Success rate: %d/%d (%.1f%%)\n', n_success, total_runs, 100*n_success/total_runs);
        end

        function save_summary(obj)
            % SAVE_SUMMARY Save analysis summary to disk

            summary_file = fullfile(obj.output_dir, 'sensitivity_analysis_summary.mat');

            summary_data = struct();
            summary_data.conditions = obj.conditions;
            summary_data.param_ranges = obj.param_ranges;
            summary_data.n_levels = obj.n_levels;
            summary_data.n_reps = obj.n_reps;
            summary_data.model_defaults = obj.model_defaults;
            summary_data.analysis_start_time = obj.analysis_start_time;
            summary_data.analysis_completed = datestr(now);

            % Compute overall statistics
            condition_names = fieldnames(obj.results);
            total_success = 0;
            total_runs = 0;

            all_conditions_summary = struct();
            for c_idx = 1:length(condition_names)
                cond_name = condition_names{c_idx};
                param_names = fieldnames(obj.results.(cond_name));

                cond_summary = struct();
                for p_idx = 1:length(param_names)
                    param_name = param_names{p_idx};
                    meta = obj.results.(cond_name).(param_name).metadata;
                    cond_summary.(param_name) = meta;
                    total_success = total_success + meta.n_success;
                    total_runs = total_runs + (meta.n_success + meta.n_failed);
                end
                all_conditions_summary.(cond_name) = cond_summary;
            end

            summary_data.all_conditions_summary = all_conditions_summary;
            if total_runs > 0
                summary_data.overall_success_rate = total_success / total_runs;
            else
                summary_data.overall_success_rate = 0;
            end

            save(summary_file, 'summary_data', '-v7.3');
            fprintf('Summary saved to: %s\n', summary_file);
        end

        function plot_sensitivity_heatmap(obj, ax, results_reshaped, metadata, ...
                hist_bins, metric, y_label, param_name, cond_name, show_ylabel, show_xlabel)
            % PLOT_SENSITIVITY_HEATMAP Create a single heatmap subplot

            n_levels_local = metadata.n_levels;
            n_reps_local = metadata.n_reps;
            param_levels = metadata.param_levels;

            num_bins = length(hist_bins) - 1;
            histogram_matrix = zeros(num_bins, n_levels_local);

            % Build histogram for each level
            for level_idx = 1:n_levels_local
                values_level = [];

                for rep_idx = 1:n_reps_local
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
                end
            end

            % Normalize by n_reps for probability
            histogram_matrix = histogram_matrix / n_reps_local;

            % Compute y-coordinates for plotting
            finite_edges = hist_bins(~isinf(hist_bins));
            if length(finite_edges) >= 2
                step_size = finite_edges(2) - finite_edges(1);
                y_coords = zeros(num_bins, 1);
                y_coords(1) = finite_edges(1) - step_size/2;
                for k = 2:length(finite_edges)
                    y_coords(k) = (finite_edges(k-1) + finite_edges(k)) / 2;
                end
                y_coords(end) = finite_edges(end) + step_size/2;
            else
                y_coords = 1:num_bins;
            end

            % Plot
            imagesc(ax, param_levels, y_coords, histogram_matrix);
            hold(ax, 'on');
            yline(ax, 0, '--', 'Color', [0 0.7 0], 'LineWidth', 2, 'Alpha', 0.5);
            hold(ax, 'off');

            colormap(ax, flipud(gray));
            caxis(ax, [0, 1]);
            axis(ax, 'xy');
            box(ax, 'on');

            % Labels
            if show_ylabel
                ylabel(ax, y_label, 'Interpreter', 'latex', 'FontSize', 14);

                % Add condition name on left
                yl = ylim(ax);
                xl = xlim(ax);
                text(ax, xl(1) - 0.3*(xl(2)-xl(1)), mean(yl), ...
                    strrep(cond_name, '_', ' '), ...
                    'Rotation', 90, 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', 'FontSize', 12);
            end

            if show_xlabel
                xlabel(ax, strrep(param_name, '_', '\_'), 'FontSize', 12);
            end

            % Limit x-ticks
            xticks(ax, linspace(param_levels(1), param_levels(end), 3));
        end
    end
end

