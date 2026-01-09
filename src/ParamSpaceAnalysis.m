classdef ParamSpaceAnalysis < handle
    % PARAMSPACEANALYSIS Multi-dimensional parameter grid search for SRNN
    %
    % This class performs parameter space exploration by generating all
    % combinations of specified parameters, simulating each configuration
    % across multiple adaptation conditions, and analyzing the results.
    %
    % Key features:
    %   - Multi-dimensional grid (not 1D sweeps like SensitivityAnalysis)
    %   - Same network W used across all 4 conditions for fair comparison
    %   - Randomized execution order for representative early-stopping
    %   - Batched parfor with checkpoint files for resume capability
    %   - Histogram-based visualization of results
    %
    % Usage:
    %   psa = ParamSpaceAnalysis('n_levels', 4, 'note', 'demo');
    %   psa.add_grid_parameter('level_of_chaos', [0.5, 3.0]);
    %   psa.add_grid_parameter('n', [50, 200]);
    %   psa.run();
    %   psa.plot('metric', 'LLE');
    %
    % See also: SRNNModel, SensitivityAnalysis
    
    %% Configuration Properties
    properties
        grid_params = {}            % Cell array of parameter names for grid
        param_ranges = struct()     % Struct: param_name -> [min, max]
        n_levels = 8                % Number of levels per grid parameter
        conditions                  % Cell array of condition structs
        integer_params = {'n', 'indegree', 'n_a_E', 'n_a_I', 'n_b_E', 'n_b_I'}
    end
    
    %% Model Default Properties
    properties
        model_defaults = struct()   % Default SRNNModel properties
        verbose = true              % Print progress during execution
    end
    
    %% Execution Properties
    properties
        batch_size = 25             % Number of configs per batch
        output_dir                  % Base directory for saving results
        note = ''                   % Optional note for folder naming
    end
    
    %% Results (SetAccess = private)
    properties (SetAccess = private)
        results = struct()          % Stored results after run
        has_run = false             % Flag indicating if analysis has run
        analysis_start_time         % Timestamp when analysis started
        param_vectors               % Cell array of parameter level vectors
        all_configs                 % Cell array of all config structs
        shuffled_indices            % Randomized order for execution
        num_combinations            % Total number of grid points
    end
    
    %% Constructor
    methods
        function obj = ParamSpaceAnalysis(varargin)
            % PARAMSPACEANALYSIS Constructor with name-value pairs
            %
            % Usage:
            %   psa = ParamSpaceAnalysis()
            %   psa = ParamSpaceAnalysis('n_levels', 8, 'batch_size', 50)
            
            % Set default conditions
            obj.set_default_conditions();
            
            % Parse name-value pairs
            for i = 1:2:length(varargin)
                if isprop(obj, varargin{i})
                    obj.(varargin{i}) = varargin{i+1};
                else
                    warning('ParamSpaceAnalysis:UnknownProperty', ...
                        'Unknown property: %s', varargin{i});
                end
            end
            
            % Set default output directory
            if isempty(obj.output_dir)
                obj.output_dir = pwd;
            end
        end
    end
    
    %% Public Methods
    methods
        function add_grid_parameter(obj, param_name, param_range)
            % ADD_GRID_PARAMETER Add a parameter to the multi-dimensional grid
            %
            % Usage:
            %   psa.add_grid_parameter('level_of_chaos', [0.5, 3.0])
            
            if ~ischar(param_name) && ~isstring(param_name)
                error('ParamSpaceAnalysis:InvalidInput', ...
                    'param_name must be a string or char array');
            end
            
            if ~isnumeric(param_range) || length(param_range) ~= 2
                error('ParamSpaceAnalysis:InvalidInput', ...
                    'param_range must be a 2-element numeric array [min, max]');
            end
            
            if param_range(2) <= param_range(1)
                error('ParamSpaceAnalysis:InvalidInput', ...
                    'param_range(2) must be > param_range(1)');
            end
            
            % Add to grid_params if not already present
            if ~ismember(param_name, obj.grid_params)
                obj.grid_params{end+1} = param_name;
            end
            
            obj.param_ranges.(param_name) = param_range;
            
            if obj.verbose
                fprintf('Added grid parameter: %s, range: [%.3g, %.3g]\n', ...
                    param_name, param_range(1), param_range(2));
            end
        end
        
        function remove_grid_parameter(obj, param_name)
            % REMOVE_GRID_PARAMETER Remove a parameter from the grid
            
            idx = find(strcmp(obj.grid_params, param_name));
            if ~isempty(idx)
                obj.grid_params(idx) = [];
                if isfield(obj.param_ranges, param_name)
                    obj.param_ranges = rmfield(obj.param_ranges, param_name);
                end
                if obj.verbose
                    fprintf('Removed grid parameter: %s\n', param_name);
                end
            else
                warning('ParamSpaceAnalysis:ParamNotFound', ...
                    'Parameter %s not found in grid', param_name);
            end
        end
        
        function set_conditions(obj, conditions_cell)
            % SET_CONDITIONS Set custom conditions for the analysis
            %
            % Usage:
            %   psa.set_conditions({
            %       struct('name', 'no_adaptation', 'n_a_E', 0, 'n_b_E', 0),
            %       struct('name', 'sfa_only',      'n_a_E', 3, 'n_b_E', 0)
            %   })
            
            obj.conditions = conditions_cell;
            if obj.verbose
                fprintf('Set %d custom conditions\n', length(conditions_cell));
            end
        end
        
        function run(obj)
            % RUN Execute the full parameter space analysis
            %
            % This method:
            %   1. Generates the multi-dimensional parameter grid
            %   2. Randomizes execution order
            %   3. Runs batched parfor with checkpoint files
            %   4. Consolidates results into per-condition MAT files
            
            % Validate
            if isempty(obj.grid_params)
                error('ParamSpaceAnalysis:NoParameters', ...
                    'No grid parameters defined. Use add_grid_parameter() first.');
            end
            
            if isempty(obj.conditions)
                error('ParamSpaceAnalysis:NoConditions', ...
                    'No conditions defined.');
            end
            
            % Create timestamped output directory
            obj.analysis_start_time = datetime('now');
            dt_str = lower(strrep(datestr(obj.analysis_start_time, 'mmm_dd_yy_HH_MM_AM'), ':', '_'));
            
            if ~isempty(obj.note)
                folder_name = sprintf('param_space_%s_nLevs_%d_%s', ...
                    obj.note, obj.n_levels, dt_str);
            else
                folder_name = sprintf('param_space_nLevs_%d_%s', ...
                    obj.n_levels, dt_str);
            end
            
            obj.output_dir = fullfile(obj.output_dir, folder_name);
            if ~exist(obj.output_dir, 'dir')
                mkdir(obj.output_dir);
            end
            
            % Generate parameter grid
            obj.generate_grid();
            
            % Print summary
            fprintf('\n========================================\n');
            fprintf('=== SRNN Parameter Space Analysis ===\n');
            fprintf('========================================\n');
            fprintf('Grid parameters: %s\n', strjoin(obj.grid_params, ', '));
            fprintf('Levels per parameter: %d\n', obj.n_levels);
            fprintf('Total grid combinations: %d\n', obj.num_combinations);
            fprintf('Conditions: %s\n', strjoin(cellfun(@(c) c.name, obj.conditions, 'UniformOutput', false), ', '));
            fprintf('Batch size: %d\n', obj.batch_size);
            fprintf('Output directory: %s\n', obj.output_dir);
            fprintf('========================================\n\n');
            
            % Create temp directory for batch results
            temp_dir = fullfile(obj.output_dir, 'temp_batches');
            if ~exist(temp_dir, 'dir')
                mkdir(temp_dir);
            end
            
            % Create condition output directories
            for c_idx = 1:length(obj.conditions)
                cond_dir = fullfile(obj.output_dir, obj.conditions{c_idx}.name);
                if ~exist(cond_dir, 'dir')
                    mkdir(cond_dir);
                end
            end
            
            overall_start = tic;
            
            % Run batched simulation
            obj.run_batched_simulation(temp_dir);
            
            % Consolidate batch results
            obj.consolidate_results(temp_dir);
            
            overall_elapsed = toc(overall_start);
            fprintf('\n========================================\n');
            fprintf('=== Analysis Complete ===\n');
            fprintf('Total time: %.2f hours\n', overall_elapsed/3600);
            fprintf('========================================\n');
            
            obj.has_run = true;
            
            % Save summary
            obj.save_summary();
        end
        
        function plot(obj, varargin)
            % PLOT Generate histogram plots of metrics across parameter space
            %
            % Usage:
            %   psa.plot()
            %   psa.plot('metric', 'LLE')
            %   psa.plot('metric', 'mean_rate')
            
            if ~obj.has_run && isempty(fieldnames(obj.results))
                error('ParamSpaceAnalysis:NotRun', ...
                    'Analysis has not been run yet. Call run() first.');
            end
            
            % Parse arguments
            metric = 'LLE';
            for i = 1:2:length(varargin)
                if strcmpi(varargin{i}, 'metric')
                    metric = varargin{i+1};
                end
            end
            
            % Define histogram bins
            if strcmpi(metric, 'LLE')
                hist_range = [-1.5, 1.5];
                n_bins = 25;
                y_label = 'LLE (\lambda_1)';
            elseif strcmpi(metric, 'mean_rate')
                hist_range = [0, 10];
                n_bins = 25;
                y_label = 'Mean Firing Rate (Hz)';
            else
                hist_range = [-10, 10];
                n_bins = 25;
                y_label = metric;
            end
            
            hist_bins = [linspace(hist_range(1), hist_range(2), n_bins + 1), inf];
            if strcmpi(metric, 'LLE')
                hist_bins = [-inf, hist_bins];
            end
            
            % Get condition names
            condition_names = cellfun(@(c) c.name, obj.conditions, 'UniformOutput', false);
            num_conditions = length(condition_names);
            
            % Readable titles
            condition_titles = containers.Map(...
                {'no_adaptation', 'sfa_only', 'std_only', 'sfa_and_std'}, ...
                {'No Adaptation', 'SFA Only', 'STD Only', 'SFA + STD'});
            
            % Create figure
            fig = figure('Name', sprintf('%s Distribution', metric), ...
                'Position', [100, 100, 300 * num_conditions, 300]);
            
            for c_idx = 1:num_conditions
                cond_name = condition_names{c_idx};
                ax = subplot(1, num_conditions, c_idx);
                
                % Extract metric values
                if isfield(obj.results, cond_name)
                    results_cell = obj.results.(cond_name);
                    values = [];
                    for k = 1:length(results_cell)
                        res = results_cell{k};
                        if isstruct(res) && isfield(res, 'success') && res.success
                            if isfield(res, metric) && ~isnan(res.(metric))
                                values(end+1) = res.(metric);
                            end
                        end
                    end
                else
                    values = [];
                end
                
                % Plot histogram
                if ~isempty(values)
                    [counts, edges] = histcounts(values, hist_bins);
                    % Normalize to probability
                    prob = counts / sum(counts);
                    
                    % Create finite edges for plotting
                    finite_edges = edges;
                    step = (hist_range(2) - hist_range(1)) / n_bins;
                    if isinf(finite_edges(1))
                        finite_edges(1) = hist_range(1) - step;
                    end
                    if isinf(finite_edges(end))
                        finite_edges(end) = hist_range(2) + step;
                    end
                    
                    histogram('BinEdges', finite_edges, 'BinCounts', prob, ...
                        'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
                    
                    hold on;
                    if strcmpi(metric, 'LLE')
                        xline(0, '--', 'Color', [0 0.7 0], 'LineWidth', 2);
                    end
                    hold off;
                end
                
                % Labels
                if condition_titles.isKey(cond_name)
                    title(condition_titles(cond_name), 'FontSize', 14);
                else
                    title(strrep(cond_name, '_', ' '), 'FontSize', 14);
                end
                
                if c_idx == 1
                    ylabel('Probability', 'FontSize', 12);
                end
                xlabel(y_label, 'FontSize', 12);
                box off;
            end
            
            % Link y-axes
            ax_handles = findobj(fig, 'Type', 'Axes');
            linkaxes(ax_handles, 'y');
            
            % Save figure
            fig_dir = fullfile(obj.output_dir, 'figures');
            if ~exist(fig_dir, 'dir')
                mkdir(fig_dir);
            end
            saveas(fig, fullfile(fig_dir, sprintf('%s_distribution.png', metric)));
            saveas(fig, fullfile(fig_dir, sprintf('%s_distribution.fig', metric)));
            fprintf('Figure saved to: %s\n', fig_dir);
        end
        
        function load_results(obj, results_dir)
            % LOAD_RESULTS Load results from a previous run
            %
            % Usage:
            %   psa.load_results('/path/to/param_space_...')
            
            obj.output_dir = results_dir;
            
            % Load summary
            summary_file = fullfile(results_dir, 'param_space_summary.mat');
            if exist(summary_file, 'file')
                loaded = load(summary_file);
                if isfield(loaded, 'summary_data')
                    obj.grid_params = loaded.summary_data.grid_params;
                    obj.param_ranges = loaded.summary_data.param_ranges;
                    obj.n_levels = loaded.summary_data.n_levels;
                    obj.conditions = loaded.summary_data.conditions;
                end
            end
            
            % Load per-condition results
            for c_idx = 1:length(obj.conditions)
                cond_name = obj.conditions{c_idx}.name;
                results_file = fullfile(results_dir, cond_name, ...
                    sprintf('param_space_results_%s.mat', cond_name));
                if exist(results_file, 'file')
                    loaded = load(results_file, 'results');
                    obj.results.(cond_name) = loaded.results;
                    fprintf('Loaded %d results for condition: %s\n', ...
                        length(loaded.results), cond_name);
                end
            end
            
            obj.has_run = true;
        end
    end
    
    %% Private Methods
    methods (Access = private)
        function set_default_conditions(obj)
            % SET_DEFAULT_CONDITIONS Initialize the four adaptation conditions
            
            obj.conditions = { ...
                struct('name', 'no_adaptation', 'n_a_E', 0, 'n_b_E', 0), ...
                struct('name', 'sfa_only',      'n_a_E', 3, 'n_b_E', 0), ...
                struct('name', 'std_only',      'n_a_E', 0, 'n_b_E', 1), ...
                struct('name', 'sfa_and_std',   'n_a_E', 3, 'n_b_E', 1) ...
            };
        end
        
        function generate_grid(obj)
            % GENERATE_GRID Create the multi-dimensional parameter grid
            
            n_params = length(obj.grid_params);
            obj.param_vectors = cell(1, n_params);
            
            for i = 1:n_params
                param_name = obj.grid_params{i};
                param_range = obj.param_ranges.(param_name);
                
                if ismember(param_name, obj.integer_params)
                    obj.param_vectors{i} = round(linspace(param_range(1), param_range(2), obj.n_levels));
                else
                    obj.param_vectors{i} = linspace(param_range(1), param_range(2), obj.n_levels);
                end
            end
            
            % Generate all combinations using ndgrid
            grid_cells = cell(size(obj.param_vectors));
            [grid_cells{:}] = ndgrid(obj.param_vectors{:});
            
            obj.num_combinations = numel(grid_cells{1});
            obj.all_configs = cell(obj.num_combinations, 1);
            
            for i = 1:obj.num_combinations
                config = struct();
                for j = 1:n_params
                    config.(obj.grid_params{j}) = grid_cells{j}(i);
                end
                obj.all_configs{i} = config;
            end
            
            % Randomize order for representative early-stopping
            rng('shuffle');
            obj.shuffled_indices = randperm(obj.num_combinations);
            
            fprintf('Generated %d parameter combinations (randomized order)\n', obj.num_combinations);
        end
        
        function run_batched_simulation(obj, temp_dir)
            % RUN_BATCHED_SIMULATION Execute simulations in batches with checkpoints
            
            num_batches = ceil(obj.num_combinations / obj.batch_size);
            conditions_local = obj.conditions;
            num_conditions = length(conditions_local);
            
            fprintf('Running %d combinations in %d batches...\n', obj.num_combinations, num_batches);
            
            for batch_idx = 1:num_batches
                batch_file = fullfile(temp_dir, sprintf('batch_%d.mat', batch_idx));
                
                % Skip if batch already completed (resume capability)
                if exist(batch_file, 'file')
                    fprintf('Batch %d/%d already completed. Skipping.\n', batch_idx, num_batches);
                    continue;
                end
                
                start_idx = (batch_idx - 1) * obj.batch_size + 1;
                end_idx = min(batch_idx * obj.batch_size, obj.num_combinations);
                batch_indices = obj.shuffled_indices(start_idx:end_idx);
                current_batch_size = length(batch_indices);
                
                fprintf('\n--- Batch %d/%d (configs %d-%d) ---\n', ...
                    batch_idx, num_batches, start_idx, end_idx);
                
                % Create jobs: each config runs ALL conditions with SAME network seed
                total_jobs = current_batch_size * num_conditions;
                jobs = cell(total_jobs, 1);
                job_idx = 1;
                
                for k = 1:current_batch_size
                    config_idx = batch_indices(k);
                    config = obj.all_configs{config_idx};
                    
                    % Same network seed for all conditions of this config
                    network_seed = config_idx;
                    
                    for c_idx = 1:num_conditions
                        job = struct();
                        job.config = config;
                        job.config_idx = config_idx;
                        job.condition = conditions_local{c_idx};
                        job.condition_idx = c_idx;
                        job.network_seed = network_seed;
                        jobs{job_idx} = job;
                        job_idx = job_idx + 1;
                    end
                end
                
                % Extract values for parfor
                model_defaults_local = obj.model_defaults;
                grid_params_local = obj.grid_params;
                verbose_local = obj.verbose;
                
                % Run parfor
                parallel_results = cell(total_jobs, 1);
                batch_start = tic;
                
                parfor j = 1:total_jobs
                    job = jobs{j};
                    run_start = tic;
                    
                    try
                        % Build model arguments
                        model_args = { ...
                            'n_a_E', job.condition.n_a_E, ...
                            'n_b_E', job.condition.n_b_E, ...
                            'rng_seeds', [job.network_seed, job.network_seed + 1] ...
                        };
                        
                        % Add grid parameters
                        for p_idx = 1:length(grid_params_local)
                            pname = grid_params_local{p_idx};
                            model_args = [model_args, {pname, job.config.(pname)}];
                        end
                        
                        % Add model defaults (don't override grid params or condition params)
                        default_fields = fieldnames(model_defaults_local);
                        for d_idx = 1:length(default_fields)
                            fname = default_fields{d_idx};
                            if ~ismember(fname, grid_params_local) && ...
                               ~strcmp(fname, 'n_a_E') && ~strcmp(fname, 'n_b_E')
                                model_args = [model_args, {fname, model_defaults_local.(fname)}];
                            end
                        end
                        
                        % Create and run model
                        model = SRNNModel(model_args{:});
                        model.build();
                        model.run();
                        
                        % Extract results
                        result = struct();
                        result.success = true;
                        result.config = job.config;
                        result.config_idx = job.config_idx;
                        result.condition_name = job.condition.name;
                        result.network_seed = job.network_seed;
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
                        
                        parallel_results{j} = result;
                        
                    catch ME
                        result = struct();
                        result.success = false;
                        result.error_message = ME.message;
                        result.config = job.config;
                        result.config_idx = job.config_idx;
                        result.condition_name = job.condition.name;
                        result.network_seed = job.network_seed;
                        result.run_duration = toc(run_start);
                        result.LLE = NaN;
                        result.mean_rate = NaN;
                        
                        parallel_results{j} = result;
                        
                        if verbose_local
                            fprintf('  ERROR job %d: %s\n', j, ME.message);
                        end
                    end
                end
                
                batch_elapsed = toc(batch_start);
                
                % Organize results by condition
                batch_results = struct();
                for c_idx = 1:num_conditions
                    cond_name = conditions_local{c_idx}.name;
                    batch_results.(cond_name) = {};
                end
                
                for j = 1:total_jobs
                    res = parallel_results{j};
                    cond_name = res.condition_name;
                    batch_results.(cond_name){end+1} = res;
                end
                
                % Save batch checkpoint
                save(batch_file, 'batch_results', 'batch_indices', '-v7.3');
                
                % Count successes
                n_success = sum(cellfun(@(r) r.success, parallel_results));
                fprintf('Batch %d completed in %.1f min (%d/%d successful)\n', ...
                    batch_idx, batch_elapsed/60, n_success, total_jobs);
            end
        end
        
        function consolidate_results(obj, temp_dir)
            % CONSOLIDATE_RESULTS Merge batch files into per-condition MAT files
            
            fprintf('\nConsolidating batch results...\n');
            
            num_batches = ceil(obj.num_combinations / obj.batch_size);
            
            % Initialize results storage
            for c_idx = 1:length(obj.conditions)
                cond_name = obj.conditions{c_idx}.name;
                obj.results.(cond_name) = cell(obj.num_combinations, 1);
            end
            
            % Load and merge batches
            all_found = true;
            for batch_idx = 1:num_batches
                batch_file = fullfile(temp_dir, sprintf('batch_%d.mat', batch_idx));
                
                if exist(batch_file, 'file')
                    loaded = load(batch_file);
                    batch_results = loaded.batch_results;
                    batch_indices = loaded.batch_indices;
                    
                    for c_idx = 1:length(obj.conditions)
                        cond_name = obj.conditions{c_idx}.name;
                        cond_results = batch_results.(cond_name);
                        
                        for k = 1:length(cond_results)
                            res = cond_results{k};
                            if isstruct(res) && isfield(res, 'config_idx')
                                obj.results.(cond_name){res.config_idx} = res;
                            end
                        end
                    end
                else
                    fprintf('Warning: Batch file %d not found\n', batch_idx);
                    all_found = false;
                end
            end
            
            % Save per-condition results
            for c_idx = 1:length(obj.conditions)
                cond_name = obj.conditions{c_idx}.name;
                results = obj.results.(cond_name);
                
                save_file = fullfile(obj.output_dir, cond_name, ...
                    sprintf('param_space_results_%s.mat', cond_name));
                save(save_file, 'results', '-v7.3');
                
                % Count successes
                n_success = sum(cellfun(@(r) isstruct(r) && isfield(r, 'success') && r.success, results));
                fprintf('Condition %s: %d/%d successful, saved to %s\n', ...
                    cond_name, n_success, obj.num_combinations, save_file);
            end
            
            % Clean up temp directory if all successful
            if all_found
                rmdir(temp_dir, 's');
                fprintf('Temp directory cleaned up.\n');
            else
                fprintf('Temp directory retained due to missing batches.\n');
            end
        end
        
        function save_summary(obj)
            % SAVE_SUMMARY Save analysis summary to disk
            
            summary_file = fullfile(obj.output_dir, 'param_space_summary.mat');
            
            summary_data = struct();
            summary_data.grid_params = obj.grid_params;
            summary_data.param_ranges = obj.param_ranges;
            summary_data.param_vectors = obj.param_vectors;
            summary_data.n_levels = obj.n_levels;
            summary_data.num_combinations = obj.num_combinations;
            summary_data.conditions = obj.conditions;
            summary_data.model_defaults = obj.model_defaults;
            summary_data.shuffled_indices = obj.shuffled_indices;
            summary_data.analysis_start_time = obj.analysis_start_time;
            summary_data.analysis_completed = datestr(now);
            
            % Compute statistics per condition
            for c_idx = 1:length(obj.conditions)
                cond_name = obj.conditions{c_idx}.name;
                if isfield(obj.results, cond_name)
                    results = obj.results.(cond_name);
                    n_success = sum(cellfun(@(r) isstruct(r) && isfield(r, 'success') && r.success, results));
                    summary_data.stats.(cond_name).n_success = n_success;
                    summary_data.stats.(cond_name).n_total = length(results);
                    summary_data.stats.(cond_name).success_rate = n_success / length(results);
                end
            end
            
            save(summary_file, 'summary_data', '-v7.3');
            fprintf('Summary saved to: %s\n', summary_file);
        end
    end
end
