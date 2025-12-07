classdef PairedPulseMIAnalysis < handle
    % PAIREDPULSEMIANALYSIS Paired-pulse mutual information analysis
    %
    % This class performs paired-pulse stimulation experiments to measure
    % mutual information (MI) between stimulus amplitude and neural response
    % across the four adaptation conditions. The same network connectivity
    % is used across all conditions for fair comparison.
    %
    % The paired-pulse protocol:
    %   - Pulse 1: Variable amplitude (randomized across trials)
    %   - Pulse 2: Fixed amplitude (probe pulse)
    %   - MI measures how well the response during pulse 2 encodes the
    %     amplitude of pulse 1 (memory/history dependence)
    %
    % Usage:
    %   pp = PairedPulseMIAnalysis('n_networks', 10, 'note', 'demo');
    %   pp.run();
    %   pp.plot();
    %
    % See also: SRNNModel, ParamSpaceAnalysis
    
    %% Stimulus Configuration Properties
    properties
        % Paired-pulse timing
        pulse1_width = 2            % Width of first pulse (s)
        pulse2_width = 2            % Width of second pulse (s)
        inter_pulse_interval = 3    % Time from pulse 1 start to pulse 2 start (s)
        repeat_interval = 15        % Time between pair starts (s)
        
        % Pulse amplitudes
        pulse2_amp = 2.5            % Fixed amplitude of second pulse
        n_amp_levels = 5            % Number of discrete amplitude levels for pulse 1
        amp_scale = 0.5             % Scaling factor for pulse 1 amplitudes
        
        % Baseline
        DC_level = 0.1              % DC baseline input level
        ramp_duration = 10          % Ramp-up duration to DC (s)
        
        % Stimulus channel
        stim_channel = 1            % Which neuron receives the stimulus
    end
    
    %% Simulation Configuration Properties
    properties
        n_networks = 20             % Number of different network realizations
        conditions                  % Cell array of condition structs
        model_defaults = struct()   % Default SRNNModel properties
        verbose = true              % Print progress during execution
    end
    
    %% MI Analysis Properties
    properties
        delay_step_samples = 10     % Step size for delay vector (samples)
        n_bins_in = 5               % Number of bins for input discretization
        n_bins_out = 5              % Number of bins for output discretization
    end
    
    %% Output Properties
    properties
        output_dir                  % Base directory for saving results
        note = ''                   % Optional note for folder naming
    end
    
    %% Grid Parameter Properties (for parameter space exploration)
    properties
        grid_params = {}            % Cell array of parameter names for grid search
        param_ranges = struct()     % Struct with min/max ranges for each grid param
        n_levels = 5                % Number of levels per grid parameter
        batch_size = 10             % Networks per batch for checkpointing
    end
    
    %% Results (SetAccess = private)
    properties (SetAccess = private)
        results = struct()          % Stored results after run
        has_run = false             % Flag indicating if analysis has run
        analysis_start_time         % Timestamp when analysis started
        all_configs = {}            % All generated parameter configurations
        shuffled_indices = []       % Randomized execution order
    end
    
    %% Constructor
    methods
        function obj = PairedPulseMIAnalysis(varargin)
            % PAIREDPULSEMIANALYSIS Constructor with name-value pairs
            
            % Set default conditions
            obj.set_default_conditions();
            
            % Set model defaults
            obj.set_model_defaults();
            
            % Parse name-value pairs
            for i = 1:2:length(varargin)
                if isprop(obj, varargin{i})
                    obj.(varargin{i}) = varargin{i+1};
                else
                    warning('PairedPulseMIAnalysis:UnknownProperty', ...
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
        function set_conditions(obj, conditions_cell)
            % SET_CONDITIONS Set custom conditions
            obj.conditions = conditions_cell;
            if obj.verbose
                fprintf('Set %d custom conditions\n', length(conditions_cell));
            end
        end
        
        function add_grid_parameter(obj, param_name, param_range)
            % ADD_GRID_PARAMETER Add a parameter to the grid search
            %
            % Usage:
            %   pp.add_grid_parameter('level_of_chaos', [0.5, 2.5]);
            %   pp.add_grid_parameter('n', [50, 200]);
            %
            % All combinations of grid parameters will be tested.
            
            if ~ismember(param_name, obj.grid_params)
                obj.grid_params{end+1} = param_name;
            end
            obj.param_ranges.(param_name) = param_range;
            
            if obj.verbose
                fprintf('Added grid parameter: %s = [%.2f, %.2f] (%d levels)\n', ...
                    param_name, param_range(1), param_range(2), obj.n_levels);
            end
        end
        
        function run(obj)
            % RUN Execute the paired-pulse MI analysis
            
            if isempty(obj.conditions)
                error('PairedPulseMIAnalysis:NoConditions', 'No conditions defined.');
            end
            
            % Create timestamped output directory
            obj.analysis_start_time = datetime('now');
            dt_str = lower(strrep(datestr(obj.analysis_start_time, 'mmm_dd_yy_HH_MM_AM'), ':', '_'));
            
            if ~isempty(obj.note)
                folder_name = sprintf('paired_pulse_MI_%s_nNets_%d_%s', ...
                    obj.note, obj.n_networks, dt_str);
            else
                folder_name = sprintf('paired_pulse_MI_nNets_%d_%s', ...
                    obj.n_networks, dt_str);
            end
            
            obj.output_dir = fullfile(obj.output_dir, folder_name);
            if ~exist(obj.output_dir, 'dir')
                mkdir(obj.output_dir);
            end
            
            % Print summary
            fprintf('\n========================================\n');
            fprintf('=== Paired-Pulse MI Analysis ===\n');
            fprintf('========================================\n');
            fprintf('Number of networks: %d\n', obj.n_networks);
            fprintf('Amplitude levels: %d\n', obj.n_amp_levels);
            fprintf('Conditions: %s\n', strjoin(cellfun(@(c) c.name, obj.conditions, 'UniformOutput', false), ', '));
            fprintf('Output directory: %s\n', obj.output_dir);
            fprintf('========================================\n\n');
            
            overall_start = tic;
            
            % Initialize results
            n_conditions = length(obj.conditions);
            for c_idx = 1:n_conditions
                cond_name = obj.conditions{c_idx}.name;
                obj.results.(cond_name) = struct();
                obj.results.(cond_name).MI_vs_delay_all = {};
                obj.results.(cond_name).delay_vec_sec = [];
                obj.results.(cond_name).network_results = {};
            end
            
            % Run simulations
            for net_idx = 1:obj.n_networks
                fprintf('\n--- Network %d/%d ---\n', net_idx, obj.n_networks);
                
                % Same network seed for all conditions
                network_seed = net_idx;
                
                for c_idx = 1:n_conditions
                    condition = obj.conditions{c_idx};
                    
                    try
                        result = obj.run_single_network(network_seed, condition);
                        result.success = true;
                    catch ME
                        result = struct();
                        result.success = false;
                        result.error_message = ME.message;
                        result.MI_vs_delay = [];
                        result.delay_vec_sec = [];
                        if obj.verbose
                            fprintf('  ERROR (%s): %s\n', condition.name, ME.message);
                        end
                    end
                    
                    result.network_seed = network_seed;
                    result.condition = condition;
                    
                    % Store results
                    cond_name = condition.name;
                    obj.results.(cond_name).network_results{end+1} = result;
                    if result.success && ~isempty(result.MI_vs_delay)
                        obj.results.(cond_name).MI_vs_delay_all{end+1} = result.MI_vs_delay;
                        if isempty(obj.results.(cond_name).delay_vec_sec)
                            obj.results.(cond_name).delay_vec_sec = result.delay_vec_sec;
                        end
                    end
                end
            end
            
            overall_elapsed = toc(overall_start);
            fprintf('\n========================================\n');
            fprintf('=== Analysis Complete ===\n');
            fprintf('Total time: %.2f minutes\n', overall_elapsed/60);
            fprintf('========================================\n');
            
            obj.has_run = true;
            
            % Save results
            obj.save_results();
        end
        
        function plot(obj, varargin)
            % PLOT Create comparison plot of MI vs delay across conditions
            
            if ~obj.has_run && isempty(fieldnames(obj.results))
                error('PairedPulseMIAnalysis:NotRun', 'Analysis has not been run yet.');
            end
            
            % Parse optional arguments
            show_individual = false;
            for i = 1:2:length(varargin)
                if strcmpi(varargin{i}, 'show_individual')
                    show_individual = varargin{i+1};
                end
            end
            
            % Colors for conditions
            colors = containers.Map(...
                {'no_adaptation', 'sfa_only', 'std_only', 'sfa_and_std'}, ...
                {[0 0 1], [1 0 0], [0 0.7 0], [0.7 0 0.7]});
            
            % Readable titles
            condition_titles = containers.Map(...
                {'no_adaptation', 'sfa_only', 'std_only', 'sfa_and_std'}, ...
                {'No Adaptation', 'SFA Only', 'STD Only', 'SFA + STD'});
            
            % Create figure
            fig = figure('Name', 'Paired-Pulse MI Comparison', 'Position', [100, 100, 800, 600]);
            hold on;
            
            condition_names = fieldnames(obj.results);
            legend_entries = {};
            
            for c_idx = 1:length(condition_names)
                cond_name = condition_names{c_idx};
                MI_all = obj.results.(cond_name).MI_vs_delay_all;
                delay_vec = obj.results.(cond_name).delay_vec_sec;
                
                if isempty(MI_all) || isempty(delay_vec)
                    continue;
                end
                
                % Convert to matrix
                MI_matrix = vertcat(MI_all{:});
                
                % Compute mean and std
                MI_mean = mean(MI_matrix, 1, 'omitnan');
                MI_std = std(MI_matrix, 0, 1, 'omitnan');
                
                % Get color
                if colors.isKey(cond_name)
                    color = colors(cond_name);
                else
                    color = [0.5 0.5 0.5];
                end
                
                % Plot individual traces (optional)
                if show_individual
                    for k = 1:size(MI_matrix, 1)
                        plot(delay_vec, MI_matrix(k, :), '-', 'Color', [color 0.2], 'LineWidth', 0.5);
                    end
                end
                
                % Plot shaded error region
                fill([delay_vec, fliplr(delay_vec)], ...
                     [MI_mean + MI_std, fliplr(MI_mean - MI_std)], ...
                     color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                
                % Plot mean
                if condition_titles.isKey(cond_name)
                    disp_name = condition_titles(cond_name);
                else
                    disp_name = strrep(cond_name, '_', ' ');
                end
                plot(delay_vec, MI_mean, '-', 'Color', color, 'LineWidth', 2.5, ...
                    'DisplayName', disp_name);
                legend_entries{end+1} = disp_name;
            end
            
            xlabel('Delay from Pulse 2 Start (s)', 'FontSize', 14);
            ylabel('Mutual Information (bits)', 'FontSize', 14);
            title('Paired-Pulse Mutual Information', 'FontSize', 16);
            legend(legend_entries, 'Location', 'best');
            grid on;
            box on;
            set(gca, 'FontSize', 12);
            xlim([0, obj.pulse2_width]);
            
            hold off;
            
            % Save figure
            fig_dir = fullfile(obj.output_dir, 'figures');
            if ~exist(fig_dir, 'dir')
                mkdir(fig_dir);
            end
            saveas(fig, fullfile(fig_dir, 'MI_vs_delay_comparison.png'));
            saveas(fig, fullfile(fig_dir, 'MI_vs_delay_comparison.fig'));
            fprintf('Figure saved to: %s\n', fig_dir);
            
            % Print summary
            fprintf('\n=== MI Summary ===\n');
            for c_idx = 1:length(condition_names)
                cond_name = condition_names{c_idx};
                MI_all = obj.results.(cond_name).MI_vs_delay_all;
                if ~isempty(MI_all)
                    MI_matrix = vertcat(MI_all{:});
                    max_MI = max(mean(MI_matrix, 1, 'omitnan'));
                    fprintf('%s: Max Mean MI = %.4f bits\n', cond_name, max_MI);
                end
            end
        end
        
        function load_results(obj, results_dir)
            % LOAD_RESULTS Load results from a previous run
            
            obj.output_dir = results_dir;
            results_file = fullfile(results_dir, 'paired_pulse_MI_results.mat');
            if exist(results_file, 'file')
                loaded = load(results_file);
                obj.results = loaded.results;
                if isfield(loaded, 'config')
                    obj.n_networks = loaded.config.n_networks;
                    obj.conditions = loaded.config.conditions;
                end
                obj.has_run = true;
                fprintf('Loaded results from: %s\n', results_file);
            else
                error('Results file not found: %s', results_file);
            end
        end
        
        function plot_debug(obj, network_seed, condition_name)
            % PLOT_DEBUG Run a single simulation and visualize stimulus + response
            %
            % Usage:
            %   pp.plot_debug(1, 'no_adaptation')
            %
            % This helps diagnose why MI might be zero by showing:
            %   1. The external stimulus (paired pulses)
            %   2. The firing rate of the stimulated neuron
            %   3. The stimulus-response relationship
            
            if nargin < 3
                condition_name = 'no_adaptation';
            end
            if nargin < 2
                network_seed = 1;
            end
            
            % Find matching condition
            condition = [];
            for c_idx = 1:length(obj.conditions)
                if strcmp(obj.conditions{c_idx}.name, condition_name)
                    condition = obj.conditions{c_idx};
                    break;
                end
            end
            if isempty(condition)
                error('Condition not found: %s', condition_name);
            end
            
            fprintf('Running debug simulation for %s (seed=%d)...\n', condition_name, network_seed);
            
            % Build model arguments
            model_args = { ...
                'n_a_E', condition.n_a_E, ...
                'n_b_E', condition.n_b_E, ...
                'rng_seeds', [network_seed, network_seed + 1] ...
            };
            
            % Add model defaults
            default_fields = fieldnames(obj.model_defaults);
            for d_idx = 1:length(default_fields)
                fname = default_fields{d_idx};
                if ~strcmp(fname, 'n_a_E') && ~strcmp(fname, 'n_b_E')
                    model_args = [model_args, {fname, obj.model_defaults.(fname)}];
                end
            end
            
            % Create custom input configuration for paired-pulse
            input_config = obj.create_paired_pulse_config(network_seed);
            model_args = [model_args, {'input_config', input_config}];
            
            % Create and run model
            model = SRNNModel(model_args{:});
            model.build();
            model.run();
            
            % Extract data (use x instead of r to avoid saturation)
            t = model.plot_data.t;
            x_all = [model.plot_data.x.E; model.plot_data.x.I];
            x_stim = x_all(obj.stim_channel, :);
            
            % Get external input for the stim channel
            u_stim = model.plot_data.u.E(obj.stim_channel, :);
            
            % Create figure with 3 subplots
            figure('Name', sprintf('Debug: %s (seed=%d)', condition_name, network_seed), ...
                'Position', [100, 100, 1200, 800]);
            
            % Subplot 1: Stimulus
            subplot(3, 1, 1);
            plot(t, u_stim, 'b-', 'LineWidth', 1);
            xlabel('Time (s)');
            ylabel('Stimulus');
            title(sprintf('External Input to Channel %d', obj.stim_channel));
            xlim([0, min(60, max(t))]);  % Show first 60 seconds
            grid on;
            
            % Subplot 2: Dendritic potential
            subplot(3, 1, 2);
            plot(t, x_stim, 'r-', 'LineWidth', 1);
            xlabel('Time (s)');
            ylabel('Dendritic Potential');
            title(sprintf('Dendritic Potential of Channel %d', obj.stim_channel));
            xlim([0, min(60, max(t))]);
            grid on;
            
            % Subplot 3: Stimulus-response relationship at pulse 2
            subplot(3, 1, 3);
            pulse2_start_times = input_config.pair_start_times + input_config.inter_pulse_interval;
            pulse1_amps = input_config.pulse1_amps;
            
            % Get firing rate at middle of pulse 2
            mid_delay_sec = input_config.pulse2_width / 2;
            fs_deci = obj.model_defaults.fs / model.plot_deci;
            
            x_at_pulse2 = zeros(size(pulse2_start_times));
            for i = 1:length(pulse2_start_times)
                [~, idx] = min(abs(t - (pulse2_start_times(i) + mid_delay_sec)));
                if idx <= length(x_stim)
                    x_at_pulse2(i) = x_stim(idx);
                else
                    x_at_pulse2(i) = NaN;
                end
            end
            
            scatter(pulse1_amps, x_at_pulse2, 100, 'filled', 'MarkerFaceAlpha', 0.7);
            xlabel('Pulse 1 Amplitude');
            ylabel('Dendritic Potential at Pulse 2 Midpoint');
            title('Stimulus-Response Relationship');
            grid on;
            
            % Add correlation
            valid_idx = ~isnan(x_at_pulse2);
            if sum(valid_idx) > 2
                r_corr = corrcoef(pulse1_amps(valid_idx), x_at_pulse2(valid_idx));
                title(sprintf('Stimulus-Response Relationship (r = %.3f)', r_corr(1,2)));
            end
            
            % Print diagnostic info
            fprintf('\n=== Debug Info ===\n');
            fprintf('Simulation time: [%.1f, %.1f] s\n', obj.model_defaults.T_range(1), obj.model_defaults.T_range(2));
            fprintf('Number of pulse pairs: %d\n', length(pulse2_start_times));
            fprintf('Pulse 1 amplitudes: %s\n', mat2str(pulse1_amps, 2));
            fprintf('Dendritic potential at pulse 2: %s\n', mat2str(x_at_pulse2, 3));
            fprintf('Stimulus range: [%.2f, %.2f]\n', min(u_stim), max(u_stim));
            fprintf('Dendritic potential range: [%.2f, %.2f]\n', min(x_stim), max(x_stim));
            fprintf('==================\n');
        end
        
        function plot_histograms(obj, varargin)
            % PLOT_HISTOGRAMS Plot LLE and mean rate distributions per condition
            %
            % Usage:
            %   pp.plot_histograms();
            %   pp.plot_histograms('lle_range', [-1.5, 1.5], 'rate_range', [0, 5]);
            
            if ~obj.has_run
                error('Analysis has not been run yet.');
            end
            
            % Parse options
            p = inputParser;
            addParameter(p, 'lle_range', [-1.5, 1.5]);
            addParameter(p, 'rate_range', [0, 5]);
            addParameter(p, 'n_bins', 25);
            parse(p, varargin{:});
            
            lle_range = p.Results.lle_range;
            rate_range = p.Results.rate_range;
            n_bins = p.Results.n_bins;
            
            condition_names = fieldnames(obj.results);
            n_conditions = length(condition_names);
            
            condition_titles = containers.Map(...
                {'no_adaptation', 'sfa_only', 'std_only', 'sfa_and_std'}, ...
                {'No Adaptation', 'SFA Only', 'STD Only', 'SFA + STD'});
            
            fig = figure('Name', 'LLE and Rate Distributions', 'Position', [100, 100, 650, 1020]);
            
            ax_lle = gobjects(n_conditions, 1);
            ax_rate = gobjects(n_conditions, 1);
            
            for i = 1:n_conditions
                cond_name = condition_names{i};
                results_list = obj.results.(cond_name).network_results;
                
                % Extract LLE and mean_rate values
                lles = [];
                rates = [];
                for k = 1:length(results_list)
                    r = results_list{k};
                    if isfield(r, 'LLE') && ~isempty(r.LLE)
                        lles(end+1) = r.LLE;
                    end
                    if isfield(r, 'mean_rate') && ~isempty(r.mean_rate)
                        rates(end+1) = r.mean_rate;
                    end
                end
                
                % LLE histogram
                ax_lle(i) = subplot(n_conditions, 2, 2*i-1);
                histogram(lles, linspace(lle_range(1), lle_range(2), n_bins), ...
                    'Normalization', 'probability', 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
                if i == n_conditions, xlabel('LLE (\lambda_1)', 'FontSize', 14); end
                xlim(lle_range);
                box off;
                
                if condition_titles.isKey(cond_name)
                    title_str = condition_titles(cond_name);
                else
                    title_str = strrep(cond_name, '_', ' ');
                end
                text(-0.22, 0.5, title_str, 'Units', 'normalized', 'Rotation', 90, ...
                    'FontSize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                
                % Rate histogram
                ax_rate(i) = subplot(n_conditions, 2, 2*i);
                histogram(rates, linspace(rate_range(1), rate_range(2), n_bins), ...
                    'Normalization', 'probability', 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
                if i == n_conditions, xlabel('Mean Firing Rate (Hz)', 'FontSize', 14); end
                xlim(rate_range);
                box off;
            end
            
            linkaxes([ax_lle; ax_rate], 'y');
            
            % Save figure
            fig_dir = fullfile(obj.output_dir, 'figures');
            if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
            saveas(fig, fullfile(fig_dir, 'LLE_rate_distributions.png'));
            saveas(fig, fullfile(fig_dir, 'LLE_rate_distributions.fig'));
            fprintf('Histograms saved to: %s\n', fig_dir);
        end
        
        function plot_mi_imagesc(obj, varargin)
            % PLOT_MI_IMAGESC Create LLE vs delay heatmap of mean MI
            %
            % Usage:
            %   pp.plot_mi_imagesc();
            %   pp.plot_mi_imagesc('lle_range', [-2, 1], 'n_lle_bins', 20);
            
            if ~obj.has_run
                error('Analysis has not been run yet.');
            end
            
            % Parse options
            p = inputParser;
            addParameter(p, 'lle_range', [-2, 1]);
            addParameter(p, 'n_lle_bins', 15);
            parse(p, varargin{:});
            
            lle_range = p.Results.lle_range;
            n_lle_bins = p.Results.n_lle_bins;
            
            % Pool data across all conditions
            condition_names = fieldnames(obj.results);
            pooled_mi = [];
            pooled_lle = [];
            template_delay_vec = [];
            
            for c_idx = 1:length(condition_names)
                cond_name = condition_names{c_idx};
                results_list = obj.results.(cond_name).network_results;
                delay_vec = obj.results.(cond_name).delay_vec_sec;
                
                if isempty(template_delay_vec) && ~isempty(delay_vec)
                    template_delay_vec = delay_vec;
                end
                
                for k = 1:length(results_list)
                    r = results_list{k};
                    if ~isfield(r, 'success') || ~r.success
                        continue;
                    end
                    if isfield(r, 'MI_vs_delay') && isfield(r, 'LLE')
                        mi_vec = r.MI_vs_delay(:)';
                        lle_val = r.LLE;
                        if ~isempty(mi_vec) && length(mi_vec) == length(template_delay_vec)
                            pooled_mi(end+1, :) = mi_vec;
                            pooled_lle(end+1, 1) = lle_val;
                        end
                    end
                end
            end
            
            if isempty(pooled_mi) || isempty(template_delay_vec)
                warning('No MI data available for imagesc plot.');
                return;
            end
            
            % Create binned heatmap
            lle_edges = linspace(lle_range(1), lle_range(2), n_lle_bins + 1);
            lle_centers = (lle_edges(1:end-1) + lle_edges(2:end)) / 2;
            n_delays = length(template_delay_vec);
            mi_img = NaN(n_lle_bins, n_delays);
            
            for b = 1:n_lle_bins
                in_bin = pooled_lle >= lle_edges(b) & pooled_lle < lle_edges(b+1);
                if any(in_bin)
                    mi_img(b, :) = mean(pooled_mi(in_bin, :), 1, 'omitnan');
                end
            end
            
            fig = figure('Name', 'MI vs Delay by LLE', 'Position', [643, 500, 404, 364]);
            imagesc(template_delay_vec, lle_centers, mi_img);
            axis xy;
            colormap(hot);
            xlabel('Delay (s)', 'FontSize', 14);
            ylabel('LLE', 'FontSize', 14);
            title('Mean MI vs Delay and LLE');
            colorbar;
            ylabel(colorbar, 'MI (bits)');
            
            % Save
            fig_dir = fullfile(obj.output_dir, 'figures');
            if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
            saveas(fig, fullfile(fig_dir, 'MI_vs_delay_LLE_imagesc.png'));
            saveas(fig, fullfile(fig_dir, 'MI_vs_delay_LLE_imagesc.fig'));
            fprintf('Imagesc plot saved to: %s\n', fig_dir);
        end
        
        function plot_mi_vs_lle(obj, varargin)
            % PLOT_MI_VS_LLE Plot LLE vs mean MI with bootstrap CIs
            %
            % Usage:
            %   pp.plot_mi_vs_lle();
            %   pp.plot_mi_vs_lle('delay_window', [0.5, 1.5], 'n_bootstrap', 1000);
            
            if ~obj.has_run
                error('Analysis has not been run yet.');
            end
            
            % Parse options
            p = inputParser;
            addParameter(p, 'delay_window', [0.5, 1.5]);  % seconds
            addParameter(p, 'lle_range', [-2, 1]);
            addParameter(p, 'n_lle_bins', 12);
            addParameter(p, 'n_bootstrap', 1000);
            parse(p, varargin{:});
            
            delay_window = p.Results.delay_window;
            lle_range = p.Results.lle_range;
            n_lle_bins = p.Results.n_lle_bins;
            n_bootstrap = p.Results.n_bootstrap;
            
            % Pool data
            condition_names = fieldnames(obj.results);
            pooled_mi_avg = [];
            pooled_lle = [];
            
            for c_idx = 1:length(condition_names)
                cond_name = condition_names{c_idx};
                results_list = obj.results.(cond_name).network_results;
                delay_vec = obj.results.(cond_name).delay_vec_sec;
                
                if isempty(delay_vec)
                    continue;
                end
                
                delay_indices = find(delay_vec >= delay_window(1) & delay_vec <= delay_window(2));
                
                for k = 1:length(results_list)
                    r = results_list{k};
                    if ~isfield(r, 'success') || ~r.success
                        continue;
                    end
                    if isfield(r, 'MI_vs_delay') && isfield(r, 'LLE')
                        mi_vec = r.MI_vs_delay;
                        if ~isempty(mi_vec) && ~isempty(delay_indices)
                            mi_avg = mean(mi_vec(delay_indices), 'omitnan');
                            pooled_mi_avg(end+1) = mi_avg;
                            pooled_lle(end+1) = r.LLE;
                        end
                    end
                end
            end
            
            if isempty(pooled_mi_avg)
                warning('No MI data available for LLE vs MI plot.');
                return;
            end
            
            % Bin by LLE and compute mean + bootstrap CI
            lle_edges = linspace(lle_range(1), lle_range(2), n_lle_bins + 1);
            lle_centers = (lle_edges(1:end-1) + lle_edges(2:end)) / 2;
            
            mi_mean = NaN(1, n_lle_bins);
            mi_ci_lower = NaN(1, n_lle_bins);
            mi_ci_upper = NaN(1, n_lle_bins);
            
            for b = 1:n_lle_bins
                in_bin = pooled_lle >= lle_edges(b) & pooled_lle < lle_edges(b+1);
                mi_in_bin = pooled_mi_avg(in_bin);
                if length(mi_in_bin) >= 3
                    mi_mean(b) = mean(mi_in_bin, 'omitnan');
                    try
                        ci = bootci(n_bootstrap, @mean, mi_in_bin);
                        mi_ci_lower(b) = ci(1);
                        mi_ci_upper(b) = ci(2);
                    catch
                        % If bootstrap fails, just use std
                        mi_ci_lower(b) = mi_mean(b) - std(mi_in_bin);
                        mi_ci_upper(b) = mi_mean(b) + std(mi_in_bin);
                    end
                end
            end
            
            fig = figure('Name', 'LLE vs MI', 'Position', [700, 400, 380, 300]);
            y_neg = mi_mean - mi_ci_lower;
            y_pos = mi_ci_upper - mi_mean;
            errorbar(lle_centers, mi_mean, y_neg, y_pos, 'ko-', 'LineWidth', 1.5, 'CapSize', 4);
            xlabel('LLE (\lambda_1)', 'FontSize', 14);
            ylabel('Mean MI (bits)', 'FontSize', 14);
            title(sprintf('MI (delay %.1f-%.1f s) vs LLE', delay_window(1), delay_window(2)));
            box off;
            grid on;
            
            % Save
            fig_dir = fullfile(obj.output_dir, 'figures');
            if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
            saveas(fig, fullfile(fig_dir, 'MI_vs_LLE_slice.png'));
            saveas(fig, fullfile(fig_dir, 'MI_vs_LLE_slice.fig'));
            fprintf('LLE vs MI plot saved to: %s\n', fig_dir);
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
        
        function set_model_defaults(obj)
            % SET_MODEL_DEFAULTS Set default model parameters
            obj.model_defaults.n = 50;
            obj.model_defaults.fs = 200;
            obj.model_defaults.T_range = [0, 300];  % 5 minutes
            obj.model_defaults.level_of_chaos = 1.5;
            obj.model_defaults.lya_method = 'benettin';  % Enable LLE for plots
        end
        
        function result = run_single_network(obj, network_seed, condition, config)
            % RUN_SINGLE_NETWORK Run simulation for one network and condition
            %
            % Arguments:
            %   network_seed - Seed for network W matrix
            %   condition - Struct with name, n_a_E, n_b_E
            %   config (optional) - Struct with grid parameter values
            
            if nargin < 4
                config = struct();
            end
            
            if obj.verbose
                fprintf('  Running %s... ', condition.name);
            end
            
            run_start = tic;
            
            % Build model arguments
            model_args = { ...
                'n_a_E', condition.n_a_E, ...
                'n_b_E', condition.n_b_E, ...
                'rng_seeds', [network_seed, network_seed + 1] ...
            };
            
            % Add model defaults
            default_fields = fieldnames(obj.model_defaults);
            for d_idx = 1:length(default_fields)
                fname = default_fields{d_idx};
                if ~strcmp(fname, 'n_a_E') && ~strcmp(fname, 'n_b_E')
                    model_args = [model_args, {fname, obj.model_defaults.(fname)}];
                end
            end
            
            % Apply grid config parameters (override defaults)
            config_fields = fieldnames(config);
            for c_idx = 1:length(config_fields)
                fname = config_fields{c_idx};
                model_args = [model_args, {fname, config.(fname)}];
            end
            
            % Create custom input configuration for paired-pulse
            input_config = obj.create_paired_pulse_config(network_seed);
            model_args = [model_args, {'input_config', input_config}];
            
            % Create and run model
            model = SRNNModel(model_args{:});
            model.build();
            model.run();
            
            % Extract dendritic potential (use x instead of r to avoid saturation)
            if ~isempty(model.plot_data) && isfield(model.plot_data, 'x')
                t = model.plot_data.t;
                x_all = [model.plot_data.x.E; model.plot_data.x.I];
                x_stim = x_all(obj.stim_channel, :);  % Response of stimulated neuron
                fs = obj.model_defaults.fs / model.plot_deci;  % Decimated sampling rate
            else
                error('No dendritic potential data available');
            end
            
            % Compute MI vs delay
            [MI_vs_delay, delay_vec_sec] = obj.compute_MI_vs_delay(...
                input_config, x_stim, t, fs);
            
            % Extract LLE if computed
            LLE_val = NaN;
            if ~isempty(model.lya_results) && isfield(model.lya_results, 'LLE')
                LLE_val = model.lya_results.LLE;
            end
            
            % Compute mean firing rate
            r_all = [model.plot_data.r.E; model.plot_data.r.I];
            mean_rate = mean(r_all(:), 'omitnan');
            
            result = struct();
            result.MI_vs_delay = MI_vs_delay;
            result.delay_vec_sec = delay_vec_sec;
            result.LLE = LLE_val;
            result.mean_rate = mean_rate;
            result.config = config;
            result.run_duration = toc(run_start);
            
            if obj.verbose
                fprintf('done (%.1fs, LLE=%.3f, max MI=%.4f)\n', ...
                    result.run_duration, LLE_val, max(MI_vs_delay));
            end
        end
        
        function input_config = create_paired_pulse_config(obj, seed)
            % CREATE_PAIRED_PULSE_CONFIG Create paired-pulse stimulus configuration
            
            rng(seed + 1000, 'twister');  % Separate seed for stimulus
            
            T_end = obj.model_defaults.T_range(2);
            
            % Determine pair start times
            pair_start_times = 0:obj.repeat_interval:T_end;
            pair_start_times = pair_start_times(...
                pair_start_times + obj.inter_pulse_interval + obj.pulse2_width <= T_end);
            
            n_pairs = length(pair_start_times);
            
            % Generate random amplitudes for pulse 1
            pulse1_amps = obj.amp_scale * randi([1, obj.n_amp_levels], 1, n_pairs);
            
            % Create config struct
            input_config = struct();
            input_config.type = 'paired_pulse';
            input_config.pair_start_times = pair_start_times;
            input_config.pulse1_amps = pulse1_amps;
            input_config.pulse2_amp = obj.pulse2_amp;
            input_config.pulse1_width = obj.pulse1_width;
            input_config.pulse2_width = obj.pulse2_width;
            input_config.inter_pulse_interval = obj.inter_pulse_interval;
            input_config.stim_channel = obj.stim_channel;
            input_config.DC_level = obj.DC_level;
            input_config.ramp_duration = obj.ramp_duration;
            input_config.intrinsic_drive = [];
            
            % Custom generator function (rng_seed argument is passed but not used)
            input_config.generator = @(params, T, fs, rng_seed_unused, cfg) ...
                obj.generate_paired_pulse_stimulus(params, T, fs, cfg);
        end
        
        function [u_ex, t] = generate_paired_pulse_stimulus(obj, params, T_end, fs, cfg)
            % GENERATE_PAIRED_PULSE_STIMULUS Generate the stimulus time series
            
            n = params.n;
            dt = 1 / fs;
            t = (0:dt:T_end)';
            nt = length(t);
            
            u_ex = zeros(n, nt);
            
            % Add paired pulses
            for i = 1:length(cfg.pair_start_times)
                % Pulse 1
                p1_start = cfg.pair_start_times(i);
                p1_end = p1_start + cfg.pulse1_width;
                p1_indices = t >= p1_start & t < p1_end;
                u_ex(cfg.stim_channel, p1_indices) = cfg.pulse1_amps(i);
                
                % Pulse 2
                p2_start = p1_start + cfg.inter_pulse_interval;
                p2_end = p2_start + cfg.pulse2_width;
                p2_indices = t >= p2_start & t < p2_end;
                u_ex(cfg.stim_channel, p2_indices) = cfg.pulse2_amp;
            end
            
            % Add DC with ramp
            ramp_end_time = cfg.ramp_duration;
            ramp_indices = t <= ramp_end_time;
            num_ramp_points = sum(ramp_indices);
            
            dc_profile = ones(1, nt) * cfg.DC_level;
            dc_profile(ramp_indices) = linspace(0, cfg.DC_level, num_ramp_points);
            
            u_ex = u_ex + dc_profile;
        end
        
        function [MI_vs_delay, delay_vec_sec] = compute_MI_vs_delay(obj, cfg, r_stim, t, fs)
            % COMPUTE_MI_VS_DELAY Compute mutual information as a function of delay
            
            % Get pulse 2 start indices
            pulse2_start_times = cfg.pair_start_times + cfg.inter_pulse_interval;
            pulse2_start_inds = zeros(size(pulse2_start_times));
            for i = 1:length(pulse2_start_times)
                [~, pulse2_start_inds(i)] = min(abs(t - pulse2_start_times(i)));
            end
            
            % Create delay vector
            max_delay_samples = round(fs * cfg.pulse2_width);
            delay_vec_samples = 1:obj.delay_step_samples:max_delay_samples;
            delay_vec_sec = delay_vec_samples / fs;
            
            MI_vs_delay = zeros(size(delay_vec_samples));
            
            for i_delay = 1:length(delay_vec_samples)
                delay_samples = delay_vec_samples(i_delay);
                
                % Read spike rate at delayed time points
                readout_inds = pulse2_start_inds + delay_samples;
                
                % Make sure we don't go beyond time series
                valid_pairs = readout_inds <= length(r_stim);
                
                if sum(valid_pairs) > 1
                    data_in = cfg.pulse1_amps(valid_pairs);
                    data_out = r_stim(readout_inds(valid_pairs));
                    
                    [MI_delayed, ~] = mutual_info_SISO(data_in, data_out, ...
                        obj.n_bins_in, obj.n_bins_out);
                    MI_vs_delay(i_delay) = MI_delayed;
                else
                    MI_vs_delay(i_delay) = NaN;
                end
            end
        end
        
        function save_results(obj)
            % SAVE_RESULTS Save analysis results to disk
            
            results = obj.results;
            config = struct();
            config.n_networks = obj.n_networks;
            config.conditions = obj.conditions;
            config.model_defaults = obj.model_defaults;
            config.pulse1_width = obj.pulse1_width;
            config.pulse2_width = obj.pulse2_width;
            config.inter_pulse_interval = obj.inter_pulse_interval;
            config.repeat_interval = obj.repeat_interval;
            config.pulse2_amp = obj.pulse2_amp;
            config.n_amp_levels = obj.n_amp_levels;
            config.analysis_start_time = obj.analysis_start_time;
            config.analysis_completed = datestr(now);
            
            save_file = fullfile(obj.output_dir, 'paired_pulse_MI_results.mat');
            save(save_file, 'results', 'config', '-v7.3');
            fprintf('Results saved to: %s\n', save_file);
        end
    end
end
