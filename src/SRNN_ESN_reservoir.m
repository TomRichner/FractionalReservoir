classdef SRNN_ESN_reservoir < SRNNModel
    % SRNN_ESN_RESERVOIR Echo State Network reservoir with memory capacity measurement
    %
    % This class extends SRNNModel to provide Echo State Network (ESN)
    % functionality with the ability to measure memory capacity using the
    % protocol described in Memory_capacity_protocol.md.
    %
    % Usage:
    %   esn = SRNN_ESN_reservoir('n', 100, 'level_of_chaos', 1.8);
    %   esn.build();
    %   MC = esn.run_memory_capacity();
    %
    % Memory Capacity Protocol:
    %   - Drive reservoir with scalar random input u(t) ~ U(0,1)
    %   - For each delay d, train a linear readout to reconstruct u(t-d)
    %   - Compute R^2_d between true delayed input and readout output
    %   - Memory capacity = sum of R^2_d for d = 1 to d_max
    %
    % See also: SRNNModel, Memory_capacity_protocol.md
    
    %% ESN Input Properties
    properties
        f_in = 0.3              % Fraction of neurons receiving input
        sigma_in = 0.5          % Input scaling parameter
        W_in                    % Input weight vector (n x 1)
        rng_seed_input = 3      % RNG seed for input weight generation
        input_type = 'white'    % Input type: 'white' or 'bandlimited'
        u_f_cutoff = []         % Cutoff frequency for bandlimited input (Hz)
                                % If empty, defaults to 1/(2*pi*tau_d)
    end
    
    %% Memory Capacity Protocol Properties
    properties
        T_wash = 1000           % Washout samples (discard transients)
        T_train = 5000          % Training samples
        T_test = 5000           % Test samples
        d_max = 70              % Maximum delay for memory capacity
        eta = 1e-7              % Ridge regression regularization
        dt_sample               % Sampling time step (computed from fs)
    end
    
    %% Memory Capacity Results
    properties (SetAccess = private)
        mc_results              % Struct with memory capacity results
    end
    
    %% Constructor
    methods
        function obj = SRNN_ESN_reservoir(varargin)
            % SRNN_ESN_RESERVOIR Constructor with name-value pairs
            %
            % Usage:
            %   esn = SRNN_ESN_reservoir()  % All defaults
            %   esn = SRNN_ESN_reservoir('n', 200, 'd_max', 100)
            
            % Call superclass constructor first (MATLAB requirement)
            % SRNNModel will ignore unknown properties with a warning
            obj = obj@SRNNModel(varargin{:});
            
            % Define ESN-specific property names (not in SRNNModel)
            esn_props = {'f_in', 'sigma_in', 'W_in', 'rng_seed_input', ...
                         'T_wash', 'T_train', 'T_test', 'd_max', 'eta', 'dt_sample', ...
                         'input_type', 'u_f_cutoff'};
            
            % Parse ESN-specific name-value pairs
            for i = 1:2:length(varargin)
                if ismember(varargin{i}, esn_props)
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
            
            % Set default sampling time step
            obj.dt_sample = 1 / obj.fs;
        end
    end
    
    %% Public Methods
    methods
        function build(obj)
            % BUILD Initialize the ESN reservoir
            %
            % Extends SRNNModel.build() to also generate input weights.
            
            % Call superclass build
            build@SRNNModel(obj);
            
            % Generate input weights
            obj.generate_input_weights();
            
            fprintf('ESN reservoir built. Ready for memory capacity measurement.\n');
        end
        
        function generate_input_weights(obj)
            % GENERATE_INPUT_WEIGHTS Create the input weight vector W_in
            %
            % Creates a sparse input weight vector where f_in fraction of
            % neurons receive input, with weights drawn uniformly from
            % [-sigma_in/2, sigma_in/2].
            
            rng(obj.rng_seed_input);
            
            % Initialize input weights to zero
            obj.W_in = zeros(obj.n, 1);
            
            % Select neurons to receive input
            n_input = round(obj.f_in * obj.n);
            input_neurons = randperm(obj.n, n_input);
            
            % Assign random weights
            obj.W_in(input_neurons) = obj.sigma_in * (rand(n_input, 1) - 0.5);
            
            fprintf('Input weights generated: %d neurons receive input (%.1f%%)\n', ...
                n_input, 100 * n_input / obj.n);
        end
        
        function [MC, R2_d, mc_results] = run_memory_capacity(obj, varargin)
            % RUN_MEMORY_CAPACITY Measure memory capacity of the reservoir
            %
            % [MC, R2_d, results] = run_memory_capacity()
            % [MC, R2_d, results] = run_memory_capacity('verbose', true)
            %
            % Outputs:
            %   MC        - Total memory capacity (sum of R^2_d)
            %   R2_d      - R^2 for each delay d (1 x d_max)
            %   results   - Struct with detailed results
            %
            % Protocol:
            %   1. Generate scalar random input sequence u ~ U(0,1)
            %   2. Run reservoir with piecewise-constant input
            %   3. Collect reservoir states (firing rates)
            %   4. Discard washout period
            %   5. Train linear readouts for each delay via ridge regression
            %   6. Compute R^2_d on test set
            %   7. Sum R^2_d to get memory capacity
            
            if ~obj.is_built
                error('SRNN_ESN_reservoir:NotBuilt', ...
                    'Reservoir must be built first. Call build().');
            end
            
            % Parse optional arguments
            verbose = true;
            for i = 1:2:length(varargin)
                if strcmpi(varargin{i}, 'verbose')
                    verbose = varargin{i+1};
                end
            end
            
            T_total = obj.T_wash + obj.T_train + obj.T_test;
            
            if verbose
                fprintf('Running memory capacity measurement...\n');
                fprintf('  Washout: %d, Train: %d, Test: %d samples\n', ...
                    obj.T_wash, obj.T_train, obj.T_test);
                fprintf('  Max delay: %d\n', obj.d_max);
            end
            
            %% Step 1: Generate scalar random input sequence
            rng(obj.rng_seeds(2));  % Use stimulus seed for reproducibility
            
            if strcmpi(obj.input_type, 'bandlimited')
                % Generate zero-mean white noise for filtering
                u_raw = rand(T_total, 1) - 0.5;
                
                % Determine cutoff frequency
                if isempty(obj.u_f_cutoff)
                    f_cut = 1 / (2 * pi * obj.tau_d);
                else
                    f_cut = obj.u_f_cutoff;
                end
                
                % Design 4th-order Butterworth low-pass filter
                f_nyq = obj.fs / 2;
                [b_filt, a_filt] = butter(3, f_cut / f_nyq, 'low');
                
                % Apply zero-phase filtering
                u_filtered = filtfilt(b_filt, a_filt, u_raw);
                
                % Rescale to [0, 1]
                u_scalar = (u_filtered - min(u_filtered)) / (max(u_filtered) - min(u_filtered));
                
                if verbose
                    fprintf('  Using bandlimited input (f_cutoff = %.2f Hz)\n', f_cut);
                end
            else
                u_scalar = rand(T_total, 1);  % U(0,1) white noise
            end
            
            %% Step 2: Run reservoir and collect states
            if verbose
                fprintf('  Running reservoir simulation...\n');
            end
            
            obj.run_reservoir_esn(u_scalar);
            t_all = obj.t_out;
            
            % Unpack states using standard utility
            [x_all, a_all, b_all, r_all, br_all] = unpack_and_compute_states(obj.S_out, obj.cached_params);
            
            % Combine E and I firing rates for training (need n x T matrix)
            R_all = [r_all.E; r_all.I];  % n x T
            
            %% Step 3: Discard washout and split into train/test
            T_eff = obj.T_train + obj.T_test;
            
            % Reservoir states after washout
            R_eff = R_all(:, (obj.T_wash + 1):end);  % n x T_eff
            u_eff = u_scalar((obj.T_wash + 1):end);  % T_eff x 1
            
            % Split into training and test
            R_train = R_eff(:, 1:obj.T_train)';          % T_train x n
            R_test = R_eff(:, (obj.T_train + 1):end)';   % T_test x n
            u_train = u_eff(1:obj.T_train);              % T_train x 1
            u_test = u_eff((obj.T_train + 1):end);       % T_test x 1
            
            %% Step 4-5: Train readouts and compute R^2 for each delay
            if verbose
                fprintf('  Training readouts for %d delays...\n', obj.d_max);
            end
            
            R2_d = zeros(1, obj.d_max);
            weights_all = zeros(obj.n, obj.d_max);
            
            % Store predictions for selected delays (for plotting)
            predictions = struct();
            
            for d = 1:obj.d_max
                % Build delayed targets for training
                % Target: u(t-d) where t is current time
                % For training indices (d+1):T_train, target is u(1:(T_train-d))
                train_indices = (d + 1):obj.T_train;
                target_indices = 1:(obj.T_train - d);
                
                if length(train_indices) < 10
                    warning('SRNN_ESN_reservoir:InsufficientData', ...
                        'Delay %d has insufficient training data.', d);
                    continue;
                end
                
                X_train = R_train(train_indices, :);
                y_train = u_train(target_indices);
                
                % Train linear readout with ridge regression
                w_d = obj.train_linear_readout(X_train, y_train, obj.eta);
                weights_all(:, d) = w_d;
                
                % Build delayed targets for test
                % For test indices (d+1):T_test, target is u_test(1:(T_test-d))
                test_indices = (d + 1):obj.T_test;
                test_target_indices = 1:(obj.T_test - d);
                
                if length(test_indices) < 10
                    continue;
                end
                
                X_test = R_test(test_indices, :);
                y_test_true = u_test(test_target_indices);
                
                % Compute predictions
                y_test_pred = X_test * w_d;
                
                % Store predictions for this delay
                predictions(d).y_true = y_test_true;
                predictions(d).y_pred = y_test_pred;
                predictions(d).t_indices = test_indices;
                
                % Compute R^2
                R2_d(d) = obj.compute_R2(y_test_true, y_test_pred);
                
                if verbose && mod(d, 10) == 0
                    fprintf('    Delay %d: R^2 = %.4f\n', d, R2_d(d));
                end
            end
            
            %% Step 6: Compute total memory capacity
            MC = sum(R2_d);
            
            if verbose
                fprintf('  Total Memory Capacity: %.4f\n', MC);
            end
            
            %% Store results including time series for plotting
            mc_results = struct();
            mc_results.MC = MC;
            mc_results.R2_d = R2_d;
            mc_results.d = 1:obj.d_max;
            mc_results.weights = weights_all;
            mc_results.T_wash = obj.T_wash;
            mc_results.T_train = obj.T_train;
            mc_results.T_test = obj.T_test;
            mc_results.eta = obj.eta;
            mc_results.u_scalar = u_scalar;
            mc_results.predictions = predictions;
            
            % Store time series data for test period (for plotting)
            % Use decimated data following SRNNModel pattern
            test_start_idx = obj.T_wash + obj.T_train + 1;
            mc_results.t_test = t_all(test_start_idx:end);
            mc_results.u_test = u_test;
            
            % Store unpacked states for test period
            mc_results.x = x_all;
            mc_results.a = a_all;
            mc_results.b = b_all;
            mc_results.r = r_all;
            mc_results.br = br_all;
            
            %% Step 7: Compute Lyapunov exponents using parent class method
            if ~strcmpi(obj.lya_method, 'none')
                obj.lya_T_interval = [t_all(obj.T_wash + 1), t_all(end)];  % After washout
                obj.compute_lyapunov();
            end
            
            mc_results.lya_results = obj.lya_results;
            
            obj.mc_results = mc_results;
        end
        
        function run_reservoir_esn(obj, u_scalar)
            % RUN_RESERVOIR_ESN Run reservoir with single ODE integration
            %
            % run_reservoir_esn(u_scalar)
            %
            % Runs the reservoir simulation using a single ODE solver call,
            % matching the pattern used in SRNNModel.run().
            %
            % Inputs:
            %   u_scalar - Scalar input sequence (T x 1)
            %
            % Results are stored in inherited properties:
            %   obj.S_out - State trajectory (T x N_sys_eqs)
            %   obj.t_out - Time vector (T x 1)
            
            % Generate ESN stimulus (sets u_ex, t_ex, u_interpolant, S0)
            obj.generate_esn_stimulus(u_scalar);
            
            params = obj.cached_params;
            dt = 1 / obj.fs;
            
            % Set up ODE options (matching SRNNModel pattern)
            if isempty(obj.ode_opts)
                obj.ode_opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-5, 'MaxStep', dt);
            end
            
            % Define RHS function using static method (avoids OOP overhead)
            params.u_interpolant = obj.u_interpolant;
            rhs = @(t, S) SRNNModel.dynamics_fast(t, S, params);
            
            % Integrate entire trajectory at once
            fprintf('  Integrating ESN dynamics...\n');
            tic
            [obj.t_out, obj.S_out] = obj.ode_solver(rhs, obj.t_ex, obj.S0, obj.ode_opts);
            integration_time = toc;
            fprintf('  Integration complete in %.2f seconds.\n', integration_time);
        end
        
        function [fig_handle, ax_handles] = plot_memory_capacity(obj, varargin)
            % PLOT_MEMORY_CAPACITY Plot memory capacity results
            %
            % [fig, axes] = plot_memory_capacity()
            %
            % Creates a figure showing:
            %   - R^2_d vs delay d
            %   - Total MC value
            
            if isempty(obj.mc_results)
                error('SRNN_ESN_reservoir:NoResults', ...
                    'No memory capacity results. Run run_memory_capacity() first.');
            end
            
            fig_handle = figure('Position', [100, 100, 800, 400]);
            
            % Plot R^2 vs delay
            ax1 = subplot(1, 2, 1);
            bar(obj.mc_results.d, obj.mc_results.R2_d, 'FaceColor', [0.3, 0.6, 0.9]);
            xlabel('Delay d (samples)');
            ylabel('R^2_d');
            title('Memory Capacity by Delay');
            grid on;
            ax_handles(1) = ax1;
            
            % Plot cumulative MC
            ax2 = subplot(1, 2, 2);
            cumMC = cumsum(obj.mc_results.R2_d);
            plot(obj.mc_results.d, cumMC, 'b-', 'LineWidth', 2);
            xlabel('Delay d (samples)');
            ylabel('Cumulative MC');
            title(sprintf('Cumulative Memory Capacity (Total: %.2f)', obj.mc_results.MC));
            grid on;
            ax_handles(2) = ax2;
            
            sgtitle(sprintf('Memory Capacity Analysis (MC = %.2f)', obj.mc_results.MC));
        end
        
        function [fig_handle, ax_handles] = plot_esn_timeseries(obj, delays_to_plot, varargin)
            % PLOT_ESN_TIMESERIES Plot comprehensive time series for ESN reservoir
            %
            % [fig, axes] = plot_esn_timeseries(delays_to_plot)
            % [fig, axes] = plot_esn_timeseries(delays_to_plot, 'title', 'My Title')
            %
            % Inputs:
            %   delays_to_plot - Vector of delays to show (e.g., [1, 10, 30, 50])
            %
            % Creates a figure with:
            %   - u(t): Scalar input
            %   - x(t): Dendritic states (line plots)
            %   - r(t): Firing rates (line plots)
            %   - br(t): Synaptic output (if STD enabled)
            %   - a(t): Adaptation states (if SFA enabled)
            %   - b(t): STD states (if STD enabled)
            %   - Lyapunov exponent (if computed)
            %   - For each delay d: u(t-d) vs y_d(t) overlay
            
            if isempty(obj.mc_results) || ~isfield(obj.mc_results, 'x')
                error('SRNN_ESN_reservoir:NoResults', ...
                    'No time series data. Run run_memory_capacity() first.');
            end
            
            if nargin < 2 || isempty(delays_to_plot)
                delays_to_plot = [1, 10, 30, 50];
            end
            
            % Filter delays that exist in results
            delays_to_plot = delays_to_plot(delays_to_plot <= obj.d_max);
            n_delays = length(delays_to_plot);
            
            % Parse optional arguments
            fig_title = '';
            for i = 1:2:length(varargin)
                if strcmpi(varargin{i}, 'title')
                    fig_title = varargin{i+1};
                end
            end
            
            % Get stored data (using new format from unpack_and_compute_states)
            t = obj.mc_results.t_test;
            u_scalar = obj.mc_results.u_test;
            x = obj.mc_results.x;      % Struct with .E and .I
            r = obj.mc_results.r;      % Struct with .E and .I
            br = obj.mc_results.br;    % Struct with .E and .I
            a = obj.mc_results.a;      % Struct with .E and .I
            b = obj.mc_results.b;      % Struct with .E and .I
            params = obj.cached_params;
            
            % Slice to test period only
            test_start_idx = obj.T_wash + obj.T_train + 1;
            test_end_idx = size(x.E, 2);
            test_indices = test_start_idx:test_end_idx;
            
            % Determine which subplots are needed
            has_adaptation = params.n_a_E > 0 || params.n_a_I > 0;
            has_std = params.n_b_E > 0 || params.n_b_I > 0;
            
            % Check for Lyapunov results
            has_lyapunov = isfield(obj.mc_results, 'lya_results') && ...
                           ~isempty(obj.mc_results.lya_results) && ...
                           isfield(obj.mc_results.lya_results, 'LLE');
            
            % Calculate number of base plots
            n_base_plots = 3;  % u(t), x(t), r(t) always present
            if has_std && ~isempty(br.E)
                n_base_plots = n_base_plots + 1;  % br(t)
            end
            if has_adaptation && (~isempty(a.E) || ~isempty(a.I))
                n_base_plots = n_base_plots + 1;  % a(t)
            end
            if has_std
                n_base_plots = n_base_plots + 1;  % b(t)
            end
            if has_lyapunov
                n_base_plots = n_base_plots + 1;  % Lyapunov exponent
            end
            
            n_total_plots = n_base_plots + n_delays;
            
            % Create figure
            fig_height = min(100 + 80 * n_total_plots, 1000);
            fig_handle = figure('Position', [100, 100, 1000, fig_height]);
            tiledlayout(n_total_plots, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
            
            ax_handles = [];
            
            % Get colormaps
            cmap_I = inhibitory_colormap(8);
            cmap_E = excitatory_colormap(8);
            
            %% Plot 1: Scalar input u(t)
            ax_handles(end+1) = nexttile;
            plot(t, u_scalar, 'k-', 'LineWidth', 0.5);
            ylabel('u(t)');
            set(gca, 'XTickLabel', []);
            grid on;
            
            %% Plot 2: Dendritic states x(t) as line plots
            ax_handles(end+1) = nexttile;
            % Use test period indices
            x_E_test = x.E(:, test_indices);
            x_I_test = x.I(:, test_indices);
            % Plot I neurons first (background)
            plot_lines_with_colormap(t, x_I_test, cmap_I);
            hold on;
            % Plot E neurons on top
            plot_lines_with_colormap(t, x_E_test, cmap_E);
            hold off;
            ylabel('dendrite');
            set(gca, 'XTickLabel', []);
            
            %% Plot 3: Firing rates r(t) as line plots
            ax_handles(end+1) = nexttile;
            r_E_test = r.E(:, test_indices);
            r_I_test = r.I(:, test_indices);
            % Plot I neurons first (background)
            plot_lines_with_colormap(t, r_I_test, cmap_I);
            hold on;
            % Plot E neurons on top
            plot_lines_with_colormap(t, r_E_test, cmap_E);
            hold off;
            ylabel('firing rate');
            ylim([0, 1]);
            yticks([0, 1]);
            set(gca, 'XTickLabel', []);
            
            %% Plot 4 (conditional): Synaptic output br(t)
            if has_std && ~isempty(br.E)
                ax_handles(end+1) = nexttile;
                br_E_test = br.E(:, test_indices);
                br_I_test = br.I(:, test_indices);
                % Plot I neurons first (background)
                plot_lines_with_colormap(t, br_I_test, cmap_I);
                hold on;
                % Plot E neurons on top
                plot_lines_with_colormap(t, br_E_test, cmap_E);
                hold off;
                ylabel('synaptic output');
                ylim([0, 1]);
                yticks([0, 1]);
                set(gca, 'XTickLabel', []);
            end
            
            %% Plot 5 (conditional): Adaptation states a(t)
            if has_adaptation && (~isempty(a.E) || ~isempty(a.I))
                ax_handles(end+1) = nexttile;
                has_plotted = false;
                % Plot I adaptation first (background)
                if ~isempty(a.I) && params.n_a_I > 0
                    % Sum across timescales: (n_I x n_a_I x T) -> (n_I x T)
                    a_I_sum = squeeze(sum(a.I(:, :, test_indices), 2));
                    if size(a_I_sum, 2) == 1
                        a_I_sum = a_I_sum';  % Ensure n_I x T
                    end
                    plot_lines_with_colormap(t, a_I_sum, cmap_I);
                    has_plotted = true;
                end
                % Plot E adaptation on top
                if ~isempty(a.E) && params.n_a_E > 0
                    if has_plotted
                        hold on;
                    end
                    % Sum across timescales: (n_E x n_a_E x T) -> (n_E x T)
                    a_E_sum = squeeze(sum(a.E(:, :, test_indices), 2));
                    if size(a_E_sum, 2) == 1
                        a_E_sum = a_E_sum';  % Ensure n_E x T
                    end
                    plot_lines_with_colormap(t, a_E_sum, cmap_E);
                end
                hold off;
                ylabel('adaptation');
                set(gca, 'XTickLabel', []);
                grid on;
            end
            
            %% Plot 6 (conditional): STD states b(t)
            if has_std
                ax_handles(end+1) = nexttile;
                has_plotted = false;
                % Plot I STD first (background)
                if ~isempty(b.I) && params.n_b_I > 0
                    b_I_test = b.I(:, test_indices);
                    if ~all(b_I_test(:) == 1)  % Check if actual STD dynamics
                        plot_lines_with_colormap(t, b_I_test, cmap_I);
                        has_plotted = true;
                    end
                end
                % Plot E STD on top
                if ~isempty(b.E) && params.n_b_E > 0
                    b_E_test = b.E(:, test_indices);
                    if ~all(b_E_test(:) == 1)  % Check if actual STD dynamics
                        if has_plotted
                            hold on;
                        end
                        plot_lines_with_colormap(t, b_E_test, cmap_E);
                    end
                end
                hold off;
                ylabel('depression');
                ylim([0, 1]);
                yticks([0, 1]);
                set(gca, 'XTickLabel', []);
                grid on;
            end
            
            %% Plot (conditional): Lyapunov exponent
            if has_lyapunov
                ax_handles(end+1) = nexttile;
                plot_lyapunov(obj.mc_results.lya_results, 'benettin', {'filtered', 'local', 'EOC', 'value'});
                set(gca, 'XTickLabel', []);
            end
            
            %% Delay reconstruction plots
            colors_pred = lines(n_delays);
            
            for i = 1:n_delays
                d = delays_to_plot(i);
                ax_handles(end+1) = nexttile;
                
                if d <= length(obj.mc_results.predictions) && ...
                        isfield(obj.mc_results.predictions(d), 'y_true') && ...
                        ~isempty(obj.mc_results.predictions(d).y_true)
                    
                    y_true = obj.mc_results.predictions(d).y_true;
                    y_pred = obj.mc_results.predictions(d).y_pred;
                    t_indices = obj.mc_results.predictions(d).t_indices;
                    t_delay = t(t_indices);
                    
                    % Plot true delayed input
                    plot(t_delay, y_true, 'k-', 'LineWidth', 0.8, 'DisplayName', sprintf('u(t-%d)', d));
                    hold on;
                    % Plot prediction
                    plot(t_delay, y_pred, '-', 'Color', colors_pred(i,:), 'LineWidth', 1.2, ...
                        'DisplayName', sprintf('\\hat{y}_{%d}(t)', d));
                    hold off;
                    
                    R2 = obj.mc_results.R2_d(d);
                    title(sprintf('Delay d=%d: R^2=%.3f', d, R2), 'FontWeight', 'normal');
                    ylabel('Value');
                    legend('Location', 'best', 'FontSize', 7);
                else
                    text(0.5, 0.5, sprintf('No data for delay %d', d), ...
                        'HorizontalAlignment', 'center');
                    title(sprintf('Delay d=%d', d), 'FontWeight', 'normal');
                end
                
                if i < n_delays
                    set(gca, 'XTickLabel', []);
                else
                    xlabel('Time (s)');
                end
                grid on;
            end
            
            % Link all axes
            linkaxes(ax_handles, 'x');
            
            % Limit x-axis to test period
            xlim([t(1), t(end)]);
            
            % Add figure title
            if ~isempty(fig_title)
                sgtitle(fig_title);
            else
                sgtitle(sprintf('ESN Time Series (MC = %.2f)', obj.mc_results.MC));
            end
        end
        
        function reset(obj)
            % RESET Clear built state and memory capacity results
            
            reset@SRNNModel(obj);
            obj.W_in = [];
            obj.mc_results = [];
            fprintf('ESN reservoir reset.\n');
        end
    end
    
    %% Private Methods
    methods (Access = private)
        function generate_esn_stimulus(obj, u_scalar)
            % GENERATE_ESN_STIMULUS Generate ESN input stimulus from scalar input
            %
            % Builds the full u_ex matrix, t_ex time vector, u_interpolant,
            % and initial state S0. Follows the pattern in SRNNModel.generate_stimulus().
            %
            % Inputs:
            %   u_scalar - Scalar input sequence (T x 1)
            
            T = length(u_scalar);
            dt = 1 / obj.fs;
            
            % Build time vector
            obj.t_ex = (0:(T-1))' * dt;
            
            % Map scalar input to neural input
            obj.u_ex = obj.W_in * u_scalar';  % n x T
            
            % Create interpolant (piecewise constant, matching ESN discrete input)
            obj.u_interpolant = griddedInterpolant(obj.t_ex, obj.u_ex', 'previous', 'nearest');
            obj.cached_params.u_interpolant = obj.u_interpolant;
            
            % Initialize state
            obj.S0 = initialize_state(obj.cached_params);
        end
    end
    
    %% Static Methods for Ridge Regression and R^2
    methods (Static)
        function w = train_linear_readout(X, y, eta)
            % TRAIN_LINEAR_READOUT Train a linear readout using ridge regression
            %
            % w = train_linear_readout(X, y, eta)
            %
            % Inputs:
            %   X   - Design matrix (T x n_features)
            %   y   - Target vector (T x 1)
            %   eta - Regularization parameter
            %
            % Outputs:
            %   w   - Weight vector (n_features x 1)
            %
            % Closed-form solution:
            %   w = (X'X + eta*I)^(-1) * X'y
            
            n_features = size(X, 2);
            XTX = X' * X;
            XTy = X' * y;
            w = (XTX + eta * eye(n_features)) \ XTy;
        end
        
        function R2 = compute_R2(y_true, y_pred)
            % COMPUTE_R2 Compute coefficient of determination
            %
            % R2 = compute_R2(y_true, y_pred)
            %
            % Computes R^2 as per the memory capacity protocol:
            %   R^2 = cov(y_true, y_pred)^2 / (var(y_true) * var(y_pred))
            %
            % This is the squared correlation coefficient.
            
            % Ensure column vectors
            y_true = y_true(:);
            y_pred = y_pred(:);
            
            % Handle degenerate cases
            var_true = var(y_true);
            var_pred = var(y_pred);
            
            if var_true < 1e-12 || var_pred < 1e-12
                R2 = 0;
                return;
            end
            
            % Compute squared correlation
            cov_matrix = cov(y_true, y_pred);
            cov_val = cov_matrix(1, 2);
            
            R2 = (cov_val^2) / (var_true * var_pred);
            
            % Clamp to [0, 1]
            R2 = max(0, min(1, R2));
        end
    end
end
