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
                         'T_wash', 'T_train', 'T_test', 'd_max', 'eta', 'dt_sample'};
            
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
            u_scalar = rand(T_total, 1);  % U(0,1) input
            
            %% Step 2: Run reservoir and collect states
            if verbose
                fprintf('  Running reservoir simulation...\n');
            end
            
            [R_all, t_all] = obj.run_reservoir_esn(u_scalar);
            
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
            
            %% Store results
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
            
            obj.mc_results = mc_results;
        end
        
        function [R_all, t_all] = run_reservoir_esn(obj, u_scalar)
            % RUN_RESERVOIR_ESN Run the reservoir with scalar input sequence
            %
            % [R_all, t_all] = run_reservoir_esn(u_scalar)
            %
            % Inputs:
            %   u_scalar - Scalar input sequence (T x 1)
            %
            % Outputs:
            %   R_all    - Reservoir states (firing rates) (n x T)
            %   t_all    - Time vector (T x 1)
            
            T = length(u_scalar);
            dt = obj.dt_sample;
            
            % Initialize storage
            R_all = zeros(obj.n, T);
            t_all = (0:(T-1))' * dt;
            
            % Initialize state
            params = obj.cached_params;
            S = initialize_state(params);
            
            % Set up ODE options if not already set
            if isempty(obj.ode_opts)
                obj.ode_opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-5, 'MaxStep', dt);
            end
            
            % Map scalar input to neural input
            u_neural_func = @(u) obj.W_in * u;
            
            % Create a temporary interpolant for each step
            params.u_interpolant = [];  % Will be set each step
            
            % Simulation loop
            for t_idx = 1:T
                % Current and next time
                t_now = (t_idx - 1) * dt;
                t_next = t_idx * dt;
                
                % Current input (piecewise constant)
                u_current = u_neural_func(u_scalar(t_idx));
                
                % Create interpolant for this step (constant value)
                % Update the object's u_interpolant so obj.dynamics() can use it
                t_interp = [t_now, t_next];
                u_interp = [u_current, u_current];  % n x 2
                obj.u_interpolant = griddedInterpolant(t_interp, u_interp', 'previous', 'nearest');
                
                % Also update cached_params for consistency
                obj.cached_params.u_interpolant = obj.u_interpolant;
                
                % Define RHS using the public dynamics method
                rhs = @(t, S) obj.dynamics(t, S);
                
                % Integrate one step
                [~, S_out] = obj.ode_solver(rhs, [t_now, t_next], S, obj.ode_opts);
                S = S_out(end, :)';
                
                % Extract firing rates from state
                r = obj.extract_firing_rates(S);
                R_all(:, t_idx) = r;
            end
        end
        
        function r = extract_firing_rates(obj, S)
            % EXTRACT_FIRING_RATES Extract firing rates from state vector
            %
            % r = extract_firing_rates(S)
            %
            % Extracts the firing rates r = phi(x_eff) from the state vector,
            % accounting for adaptation variables.
            
            params = obj.cached_params;
            n = params.n;
            n_E = params.n_E;
            n_I = params.n_I;
            n_a_E = params.n_a_E;
            n_a_I = params.n_a_I;
            n_b_E = params.n_b_E;
            n_b_I = params.n_b_I;
            E_indices = params.E_indices;
            I_indices = params.I_indices;
            c_E = params.c_E;
            c_I = params.c_I;
            activation_fn = params.activation_function;
            
            % Unpack state variables
            current_idx = 0;
            
            % Adaptation states for E neurons
            len_a_E = n_E * n_a_E;
            if len_a_E > 0
                a_E = reshape(S(current_idx + (1:len_a_E)), n_E, n_a_E);
            else
                a_E = [];
            end
            current_idx = current_idx + len_a_E;
            
            % Adaptation states for I neurons
            len_a_I = n_I * n_a_I;
            if len_a_I > 0
                a_I = reshape(S(current_idx + (1:len_a_I)), n_I, n_a_I);
            else
                a_I = [];
            end
            current_idx = current_idx + len_a_I;
            
            % STD states for E neurons
            len_b_E = n_E * n_b_E;
            if len_b_E > 0
                b_E = S(current_idx + (1:len_b_E));
            else
                b_E = [];
            end
            current_idx = current_idx + len_b_E;
            
            % STD states for I neurons
            len_b_I = n_I * n_b_I;
            if len_b_I > 0
                b_I = S(current_idx + (1:len_b_I));
            else
                b_I = [];
            end
            current_idx = current_idx + len_b_I;
            
            % Dendritic states
            x = S(current_idx + (1:n));
            
            % Compute effective input with adaptation
            x_eff = x;
            if n_E > 0 && n_a_E > 0 && ~isempty(a_E)
                x_eff(E_indices) = x_eff(E_indices) - c_E * sum(a_E, 2);
            end
            if n_I > 0 && n_a_I > 0 && ~isempty(a_I)
                x_eff(I_indices) = x_eff(I_indices) - c_I * sum(a_I, 2);
            end
            
            % Apply STD effect for output (b * r)
            b = ones(n, 1);
            if n_b_E > 0 && ~isempty(b_E)
                b(E_indices) = b_E;
            end
            if n_b_I > 0 && ~isempty(b_I)
                b(I_indices) = b_I;
            end
            
            % Compute firing rates
            r = b .* activation_fn(x_eff);
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
        
        function reset(obj)
            % RESET Clear built state and memory capacity results
            
            reset@SRNNModel(obj);
            obj.W_in = [];
            obj.mc_results = [];
            fprintf('ESN reservoir reset.\n');
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
