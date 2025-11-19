classdef SRNN_ESN < handle
    % SRNN_ESN: Echo State Network wrapper for the fractional SRNN reservoir
    %
    % This class implements a complete Reservoir Computing / Echo State Network
    % paradigm using the fractional-order SRNN reservoir (SRNN_reservoir.m) as
    % the dynamic core and adding a trainable linear readout layer.
    %
    % Key features:
    %   - Wraps the fractional SRNN reservoir dynamics
    %   - Implements ridge regression for readout training
    %   - Handles time series data with proper temporal splitting and washout
    %   - Configurable feature extraction from reservoir states
    %   - Optional inclusion of raw input as readout features
    %
    % Example usage:
    %   esn = SRNN_ESN(params);
    %   metrics = esn.trainReadout(U_train, Y_train, train_options);
    %   Y_pred = esn.predict(U_test);
    
    properties
        % Reservoir architecture parameters
        n              % Total number of neurons
        n_E            % Number of excitatory neurons
        n_I            % Number of inhibitory neurons
        E_indices      % Indices of excitatory neurons
        I_indices      % Indices of inhibitory neurons
        
        % Adaptation and STD parameters
        n_a_E          % Number of adaptation timescales for E neurons
        n_a_I          % Number of adaptation timescales for I neurons
        n_b_E          % STD flag for E neurons (0 or 1)
        n_b_I          % STD flag for I neurons (0 or 1)
        
        % Network connectivity and dynamics
        W              % Recurrent weight matrix (n x n)
        W_in           % Input weight matrix (n x n_inputs)
        tau_d          % Dendritic time constant
        tau_a_E        % Adaptation time constants for E neurons
        tau_a_I        % Adaptation time constants for I neurons
        tau_b_E_rec    % STD recovery time constant for E neurons
        tau_b_E_rel    % STD release time constant for E neurons
        tau_b_I_rec    % STD recovery time constant for I neurons
        tau_b_I_rel    % STD release time constant for I neurons
        c_E            % Adaptation scaling for E neurons
        c_I            % Adaptation scaling for I neurons
        activation_function  % Nonlinearity (function handle)
        
        % Reservoir state
        S              % Current state vector [a_E(:); a_I(:); b_E(:); b_I(:); x(:)]
        S0             % Initial state (for reset)
        
        % Readout layer
        W_out          % Readout weight matrix
        b_out          % Readout bias vector
        
        % Configuration
        which_states   % Which states to use as features ('x', 'r', 'all')
        include_input  % Whether to include raw input in features
        lambda         % Ridge regression regularization parameter
        
        % Metadata
        n_inputs       % Number of input dimensions
        n_outputs      % Number of output dimensions
        is_trained     % Flag indicating if readout has been trained
    end
    
    methods
        function obj = SRNN_ESN(params)
            % SRNN_ESN: Constructor for the ESN class
            %
            % Input:
            %   params - struct with fields:
            %     Required:
            %       n - total number of neurons
            %       n_E - number of excitatory neurons
            %       n_I - number of inhibitory neurons
            %       W - recurrent weight matrix (n x n)
            %       W_in - input weight matrix (n x n_inputs)
            %       tau_d - dendritic time constant
            %       activation_function - nonlinearity (function handle)
            %     Optional:
            %       n_a_E, n_a_I - adaptation timescales (default: 0)
            %       tau_a_E, tau_a_I - adaptation time constants
            %       n_b_E, n_b_I - STD flags (default: 0)
            %       tau_b_E_rec, tau_b_E_rel, tau_b_I_rec, tau_b_I_rel
            %       c_E, c_I - adaptation scaling (default: 1.0)
            %       which_states - 'x' (default), 'r', or 'all'
            %       include_input - true/false (default: false)
            %       lambda - regularization (default: 1e-6)
            
            % Required parameters
            obj.n = params.n;
            obj.n_E = params.n_E;
            obj.n_I = params.n_I;
            obj.W = params.W;
            obj.W_in = params.W_in;
            obj.tau_d = params.tau_d;
            obj.activation_function = params.activation_function;
            
            % Compute indices
            obj.E_indices = 1:obj.n_E;
            obj.I_indices = (obj.n_E + 1):(obj.n_E + obj.n_I);
            
            % Input dimension
            obj.n_inputs = size(obj.W_in, 2);
            
            % Optional adaptation parameters
            obj.n_a_E = getFieldOrDefault(params, 'n_a_E', 0);
            obj.n_a_I = getFieldOrDefault(params, 'n_a_I', 0);
            obj.tau_a_E = getFieldOrDefault(params, 'tau_a_E', []);
            obj.tau_a_I = getFieldOrDefault(params, 'tau_a_I', []);
            
            % Optional STD parameters
            obj.n_b_E = getFieldOrDefault(params, 'n_b_E', 0);
            obj.n_b_I = getFieldOrDefault(params, 'n_b_I', 0);
            obj.tau_b_E_rec = getFieldOrDefault(params, 'tau_b_E_rec', inf);
            obj.tau_b_E_rel = getFieldOrDefault(params, 'tau_b_E_rel', inf);
            obj.tau_b_I_rec = getFieldOrDefault(params, 'tau_b_I_rec', inf);
            obj.tau_b_I_rel = getFieldOrDefault(params, 'tau_b_I_rel', inf);
            
            % Adaptation scaling
            obj.c_E = getFieldOrDefault(params, 'c_E', 1.0);
            obj.c_I = getFieldOrDefault(params, 'c_I', 1.0);
            
            % Configuration
            obj.which_states = getFieldOrDefault(params, 'which_states', 'x');
            obj.include_input = getFieldOrDefault(params, 'include_input', false);
            obj.lambda = getFieldOrDefault(params, 'lambda', 1e-6);
            
            % Initialize state
            obj.resetState();
            
            % Readout not yet trained
            obj.is_trained = false;
            obj.W_out = [];
            obj.b_out = [];
            obj.n_outputs = 0;
        end
        
        function metrics = trainReadout(obj, U, Y, options)
            % trainReadout: Train the linear readout layer using ridge regression
            %
            % Inputs:
            %   U - Input time series (n_timesteps x n_inputs)
            %   Y - Target time series (n_timesteps x n_outputs)
            %   options - struct with fields:
            %     train_ratio - fraction for training (default: 0.6)
            %     val_ratio - fraction for validation (default: 0.2)
            %     washout_steps - number of initial steps to ignore (default: 100)
            %     lambda - ridge regularization parameter (default: obj.lambda)
            %
            % Output:
            %   metrics - struct with train_mse, val_mse, train_nrmse, val_nrmse
            
            % Parse options
            if nargin < 4
                options = struct();
            end
            train_ratio = getFieldOrDefault(options, 'train_ratio', 0.6);
            val_ratio = getFieldOrDefault(options, 'val_ratio', 0.2);
            washout_steps = getFieldOrDefault(options, 'washout_steps', 100);
            lambda = getFieldOrDefault(options, 'lambda', obj.lambda);
            
            % Update lambda
            obj.lambda = lambda;
            
            % Get dimensions
            n_timesteps = size(U, 1);
            obj.n_outputs = size(Y, 2);
            
            % Compute split indices (temporal order preserved)
            n_train = floor(n_timesteps * train_ratio);
            n_val = floor(n_timesteps * val_ratio);
            
            % Split data temporally
            U_train = U(1:n_train, :);
            Y_train = Y(1:n_train, :);
            U_val = U(n_train+1:n_train+n_val, :);
            Y_val = Y(n_train+1:n_train+n_val, :);
            
            fprintf('Training readout layer...\n');
            fprintf('  Total samples: %d\n', n_timesteps);
            fprintf('  Training samples: %d (after washout: %d)\n', n_train, n_train - washout_steps);
            fprintf('  Validation samples: %d\n', n_val);
            fprintf('  Washout steps: %d\n', washout_steps);
            fprintf('  Lambda (regularization): %.2e\n', lambda);
            
            % Reset reservoir and run over training data
            obj.resetState();
            [X_train, ~] = obj.runReservoir(U_train);
            
            % Apply washout: remove first washout_steps samples
            if washout_steps >= n_train
                error('Washout steps (%d) must be less than training samples (%d)', ...
                      washout_steps, n_train);
            end
            X_train = X_train(washout_steps+1:end, :);
            Y_train_use = Y_train(washout_steps+1:end, :);
            
            % Train readout using ridge regression
            % W_out = (X'*X + lambda*I) \ (X'*Y)
            n_features = size(X_train, 2);
            W_ridge = X_train' * X_train + lambda * eye(n_features);
            obj.W_out = W_ridge \ (X_train' * Y_train_use);
            
            % Compute bias as mean of residuals
            Y_train_pred = X_train * obj.W_out;
            obj.b_out = mean(Y_train_use - Y_train_pred, 1)';
            
            % Evaluate on training set (after washout)
            Y_train_pred = Y_train_pred + obj.b_out';
            train_metrics = compute_metrics(Y_train_pred, Y_train_use);
            
            % Evaluate on validation set
            obj.resetState();
            [X_val, ~] = obj.runReservoir(U_val);
            Y_val_pred = X_val * obj.W_out + obj.b_out';
            val_metrics = compute_metrics(Y_val_pred, Y_val);
            
            % Mark as trained
            obj.is_trained = true;
            
            % Compile metrics
            metrics = struct();
            metrics.train_mse = train_metrics.mse;
            metrics.train_nrmse = train_metrics.nrmse;
            metrics.val_mse = val_metrics.mse;
            metrics.val_nrmse = val_metrics.nrmse;
            
            fprintf('  Training MSE: %.6f, NRMSE: %.4f\n', metrics.train_mse, metrics.train_nrmse);
            fprintf('  Validation MSE: %.6f, NRMSE: %.4f\n', metrics.val_mse, metrics.val_nrmse);
        end
        
        function [Y_pred, X_features] = predict(obj, U)
            % predict: Generate predictions using the trained readout
            %
            % Input:
            %   U - Input time series (n_timesteps x n_inputs)
            %
            % Output:
            %   Y_pred - Predicted outputs (n_timesteps x n_outputs)
            %   X_features - Extracted features (optional, for analysis)
            
            if ~obj.is_trained
                error('SRNN_ESN:NotTrained', ...
                      'Readout layer has not been trained. Call trainReadout() first.');
            end
            
            % Run reservoir over input sequence
            [X_features, ~] = obj.runReservoir(U);
            
            % Apply readout
            Y_pred = X_features * obj.W_out + obj.b_out';
        end
        
        function [Y_gen, X_features] = generateAutonomous(obj, initial_data, n_steps, options)
            % generateAutonomous: Generate predictions in closed-loop (generative mode)
            %
            % This method tests if the reservoir has truly learned the dynamics
            % by autonomously generating predictions where each predicted output
            % is fed back as the next input.
            %
            % Inputs:
            %   initial_data - Initial input sequence for washout (n_washout x n_inputs)
            %                  OR struct with field 'input' containing the sequence
            %   n_steps - Number of steps to generate autonomously
            %   options - struct with optional fields:
            %     horizon - prediction horizon (default: 1 for true closed-loop)
            %     washout_steps - initial steps from initial_data to use (default: all)
            %     return_features - whether to return reservoir features (default: false)
            %
            % Outputs:
            %   Y_gen - Generated output sequence (n_steps x n_outputs)
            %   X_features - Reservoir features during generation (optional)
            
            if ~obj.is_trained
                error('SRNN_ESN:NotTrained', ...
                      'Readout layer has not been trained. Call trainReadout() first.');
            end
            
            % Parse options
            if nargin < 4
                options = struct();
            end
            horizon = getFieldOrDefault(options, 'horizon', 1);
            washout_steps = getFieldOrDefault(options, 'washout_steps', size(initial_data, 1));
            return_features = getFieldOrDefault(options, 'return_features', false);
            
            % Validate inputs
            if washout_steps > size(initial_data, 1)
                warning('SRNN_ESN:WashoutTooLarge', ...
                        'washout_steps (%d) exceeds initial_data length (%d). Using all data.', ...
                        washout_steps, size(initial_data, 1));
                washout_steps = size(initial_data, 1);
            end
            
            fprintf('Generating autonomous predictions...\n');
            fprintf('  Washout steps: %d\n', washout_steps);
            fprintf('  Generation steps: %d\n', n_steps);
            fprintf('  Prediction horizon: %d\n', horizon);
            
            % Phase 1: Washout with initial data
            U_washout = initial_data(1:washout_steps, :);
            [X_washout, ~] = obj.runReservoir(U_washout);
            
            % Phase 2: Autonomous generation
            Y_gen = zeros(n_steps, obj.n_outputs);
            if return_features
                % Pre-allocate based on feature dimension
                n_features = size(obj.W_out, 1);
                X_features = zeros(n_steps, n_features);
            else
                X_features = [];
            end
            
            % Get initial prediction from last washout step
            current_input = X_washout(end, :) * obj.W_out + obj.b_out';
            
            % Prepare parameters for single-step reservoir dynamics
            dt = 1.0; this is in seconds, so 1.0 is a 1 second step, very big
            
            for t = 1:n_steps
                % Create input for this time step (reshape to column vector)
                u_t = (obj.W_in * current_input')';  % (1 x n)
                
                % Integrate reservoir for one time step
                % Need at least 2 time points for ODE solver
                t_span = [0, dt];
                
                % Create interpolation function for constant input over this step
                t_ex = [0, dt];
                u_ex = [u_t; u_t]';  % (n x 2) - constant input over time step
                
                % Pack parameters for SRNN_reservoir
                params = struct();
                params.n = obj.n;
                params.n_E = obj.n_E;
                params.n_I = obj.n_I;
                params.E_indices = obj.E_indices;
                params.I_indices = obj.I_indices;
                params.n_a_E = obj.n_a_E;
                params.n_a_I = obj.n_a_I;
                params.n_b_E = obj.n_b_E;
                params.n_b_I = obj.n_b_I;
                params.W = obj.W;
                params.tau_d = obj.tau_d;
                params.tau_a_E = obj.tau_a_E;
                params.tau_a_I = obj.tau_a_I;
                params.tau_b_E_rec = obj.tau_b_E_rec;
                params.tau_b_E_rel = obj.tau_b_E_rel;
                params.tau_b_I_rec = obj.tau_b_I_rec;
                params.tau_b_I_rel = obj.tau_b_I_rel;
                params.c_E = obj.c_E;
                params.c_I = obj.c_I;
                params.activation_function = obj.activation_function;
                
                % Solve ODE for one time step
                odefun = @(time, S) SRNN_reservoir(time, S, t_ex, u_ex, params);
                options_ode = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
                [~, S_history] = ode23s(odefun, t_span, obj.S, options_ode);
                
                % Update state
                obj.S = S_history(end, :)';
                
                % Extract features from final state
                X_t = obj.extractFeatures(S_history(end, :), current_input);
                
                % Generate prediction
                Y_t = X_t * obj.W_out + obj.b_out';
                Y_gen(t, :) = Y_t;
                
                if return_features
                    X_features(t, :) = X_t;
                end
                
                % Feedback: use prediction as next input
                current_input = Y_t;
            end
            
            fprintf('  Autonomous generation complete!\n');
        end
        
        function resetState(obj)
            % resetState: Reset reservoir to initial conditions
            
            % Compute total state dimension
            len_a_E = obj.n_E * obj.n_a_E;
            len_a_I = obj.n_I * obj.n_a_I;
            len_b_E = obj.n_E * obj.n_b_E;
            len_b_I = obj.n_I * obj.n_b_I;
            len_x = obj.n;
            
            % Initialize state to zeros (or small random values)
            obj.S0 = zeros(len_a_E + len_a_I + len_b_E + len_b_I + len_x, 1);
            
            % Initialize b states to 1 (no depression initially)
            if obj.n_b_E > 0
                idx_start = len_a_E + len_a_I + 1;
                idx_end = idx_start + len_b_E - 1;
                obj.S0(idx_start:idx_end) = 1;
            end
            if obj.n_b_I > 0
                idx_start = len_a_E + len_a_I + len_b_E + 1;
                idx_end = idx_start + len_b_I - 1;
                obj.S0(idx_start:idx_end) = 1;
            end
            
            obj.S = obj.S0;
        end
        
        function S = getState(obj)
            % getState: Get current reservoir state
            S = obj.S;
        end
        
        function setState(obj, S)
            % setState: Set reservoir state
            obj.S = S;
        end
        
        function setHyperparams(obj, hyperparams)
            % setHyperparams: Update hyperparameters
            %
            % Input:
            %   hyperparams - struct with fields to update:
            %     lambda - regularization parameter
            %     which_states - feature extraction mode
            %     include_input - whether to include raw input
            
            if isfield(hyperparams, 'lambda')
                obj.lambda = hyperparams.lambda;
            end
            if isfield(hyperparams, 'which_states')
                obj.which_states = hyperparams.which_states;
            end
            if isfield(hyperparams, 'include_input')
                obj.include_input = hyperparams.include_input;
            end
        end
        
    end
    
    methods (Access = private)
        
        function [X_features, S_history] = runReservoir(obj, U)
            % runReservoir: Simulate reservoir dynamics over input sequence
            %
            % Input:
            %   U - Input time series (n_timesteps x n_inputs)
            %
            % Output:
            %   X_features - Extracted features (n_timesteps x n_features)
            %   S_history - Full state history (for analysis)
            
            n_timesteps = size(U, 1);
            dt = 1.0; % Time step (can be made configurable), this is in seconds, so 1.0 is a 1 second step, very big
            
            % Time vectors for ODE solver
            % Ensure t_span has at least 2 elements for ODE solver
            if n_timesteps == 1
                t_span = [0, dt];
                t_ex = [0, dt];
            else
                t_span = 0:dt:(n_timesteps-1)*dt;
                t_ex = t_span; % Time points for input interpolation
            end
            
            % Prepare external input: u_ex = W_in * U'
            % u_ex should be (n x n_timesteps) for interpolation in SRNN_reservoir
            u_ex = obj.W_in * U'; % (n x n_timesteps)
            
            % For single timestep, replicate input to match t_ex length
            if n_timesteps == 1
                u_ex = [u_ex, u_ex]; % Replicate to (n x 2)
            end
            
            % Pack parameters for SRNN_reservoir
            params = struct();
            params.n = obj.n;
            params.n_E = obj.n_E;
            params.n_I = obj.n_I;
            params.E_indices = obj.E_indices;
            params.I_indices = obj.I_indices;
            params.n_a_E = obj.n_a_E;
            params.n_a_I = obj.n_a_I;
            params.n_b_E = obj.n_b_E;
            params.n_b_I = obj.n_b_I;
            params.W = obj.W;
            params.tau_d = obj.tau_d;
            params.tau_a_E = obj.tau_a_E;
            params.tau_a_I = obj.tau_a_I;
            params.tau_b_E_rec = obj.tau_b_E_rec;
            params.tau_b_E_rel = obj.tau_b_E_rel;
            params.tau_b_I_rec = obj.tau_b_I_rec;
            params.tau_b_I_rel = obj.tau_b_I_rel;
            params.c_E = obj.c_E;
            params.c_I = obj.c_I;
            params.activation_function = obj.activation_function;
            
            % Solve ODE
            odefun = @(t, S) SRNN_reservoir(t, S, t_ex, u_ex, params);
            
            % Use ode23s (stiff solver) for fractional dynamics
            options_ode = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
            [~, S_history] = ode23s(odefun, t_span, obj.S, options_ode); % 
            
            % Update current state
            obj.S = S_history(end, :)';
            
            % For single timestep input, only extract final state
            if n_timesteps == 1
                S_history = S_history(end, :); % Only keep final state
            end
            
            % Extract features from state history
            X_features = obj.extractFeatures(S_history, U);
        end
        
        function X = extractFeatures(obj, S_history, U)
            % extractFeatures: Extract configured state variables as features
            %
            % Input:
            %   S_history - State trajectory (n_timesteps x n_state_vars)
            %   U - Original input time series (n_timesteps x n_inputs)
            %
            % Output:
            %   X - Feature matrix (n_timesteps x n_features)
            
            n_timesteps = size(S_history, 1);
            
            % Compute state dimensions
            len_a_E = obj.n_E * obj.n_a_E;
            len_a_I = obj.n_I * obj.n_a_I;
            len_b_E = obj.n_E * obj.n_b_E;
            len_b_I = obj.n_I * obj.n_b_I;
            
            % Extract dendritic states x
            x_start_idx = len_a_E + len_a_I + len_b_E + len_b_I + 1;
            x_end_idx = x_start_idx + obj.n - 1;
            x_history = S_history(:, x_start_idx:x_end_idx); % (n_timesteps x n)
            
            % Extract features based on configuration
            switch obj.which_states
                case 'x'
                    % Use dendritic states only (default)
                    X = x_history;
                    
                case 'r'
                    % Compute firing rates from x
                    % Need to recompute r = b .* activation_function(x_eff)
                    X = zeros(n_timesteps, obj.n);
                    for t = 1:n_timesteps
                        S_t = S_history(t, :)';
                        [r_t, ~] = obj.computeRates(S_t);
                        X(t, :) = r_t';
                    end
                    
                case 'all'
                    % Use all state variables
                    X = S_history;
                    
                otherwise
                    error('SRNN_ESN:InvalidConfig', ...
                          'which_states must be ''x'', ''r'', or ''all''');
            end
            
            % Optionally include raw input
            if obj.include_input
                X = [X, U];
            end
        end
        
        function [r, x_eff] = computeRates(obj, S)
            % computeRates: Compute firing rates from state vector
            % (Helper function to extract r when needed)
            
            % Unpack state variables (following SRNN_reservoir.m structure)
            current_idx = 0;
            
            len_a_E = obj.n_E * obj.n_a_E;
            if len_a_E > 0
                a_E = reshape(S(current_idx + (1:len_a_E)), obj.n_E, obj.n_a_E);
            else
                a_E = [];
            end
            current_idx = current_idx + len_a_E;
            
            len_a_I = obj.n_I * obj.n_a_I;
            if len_a_I > 0
                a_I = reshape(S(current_idx + (1:len_a_I)), obj.n_I, obj.n_a_I);
            else
                a_I = [];
            end
            current_idx = current_idx + len_a_I;
            
            len_b_E = obj.n_E * obj.n_b_E;
            if len_b_E > 0
                b_E = S(current_idx + (1:len_b_E));
            else
                b_E = [];
            end
            current_idx = current_idx + len_b_E;
            
            len_b_I = obj.n_I * obj.n_b_I;
            if len_b_I > 0
                b_I = S(current_idx + (1:len_b_I));
            else
                b_I = [];
            end
            current_idx = current_idx + len_b_I;
            
            x = S(current_idx + (1:obj.n));
            
            % Compute effective x with adaptation
            x_eff = x;
            if obj.n_E > 0 && obj.n_a_E > 0 && ~isempty(a_E)
                x_eff(obj.E_indices) = x_eff(obj.E_indices) - obj.c_E * sum(a_E, 2);
            end
            if obj.n_I > 0 && obj.n_a_I > 0 && ~isempty(a_I)
                x_eff(obj.I_indices) = x_eff(obj.I_indices) - obj.c_I * sum(a_I, 2);
            end
            
            % Apply STD
            b = ones(obj.n, 1);
            if obj.n_b_E > 0 && ~isempty(b_E)
                b(obj.E_indices) = b_E;
            end
            if obj.n_b_I > 0 && ~isempty(b_I)
                b(obj.I_indices) = b_I;
            end
            
            % Compute firing rates
            r = b .* obj.activation_function(x_eff);
        end
        
    end
end

function value = getFieldOrDefault(s, field, default_value)
    % Helper function to get field from struct or return default
    if isfield(s, field)
        value = s.(field);
    else
        value = default_value;
    end
end

