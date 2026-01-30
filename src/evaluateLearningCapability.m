function learning_metrics = evaluateLearningCapability(esn, params, options)
%EVALUATELEARNINGCAPABILITY Evaluate reservoir learning capability
%
% This function evaluates how well a reservoir can learn to perform
% various tasks, which is essential for assessing the practical utility
% of different parameter configurations.
%
% Inputs:
%   esn - SRNN_ESN object (already initialized)
%   params - Parameter struct with reservoir configuration
%   options - struct with optional fields:
%     task - 'narma10', 'prediction', 'memory' (default: 'prediction')
%     n_train - number of training samples (default: 2000)
%     n_test - number of test samples (default: 500)
%     washout - washout steps (default: 100)
%     lambda - regularization parameter (default: 1e-6)
%     verbose - print progress (default: false)
%
% Outputs:
%   learning_metrics - struct with:
%     train_mse - Training mean squared error
%     test_mse - Test mean squared error
%     train_nrmse - Training normalized RMSE
%     test_nrmse - Test normalized RMSE  
%     R2 - Coefficient of determination on test set
%     memory_capacity - (if task='memory') total memory capacity
%     task - Name of the task performed
%     success - Whether learning evaluation completed successfully

    % Default options
    if nargin < 3
        options = struct();
    end
    
    task = getFieldOrDefault(options, 'task', 'prediction');
    n_train = getFieldOrDefault(options, 'n_train', 2000);
    n_test = getFieldOrDefault(options, 'n_test', 500);
    washout = getFieldOrDefault(options, 'washout', 100);
    lambda = getFieldOrDefault(options, 'lambda', 1e-6);
    verbose = getFieldOrDefault(options, 'verbose', false);
    
    % Initialize output
    learning_metrics = struct();
    learning_metrics.task = task;
    learning_metrics.success = false;
    learning_metrics.train_mse = NaN;
    learning_metrics.test_mse = NaN;
    learning_metrics.train_nrmse = NaN;
    learning_metrics.test_nrmse = NaN;
    learning_metrics.R2 = NaN;
    learning_metrics.memory_capacity = NaN;
    
    try
        dt = getFieldOrDefault(params, 'dt', 0.1);
        
        switch lower(task)
            case 'prediction'
                % Time series prediction task (1-step ahead)
                [U_train, Y_train, U_test, Y_test] = generatePredictionTask(n_train, n_test, dt);
                
            case 'narma10'
                % NARMA-10 task (nonlinear benchmark)
                [U_train, Y_train, U_test, Y_test] = generateNARMA10Task(n_train, n_test);
                
            case 'memory'
                % Memory capacity task
                [U_train, Y_train, U_test, Y_test, delays] = generateMemoryTask(n_train, n_test, 50);
                
            otherwise
                warning('evaluateLearningCapability:UnknownTask', ...
                    'Unknown task: %s. Using prediction task.', task);
                [U_train, Y_train, U_test, Y_test] = generatePredictionTask(n_train, n_test, dt);
        end
        
        % Clear persistent variables to avoid dimension issues
        clear SRNN_reservoir SRNN_reservoir_DDE;
        
        % Reset reservoir state
        esn.resetState();
        
        % Run reservoir on training data
        [X_train, ~] = esn.runReservoir(U_train);
        
        % Apply washout
        if washout >= size(X_train, 1)
            washout = max(1, floor(size(X_train, 1) / 10));
        end
        X_train_use = X_train(washout+1:end, :);
        Y_train_use = Y_train(washout+1:end, :);
        
        % Train readout using ridge regression
        n_features = size(X_train_use, 2);
        XTX = X_train_use' * X_train_use;
        XTY = X_train_use' * Y_train_use;
        W_out = (XTX + lambda * eye(n_features)) \ XTY;
        
        % Compute training predictions
        Y_train_pred = X_train_use * W_out;
        
        % Compute training metrics
        learning_metrics.train_mse = mean((Y_train_use - Y_train_pred).^2, 'all');
        train_std = std(Y_train_use(:));
        if train_std > 1e-10
            learning_metrics.train_nrmse = sqrt(learning_metrics.train_mse) / train_std;
        else
            learning_metrics.train_nrmse = NaN;
        end
        
        % Clear and run on test data
        clear SRNN_reservoir SRNN_reservoir_DDE;
        esn.resetState();
        [X_test, ~] = esn.runReservoir(U_test);
        
        % Compute test predictions
        Y_test_pred = X_test * W_out;
        
        % Compute test metrics
        learning_metrics.test_mse = mean((Y_test - Y_test_pred).^2, 'all');
        test_std = std(Y_test(:));
        if test_std > 1e-10
            learning_metrics.test_nrmse = sqrt(learning_metrics.test_mse) / test_std;
        else
            learning_metrics.test_nrmse = NaN;
        end
        
        % Compute R² (coefficient of determination)
        SS_res = sum((Y_test - Y_test_pred).^2);
        SS_tot = sum((Y_test - mean(Y_test)).^2);
        if SS_tot > 1e-10
            learning_metrics.R2 = 1 - SS_res / SS_tot;
        else
            learning_metrics.R2 = NaN;
        end
        
        % For memory task, compute memory capacity
        if strcmpi(task, 'memory')
            learning_metrics.memory_capacity = computeMemoryCapacity(...
                X_train_use, Y_train_use, X_test, Y_test, delays, lambda);
        end
        
        learning_metrics.success = true;
        
        if verbose
            fprintf('Learning evaluation (%s task):\n', task);
            fprintf('  Train MSE: %.6f, NRMSE: %.4f\n', ...
                learning_metrics.train_mse, learning_metrics.train_nrmse);
            fprintf('  Test MSE: %.6f, NRMSE: %.4f, R²: %.4f\n', ...
                learning_metrics.test_mse, learning_metrics.test_nrmse, learning_metrics.R2);
        end
        
    catch ME
        learning_metrics.success = false;
        learning_metrics.error_message = ME.message;
        if verbose
            warning('evaluateLearningCapability:Failed', ...
                'Learning evaluation failed: %s', ME.message);
        end
    end
end

%% Task Generation Functions

function [U_train, Y_train, U_test, Y_test] = generatePredictionTask(n_train, n_test, dt)
    %GENERATEPREDICTIONTASK Generate a time series prediction task
    % Uses Mackey-Glass chaotic time series for prediction
    
    total_samples = n_train + n_test + 1000;  % Extra for washout
    
    % Try to use Mackey-Glass, fallback to sine mixture
    try
        [~, mg] = generate_mackey_glass('tau', 17, 'n_samples', total_samples, ...
            'dt', dt, 'discard', 500);
        signal = mg(:);
    catch
        % Fallback: mixture of sines (still challenging)
        t = (0:total_samples-1)' * dt;
        signal = 0.5 * sin(0.02 * t) + 0.3 * sin(0.05 * t + 0.5) + ...
                 0.2 * sin(0.11 * t + 1.0) + 0.1 * randn(size(t));
        signal = (signal - mean(signal)) / std(signal) * 0.5;  % Normalize
    end
    
    % Task: predict next value
    U = signal(1:end-1);
    Y = signal(2:end);
    
    % Split into train/test
    U_train = U(1:n_train);
    Y_train = Y(1:n_train);
    U_test = U(n_train+1:n_train+n_test);
    Y_test = Y(n_train+1:n_train+n_test);
end

function [U_train, Y_train, U_test, Y_test] = generateNARMA10Task(n_train, n_test)
    %GENERATENARMA10TASK Generate NARMA-10 benchmark task
    % NARMA-10 is a standard nonlinear benchmark for reservoir computing
    
    total_samples = n_train + n_test + 200;
    
    % Generate uniform random input
    u = 0.5 * rand(total_samples, 1);
    
    % Generate NARMA-10 output
    y = zeros(total_samples, 1);
    
    for t = 11:total_samples
        y(t) = 0.3 * y(t-1) + 0.05 * y(t-1) * sum(y(t-10:t-1)) + ...
               1.5 * u(t-1) * u(t-10) + 0.1;
    end
    
    % Remove initial transient
    u = u(101:end);
    y = y(101:end);
    
    % Split into train/test
    U_train = u(1:n_train);
    Y_train = y(1:n_train);
    U_test = u(n_train+1:n_train+n_test);
    Y_test = y(n_train+1:n_train+n_test);
end

function [U_train, Y_train, U_test, Y_test, delays] = generateMemoryTask(n_train, n_test, max_delay)
    %GENERATEMEMORYTASK Generate memory capacity task
    % Task: reconstruct delayed versions of input
    
    total_samples = n_train + n_test + max_delay + 100;
    
    % Generate random input
    u = rand(total_samples, 1) - 0.5;  % Zero-mean uniform
    
    % Create targets for each delay
    delays = 1:max_delay;
    n_delays = length(delays);
    
    Y = zeros(total_samples - max_delay, n_delays);
    for d = 1:n_delays
        Y(:, d) = u(max_delay - delays(d) + 1:end - delays(d));
    end
    
    % Input is aligned with targets
    U = u(max_delay + 1:end);
    
    % Split into train/test
    U_train = U(1:n_train);
    Y_train = Y(1:n_train, :);
    U_test = U(n_train+1:n_train+n_test);
    Y_test = Y(n_train+1:n_train+n_test, :);
end

function MC = computeMemoryCapacity(X_train, Y_train, X_test, Y_test, delays, lambda)
    %COMPUTEMEMORYCAPACITY Compute total memory capacity
    
    n_delays = length(delays);
    R2_delays = zeros(n_delays, 1);
    n_features = size(X_train, 2);
    
    for d = 1:n_delays
        % Train readout for this delay
        y_train = Y_train(:, d);
        y_test = Y_test(:, d);
        
        % Ridge regression
        XTX = X_train' * X_train;
        XTy = X_train' * y_train;
        w = (XTX + lambda * eye(n_features)) \ XTy;
        
        % Predict
        y_pred = X_test * w;
        
        % Compute R² (squared correlation)
        var_true = var(y_test);
        var_pred = var(y_pred);
        
        if var_true > 1e-12 && var_pred > 1e-12
            cov_matrix = cov(y_test, y_pred);
            cov_val = cov_matrix(1, 2);
            R2_delays(d) = (cov_val^2) / (var_true * var_pred);
            R2_delays(d) = max(0, min(1, R2_delays(d)));
        else
            R2_delays(d) = 0;
        end
    end
    
    % Total memory capacity is sum of R² values
    MC = sum(R2_delays);
end

function val = getFieldOrDefault(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end

