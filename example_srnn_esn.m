%% Example: Fractional Reservoir Echo State Network for Mackey-Glass Prediction
% 
% This script demonstrates the complete Reservoir Computing / Echo State Network
% paradigm using the fractional-order SRNN reservoir with spike-frequency
% adaptation and short-term synaptic depression.
%
% Task: One-step-ahead prediction of the chaotic Mackey-Glass time series
%
% The script:
%   1. Generates Mackey-Glass time series data
%   2. Creates a fractional SRNN reservoir with ESN readout
%   3. Trains the readout layer using ridge regression
%   4. Evaluates on validation and test sets
%   5. Visualizes predictions and performance

clear; close all; clc;

fprintf('=================================================================\n');
fprintf('  Fractional Reservoir ESN: Mackey-Glass Time Series Prediction\n');
fprintf('=================================================================\n\n');

%% 1. Generate Mackey-Glass Time Series
fprintf('Step 1: Generating Mackey-Glass time series...\n');

% Generate chaotic time series (tau=17 for chaotic regime)
[t, x_mg] = generate_mackey_glass('tau', 17, ...
                                   'n_samples', 5000, ...
                                   'dt', 1.0, ...
                                   'discard', 1000);

fprintf('  Generated %d samples\n', length(x_mg));
fprintf('  Time range: [%.1f, %.1f]\n', t(1), t(end));
fprintf('  Value range: [%.3f, %.3f]\n\n', min(x_mg), max(x_mg));

%% 2. Create Prediction Task
fprintf('Step 2: Creating prediction task...\n');

% One-step-ahead prediction: predict x(t+1) from x(t)
prediction_horizon = 1;

% Input: x(t)
U = x_mg(1:end-prediction_horizon);

% Target: x(t+prediction_horizon)
Y = x_mg(1+prediction_horizon:end);

fprintf('  Task: Predict x(t+%d) from x(t)\n', prediction_horizon);
fprintf('  Input dimension: %d\n', size(U, 2));
fprintf('  Output dimension: %d\n\n', size(Y, 2));

%% 3. Initialize Fractional SRNN-ESN
fprintf('Step 3: Initializing fractional reservoir ESN...\n');

% Reservoir architecture
n = 200;                    % Total number of neurons
fraction_E = 0.8;           % 80% excitatory
n_E = round(n * fraction_E);
n_I = n - n_E;

fprintf('  Reservoir size: %d neurons (%d E, %d I)\n', n, n_E, n_I);

% Fractional-order parameters: Multiple adaptation timescales
n_a_E = 3;  % 3 adaptation timescales for E neurons
n_a_I = 2;  % 2 adaptation timescales for I neurons
tau_a_E = [100, 500, 2000];  % Fast, medium, slow adaptation
tau_a_I = [50, 300];

fprintf('  Adaptation timescales (E): [%s]\n', num2str(tau_a_E));
fprintf('  Adaptation timescales (I): [%s]\n', num2str(tau_a_I));

% Short-term depression parameters
n_b_E = 1;  % Enable STD for E neurons
n_b_I = 0;  % No STD for I neurons
tau_b_E_rec = 800;   % Recovery time constant
tau_b_E_rel = 100;   % Release time constant

fprintf('  STD enabled for E neurons (tau_rec=%.0f, tau_rel=%.0f)\n', ...
        tau_b_E_rec, tau_b_E_rel);

% Dendritic time constant
tau_d = 10;

% Adaptation scaling factors
c_E = 0.5;  % Moderate adaptation effect for E neurons
c_I = 0.3;  % Weaker adaptation for I neurons

% Create recurrent weight matrix W with spectral radius scaling
rng(42);  % For reproducibility
W = randn(n, n);

% Dale's law: E neurons have positive weights, I neurons negative
W(:, 1:n_E) = abs(W(:, 1:n_E));      % E columns positive
W(:, n_E+1:end) = -abs(W(:, n_E+1:end));  % I columns negative

% Scale to desired spectral radius
spectral_radius = 1.2;
rho = max(abs(eig(W)));
W = (spectral_radius / rho) * W;

fprintf('  Spectral radius: %.2f\n', spectral_radius);

% Create input weight matrix with uniform random scaling
input_scaling = 0.5;
n_inputs = size(U, 2);
W_in = (2 * rand(n, n_inputs) - 1) * input_scaling;

fprintf('  Input scaling: %.2f\n', input_scaling);

% Activation function (ReLU-like with threshold)
activation_function = @(x) max(0, x);

% Pack parameters
params = struct();
params.n = n;
params.n_E = n_E;
params.n_I = n_I;
params.W = W;
params.W_in = W_in;
params.tau_d = tau_d;
params.n_a_E = n_a_E;
params.n_a_I = n_a_I;
params.tau_a_E = tau_a_E;
params.tau_a_I = tau_a_I;
params.n_b_E = n_b_E;
params.n_b_I = n_b_I;
params.tau_b_E_rec = tau_b_E_rec;
params.tau_b_E_rel = tau_b_E_rel;
params.tau_b_I_rec = inf;  % No STD for I neurons
params.tau_b_I_rel = inf;
params.c_E = c_E;
params.c_I = c_I;
params.activation_function = activation_function;
params.which_states = 'x';  % Use dendritic states as features (default)
params.include_input = false;  % Don't include raw input in features
params.lambda = 1e-6;  % Ridge regularization

% Create ESN
esn = SRNN_ESN(params);

fprintf('  Feature extraction: using dendritic states (x)\n');
fprintf('  Include raw input: %s\n\n', mat2str(params.include_input));

%% 4. Train the Readout Layer
fprintf('Step 4: Training readout layer...\n');

% Training configuration
train_options = struct();
train_options.train_ratio = 0.6;      % 60% for training
train_options.val_ratio = 0.2;        % 20% for validation
train_options.washout_steps = 200;    % Washout period
train_options.lambda = 1e-6;          % Ridge regularization

fprintf('  Data split: %.0f%% train, %.0f%% val, %.0f%% test\n', ...
        train_options.train_ratio*100, ...
        train_options.val_ratio*100, ...
        (1 - train_options.train_ratio - train_options.val_ratio)*100);

% Train
tic;
metrics = esn.trainReadout(U, Y, train_options);
train_time = toc;

fprintf('  Training completed in %.2f seconds\n\n', train_time);

%% 5. Evaluate on Test Set
fprintf('Step 5: Evaluating on test set...\n');

% Compute split indices
n_total = length(U);
n_train = floor(n_total * train_options.train_ratio);
n_val = floor(n_total * train_options.val_ratio);
n_test = n_total - n_train - n_val;

% Test data
U_test = U(n_train+n_val+1:end, :);
Y_test = Y(n_train+n_val+1:end, :);

fprintf('  Test samples: %d\n', n_test);

% Reset and predict
esn.resetState();
tic;
Y_test_pred = esn.predict(U_test);
test_time = toc;

% Compute test metrics
test_metrics = compute_metrics(Y_test_pred, Y_test);

fprintf('  Test MSE: %.6f\n', test_metrics.mse);
fprintf('  Test NRMSE: %.4f\n', test_metrics.nrmse);
fprintf('  Prediction completed in %.2f seconds\n\n', test_time);

%% 6. Visualization
fprintf('Step 6: Creating visualizations...\n');

% Prepare data for plotting
n_plot_samples = min(500, n_test);  % Plot up to 500 samples

% Training data (after washout)
washout = train_options.washout_steps;
U_train_full = U(1:n_train, :);
Y_train_full = Y(1:n_train, :);
esn.resetState();
Y_train_pred = esn.predict(U_train_full);
Y_train_pred = Y_train_pred(washout+1:end, :);
Y_train_full = Y_train_full(washout+1:end, :);
t_train = (washout:n_train-1)';

% Validation data
U_val = U(n_train+1:n_train+n_val, :);
Y_val = Y(n_train+1:n_train+n_val, :);
esn.resetState();
Y_val_pred = esn.predict(U_val);
t_val = (n_train:n_train+n_val-1)';

% Test data (use subset for clearer visualization)
t_test = (n_train+n_val:n_train+n_val+n_plot_samples-1)';
Y_test_plot = Y_test(1:n_plot_samples, :);
Y_test_pred_plot = Y_test_pred(1:n_plot_samples, :);

% Create figure
figure('Position', [100, 100, 1200, 800]);

% Subplot 1: Training predictions
subplot(3, 2, 1);
plot(t_train(1:min(300, end)), Y_train_full(1:min(300, end)), 'b-', 'LineWidth', 1.5);
hold on;
plot(t_train(1:min(300, end)), Y_train_pred(1:min(300, end)), 'r--', 'LineWidth', 1);
xlabel('Time Step');
ylabel('Value');
title('Training Set Predictions');
legend('Target', 'Predicted', 'Location', 'best');
grid on;

% Subplot 2: Training error
subplot(3, 2, 2);
plot(t_train, abs(Y_train_full - Y_train_pred), 'k-', 'LineWidth', 0.5);
xlabel('Time Step');
ylabel('Absolute Error');
title(sprintf('Training Error (MSE=%.6f)', metrics.train_mse));
grid on;

% Subplot 3: Validation predictions
subplot(3, 2, 3);
plot(t_val, Y_val, 'b-', 'LineWidth', 1.5);
hold on;
plot(t_val, Y_val_pred, 'r--', 'LineWidth', 1);
xlabel('Time Step');
ylabel('Value');
title('Validation Set Predictions');
legend('Target', 'Predicted', 'Location', 'best');
grid on;

% Subplot 4: Validation error
subplot(3, 2, 4);
plot(t_val, abs(Y_val - Y_val_pred), 'k-', 'LineWidth', 0.5);
xlabel('Time Step');
ylabel('Absolute Error');
title(sprintf('Validation Error (MSE=%.6f)', metrics.val_mse));
grid on;

% Subplot 5: Test predictions
subplot(3, 2, 5);
plot(t_test, Y_test_plot, 'b-', 'LineWidth', 1.5);
hold on;
plot(t_test, Y_test_pred_plot, 'r--', 'LineWidth', 1);
xlabel('Time Step');
ylabel('Value');
title('Test Set Predictions');
legend('Target', 'Predicted', 'Location', 'best');
grid on;

% Subplot 6: Test error
subplot(3, 2, 6);
plot(t_test, abs(Y_test_plot - Y_test_pred_plot), 'k-', 'LineWidth', 0.5);
xlabel('Time Step');
ylabel('Absolute Error');
title(sprintf('Test Error (MSE=%.6f)', test_metrics.mse));
grid on;

% Overall title
sgtitle('Fractional SRNN-ESN: Mackey-Glass Time Series Prediction', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Adjust layout
set(gcf, 'Color', 'w');

fprintf('  Visualization complete!\n\n');

%% 7. Summary
fprintf('=================================================================\n');
fprintf('  SUMMARY\n');
fprintf('=================================================================\n');
fprintf('Reservoir Configuration:\n');
fprintf('  - Size: %d neurons (%d E, %d I)\n', n, n_E, n_I);
fprintf('  - Spectral radius: %.2f\n', spectral_radius);
fprintf('  - Adaptation timescales (E): %d\n', n_a_E);
fprintf('  - Adaptation timescales (I): %d\n', n_a_I);
fprintf('  - STD enabled: E neurons only\n\n');

fprintf('Performance Metrics:\n');
fprintf('  Training:   MSE = %.6f, NRMSE = %.4f\n', metrics.train_mse, metrics.train_nrmse);
fprintf('  Validation: MSE = %.6f, NRMSE = %.4f\n', metrics.val_mse, metrics.val_nrmse);
fprintf('  Test:       MSE = %.6f, NRMSE = %.4f\n', test_metrics.mse, test_metrics.nrmse);
fprintf('=================================================================\n');

% Additional scatter plot for test predictions
figure('Position', [150, 150, 600, 500]);
scatter(Y_test, Y_test_pred, 10, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot([min(Y_test), max(Y_test)], [min(Y_test), max(Y_test)], 'r--', 'LineWidth', 2);
xlabel('Target Value');
ylabel('Predicted Value');
title(sprintf('Test Set: Target vs Predicted (NRMSE=%.4f)', test_metrics.nrmse));
grid on;
axis equal;
xlim([min(Y_test), max(Y_test)]);
ylim([min(Y_test), max(Y_test)]);
legend('Predictions', 'Perfect Fit', 'Location', 'best');
set(gcf, 'Color', 'w');

fprintf('\nExample completed successfully!\n');

