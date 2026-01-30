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

clear; clc;

fprintf('=================================================================\n');
fprintf('  Fractional Reservoir ESN: Mackey-Glass Time Series Prediction\n');
fprintf('=================================================================\n\n');

%% 1. Generate Mackey-Glass Time Series
fprintf('Step 1: Generating Mackey-Glass time series...\n');

% Generate chaotic time series (tau=17 for chaotic regime)
tau_mackey = 17;
dt = 1;  % Time step (seconds) - used for both data sampling and reservoir integration
[t, x_mg] = generate_mackey_glass('tau', tau_mackey, ...
                                   'n_samples', 9000, ...
                                   'dt', dt, ...
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
n = 100;                    % Total number of neurons % often the best looking chaos starts art n = 100
fraction_E = 0.5;           % 80% excitatory
n_E = round(n * fraction_E);
n_I = n - n_E;

fprintf('  Reservoir size: %d neurons (%d E, %d I)\n', n, n_E, n_I);

% Fractional-order parameters: Multiple adaptation timescales
n_a_E = 3;  % 3 adaptation timescales for E neurons
n_a_I = 2;  % 2 adaptation timescales for I neurons % this makes I neurons have some adaptation, which may be good
tau_a_E = logspace(log10(0.25), log10(25), n_a_E);  % Fast, medium, slow adaptation,  if n_a_E == 1, then only 25 is used.  
tau_a_I = logspace(log10(0.25), log10(25), n_a_I);

fprintf('  Adaptation timescales (E): [%s]\n', num2str(tau_a_E));
fprintf('  Adaptation timescales (I): [%s]\n', num2str(tau_a_I));

% Short-term depression parameters
n_b_E = 0;  % Enable STD for E neurons  % can set to zero for less stability
n_b_I = 0;  % No STD for I neurons
tau_b_E_rec = 2;   % Recovery time constant
tau_b_E_rel = 0.5;   % Release time constant

if n_b_E ~= 0
    fprintf('  STD enabled for E neurons (tau_rec=%.0f, tau_rel=%.0f)\n', ...
        tau_b_E_rec, tau_b_E_rel);
else
    fprintf('  STD disabled for E neurons\n');
end

if n_b_I ~=0
    fprintf('  STD enabled for I neurons (tau_rec=%.0f, tau_rel=%.0f)\n', ...
            tau_b_I_rec, tau_b_I_rel);
else
    fprintf('  STD disabled for I neurons\n');
end

% Dendritic time constant
tau_d = 1;

% Adaptation scaling factors
c_E = 0.1*1/3;  % Moderate adaptation effect for E neurons % you can make this lower to allow more "free" dynamics
c_I = 0.1;  % Weaker adaptation for I neurons

% Create recurrent weight matrix W with spectral radius scaling
  % For reproducibility, % because 42
W = randn(n, n);

% Dale's law: E neurons have positive weights, I neurons negative
W(:, 1:n_E) = abs(W(:, 1:n_E));      % E columns positive
W(:, n_E+1:end) = -abs(W(:, n_E+1:end));  % I columns negative 
% with this construction of W, the random matrix theory of expected spectral radius won't necessarily work
% but since we are computing the actual abscissa of W, it all should work out fine.
% my create_W_matrix(params) might improve on this
W = W-mean(W,2); % this ensures that all rows sum to zero 
% which reduces the outlier eigs of W, and therefore "richer" dynamics which are not dominated by a single real eig or pair of eigs
% Create input weight matrix with uniform random scaling
input_scaling = 0.5;
n_inputs = size(U, 2);
W_in = (2 * rand(n, n_inputs) - 1) * input_scaling;
W_in(rand(n, n_inputs) > 0.8) = 0;  % Make W_in sparse with X% of neurons getting input

fprintf('  Input scaling: %.2f\n', input_scaling);

% Activation function (ReLU-like with threshold)
a_0=0.5; % bias of the activation function.
%activation_function = @(x) min(max(0, x-a_0),1); % hard sigmoid with bias to middle. Sigmoid is much better than relu for chaos, and then you don't need STD to ensure stability
activation_function = @(x) tanh(x)+1;
% if the sigmoid is biased to the center, then the level_of_chaos is true.  
% if not biased to the center, then actual level of chaos is lower due to rectification or saturation

%level of chaos
level_of_chaos = 1.5; % this can go much higher up if adaptation and depression are used.  1.0 assumes no adaptation or depression
W_eigs = eig(W);
abscissa_0 = max(real(W_eigs)); 
gamma = 1 / abscissa_0; % calcualtes the adjustment to make 

% Pack parameters
params = struct();
params.n = n;
params.n_E = n_E;
params.n_I = n_I;
params.W = level_of_chaos * gamma * W; % this looks good.  gamma corrects W to have absciss of 1, wand then level_of_chaos pushes it into chaos
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
params.dt = dt;  % Time step for ODE integration (matches Mackey-Glass sampling rate)

% Create ESN
esn = SRNN_ESN(params);

fprintf('  Feature extraction: using dendritic states (x)\n');
fprintf('  Include raw input: %s\n', mat2str(params.include_input));
fprintf('  Integration time step (dt): %.2f seconds\n\n', params.dt);

%% 4. Train the Readout Layer
fprintf('Step 4: Training readout layer...\n');

% Training configuration
train_options = struct();
train_options.train_ratio = 0.4;      % 60% for training
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
figure;

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

%% 6b. 3D Attractor Visualizations (Delay Embedding)
fprintf('Step 6b: Creating 3D attractor visualizations...\n');

% Delay for embedding (matches Mackey-Glass generation parameter tau)
tau_embed = tau_mackey;

% Helper function to create delay embedding coordinates
create_embedding = @(signal, tau) [...
    signal(1:end-2*tau), ...
    signal(1+tau:end-tau), ...
    signal(1+2*tau:end)];

% Determine number of samples to plot (avoid cluttered visualization)
n_samples_3d_train = min(1000, length(Y_train_full) - 2*tau_embed);
n_samples_3d_val = min(1000, length(Y_val) - 2*tau_embed);
n_samples_3d_test = min(1000, length(Y_test) - 2*tau_embed);

% Create embeddings for training set
Y_train_embed = create_embedding(Y_train_full(1:n_samples_3d_train+2*tau_embed), tau_embed);
Y_train_pred_embed = create_embedding(Y_train_pred(1:n_samples_3d_train+2*tau_embed), tau_embed);

% Create embeddings for validation set
Y_val_embed = create_embedding(Y_val(1:n_samples_3d_val+2*tau_embed), tau_embed);
Y_val_pred_embed = create_embedding(Y_val_pred(1:n_samples_3d_val+2*tau_embed), tau_embed);

% Create embeddings for test set
Y_test_embed = create_embedding(Y_test(1:n_samples_3d_test+2*tau_embed), tau_embed);
Y_test_pred_embed = create_embedding(Y_test_pred(1:n_samples_3d_test+2*tau_embed), tau_embed);

% Create figure with 3D attractor plots
figure;

% Subplot 1: Training Set Attractor
subplot(1, 3, 1);
plot3(Y_train_embed(:,1), Y_train_embed(:,2), Y_train_embed(:,3), ...
      'b-', 'LineWidth', 1.2);
hold on;
plot3(Y_train_pred_embed(:,1), Y_train_pred_embed(:,2), Y_train_pred_embed(:,3), ...
      'r-', 'LineWidth', 1.0, 'Color', [1, 0.5, 0]);
xlabel(sprintf('x(t)'));
ylabel(sprintf('x(t+%d)', tau_embed));
zlabel(sprintf('x(t+%d)', 2*tau_embed));
title('Training Set Attractor');
legend('Ground Truth', 'Predicted', 'Location', 'best');
grid on;
view(-37.5, 30);  % Good viewing angle for 3D structure
axis tight;

% Subplot 2: Validation Set Attractor
subplot(1, 3, 2);
plot3(Y_val_embed(:,1), Y_val_embed(:,2), Y_val_embed(:,3), ...
      'b-', 'LineWidth', 1.2);
hold on;
plot3(Y_val_pred_embed(:,1), Y_val_pred_embed(:,2), Y_val_pred_embed(:,3), ...
      'r-', 'LineWidth', 1.0, 'Color', [1, 0.5, 0]);
xlabel(sprintf('x(t)'));
ylabel(sprintf('x(t+%d)', tau_embed));
zlabel(sprintf('x(t+%d)', 2*tau_embed));
title('Validation Set Attractor');
legend('Ground Truth', 'Predicted', 'Location', 'best');
grid on;
view(-37.5, 30);
axis tight;

% Subplot 3: Test Set Attractor
subplot(1, 3, 3);
plot3(Y_test_embed(:,1), Y_test_embed(:,2), Y_test_embed(:,3), ...
      'b-', 'LineWidth', 1.2);
hold on;
plot3(Y_test_pred_embed(:,1), Y_test_pred_embed(:,2), Y_test_pred_embed(:,3), ...
      'r-', 'LineWidth', 1.0, 'Color', [1, 0.5, 0]);
xlabel(sprintf('x(t)'));
ylabel(sprintf('x(t+%d)', tau_embed));
zlabel(sprintf('x(t+%d)', 2*tau_embed));
title('Test Set Attractor');
legend('Ground Truth', 'Predicted', 'Location', 'best');
grid on;
view(-37.5, 30);
axis tight;

% Overall title
sgtitle(sprintf('3D Attractor Reconstruction (Delay Embedding, \\tau=%d)', tau_embed), ...
        'FontSize', 14, 'FontWeight', 'bold');

% Adjust layout
set(gcf, 'Color', 'w');

fprintf('  3D attractor visualization complete!\n');
fprintf('  Embedding delay: %d time steps\n', tau_embed);
fprintf('  Points plotted - Train: %d, Val: %d, Test: %d\n\n', ...
        n_samples_3d_train, n_samples_3d_val, n_samples_3d_test);

%% 7. Generative Mode (Closed-Loop Autonomous Prediction)
fprintf('=================================================================\n');
fprintf('Step 7: Testing generative mode (closed-loop prediction)...\n');
fprintf('=================================================================\n\n');

% 7a. Configuration
fprintf('Step 7a: Configuring generative mode parameters...\n');

gen_horizon = 1;  % True closed-loop: feed each prediction back immediately
gen_length = n_test;  % Generate same length as test set (default)
gen_washout = 100;  % Washout steps before generation

fprintf('  Prediction horizon: %d\n', gen_horizon);
fprintf('  Generation length: %d steps\n', gen_length);
fprintf('  Washout steps: %d\n\n', gen_washout);

% 7b. Generation from Validation End
fprintf('Step 7b: Generating from end of validation set...\n');

% Prepare initial data: last gen_washout steps of validation
U_init_val = U(max(1, n_train+n_val-gen_washout+1):n_train+n_val, :);

% Reset and generate
esn.resetState();
gen_options = struct();
gen_options.horizon = gen_horizon;
gen_options.washout_steps = size(U_init_val, 1);

tic;
Y_gen_from_val = esn.generateAutonomous(U_init_val, gen_length, gen_options);
gen_time_val = toc;

fprintf('  Generated %d steps in %.2f seconds\n', gen_length, gen_time_val);
fprintf('  Generation rate: %.1f steps/sec\n\n', gen_length/gen_time_val);

% 7c. Generation from Test Set Start (with ground truth comparison)
fprintf('Step 7c: Generating from test set start (with ground truth)...\n');

% Prepare initial data: first gen_washout steps of test set
U_init_test = U_test(1:min(gen_washout, size(U_test, 1)), :);
gen_length_test = min(gen_length, n_test - gen_washout);

% Reset and generate
esn.resetState();
gen_options_test = struct();
gen_options_test.horizon = gen_horizon;
gen_options_test.washout_steps = size(U_init_test, 1);

tic;
Y_gen_from_test = esn.generateAutonomous(U_init_test, gen_length_test, gen_options_test);
gen_time_test = toc;

% Get corresponding ground truth (after washout)
Y_test_true = Y_test(gen_washout+1:gen_washout+gen_length_test, :);

fprintf('  Generated %d steps in %.2f seconds\n', gen_length_test, gen_time_test);
fprintf('  Generation rate: %.1f steps/sec\n\n', gen_length_test/gen_time_test);

% 7d. Divergence Analysis
fprintf('Step 7d: Analyzing divergence from ground truth...\n');

% Compute point-wise divergence
divergence = abs(Y_gen_from_test - Y_test_true);
mean_divergence = mean(divergence);
std_divergence = std(divergence);

% Compute correlation
correlation = corr(Y_gen_from_test(:), Y_test_true(:));

% Find time to critical divergence (e.g., when error exceeds 2*std of signal)
critical_threshold = 2 * std(Y_test_true(:));
critical_indices = find(divergence > critical_threshold, 1, 'first');
if isempty(critical_indices)
    time_to_critical = gen_length_test;
else
    time_to_critical = critical_indices;
end

fprintf('  Mean divergence: %.6f\n', mean_divergence);
fprintf('  Std divergence: %.6f\n', std_divergence);
fprintf('  Correlation: %.6f\n', correlation);
fprintf('  Time to critical divergence: %d steps (%.1f%%)\n\n', ...
        time_to_critical, 100*time_to_critical/gen_length_test);

% 7e. Attractor Comparison
fprintf('Step 7e: Creating attractor comparison...\n');

% Create delay embeddings for generated sequences
n_samples_attractor = min(1000, gen_length_test - 2*tau_embed);

Y_gen_embed = create_embedding(Y_gen_from_test(1:n_samples_attractor+2*tau_embed), tau_embed);
Y_true_embed = create_embedding(Y_test_true(1:n_samples_attractor+2*tau_embed), tau_embed);

fprintf('  Created delay embeddings with tau=%d\n', tau_embed);
fprintf('  Using %d samples for attractor visualization\n\n', n_samples_attractor);

% 7f. Lyapunov Analysis
fprintf('Step 7f: Computing Lyapunov exponents...\n');

% Prepare parameters for Lyapunov analysis
% Note: We need to prepare the data in a format suitable for benettin_algorithm
% This is a simplified approach - for full analysis, we'd need the actual dynamics

% For time series data, we estimate LLE using nearest-neighbor divergence
% or use the benettin algorithm if we have access to the dynamics
% Here we'll compute a simplified metric based on local divergence rates

% Compute local exponential divergence rate
window_size = 50;  % Window for local divergence estimation
n_windows = floor(gen_length_test / window_size);

local_div_rates = zeros(n_windows, 1);
for w = 1:n_windows
    idx_start = (w-1)*window_size + 1;
    idx_end = min(w*window_size, gen_length_test);
    
    local_error = norm(Y_gen_from_test(idx_start:idx_end) - Y_test_true(idx_start:idx_end));
    initial_sep = norm(Y_gen_from_test(idx_start) - Y_test_true(idx_start)) + 1e-10;
    
    local_div_rates(w) = log(local_error / initial_sep) / (window_size);
end

mean_div_rate = mean(local_div_rates(isfinite(local_div_rates)));

fprintf('  Mean local divergence rate: %.6f\n', mean_div_rate);
fprintf('  (Positive = trajectories diverging, typical for chaotic systems)\n\n');

% 7g. Spectral Analysis
fprintf('Step 7g: Performing spectral analysis...\n');

% Compute power spectral density using FFT
n_fft = 2^nextpow2(gen_length_test);
fs_signal = 1/dt;  % Sampling frequency in Hz (inverse of time step)

% True signal spectrum
Y_test_fft = fft(Y_test_true, n_fft);
psd_true = abs(Y_test_fft(1:n_fft/2+1)).^2 / gen_length_test;
freq = fs_signal * (0:(n_fft/2)) / n_fft;

% Generated signal spectrum
Y_gen_fft = fft(Y_gen_from_test, n_fft);
psd_gen = abs(Y_gen_fft(1:n_fft/2+1)).^2 / gen_length_test;

% Find dominant frequencies
[~, idx_true] = max(psd_true(2:end));  % Skip DC component
[~, idx_gen] = max(psd_gen(2:end));
dominant_freq_true = freq(idx_true+1);
dominant_freq_gen = freq(idx_gen+1);

fprintf('  Dominant frequency (true): %.6f Hz\n', dominant_freq_true);
fprintf('  Dominant frequency (generated): %.6f Hz\n', dominant_freq_gen);
fprintf('  Frequency error: %.6f Hz (%.2f%%)\n\n', ...
        abs(dominant_freq_true - dominant_freq_gen), ...
        100*abs(dominant_freq_true - dominant_freq_gen)/dominant_freq_true);

% Visualization: Generative Mode Results
fprintf('Step 7h: Creating generative mode visualizations...\n');

% Figure 1: Time Series Comparison
figure;

% Subplot 1: Generated vs True (full sequence)
subplot(3, 2, 1);
t_gen = (1:gen_length_test)';
plot_samples = min(500, gen_length_test);
plot(t_gen(1:plot_samples), Y_test_true(1:plot_samples), 'b-', 'LineWidth', 1.5);
hold on;
plot(t_gen(1:plot_samples), Y_gen_from_test(1:plot_samples), 'r--', 'LineWidth', 1);
xlabel('Time Step');
ylabel('Value');
title('Generated vs True (First 500 steps)');
legend('Ground Truth', 'Generated', 'Location', 'best');
grid on;

% Subplot 2: Divergence over time
subplot(3, 2, 2);
plot(t_gen, divergence, 'k-', 'LineWidth', 1);
hold on;
yline(critical_threshold, 'r--', 'LineWidth', 1.5, 'Label', 'Critical Threshold');
if ~isempty(critical_indices)
    xline(time_to_critical, 'r:', 'LineWidth', 1.5);
end
xlabel('Time Step');
ylabel('Absolute Error');
title(sprintf('Divergence (Mean=%.4f)', mean_divergence));
grid on;

% Subplot 3: Correlation over time (rolling window)
subplot(3, 2, 3);
window = 100;
rolling_corr = zeros(gen_length_test - window + 1, 1);
for i = 1:(gen_length_test - window + 1)
    rolling_corr(i) = corr(Y_gen_from_test(i:i+window-1), Y_test_true(i:i+window-1));
end
plot(window:gen_length_test, rolling_corr, 'b-', 'LineWidth', 1);
xlabel('Time Step');
ylabel('Correlation');
title('Rolling Correlation (window=100)');
ylim([0, 1]);
grid on;

% Subplot 4: Local divergence rates
subplot(3, 2, 4);
t_windows = (0:n_windows-1) * window_size + window_size/2;
plot(t_windows, local_div_rates, 'k-', 'LineWidth', 1);
hold on;
yline(0, 'r--', 'LineWidth', 1);
xlabel('Time Step');
ylabel('Local Divergence Rate');
title(sprintf('Local Divergence Rate (Mean=%.4f)', mean_div_rate));
grid on;

% Subplot 5: Scatter plot
subplot(3, 2, 5);
scatter(Y_test_true, Y_gen_from_test, 10, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot([min(Y_test_true), max(Y_test_true)], [min(Y_test_true), max(Y_test_true)], ...
     'r--', 'LineWidth', 2);
xlabel('True Value');
ylabel('Generated Value');
title(sprintf('Generated vs True (r=%.4f)', correlation));
axis equal;
grid on;

% Subplot 6: Power Spectral Density
subplot(3, 2, 6);
semilogy(freq, psd_true, 'b-', 'LineWidth', 1.5);
hold on;
semilogy(freq, psd_gen, 'r--', 'LineWidth', 1);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density');
title('Spectral Comparison');
legend('Ground Truth', 'Generated', 'Location', 'best');
xlim([0, 0.1]);  % Focus on low frequencies
grid on;

sgtitle('Generative Mode: Performance Analysis', 'FontSize', 14, 'FontWeight', 'bold');
set(gcf, 'Color', 'w');

% Figure 2: 3D Attractor Comparison
figure;

subplot(1, 2, 1);
plot3(Y_true_embed(:,1), Y_true_embed(:,2), Y_true_embed(:,3), ...
      'b-', 'LineWidth', 1.2);
hold on;
plot3(Y_gen_embed(:,1), Y_gen_embed(:,2), Y_gen_embed(:,3), ...
      'r--', 'LineWidth', 1.0);
xlabel(sprintf('x(t)'));
ylabel(sprintf('x(t+%d)', tau_embed));
zlabel(sprintf('x(t+%d)', 2*tau_embed));
title('3D Attractor: Generated vs True');
legend('Ground Truth', 'Generated', 'Location', 'best');
grid on;
view(-37.5, 30);
axis tight;

% Overlay plot for better comparison
subplot(1, 2, 2);
plot3(Y_true_embed(:,1), Y_true_embed(:,2), Y_true_embed(:,3), ...
      'b-', 'LineWidth', 1.5, 'Color', [0, 0, 1, 0.6]);
hold on;
plot3(Y_gen_embed(:,1), Y_gen_embed(:,2), Y_gen_embed(:,3), ...
      'r-', 'LineWidth', 1.5, 'Color', [1, 0, 0, 0.6]);
xlabel(sprintf('x(t)'));
ylabel(sprintf('x(t+%d)', tau_embed));
zlabel(sprintf('x(t+%d)', 2*tau_embed));
title('3D Attractor Overlay');
legend('Ground Truth', 'Generated', 'Location', 'best');
grid on;
view(-37.5, 30);
axis tight;

sgtitle(sprintf('Attractor Reconstruction (Delay Embedding, \\tau=%d)', tau_embed), ...
        'FontSize', 14, 'FontWeight', 'bold');
set(gcf, 'Color', 'w');

fprintf('  Generative mode visualization complete!\n\n');

%% 8. Summary
fprintf('=================================================================\n');
fprintf('  SUMMARY\n');
fprintf('=================================================================\n');
fprintf('Reservoir Configuration:\n');
fprintf('  - Size: %d neurons (%d E, %d I)\n', n, n_E, n_I);
fprintf('  - Adaptation timescales (E): %d\n', n_a_E);
fprintf('  - Adaptation timescales (I): %d\n', n_a_I);
fprintf('  - STD enabled: E neurons only\n\n');

fprintf('Performance Metrics (Teacher Forcing):\n');
fprintf('  Training:   MSE = %.6f, NRMSE = %.4f\n', metrics.train_mse, metrics.train_nrmse);
fprintf('  Validation: MSE = %.6f, NRMSE = %.4f\n', metrics.val_mse, metrics.val_nrmse);
fprintf('  Test:       MSE = %.6f, NRMSE = %.4f\n\n', test_metrics.mse, test_metrics.nrmse);

fprintf('Generative Mode Performance (Closed-Loop):\n');
fprintf('  Generation length: %d steps (horizon=%d)\n', gen_length_test, gen_horizon);
fprintf('  Mean divergence: %.6f\n', mean_divergence);
fprintf('  Std divergence: %.6f\n', std_divergence);
fprintf('  Correlation with ground truth: %.6f\n', correlation);
fprintf('  Time to critical divergence: %d steps (%.1f%% of sequence)\n', ...
        time_to_critical, 100*time_to_critical/gen_length_test);
fprintf('  Mean local divergence rate: %.6f\n', mean_div_rate);
fprintf('  Dominant frequency (true): %.6f Hz\n', dominant_freq_true);
fprintf('  Dominant frequency (generated): %.6f Hz\n', dominant_freq_gen);
fprintf('  Spectral error: %.2f%%\n', ...
        100*abs(dominant_freq_true - dominant_freq_gen)/dominant_freq_true);
fprintf('=================================================================\n');

% Additional scatter plot for test predictions
figure;
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

