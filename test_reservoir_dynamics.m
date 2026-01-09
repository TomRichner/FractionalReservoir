%% Test Script: Reservoir Dynamics Visualization
%
% This script tests the dynamical behavior of the fractional SRNN reservoir
% by driving it with different input signals and visualizing neuron states.
%
% Purpose: Explore and tune reservoir parameters to observe different
%          dynamical regimes (stable, chaotic, oscillatory, etc.)
%
% Usage: Adjust parameters in Section 1, then run to see dynamics

clear all; 
clc;

fprintf('=================================================================\n');
fprintf('  Testing Reservoir Dynamics - Parameter Exploration\n');
fprintf('=================================================================\n\n');

%% 1. ADJUSTABLE PARAMETERS
fprintf('Step 1: Setting up parameters...\n');

% === Reservoir Architecture ===
n = 210;                      % Total number of neurons (try: 20-200)
fraction_E = 0.5;            % Fraction of excitatory neurons (try: 0.3-0.5)

% === Time Parameters ===
dt = 0.1;                    % Integration time step (seconds)
n_steps = 300;               % Number of time steps to stimulate

% === Adaptation Parameters ===
n_a_E = 3;                   % Number of adaptation timescales for E (try: 0-3)
n_a_I = 1;                   % Number of adaptation timescales for I (try: 0-3)
tau_a_E = logspace(log10(0.25), log10(25), n_a_E);  % E adaptation timescales
tau_a_I = logspace(log10(0.25), log10(25), n_a_I);  % I adaptation timescales
c_E = 0.1/7;                 % Adaptation strength for E (try: 0.0-0.5)
c_I = 0.1/4;                   % Adaptation strength for I (try: 0.0-0.5)

% === Short-term Depression Parameters ===
n_b_E = 1;                   % Enable STD for E neurons (0=off, 1=on)
n_b_I = 1;                   % Enable STD for I neurons (0=off, 1=on)
tau_b_E_rec = 0.6;           % Recovery time constant
tau_b_E_rel = 0.1;           % Release time constant
tau_b_I_rec = 0.4;           
tau_b_I_rel = 0.5;  

% === Network Dynamics ===
tau_d = 0.55;                % Dendritic time constant
level_of_chaos = 3.7;        % Chaos level (try: 0.5-2.5) - higher = more chaotic
lags = 0.03; %or 0.04
% === Activation Function ===
%a_0 = 0.000000001;                   % Activation bias (0.5 = centered)
%activation_function = @(x) min(max(0, x-a_0), 1);  % Hard sigmoid
%activation_function = @(x) 0.5*tanh(2*x)+1;
S_a = 0.85;
S_c = 0.4;
activation_function = @(x) piecewiseSigmoid(x, S_a, S_c);

% === Input Configuration ===
input_scaling = 0.75;         % Input weight scaling (try: 0.1-1.0)
input_type = 'pulse';         % Options: 'sine', 'pulse', 'step', 'noise', 'mackey'
input_amplitude = 1;       % Input amplitude

fprintf('  Reservoir: %d neurons (%.0f%% E, %.0f%% I)\n', ...
        n, fraction_E*100, (1-fraction_E)*100);
fprintf('  Adaptation: E(%d scales), I(%d scales)\n', n_a_E, n_a_I);
fprintf('  Chaos level: %.2f\n', level_of_chaos);
fprintf('  Input: %s, scaling=%.2f\n\n', input_type, input_scaling);

%% 2. Generate Input Signal
fprintf('Step 2: Generating input signal...\n');

t = (0:n_steps-1)' * dt;

% Define silence duration in seconds
N = 10; % Duration before stimulus in seconds
M = 200; % Duration after stimulus in seconds
n_silence_steps_start = round(N / dt);
n_silence_steps_end = round(M / dt);
total_steps = n_steps + n_silence_steps_start + n_silence_steps_end;

% Initialize input signal with silence
U = zeros(total_steps, 1);

% Generate input signal based on the specified type
switch lower(input_type)
    case 'sine'
        % Sinusoidal input
        input_frequency = 0.001;      % Hz
        U(n_silence_steps_start+1:n_silence_steps_start+n_steps) = ...
            input_amplitude * sin(2*pi*input_frequency*t);
        
    case 'pulse'
        % Periodic pulses
        pulse_period = 20;
        pulse_width = 5;
        pulse_signal = zeros(n_steps, 1);
        pulse_signal(mod(0:n_steps-1, pulse_period) < pulse_width) = input_amplitude;
        U(n_silence_steps_start+1:n_silence_steps_start+n_steps) = pulse_signal;
        
    case 'step'
        % Step input at quarter-way
        step_signal = zeros(n_steps, 1);
        step_signal(round(n_steps/4):end) = input_amplitude;
        U(n_silence_steps_start+1:n_silence_steps_start+n_steps) = step_signal;
        
    case 'noise'
        % White noise
        U(n_silence_steps_start+1:n_silence_steps_start+n_steps) = ...
            input_amplitude * randn(n_steps, 1);
        
    case 'mackey'
        % Mackey-Glass chaotic signal
        [~, mg] = generate_mackey_glass('tau', 17, 'n_samples', n_steps+1000, ...
                                         'dt', dt, 'discard', 1000);
        U(n_silence_steps_start+1:n_silence_steps_start+n_steps) = ...
            input_amplitude * mg(1:n_steps);
        
    otherwise
        error('Unknown input_type: %s', input_type);
end

% Make it a column vector and proper dimensionality
U = U(:);
t = (0:length(U)-1)' * dt;
n_inputs = size(U, 2);
fprintf('  Generated %d steps of %s input\n', total_steps, input_type);
fprintf('  Time range: [%.2f, %.2f] seconds\n\n', t(1), t(end));

%% 3. Build Reservoir
fprintf('Step 3: Building reservoir...\n');

% Network structure
n_E = round(n * fraction_E);
n_I = n - n_E;

% Create recurrent weight matrix
rng(42);  % Reproducibility
W = randn(n, n);

% Apply Dale's law
W(:, 1:n_E) = abs(W(:, 1:n_E));           % E columns positive
W(:, n_E+1:end) = -abs(W(:, n_E+1:end));  % I columns negative

% Center rows (helps create richer dynamics)
W = W - mean(W, 2);


% Apply chaos scaling
W_eigs = eig(W);
abscissa_0 = max(real(W_eigs));
gamma = 1 / abscissa_0;
W = level_of_chaos * gamma * W;

fprintf('  Initial abscissa: %.3f\n', abscissa_0);

% Create input weight matrix
W_in = (2 * rand(n, n_inputs) - 1) * input_scaling;
W_in(rand(n, n_inputs) > 0.2) = 0;  % Make W_in sparse with X% of neurons getting input

fprintf('  Input weights: %d x %d\n\n', n, n_inputs);
%% 4. Pack Parameters and Create ESN
fprintf('Step 4: Initializing SRNN-ESN...\n');

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
params.tau_b_I_rec = tau_b_I_rec;
params.tau_b_I_rel = tau_b_I_rel;
params.c_E = c_E;
params.c_I = c_I;
params.lags = lags;
params.activation_function = activation_function;
params.activation_function_derivative = @(x) piecewiseSigmoidDerivative(x, S_a, S_c);
params.E_indices = 1:n_E;
params.I_indices = (n_E+1):n;
params.which_states = 'x';
params.include_input = false;
params.lambda = 1e-6;
params.dt = dt;

esn = SRNN_ESN(params);

fprintf('  SRNN-ESN initialized\n\n');

%% 5. Drive Reservoir with Input (No Training)
fprintf('Step 5: Driving reservoir with input...\n');

% Reset state
esn.resetState();

% Run reservoir and get full state history
[~, S_history] = esn.runReservoir(U);

% Extract dendritic states x from state history
len_a_E = n_E * n_a_E;
len_a_I = n_I * n_a_I;
len_b_E = n_E * n_b_E;
len_b_I = n_I * n_b_I;
x_start_idx = len_a_E + len_a_I + len_b_E + len_b_I + 1;
x_end_idx = x_start_idx + n - 1;
states_x = S_history(:, x_start_idx:x_end_idx); % (n_steps x n)

% Compute firing rates r from states
% r = activation_function(x_eff) where x_eff = x - adaptation + depression
states_r = zeros(total_steps, n);

for i = 1:total_steps
    % Extract state components for this timestep
    current_idx = 1;
    
    % Extract adaptation variables a_E
    if n_a_E > 0
        len_a_E = n_E * n_a_E;
        a_E = reshape(S_history(i, current_idx:current_idx+len_a_E-1), n_a_E, n_E)';
        current_idx = current_idx + len_a_E;
    else
        a_E = zeros(n_E, 0);
    end
    
    % Extract adaptation variables a_I
    if n_a_I > 0
        len_a_I = n_I * n_a_I;
        a_I = reshape(S_history(i, current_idx:current_idx+len_a_I-1), n_a_I, n_I)';
        current_idx = current_idx + len_a_I;
    else
        a_I = zeros(n_I, 0);
    end
    
    % Skip depression variables (not needed for rate computation in this model)
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
    
    % Apply activation function to get firing rates
    states_r(i, :) = activation_function(x_eff);
end

fprintf('  Simulation complete: %d time steps\n\n', total_steps);

%% 6. Analyze Dynamics
fprintf('Step 6: Analyzing dynamics...\n');

% Compute statistics
mean_activity = mean(states_r, 1);
std_activity = std(states_r, 1);
max_activity = max(states_r, [], 1);

fprintf('  Mean firing rate: %.4f (std=%.4f)\n', mean(mean_activity), mean(std_activity));
fprintf('  Max firing rate: %.4f\n', max(max_activity));

% Compute dimensionality (participation ratio)
C = cov(states_r);
eigvals = eig(C);
eigvals = eigvals(eigvals > 1e-10);  % Remove numerical zeros
participation_ratio = sum(eigvals)^2 / sum(eigvals.^2);

fprintf('  Effective dimensionality: %.2f / %d\n', participation_ratio, n);

% Compute autocorrelation of mean population activity
recovery_start_idx = n_silence_steps_start+n_steps;
pop_activity = mean(states_r(recovery_start_idx:end,:),2);
[acf, corr_lags] = autocorr(pop_activity, NumLags=2000);
acf = acf(corr_lags >= 0);

fprintf('  Population activity range: [%.4f, %.4f]\n\n', min(pop_activity), max(pop_activity));

%% 7. Visualization
fprintf('Step 7: Creating visualizations...\n');

figure('Color', 'w');

% Subplot 1: Input Signal
subplot(3, 3, 1);
plot(t, U, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Input');
title(sprintf('Input Signal (%s)', input_type));
grid on;

% Subplot 2: Neuron States (Raster-like heatmap)
subplot(3, 3, 2);
imagesc(t, 1:n, states_r');
colormap(hot);
colorbar;
xlabel('Time (s)');
ylabel('Neuron Index');
title('Firing Rates (All Neurons)');
set(gca, 'YDir', 'normal');

% Subplot 3: Neuron Traces
subplot(3, 3, 3);
plot(t, states_r(:, 1:n), 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Dendritic State');
title('Neurons Activity');
grid on;
%ylim([-1, 1]);

% Subplot 4: E vs I Population Activity
subplot(3, 3, 4);
E_pop = mean(states_r(:, 1:n_E), 2);
I_pop = mean(states_r(:, n_E+1:end), 2);
plot(t, E_pop, 'r-', 'LineWidth', 1.5);
hold on;
plot(t, I_pop, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Mean State');
title('E/I Population Activity');
legend('Excitatory', 'Inhibitory', 'Location', 'best');
grid on;

% Subplot 5: Mean Activity per Neuron
subplot(3, 3, 5);
bar(1:n, mean_activity);
hold on;
xline(n_E+0.5, 'r--', 'LineWidth', 2);
xlabel('Neuron Index');
ylabel('Mean State');
title('Average Activity per Neuron');
legend('', 'E/I boundary');
grid on;

% Subplot 6: Population Autocorrelation
subplot(3, 3, 6);
plot(corr_lags(corr_lags>=0)*dt, acf, 'k-', 'LineWidth', 1.5);
xlabel('Lag (s)');
ylabel('Autocorrelation');
title('Population Activity Autocorrelation');
grid on;

% Subplot 7: Phase Space (3D)
[coeff,score,latent] = pca(states_r(recovery_start_idx:end, :),'NumComponents',3);
subplot(3, 3, 7);
plot3(score(:, 1), score(:, 2), score(:, 3), 'b-', 'LineWidth', 0.5);
title('3D Phase Space - Recovery Phase');
grid on;
view(-37.5, 30);

% Subplot 8: State Distribution
subplot(3, 3, 8);
histogram(states_r(:), 30, 'Normalization', 'probability');
xlabel('Firing Rate');
ylabel('Probability');
title('Distribution of Firing Rates');
grid on;

% Subplot 9: Activity Covariance Matrix
subplot(3, 3, 9);
imagesc(corr(states_r));
colormap(gca, 'jet');
colorbar;
xlabel('Neuron Index');
ylabel('Neuron Index');
title('Correlation Matrix');
axis square;

% Overall title
sgtitle(sprintf('Reservoir Dynamics Test | n=%d, chaos=%.2f, input=%s', ...
                n, level_of_chaos, input_type), ...
        'FontSize', 14, 'FontWeight', 'bold');

fprintf('  Visualization complete!\n\n');

%% 8. Jacobian Analysis at Multiple Timepoints
fprintf('Step 8: Computing Jacobian at 10 timepoints...\n');

% Select 10 evenly-spaced timepoints across the simulation
n_jacobian_times = 10;
J_times = round(linspace(1, total_steps, n_jacobian_times));

% Compute Jacobians at selected timepoints
J_array = compute_Jacobian_at_indices(S_history, J_times, params);

% Compute eigenvalues at each timepoint
n_eigs = size(J_array, 1);
all_eigs = zeros(n_eigs, n_jacobian_times);
spectral_abscissa = zeros(n_jacobian_times, 1);

for i = 1:n_jacobian_times
    eigs_i = eig(J_array(:,:,i));
    all_eigs(:,i) = eigs_i;
    spectral_abscissa(i) = max(real(eigs_i));
end

fprintf('  Jacobian dimensions: %d x %d\n', size(J_array,1), size(J_array,2));
fprintf('  Spectral abscissa range: [%.4f, %.4f]\n', min(spectral_abscissa), max(spectral_abscissa));

% Create Jacobian visualization figure
figure('Color', 'w', 'Name', 'Jacobian Analysis');

% Subplot 1: Eigenvalue spectra for all timepoints
subplot(1, 3, 1);
colors = parula(n_jacobian_times);
hold on;
for i = 1:n_jacobian_times
    scatter(real(all_eigs(:,i)), imag(all_eigs(:,i)), 15, colors(i,:), 'filled', 'MarkerFaceAlpha', 0.5);
end
xline(0, 'k--', 'LineWidth', 1);
xlabel('Real Part');
ylabel('Imaginary Part');
title('Eigenvalue Spectra (10 Timepoints)');
colormap(gca, parula);
cb = colorbar;
cb.Label.String = 'Time (s)';
clim([t(J_times(1)), t(J_times(end))]);
grid on;
axis equal;

% Subplot 2: Spectral abscissa over time
subplot(1, 3, 2);
plot(t(J_times), spectral_abscissa, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
hold on;
yline(0, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Spectral Abscissa');
title('Max Real Eigenvalue vs Time');
legend('Spectral Abscissa', 'Stability Boundary', 'Location', 'best');
grid on;

% Subplot 3: Jacobian heatmap at middle timepoint
subplot(1, 3, 3);
mid_idx = round(n_jacobian_times / 2);
J_mid = J_array(:,:,mid_idx);
imagesc(J_mid);
colormap(gca, bluewhitered_colormap(256));
colorbar;
xlabel('State Index');
ylabel('State Index');
title(sprintf('Jacobian at t=%.2fs', t(J_times(mid_idx))));
axis square;

sgtitle('Jacobian Analysis', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('  Jacobian analysis complete!\n\n');

%% 9. Summary
fprintf('=================================================================\n');
fprintf('  SUMMARY\n');
fprintf('=================================================================\n');
fprintf('Configuration:\n');
fprintf('  - Neurons: %d (%d E, %d I)\n', n, n_E, n_I);
fprintf('  - Chaos level: %.2f\n', level_of_chaos);
fprintf('  - Adaptation: E(%d), I(%d)\n', n_a_E, n_a_I);
fprintf('  - Input: %s (scaling=%.2f)\n', input_type, input_scaling);
fprintf('  - Simulation: %d steps (%.2f seconds)\n\n', total_steps, t(end));

fprintf('Dynamics:\n');
fprintf('  - Mean firing rate: %.4f\n', mean(mean_activity));
fprintf('  - Effective dimensionality: %.2f / %d (%.1f%%)\n', ...
        participation_ratio, n, 100*participation_ratio/n);
fprintf('  - Population activity range: [%.4f, %.4f]\n', ...
        min(pop_activity), max(pop_activity));
fprintf('=================================================================\n\n');

fprintf('Test script completed!\n');

%% 10. SPEECH RECOGNITION CLASSIFICATION TASK
% =========================================================================
% This section implements a speech command classification task using the
% fractional ESN. We use the Google Speech Commands dataset with raw audio
% waveforms (no spectrograms) to classify spoken words.
% =========================================================================

fprintf('\n=================================================================\n');
fprintf('  SECTION 10: Speech Recognition Classification\n');
fprintf('=================================================================\n\n');

% === 10.1 Dataset Download and Import ===
fprintf('Step 10.1: Downloading and importing Google Speech Commands dataset...\n');

% Dataset URL and local storage
datasetURL = 'https://storage.googleapis.com/download.tensorflow.org/data/speech_commands_v0.02.tar.gz';
datasetFolder = fullfile(tempdir, 'google_speech_commands');

% Download and extract if not already present
if ~exist(fullfile(datasetFolder, 'yes'), 'dir')
    fprintf('  Downloading dataset (this may take a few minutes)...\n');
    
    % Create folder
    if ~exist(datasetFolder, 'dir')
        mkdir(datasetFolder);
    end
    
    % Download
    tarFile = fullfile(tempdir, 'speech_commands.tar.gz');
    if ~exist(tarFile, 'file')
        websave(tarFile, datasetURL);
    end
    
    % Extract
    fprintf('  Extracting dataset...\n');
    untar(tarFile, datasetFolder);
    fprintf('  Dataset extracted to: %s\n', datasetFolder);
else
    fprintf('  Dataset already exists at: %s\n', datasetFolder);
end

% Define target classes (small subset for faster training)
targetClasses = categorical({'yes', 'no', 'up', 'down', 'stop', 'go'});
numClasses = numel(targetClasses);

fprintf('  Target classes: %s\n', strjoin(string(targetClasses), ', '));

% Create audioDatastore for each class
fprintf('  Loading audio files...\n');
ads = audioDatastore(datasetFolder, ...
    'IncludeSubfolders', true, ...
    'LabelSource', 'foldernames');

% Filter to only target classes
isTargetClass = ismember(ads.Labels, targetClasses);
ads = subset(ads, isTargetClass);

% Remove unused categories from the labels (important for correct counting)
ads.Labels = removecats(ads.Labels);

fprintf('  Total samples after filtering: %d\n', numel(ads.Files));

% Display class distribution (now only shows target classes)
uniqueLabels = categories(ads.Labels);
labelCounts = zeros(numel(uniqueLabels), 1);
fprintf('  Class distribution:\n');
for i = 1:numel(uniqueLabels)
    labelCounts(i) = sum(ads.Labels == uniqueLabels{i});
    fprintf('    %s: %d samples\n', uniqueLabels{i}, labelCounts(i));
end

% === 10.2 Audio Preprocessing ===
fprintf('\nStep 10.2: Preprocessing audio data...\n');

% Audio parameters - use frame-based features to match reservoir dynamics
originalFs = 16000;    % Original sample rate (16 kHz)
intermediateSamples = 16000;  % Keep 1 second of audio at original rate first

% Frame-based processing to match reservoir timescale (tau_d = 0.55s, dt = 0.1s)
nFrames = 10;          % Number of frames (matches ~1 second with dt=0.1)
dt_speech = dt;        % Use same dt as reservoir (0.1s)

fprintf('  Original sample rate: %d Hz\n', originalFs);
fprintf('  Frame-based processing: %d frames (dt=%.2fs, total=%.1fs)\n', nFrames, dt_speech, nFrames*dt_speech);
fprintf('  This matches reservoir dynamics (tau_d=%.2fs)\n', tau_d);

% Balance the dataset by subsampling majority classes
fprintf('  Balancing dataset...\n');
minSamplesPerClass = min(labelCounts);
maxSamplesPerClass = min(minSamplesPerClass, 500);  % Limit for faster training

fprintf('  Using %d samples per class\n', maxSamplesPerClass);

% Create balanced subset
balancedIndices = [];
for i = 1:numel(uniqueLabels)
    classIndices = find(ads.Labels == uniqueLabels{i});
    % Randomly select maxSamplesPerClass samples
    rng(42);  % Reproducibility
    nToSelect = min(maxSamplesPerClass, numel(classIndices));
    selectedIndices = classIndices(randperm(numel(classIndices), nToSelect));
    balancedIndices = [balancedIndices; selectedIndices];
end

% Shuffle the balanced indices
rng(42);
balancedIndices = balancedIndices(randperm(numel(balancedIndices)));
adsBalanced = subset(ads, balancedIndices);

fprintf('  Balanced dataset size: %d samples\n', numel(adsBalanced.Files));

% Load and preprocess all audio files into frame-based features
fprintf('  Loading and preprocessing audio files...\n');
numSamples = numel(adsBalanced.Files);
audioFrames = zeros(numSamples, nFrames);  % Frame-based RMS amplitude
audioRaw = cell(numSamples, 1);  % Keep raw audio for visualization
labels = adsBalanced.Labels;

for i = 1:numSamples
    if mod(i, 500) == 0
        fprintf('    Processed %d/%d samples...\n', i, numSamples);
    end
    
    % Read audio file
    [audio, fs] = audioread(adsBalanced.Files{i});
    
    % Convert to mono if stereo
    if size(audio, 2) > 1
        audio = mean(audio, 2);
    end
    
    % Resample to intermediate rate if needed
    if fs ~= originalFs
        audio = resample(audio, originalFs, fs);
    end
    
    % Pad or truncate to 1 second
    if length(audio) < intermediateSamples
        audio = [audio; zeros(intermediateSamples - length(audio), 1)];
    elseif length(audio) > intermediateSamples
        audio = audio(1:intermediateSamples);
    end
    
    % Store raw audio for visualization
    audioRaw{i} = audio;
    
    % Compute frame-based RMS amplitude (envelope)
    frameSize = floor(intermediateSamples / nFrames);
    for f = 1:nFrames
        frameStart = (f-1) * frameSize + 1;
        frameEnd = min(f * frameSize, intermediateSamples);
        frameData = audio(frameStart:frameEnd);
        % RMS amplitude of frame
        audioFrames(i, f) = sqrt(mean(frameData.^2));
    end
    
    % NOTE: Per-sample normalization removed - using global normalization instead
    % This preserves relative differences between samples
end

fprintf('  Audio preprocessing complete!\n');

% Global normalization (preserves relative differences between samples)
globalMax = max(audioFrames(:));
globalMin = min(audioFrames(:));
fprintf('  RMS range before normalization: [%.4f, %.4f]\n', globalMin, globalMax);
if globalMax > 0
    audioFrames = audioFrames / globalMax;  % Scale to [0, 1] globally
end
fprintf('  Applied global normalization (preserves inter-sample differences)\n');

fprintf('  Frame-based data shape: %d samples x %d frames\n', size(audioFrames, 1), size(audioFrames, 2));

% Split into train/validation/test sets
fprintf('\n  Splitting dataset...\n');
trainRatio = 0.7;
valRatio = 0.15;
testRatio = 0.15;

nTrain = floor(numSamples * trainRatio);
nVal = floor(numSamples * valRatio);
nTest = numSamples - nTrain - nVal;

% Create random permutation for splitting
rng(42);
shuffleIdx = randperm(numSamples);

trainIdx = shuffleIdx(1:nTrain);
valIdx = shuffleIdx(nTrain+1:nTrain+nVal);
testIdx = shuffleIdx(nTrain+nVal+1:end);

audioTrain = audioFrames(trainIdx, :);
audioVal = audioFrames(valIdx, :);
audioTest = audioFrames(testIdx, :);

% Also keep raw audio for visualization
audioRawTrain = audioRaw(trainIdx);
audioRawVal = audioRaw(valIdx);
audioRawTest = audioRaw(testIdx);

labelsTrain = labels(trainIdx);
labelsVal = labels(valIdx);
labelsTest = labels(testIdx);

fprintf('  Train set: %d samples\n', nTrain);
fprintf('  Validation set: %d samples\n', nVal);
fprintf('  Test set: %d samples\n', nTest);

% Create one-hot encoded targets
fprintf('  Creating one-hot encoded targets...\n');
classNames = categories(targetClasses);

YTrain = zeros(nTrain, numClasses);
YVal = zeros(nVal, numClasses);
YTest = zeros(nTest, numClasses);

for i = 1:nTrain
    classIdx = find(strcmp(classNames, char(labelsTrain(i))));
    YTrain(i, classIdx) = 1;
end
for i = 1:nVal
    classIdx = find(strcmp(classNames, char(labelsVal(i))));
    YVal(i, classIdx) = 1;
end
for i = 1:nTest
    classIdx = find(strcmp(classNames, char(labelsTest(i))));
    YTest(i, classIdx) = 1;
end

% === 10.3 Build ESN for Classification ===
fprintf('\nStep 10.3: Building ESN for speech classification...\n');

% Reuse reservoir parameters from Sections 1-4
fprintf('  Reusing reservoir parameters from previous sections...\n');
fprintf('  Network: %d neurons (%d E, %d I)\n', n, n_E, n_I);
fprintf('  Adaptation: E(%d scales), I(%d scales)\n', n_a_E, n_a_I);
fprintf('  STD: E(%d), I(%d)\n', n_b_E, n_b_I);

% Use same dt as main reservoir (0.1s) - now matches frame-based audio
% dt_speech was already set in Section 10.2 as dt_speech = dt

% Create new input weight matrix with STRONGER scaling to drive reservoir
% This ensures the reservoir responds differently to different inputs
rng(123);
input_scaling_speech = 3.0;  % Strong input scaling to overcome autonomous dynamics
W_in_speech = (2 * rand(n, 1) - 1) * input_scaling_speech;
fprintf('  Input scaling for speech: %.1f (stronger to drive reservoir)\n', input_scaling_speech);
lambda_speech = 1e-6;
% Create ESN parameters struct reusing existing reservoir params
params_speech = struct();
params_speech.n = n;
params_speech.n_E = n_E;
params_speech.n_I = n_I;
params_speech.W = W;                             % Reuse recurrent weights
params_speech.W_in = W_in_speech;                % New input weights for 1D audio
params_speech.tau_d = tau_d;                     % Reuse dendritic time constant
params_speech.n_a_E = n_a_E;                     % Reuse adaptation config
params_speech.n_a_I = n_a_I;
params_speech.tau_a_E = tau_a_E;
params_speech.tau_a_I = tau_a_I;
params_speech.n_b_E = n_b_E;                     % Reuse STD config
params_speech.n_b_I = n_b_I;
params_speech.tau_b_E_rec = tau_b_E_rec;
params_speech.tau_b_E_rel = tau_b_E_rel;
params_speech.tau_b_I_rec = tau_b_I_rec;
params_speech.tau_b_I_rel = tau_b_I_rel;
params_speech.c_E = c_E;                         % Reuse adaptation scaling
params_speech.c_I = c_I;
params_speech.lags = lags;                       % Reuse synaptic delays
params_speech.activation_function = activation_function;  % Reuse activation
params_speech.which_states = 'x';
params_speech.include_input = false;
params_speech.lambda = lambda_speech;
params_speech.dt = dt_speech;                    % Time step for audio

% Create ESN instance for classification
esn_speech = SRNN_ESN(params_speech);

fprintf('  ESN created with existing reservoir configuration\n');
fprintf('  Chaos level: %.2f\n', level_of_chaos);
fprintf('  Input scaling: %.2f\n', input_scaling);
fprintf('  Audio dt: %.6f s (%.0f Hz)\n', dt_speech, 1/dt_speech);

% === 10.4 Extract Reservoir Features (Flattened) ===
fprintf('\nStep 10.4: Extracting reservoir features...\n');

% Function to process audio through reservoir and flatten all states
function [features, stateHistory] = processAudioThroughReservoir(esn, audioSample, nNeurons)
    % Reset reservoir state
    esn.resetState();
    
    % Reshape audio to column vector (n_timesteps x 1)
    U = audioSample(:);
    nTimesteps = length(U);
    
    % Run reservoir and get full state history
    [X, ~] = esn.runReservoir(U);
    
    % Store state history for visualization (T x N)
    stateHistory = X;
    
    % Flatten all states into single vector: (T x N) -> (T*N,)
    % This captures the full temporal dynamics as a "block"
    features = reshape(X', 1, []);  % Row vector of T*N features
end

% Feature dimension = T * N (flattened reservoir states)
nTimesteps = nFrames;  % Number of timesteps = number of audio frames
nFeatures = nTimesteps * n;  % T * N = 10 * 210 = 2100 features
fprintf('  Feature extraction: flattened reservoir states (T*N)\n');
fprintf('  Timesteps: %d, Neurons: %d\n', nTimesteps, n);
fprintf('  Feature dimension: %d (%d x %d)\n', nFeatures, nTimesteps, n);

% Extract features for training set
fprintf('  Processing training set (%d samples)...\n', nTrain);
XTrain = zeros(nTrain, nFeatures);
stateHistoryTrain = cell(nTrain, 1);  % Store for visualization

tic;
for i = 1:nTrain
    if mod(i, 100) == 0
        elapsed = toc;
        eta = elapsed / i * (nTrain - i);
        fprintf('    Sample %d/%d (ETA: %.1f sec)...\n', i, nTrain, eta);
    end
    [XTrain(i, :), stateHistoryTrain{i}] = processAudioThroughReservoir(esn_speech, audioTrain(i, :), n);
end
fprintf('  Training features extracted in %.1f seconds\n', toc);

% Extract features for validation set
fprintf('  Processing validation set (%d samples)...\n', nVal);
XVal = zeros(nVal, nFeatures);
stateHistoryVal = cell(nVal, 1);

tic;
for i = 1:nVal
    [XVal(i, :), stateHistoryVal{i}] = processAudioThroughReservoir(esn_speech, audioVal(i, :), n);
end
fprintf('  Validation features extracted in %.1f seconds\n', toc);

% Extract features for test set
fprintf('  Processing test set (%d samples)...\n', nTest);
XTest = zeros(nTest, nFeatures);
stateHistoryTest = cell(nTest, 1);  % Store for visualization

tic;
for i = 1:nTest
    [XTest(i, :), stateHistoryTest{i}] = processAudioThroughReservoir(esn_speech, audioTest(i, :), n);
end
fprintf('  Test features extracted in %.1f seconds\n', toc);

% === 10.5 Train Two-Layer Network with SGD ===
fprintf('\nStep 10.5: Training two-layer readout network with SGD...\n');

% Network architecture
nHidden = floor(n / 2);  % Hidden layer size = N/2 = 105
nOutput = numClasses;    % Output layer size = 6 classes

fprintf('  Architecture: %d -> %d (ReLU) -> %d (Softmax)\n', nFeatures, nHidden, nOutput);

% Activation functions
relu = @(x) max(0, x);
relu_deriv = @(x) double(x > 0);

softmax = @(z) exp(z - max(z, [], 2)) ./ sum(exp(z - max(z, [], 2)), 2);

% Cross-entropy loss
cross_entropy = @(y_pred, y_true) -sum(y_true .* log(y_pred + 1e-10), 2);

% Initialize weights (Xavier initialization)
rng(42);
W1 = randn(nFeatures, nHidden) * sqrt(2 / nFeatures);
b1 = zeros(1, nHidden);
W2 = randn(nHidden, nOutput) * sqrt(2 / nHidden);
b2 = zeros(1, nOutput);

% SGD hyperparameters with momentum
learningRate = 0.1;    % Increased learning rate
momentum = 0.9;        % Momentum coefficient
nEpochs = 100;         % More epochs for convergence
batchSize = 32;        % Mini-batch size

fprintf('  Learning rate: %.4f\n', learningRate);
fprintf('  Momentum: %.2f\n', momentum);
fprintf('  Epochs: %d\n', nEpochs);
fprintf('  Batch size: %d\n', batchSize);

% Normalize features (important for SGD stability)
XTrain_mean = mean(XTrain, 1);
XTrain_std = std(XTrain, 0, 1) + 1e-8;
XTrain_norm = (XTrain - XTrain_mean) ./ XTrain_std;
XVal_norm = (XVal - XTrain_mean) ./ XTrain_std;
XTest_norm = (XTest - XTrain_mean) ./ XTrain_std;

% Initialize velocity for momentum
vW1 = zeros(size(W1));
vb1 = zeros(size(b1));
vW2 = zeros(size(W2));
vb2 = zeros(size(b2));

% Training history
trainLossHistory = zeros(nEpochs, 1);
valLossHistory = zeros(nEpochs, 1);
trainAccHistory = zeros(nEpochs, 1);
valAccHistory = zeros(nEpochs, 1);

fprintf('  Training...\n');
tic;

for epoch = 1:nEpochs
    % Shuffle training data
    shuffleIdx = randperm(nTrain);
    epochLoss = 0;
    nBatches = ceil(nTrain / batchSize);
    
    for batch = 1:nBatches
        % Get batch indices
        batchStart = (batch - 1) * batchSize + 1;
        batchEnd = min(batch * batchSize, nTrain);
        batchIdx = shuffleIdx(batchStart:batchEnd);
        batchSizeActual = length(batchIdx);
        
        % Get batch data
        X_batch = XTrain_norm(batchIdx, :);
        Y_batch = YTrain(batchIdx, :);
        
        % Forward pass
        z1 = X_batch * W1 + b1;           % Pre-activation hidden
        h = relu(z1);                      % Hidden activations (ReLU)
        z2 = h * W2 + b2;                  % Pre-activation output
        y_pred = softmax(z2);              % Output probabilities (Softmax)
        
        % Compute loss
        batchLoss = mean(cross_entropy(y_pred, Y_batch));
        epochLoss = epochLoss + batchLoss * batchSizeActual;
        
        % Backward pass (gradient computation)
        % Output layer gradient (softmax + cross-entropy simplifies to y_pred - y_true)
        dz2 = (y_pred - Y_batch) / batchSizeActual;
        dW2 = h' * dz2;
        db2 = sum(dz2, 1);
        
        % Hidden layer gradient
        dh = dz2 * W2';
        dz1 = dh .* relu_deriv(z1);
        dW1 = X_batch' * dz1;
        db1 = sum(dz1, 1);
        
        % SGD with momentum update
        vW1 = momentum * vW1 - learningRate * dW1;
        vb1 = momentum * vb1 - learningRate * db1;
        vW2 = momentum * vW2 - learningRate * dW2;
        vb2 = momentum * vb2 - learningRate * db2;
        
        W1 = W1 + vW1;
        b1 = b1 + vb1;
        W2 = W2 + vW2;
        b2 = b2 + vb2;
    end
    
    % Compute epoch metrics
    trainLossHistory(epoch) = epochLoss / nTrain;
    
    % Training accuracy
    h_train = relu(XTrain_norm * W1 + b1);
    y_train_pred = softmax(h_train * W2 + b2);
    [~, predIdx] = max(y_train_pred, [], 2);
    [~, trueIdx] = max(YTrain, [], 2);
    trainAccHistory(epoch) = mean(predIdx == trueIdx) * 100;
    
    % Validation loss and accuracy
    h_val = relu(XVal_norm * W1 + b1);
    y_val_pred = softmax(h_val * W2 + b2);
    valLossHistory(epoch) = mean(cross_entropy(y_val_pred, YVal));
    [~, predIdx] = max(y_val_pred, [], 2);
    [~, trueIdx] = max(YVal, [], 2);
    valAccHistory(epoch) = mean(predIdx == trueIdx) * 100;
    
    % Print progress every 10 epochs
    if mod(epoch, 10) == 0 || epoch == 1
        fprintf('    Epoch %d/%d: Train Loss=%.4f, Train Acc=%.1f%%, Val Acc=%.1f%%\n', ...
                epoch, nEpochs, trainLossHistory(epoch), trainAccHistory(epoch), valAccHistory(epoch));
    end
end

fprintf('  Training completed in %.1f seconds\n', toc);
fprintf('  Final Train Accuracy: %.2f%%\n', trainAccHistory(end));
fprintf('  Final Val Accuracy: %.2f%%\n', valAccHistory(end));

% === 10.6 Evaluate on All Sets ===
fprintf('\nStep 10.6: Evaluating classification performance...\n');

% Function to compute predictions with two-layer network
function [predLabels, accuracy, probs] = evaluateTwoLayerNetwork(X_norm, Y, W1, b1, W2, b2, classNames)
    % Forward pass
    h = max(0, X_norm * W1 + b1);  % ReLU
    z = h * W2 + b2;
    probs = exp(z - max(z, [], 2)) ./ sum(exp(z - max(z, [], 2)), 2);  % Softmax
    
    % Get predicted class (argmax)
    [~, predIdx] = max(probs, [], 2);
    [~, trueIdx] = max(Y, [], 2);
    
    % Convert to categorical labels
    predLabels = categorical(classNames(predIdx));
    
    % Compute accuracy
    accuracy = mean(predIdx == trueIdx) * 100;
end

% Evaluate on all sets
[predTrain, accTrain, ~] = evaluateTwoLayerNetwork(XTrain_norm, YTrain, W1, b1, W2, b2, classNames);
fprintf('  Training Accuracy: %.2f%%\n', accTrain);

[predVal, accVal, ~] = evaluateTwoLayerNetwork(XVal_norm, YVal, W1, b1, W2, b2, classNames);
fprintf('  Validation Accuracy: %.2f%%\n', accVal);

[predTest, accTest, probsTest] = evaluateTwoLayerNetwork(XTest_norm, YTest, W1, b1, W2, b2, classNames);
fprintf('  Test Accuracy: %.2f%%\n', accTest);

% Get true labels for test set
[~, trueTestIdx] = max(YTest, [], 2);
trueTest = categorical(classNames(trueTestIdx));

% Store raw output for visualization (use probabilities)
rawOutputTest = probsTest;

% === 10.7 Visualization ===
fprintf('\nStep 10.7: Generating classification plots...\n');

% Create new figure for classification results
figure('Color', 'w', 'Name', 'Speech Classification Results');

% Subplot 1: Training Loss Curves
subplot(2, 4, 1);
plot(1:nEpochs, trainLossHistory, 'b-', 'LineWidth', 1.5);
hold on;
plot(1:nEpochs, valLossHistory, 'r-', 'LineWidth', 1.5);
hold off;
xlabel('Epoch');
ylabel('Cross-Entropy Loss');
title('Training Progress: Loss');
legend('Train', 'Validation', 'Location', 'best');
grid on;

% Subplot 2: Training Accuracy Curves
subplot(2, 4, 2);
plot(1:nEpochs, trainAccHistory, 'b-', 'LineWidth', 1.5);
hold on;
plot(1:nEpochs, valAccHistory, 'r-', 'LineWidth', 1.5);
hold off;
xlabel('Epoch');
ylabel('Accuracy (%)');
title('Training Progress: Accuracy');
legend('Train', 'Validation', 'Location', 'best');
ylim([0 100]);
grid on;

% Subplot 3: Confusion Matrix
subplot(2, 4, 3);
cm = confusionmat(trueTest, predTest);
confusionchart(cm, classNames, ...
    'Title', 'Test Confusion Matrix', ...
    'RowSummary', 'row-normalized', ...
    'ColumnSummary', 'column-normalized');

% Subplot 4: Final Accuracy Bar Chart
subplot(2, 4, 4);
accuracies = [accTrain, accVal, accTest];
bar(1:3, accuracies, 'FaceColor', [0.2 0.6 0.8]);
set(gca, 'XTickLabel', {'Train', 'Val', 'Test'});
ylabel('Accuracy (%)');
title('Final Accuracy');
ylim([0 100]);
grid on;

% Add accuracy values on bars
for i = 1:3
    text(i, accuracies(i) + 2, sprintf('%.1f%%', accuracies(i)), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 8);
end

% Subplot 5: Per-class Accuracy
subplot(2, 4, 5);
perClassAcc = zeros(numClasses, 1);
for i = 1:numClasses
    classIdxTrue = trueTestIdx == i;
    [~, predIdx] = max(rawOutputTest, [], 2);
    perClassAcc(i) = mean(predIdx(classIdxTrue) == i) * 100;
end
bar(1:numClasses, perClassAcc, 'FaceColor', [0.8 0.4 0.2]);
set(gca, 'XTickLabel', classNames);
ylabel('Accuracy (%)');
title('Per-Class Test Accuracy');
ylim([0 100]);
grid on;
xtickangle(45);

% Subplot 6: Frame-based Audio Features
subplot(2, 4, 6);
nSamplesToShow = 6;
rng(123);
sampleIdx = randperm(nTest, min(nSamplesToShow, nTest));
t_frames = (0:nFrames-1) * dt_speech;

hold on;
colors = lines(nSamplesToShow);
for i = 1:nSamplesToShow
    idx = sampleIdx(i);
    offset = (i-1) * 1.2;
    plot(t_frames, audioTest(idx, :) + offset, '-o', 'Color', colors(i,:), 'LineWidth', 1.5, 'MarkerSize', 4);
    
    % Add label
    trueLabel = char(trueTest(idx));
    predLabel = char(predTest(idx));
    if strcmp(trueLabel, predLabel)
        labelStr = sprintf('%s (OK)', trueLabel);
        textColor = [0 0.5 0];
    else
        labelStr = sprintf('%s->%s', trueLabel, predLabel);
        textColor = [0.8 0 0];
    end
    text(t_frames(end) + 0.02, offset + 0.5, labelStr, ...
         'Color', textColor, 'FontSize', 7, 'FontWeight', 'bold');
end
hold off;
xlabel('Time (s)');
ylabel('RMS Amplitude');
title('Sample Predictions');
xlim([0, max(t_frames) + 0.25]);
grid on;

% Subplot 7: Feature Distribution (Flattened States)
subplot(2, 4, 7);
histogram(XTest_norm(:), 50, 'Normalization', 'probability', ...
          'FaceColor', [0.4 0.4 0.8], 'EdgeColor', 'none');
xlabel('Normalized Feature Value');
ylabel('Probability');
title('Flattened Feature Dist.');
grid on;

% Subplot 8: Network Weights Visualization
subplot(2, 4, 8);
% Show W2 weights (hidden to output) as heatmap
imagesc(W2');
colormap(gca, bluewhitered_colormap(256));
colorbar;
xlabel('Hidden Unit');
ylabel('Class');
set(gca, 'YTick', 1:numClasses, 'YTickLabel', classNames);
title('Output Layer Weights (W2)');

% Overall title
sgtitle(sprintf('Fractional ESN Speech Classification (Two-Layer SGD) | Test Accuracy: %.1f%%', accTest), ...
        'FontSize', 12, 'FontWeight', 'bold');

fprintf('  Classification figure complete!\n');

% === 10.7b Neuron Dynamics Visualization ===
fprintf('  Generating neuron dynamics visualization...\n');

% Create separate figure for neuron dynamics (one example per word)
figure('Color', 'w', 'Name', 'Reservoir Dynamics per Word');

% Find one example of each class from test set
exampleIdx = zeros(numClasses, 1);
for i = 1:numClasses
    classIndices = find(trueTestIdx == i);
    if ~isempty(classIndices)
        exampleIdx(i) = classIndices(1);  % Take first example of each class
    end
end

% Time vector for reservoir states
t_reservoir = (0:nFrames-1) * dt_speech;

% Plot neuron states for each word class
for i = 1:numClasses
    subplot(2, 3, i);
    
    idx = exampleIdx(i);
    stateHistory = stateHistoryTest{idx};  % (nFrames x n)
    
    % Plot as heatmap
    imagesc(t_reservoir, 1:n, stateHistory');
    colormap(gca, 'jet');
    colorbar;
    
    xlabel('Time (s)');
    ylabel('Neuron Index');
    
    % Get prediction for this example
    predLabel = char(predTest(idx));
    if strcmp(classNames{i}, predLabel)
        titleColor = [0 0.5 0];
        titleStr = sprintf('"%s" (correct)', classNames{i});
    else
        titleColor = [0.8 0 0];
        titleStr = sprintf('"%s" -> "%s"', classNames{i}, predLabel);
    end
    title(titleStr, 'Color', titleColor, 'FontWeight', 'bold');
    
    % Add E/I boundary line
    hold on;
    yline(n_E + 0.5, 'w--', 'LineWidth', 1.5);
    hold off;
end

sgtitle('Reservoir Neuron Dynamics per Word Class (E|I boundary in white)', ...
        'FontSize', 14, 'FontWeight', 'bold');

fprintf('  Neuron dynamics figure complete!\n\n');

% === 10.8 Classification Summary ===
fprintf('=================================================================\n');
fprintf('  SPEECH CLASSIFICATION SUMMARY\n');
fprintf('=================================================================\n');
fprintf('Dataset:\n');
fprintf('  - Classes: %s\n', strjoin(classNames', ', '));
fprintf('  - Samples per class: %d\n', maxSamplesPerClass);
fprintf('  - Total samples: %d (Train: %d, Val: %d, Test: %d)\n', ...
        numSamples, nTrain, nVal, nTest);
fprintf('\nReservoir Configuration:\n');
fprintf('  - Neurons: %d (%d E, %d I)\n', n, n_E, n_I);
fprintf('  - Chaos level: %.2f\n', level_of_chaos);
fprintf('  - Adaptation: E(%d scales), I(%d scales)\n', n_a_E, n_a_I);
fprintf('  - STD: E(%d), I(%d)\n', n_b_E, n_b_I);
fprintf('  - Speech input scaling: %.2f\n', input_scaling_speech);
fprintf('\nTwo-Layer Readout Network:\n');
fprintf('  - Architecture: %d -> %d (ReLU) -> %d (Softmax)\n', nFeatures, nHidden, nOutput);
fprintf('  - Training: SGD with momentum + cross-entropy loss\n');
fprintf('  - Learning rate: %.4f, Momentum: %.2f\n', learningRate, momentum);
fprintf('  - Epochs: %d, Batch size: %d\n', nEpochs, batchSize);
fprintf('\nResults:\n');
fprintf('  - Training Accuracy: %.2f%%\n', accTrain);
fprintf('  - Validation Accuracy: %.2f%%\n', accVal);
fprintf('  - Test Accuracy: %.2f%%\n', accTest);
fprintf('\nPer-class Accuracy:\n');
for i = 1:numClasses
    fprintf('  - %s: %.1f%%\n', classNames{i}, perClassAcc(i));
end
fprintf('=================================================================\n\n');

fprintf('Section 10 completed: Speech Recognition Classification!\n');

