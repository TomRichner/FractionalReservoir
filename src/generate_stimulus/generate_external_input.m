function [u_ex, t_ex] = generate_external_input(params, T, fs, rng_seed, input_config)
% generate_external_input - Generate external input for SRNN simulation
%
% Syntax:
%   [u_ex, t_ex] = generate_external_input(params, T, fs, rng_seed, input_config)
%
% Description:
%   Creates a random sparse step function as external input to the network.
%   The input consists of steps with random amplitudes applied to a sparse
%   subset of neurons. Optionally adds constant intrinsic drive.
%
%   If input_config contains a 'generator' field with a function handle,
%   that function is called instead to generate custom stimulus patterns.
%
% Inputs:
%   params       - Struct containing:
%                  .n - Total number of neurons
%   T            - Duration of simulation (seconds)
%   fs           - Sampling frequency (Hz)
%   rng_seed     - Random seed for reproducibility
%   input_config - Struct containing either:
%                  (Standard mode):
%                  .n_steps          - Number of steps in the input
%                  .step_density     - Fraction of neurons receiving input (0-1)
%                                      [used if step_density_E/I not specified]
%                  .step_density_E   - (Optional) Fraction of E neurons receiving input
%                  .step_density_I   - (Optional) Fraction of I neurons receiving input
%                  .positive_only    - (Optional) If true, generate only positive amplitudes
%                  .amp              - Amplitude scaling factor
%                  .no_stim_pattern  - Logical array (1 x n_steps) where true = no stim
%                  .intrinsic_drive  - Constant drive added to all neurons (n x 1)
%
%                  (Custom mode):
%                  .generator        - Function handle: @(params, T, fs, rng_seed, cfg)
%                                      Must return [u_ex, t_ex]
%
% Outputs:
%   u_ex - n x nt external input matrix
%   t_ex - nt x 1 time vector
%
% Example (standard mode):
%   params.n = 100;
%   input_config.n_steps = 7;
%   input_config.step_density = 0.2;
%   input_config.amp = 1;
%   input_config.no_stim_pattern = false(1, 7);
%   input_config.no_stim_pattern(1:2:end) = true;
%   input_config.intrinsic_drive = zeros(100, 1);
%   [u_ex, t_ex] = generate_external_input(params, 150, 100, 2, input_config);
%
% Example (custom generator):
%   input_config.generator = @my_custom_stimulus_generator;
%   [u_ex, t_ex] = generate_external_input(params, 150, 100, 2, input_config);

% Check for custom generator function
if isfield(input_config, 'generator') && isa(input_config.generator, 'function_handle')
    [u_ex, t_ex] = input_config.generator(params, T, fs, rng_seed, input_config);
    return;
end

% Standard sparse step stimulus generation
rng(rng_seed);  % Set random seed for reproducibility

% Create time vector
dt = 1 / fs;
t_ex = (0:dt:T)';
nt = length(t_ex);

% Extract configuration parameters
n_steps = input_config.n_steps;
step_density = input_config.step_density;
amp = input_config.amp;
no_stim_pattern = input_config.no_stim_pattern;
intrinsic_drive = input_config.intrinsic_drive;

% Calculate step period and length
step_period = fix(T / n_steps);
step_length = round(step_period * fs);  % Number of time points per step

% Precompute random sparse step matrix (vectorized)
% Generate all random amplitudes at once: n x n_steps
if isfield(input_config, 'positive_only') && input_config.positive_only
    random_sparse_step = amp * abs(randn(params.n, n_steps));
else
    random_sparse_step = amp * randn(params.n, n_steps);
end

% Create sparse mask with separate E/I densities if specified
sparse_mask = false(params.n, n_steps);

% Handle E neurons
if isfield(input_config, 'step_density_E')
    density_E = input_config.step_density_E;
else
    density_E = step_density;  % Fall back to uniform step_density
end

% Handle I neurons
if isfield(input_config, 'step_density_I')
    density_I = input_config.step_density_I;
else
    density_I = step_density;  % Fall back to uniform step_density
end

% Get E/I indices (use params if available, otherwise compute from f)
if isfield(params, 'E_indices') && isfield(params, 'I_indices')
    E_idx = params.E_indices;
    I_idx = params.I_indices;
else
    % Fallback: compute from f if provided, otherwise assume 50% E
    f_val = 0.5;
    if isfield(params, 'f')
        f_val = params.f;
    end
    n_E = round(params.n * f_val);
    E_idx = 1:n_E;
    I_idx = (n_E + 1):params.n;
end

% Apply masks separately to E and I populations
sparse_mask(E_idx, :) = rand(length(E_idx), n_steps) < density_E;
sparse_mask(I_idx, :) = rand(length(I_idx), n_steps) < density_I;

random_sparse_step = random_sparse_step .* sparse_mask;

% Apply no-stimulation pattern (zero out specified steps)
random_sparse_step(:, no_stim_pattern) = 0;

% Create u_ex using the precomputed random sparse steps
u_ex = zeros(params.n, nt);
for step_idx = 1:n_steps
    % Determine time indices for this step
    start_idx = (step_idx - 1) * step_length + 1;
    end_idx = min(step_idx * step_length, nt);

    if start_idx > nt
        break;
    end

    % Apply the step to all neurons (broadcasting the column vector)
    u_ex(:, start_idx:end_idx) = repmat(random_sparse_step(:, step_idx), 1, end_idx - start_idx + 1);
end

% Add intrinsic drive to neurons
u_ex = u_ex + intrinsic_drive;
end


