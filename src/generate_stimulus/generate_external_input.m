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
    random_sparse_step = amp * randn(params.n, n_steps);
    
    % Make it sparse based on step_density
    sparse_mask = rand(params.n, n_steps) < step_density;
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


