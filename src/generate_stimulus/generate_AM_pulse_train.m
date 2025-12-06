function [u_ex, t_ex] = generate_AM_pulse_train(params, T, fs, rng_seed, pulse_config)
% GENERATE_AM_PULSE_TRAIN Generates amplitude-modulated pulse train input.
%
% Syntax:
%   [u_ex, t_ex] = generate_AM_pulse_train(params, T, fs, rng_seed, pulse_config)
%
% Inputs:
%   params       - Struct with .n (number of neurons)
%   T            - Duration (seconds)
%   fs           - Sampling frequency (Hz)
%   rng_seed     - Random seed
%   pulse_config - Struct with:
%                  .amp_max      - Maximum pulse amplitude (default 3)
%                  .width_ms     - Pulse width in ms (default 0.5)
%                  .neuron_frac  - Fraction of neurons to stimulate (default 0.2)
%                  .pulse_period - Time between pulses (default 0.2s)
%                  .mod_freq     - Modulation frequency (default 0.5 Hz)
%                  .mod_depth    - Modulation depth (default 0.5)
%                                  Amplitude oscillates between amp_max and (1-mod_depth)*amp_max
%
% Description:
%   Generates a pulse train with pulses every 0.2s.
%   The amplitude is modulated sinusoidally at 0.5 Hz.
%   The peak of the modulation aligns with a pulse time.
%   The same random set of neurons is stimulated for all pulses.

    if nargin < 5
        pulse_config = struct();
    end
    
    % Default configuration
    amp_max = 3;
    if isfield(pulse_config, 'amp_max'), amp_max = pulse_config.amp_max; end
    
    width_ms = 0.5;
    if isfield(pulse_config, 'width_ms'), width_ms = pulse_config.width_ms; end
    
    neuron_frac = 0.2;
    if isfield(pulse_config, 'neuron_frac'), neuron_frac = pulse_config.neuron_frac; end
    
    pulse_period = 0.2;
    if isfield(pulse_config, 'pulse_period'), pulse_period = pulse_config.pulse_period; end
    
    mod_freq = 0.5;
    if isfield(pulse_config, 'mod_freq'), mod_freq = pulse_config.mod_freq; end
    
    mod_depth = 0.5;
    if isfield(pulse_config, 'mod_depth'), mod_depth = pulse_config.mod_depth; end

    rng(rng_seed);

    dt = 1/fs;
    t_ex = (0:dt:T)';
    nt = length(t_ex);
    u_ex = zeros(params.n, nt);

    % Select neurons to stimulate
    n_stim = round(params.n * neuron_frac);
    stim_indices = randperm(params.n, n_stim);

    % Pulse width in samples
    width_samples = max(1, round(width_ms / 1000 * fs));
    
    % Generate pulses
    % Pulse times: 0, 0.2, 0.4, ...
    pulse_times = 0:pulse_period:T;
    
    for i = 1:length(pulse_times)
        t_pulse = pulse_times(i);
        
        % Calculate amplitude modulation
        % We want peak at t=0 (aligned with first pulse)
        % Formula: A(t) = A_max * (1 - depth/2 * (1 - cos(2*pi*f*t)))
        % When cos=1 (t=0), A = A_max * (1 - depth/2 * 0) = A_max
        % When cos=-1, A = A_max * (1 - depth/2 * 2) = A_max * (1 - depth)
        % This oscillates between A_max and A_max*(1-depth).
        
        mod_val = cos(2 * pi * mod_freq * t_pulse);
        current_amp = amp_max * (1 - (mod_depth/2) * (1 - mod_val));
        
        % Apply pulse
        idx_start = round(t_pulse * fs) + 1;
        idx_end = idx_start + width_samples - 1;
        
        if idx_end <= nt
            u_ex(stim_indices, idx_start:idx_end) = current_amp;
        end
    end

end

