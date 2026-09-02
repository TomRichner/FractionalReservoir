function [u_ex, t_ex] = generate_paired_pulse_input(params, T, fs, rng_seed, pulse_config)
% GENERATE_PAIRED_PULSE_INPUT Generates paired-pulse input for SRNN simulation.
%
% Syntax:
%   [u_ex, t_ex] = generate_paired_pulse_input(params, T, fs, rng_seed, pulse_config)
%
% Inputs:
%   params       - Struct with .n (number of neurons)
%   T            - Duration (seconds)
%   fs           - Sampling frequency (Hz)
%   rng_seed     - Random seed
%   pulse_config - Struct with:
%                  .amp          - Pulse amplitude (default 3)
%                  .width_ms     - Pulse width in ms (default 0.5)
%                  .neuron_frac  - Fraction of neurons to stimulate (default 0.2)
%
% Description:
%   Generates a repeating 10s cycle of stimulation:
%   - t=1.25s: Single pulse
%   - t=3.75s: Paired pulse, ISI = 2 ms
%   - t=6.25s: Paired pulse, ISI = 10 ms
%   - t=8.75s: Paired pulse, ISI = 100 ms
%
%   The same random set of neurons is stimulated for all pulses.

    if nargin < 5
        pulse_config = struct();
    end
    
    % Default configuration
    amp = 3;
    if isfield(pulse_config, 'amp'), amp = pulse_config.amp; end
    
    width_ms = 0.5;
    if isfield(pulse_config, 'width_ms'), width_ms = pulse_config.width_ms; end
    
    neuron_frac = 0.2;
    if isfield(pulse_config, 'neuron_frac'), neuron_frac = pulse_config.neuron_frac; end

    rng(rng_seed);

    dt = 1/fs;
    t_ex = (0:dt:T)';
    nt = length(t_ex);
    u_ex = zeros(params.n, nt);

    % Select neurons to stimulate
    n_stim = round(params.n * neuron_frac);
    stim_indices = randperm(params.n, n_stim);

    % Cycle parameters
    cycle_len_sec = 10;
    cycle_len_samples = round(cycle_len_sec * fs);
    
    % Pulse width in samples
    width_samples = max(1, round(width_ms / 1000 * fs));
    
    % Base times for events within a cycle (seconds)
    event_times = [1.25, 3.75, 6.25, 8.75];
    
    % ISIs for the events (seconds). 0 means single pulse.
    isis = [0, 0.002, 0.010, 0.100];

    % Create one cycle template
    cycle_template = zeros(1, cycle_len_samples);
    
    for i = 1:length(event_times)
        t_start = event_times(i);
        isi = isis(i);
        
        % First pulse
        idx1 = round(t_start * fs) + 1;
        if idx1 + width_samples - 1 <= cycle_len_samples
            cycle_template(idx1 : idx1 + width_samples - 1) = amp;
        end
        
        % Second pulse (if ISI > 0)
        if isi > 0
            t_start_2 = t_start + isi;
            idx2 = round(t_start_2 * fs) + 1;
            if idx2 + width_samples - 1 <= cycle_len_samples
                cycle_template(idx2 : idx2 + width_samples - 1) = amp;
            end
        end
    end

    % Replicate cycle
    num_cycles = ceil(nt / cycle_len_samples);
    full_stim_trace = repmat(cycle_template, 1, num_cycles);
    full_stim_trace = full_stim_trace(1:nt); % Trim to actual length

    % Apply to selected neurons
    u_ex(stim_indices, :) = repmat(full_stim_trace, n_stim, 1);

end

