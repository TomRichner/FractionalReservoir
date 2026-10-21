function [u_ex, t_ex, info] = generate_am_pulse_train(params, T, fs, rng_seed, config)
% GENERATE_AM_PULSE_TRAIN Pulse train whose amplitude is sinusoidally modulated.
%
%   [u_ex, t_ex]       = GENERATE_AM_PULSE_TRAIN(params, T, fs, rng_seed, config)
%   [u_ex, t_ex, info] = GENERATE_AM_PULSE_TRAIN(...)
%
% The human sine-stimulation protocol (brain-stim-analysis, CLAUDE.md
% "Stimuli"): rectangular pulses every pulse_period (5 Hz) whose amplitude
% follows
%
%   A(t) = amp * (1 - (mod_depth/2) * (1 + cos(2 pi (t - t_first) / mod_period)))
%
% so the FIRST pulse sits at the envelope's trough, (1 - mod_depth) of the
% peak, and the peak comes half a modulation period later. With the default
% mod_depth 0.3 that is the experiment's 0.85 - 0.15 cos(...) times the
% peak: 70 % to 100 %, relative modulation 0.15/0.85 = 0.176.
%
% Usable as input_config.generator on SRNNCellTypePairs (the model calls it
% with two outputs; the third is for the analysis that has to know where
% the pulses are). Stimulus time 0 is t = 0 here; a negative T_range(1)
% prepends zeros for settling.
%
% Inputs:
%   params   - model params: .n (and what stim_target_indices needs)
%   T        - duration (s); t_ex = (0:1/fs:T)'
%   fs       - sampling rate (Hz)
%   rng_seed - seed of the driven subset (stim_target_indices)
%   config   - struct (the model's input_config), optional fields:
%              .amp           peak pulse amplitude, input units (default 1)
%              .width_s       pulse width in SECONDS (default 1e-3); rounded
%                             to whole samples, at least one
%              .pulse_period  s between pulses (default 0.2: 5 Hz)
%              .mod_period    s per modulation cycle (default 2)
%              .mod_depth     peak-to-trough depth relative to the peak
%                             (default 0.3)
%              .t_first       time of the first pulse, s (default 0)
%              .t_stop        no pulse at or after this time (default Inf)
%              .neuron_frac, .target_types   the driven subset, see
%                             stim_target_indices
%
% Outputs:
%   u_ex - n x numel(t_ex) input
%   t_ex - (0:1/fs:T)'
%   info - struct: pulse_times (s), pulse_amps, pulse_phase (radians of
%          the envelope, 0 = trough, pi = peak), width_samples, stim_idx,
%          and the resolved settings
%
% Note on the pulse shape: the model reads u through a LINEAR interpolant,
% so a pulse of w samples is a trapezoid whose integral is exactly
% w * (1/fs) * amp -- the charge is right even though the edges are ramps.
%
% Revived from generate_AM_pulse_train (deleted in 6effc85), which put the
% envelope PEAK on the first pulse and took the width in samples.
%
% See also: generate_paired_pulse_input, stim_target_indices, run_stim_evoked

amp          = get_field(config, 'amp', 1);
width_s      = get_field(config, 'width_s', 1e-3);
pulse_period = get_field(config, 'pulse_period', 0.2);
mod_period   = get_field(config, 'mod_period', 2);
mod_depth    = get_field(config, 'mod_depth', 0.3);
t_first      = get_field(config, 't_first', 0);
t_stop       = get_field(config, 't_stop', Inf);

dt   = 1 / fs;
t_ex = (0:dt:T)';
nt   = numel(t_ex);
u_ex = zeros(params.n, nt);

stim_idx      = stim_target_indices(params, rng_seed, config);
width_samples = max(1, round(width_s * fs));

pulse_times = t_first:pulse_period:T;
pulse_times = pulse_times(pulse_times < t_stop);
phase       = 2 * pi * (pulse_times - t_first) / mod_period;
pulse_amps  = amp * (1 - (mod_depth / 2) * (1 + cos(phase)));

keep = false(size(pulse_times));
for k = 1:numel(pulse_times)
    i0 = round(pulse_times(k) * fs) + 1;
    i1 = i0 + width_samples - 1;
    if i1 <= nt
        u_ex(stim_idx, i0:i1) = pulse_amps(k);
        keep(k) = true;
    end
end

info = struct('pulse_times', pulse_times(keep), 'pulse_amps', pulse_amps(keep), ...
    'pulse_phase', mod(phase(keep), 2 * pi), 'width_samples', width_samples, ...
    'stim_idx', stim_idx, 'amp', amp, 'width_s', width_s, ...
    'pulse_period', pulse_period, 'mod_period', mod_period, ...
    'mod_depth', mod_depth, 't_first', t_first);
end

function v = get_field(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end
