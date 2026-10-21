function cfg = sine_paired_stim_config()
% SINE_PAIRED_STIM_CONFIG The sine- and paired-pulse stimulation experiment on the paper's network (TR, 2026-10-20).
%
%   cfg = sine_paired_stim_config();
%   run_stim_evoked(cfg);  fig_sine_paired_stim('run_dir', cfg.run_dir, 'out_dir', ...)
%
% Or press Run on sine_paired_stim_run.m, which does both.
%
% WHAT IT IS FOR. Tom's Minnesota Neuromodulation poster: in human SEEG the
% cortico-cortical evoked potential to a 5 Hz pulse train whose amplitude is
% sinusoidally modulated (2, 5, 10 s periods) LEADS the stimulus envelope by
% about 10 degrees at every period (medians 9.0, 6.6, 3.9 deg), and the
% second response of a 100 ms pair is suppressed far below the linear
% prediction. This bundle drives the manuscript's multiple-timescale
% network with those very stimuli, against single-timescale and no
% adaptation, to show how SFA and STD produce both.
%
% Every setting is stated here; no other config is called. Every value of
% the network is the manuscript preset's, unchanged (the mu7revisedMedium
% network: n = 500, piecewise activation, S_c 0.2 +/- 0.1, SFA on E and I
% c = 0.5 with the [0.25 1.58 10] s ladder and tau_a_spread 0.25, dual STD on
% every route, sigma_u_noise 0.025, so the integrator is sra1).
%
% DEVIATIONS FROM THE MANUSCRIPT RUNS, all in the stimulus and the timing:
%   fs 2000 Hz (the paper runs at 400 Hz). A pulse must be shorter than the
%       2 ms pair interval and is 1 ms = 2 samples here; at 400 Hz a pulse is
%       one 2.5-ms sample and the 2 ms pair cannot exist.
%   input_config: the stimulus generators below instead of the preset's
%       (empty) step pattern.
%   Only E neurons are driven (stim.target_types {'E'}, half of them). In
%       this network E and I are statistically identical except for the sign
%       of their output, so driving both alike makes the evoked recurrent
%       input to every other neuron cancel to first order (a pilot on
%       2026-10-20: RMS 10-60 ms of the undriven neurons' mean x 0.006
%       driving E and I, 0.061 driving E only, at the same amplitude).
%   The CCEP-like readout is the mean x of the UNDRIVEN neurons (the
%       recording site reached only through the network); the mean over all
%       neurons is also reported, but it is dominated by the driven
%       neurons' direct, linear jump in x, which no adaptation acts on.
%   amp 40 with 1 ms pulses: a 0.4 jump in a driven neuron's x at the peak
%       of the envelope (amp * width / tau_d), against the setpoint 0.2,
%       x's own noise std about 0.06 and the linear range of the piecewise
%       sigmoid (+/- 0.4 about S_c). Pilot: amp 100 already saturates the
%       undriven response (RMS 0.074 against 0.061 at 40).
%
% Durations: 10 s of noise alone before the first pulse (T_range(1) = -10),
% then 20 / 10 / 6 cycles at 2 / 5 / 10 s; the first discard_s = 10 s of
% stimulation (the slow adaptation building up to the train's mean) is left
% out of every average, leaving 15 / 8 / 5 cycles. Paired pulse: 5 cycles of
% the 10 s protocol, the first left out.
%
% Lyapunov estimation off; no parallel pool (other work shares the machine).

cfg = struct();

%% The network
cfg.preset_name = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0p025_dualStd_3cond_mu7revisedMedium';
cfg.conditions  = {};          % {} -> every condition of the preset
cfg.rng_seeds   = [1 2];       % W / setpoints / ladders; driven subset (seed 2). Same for every condition
cfg.sigma_u_noise = [];        % [] -> the preset's 0.025

%% Timing
cfg.fs           = 2000;       % Hz; see the deviations above
cfg.readout_deci = 2;          % saved readouts at 1000 Hz
cfg.settle_s     = 10;         % noise alone before the first pulse
cfg.discard_s    = 10;         % stimulation time left out of the averages

%% The protocols
cfg.sine_periods  = [2 5 10];  % s
cfg.sine_cycles   = [20 10 6]; % cycles simulated per period
cfg.paired_cycles = 5;         % 10-s paired-pulse cycles

cfg.stim = struct( ...
    'amp',          40, ...    % peak pulse amplitude (input units; u_ex_scale stays 1)
    'width_s',      1e-3, ...  % 1 ms pulses
    'pulse_period', 0.2, ...   % 5 Hz
    'mod_depth',    0.3, ...   % 70 % .. 100 % of the peak, first pulse at the trough
    'neuron_frac',  0.5, ...   % half of the target pool
    'target_types', {{'E'}});  % E neurons only; see the deviations above

%% Where things land
cfg.run_dir  = 'data/sine_paired_stim';   % must be empty; run_stim_evoked refuses otherwise
cfg.fig_root = 'figs/sine_paired_stim';
cfg.verbose  = 'minimal';
end
