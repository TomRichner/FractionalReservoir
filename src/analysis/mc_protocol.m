function cfg = mc_protocol(cfg)
% MC_PROTOCOL The memory-capacity experiment, identical in every run mode.
%
%   cfg = MC_PROTOCOL()      the protocol as a struct
%   cfg = MC_PROTOCOL(cfg)   the same fields added to an existing struct
%
% THE ONE PLACE these values are written. run_memory_capacity (the paired
% trials at the preset network) and run_mc_gain_sweep (memory capacity at
% every point of a gain-by-bias sweep) both take them from here, so the two
% stages measure the same task. mc_delay_axis gives the delays it scores.
%
% Durations, in seconds and Hz. None of it is worth trading for speed:
%   T_train_sec  decides whether MC means anything at all. The ridge readout
%                is fit on T_train_sec / T_hold hold-samples and needs more
%                of them than it has features (the network's n). 600 s gives
%                2000, well posed for n = 300 and n = 500.
%   d_max_sec    caps the horizon that can be measured; observed horizons
%                reach 5.5 to 5.7 s, so 15 s leaves room.
%
% See also: mc_delay_axis, run_memory_capacity, run_mc_gain_sweep, SRNN_ESN_reservoir

if nargin < 1
    cfg = struct();
end

cfg.T_wash_sec  = 10;     % washout before anything is recorded
cfg.T_train_sec = 600;    % readout training window -> 2000 hold-samples
cfg.T_test_sec  = 150;    % readout test window
cfg.d_max_sec   = 15;     % longest delay scored
cfg.fs          = 200;

% sample_hold: i.i.d. values held for T_hold, so MC is measured in HOLD UNITS.
% That is the fair choice for a low-pass reservoir and is free of the
% autocorrelation inflation that white input produces.
cfg.input_type     = 'sample_hold';
cfg.T_hold         = 0.3;      % s; sets the MC delay increment
cfg.u_f_cutoff     = 5;        % only used by 'bandlimited'
cfg.u_alpha        = 1;        % only used by 'one_over_f'
% 'synaptic' reads out the synaptic output (the rate times the route's
% depression and facilitation factor), exposing the STD state to the linear
% readout. Unchanged for the no-STD conditions, where it equals the rate.
cfg.readout_signal = 'synaptic';
cfg.R2_threshold_for_horizon = 0.10;
cfg.seed_net_base  = 3000;
cfg.seed_stim_base = 4000;
end
