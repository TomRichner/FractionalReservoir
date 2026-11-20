function [m, info] = small_network_model(condition_name, opts)
%SMALL_NETWORK_MODEL The 80-neuron version of the mu7chain network under one condition, built and not run.
%
%   m = small_network_model(condition_name)
%   [m, info] = small_network_model(condition_name, Name, Value, ...)
%
% THE ONE PLACE the 80-neuron network is defined (lead's decision,
% 2026-11-19). Later scripts call this and do not re-derive it.
%
% The network is the mu7chain bundle's network
% (celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain)
% with ONLY its size changed, by the rule of the repository's 100-neuron
% regression preset (regression_pairs_noStim_noise0p025_dualStd_3cond):
%
%   n                 500 -> 80
%   indegree          100 -> 16    (connection probability alpha stays 0.2;
%                                   the mask is Bernoulli)
%   F_ref_n           500 -> 80    (F_tracks_network stays false), so the
%   F_ref_indegree    100 -> 16    weight unit is the small network's own:
%                                  F = 1/sqrt(80 * 0.2 * 1.8) = 0.18634
%
% Everything else is the preset's: mu_tilde_relative [7 -7; 7 -7],
% sigma_tilde_relative 1.5, f = [0.5 0.5], per-neuron setpoints, SFA and
% depression spreads, noise 0.025, no stimulus. Consequences: the
% theoretical bulk radius R = 4.799 is the same as at 500 neurons (with F
% the network's own, R depends only on alpha and the tildes), while the
% mean (rank-one, E minus I) part of W has norm about N alpha F 7 = 20.9
% against 52.2 at 500 neurons.
%
% condition_name   one of the preset's conditions:
%                  'no_adaptation', 'sfa1_std1' (one timescale),
%                  'sfa3_stdc2' (chained multiple-timescale), 'sfa3_std2'
%                  (product multiple-timescale; not used by
%                  explore_small_network)
%
% Name-value options
%   seed         network seed pair rng_seeds (default [7 8], the mu7chain
%                examples' seed): W, the setpoints, the spreads and x0
%                follow the first entry, the input the second
%   T_range      simulated interval in s (default [0 15])
%   fs           sampling rate in Hz (default 500)
%   lya_method   default 'none'; 'topk' with the lya_* options below
%   lya_K, lya_dt, lya_T_interval, lya_warmup
%                Lyapunov settings, used only when lya_method is not
%                'none' (defaults 20, 0.05 s, [5 15] s, 5 s); lya_K_auto
%                is always false here, so exactly lya_K exponents are taken
%   model_args   cell of further name-value model arguments, applied LAST,
%                so they override everything above, the size included
%
% Outputs
%   m      the SRNNCellTypePairs model, built (W in place), not run;
%          call m.run()
%   info   struct: preset, condition, and params (every value this
%          function set on top of the preset)
%
% See also: explore_small_network, build_from_preset, srnn_param_preset

arguments
    condition_name (1,:) char {mustBeMember(condition_name, ...
        {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2', 'sfa3_std2'})}
    opts.seed (1,2) double = [7 8]
    opts.T_range (1,2) double = [0 15]
    opts.fs (1,1) double {mustBePositive} = 500
    opts.lya_method (1,:) char = 'none'
    opts.lya_K (1,1) double {mustBeInteger, mustBePositive} = 20
    opts.lya_dt (1,1) double {mustBePositive} = 0.05
    opts.lya_T_interval (1,2) double = [5 15]
    opts.lya_warmup (1,1) double {mustBeNonnegative} = 5
    opts.model_args cell = {}
end

preset_name = ['celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_' ...
    'stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain'];

% The size override: the only change to the preset's network.
size_args = {'n', 80, 'indegree', 16, 'F_ref_n', 80, 'F_ref_indegree', 16};

% The run settings. The preset is stochastic (noise 0.025), so the
% integrator is the stochastic sra1, as build_from_preset would pick.
run_args = {'rng_seeds', opts.seed, 'T_range', opts.T_range, ...
    'fs', opts.fs, 'ode_solver', 'sra1', 'lya_method', opts.lya_method};
if ~strcmpi(opts.lya_method, 'none')
    run_args = [run_args, {'lya_K', opts.lya_K, 'lya_K_auto', false, ...
        'lya_dt', opts.lya_dt, 'lya_T_interval', opts.lya_T_interval, ...
        'lya_warmup', opts.lya_warmup}];
end

m = build_from_preset(preset_name, condition_name, ...
    size_args{:}, run_args{:}, opts.model_args{:});

params = struct( ...
    'n', 80, 'indegree', 16, 'F_ref_n', 80, 'F_ref_indegree', 16, ...
    'rng_seeds', opts.seed, 'T_range', opts.T_range, 'fs', opts.fs, ...
    'ode_solver', 'sra1', 'lya_method', opts.lya_method, ...
    'lya_K', opts.lya_K, 'lya_K_auto', false, 'lya_dt', opts.lya_dt, ...
    'lya_T_interval', opts.lya_T_interval, 'lya_warmup', opts.lya_warmup, ...
    'model_args', {opts.model_args});
info = struct('preset', preset_name, 'condition', condition_name, ...
    'params', params);
end
