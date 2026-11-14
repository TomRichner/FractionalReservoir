function opts = mc_gain_sweep_options(opts)
% MC_GAIN_SWEEP_OPTIONS The options of one gain sweep of memory capacity, defaulted and checked.
%
%   opts = MC_GAIN_SWEEP_OPTIONS(name, value, ...)
%
% run_mc_gain_sweep's options, in the one place that names them and gives
% their defaults. run_mc_gain_sweep passes its arguments through here, and
% so does run_mc_gain_sweep_set for each run of a set, so a single run and
% a run of a set cannot differ in what an option means. What each option
% does is documented in run_mc_gain_sweep's header.
%
% See also: run_mc_gain_sweep, run_mc_gain_sweep_set

arguments
    opts.preset_name  (1,:) char   = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0_chainStd_4cond_mu7chain'
    opts.run_mode     (1,:) char   = 'essential'
    opts.gain_levels  (1,:) double = []
    opts.second_axis  (1,:) char   = 'S_c'
    opts.second_levels (1,:) double = []
    opts.bias_levels  (1,:) double = []     % second_levels, for the default axis
    opts.seeds        (1,:) double = []
    opts.conditions   (1,:) cell   = {}
    opts.output_dir   (1,:) char   = ''
    opts.ode_solver   (1,:) char   = 'sra1'
    opts.verbose                   = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (or a logical); see verbose_level
    opts.model_args   (1,:) cell   = {}
    opts.protocol     (1,1) struct = struct()
    opts.lya_K        (1,1) double {mustBeInteger, mustBePositive} = 5
    opts.lya_T_sec    (1,1) double {mustBePositive} = 30
    opts.lya_acc_sec  (1,1) double {mustBePositive} = 20
    opts.segment_len  (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(opts.segment_len, 2)} = 600
    opts.retry_failed (1,1) logical = true
    opts.use_parallel (1,1) logical = true
end
end
