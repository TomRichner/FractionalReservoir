function [m, groups, info] = run13_condition_model(condition_name, opts)
%RUN13_CONDITION_MODEL The run13 IED-like network under one adaptation condition, built and not run.
%
%   [m, groups, info] = run13_condition_model(condition_name)
%   [m, groups, info] = run13_condition_model(condition_name, Name, Value, ...)
%
% THE ONE PLACE the three conditions of the run13 exploration are defined
% (pilot 0a, 2026-11-12). Follow-up pilots call this and do not re-derive
% the conditions.
%
% condition_name   'original'     run13 as ied_build_network builds it: one
%                                 SFA timescale per neuron, one depression
%                                 pool on E outputs (E->E, E->I), none on I
%                  'chain_mts'    three SFA timescales (same total c),
%                                 depression on the same routes as the
%                                 chained vesicle pool (std.model = 'chain')
%                  'product_mts'  three SFA timescales, depression on the
%                                 same routes as the product of two
%                                 independent pools
%
% Name-value options
%   duration      s observed after the 10-s settling interval (default 100);
%                 T_range = [-10 duration], T_plot = [0 duration]
%   seed          network seed (default 42, run13's). W, the setpoints, the
%                 SFA ladders and x0 follow it
%   noise_seed    Wiener seed (default 224779 + (seed - 42), the offset rule
%                 of ied_stim_experiment; 224779 is run13's)
%   input_config  [] (default) = run13's: zero input. Otherwise a struct
%                 passed to ied_build_network as is, which copies its fields
%                 into a fresh input_config; give a complete one, for example
%                 a dc_staircase_stimulus struct as ied_stim_experiment does
%   lya_method    default 'none' (the trajectory does not depend on it)
%   model_args    cell of further name-value model arguments, applied LAST,
%                 so they override everything above, the condition's tau_a
%                 and synapse_config included
%
% Outputs
%   m       the IEDExplorationNetwork, built, run13's modular W in place;
%           call m.run()
%   groups  n x 1 module label (1..5) of each neuron
%   info    struct: condition, connectivity (ied_connectivity's summary),
%           cfg (ied_base_config with duration, seed and noise seed set),
%           and params, every value this function chose (see below)
%
% THE CONSTANTS (Tom, 2026-11-12: run13's value where run13 has the
% parameter; the mu7chain bundle's current chained-MTS default where only
% the MTS form has it):
%   SFA ladder   original E 0.40 s, I 0.75 s (run13, one rung, log-normal
%                per-neuron spread [0.70 0.50]); MTS E [0.40 1.58 10] s,
%                I [0.75 1.58 10] s: run13's rung is the fastest, 1.58 and
%                10 s are the bundle's log_ladder(0.25, 10, 3) slow rungs;
%                same spread, applied to the ladder's ends per neuron
%   c            [0.5 0.5], the total SFA budget, in every condition
%   original     tau_rec 2 s, tau_rel 0.25 s (run13)
%   chain        tau_rec [2 5] s, tau_rel 0.25 s, n2 2 (releasable pool
%                run13's; reserve refill 5 s and n2 2 the bundle's)
%   product      tau_rec [2 4] s, tau_rel [0.25 0.5] s (first pool run13's,
%                second the bundle's)
% Depression on E outputs only; no depression spread (run13 has none).
%
% See also: example_run13_chain_mts, ied_build_network, ied_base_config

arguments
    condition_name (1,:) char {mustBeMember(condition_name, {'original', 'chain_mts', 'product_mts'})}
    opts.duration (1,1) double {mustBePositive} = 100
    opts.seed (1,1) double = 42
    opts.noise_seed double = []
    opts.input_config = []
    opts.lya_method (1,:) char = 'none'
    opts.model_args cell = {}
end

cfg = ied_base_config();
seed0 = cfg.seed;
cfg.duration = opts.duration;
cfg.seed = opts.seed;
if isempty(opts.noise_seed)
    cfg.noise_seed = cfg.noise_seed + (opts.seed - seed0);
else
    cfg.noise_seed = opts.noise_seed;
end
assert(~cfg.std_on_I, 'run13 depresses E outputs only.');

% The constants of the MTS forms (bundle defaults, mu7chain preset).
ladder     = log_ladder(0.25, 10, 3);
sfa_slow   = ladder(2:3);
chain_rec2 = 5;   chain_n2  = 2;
prod2_rec  = 4;   prod2_rel = 0.5;
assert(all(sfa_slow(1) > cfg.tau_a) && chain_rec2 > cfg.std_rec ...
    && ~isequal([prod2_rec prod2_rel], [cfg.std_rec cfg.std_rel]), ...
    'A bundle default does not make sense beside run13''s own value.');

std1 = struct('tau_rec', cfg.std_rec, 'tau_rel', cfg.std_rel);
switch condition_name
    case 'original'
        tau_a = {cfg.tau_a(1), cfg.tau_a(2)};
        std_cfg = std1;
        cond_args = {};      % exactly ied_build_network's construction
    case 'chain_mts'
        tau_a = {[cfg.tau_a(1) sfa_slow], [cfg.tau_a(2) sfa_slow]};
        std_cfg = struct('model', 'chain', 'tau_rec', [cfg.std_rec chain_rec2], ...
                         'tau_rel', cfg.std_rel, 'n2', chain_n2);
        cond_args = {'tau_a', tau_a, 'synapse_config', e_routes(std_cfg)};
    case 'product_mts'
        tau_a = {[cfg.tau_a(1) sfa_slow], [cfg.tau_a(2) sfa_slow]};
        std_cfg = struct('tau_rec', [cfg.std_rec prod2_rec], 'tau_rel', [cfg.std_rel prod2_rel]);
        cond_args = {'tau_a', tau_a, 'synapse_config', e_routes(std_cfg)};
end

% run13's input: run_ied_experiment's, i.e. its base preset's input_config
% with the drive zeroed and every step off. (ied_build_network's own
% default, a struct with only those two fields, no longer builds: the
% class's step generator now requires n_steps.)
if isempty(opts.input_config)
    preset_args = srnn_param_preset('celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0p025_dualStd_3cond_mu7revisedMedium');
    input_config = preset_args.input_config;
    input_config.intrinsic_drive = 0;
    input_config.no_stim_pattern = true(1, 3);
else
    input_config = opts.input_config;
end

[m, groups, connectivity] = ied_build_network(cfg, input_config, ...
    'lya_method', opts.lya_method, cond_args{:}, opts.model_args{:});

params = struct( ...
    'condition', condition_name, ...
    'n', cfg.n, 'f', cfg.f, 'groups', cfg.groups, ...
    'p_within', cfg.p_within, 'p_between', cfg.p_between, ...
    'wE', cfg.wE, 'wI', cfg.wI, 'w_sd', cfg.w_sd, 'between_scale', cfg.between_scale, ...
    'S_c_mean', cfg.Sc, 'S_c_sd', cfg.Sc_sd, ...
    'c', cfg.c, 'tau_a', {tau_a}, 'tau_a_spread', cfg.tau_spread, ...
    'tau_a_run13', cfg.tau_a, 'std_routes', 'E->E, E->I', 'std', std_cfg, ...
    'std_run13', std1, ...
    'sigma_u_noise', cfg.noise, 'fs', cfg.fs, 'ode_solver', 'sra1', ...
    'warmup', cfg.warmup, 'duration', cfg.duration, ...
    'seed', cfg.seed, 'rng_seeds', [cfg.seed cfg.seed + 1], 'noise_seed', cfg.noise_seed, ...
    'input_config', input_config, 'model_args', {opts.model_args}, ...
    'tau_a_clipped', m.tau_a_clipped);
info = struct('condition', condition_name, 'cfg', cfg, ...
    'connectivity', connectivity, 'params', params);
end

function sc = e_routes(std_cfg)
% Depression on E outputs only, as run13.
sc = struct();
sc.E.E.std = std_cfg;  sc.E.I.std = std_cfg;
sc.I.E = struct();     sc.I.I = struct();
end
