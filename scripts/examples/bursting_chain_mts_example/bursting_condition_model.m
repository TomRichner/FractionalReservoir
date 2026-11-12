function [model, values] = bursting_condition_model(condition_name, opts)
% BURSTING_CONDITION_MODEL The bursting network under one of its three named adaptation conditions.
%
%   [model, values] = BURSTING_CONDITION_MODEL(condition_name)
%   [model, values] = BURSTING_CONDITION_MODEL(condition_name, name, value, ...)
%
% THE ONE PLACE the three conditions of the block-B bursting pilots are
% defined. A follow-up pilot (stimulation, ablation, DC versus pulses) calls
% this and does not re-derive them.
%
% condition_name   'original'     the preset's own fully adapted condition
%                                 (bursting_pairs: sfa3_std1)
%                  'chain_mts'    the same, with depression as the chained
%                                 vesicle pool on the same routes
%                  'product_mts'  the same, with depression as the product of
%                                 two independent pools on the same routes
%
% Name-value options
%   'preset_name'  'bursting_pairs'   the network (see srnn_param_preset)
%   'rng_seeds'    [42 42]            W, setpoints, stimulus, x0, noise path
%   'T_settle'     10                 s simulated before t = 0
%   'T_show'       100                s after t = 0; T_range = [-T_settle, T_show]
%   'input_config' []                 [] = no external input at all (one
%                                     step that is never on, intrinsic_drive
%                                     0); otherwise handed to the model as is
%   'plot_deci'    4                  400 Hz integration -> 100 Hz plot_data
%   'model_args'   {}                 further name-value pairs for the
%                                     constructor, applied LAST, so they win
%                                     over the preset AND the condition (a
%                                     stimulus, another synapse_config,
%                                     another tau_a, sigma_u_noise, ...)
%   'chain_std'    the chain's depression struct (default below)
%   'prod_std'     the product's depression struct (default below)
%   'build'        true               false returns the model unbuilt
%
% Returns the BUILT (not run) model, and `values`, a struct of every value
% the condition was built with: the preset's fields, the condition's tau_a
% and synapse_config, the integrator, seeds, T_range, input_config and
% model_args.
%
% THE CONSTANTS (Tom, 2026-11-12: parameters the original has keep the
% original's values; parameters only the MTS forms have take the current
% defaults of the mu7chain bundle's chained MTS condition):
%   original     SFA tau_a{E} = [0.25 1.581 10] s (log_ladder(0.25, 10, 3)),
%                c = [0.5 0] (total budget, divided by K = 3; no SFA on I);
%                STD on E->E and E->I: tau_rec 1 s, tau_rel 0.25 s.
%   chain_mts    SFA as the original (it already has three rungs).
%                E->E, E->I: model 'chain', tau_rec [1 5] s (releasable
%                refill = original's tau_rec; reserve refill 5 s = mu7chain
%                default), tau_rel 0.25 s (original's), n2 = 2 (mu7chain
%                default).
%   product_mts  SFA as the original. E->E, E->I: tau_rec [1 4] s,
%                tau_rel [0.25 0.5] s (first pool = original's; second
%                pool = the mu7chain bundle's product second pool).
% Every other value (n 50, indegree 10, f [0.7 0.3], mu_tilde_relative
% [3 -2; 3 -2], sigma 1, level_of_chaos 1, piecewise S_a 0.9 S_c 0.5,
% tau_d 0.1, no noise) is the preset's, untouched.
%
% See also: example_bursting_chain_mts, bursting_chain_mts_figure,
%           build_from_preset, srnn_param_preset

arguments
    condition_name (1,:) char {mustBeMember(condition_name, {'original', 'chain_mts', 'product_mts'})}
    opts.preset_name  (1,:) char   = 'bursting_pairs'
    opts.rng_seeds    (1,2) double = [42 42]
    opts.T_settle     (1,1) double = 10
    opts.T_show       (1,1) double = 100
    opts.input_config              = []
    opts.plot_deci    (1,1) double = 4
    opts.model_args   (1,:) cell   = {}
    opts.chain_std    (1,1) struct = struct('model', 'chain', 'tau_rec', [1 5], 'tau_rel', 0.25, 'n2', 2)
    opts.prod_std     (1,1) struct = struct('tau_rec', [1 4], 'tau_rel', [0.25 0.5])
    opts.build        (1,1) logical = true
end

input_config = opts.input_config;
if isempty(input_config)
    input_config = struct('n_steps', 1, 'step_density', struct(), 'amp', 0, ...
        'no_stim_pattern', true, 'intrinsic_drive', 0, 'positive_only', false);
end

[preset, model_class, conds] = srnn_param_preset(opts.preset_name);
orig_name = full_adaptation_condition(conds);
orig = conds{strcmp(cellfun(@(c) c.name, conds, 'UniformOutput', false), orig_name)};

% The MTS conditions replace the depression on exactly the routes the
% original depresses, and nowhere else.
sc = orig.synapse_config;
switch condition_name
    case 'original'
        new_std = [];
    case 'chain_mts'
        new_std = opts.chain_std;
    case 'product_mts'
        new_std = opts.prod_std;
end
if ~isempty(new_std)
    pres = fieldnames(sc);
    for i = 1:numel(pres)
        posts = fieldnames(sc.(pres{i}));
        for j = 1:numel(posts)
            if isfield(sc.(pres{i}).(posts{j}), 'std')
                sc.(pres{i}).(posts{j}).std = new_std;
            end
        end
    end
end

% The integrator rule is build_from_preset's (sra1 if the preset is noisy,
% rk4 otherwise); a model_args 'ode_solver' or 'sigma_u_noise' wins.
if isfield(preset, 'sigma_u_noise') && any(preset.sigma_u_noise(:) > 0)
    solver = 'sra1';
else
    solver = 'rk4';
end
T_range = [-opts.T_settle, opts.T_show];
cond_args = {'tau_a', orig.tau_a, 'synapse_config', sc};
run_args  = {'ode_solver', solver, 'T_range', T_range, 'rng_seeds', opts.rng_seeds, ...
             'input_config', input_config, 'lya_method', 'none', ...
             'plot_deci', opts.plot_deci, 'verbose', 'near-none'};
args = [struct2namevalue(preset), cond_args, run_args, opts.model_args];

model = feval(model_class, args{:});
if opts.build
    model.build();
end

values = struct();
values.condition_name     = condition_name;
values.preset_name        = opts.preset_name;
values.preset_condition   = orig_name;
values.model_class        = model_class;
values.preset             = preset;
values.tau_a              = orig.tau_a;
values.synapse_config     = sc;
values.ode_solver         = model.ode_solver;
values.fs                 = model.fs;
values.rng_seeds          = opts.rng_seeds;
values.T_range            = T_range;
values.input_config       = input_config;
values.model_args         = opts.model_args;
values.sigma_u_noise      = model.sigma_u_noise;
values.c                  = model.c;
end
