function [model, params, ws, net] = cx_network(plan, entry, T)
%CX_NETWORK One network of the concept_examples stage, built (not run): the model, its params and weight statistics.
%
%   [model, params, ws, net] = cx_network(plan, entry)        structure only (T_range [0 1])
%   [model, params, ws, net] = cx_network(plan, entry, T)     ready to run over [0 T]
%
%   plan    concept_examples_stage_plan(cfg)
%   entry   one element of plan.examples or plan.gallery
%   model   the built OspecNet (an SRNNCellTypePairs) of the preset's
%           'no_adaptation' condition, with entry.level_of_chaos, the
%           network name's block means and added structure (ospec_network),
%           and rng_seeds entry.seed
%   params  model.get_params()
%   ws      ospec_weight_stats(params)
%   net     ospec_network's struct for the name
%
% The build is ospec_build_model's (the global stream set to the twister
% first, so a seed names one network in a client and on a pool worker).
% With T given, the Lyapunov settings follow plan.lya: Benettin over
% [plan.lya_from, T], warmup min(5, lya_from).
%
% See also: run_concept_examples_stage, ospec_build_model, ospec_network

arguments
    plan  (1,1) struct
    entry (1,1) struct
    T double = []
end
p = srnn_param_preset(plan.preset_name);
net = ospec_network(entry.name, p.mu_tilde_relative);
size_args = {};
if plan.n_override > 0 && plan.n_override ~= p.n
    % The outlier stages' size rule (ostruct_structural_check): alpha 0.2
    % and the smaller network's own weight unit.
    k_in = max(2, round(0.2 * plan.n_override));
    size_args = {'n', plan.n_override, 'indegree', k_in, 'F_ref_n', plan.n_override, 'F_ref_indegree', k_in};
end
if isempty(T)
    T_range = [0 1];
    lya_args = {'lya_method', 'none'};
else
    T_range = [0 T];
    if plan.lya
        lya_args = {'lya_method', 'benettin', 'lya_T_interval', [plan.lya_from T], ...
            'lya_warmup', min(5, plan.lya_from)};
    else
        lya_args = {'lya_method', 'none'};
    end
end
B = struct('preset_name', plan.preset_name, 'condition', plan.condition, ...
    'net_args', {[net.net_args, {'level_of_chaos', entry.level_of_chaos}]}, 'size_args', {size_args}, ...
    'added_structure', net.added_structure, 'rng_seeds', entry.seed, 'T_range', T_range, 'fs', plan.fs, ...
    'lya_args', {lya_args}, 'ode_solver', plan.ode_solver);
model = ospec_build_model(B);
params = model.get_params();
ws = ospec_weight_stats(params);
end
