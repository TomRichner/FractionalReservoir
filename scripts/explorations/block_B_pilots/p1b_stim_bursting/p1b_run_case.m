function C = p1b_run_case(condition_name, stim, S)
% P1B_RUN_CASE Run one adaptation condition under one stimulus case; return its traces and measures.
%
%   C = P1B_RUN_CASE(condition_name, stim, S)
%
% condition_name  'original' | 'chain_mts' | 'product_mts', built by
%                 bursting_condition_model (the one definition)
% stim            one stimulus case (see p1b_stim_bursting): name, label,
%                 input_config
% S               the settings struct of p1b_stim_bursting
%
% C holds, over [0, T_show] at 100 Hz: t, mean_r_E, mean_r_I, the per-neuron
% rates r_E and r_I (single), mean_adaptation (mean (c/K) sum_k a_k over E),
% mean_depression (mean depression factor over every depressed route, from
% SRNNCellTypePairs.std_factors, i.e. std_route_factor), u_mean (input
% averaged over all neurons), u_target (input averaged over the stimulated
% neurons), the bursts (p1b_detect_bursts over the whole [0, T_show]),
% measures (p1b_measures) and values (every value the model was built with).

[model, values] = bursting_condition_model(condition_name, ...
    'rng_seeds', S.rng_seeds, 'T_settle', S.T_settle, 'T_show', S.T_show, ...
    'input_config', stim.input_config, 'plot_deci', S.plot_deci);
model.run();

p      = model.plot_data;
params = model.get_params();
keep   = p.t >= 0 & p.t <= S.T_show;
C = struct();
C.condition = condition_name;
C.stim_name = stim.name;
C.stim_label = stim.label;
C.t   = reshape(p.t(keep), 1, []);
rE = p.r.E(:, keep); rI = p.r.I(:, keep);
C.r_E = single(rE);
C.r_I = single(rI);
C.mean_r_E = mean(rE, 1);
C.mean_r_I = mean(rI, 1);

K = size(p.a.E, 2);
C.mean_adaptation = squeeze(mean(model.c(1) / K * sum(p.a.E(:, :, keep), 2), 1))';

factor = SRNNCellTypePairs.std_factors(p.b, params);
names = model.cell_type_names;
acc = []; routes = {};
for i = 1:numel(names)
    for j = 1:numel(names)
        v = factor.(names{i}).(names{j});
        if isempty(v), continue; end
        acc = [acc; v(:, keep)]; %#ok<AGROW>
        routes{end+1} = sprintf('%s->%s', names{i}, names{j}); %#ok<AGROW>
    end
end
C.mean_depression   = mean(acc, 1);
C.depression_routes = routes;

u = [p.u.E; p.u.I];
u = u(:, keep);
C.u_mean = mean(u, 1);
target = stim.target;
if isempty(target), target = 1:model.n; end
C.target = target;
C.u_target = mean(u(target, :), 1);

C.bursts   = p1b_detect_bursts(C.t, C.mean_r_E, S.burst);
C.measures = p1b_measures(C, S);
values = rmfield(values, 'input_config');   % function handle; described in stim
values.stim = rmfield(stim, 'input_config');
values.stim.input_config_fields = rmfield(stim.input_config, 'generator');
values.stim.generator = func2str(stim.input_config.generator);
C.values = values;
end
