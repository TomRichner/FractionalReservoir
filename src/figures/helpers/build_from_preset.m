function [model, info] = build_from_preset(preset_name, condition_name, varargin)
% BUILD_FROM_PRESET Construct and build a model from a preset + one condition.
%
%   model = BUILD_FROM_PRESET('celltype_pairs_Sc0p2_noise0p025_dualStd_4cond')
%   model = BUILD_FROM_PRESET(p, 'sfa_only')
%   model = BUILD_FROM_PRESET(p, 'sfa_and_std', 'n', 50, 'T_range', [0 20])
%
% Trailing name-value pairs are applied LAST, so they win over both the preset
% and the condition -- the same precedence run_single_job uses for a swept
% parameter.
%
% Returns the BUILT model, plus an info struct recording model_class, the
% condition used, and the integrator that was selected.
%
% WHY THIS EXISTS -- three traps, each of which bit a figure during the refactor:
%
% 1. THE INTEGRATOR. A preset carrying sigma_u_noise > 0 REQUIRES a stochastic
%    scheme; construction fails outright with
%      "sigma_u_noise = 0.025 requires a stochastic integrator"
%    against the class default 'ode45'. In the sweep pipeline
%    analysis_run_config picks one, but a figure does not go through it, so each
%    figure was re-deriving the rule. Selected here once: sra1 when the preset is
%    stochastic ('sra1' costs the same two drift evaluations as heun but is
%    strong order 1.5), rk4 otherwise.
%
% 2. ARGUMENT EXPANSION. The preset and the condition are separate structs, and
%
%      feval(cls, struct2namevalue(preset), struct2namevalue(cond))
%
%    passes TWO CELLS as two arguments -- the constructor sees a cell where it
%    expects a property name and fails deep inside its parser. They must be
%    concatenated into one flat cell and expanded with {:}.
%
% 3. THE CONDITION CARRIES THE ADAPTATION COUNTS. n_a / n_a_E / n_b_E and
%    synapse_config are condition-owned and never appear in a preset, so a model
%    built from the preset ALONE has no adaptation at all. Naming a condition is
%    not optional decoration.
%
% See also: srnn_param_preset, srnn_adaptation_conditions, struct2namevalue

arguments
    preset_name    (1,:) char
    condition_name (1,:) char = 'sfa_and_std'
end
arguments (Repeating)
    varargin
end

[preset, model_class, conditions] = srnn_param_preset(preset_name);

names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
ci = find(strcmp(names, condition_name), 1);
if isempty(ci)
    error('build_from_preset:NoSuchCondition', ...
        'Preset ''%s'' has no condition ''%s''. Available: %s.', ...
        preset_name, condition_name, strjoin(names, ', '));
end
cond = rmfield(conditions{ci}, 'name');

% Integrator, from the preset's noise. See trap 1 above.
if isfield(preset, 'sigma_u_noise') && any(preset.sigma_u_noise(:) > 0)
    solver = 'sra1';
else
    solver = 'rk4';
end

args = [struct2namevalue(preset), struct2namevalue(cond), ...
        {'ode_solver', solver}, varargin];

model = feval(model_class, args{:});
model.build();

info = struct('model_class', model_class, 'preset_name', preset_name, ...
              'condition', condition_name, 'ode_solver', solver);
end
