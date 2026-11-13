function [gain, bias] = mc_gain_sweep_preset_point(preset_name, model_args)
% MC_GAIN_SWEEP_PRESET_POINT A network's own gain and bias: where it sits on the gain-by-bias map.
%
%   [gain, bias] = MC_GAIN_SWEEP_PRESET_POINT(preset_name)
%   [gain, bias] = MC_GAIN_SWEEP_PRESET_POINT(preset_name, model_args)
%
%   preset_name  a SRNNCellTypePairs preset
%   model_args   name-value overrides of the preset (cell), as
%                run_mc_gain_sweep takes them; default none
%   gain         level_of_chaos of that network
%   bias         S_c, the centre of its setpoint draw
%
% Read off a model constructed (not built) from the preset, the overrides
% and the preset's first condition, so the class's defaults and its
% last-write-wins rule apply and are not restated. The sweep records the
% point in its settings, the reader of a run in progress recovers it from a
% job file's signature, and the bundle's stage runs each gain-only network
% at this bias.
%
% See also: run_mc_gain_sweep, read_mc_gain_sweep_jobs, mc_gain_sweep_stage_plan

arguments
    preset_name (1,:) char
    model_args  (1,:) cell = {}
end

[preset, model_class, conditions] = srnn_param_preset(preset_name);
if ~strcmp(model_class, 'SRNNCellTypePairs')
    error('mc_gain_sweep_preset_point:BadModelClass', ...
        'Preset ''%s'' is written for %s; the gain-by-bias sweep runs SRNNCellTypePairs presets.', ...
        preset_name, model_class);
end
condition_args = struct2namevalue(rmfield(conditions{1}, 'name'));
preset_args = struct2namevalue(preset);
model = SRNNCellTypePairs(preset_args{:}, model_args{:}, condition_args{:});
gain = model.level_of_chaos;
bias = model.S_c;
end
