function [first, second] = mc_gain_sweep_preset_point(preset_name, model_args, second_axis, first_axis)
% MC_GAIN_SWEEP_PRESET_POINT A network's own values on the two axes of a sweep: where it sits on the sweep's map.
%
%   [gain, bias]   = MC_GAIN_SWEEP_PRESET_POINT(preset_name)
%   [gain, bias]   = MC_GAIN_SWEEP_PRESET_POINT(preset_name, model_args)
%   [gain, second] = MC_GAIN_SWEEP_PRESET_POINT(preset_name, model_args, second_axis)
%   [first, second] = MC_GAIN_SWEEP_PRESET_POINT(preset_name, model_args, second_axis, first_axis)
%
%   preset_name  a SRNNCellTypePairs preset
%   model_args   name-value overrides of the preset (cell), as
%                run_mc_gain_sweep takes them; default none
%   second_axis  the model property of the sweep's second axis; default
%                'S_c', the bias
%   first_axis   the model property of its first axis; default
%                'level_of_chaos', the gain
%   first, second  that network's value of each property: for S_c, the
%                center of its setpoint draw. A property the model holds as
%                one value per cell type (self_feedback) is returned as the
%                one value the types share; a network whose types differ
%                has no single point on the axis, which is an error.
%
% Read off a model constructed (not built) from the preset, the overrides
% and the preset's first condition, so the class's defaults and its
% last-write-wins rule apply and are not restated. The sweep records the
% point in its settings, the reader of a run in progress recovers it from a
% job file's signature, and the bundle's stage runs each variant network at
% its value of the second axis.
%
% See also: run_mc_gain_sweep, read_mc_gain_sweep_jobs, mc_gain_sweep_stage_plan,
%           mc_sweep_axis

arguments
    preset_name (1,:) char
    model_args  (1,:) cell = {}
    second_axis (1,:) char = 'S_c'
    first_axis  (1,:) char = 'level_of_chaos'
end

[preset, model_class, conditions] = srnn_param_preset(preset_name);
if ~strcmp(model_class, 'SRNNCellTypePairs')
    error('mc_gain_sweep_preset_point:BadModelClass', ...
        'Preset ''%s'' is written for %s; the gain sweep of memory capacity runs SRNNCellTypePairs presets.', ...
        preset_name, model_class);
end
condition_args = struct2namevalue(rmfield(conditions{1}, 'name'));
preset_args = struct2namevalue(preset);
model = SRNNCellTypePairs(preset_args{:}, model_args{:}, condition_args{:});
first = one_value(model, first_axis, preset_name);
second = one_value(model, second_axis, preset_name);
end

%% ------------------------------------------------------------------------
function v = one_value(model, property, preset_name)
v = unique(model.(property));
if ~(isnumeric(v) && isscalar(v))
    error('mc_gain_sweep_preset_point:NotOneValue', ...
        'Preset ''%s'' has no single value of %s (it holds %s), so it has no point on that axis.', ...
        preset_name, property, mat2str(model.(property)));
end
end
