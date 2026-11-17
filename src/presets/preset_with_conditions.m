function [preset, model_class, conditions] = preset_with_conditions(preset_name, requested)
% PRESET_WITH_CONDITIONS A preset's three outputs, its conditions narrowed to the ones a run simulates (FD-030).
%
%   [preset, model_class, conditions] = PRESET_WITH_CONDITIONS(preset_name)
%   [preset, model_class, conditions] = PRESET_WITH_CONDITIONS(preset_name, requested)
%
%   preset_name  srnn_param_preset key.
%   requested    cell of condition names to simulate, or empty ({}, '', [])
%                for every condition the preset states. Default empty.
%
%   preset, model_class   srnn_param_preset's first two outputs, unchanged.
%   conditions            its third output, the requested conditions only, in
%                         the requested order; with an empty request, all of
%                         them in the preset's order, exactly as
%                         srnn_param_preset returns them.
%
% WHY. A bundle may simulate fewer conditions than its preset states (Tom,
% 2026-11-16, FD-030: product MTS is left out of the mu7chain runs from
% 'fast' upward). The preset is not edited and no second preset is added:
% the bundle names the conditions it runs, and every stage takes that list
% as its 'conditions' argument and resolves it here, in place of calling
% srnn_param_preset directly. A stage given no list behaves exactly as
% before.
%
% THE RULES are select_conditions' (one statement for figures and runs): a
% name the preset does not state is an error, raised here as
% 'preset_with_conditions:NoSuchCondition' with the preset named; a
% duplicate is an error; names match exactly.
%
% See also: srnn_param_preset, select_conditions, run_all_paper_analyses

arguments
    preset_name (1,:) char
    requested = {}
end

[preset, model_class, conditions] = srnn_param_preset(preset_name);
available = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
try
    [~, idx] = select_conditions(available, requested);
catch ME
    if strcmp(ME.identifier, 'select_conditions:NoSuchCondition')
        missing = setdiff(cellstr(requested), available, 'stable');
        error('preset_with_conditions:NoSuchCondition', ...
            'Preset ''%s'' states no condition named %s. Its conditions: %s.', ...
            preset_name, strjoin(missing, ', '), strjoin(available, ', '));
    end
    rethrow(ME);
end
conditions = conditions(idx);
end
