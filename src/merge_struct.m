function s = merge_struct(base, over)
% MERGE_STRUCT Copy OVER's fields onto BASE, with OVER winning.
%
% Usage:
%   s = merge_struct(srnn_param_preset('overconnected'), cfg.model)
%
% Used to layer a parameter preset (the base) under the run_mode timing
% settings (which win). No collision checking is done here on purpose:
% ParamSpaceAnalysis2.validate_model_defaults reports anything unusable at
% run(), which is the single place that knows the grid and the conditions.
%
% See also: srnn_param_preset, analysis_run_config, ParamSpaceAnalysis2

if ~isstruct(base) || ~isscalar(base)
    error('merge_struct:InvalidInput', 'base must be a scalar struct.');
end
if ~isstruct(over) || ~isscalar(over)
    error('merge_struct:InvalidInput', 'over must be a scalar struct.');
end

s = base;
f = fieldnames(over);
for i = 1:numel(f)
    s.(f{i}) = over.(f{i});
end
end
