function A = mc_sweep_second_axis(source)
% MC_SWEEP_SECOND_AXIS The second axis of a gain sweep of memory capacity: its model property and the names it goes by.
%
%   A = MC_SWEEP_SECOND_AXIS()              % the default axis: the bias, S_c
%   A = MC_SWEEP_SECOND_AXIS('self_feedback')
%   A = MC_SWEEP_SECOND_AXIS(s)             % a run's settings or a job's signature
%
% run_mc_gain_sweep sweeps the gain (level_of_chaos) against ONE other
% scalar model property. Which property is the run's choice (its option
% second_axis); this function is the one place that turns the property into
% the names the stage, its job files, its table and its readers use, so they
% cannot disagree.
%
%   source   the property's name; or a struct that may carry the field
%            second_axis (a run's settings, a job file's signature). Empty,
%            or a struct without the field, is the default axis: results and
%            job files written before the axis could be chosen carry no such
%            field, and their second axis is the bias.
%
%   A.property     the model property swept ('S_c', 'self_feedback', ...)
%   A.column       the name of the axis in a result: the column of the
%                  results table, the field of a job struct, the stem of the
%                  settings fields <column>_levels and preset_<column>, and
%                  the quantity's name in mc_sweep_quantity
%   A.file_token   what precedes the value in a job file's name
%   A.levels_option  run_mc_gain_sweep's option that carries its levels
%   A.is_default   true for the bias
%
% THE BIAS KEEPS ITS OWN NAMES: the property S_c is the column 'bias', the
% file token 'b' and the option 'bias_levels', as in every result, job file
% and bundle made while the bias was the only second axis, so they read,
% resume and plan unchanged. Any other property goes by its own name as
% column and file token, and its levels are the option 'second_levels'.
%
% See also: run_mc_gain_sweep, mc_gain_sweep_job_rows, read_mc_gain_sweep,
%           mc_sweep_quantity

arguments
    source = ''
end

DEFAULT_PROPERTY = 'S_c';

property = DEFAULT_PROPERTY;
if isstruct(source)
    if isfield(source, 'second_axis') && ~isempty(source.second_axis)
        property = char(source.second_axis);
    end
elseif ~isempty(source)
    property = char(source);
end
if ~isvarname(property)
    error('mc_sweep_second_axis:BadProperty', ...
        'The second axis is named by a model property; got ''%s''.', property);
end

A = struct('property', property, 'column', property, 'file_token', property, ...
    'levels_option', 'second_levels', 'is_default', strcmp(property, DEFAULT_PROPERTY));
if A.is_default
    A.column = 'bias';
    A.file_token = 'b';
    A.levels_option = 'bias_levels';
end
end
