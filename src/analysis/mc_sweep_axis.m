function A = mc_sweep_axis(position, source)
% MC_SWEEP_AXIS One axis of a block C sweep of memory capacity: its model property and the names it goes by.
%
%   A = MC_SWEEP_AXIS('first')                    % the default first axis: the gain
%   A = MC_SWEEP_AXIS('second')                   % the default second axis: the bias, S_c
%   A = MC_SWEEP_AXIS('second', 'self_coupling')
%   A = MC_SWEEP_AXIS(position, s)                % a run's settings or a job's signature
%
% run_mc_gain_sweep sweeps one scalar model property against another on a
% grid. Which two is the run's choice (its options first_axis and
% second_axis); THIS FUNCTION IS THE ONE PLACE that turns an axis's property
% into the names the stage, its job files, its table, its settings, its
% readers, the statistics and the figures use, so they cannot disagree. The
% two axes are the same kind of thing and are described by the same fields;
% only their defaults differ.
%
%   position  'first' (the columns of the maps, the x axis) or 'second' (the
%             rows, the y axis)
%   source    the property's name; or a struct that may carry the field
%             first_axis / second_axis (a run's settings, a job file's
%             signature). Empty, or a struct without the field, is the
%             default axis: results and job files written before an axis
%             could be chosen carry no such field, and their axis is the
%             default.
%
%   A.position       'first' or 'second'
%   A.property       the model property swept ('level_of_chaos', 'S_c',
%                    'self_coupling', 'sigma_u_noise', ...)
%   A.column         the axis's name in a result: the column of the results
%                    table, the field of a job struct, the stem of the
%                    settings fields <column>_levels and preset_<column>,
%                    and the quantity's name in mc_sweep_quantity
%   A.file_token     what precedes the value in a job file's name
%   A.option         the stem of run_mc_gain_sweep's and mc_gain_sweep_grid's
%                    options for this axis: <option>_levels, and for the
%                    first axis <option>_range and <option>_spacing
%   A.generic_option the stem the options have on any axis: 'first' or
%                    'second' (first_levels, second_range, ...)
%   A.levels_option  <option>_levels, the option that carries its levels
%   A.signature_field  the field of a job's signature and a run's settings
%                    that names a non-default axis: first_axis / second_axis
%   A.is_default     true for the default property
%
% THE DEFAULT AXES KEEP THEIR OWN NAMES, as in every result, job file and
% bundle made while they were the only axes, so those read, resume and plan
% unchanged:
%   first   level_of_chaos: column 'gain', file token 'g', options gain_levels,
%           gain_range, gain_spacing, settings gain_levels and preset_gain
%   second  S_c: column 'bias', file token 'b', option bias_levels, settings
%           bias_levels and preset_bias
% ANY OTHER PROPERTY goes by its own name as column and file token, and its
% options are the generic ones (first_levels, first_range, first_spacing;
% second_levels, second_range). The generic options are also accepted on a
% default axis (second_levels for the bias, as before), so a caller that
% names its axes need not know which is a default.
% A job's signature carries first_axis / second_axis ONLY for a non-default
% axis, so the signature of a default run is the one its files have always
% carried.
%
% See also: run_mc_gain_sweep, mc_gain_sweep_grid, mc_gain_sweep_job_rows,
%           read_mc_gain_sweep, mc_sweep_quantity

arguments
    position (1,:) char {mustBeMember(position, {'first', 'second'})}
    source = ''
end

switch position
    case 'first'
        defaults = struct('property', 'level_of_chaos', 'column', 'gain', 'file_token', 'g', ...
            'option', 'gain');
    case 'second'
        defaults = struct('property', 'S_c', 'column', 'bias', 'file_token', 'b', ...
            'option', 'bias');
end
signature_field = [position '_axis'];

property = defaults.property;
if isstruct(source)
    if isfield(source, signature_field) && ~isempty(source.(signature_field))
        property = char(source.(signature_field));
    end
elseif ~isempty(source)
    property = char(source);
end
if ~isvarname(property)
    error('mc_sweep_axis:BadProperty', ...
        'The %s axis is named by a model property; got ''%s''.', position, property);
end

A = struct('position', position, 'property', property, 'column', property, ...
    'file_token', property, 'option', position, 'generic_option', position, ...
    'levels_option', '', 'signature_field', signature_field, ...
    'is_default', strcmp(property, defaults.property));
if A.is_default
    A.column = defaults.column;
    A.file_token = defaults.file_token;
    A.option = defaults.option;
end
A.levels_option = [A.option '_levels'];
end
