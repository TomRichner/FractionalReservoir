function validate_preset_conditions(conditions, model_class, preset_name)
% VALIDATE_PRESET_CONDITIONS Check a preset's condition set before it is used.
%
%   VALIDATE_PRESET_CONDITIONS(conditions, model_class, preset_name)
%
% Called at the end of srnn_param_preset, so a malformed condition set fails in
% the file you are editing rather than deep in a sweep.
%
% WHY THIS EXISTS. Conditions were the only unvalidated layer of a run.
% model_defaults gets ParamSpaceAnalysis2.validate_model_defaults, which lists
% every problem at once before the output directory is even created; grid axes
% are rejected up front too. A condition field, by contrast, was checked only by
% the model constructor -- which runs inside a parfor worker, mid-sweep, after
% the run directory exists. Now that presets write their conditions out by hand
% rather than receiving them from a shared builder, that gap is where typos land.
%
% DELIBERATELY GENERAL. It knows nothing about tau_a, synapse_config or
% adaptation. run_single_job is equally general -- it appends every field except
% 'name' as a constructor argument -- and a validator narrower than the contract
% it guards would block conditions the pipeline supports. In particular it does
% NOT require a Pairs preset to state tau_a/synapse_config; a preset that omits
% both is accepted, and is visible as such when you read it.
%
% THE FIVE CHECKS
%
%   1. cell array of SCALAR structs. Catches the one real trap in writing
%      conditions by hand: struct('tau_a', taus) with a CELL taus silently builds
%      a 1xN struct ARRAY rather than a scalar struct with a cell field. The
%      double-brace form struct('tau_a', {taus}) is correct but easy to drop, so
%      the preset uses dot assignment instead, which cannot express the mistake.
%   2. 'name' present, non-empty char, and unique across the set. It is the
%      output subdirectory name (ParamSpaceAnalysis2.run) and the key every
%      plotter and title lookup uses, so a duplicate silently merges two regimes.
%   3. ALL CONDITIONS SET THE SAME FIELDS. Only the fields a condition actually
%      sets are condition-owned (see ParamSpaceAnalysis2.condition_set_fields), so
%      one condition omitting a field does not inherit its siblings' value -- it
%      inherits model_defaults, silently running a different experiment under a
%      name that says otherwise. State a mechanism as explicitly off rather than
%      leaving it out.
%   4. every non-'name' field is publicly settable on model_class. Same
%      classification validate_model_defaults uses, and the same all-at-once
%      reporting; Dependent and non-public properties get their own message
%      because "not a property" would be wrong and unhelpful for those.
%   5. every name is a key of srnn_condition_titles. A name that is not there
%      does not fail -- it produces a raw snake_case panel title in a figure,
%      found hours later when the figures are read. This is a "the lookup must
%      exist" check only; the titles file remains the source of the text, which
%      is what lets a saved run be retitled without recomputing it.
%
% Uses meta.class directly rather than ParamSpaceAnalysis2.srnn_property_info,
% which is a PRIVATE static and unreachable from here. That also keeps
% src/presets from depending on src/analysis.
%
% See also: srnn_param_preset, srnn_condition_titles,
%           ParamSpaceAnalysis2/validate_model_defaults

arguments
    conditions
    model_class (1,:) char
    preset_name (1,:) char
end

%% 1. Shape
if ~iscell(conditions)
    error('validate_preset_conditions:NotACell', ...
        'Preset ''%s'': conditions must be a CELL ARRAY of structs, got %s.', ...
        preset_name, class(conditions));
end
if isempty(conditions)
    error('validate_preset_conditions:NoConditions', ...
        'Preset ''%s'': conditions is empty; every preset defines at least one.', ...
        preset_name);
end

problems = {};
for i = 1:numel(conditions)
    c = conditions{i};
    if ~isstruct(c)
        problems{end+1} = sprintf('  condition %d is a %s, not a struct.', ...
            i, class(c)); %#ok<AGROW>
    elseif ~isscalar(c)
        % The struct('field', <cell>) trap -- see the header.
        problems{end+1} = sprintf(['  condition %d is a %dx%d struct ARRAY, not a ' ...
            'scalar struct. struct(''f'', c) with a CELL c builds one element per ' ...
            'cell entry; use dot assignment, or struct(''f'', {c}).'], ...
            i, size(c, 1), size(c, 2)); %#ok<AGROW>
    end
end
raise(problems, preset_name, 'validate_preset_conditions:BadShape', ...
    'condition(s) that are not scalar structs');

%% 2. Names
names = cell(1, numel(conditions));
for i = 1:numel(conditions)
    c = conditions{i};
    if ~isfield(c, 'name')
        problems{end+1} = sprintf('  condition %d has no ''name'' field.', i); %#ok<AGROW>
        names{i} = '';
    elseif ~ischar(c.name) || isempty(c.name)
        problems{end+1} = sprintf(['  condition %d: name must be a non-empty ' ...
            'char row, got %s.'], i, class(c.name)); %#ok<AGROW>
        names{i} = '';
    else
        names{i} = c.name;
    end
end
raise(problems, preset_name, 'validate_preset_conditions:BadName', ...
    'condition(s) with an unusable name');

[uniq, ~, idx] = unique(names);
counts = accumarray(idx(:), 1);
dupes = uniq(counts > 1);
if ~isempty(dupes)
    error('validate_preset_conditions:DuplicateName', ...
        ['Preset ''%s'': duplicate condition name(s): %s.\n' ...
         'The name is the output subdirectory and the key every plotter and ' ...
         'title lookup uses, so duplicates silently merge two regimes.'], ...
        preset_name, strjoin(dupes, ', '));
end

%% 3. Every condition sets the same fields
field_sets = cellfun(@(c) sort(setdiff(fieldnames(c), {'name'}))', ...
    conditions, 'UniformOutput', false);
reference = field_sets{1};
for i = 2:numel(field_sets)
    if isequal(field_sets{i}, reference); continue; end
    missing = setdiff(reference, field_sets{i});
    extra   = setdiff(field_sets{i}, reference);
    detail = '';
    if ~isempty(missing)
        detail = sprintf('%s missing %s;', detail, strjoin(missing, ', '));
    end
    if ~isempty(extra)
        detail = sprintf('%s extra %s;', detail, strjoin(extra, ', '));
    end
    problems{end+1} = sprintf('  ''%s'' vs ''%s'':%s', ...
        names{i}, names{1}, detail); %#ok<AGROW>
end
if ~isempty(problems)
    error('validate_preset_conditions:FieldSetMismatch', ...
        ['Preset ''%s'': conditions do not all set the same fields:\n%s\n\n' ...
         'Only the fields a condition ACTUALLY SETS are condition-owned, so an ' ...
         'omitted field does not fall back to a sibling condition -- it falls ' ...
         'back to model_defaults, running a different experiment under a name ' ...
         'that says otherwise. State a mechanism as explicitly off (tau_a = ' ...
         '{zeros(1,0), ...}, synapse_config = struct()) rather than omitting it.'], ...
        preset_name, strjoin(problems, '\n'));
end

%% 4. Fields are settable on the model class
info = property_info(model_class);
for i = 1:numel(reference)
    name = reference{i};
    if ismember(name, info.settable); continue; end
    if ismember(name, info.dependent)
        problems{end+1} = sprintf(['  ''%s'' is a Dependent (computed) property ' ...
            'of %s and cannot be set; set the properties it is derived from ' ...
            'instead.%s'], name, model_class, suggest(name, info.settable)); %#ok<AGROW>
    elseif ismember(name, info.nonpublic)
        problems{end+1} = sprintf(['  ''%s'' is not publicly settable on %s ' ...
            '(it is computed during build/run).'], name, model_class); %#ok<AGROW>
    else
        problems{end+1} = sprintf('  ''%s'' is not a property of %s.%s', ...
            name, model_class, suggest(name, info.settable)); %#ok<AGROW>
    end
end
if ~isempty(problems)
    error('validate_preset_conditions:UnknownField', ...
        ['Preset ''%s'': %d condition field(s) cannot be set on %s:\n%s\n\n' ...
         'Every field except ''name'' is passed to the constructor by ' ...
         'ParamSpaceAnalysis2.run_single_job.'], ...
        preset_name, numel(problems), model_class, strjoin(problems, '\n'));
end

%% 5. Every name has a title
titles = srnn_condition_titles();
for i = 1:numel(names)
    if ~titles.isKey(names{i})
        problems{end+1} = sprintf('  ''%s''', names{i}); %#ok<AGROW>
    end
end
if ~isempty(problems)
    error('validate_preset_conditions:NoTitle', ...
        ['Preset ''%s'': condition name(s) with no entry in ' ...
         'srnn_condition_titles:\n%s\n\n' ...
         'Without one the plotters fall back to the raw snake_case name as a ' ...
         'panel title, which is found only when the figures are read. Add the ' ...
         'name and its display text there.'], ...
        preset_name, strjoin(problems, '\n'));
end
end

%% ------------------------------------------------------------------------
function raise(problems, preset_name, id, what)
% Report every problem of one kind at once, as validate_model_defaults does: a
% malformed set is usually malformed in the same way several times over, and
% fixing them one error per run is needless.
if isempty(problems); return; end
error(id, 'Preset ''%s'': %d %s:\n%s', ...
    preset_name, numel(problems), what, strjoin(problems, '\n'));
end

function info = property_info(class_name)
% Mirrors ParamSpaceAnalysis2.srnn_property_info, which is a private static and
% so unreachable from here. Kept identical in its classification rules -- notably
% the ischar guard on SetAccess, which treats SetAccess = ?SomeClass as
% non-public -- so a name accepted here is accepted there.
%
% Not cached: this runs once per srnn_param_preset call, not once per result the
% way effective_param does, so the persistent that justifies caching there would
% only add the stale-after-classdef-edit hazard here.
mc = meta.class.fromName(class_name);
if isempty(mc)
    error('validate_preset_conditions:UnknownModelClass', ...
        '''%s'' is not a class on the path.', class_name);
end
props = mc.PropertyList;
names = {props.Name};
is_dependent = [props.Dependent];
has_setter   = ~cellfun(@isempty, {props.SetMethod});
is_public    = cellfun(@(a) ischar(a) && strcmp(a, 'public'), {props.SetAccess});

info = struct();
info.settable  = names(is_public & (~is_dependent | has_setter));
info.dependent = names(is_dependent & ~has_setter);
info.nonpublic = names(~is_public & ~is_dependent);
end

function s = suggest(name, candidates)
% Trailing " Did you mean ...?" text, matching suggest_property's behaviour.
hit = candidates(strcmpi(candidates, name));
if isempty(hit)
    hit = candidates(startsWith(candidates, name, 'IgnoreCase', true) | ...
        contains(candidates, name, 'IgnoreCase', true));
end
if isempty(hit)
    s = '';
elseif isscalar(hit)
    s = sprintf(' Did you mean ''%s''?', hit{1});
else
    s = sprintf(' Did you mean one of: %s?', ...
        strjoin(hit(1:min(5, numel(hit))), ', '));
end
end
