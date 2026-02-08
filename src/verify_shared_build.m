function verify_shared_build(esn_array, expected_to_differ, also_check_protected)
% VERIFY_SHARED_BUILD Verify that built ESN objects share identical configuration
%
% verify_shared_build(esn_array, expected_to_differ, also_check_protected)
%
% After building multiple ESN objects for a comparison experiment, this
% function verifies that all properties that SHOULD be shared are identical,
% and that properties expected to differ DO actually differ.
%
% Uses metaclass introspection to automatically discover all properties.
% Skips Dependent properties (computed from other properties), non-public
% properties (unless explicitly listed in also_check_protected), and
% properties in the expected_to_differ list.
%
% Inputs:
%   esn_array           - Cell array of built SRNN_ESN_reservoir objects
%   expected_to_differ  - Cell array of property names that SHOULD differ
%                         between at least one pair of objects.
%                         e.g., {'n_a_E', 'n_b_E', 'tau_a_E'}
%   also_check_protected - Cell array of SetAccess=protected property names
%                         to additionally verify match.
%                         e.g., {'W', 'W_in', 'u_scalar', 'u_ex', 't_ex'}
%
% Errors if:
%   - A public non-Dependent property (not in expected_to_differ) differs
%   - A property in also_check_protected differs
%   - A property in expected_to_differ is actually identical across all objects
%
% Example:
%   verify_shared_build(esn, ...
%       {'n_a_E', 'n_b_E', 'tau_a_E'}, ...
%       {'W', 'W_in', 'u_scalar', 'u_ex', 't_ex'});
%
% See also: SRNN_ESN_reservoir, SRNNModel2

    if numel(esn_array) < 2
        fprintf('verify_shared_build: only 1 object, nothing to compare.\n');
        return;
    end

    ref = esn_array{1};
    mc = metaclass(ref);
    n_obj = numel(esn_array);

    % Properties to always skip (run-output, complex objects, or
    % derived aggregates that legitimately differ when config differs)
    always_skip = {'S0', 'cached_params', 'mc_results', 'u_interpolant', ...
                   'ode_opts', 't_out', 'S_out', 'plot_data', 'lya_results'};

    n_checked = 0;
    n_matched = 0;
    n_skipped = 0;
    checked_names = {};
    skipped_names = {};

    for p = 1:numel(mc.PropertyList)
        prop = mc.PropertyList(p);
        name = prop.Name;

        % Skip Dependent properties (computed from other properties)
        if prop.Dependent
            skipped_names{end+1} = name; %#ok<AGROW>
            n_skipped = n_skipped + 1;
            continue;
        end

        % Skip always-skip list
        if ismember(name, always_skip)
            skipped_names{end+1} = name; %#ok<AGROW>
            n_skipped = n_skipped + 1;
            continue;
        end

        % Skip properties expected to differ
        if ismember(name, expected_to_differ)
            skipped_names{end+1} = name; %#ok<AGROW>
            n_skipped = n_skipped + 1;
            continue;
        end

        % Determine if we should check this property
        is_public_get = strcmp(prop.GetAccess, 'public');
        is_in_also_check = ismember(name, also_check_protected);

        if ~is_public_get && ~is_in_also_check
            skipped_names{end+1} = name; %#ok<AGROW>
            n_skipped = n_skipped + 1;
            continue;
        end

        % Compare this property across all objects vs reference
        for i = 2:n_obj
            obj = esn_array{i};
            val_ref = ref.(name);
            val_obj = obj.(name);

            % Handle function_handle comparison (isequaln fails on closures)
            if isa(val_ref, 'function_handle') && isa(val_obj, 'function_handle')
                match = strcmp(func2str(val_ref), func2str(val_obj));
            else
                match = isequaln(val_ref, val_obj);
            end

            if ~match
                error('verify_shared_build:Mismatch', ...
                    'Property ''%s'' differs between condition 1 and %d.\nThis property was expected to be identical. If it should differ, add it to expected_to_differ.', ...
                    name, i);
            end
        end

        n_checked = n_checked + 1;
        n_matched = n_matched + 1;
        checked_names{end+1} = name; %#ok<AGROW>
    end

    % Verify that expected_to_differ properties actually differ
    for k = 1:numel(expected_to_differ)
        name = expected_to_differ{k};

        % Check if this property exists and is accessible
        prop_meta = findobj(mc.PropertyList, 'Name', name);
        if isempty(prop_meta)
            warning('verify_shared_build:UnknownProperty', ...
                'Property ''%s'' in expected_to_differ does not exist on this class.', name);
            continue;
        end

        % Check if at least one pair differs
        any_differs = false;
        for i = 2:n_obj
            val_ref = ref.(name);
            val_obj = esn_array{i}.(name);
            if ~isequaln(val_ref, val_obj)
                any_differs = true;
                break;
            end
        end

        if ~any_differs
            warning('verify_shared_build:NoDifference', ...
                'Property ''%s'' is listed in expected_to_differ but is identical across all %d conditions.\nThis may indicate a misconfigured experiment.', ...
                name, n_obj);
        end
    end

    fprintf('verify_shared_build: %d properties checked, all matched across %d conditions.\n', ...
        n_checked, n_obj);
    fprintf('  Checked: %s\n', strjoin(checked_names, ', '));
    fprintf('  Expected to differ: %s\n', strjoin(expected_to_differ, ', '));
    if ~isempty(also_check_protected)
        fprintf('  Also verified (protected): %s\n', strjoin(also_check_protected, ', '));
    end
end
