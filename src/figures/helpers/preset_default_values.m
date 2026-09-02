function default_value = preset_default_values(data_root, params)
% PRESET_DEFAULT_VALUES Default value of each swept parameter, from the preset.
%
%   default_value = PRESET_DEFAULT_VALUES(data_root, params)
%
% Inputs:
%   data_root : a run_all_<dt> folder (holds run_manifest.mat and the
%               1D_sensitivity_* sweep subfolders)
%   params    : cellstr of parameter names to resolve
%
% Output:
%   default_value : containers.Map name -> scalar default. Parameters that
%                   cannot be resolved to a finite scalar are OMITTED (with a
%                   warning), so callers should test isKey before using one.
%
% Answers "what was this parameter before the sweep moved it" -- the 0% point of
% a percent axis and the position of a default marker.
%
% The run's own resolved_defaults cannot answer it: ParamSpaceAnalysis2
% deliberately EXCLUDES grid axes from that struct, and every parameter asked
% about here is the axis of one of the sweeps. So the source is the preset named
% in the run's own run_manifest.mat, read fresh each time -- re-pointing
% data_root at a run built from a different preset then moves the reference
% automatically instead of silently reusing another network's numbers.
%
% NOTHING IS ASSUMED ABOUT HOW A NAME MAPS TO A PRESET FIELD. Several of these
% are Dependent scalar aliases rather than stored properties: f_E aliases f(1),
% mu_EE_relative aliases mu_tilde_relative(1,1) indexed (post, pre), and a 1 x C
% mu_tilde_relative is a presynaptic row broadcast down the columns. Rather than
% reimplement any of that here -- where it could drift out of step with the
% class -- a model is CONSTRUCTED from the preset and each value is read off it,
% so the class's own getters do the resolving. The model class also comes from
% the manifest, i.e. the class that actually ran.
%
% Falls back to the bare class default, then omits the parameter with a warning.
%
% Extracted verbatim from the local subfunction of
% Fig_sensitivity_analysis_allStd.m so fig_sensitivity_medians can share it.
%
% See also: apply_percent_axis, mark_default_value, srnn_param_preset

    default_value = containers.Map('KeyType', 'char', 'ValueType', 'double');

    manifest_file = fullfile(data_root, 'run_manifest.mat');
    if ~isfile(manifest_file)
        warning('preset_default_values:NoManifest', ...
            ['No run_manifest.mat in %s, so the preset defaults are unknown; ' ...
             'percent axes fall back to raw units and no default markers are ' ...
             'drawn.'], data_root);
        return;
    end
    M = load(manifest_file);
    preset_name = M.run_manifest.preset_name;

    [preset_defaults, preset_class] = srnn_param_preset(preset_name);
    if isfield(M.run_manifest, 'model_class')
        model_class = M.run_manifest.model_class;   % what actually ran
    else
        model_class = preset_class;
    end

    % Build the model the preset describes. Construction only sets properties --
    % build()/run() are never called -- so this is cheap.
    %
    % A preset alone is not always constructible: a preset carrying
    % sigma_u_noise > 0 does NOT carry ode_solver (that comes from
    % analysis_run_config), and the constructor rejects noise on a deterministic
    % solver. So the run's own stored resolved_defaults is layered UNDER the
    % preset -- it supplies exactly those run-mode settings, as they actually
    % ran, while the preset keeps the final say on the physics. resolved_defaults
    % omits each sweep's own grid axis, so it cannot contaminate the values being
    % read back here.
    run_defaults = stored_resolved_defaults(data_root);
    if isempty(fieldnames(run_defaults))
        construct_from = preset_defaults;
    else
        construct_from = merge_struct(run_defaults, preset_defaults);
    end

    nv_merged = namedargs2cell_local(construct_from);
    nv_preset = namedargs2cell_local(preset_defaults);
    probe = [];
    try
        probe = feval(model_class, nv_merged{:});
    catch
        try
            probe = feval(model_class, nv_preset{:});
        catch probe_err
            warning('preset_default_values:PresetProbeFailed', ...
                ['Could not construct %s from preset ''%s'' (%s), so parameter ' ...
                 'defaults fall back to bare class defaults.'], ...
                model_class, preset_name, probe_err.message);
        end
    end

    for k = 1:numel(params)
        name = params{k};
        if ~isempty(probe) && isprop(probe, name)
            val = probe.(name);
        else
            try
                val = ParamSpaceAnalysis2.class_default(name);
            catch
                val = [];
            end
        end
        if isscalar(val) && isnumeric(val) && isfinite(val)
            default_value(name) = val;
        else
            warning('preset_default_values:NoDefault', ...
                ['No scalar default resolved for ''%s'' from preset ''%s''; ' ...
                 'that axis stays in raw units with no default marker.'], ...
                name, preset_name);
        end
    end
end

function rd = stored_resolved_defaults(data_root)
% STORED_RESOLVED_DEFAULTS The parameter set a sweep in this run actually used.
%
% Read from any one sweep's param_space_summary.mat, which ParamSpaceAnalysis2
% documents as the run's metadata record. Only run-mode settings are wanted from
% it (ode_solver above all), so which sweep it comes from does not matter --
% those are identical across the sweeps of one run_all.
%
% Returns an empty struct rather than erroring: this is a convenience layer, and
% the caller falls back to the preset alone.
    rd = struct();
    listing = dir(fullfile(data_root, '1D_sensitivity_*'));
    listing = listing([listing.isdir]);
    for k = 1:numel(listing)
        f = fullfile(listing(k).folder, listing(k).name, 'param_space_summary.mat');
        if ~isfile(f); continue; end
        try
            S = load(f, 'summary_data');
            if isfield(S, 'summary_data') && isfield(S.summary_data, 'resolved_defaults')
                rd = S.summary_data.resolved_defaults;
                return;
            end
        catch
            % try the next sweep
        end
    end
end

function nv = namedargs2cell_local(s)
% Flatten a struct of model overrides into a name-value cell, the same shape
% ParamSpaceAnalysis2 passes to the model constructor.
    f = fieldnames(s);
    nv = cell(1, 2*numel(f));
    nv(1:2:end) = f;
    for i = 1:numel(f)
        nv{2*i} = s.(f{i});
    end
end
