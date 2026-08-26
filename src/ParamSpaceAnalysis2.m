classdef ParamSpaceAnalysis2 < handle
    % PARAMSPACEANALYSIS Multi-dimensional parameter grid search for SRNN
    %
    % This class performs parameter space exploration by generating all
    % combinations of specified parameters, simulating each configuration
    % across multiple adaptation conditions, and analyzing the results.
    %
    % Key features:
    %   - Multi-dimensional grid (not 1D sweeps like SensitivityAnalysis)
    %   - Same network W used across all 4 conditions for fair comparison
    %   - Randomized execution order for representative early-stopping
    %   - Batched parfor with checkpoint files for resume capability
    %   - Histogram-based visualization of results
    %
    % Usage:
    %   psa = ParamSpaceAnalysis2('n_levels', 4, 'note', 'demo');
    %   psa.add_grid_parameter('level_of_chaos', [0.5, 3.0]);
    %   psa.add_grid_parameter('n', [50, 200]);
    %   psa.run();
    %   psa.plot('metric', 'LLE');
    %
    % See also: SRNNModel2, SensitivityAnalysis

    %% Configuration Properties
    properties
        grid_params = {}            % Cell array of parameter names for grid
        param_ranges = struct()     % Struct: param_name -> [min, max]
        n_levels = 8                % Number of levels per grid parameter
        conditions                  % Cell array of condition structs
        integer_params = {'n', 'indegree', 'n_a_E', 'n_a_I', 'n_b_E', 'n_b_I'}
        reps                        % Repetition index for multiple runs per grid point
        explicit_vectors = struct() % Struct: param_name -> explicit vector (when length > 2)
        vector_param_config = struct() % Struct: param_name -> config for vector params
        randomize_order = true         % Whether to randomize execution order (false for sensitivity)
        network_seed_offset = 0        % Added to every per-config network seed; set per run (e.g. run_index*1e6) so repeated runs of the same config draw independent networks for pooling

        % Fraction of the full grid to actually simulate, in (0, 1].
        %
        % A high-dimensional grid costs n_levels^n_params, which grows past
        % what an overnight run can hold long before the QUESTION stops being
        % answerable: the marginals and the pooled histograms are estimated
        % from a sample, not enumerated. subset_fraction = 0.15 runs a random
        % 15% of the grid points and leaves the rest empty.
        %
        % Gaps are safe by construction and this is worth stating, because it
        % is the whole reason no other code needed changing:
        %   - results are stored by GRID POSITION (config_idx), never by
        %     execution order, so a partial run puts every result at its true
        %     coordinates and leaves holes;
        %   - the network seed is config_idx*100 + offset, tied to position
        %     rather than order, so a grid point draws the SAME network whether
        %     it was the 3rd or the 300th thing run -- a subset is a subset of
        %     the same experiment, not a different one;
        %   - every pooling/plotting path filters on isstruct(res) && res.success,
        %     so empty slots are skipped exactly like failed ones already were.
        %
        % REQUIRES randomize_order = true, enforced in generate_grid. A subset
        % of a SEQUENTIAL order is the first K configs in ndgrid order, i.e. a
        % systematic corner of the grid (smallest n, smallest f, weakest gain),
        % which is not a sample of anything. That is why the sensitivity sweeps
        % -- which need randomize_order = false for their ordered axes -- cannot
        % use this, and do not need to: they are 1-D.
        subset_fraction = 1            % Fraction of grid points to run, (0, 1]. Needs randomize_order.
    end

    %% Model Default Properties
    properties
        % Which model class to sweep. Any class constructible from name-value
        % pairs with build()/run() methods works -- PSA needs a constructor and
        % a property list, not a common base class, so the two model classes are
        % duck-typed rather than sharing a hierarchy.
        model_class = 'SRNNModel2'  % e.g. 'SRNNModel2' | 'SRNNCellTypePairs'
        model_defaults = struct()   % Default properties of model_class
        verbose = true              % Print progress during execution
    end

    %% Execution Properties
    properties
        batch_size = 25             % Number of configs per batch
        output_dir                  % Base directory for saving results
        note = ''                   % Optional note for folder naming
        folder_prefix = 'param_space'  % Prefix for output folder name
        store_local_lya = false     % Whether to store decimated local Lyapunov time series
        store_local_lya_dt = 0.1    % Time resolution for stored local_lya (seconds)
        use_parallel = true         % Whether to use parfor (requires Parallel Computing Toolbox)
    end

    %% Results (SetAccess = private)
    properties (SetAccess = private)
        results = struct()          % Stored results after run
        has_run = false             % Flag indicating if analysis has run
        analysis_start_time         % Timestamp when analysis started
        param_vectors               % Cell array of parameter level vectors
        all_configs                 % Cell array of all config structs
        shuffled_indices            % Randomized/sequential order for execution
        num_combinations            % Total number of grid points
        vector_param_lookup = struct() % Struct: param_name -> cell array of pre-generated vectors
        resolved_defaults = struct()   % Every SRNNModel2 parameter as actually used (set by run); empty for runs predating this field
    end

    %% Constructor
    methods
        function obj = ParamSpaceAnalysis2(varargin)
            % PARAMSPACEANALYSIS Constructor with name-value pairs
            %
            % Usage:
            %   psa = ParamSpaceAnalysis2()
            %   psa = ParamSpaceAnalysis2('n_levels', 8, 'batch_size', 50)

            % Set default conditions
            obj.set_default_conditions();

            % Parse name-value pairs
            for i = 1:2:length(varargin)
                if isprop(obj, varargin{i})
                    obj.(varargin{i}) = varargin{i+1};
                else
                    warning('ParamSpaceAnalysis2:UnknownProperty', ...
                        'Unknown property: %s', varargin{i});
                end
            end

            % Set default output directory
            if isempty(obj.output_dir)
                % Default to 'data/param_space' in the project root
                src_path = fileparts(mfilename('fullpath'));
                project_root = fileparts(src_path);
                obj.output_dir = fullfile(project_root, 'data', 'param_space');
            end
        end
    end

    %% Public Methods
    methods
        function add_grid_parameter(obj, param_name, param_range)
            % ADD_GRID_PARAMETER Add a parameter to the multi-dimensional grid
            %
            % Usage:
            %   psa.add_grid_parameter('level_of_chaos', [0.5, 3.0])

            if ~ischar(param_name) && ~isstring(param_name)
                error('ParamSpaceAnalysis2:InvalidInput', ...
                    'param_name must be a string or char array');
            end

            if ~isnumeric(param_range) || length(param_range) < 2
                error('ParamSpaceAnalysis2:InvalidInput', ...
                    'param_range must be a numeric array with at least 2 elements');
            end

            ParamSpaceAnalysis2.reject_condition_field(param_name);

            % Add to grid_params if not already present
            if ~ismember(param_name, obj.grid_params)
                obj.grid_params{end+1} = param_name;
            end

            if length(param_range) == 2
                % Range mode: [min, max] -> will use n_levels
                if param_range(2) < param_range(1)
                    error('ParamSpaceAnalysis2:InvalidInput', ...
                        'param_range(2) must be >= param_range(1)');
                end
                obj.param_ranges.(param_name) = param_range;
                % Remove from explicit_vectors if previously set
                if isfield(obj.explicit_vectors, param_name)
                    obj.explicit_vectors = rmfield(obj.explicit_vectors, param_name);
                end
                if obj.verbose
                    fprintf('Added grid parameter: %s, range: [%.3g, %.3g] (n_levels=%d)\n', ...
                        param_name, param_range(1), param_range(2), obj.n_levels);
                end
            else
                % Explicit vector mode: use values directly
                obj.explicit_vectors.(param_name) = param_range;
                % Remove from param_ranges if previously set
                if isfield(obj.param_ranges, param_name)
                    obj.param_ranges = rmfield(obj.param_ranges, param_name);
                end
                if obj.verbose
                    fprintf('Added grid parameter: %s, explicit vector with %d values\n', ...
                        param_name, length(param_range));
                end
            end
        end

        function remove_grid_parameter(obj, param_name)
            % REMOVE_GRID_PARAMETER Remove a parameter from the grid

            idx = find(strcmp(obj.grid_params, param_name));
            if ~isempty(idx)
                obj.grid_params(idx) = [];
                if isfield(obj.param_ranges, param_name)
                    obj.param_ranges = rmfield(obj.param_ranges, param_name);
                end
                if obj.verbose
                    fprintf('Removed grid parameter: %s\n', param_name);
                end
            else
                warning('ParamSpaceAnalysis2:ParamNotFound', ...
                    'Parameter %s not found in grid', param_name);
            end
        end

        function add_vector_parameter(obj, param_name, varargin)
            % ADD_VECTOR_PARAMETER Add a vector-valued parameter to the grid
            %
            % For parameters like tau_a_E where the parameter is a vector and
            % one end is varied across levels. Internally, the grid uses integer
            % indices that map to pre-generated vectors.
            %
            % Usage:
            %   psa.add_vector_parameter('tau_a_E', ...
            %       'vary_element', 'last', ...
            %       'fixed_value', 0.25, ...
            %       'vary_range', [5, 60], ...
            %       'n_elements', 3, ...
            %       'spacing', 'log', ...
            %       'level_spacing', 'linear')

            % Parse name-value pairs
            p = inputParser;
            addRequired(p, 'param_name', @(x) ischar(x) || isstring(x));
            addParameter(p, 'vary_element', 'last', @(x) ismember(x, {'first', 'last'}));
            addParameter(p, 'fixed_value', [], @isnumeric);
            addParameter(p, 'vary_range', [], @(x) isnumeric(x) && length(x) == 2);
            addParameter(p, 'n_elements', [], @(x) isnumeric(x) && x >= 2);
            addParameter(p, 'spacing', 'linear', @(x) ismember(x, {'linear', 'log'}));
            addParameter(p, 'level_spacing', 'linear', @(x) ismember(x, {'linear', 'log'}));
            parse(p, param_name, varargin{:});

            % Validate required parameters
            if isempty(p.Results.fixed_value)
                error('ParamSpaceAnalysis2:InvalidInput', 'fixed_value is required');
            end
            if isempty(p.Results.vary_range)
                error('ParamSpaceAnalysis2:InvalidInput', 'vary_range is required');
            end
            if isempty(p.Results.n_elements)
                error('ParamSpaceAnalysis2:InvalidInput', 'n_elements is required');
            end
            if p.Results.vary_range(2) < p.Results.vary_range(1)
                error('ParamSpaceAnalysis2:InvalidInput', 'vary_range(2) must be >= vary_range(1)');
            end

            ParamSpaceAnalysis2.reject_condition_field(param_name);

            % Store config
            vpc = struct();
            vpc.vary_element = p.Results.vary_element;
            vpc.fixed_value = p.Results.fixed_value;
            vpc.vary_range = p.Results.vary_range;
            vpc.n_elements = p.Results.n_elements;
            vpc.spacing = p.Results.spacing;
            vpc.level_spacing = p.Results.level_spacing;
            obj.vector_param_config.(param_name) = vpc;

            % Add to grid_params if not already present
            if ~ismember(param_name, obj.grid_params)
                obj.grid_params{end+1} = param_name;
            end

            % Remove from scalar param storage if previously set
            if isfield(obj.param_ranges, param_name)
                obj.param_ranges = rmfield(obj.param_ranges, param_name);
            end
            if isfield(obj.explicit_vectors, param_name)
                obj.explicit_vectors = rmfield(obj.explicit_vectors, param_name);
            end

            if obj.verbose
                fprintf('Added vector parameter: %s, vary_%s [%.3g, %.3g], %d elements, %s spacing, %s level_spacing\n', ...
                    param_name, vpc.vary_element, vpc.vary_range(1), vpc.vary_range(2), ...
                    vpc.n_elements, vpc.spacing, vpc.level_spacing);
            end
        end

        function set_conditions(obj, conditions_cell)
            % SET_CONDITIONS Set custom conditions for the analysis
            %
            % Usage:
            %   psa.set_conditions({
            %       struct('name', 'no_adaptation', 'n_a_E', 0, 'n_b_E', 0),
            %       struct('name', 'sfa_only',      'n_a_E', 3, 'n_b_E', 0)
            %   })

            obj.conditions = conditions_cell;
            if obj.verbose
                fprintf('Set %d custom conditions\n', length(conditions_cell));
            end
        end

        function validate_noise_settings(obj)
            % VALIDATE_NOISE_SETTINGS Reject a swept sigma_u_noise paired with a
            % deterministic integrator, before the sweep starts.
            %
            % The model would catch this too, but only when it reached the first
            % NONZERO grid point -- which for a sweep starting at sigma = 0 can
            % be an hour in, after a partial run directory has been written.
            % Runs next to validate_model_defaults, before the output directory
            % exists.

            name = 'sigma_u_noise';
            sigma_levels = [];
            if isfield(obj.model_defaults, name)
                sigma_levels = obj.model_defaults.(name)(:);
            end
            % Grid levels live in param_ranges ([min max], expanded to n_levels)
            % or in explicit_vectors; the largest value is what decides.
            if isfield(obj.param_ranges, name)
                sigma_levels = [sigma_levels; obj.param_ranges.(name)(:)];
            end
            if isfield(obj.explicit_vectors, name)
                sigma_levels = [sigma_levels; obj.explicit_vectors.(name)(:)];
            end
            if isempty(sigma_levels) || ~any(sigma_levels > 0)
                return;     % nothing asks for noise
            end

            solver = 'ode45';   % the class default, unless overridden
            if isfield(obj.model_defaults, 'ode_solver')
                solver = obj.model_defaults.ode_solver;
            end
            % Delegate so the rule and its message live in one place. Report the
            % largest requested level, which is the one that will certainly bite.
            SRNNModel2.check_noise_settings(max(sigma_levels), solver, 'ParamSpaceAnalysis2');
        end

        function validate_model_defaults(obj)
            % VALIDATE_MODEL_DEFAULTS Check model_defaults against SRNNModel2's properties
            %
            % Errors on names that can never take effect: unknown properties
            % (typos), Dependent properties, and properties that are not publicly
            % settable. All such problems are accumulated into a SINGLE error so a
            % bad override struct can be fixed in one pass.
            %
            % Warns only on names a CONDITION sets, which can therefore never take
            % effect in any run. Names shadowed by a grid axis are NOT warned
            % about: a parameter preset legitimately carries values for axes that
            % some sweeps vary and others hold fixed, so that case is expected and
            % is merely reported by run() (see report_shadowed_defaults).
            %
            % This matters because SRNNModel2's constructor only WARNS on an
            % unknown name, and a warning raised inside a parfor worker is easily
            % lost; setting a Dependent or protected property throws mid-sweep.
            % Called automatically at the top of run(); call it directly to check a
            % configuration before committing to a long sweep.

            fields = fieldnames(obj.model_defaults);
            if isempty(fields); return; end

            info = ParamSpaceAnalysis2.srnn_property_info(obj.model_class);
            condition_fields = obj.condition_set_fields();
            problems = {};

            for i = 1:numel(fields)
                name = fields{i};
                if ismember(name, info.settable)
                    % Valid property. A grid axis overriding it is expected; only a
                    % condition-set name can never take effect anywhere.
                    if ismember(name, condition_fields) && ~ismember(name, obj.grid_params)
                        warning('ParamSpaceAnalysis2:ConditionParamShadowed', ...
                            ['model_defaults.%s is set by every condition; the condition ' ...
                            'value wins and the default can never take effect.'], name);
                    end
                elseif ismember(name, info.dependent)
                    problems{end+1} = sprintf(['  ''%s'' is a Dependent (computed) property ' ...
                        'of %s and cannot be set; set the properties it is ' ...
                        'derived from instead.%s'], ...
                        name, obj.model_class, ...
                        ParamSpaceAnalysis2.suggest_property(name, info.settable)); %#ok<AGROW>
                elseif ismember(name, info.nonpublic)
                    problems{end+1} = sprintf(['  ''%s'' is not publicly settable on ' ...
                        '%s (it is computed during build/run).'], name, obj.model_class); %#ok<AGROW>
                else
                    problems{end+1} = sprintf('  ''%s'' is not a property of %s.%s', ...
                        name, obj.model_class, ...
                        ParamSpaceAnalysis2.suggest_property(name, info.settable)); %#ok<AGROW>
                end
            end

            if ~isempty(problems)
                error('ParamSpaceAnalysis2:InvalidModelDefaults', ...
                    'model_defaults contains %d field(s) that cannot take effect:\n%s', ...
                    numel(problems), strjoin(problems, '\n'));
            end
        end

        function resolve_model_defaults(obj)
            % RESOLVE_MODEL_DEFAULTS Freeze this run's full SRNNModel2 parameter set
            %
            % Records every publicly settable SRNNModel2 property as it will
            % actually be used, so a saved run describes itself and readers do
            % not have to re-derive run_single_job's precedence. Grid parameters
            % and condition fields are excluded: those vary per job and are
            % already recorded per result.
            %
            % Built by CONSTRUCTING a model rather than merging structs, so the
            % set_defaults side effects are captured too -- the activation
            % handles, input_config, and plot_deci computed from fs/plot_freq.
            %
            % The RMT parameters freeze cleanly now that they are stored as
            % _relative multipliers of F: mu_E_tilde and friends are Dependent,
            % so they drop out of srnn_property_info().settable automatically,
            % and mu_E_tilde_relative and friends are captured like any other
            % scalar. F itself is reconstructed per grid point from n/indegree
            % (or pinned, when F_tracks_network is false -- also recorded here).

            excluded = ParamSpaceAnalysis2.per_job_param_names( ...
                obj.grid_params, obj.condition_set_fields());

            model_args = {};
            default_fields = fieldnames(obj.model_defaults);
            for i = 1:numel(default_fields)
                fname = default_fields{i};
                if ~ismember(fname, excluded)
                    model_args = [model_args, {fname, obj.model_defaults.(fname)}]; %#ok<AGROW>
                end
            end

            m = feval(obj.model_class, model_args{:});
            info = ParamSpaceAnalysis2.srnn_property_info(obj.model_class);
            obj.resolved_defaults = struct();
            for i = 1:numel(info.settable)
                pname = info.settable{i};
                if ~ismember(pname, excluded)
                    obj.resolved_defaults.(pname) = m.(pname);
                end
            end
        end

        function report_shadowed_defaults(obj)
            % REPORT_SHADOWED_DEFAULTS Note which model_defaults the grid supersedes
            %
            % Informational only, and not a warning: a parameter preset carries a
            % value for every axis, and each sweep varies a different subset, so
            % this situation is normal. Printing it once keeps the fact visible
            % without the alarm fatigue a per-field warning would cause.
            if ~obj.verbose; return; end
            shadowed = intersect(fieldnames(obj.model_defaults), obj.grid_params);
            if ~isempty(shadowed)
                shadowed = sort(shadowed);
                fprintf('[psa] model_defaults superseded by grid axes: %s\n', ...
                    strjoin(shadowed(:)', ', '));
            end
        end

        function names = condition_set_fields(obj)
            % CONDITION_SET_FIELDS Condition fields at least one condition sets.
            %
            % Only these are genuinely supplied per job. The default conditions
            % set n_a_E and n_b_E but never n_a_I / n_b_I, so treating all four as
            % condition-owned would silently drop a model_defaults value for the
            % inhibitory ones (which is exactly the bug this replaces).
            names = {};
            for c_idx = 1:numel(obj.conditions)
                names = union(names, setdiff(fieldnames(obj.conditions{c_idx}), {'name'}));
            end
            names = names(:)';
        end

        function val = effective_param(obj, res, name)
            % EFFECTIVE_PARAM Value a parameter actually had for one result
            %
            % Resolves NAME with the same precedence run_single_job uses to build
            % the model, so readers never have to reimplement it:
            %   1. vector grid parameter -- res.config holds the level INDEX, not
            %      the value, so it is resolved through vector_param_lookup
            %   2. scalar grid parameter -- res.config
            %   3. condition field (n_a_E / n_b_E / ...), via res.condition_name
            %   4. obj.resolved_defaults, the set frozen at run time
            %   5. obj.model_defaults, then the SRNNModel2 class default -- the
            %      fallback for objects loaded from runs predating (4)
            %
            % Pass RES = [] for a parameter with no per-result value, to resolve
            % straight from model_defaults / the class default.
            %
            % Usage:
            %   f = psa.effective_param(res, 'f');
            %   T = psa.effective_param([], 'T_range');

            info = ParamSpaceAnalysis2.srnn_property_info(obj.model_class);
            if ~ismember(name, info.settable)
                error('ParamSpaceAnalysis2:UnknownModelParam', ...
                    '''%s'' is not a settable property of %s.%s', ...
                    name, obj.model_class, ...
                    ParamSpaceAnalysis2.suggest_property(name, info.settable));
            end

            % 1-2. Per-config grid value
            if isstruct(res) && isfield(res, 'config') && isfield(res.config, name)
                if isfield(obj.vector_param_lookup, name)
                    val = obj.vector_param_lookup.(name){res.config.(name)};
                elseif isfield(obj.vector_param_config, name)
                    % A VECTOR parameter whose decoder is missing. res.config
                    % holds a grid LEVEL INDEX here, not a value, so returning it
                    % would hand back 1..n_levels dressed up as seconds. Fail
                    % loudly instead. This guard is here because that once
                    % happened for real (fixed 2026-08-14) and was only caught
                    % because a tau sweep quietly plotted against 1..13.
                    error('ParamSpaceAnalysis2:MissingVectorLookup', ...
                        ['''%s'' is a vector parameter, so res.config.%s holds a ' ...
                         'grid LEVEL INDEX (%s) rather than a value, and ' ...
                         'vector_param_lookup is not populated to decode it.\n' ...
                         'Load the run with ParamSpaceAnalysis2.from_dir(dir), ' ...
                         'which restores the lookup from psa_object.mat.'], ...
                        name, name, mat2str(res.config.(name)));
                else
                    val = res.config.(name);
                end
                return;
            end

            % 3. Per-condition value
            if isstruct(res) && isfield(res, 'condition_name')
                for c_idx = 1:numel(obj.conditions)
                    cond = obj.conditions{c_idx};
                    if strcmp(cond.name, res.condition_name)
                        if isfield(cond, name)
                            val = cond.(name);
                            return;
                        end
                        break;
                    end
                end
            end

            % 4-5. Run-wide default: prefer the set frozen at run time, which
            % already folds in the class defaults. The model_defaults path is
            % the fallback for runs saved before resolved_defaults existed.
            if isfield(obj.resolved_defaults, name)
                val = obj.resolved_defaults.(name);
            elseif isfield(obj.model_defaults, name)
                val = obj.model_defaults.(name);
            else
                % Last resort: the class's own default. Some model classes have
                % required constructor arguments and cannot be default-built, in
                % which case there is no class default to report.
                try
                    val = ParamSpaceAnalysis2.class_default(name, obj.model_class);
                catch
                    error('ParamSpaceAnalysis2:NoValueForParam', ...
                        ['''%s'' is not set by the grid, the conditions, or ' ...
                        'model_defaults, and %s cannot be default-constructed to ' ...
                        'supply a class default. Set it in model_defaults.'], ...
                        name, obj.model_class);
                end
            end
        end

        function model = rebuild_model(obj, res)
            % REBUILD_MODEL Reconstruct the model a result came from, unbuilt.
            %
            %   model = psa.rebuild_model(res)
            %
            % Returns a CONSTRUCTED model -- build() and run() are the caller's
            % to invoke, so a caller that only wants the connectivity pays for
            % build() and nothing more.
            %
            % The reconstruction is exact because a job's network is pinned by
            % data the result carries: run_single_job passes
            % rng_seeds = [network_seed, network_seed + 1], and network_seed is
            % config_idx*100 + network_seed_offset, i.e. tied to the grid
            % POSITION. res.network_seed records it, so calling build() on what
            % this returns reproduces the same W the sweep simulated -- also
            % across a subsetted run, where a point's position is unchanged by
            % how many of its neighbours were skipped.
            %
            % Argument precedence deliberately MIRRORS run_single_job's, weakest
            % first, letting the constructor's last-write-wins settle collisions:
            %
            %   model_defaults  <  condition  <  grid parameters
            %
            % Keep the two in step. This reads model_defaults, not
            % resolved_defaults: resolved_defaults is model_defaults already
            % pushed through the constructor, so feeding it back would re-assign
            % every class default explicitly -- harmless in principle, but it
            % would no longer be the same call the sweep made.
            %
            % Used by fig_EI_weights_param_space, which needs each grid point's
            % actual weight matrix and cannot get it from the stored results:
            % a result records scalars (LLE, mean rate), never W.
            if ~isstruct(res) || ~isfield(res, 'config_idx')
                error('ParamSpaceAnalysis2:BadResult', ...
                    'rebuild_model needs a result struct with a config_idx.');
            end
            if ~isfield(res, 'network_seed')
                error('ParamSpaceAnalysis2:NoNetworkSeed', ...
                    ['Result %d carries no network_seed, so its network cannot ' ...
                     'be reproduced. Runs predating that field are out of reach.'], ...
                    res.config_idx);
            end

            model_args = {'rng_seeds', [res.network_seed, res.network_seed + 1]};

            default_fields = fieldnames(obj.model_defaults);
            for k = 1:numel(default_fields)
                model_args = [model_args, ...
                    {default_fields{k}, obj.model_defaults.(default_fields{k})}]; %#ok<AGROW>
            end

            % The condition this result ran under, by name.
            if isfield(res, 'condition_name')
                for c = 1:numel(obj.conditions)
                    if strcmp(obj.conditions{c}.name, res.condition_name)
                        cond_fields = setdiff(fieldnames(obj.conditions{c}), {'name'})';
                        for cf = 1:numel(cond_fields)
                            model_args = [model_args, ...
                                {cond_fields{cf}, obj.conditions{c}.(cond_fields{cf})}]; %#ok<AGROW>
                        end
                        break;
                    end
                end
            end

            % Grid parameters. A vector parameter stores a LEVEL INDEX in
            % res.config, so it goes through the lookup rather than being passed
            % as the number it appears to be.
            if isfield(res, 'config')
                gp = fieldnames(res.config);
                for k = 1:numel(gp)
                    if isfield(obj.vector_param_lookup, gp{k})
                        v = obj.vector_param_lookup.(gp{k}){res.config.(gp{k})};
                    else
                        v = res.config.(gp{k});
                    end
                    model_args = [model_args, {gp{k}, v}]; %#ok<AGROW>
                end
            end

            model = feval(obj.model_class, model_args{:});
        end

        function run(obj)
            % RUN Execute the full parameter space analysis
            %
            % This method:
            %   1. Generates the multi-dimensional parameter grid
            %   2. Randomizes execution order
            %   3. Runs batched parfor with checkpoint files
            %   4. Consolidates results into per-condition MAT files

            % Validate
            if isempty(obj.grid_params)
                error('ParamSpaceAnalysis2:NoParameters', ...
                    'No grid parameters defined. Use add_grid_parameter() first.');
            end

            if isempty(obj.conditions)
                error('ParamSpaceAnalysis2:NoConditions', ...
                    'No conditions defined.');
            end

            if ~ischar(obj.model_class) || isempty(meta.class.fromName(obj.model_class))
                error('ParamSpaceAnalysis2:UnknownModelClass', ...
                    'model_class must name a class on the path; got ''%s''.', ...
                    char(string(obj.model_class)));
            end

            % Reject a grid axis that the conditions own -- both the names this
            % class reserves and any field the configured conditions actually
            % set, since a grid value would silently override the condition.
            % Checked here as well as in add_grid_parameter/add_vector_parameter
            % because grid_params is publicly settable and could be assigned
            % directly or restored from an older saved object.
            owned = intersect(obj.grid_params, ...
                [ParamSpaceAnalysis2.condition_field_names(), obj.condition_set_fields()]);
            if ~isempty(owned)
                error('ParamSpaceAnalysis2:ConditionFieldAsGridParam', ...
                    ['%s cannot be a grid parameter: the adaptation conditions set ' ...
                    'it per job, and a grid value would silently override every ' ...
                    'condition. Vary it with set_conditions() instead.'], ...
                    strjoin(owned, ', '));
            end

            % Reject unusable model_defaults BEFORE the output directory is
            % created, so a rejected config leaves no empty dated folder behind.
            obj.validate_model_defaults();
            obj.validate_noise_settings();
            obj.report_shadowed_defaults();

            % Create timestamped output directory
            obj.analysis_start_time = datetime('now');
            dt_str = lower(datestr(obj.analysis_start_time, 'mmm_dd_yy_HH_MM'));

            if ~isempty(obj.note)
                folder_name = sprintf('%s_%s_nLevs_%d_%s', ...
                    obj.folder_prefix, obj.note, obj.n_levels, dt_str);
            else
                folder_name = sprintf('%s_nLevs_%d_%s', ...
                    obj.folder_prefix, obj.n_levels, dt_str);
            end

            obj.output_dir = fullfile(obj.output_dir, folder_name);
            if ~exist(obj.output_dir, 'dir')
                mkdir(obj.output_dir);
            end

            % Generate parameter grid
            obj.generate_grid();

            % Freeze the full parameter set this run will actually use, so the
            % saved run describes itself without a reader re-deriving it.
            obj.resolve_model_defaults();

            % Print summary
            fprintf('\n========================================\n');
            fprintf('=== SRNN Parameter Space Analysis ===\n');
            fprintf('========================================\n');
            fprintf('Grid parameters: %s\n', strjoin(obj.grid_params, ', '));
            fprintf('Levels per parameter: %d\n', obj.n_levels);
            fprintf('Total grid combinations: %d\n', obj.num_combinations);
            fprintf('Conditions: %s\n', strjoin(cellfun(@(c) c.name, obj.conditions, 'UniformOutput', false), ', '));
            fprintf('Batch size: %d\n', obj.batch_size);
            fprintf('Output directory: %s\n', obj.output_dir);
            fprintf('========================================\n\n');

            % Create temp directory for batch results
            temp_dir = fullfile(obj.output_dir, 'temp_batches');
            if ~exist(temp_dir, 'dir')
                mkdir(temp_dir);
            end

            % Create condition output directories
            for c_idx = 1:length(obj.conditions)
                cond_dir = fullfile(obj.output_dir, obj.conditions{c_idx}.name);
                if ~exist(cond_dir, 'dir')
                    mkdir(cond_dir);
                end
            end

            % Write the object BEFORE any simulation runs, results still empty.
            % psa_object.mat is the authoritative record of a run, and writing it
            % here is what makes a run that dies part-way recoverable: temp_batches
            % holds the results, and this holds the configuration needed to
            % interpret them. Overwritten with the completed object below.
            obj.save_object();

            overall_start = tic;

            % Run batched simulation
            obj.run_batched_simulation(temp_dir);

            % Consolidate batch results (internal call - skip validation/recovery)
            obj.consolidate(temp_dir);

            overall_elapsed = toc(overall_start);
            fprintf('\n========================================\n');
            fprintf('=== Analysis Complete ===\n');
            fprintf('Total time: %.2f hours\n', overall_elapsed/3600);
            fprintf('========================================\n');

            obj.has_run = true;

            % Save summary (metadata) and the completed object (authoritative)
            obj.save_summary();
            obj.save_object();
        end

        function save_object(obj)
            % SAVE_OBJECT Write psa_object.mat, the authoritative record of a run.
            %
            % ALWAYS under the variable name `psa`. This used to be each calling
            % script's job, and they disagreed: most saved `psa` but the tau
            % script saved `psa_tau_a` / `psa_tau_b`, so every reader had to know
            % about it and three of five got it wrong. from_dir selects by CLASS
            % rather than by name, so old files still load, but nothing new should
            % add to the problem.
            if isempty(obj.output_dir)
                error('ParamSpaceAnalysis2:NoOutputDir', ...
                    'output_dir must be set before saving the object.');
            end
            if ~exist(obj.output_dir, 'dir')
                mkdir(obj.output_dir);
            end
            psa = obj; %#ok<NASGU>
            save(fullfile(obj.output_dir, 'psa_object.mat'), 'psa');
        end

        function [tf, reason] = same_config(obj, other, varargin)
            % SAME_CONFIG Whether another PSA run is poolable with this one.
            %   [tf, reason] = obj.same_config(other) returns true when the two
            %   runs sweep the same grid under the same conditions and model
            %   parameters, so their results can be concatenated for combined
            %   plotting. The reps-axis LENGTH, randomize_order,
            %   network_seed_offset and subset_fraction are intentionally
            %   ignored (expected to differ across runs) -- pooling two
            %   different random subsets of the same grid is precisely what
            %   subsetting is for. `reason` explains the first mismatch found.
            %
            %   Model parameters are compared through resolved_defaults -- the
            %   FULL parameter set each run froze at run time -- so two runs with
            %   identical explicit overrides but different SRNNModel2 class
            %   defaults (say, across a commit that changed the default
            %   activation) are correctly refused.
            %
            %   Name-value:
            %     'allow_legacy' - (default false) when a run predates
            %                      resolved_defaults, fall back to comparing
            %                      model_defaults instead of refusing. Weaker:
            %                      it cannot see class defaults.
            tf = false; reason = '';

            allow_legacy = false;
            for i = 1:2:numel(varargin)
                switch lower(varargin{i})
                    case 'allow_legacy', allow_legacy = varargin{i+1};
                end
            end

            % Different model classes are never poolable, whatever the parameters.
            if ~strcmp(obj.model_class, other.model_class)
                reason = sprintf('model_class differs: %s vs %s', ...
                    obj.model_class, other.model_class);
                return;
            end

            % Grid parameter names (order-independent set).
            if ~isequal(sort(obj.grid_params(:)), sort(other.grid_params(:)))
                reason = sprintf('grid_params differ: {%s} vs {%s}', ...
                    strjoin(obj.grid_params, ','), strjoin(other.grid_params, ','));
                return;
            end

            if ~isequal(obj.n_levels, other.n_levels)
                reason = sprintf('n_levels differ: %d vs %d', obj.n_levels, other.n_levels);
                return;
            end

            % Swept (non-'reps') grid points must match exactly.
            swept = setdiff(obj.grid_params, {'reps'}, 'stable');
            for i = 1:numel(swept)
                p = swept{i};
                va = obj.param_vectors{strcmp(obj.grid_params, p)};
                vb = other.param_vectors{strcmp(other.grid_params, p)};
                if ~isequaln(va, vb)
                    reason = sprintf('grid points for ''%s'' differ', p);
                    return;
                end
            end

            % Conditions match as a set (compared by full struct content).
            if numel(obj.conditions) ~= numel(other.conditions)
                reason = 'number of conditions differ';
                return;
            end
            for i = 1:numel(obj.conditions)
                ci = obj.conditions{i};
                found = false;
                for j = 1:numel(other.conditions)
                    if isequaln(ci, other.conditions{j}); found = true; break; end
                end
                if ~found
                    reason = sprintf('condition ''%s'' not found in other run', ci.name);
                    return;
                end
            end

            % Integer params / vector configs.
            if ~isequal(sort(obj.integer_params(:)), sort(other.integer_params(:)))
                reason = 'integer_params differ'; return;
            end
            if ~isequaln(obj.explicit_vectors, other.explicit_vectors)
                reason = 'explicit_vectors differ'; return;
            end
            if ~isequaln(obj.vector_param_config, other.vector_param_config)
                reason = 'vector_param_config differ'; return;
            end

            % Model parameters: prefer the resolved set frozen at run time.
            has_this = ~isempty(fieldnames(obj.resolved_defaults));
            has_other = ~isempty(fieldnames(other.resolved_defaults));
            if has_this && has_other
                [ok, reason] = ParamSpaceAnalysis2.compare_param_structs( ...
                    obj.resolved_defaults, other.resolved_defaults, 'resolved_defaults');
                if ~ok; return; end
            elseif allow_legacy
                [ok, reason] = ParamSpaceAnalysis2.compare_param_structs( ...
                    obj.model_defaults, other.model_defaults, 'model_defaults');
                if ~ok; return; end
            else
                if ~has_this && ~has_other
                    which_side = 'both runs predate';
                elseif ~has_this
                    which_side = 'this run predates';
                else
                    which_side = 'the other run predates';
                end
                reason = sprintf(['%s resolved_defaults, so model parameters cannot ' ...
                    'be compared exactly; pass ''allow_legacy'', true to fall back ' ...
                    'to comparing model_defaults'], which_side);
                return;
            end

            tf = true;
        end

        function plot(obj, varargin)
            % PLOT Generate histogram plots of metrics across parameter space
            %
            % Usage:
            %   psa.plot()
            %   psa.plot('metric', 'LLE')
            %   psa.plot('metric', 'mean_rate')

            if ~obj.has_run && isempty(fieldnames(obj.results))
                error('ParamSpaceAnalysis2:NotRun', ...
                    'Analysis has not been run yet. Call run() first.');
            end

            % Parse arguments
            metric = 'LLE';
            pool_with = {};   % cell of other PSA objects to pool into the distribution
            for i = 1:2:length(varargin)
                if strcmpi(varargin{i}, 'metric')
                    metric = varargin{i+1};
                elseif strcmpi(varargin{i}, 'pool_with')
                    pool_with = varargin{i+1};
                end
            end
            if ~iscell(pool_with); pool_with = {pool_with}; end

            % Define histogram bins
            if strcmpi(metric, 'LLE')
                hist_range = [-1.5, 1.5];
                n_bins = 25;
                y_label = 'LLE (\lambda_1)';
            elseif strcmpi(metric, 'mean_rate')
                hist_range = [0, 1];  % Activation function caps at 1
                n_bins = 25;
                y_label = 'Mean Firing Rate';
            else
                hist_range = [-10, 10];
                n_bins = 25;
                y_label = metric;
            end

            hist_bins = [linspace(hist_range(1), hist_range(2), n_bins + 1), inf];
            if strcmpi(metric, 'LLE')
                hist_bins = [-inf, hist_bins];
            end

            % Get condition names
            condition_names = cellfun(@(c) c.name, obj.conditions, 'UniformOutput', false);
            num_conditions = length(condition_names);

            % Readable titles
            condition_titles = srnn_condition_titles();   % one source; see that file

            % Create figure
            fig = figure('Name', sprintf('%s Distribution', metric), ...
                'Position', [100, 100, 300 * num_conditions, 300]);

            for c_idx = 1:num_conditions
                cond_name = condition_names{c_idx};
                ax = subplot(1, num_conditions, c_idx);

                % Extract metric values across the grid, pooled across runs.
                values = ParamSpaceAnalysis2.collect_grid_values(obj.results, cond_name, metric);
                for pp = 1:numel(pool_with)
                    values = [values, ...
                        ParamSpaceAnalysis2.collect_grid_values(pool_with{pp}.results, cond_name, metric)]; %#ok<AGROW>
                end

                % Plot histogram
                if ~isempty(values)
                    [counts, edges] = histcounts(values, hist_bins);
                    % Normalize to probability
                    prob = counts / sum(counts);

                    % Create finite edges for plotting
                    finite_edges = edges;
                    step = (hist_range(2) - hist_range(1)) / n_bins;
                    if isinf(finite_edges(1))
                        finite_edges(1) = hist_range(1) - step;
                    end
                    if isinf(finite_edges(end))
                        finite_edges(end) = hist_range(2) + step;
                    end

                    histogram('BinEdges', finite_edges, 'BinCounts', prob, ...
                        'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);

                    hold on;
                    if strcmpi(metric, 'LLE')
                        xline(0, '--', 'Color', [0 0.7 0], 'LineWidth', 2);
                    end
                    hold off;
                end

                % Labels
                if condition_titles.isKey(cond_name)
                    title(condition_titles(cond_name), 'FontSize', 14);
                else
                    title(strrep(cond_name, '_', ' '), 'FontSize', 14);
                end

                if c_idx == 1
                    ylabel('Probability', 'FontSize', 12);
                end
                xlabel(y_label, 'FontSize', 12);
                box off;
            end

            % Link y-axes
            drawnow;  % Force figure to render before linking
            ax_handles = findobj(fig, 'Type', 'Axes');
            linkaxes(ax_handles, 'y');

            % Figure saving disabled here -- handled by calling scripts via save_some_figs_to_folder_2
            % fig_dir = fullfile(obj.output_dir, 'figures');
            % if ~exist(fig_dir, 'dir')
            %     mkdir(fig_dir);
            % end
            % saveas(fig, fullfile(fig_dir, sprintf('%s_distribution.png', metric)));
            % saveas(fig, fullfile(fig_dir, sprintf('%s_distribution.fig', metric)));
            % fprintf('Figure saved to: %s\n', fig_dir);
        end

        function plot_sensitivity(obj, varargin)
            % PLOT_SENSITIVITY Generate imagesc heatmap plots for 1D sensitivity analysis
            %
            % Creates one subplot per condition showing histogram of metric values
            % across levels of each swept parameter (non-reps grid params).
            %
            % Usage:
            %   psa.plot_sensitivity()
            %   psa.plot_sensitivity('metric', 'LLE')
            %   psa.plot_sensitivity('metric', 'LLE', 'hist_range', [-0.3, 0.1], 'n_bins', 35)

            if ~obj.has_run && isempty(fieldnames(obj.results))
                error('ParamSpaceAnalysis2:NotRun', ...
                    'Analysis has not been run yet. Call run() first.');
            end

            % Parse arguments
            metric = 'LLE';
            hist_range = [];
            n_bins = 35;
            pool_with = {};   % cell of other PSA objects to pool (more reps)
            for i = 1:2:length(varargin)
                switch lower(varargin{i})
                    case 'metric', metric = varargin{i+1};
                    case 'hist_range', hist_range = varargin{i+1};
                    case 'n_bins', n_bins = varargin{i+1};
                    case 'pool_with', pool_with = varargin{i+1};
                end
            end
            if ~iscell(pool_with); pool_with = {pool_with}; end

            % Default hist_range by metric
            if isempty(hist_range)
                if strcmpi(metric, 'LLE')
                    hist_range = [-0.3, 0.1];
                elseif strcmpi(metric, 'mean_rate')
                    hist_range = [0, 1];
                else
                    hist_range = [-1, 1];
                end
            end

            % Identify swept parameters (non-reps grid params)
            swept_params = setdiff(obj.grid_params, {'reps'}, 'stable');
            if isempty(swept_params)
                error('ParamSpaceAnalysis2:NoSweptParam', ...
                    'No non-reps grid parameter found for sensitivity plot.');
            end

            % Build histogram bins
            hist_bins = [-inf, linspace(hist_range(1), hist_range(2), n_bins), inf];

            % Get condition info
            condition_names = cellfun(@(c) c.name, obj.conditions, 'UniformOutput', false);
            num_conditions = length(condition_names);

            % Readable condition titles
            condition_titles = srnn_condition_titles();   % one source; see that file

            % Loop over each swept parameter
            for sp_idx = 1:length(swept_params)
                swept_param = swept_params{sp_idx};
                param_dim = find(strcmp(obj.grid_params, swept_param));
                reps_dim = find(strcmp(obj.grid_params, 'reps'));

                % Get x-axis values
                if isfield(obj.vector_param_config, swept_param)
                    vpc = obj.vector_param_config.(swept_param);
                    vecs = obj.vector_param_lookup.(swept_param);
                    x_values = zeros(1, length(vecs));
                    for v = 1:length(vecs)
                        if strcmp(vpc.vary_element, 'last')
                            x_values(v) = vecs{v}(end);
                        else
                            x_values(v) = vecs{v}(1);
                        end
                    end
                    x_label = sprintf('%s(%s)', strrep(swept_param, '_', '\_'), vpc.vary_element);
                else
                    x_values = obj.param_vectors{param_dim};
                    x_label = strrep(swept_param, '_', '\_');
                end

                n_levels_param = length(x_values);
                n_reps = length(obj.param_vectors{reps_dim});
                % Total reps across obj + pooled runs (sets the colour scale).
                total_reps = n_reps;
                for pp = 1:numel(pool_with)
                    rdim_pp = find(strcmp(pool_with{pp}.grid_params, 'reps'));
                    total_reps = total_reps + length(pool_with{pp}.param_vectors{rdim_pp});
                end

                % Create figure
                fig = figure('Name', sprintf('%s Sensitivity - %s', metric, swept_param), ...
                    'Position', [100, 100, 350 * num_conditions, 400]);

                for c_idx = 1:num_conditions
                    cond_name = condition_names{c_idx};
                    ax = subplot(1, num_conditions, c_idx);

                    num_hist_bins = length(hist_bins) - 1;
                    histogram_matrix = zeros(num_hist_bins, n_levels_param);
                    median_values = NaN(n_levels_param, 1);

                    for level_idx = 1:n_levels_param
                        % Rep values for this level, pooled across obj + pool_with.
                        values_level = ParamSpaceAnalysis2.collect_level_values( ...
                            obj, swept_param, level_idx, cond_name, metric);
                        for pp = 1:numel(pool_with)
                            values_level = [values_level, ...
                                ParamSpaceAnalysis2.collect_level_values( ...
                                pool_with{pp}, swept_param, level_idx, cond_name, metric)]; %#ok<AGROW>
                        end

                        if ~isempty(values_level)
                            [counts, ~] = histcounts(values_level, hist_bins);
                            histogram_matrix(:, level_idx) = counts;
                            median_values(level_idx) = median(values_level);
                        end
                    end

                    % Compute y-coordinates
                    finite_edges = hist_bins(~isinf(hist_bins));
                    step_size = finite_edges(2) - finite_edges(1);
                    y_coords = zeros(num_hist_bins, 1);
                    y_coords(1) = finite_edges(1) - step_size/2;
                    for k = 2:length(finite_edges)
                        y_coords(k) = (finite_edges(k-1) + finite_edges(k)) / 2;
                    end
                    y_coords(end) = finite_edges(end) + step_size/2;

                    % Plot
                    imagesc(ax, x_values, y_coords, histogram_matrix);
                    hold(ax, 'on');
                    yline(ax, 0, '--', 'Color', [0 0.7 0], 'LineWidth', 4, 'Alpha', 0.5);
                    plot(ax, x_values, median_values, 'b-', 'LineWidth', 4, 'Color', [0 0 1 0.55]);
                    hold(ax, 'off');

                    colormap(ax, flipud(gray));
                    caxis(ax, [0, total_reps]);
                    axis(ax, 'xy');
                    box(ax, 'on');

                    % Labels
                    xlabel(ax, x_label, 'FontSize', 14);
                    if c_idx == 1
                        if strcmpi(metric, 'LLE')
                            ylabel(ax, '$\lambda_1$', 'Interpreter', 'latex', 'FontSize', 18);
                        else
                            ylabel(ax, strrep(metric, '_', '\_'), 'FontSize', 14);
                        end
                    end

                    if condition_titles.isKey(cond_name)
                        title(ax, condition_titles(cond_name), 'FontSize', 14);
                    else
                        title(ax, strrep(cond_name, '_', ' '), 'FontSize', 14);
                    end
                end

                % Figure saving disabled here -- handled by calling scripts via save_some_figs_to_folder_2
                % fig_dir = fullfile(obj.output_dir, 'figures');
                % if ~exist(fig_dir, 'dir')
                %     mkdir(fig_dir);
                % end
                % saveas(fig, fullfile(fig_dir, sprintf('sensitivity_%s_%s.png', metric, swept_param)));
                % saveas(fig, fullfile(fig_dir, sprintf('sensitivity_%s_%s.fig', metric, swept_param)));
                % fprintf('Figure saved to: %s\n', fig_dir);
            end
        end

        function [fig_handles, f_colorbar_fig] = plot_unit_histograms(obj, varargin)
            % PLOT_UNIT_HISTOGRAMS Create unit histograms colored by f value
            %
            % Creates one figure per metric, each with one subplot per condition.
            % Points are colored by fraction excitatory (f value).
            %
            % Usage:
            %   psa.plot_unit_histograms()
            %   psa.plot_unit_histograms('Metrics', {'lle', 'br'})
            %   [figs, cb_fig] = psa.plot_unit_histograms('NormalizeMode', 'probability')
            %
            % Name-value args:
            %   Metrics       - Cell array: 'lle', 'r', 'br' (default {'lle','r'})
            %   NormalizeMode - 'count' (default) or 'probability'
            %   color_by      - SRNNModel2 parameter to colour by (default 'f')

            if ~obj.has_run && isempty(fieldnames(obj.results))
                error('ParamSpaceAnalysis2:NotRun', ...
                    'Analysis has not been run yet. Call run() first.');
            end

            % Parse name-value args
            normalize_mode = 'count';
            metrics_to_plot = {'lle', 'r'};
            color_by = 'f';
            for i = 1:2:length(varargin)
                switch lower(varargin{i})
                    case 'normalizemode', normalize_mode = varargin{i+1};
                    case 'metrics', metrics_to_plot = lower(varargin{i+1});
                    case 'color_by', color_by = varargin{i+1};
                end
            end

            % Condition info
            condition_names = cellfun(@(c) c.name, obj.conditions, 'UniformOutput', false);
            num_conditions = length(condition_names);

            condition_titles = srnn_condition_titles();   % one source; see that file

            % Metric configuration
            metric_config = struct();
            metric_config.lle = struct('field', 'LLE', 'label', '\lambda_1', 'range', [-2.3, 1.5], 'inf_both', true);
            metric_config.r = struct('field', 'mean_rate', 'label', 'Mean Firing Rate', 'range', [0, 1], 'inf_both', false);
            metric_config.br = struct('field', 'mean_synaptic_output', 'label', 'Mean Synaptic Output', 'range', [0, 1], 'inf_both', false);

            % Filter to requested metrics
            metrics = {};
            metric_labels = {};
            metric_ranges = {};
            metric_inf_both = {};
            for i = 1:length(metrics_to_plot)
                key = metrics_to_plot{i};
                if isfield(metric_config, key)
                    cfg = metric_config.(key);
                    metrics{end+1} = cfg.field; %#ok<AGROW>
                    metric_labels{end+1} = cfg.label; %#ok<AGROW>
                    metric_ranges{end+1} = cfg.range; %#ok<AGROW>
                    metric_inf_both{end+1} = cfg.inf_both; %#ok<AGROW>
                else
                    warning('ParamSpaceAnalysis2:UnknownMetric', ...
                        'Unknown metric: %s. Valid: lle, r, br', key);
                end
            end

            if isempty(metrics)
                error('ParamSpaceAnalysis2:NoValidMetrics', 'No valid metrics specified.');
            end

            n_bins = 25;

            % Precompute bin edges
            metric_edges = cell(1, length(metrics));
            for m_idx = 1:length(metrics)
                hist_range = metric_ranges{m_idx};
                edges = linspace(hist_range(1), hist_range(2), n_bins + 1);
                if metric_inf_both{m_idx}
                    edges = [-inf, edges, inf];
                else
                    edges = [edges, inf];
                end
                metric_edges{m_idx} = edges;
            end

            % Collect all colour-parameter values for global normalization
            all_f_combined = [];
            for c_idx = 1:num_conditions
                cond_name = condition_names{c_idx};
                if isfield(obj.results, cond_name)
                    results_cell = obj.results.(cond_name);
                    for k = 1:length(results_cell)
                        res = results_cell{k};
                        if isstruct(res) && isfield(res, 'success') && res.success
                            all_f_combined(end+1) = obj.effective_param(res, color_by); %#ok<AGROW>
                        end
                    end
                end
            end

            % A constant colour parameter would give a degenerate CLim, so fall
            % back to the unit range as if no values had been found at all.
            if ~isempty(all_f_combined) && max(all_f_combined) > min(all_f_combined)
                f_min = min(all_f_combined);
                f_max = max(all_f_combined);
                has_f_variation = true;
            else
                f_min = 0; f_max = 1; has_f_variation = false;
            end

            if has_f_variation
                fprintf('Coloring by %s value: [%.3f, %.3f]\n', color_by, f_min, f_max);
            end

            cmap_f = blue_gray_red_colormap(256);

            % Create figures for each metric
            fig_handles = gobjects(1, length(metrics));

            for m_idx = 1:length(metrics)
                metric = metrics{m_idx};
                metric_label = metric_labels{m_idx};
                metric_range = metric_ranges{m_idx};

                fig = figure('Name', sprintf('%s Unit Histogram', metric), ...
                    'Position', [100, 100, 300 * num_conditions, 300]);

                for c_idx = 1:num_conditions
                    cond_name = condition_names{c_idx};
                    ax = subplot(1, num_conditions, c_idx);

                    values = [];
                    f_values = [];

                    if isfield(obj.results, cond_name)
                        results_cell = obj.results.(cond_name);
                        for k = 1:length(results_cell)
                            res = results_cell{k};
                            if isstruct(res) && isfield(res, 'success') && res.success
                                if isfield(res, metric) && ~isnan(res.(metric))
                                    values(end+1) = res.(metric); %#ok<AGROW>
                                    f_values(end+1) = obj.effective_param(res, color_by); %#ok<AGROW>
                                end
                            end
                        end
                    end

                    if ~isempty(values)
                        unit_histogram_patch(values(:), f_values(:), ...
                            'BinEdges', metric_edges{m_idx}, ...
                            'SortMode', 'sorted', ...
                            'Axes', ax, ...
                            'Colormap', cmap_f, ...
                            'CLim', [f_min, f_max], ...
                            'Normalize', normalize_mode, ...
                            'EdgeColor', 'none');

                        if strcmpi(metric, 'LLE')
                            hold(ax, 'on');
                            xline(ax, 0, '--', 'Color', [0 0 0], 'LineWidth', 2);
                            hold(ax, 'off');
                        end
                    end

                    if condition_titles.isKey(cond_name)
                        title(ax, condition_titles(cond_name), 'FontWeight', 'normal');
                    else
                        title(ax, strrep(cond_name, '_', ' '), 'FontWeight', 'normal');
                    end

                    if c_idx == 1
                        if strcmpi(normalize_mode, 'probability')
                            ylabel(ax, 'Probability');
                        else
                            ylabel(ax, 'Count');
                        end
                    end
                    xlabel(ax, metric_label);
                    xlim(ax, 1.05 .* metric_range);
                end

                drawnow;
                ax_handles = findobj(fig, 'Type', 'Axes');
                linkaxes(ax_handles, 'y');
                fig_handles(m_idx) = fig;
            end

            % Create colorbar figure for f values
            f_colorbar_fig = gobjects(0);
            if has_f_variation
                f_colorbar_fig = figure('Name', 'f Value Colorbar', ...
                    'Position', [500, 200, 300, 300]);
                ax_cb = axes(f_colorbar_fig);
                gradient_img = repmat(linspace(0, 1, 256)', 1, 2);
                imagesc(ax_cb, [0 1], [f_min f_max], gradient_img);
                colormap(ax_cb, cmap_f);
                ax_cb.XTick = [];
                ax_cb.YDir = 'normal';
                ax_cb.XColor = 'none';
                ylabel(ax_cb, 'fraction excitatory', 'FontSize', 12);
                box(ax_cb, 'off');
                pbaspect(ax_cb, [0.1 1 1]);
            end

            fprintf('Unit histograms generated.\n');
        end

        function [fig, p_values] = plot_lle_by_stim_period(obj, varargin)
            % PLOT_LLE_BY_STIM_PERIOD Paired beeswarm of mean local LLE per stimulus step
            %
            % Creates a figure with one subplot per condition showing mean local
            % LLE during each stimulus step, colored by f value, with significance
            % brackets (Wilcoxon signed-rank test).
            %
            % Requires that the analysis was run with store_local_lya = true.
            %
            % Usage:
            %   psa.plot_lle_by_stim_period()
            %   psa.plot_lle_by_stim_period('transient_skip', 3)
            %   [fig, p] = psa.plot_lle_by_stim_period('periods_to_plot', [0 1 1])
            %
            % Name-value args:
            %   transient_skip  - Seconds to skip at the start of each step
            %   periods_to_plot - Logical mask over the steps to show
            %   color_by        - SRNNModel2 parameter to colour by (default 'f')

            if ~obj.has_run && isempty(fieldnames(obj.results))
                error('ParamSpaceAnalysis2:NotRun', ...
                    'Analysis has not been run yet. Call run() first.');
            end

            % Parse args
            transient_skip = 0;
            periods_to_plot_arg = [];
            color_by = 'f';
            for i = 1:2:length(varargin)
                switch lower(varargin{i})
                    case 'transient_skip', transient_skip = varargin{i+1};
                    case 'periods_to_plot', periods_to_plot_arg = logical(varargin{i+1});
                    case 'color_by', color_by = varargin{i+1};
                end
            end

            % Validate local_lya exists
            cond_names = cellfun(@(c) c.name, obj.conditions, 'UniformOutput', false);
            has_local_lya = false;
            for c_idx = 1:length(cond_names)
                cond_name = cond_names{c_idx};
                if isfield(obj.results, cond_name)
                    for k = 1:length(obj.results.(cond_name))
                        res = obj.results.(cond_name){k};
                        if isstruct(res) && isfield(res, 'local_lya') && ~isempty(res.local_lya)
                            has_local_lya = true;
                            break;
                        end
                    end
                end
                if has_local_lya, break; end
            end

            if ~has_local_lya
                error('ParamSpaceAnalysis2:NoLocalLya', ...
                    'Results do not contain local_lya. Re-run with store_local_lya = true.');
            end

            % Step configuration (run-wide, so no per-result context is needed)
            stim_config = obj.effective_param([], 'input_config');
            n_steps = stim_config.n_steps;
            no_stim_pattern = stim_config.no_stim_pattern;

            T_stim_range = obj.effective_param([], 'T_range');
            T_stim = T_stim_range(2);
            step_period = T_stim / n_steps;

            num_conditions = length(cond_names);

            % Compute mean LLE per step for each simulation
            all_means = cell(num_conditions, 1);
            all_f_values = cell(num_conditions, 1);

            for c_idx = 1:num_conditions
                cond_name = cond_names{c_idx};
                results_cell = obj.results.(cond_name);

                n_valid = 0;
                for k = 1:length(results_cell)
                    res = results_cell{k};
                    if isstruct(res) && isfield(res, 'success') && res.success && ...
                            isfield(res, 'local_lya') && ~isempty(res.local_lya)
                        n_valid = n_valid + 1;
                    end
                end

                means_matrix = NaN(n_valid, n_steps);
                f_values = NaN(n_valid, 1);
                sim_idx = 0;

                for k = 1:length(results_cell)
                    res = results_cell{k};
                    if ~isstruct(res) || ~isfield(res, 'success') || ~res.success
                        continue;
                    end
                    if ~isfield(res, 'local_lya') || isempty(res.local_lya)
                        continue;
                    end

                    sim_idx = sim_idx + 1;
                    t_lya = res.t_lya;
                    local_lya = res.local_lya;

                    f_values(sim_idx) = obj.effective_param(res, color_by);

                    valid_mask = t_lya >= 0;
                    t_lya = t_lya(valid_mask);
                    local_lya = local_lya(valid_mask);

                    for step_idx = 1:n_steps
                        step_start = (step_idx - 1) * step_period + transient_skip;
                        step_end = step_idx * step_period;
                        step_mask = t_lya >= step_start & t_lya < step_end;
                        if any(step_mask)
                            means_matrix(sim_idx, step_idx) = mean(local_lya(step_mask), 'omitnan');
                        end
                    end
                end

                all_means{c_idx} = means_matrix;
                all_f_values{c_idx} = f_values;
            end

            % Global f-value range
            all_f_combined = [];
            for c_idx = 1:num_conditions
                all_f_combined = [all_f_combined; all_f_values{c_idx}]; %#ok<AGROW>
            end
            all_f_combined = all_f_combined(~isnan(all_f_combined));

            % A constant colour parameter would give a degenerate CLim, so fall
            % back to the unit range as if no values had been found at all.
            if ~isempty(all_f_combined) && max(all_f_combined) > min(all_f_combined)
                f_min = min(all_f_combined);
                f_max = max(all_f_combined);
                has_f_variation = true;
            else
                f_min = 0; f_max = 1; has_f_variation = false;
            end

            % Create figure
            condition_titles = srnn_condition_titles();   % one source; see that file

            fig = figure('Name', 'Mean Local LLE by Step', ...
                'Position', [100, 100, 300 * num_conditions, 300]);

            cmap_f = blue_gray_red_colormap(256);
            n_colors = size(cmap_f, 1);

            p_values = NaN(num_conditions, 1);
            ax_cell = cell(num_conditions, 1);

            for c_idx = 1:num_conditions
                cond_name = cond_names{c_idx};
                means_matrix = all_means{c_idx};
                f_vals = all_f_values{c_idx};

                ax = subplot(1, num_conditions, c_idx);
                n_sims = size(means_matrix, 1);

                % Build color matrix
                sim_colors = zeros(n_sims, 3);
                for s = 1:n_sims
                    if has_f_variation && ~isnan(f_vals(s))
                        f_normalized = max(0, min(1, (f_vals(s) - f_min) / (f_max - f_min)));
                        color_idx = round(f_normalized * (n_colors - 1)) + 1;
                        sim_colors(s, :) = cmap_f(color_idx, :);
                    else
                        sim_colors(s, :) = [0.5 0.5 0.5];
                    end
                end

                % Step labels
                step_labels = cell(1, n_steps);
                for step_idx = 1:n_steps
                    if no_stim_pattern(step_idx)
                        step_labels{step_idx} = 'no-stim';
                    else
                        step_labels{step_idx} = 'stim';
                    end
                end

                % Apply periods_to_plot filter
                if isempty(periods_to_plot_arg)
                    periods_mask = true(1, n_steps);
                else
                    periods_mask = periods_to_plot_arg;
                end

                means_filtered = means_matrix(:, periods_mask);
                labels_filtered = step_labels(periods_mask);
                n_steps_plot = sum(periods_mask);

                % Stability reference line
                line(ax, [0.3, n_steps_plot + 0.7], [0, 0], ...
                    'Color', [0 0 0 0.5], 'LineStyle', '--', 'LineWidth', 1, ...
                    'HandleVisibility', 'off');
                hold(ax, 'on');

                paired_beeswarm(means_filtered, 'Axes', ax, ...
                    'Colors', sim_colors, ...
                    'MarkerSize', 0.9, ...
                    'LineWidth', 0.75, ...
                    'Labels', labels_filtered, ...
                    'Alpha', 1, ...
                    'SortStyle', 'hex', ...
                    'ShowYAxis', (c_idx == 1));

                if condition_titles.isKey(cond_name)
                    title(condition_titles(cond_name), 'FontWeight', 'normal');
                else
                    title(strrep(cond_name, '_', ' '), 'FontWeight', 'normal');
                end

                if c_idx == 1
                    ylabel('Mean Local LLE');
                else
                    ax.YAxis.Visible = 'off';
                end

                hold(ax, 'off');
                xtickangle(ax, 45);
                box(ax, 'off');

                % Significance test
                stim_means = mean(means_matrix(:, ~no_stim_pattern), 2, 'omitnan');
                no_stim_means_col = mean(means_matrix(:, no_stim_pattern), 2, 'omitnan');
                valid_pairs = ~isnan(stim_means) & ~isnan(no_stim_means_col);

                if sum(valid_pairs) >= 2
                    [p_val, ~, ~] = signrank(stim_means(valid_pairs), no_stim_means_col(valid_pairs));
                    p_values(c_idx) = p_val;
                end

                ax_cell{c_idx} = ax;
            end

            % Link y-axes
            drawnow;
            ax_handles = findobj(fig, 'Type', 'Axes');
            if length(ax_handles) >= 2
                linkaxes(ax_handles, 'y');
            end
            for i = 1:length(ax_handles)
                xlim(ax_handles(i), [0.5, n_steps_plot + 0.5]);
            end

            % Compute global y-limits
            global_ymin = inf;
            global_ymax = -inf;
            for i = 1:length(ax_handles)
                children = ax_handles(i).Children;
                for j = 1:length(children)
                    if isprop(children(j), 'YData')
                        ydata = children(j).YData;
                        ydata = ydata(isfinite(ydata));
                        if ~isempty(ydata)
                            global_ymin = min(global_ymin, min(ydata));
                            global_ymax = max(global_ymax, max(ydata));
                        end
                    end
                end
            end

            global_range = global_ymax - global_ymin;
            if global_range <= 0, global_range = 0.2; end
            global_padding = 0.05 * global_range;
            global_upper = max(global_ymax + global_padding, 0.05);
            global_lower = global_ymin - global_padding;
            global_upper_with_bracket = global_upper + 0.15 * global_range;

            if ~isempty(ax_handles)
                ylim(ax_handles(1), [global_lower, global_upper_with_bracket]);
            end

            % Add significance brackets
            bracket_y = global_upper + 0.03 * global_range;
            vert_height = 0.03 * global_range;
            text_y = bracket_y + 0.05 * global_range;

            for c_idx = 1:num_conditions
                ax = ax_cell{c_idx};
                p_val = p_values(c_idx);
                if isnan(p_val), continue; end

                hold(ax, 'on');
                x_left = 1; x_right = 2;
                bracket_x = [x_left, x_left, x_right, x_right];
                bracket_y_pts = [bracket_y - vert_height, bracket_y, bracket_y, bracket_y - vert_height];
                plot(ax, bracket_x, bracket_y_pts, 'k-', 'LineWidth', 2.5, 'HandleVisibility', 'off');

                p_str = sprintf('p = %.2g', p_val);
                text(ax, (x_left + x_right) / 2, text_y, p_str, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                    'FontSize', 9, 'HandleVisibility', 'off');
                hold(ax, 'off');
            end

            drawnow;

            % Create colorbar figure for f values
            if has_f_variation
                fig_cb = figure('Name', 'f Value Colorbar', ...
                    'Position', [500, 200, 300, 300]);
                ax_cb = axes(fig_cb);
                gradient_img = repmat(linspace(0, 1, 256)', 1, 2);
                imagesc(ax_cb, [0 1], [f_min f_max], gradient_img);
                colormap(ax_cb, cmap_f);
                ax_cb.XTick = [];
                ax_cb.YDir = 'normal';
                ax_cb.XColor = 'none';
                ylabel(ax_cb, 'fraction excitatory', 'FontSize', 12);
                box(ax_cb, 'off');
                pbaspect(ax_cb, [0.1 1 1]);
            end

            % Print stats
            fprintf('\n=== Statistical Analysis: Stim vs No-Stim ===\n');
            for c_idx = 1:num_conditions
                cond_name = cond_names{c_idx};
                means_matrix = all_means{c_idx};

                stim_means = mean(means_matrix(:, ~no_stim_pattern), 2, 'omitnan');
                no_stim_means_col = mean(means_matrix(:, no_stim_pattern), 2, 'omitnan');
                valid_pairs = ~isnan(stim_means) & ~isnan(no_stim_means_col);
                n_pairs = sum(valid_pairs);

                if condition_titles.isKey(cond_name)
                    display_name = condition_titles(cond_name);
                else
                    display_name = strrep(cond_name, '_', ' ');
                end

                if n_pairs >= 2
                    stim_valid = stim_means(valid_pairs);
                    no_stim_valid = no_stim_means_col(valid_pairs);
                    [p_value, ~, ~] = signrank(stim_valid, no_stim_valid);
                    differences = stim_valid - no_stim_valid;
                    cohens_d = mean(differences) / std(differences);

                    fprintf('%s (n=%d pairs):\n', display_name, n_pairs);
                    fprintf('  p-value: %.4g, Cohen''s d: %.4f\n', p_value, cohens_d);
                    fprintf('  Median diff (stim-nostim): %.4f\n', ...
                        median(stim_valid) - median(no_stim_valid));
                else
                    fprintf('%s: Insufficient pairs (n=%d)\n', display_name, n_pairs);
                end
            end
        end

        function load_condition_results(obj)
            % LOAD_CONDITION_RESULTS Read the per-condition result files into
            % obj.results. Configuration must already be set -- this loads
            % RESULTS ONLY.
            %
            % This was `load_results`, which also restored configuration from
            % param_space_summary.mat using its own hand-picked list of six
            % fields. That list had drifted from the twelve the summary writes
            % and the twenty-six saveobj writes, so a loaded object silently
            % lost model_class and, worse, vector_param_lookup -- leaving
            % effective_param to hand back grid LEVEL INDICES for vector
            % parameters instead of values. Both were fixed on 2026-08-14 and
            % are covered by test_psa_loaders.
            %
            % The configuration now comes from psa_object.mat via from_dir, and
            % param_space_summary.mat is metadata only. Keeping this method to
            % results alone is what stops a fourth restore list existing.

            for c_idx = 1:length(obj.conditions)
                cond_name = obj.conditions{c_idx}.name;
                results_file = fullfile(obj.output_dir, cond_name, ...
                    sprintf('param_space_results_%s.mat', cond_name));
                if exist(results_file, 'file')
                    loaded = load(results_file, 'results');
                    obj.results.(cond_name) = loaded.results;
                    if obj.verbose
                        fprintf('Loaded %d results for condition: %s\n', ...
                            length(loaded.results), cond_name);
                    end
                end
            end

            obj.has_run = true;
        end

        function consolidate(obj, temp_dir)
            % CONSOLIDATE Merge batch results into per-condition MAT files
            %
            % Usage:
            %   psa.consolidate()           % Standalone recovery after interrupted run
            %   psa.consolidate(temp_dir)   % Internal call from run() with explicit path
            %
            % When called without arguments:
            %   - Validates that output_dir is set and temp_batches/ exists
            %   - Recovers configuration from summary file or batch files
            %   - Saves summary and sets has_run = true after completion
            %
            % When called with temp_dir (internal use by run()):
            %   - Assumes configuration is already set
            %   - Does not save summary (run() handles that)

            % Determine if this is a standalone call or internal call from run()
            standalone_call = (nargin < 2 || isempty(temp_dir));

            if standalone_call
                % Standalone recovery mode - do validation and config recovery
                if isempty(obj.output_dir)
                    error('ParamSpaceAnalysis2:NoOutputDir', ...
                        'output_dir is not set. Set it to the run directory first.');
                end

                temp_dir = fullfile(obj.output_dir, 'temp_batches');

                if ~exist(temp_dir, 'dir')
                    error('ParamSpaceAnalysis2:NoTempDir', ...
                        'No temp_batches directory found in %s.\nRun may have already been consolidated, or output_dir is incorrect.', ...
                        obj.output_dir);
                end

                % Configuration comes from the object, not from a summary file.
                % run() writes psa_object.mat BEFORE batching precisely so that an
                % interrupted run still has its configuration on disk; this method
                % used to re-derive it from param_space_summary.mat with a third
                % hand-maintained field list (and, failing that, to GUESS the
                % conditions and combination count from the batch files). Both are
                % gone: from_dir loads the object and calls this, so obj is already
                % configured.
                if isempty(obj.conditions) || isempty(obj.num_combinations)
                    error('ParamSpaceAnalysis2:NotConfigured', ...
                        ['consolidate() needs a configured object (conditions and ' ...
                         'num_combinations).\nLoad the run with ' ...
                         'ParamSpaceAnalysis2.from_dir(''%s'') rather than ' ...
                         'consolidating a blank object.'], obj.output_dir);
                end

                fprintf('Consolidating results from %s...\n', temp_dir);
            end

            %% Core consolidation logic
            fprintf('\nConsolidating batch results...\n');

            % Batch count must match run_batched_simulation's, which batches
            % over shuffled_indices. That field is persisted by saveobj/loadobj,
            % so a consolidate() reached through from_dir gets the subset's
            % batch count with no extra saved state.
            n_run = numel(obj.shuffled_indices);
            num_batches = ceil(n_run / obj.batch_size);

            % Storage stays FULL-GRID length: results are indexed by config_idx,
            % so an unrun point must be an empty slot at its own coordinates,
            % not a missing one that shifts everything after it.
            for c_idx = 1:length(obj.conditions)
                cond_name = obj.conditions{c_idx}.name;
                obj.results.(cond_name) = cell(obj.num_combinations, 1);
            end

            % Load and merge batches
            all_found = true;
            for batch_idx = 1:num_batches
                batch_file = fullfile(temp_dir, sprintf('batch_%d.mat', batch_idx));

                if exist(batch_file, 'file')
                    loaded = load(batch_file);
                    batch_results = loaded.batch_results;

                    for c_idx = 1:length(obj.conditions)
                        cond_name = obj.conditions{c_idx}.name;
                        cond_results = batch_results.(cond_name);

                        for k = 1:length(cond_results)
                            res = cond_results{k};
                            if isstruct(res) && isfield(res, 'config_idx')
                                obj.results.(cond_name){res.config_idx} = res;
                            end
                        end
                    end
                else
                    fprintf('Warning: Batch file %d not found\n', batch_idx);
                    all_found = false;
                end
            end

            % Save per-condition results
            for c_idx = 1:length(obj.conditions)
                cond_name = obj.conditions{c_idx}.name;
                results = obj.results.(cond_name);

                save_file = fullfile(obj.output_dir, cond_name, ...
                    sprintf('param_space_results_%s.mat', cond_name));
                save(save_file, 'results', '-v7.3');

                % Count successes
                % Denominator is what was ATTEMPTED. Reporting against the full
                % grid would read a deliberate 15% subset as an 85% failure rate.
                n_success = sum(cellfun(@(r) isstruct(r) && isfield(r, 'success') && r.success, results));
                if n_run < obj.num_combinations
                    fprintf('Condition %s: %d/%d successful (subset of %d grid points), saved to %s\n', ...
                        cond_name, n_success, n_run, obj.num_combinations, save_file);
                else
                    fprintf('Condition %s: %d/%d successful, saved to %s\n', ...
                        cond_name, n_success, n_run, save_file);
                end
            end

            % Clean up temp directory if all successful
            if all_found
                rmdir(temp_dir, 's');
                fprintf('Temp directory cleaned up.\n');
            else
                fprintf('Temp directory retained due to missing batches.\n');
            end

            %% Finalize for standalone calls
            if standalone_call
                obj.save_summary();
                obj.has_run = true;
                fprintf('Consolidation complete. Results available in psa.results\n');
            end
        end

        function s = saveobj(obj)
            % SAVEOBJ Convert object to struct for saving with MATLAB's save()
            %
            % This method is called automatically by MATLAB's save() function.
            % Use load() to restore the object via the static loadobj() method.
            %
            % Usage:
            %   save('psa_saved.mat', 'psa');  % Automatically calls saveobj
            %   loaded = load('psa_saved.mat');
            %   psa_restored = loaded.psa;     % Automatically calls loadobj

            s = struct();

            % Configuration Properties (public)
            s.grid_params = obj.grid_params;
            s.param_ranges = obj.param_ranges;
            s.n_levels = obj.n_levels;
            s.conditions = obj.conditions;
            s.integer_params = obj.integer_params;
            s.explicit_vectors = obj.explicit_vectors;
            s.vector_param_config = obj.vector_param_config;
            s.randomize_order = obj.randomize_order;
            s.subset_fraction = obj.subset_fraction;
            s.network_seed_offset = obj.network_seed_offset;

            % Model Default Properties (public)
            s.model_class = obj.model_class;
            s.model_defaults = obj.model_defaults;
            s.verbose = obj.verbose;

            % Execution Properties (public)
            s.batch_size = obj.batch_size;
            s.output_dir = obj.output_dir;
            s.note = obj.note;
            s.folder_prefix = obj.folder_prefix;
            s.store_local_lya = obj.store_local_lya;
            s.store_local_lya_dt = obj.store_local_lya_dt;
            s.use_parallel = obj.use_parallel;

            % Results Properties (private)
            s.results = obj.results;
            s.has_run = obj.has_run;
            s.analysis_start_time = obj.analysis_start_time;
            s.param_vectors = obj.param_vectors;
            s.all_configs = obj.all_configs;
            s.shuffled_indices = obj.shuffled_indices;
            s.num_combinations = obj.num_combinations;
            s.vector_param_lookup = obj.vector_param_lookup;
            s.resolved_defaults = obj.resolved_defaults;
        end
    end

    %% Model-default introspection helpers
    methods (Static)
        function psa = from_dir(results_dir)
            % FROM_DIR Load a completed or interrupted run from its directory.
            %
            %   psa = ParamSpaceAnalysis2.from_dir('/path/to/param_space_...')
            %
            % The ONE way to read a run back off disk. psa_object.mat is the
            % authoritative artifact -- run() writes it before batching and again
            % on completion -- and this resolves the rest:
            %
            %   1. load psa_object.mat
            %   2. if it carries no results, load the per-condition result files
            %   3. if temp_batches/ is present, consolidate() first (a run that
            %      died part-way: the early object has the config, the batches
            %      have the results)
            %
            % Callers used to hand-roll that three-way choice, and each guessed
            % differently at the variable name inside psa_object.mat -- most
            % scripts saved `psa`, the tau script saved `psa_tau_a`/`psa_tau_b`.
            % Selection here is BY CLASS, not by name, so old files load without
            % anyone needing to know.
            %
            % Not named `load`: the class already defines loadobj, and `load` is a
            % core MATLAB function.
            %
            % See also: save_object, consolidate

            if nargin < 1 || isempty(results_dir) || ~ischar(results_dir) && ~isstring(results_dir)
                error('ParamSpaceAnalysis2:InvalidInput', ...
                    'from_dir requires a results directory path.');
            end
            results_dir = char(results_dir);
            if ~exist(results_dir, 'dir')
                error('ParamSpaceAnalysis2:NoSuchDir', ...
                    'Not a directory: %s', results_dir);
            end

            obj_file = fullfile(results_dir, 'psa_object.mat');
            if ~exist(obj_file, 'file')
                error('ParamSpaceAnalysis2:NoPsaObject', ...
                    ['No psa_object.mat in %s.\nSince run() writes it before ' ...
                     'batching, a directory without one predates that change or ' ...
                     'was not produced by run(). If temp_batches/ is present, ' ...
                     'set output_dir on a configured PSA and call consolidate().'], ...
                    results_dir);
            end

            % Select by CLASS rather than by variable name.
            S = load(obj_file);
            fns = fieldnames(S);
            is_psa = cellfun(@(f) isa(S.(f), 'ParamSpaceAnalysis2'), fns);
            if ~any(is_psa)
                error('ParamSpaceAnalysis2:BadPsaObject', ...
                    'No ParamSpaceAnalysis2 object found in %s (variables: %s).', ...
                    obj_file, strjoin(fns', ', '));
            end
            psa = S.(fns{find(is_psa, 1)});

            % The directory it was loaded FROM wins over the one recorded in it:
            % run directories get moved and copied between machines.
            psa.output_dir = results_dir;

            % An interrupted run: the early object has the configuration, the
            % batches have whatever finished.
            if exist(fullfile(results_dir, 'temp_batches'), 'dir')
                fprintf(['[from_dir] temp_batches/ present -- consolidating an ' ...
                    'interrupted run.\n']);
                psa.consolidate();
                return;
            end

            if isempty(fieldnames(psa.results))
                psa.load_condition_results();
            end
        end

        function val = class_default(name, class_name)
            % CLASS_DEFAULT Value a model class gives NAME with nothing overriding
            %
            % Backed by a cached throwaway model per class: the constructor runs
            % set_defaults (which fills input_config and computes plot_deci), so
            % no build() is needed. The cached object is never mutated. Using this
            % instead of a hand-copied literal is what keeps plotting fallbacks
            % from drifting when a class default changes.
            %
            % Classes whose constructor has required arguments cannot be default-
            % constructed; for those this errors, and effective_param falls back
            % before reaching here (it only consults the class default when the
            % name appears nowhere else).
            if nargin < 2, class_name = 'SRNNModel2'; end
            persistent cache
            if isempty(cache)
                cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
            end
            if ~cache.isKey(class_name) || ~isvalid(cache(class_name))
                cache(class_name) = feval(class_name);
            end
            m0 = cache(class_name);
            val = m0.(name);
        end
    end

    methods (Static, Access = private)
        function val = pooled_mean(s)
            % POOLED_MEAN Mean over every leaf array of a plot_data field,
            % ignoring NaNs. Works for any number of cell types and for either
            % nesting depth: SRNNModel2's br is keyed by type (.E, .I) while
            % SRNNCellTypePairs' synaptic_output is keyed by ROUTE
            % (.E.E, .E.I, ...), because its synapses are per-route.
            v = ParamSpaceAnalysis2.collect_leaves(s);
            v = v(~isnan(v));
            if isempty(v), val = NaN; else, val = mean(v); end
        end

        function v = collect_leaves(s)
            % COLLECT_LEAVES Every numeric element under a (possibly nested) struct.
            if ~isstruct(s)
                v = double(s(:));
                return;
            end
            f = fieldnames(s);
            parts = cell(numel(f), 1);
            for i = 1:numel(f)
                parts{i} = ParamSpaceAnalysis2.collect_leaves(s.(f{i}));
            end
            v = vertcat(parts{:});
        end

        function names = per_job_param_names(grid_params, condition_fields)
            % PER_JOB_PARAM_NAMES Parameters that vary per job, not per run
            %
            % The grid axes plus the condition fields THAT SOME CONDITION ACTUALLY
            % SETS. run_single_job supplies these from the job rather than from
            % model_defaults, so a run-wide default for any of them is ignored and
            % they are excluded from resolved_defaults.
            %
            % Passing the conditions' real field set rather than the static
            % {n_a_E, n_a_I, n_b_E, n_b_I} matters: the default conditions set only
            % n_a_E and n_b_E, so model_defaults.n_a_I DOES take effect and must
            % not be skipped or excluded.
            names = [grid_params(:)', condition_fields(:)'];
        end

        function reject_condition_field(param_name)
            % REJECT_CONDITION_FIELD Error if a name the conditions own is swept.
            if ismember(char(param_name), ParamSpaceAnalysis2.condition_field_names())
                error('ParamSpaceAnalysis2:ConditionFieldAsGridParam', ...
                    ['''%s'' cannot be a grid parameter: the adaptation conditions ' ...
                    'set it per job, and a grid value would silently override every ' ...
                    'condition. Vary it with set_conditions() instead.'], char(param_name));
            end
        end

        function names = condition_field_names()
            % CONDITION_FIELD_NAMES Names the conditions mechanism owns.
            % These may never be used as grid axes -- run_single_job appends grid
            % parameters AFTER condition parameters and the SRNNModel2 constructor
            % is last-write-wins, so a gridded n_a_E would silently override every
            % condition and collapse the adaptation regimes into one.
            names = {'n_a_E', 'n_a_I', 'n_b_E', 'n_b_I'};
        end

        function [tf, reason] = compare_param_structs(a, b, label)
            % COMPARE_PARAM_STRUCTS Field-by-field equality of two parameter structs
            %
            % Ignores properties that do not affect the dynamics (storage and
            % plotting options, per-job seeds), so pooling is not refused over a
            % different plot_freq. Function handles compare by their string form.
            % NOTE: func2str does not expose values captured by an anonymous
            % handle, so two handles differing only in a captured constant compare
            % equal -- acceptable here since pooled runs come from the same
            % script/config.
            ignore = {'rng_seeds', 'reps', 'store_full_state', ...
                'store_decimated_state', 'plot_deci', 'plot_freq', ...
                'T_plot', 'check_connectivity'};

            tf = false;
            fa = setdiff(fieldnames(a), ignore);
            fb = setdiff(fieldnames(b), ignore);
            if ~isequal(fa, fb)
                reason = sprintf('%s field sets differ', label);
                return;
            end
            for i = 1:numel(fa)
                k = fa{i};
                va = a.(k);
                vb = b.(k);
                if isa(va, 'function_handle') || isa(vb, 'function_handle')
                    if ~(isa(va, 'function_handle') && isa(vb, 'function_handle') ...
                            && strcmp(func2str(va), func2str(vb)))
                        reason = sprintf('%s.%s (function handle) differs', label, k);
                        return;
                    end
                elseif ~isequaln(va, vb)
                    reason = sprintf('%s.%s differs', label, k);
                    return;
                end
            end
            tf = true;
            reason = '';
        end

        function info = srnn_property_info(class_name)
            % SRNN_PROPERTY_INFO Classify a model class's properties for validation
            %
            % Returns a struct with cell-array fields:
            %   settable  - publicly settable (what model_defaults may contain)
            %   dependent - Dependent with no set method (computed, read-only)
            %   nonpublic - not publicly settable (protected/private)
            %
            % Cached in a persistent: the class definition does not change at
            % runtime and metaclass introspection is comparatively slow (this is
            % called once per result by effective_param, so plotting loops hit it
            % thousands of times).
            %
            % DEVELOPMENT NOTE: because it is cached, adding or renaming a model
            % property mid-session leaves this list stale, and the new name is
            % then reported as "not a property of <class>". Run `clear classes`
            % after editing a classdef -- `clear <classname>` and
            % `clear functions` do NOT clear this.
            persistent cache
            if isempty(cache)
                cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
            end
            if ~cache.isKey(class_name)
                mc = meta.class.fromName(class_name);
                if isempty(mc)
                    error('ParamSpaceAnalysis2:UnknownModelClass', ...
                        '''%s'' is not a class on the path.', class_name);
                end
                props = mc.PropertyList;
                names = {props.Name};
                is_dependent = [props.Dependent];
                has_setter = ~cellfun(@isempty, {props.SetMethod});
                % SetAccess is 'public'/'protected'/'private', or a cell array of
                % meta.class for SetAccess = ?SomeClass -- the ischar guard treats
                % the latter as non-public, which is what we want here.
                is_public = cellfun(@(a) ischar(a) && strcmp(a, 'public'), {props.SetAccess});

                info = struct();
                info.settable  = names(is_public & (~is_dependent | has_setter));
                info.dependent = names(is_dependent & ~has_setter);
                info.nonpublic = names(~is_public & ~is_dependent);
                cache(class_name) = info;
            end
            info = cache(class_name);
        end

        function s = suggest_property(name, candidates)
            % SUGGEST_PROPERTY Trailing " Did you mean ...?" text for a bad name.
            % Case-insensitive exact match first (catches tau_D -> tau_d), then
            % prefix/substring matches. Returns '' when nothing is close.
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
    end

    %% Static pooling helpers (used by plot / plot_sensitivity)
    % collect_level_values is PUBLIC: it is the one place that knows how to walk
    % the reps axis of a sweep and drop failed/NaN jobs, so a presentation script
    % that wants the same per-level rep values without the imagesc panel
    % (scripts/presentations/.../fig_sensitivity_medians) reuses it rather than
    % reimplementing the sub2ind indexing. collect_grid_values stays private.
    methods (Static)
        function vals = collect_level_values(psa, swept_param, level_idx, cond_name, metric)
            % Successful metric values across the reps axis for one swept-param
            % level and one condition. Reused to pool reps within & across runs.
            vals = [];
            if ~isfield(psa.results, cond_name); return; end
            rc = psa.results.(cond_name);
            gp = psa.grid_params;
            gsz = cellfun(@length, psa.param_vectors);
            param_dim = find(strcmp(gp, swept_param));
            reps_dim = find(strcmp(gp, 'reps'));
            n_reps = length(psa.param_vectors{reps_dim});
            for rep_idx = 1:n_reps
                subs = num2cell(ones(1, numel(gp)));
                subs{param_dim} = level_idx;
                subs{reps_dim} = rep_idx;
                lin = sub2ind(gsz, subs{:});
                if lin <= numel(rc)
                    res = rc{lin};
                    if isstruct(res) && isfield(res, 'success') && res.success ...
                            && isfield(res, metric) && ~isnan(res.(metric))
                        vals(end+1) = res.(metric); %#ok<AGROW>
                    end
                end
            end
        end
    end

    methods (Static, Access = private)
        function vals = collect_grid_values(results, cond_name, metric)
            % All successful metric values across the whole grid for one
            % condition. Reused to pool the param-space distribution across runs.
            vals = [];
            if ~isfield(results, cond_name); return; end
            rc = results.(cond_name);
            for k = 1:numel(rc)
                res = rc{k};
                if isstruct(res) && isfield(res, 'success') && res.success ...
                        && isfield(res, metric) && ~isnan(res.(metric))
                    vals(end+1) = res.(metric); %#ok<AGROW>
                end
            end
        end
    end

    %% Static Methods
    methods (Static)
        function obj = loadobj(s)
            % LOADOBJ Reconstruct object from struct when loading
            %
            % This method is called automatically by MATLAB's load() function
            % when loading a ParamSpaceAnalysis2 object that was saved with save().

            if isstruct(s)
                % Reconstruct object from struct
                obj = ParamSpaceAnalysis2();

                % Configuration Properties (public)
                if isfield(s, 'grid_params'), obj.grid_params = s.grid_params; end
                if isfield(s, 'param_ranges'), obj.param_ranges = s.param_ranges; end
                if isfield(s, 'n_levels'), obj.n_levels = s.n_levels; end
                if isfield(s, 'conditions'), obj.conditions = s.conditions; end
                if isfield(s, 'integer_params'), obj.integer_params = s.integer_params; end
                if isfield(s, 'explicit_vectors'), obj.explicit_vectors = s.explicit_vectors; end
                if isfield(s, 'vector_param_config'), obj.vector_param_config = s.vector_param_config; end
                if isfield(s, 'randomize_order'), obj.randomize_order = s.randomize_order; end
                if isfield(s, 'subset_fraction'), obj.subset_fraction = s.subset_fraction; end
                if isfield(s, 'network_seed_offset'), obj.network_seed_offset = s.network_seed_offset; end

                % Model Default Properties (public)
                if isfield(s, 'model_class'), obj.model_class = s.model_class; end
                if isfield(s, 'model_defaults'), obj.model_defaults = s.model_defaults; end
                if isfield(s, 'verbose'), obj.verbose = s.verbose; end

                % Execution Properties (public)
                if isfield(s, 'batch_size'), obj.batch_size = s.batch_size; end
                if isfield(s, 'output_dir'), obj.output_dir = s.output_dir; end
                if isfield(s, 'note'), obj.note = s.note; end
                if isfield(s, 'folder_prefix'), obj.folder_prefix = s.folder_prefix; end
                if isfield(s, 'store_local_lya'), obj.store_local_lya = s.store_local_lya; end
                if isfield(s, 'store_local_lya_dt'), obj.store_local_lya_dt = s.store_local_lya_dt; end
                if isfield(s, 'use_parallel'), obj.use_parallel = s.use_parallel; end

                % Results Properties (private) - need to set via internal assignment
                if isfield(s, 'results'), obj.results = s.results; end
                if isfield(s, 'has_run'), obj.has_run = s.has_run; end
                if isfield(s, 'analysis_start_time'), obj.analysis_start_time = s.analysis_start_time; end
                if isfield(s, 'param_vectors'), obj.param_vectors = s.param_vectors; end
                if isfield(s, 'all_configs'), obj.all_configs = s.all_configs; end
                if isfield(s, 'shuffled_indices'), obj.shuffled_indices = s.shuffled_indices; end
                if isfield(s, 'num_combinations'), obj.num_combinations = s.num_combinations; end
                if isfield(s, 'vector_param_lookup'), obj.vector_param_lookup = s.vector_param_lookup; end
                % Absent for runs saved before resolved_defaults existed; stays
                % empty, and effective_param / same_config fall back accordingly.
                if isfield(s, 'resolved_defaults'), obj.resolved_defaults = s.resolved_defaults; end
            else
                % Object was saved directly (already a ParamSpaceAnalysis2)
                obj = s;
            end
        end

        function result = run_single_job(job, model_defaults_local, grid_params_local, ...
                verbose_local, store_local_lya_local, store_local_lya_dt_local, ...
                vector_param_lookup_local, model_class_local)
            % RUN_SINGLE_JOB Execute a single simulation job
            % Extracted to allow use with both parfor and for loops

            run_start = tic;

            try
                % Build model arguments IN PRECEDENCE ORDER, weakest first, and
                % let the constructor's last-write-wins settle collisions:
                %
                %   model_defaults  <  condition  <  grid parameters
                %
                % Ordering by precedence rather than skipping duplicates matters
                % for more than tidiness. Grid parameters used to be appended
                % before model_defaults, so a swept property was assigned before
                % the defaults that give it meaning -- setting mu_EE_relative
                % before n_cellTypes, for instance. Feeding the defaults first
                % means the model is always fully configured by the time a swept
                % value lands on it.
                model_args = {'rng_seeds', [job.network_seed, job.network_seed + 1]};

                default_fields = fieldnames(model_defaults_local);
                for d_idx = 1:length(default_fields)
                    fname = default_fields{d_idx};
                    model_args = [model_args, {fname, model_defaults_local.(fname)}]; %#ok<AGROW>
                end

                % EVERY condition field except 'name' -- not a hardcoded list of
                % SRNNModel2's four adaptation counts. This is what lets a
                % condition carry a whole synapse_config struct for
                % SRNNCellTypePairs.
                cond_fields = setdiff(fieldnames(job.condition), {'name'})';
                for cf_idx = 1:numel(cond_fields)
                    field_name = cond_fields{cf_idx};
                    model_args = [model_args, {field_name, job.condition.(field_name)}]; %#ok<AGROW>
                end

                % Add grid parameters
                for p_idx = 1:length(grid_params_local)
                    pname = grid_params_local{p_idx};
                    if isfield(vector_param_lookup_local, pname)
                        % Vector parameter: look up pre-generated vector by index
                        vec_idx = job.config.(pname);
                        model_args = [model_args, {pname, vector_param_lookup_local.(pname){vec_idx}}]; %#ok<AGROW>
                    else
                        model_args = [model_args, {pname, job.config.(pname)}]; %#ok<AGROW>
                    end
                end

                % Create and run model
                model = feval(model_class_local, model_args{:});
                model.build();
                model.run();

                % Extract results
                result = struct();
                result.success = true;
                result.config = job.config;
                result.config_idx = job.config_idx;
                result.condition_name = job.condition.name;
                result.network_seed = job.network_seed;
                result.run_duration = toc(run_start);

                % Extract LLE
                if ~isempty(model.lya_results) && isfield(model.lya_results, 'LLE')
                    result.LLE = model.lya_results.LLE;
                else
                    result.LLE = NaN;
                end

                % Extract decimated local Lyapunov time series if requested
                if store_local_lya_local && ~isempty(model.lya_results)
                    if isfield(model.lya_results, 'local_lya') && isfield(model.lya_results, 't_lya')
                        current_lya_dt = model.lya_results.lya_dt;
                        deci_factor = max(1, round(store_local_lya_dt_local / current_lya_dt));

                        % Decimate local_lya and t_lya
                        result.local_lya = model.lya_results.local_lya(1:deci_factor:end);
                        result.t_lya = model.lya_results.t_lya(1:deci_factor:end);
                        result.local_lya_dt = current_lya_dt * deci_factor;
                    end
                end

                % Extract mean firing rate and mean synaptic output, pooled over
                % EVERY cell type rather than a hardcoded E/I pair: plot_data.r is
                % keyed by cell type name, so a 3-type model has .E/.PV/.SST and
                % the old .E/.I access would have silently missed a population.
                %
                % The synaptic-output field is called 'br' by SRNNModel2 (it is
                % literally b.*r) and 'synaptic_output' by SRNNCellTypePairs,
                % where facilitation makes 'br' a misnomer. Accept either.
                result.mean_rate = NaN;
                result.mean_synaptic_output = NaN;
                if ~isempty(model.plot_data) && isfield(model.plot_data, 'r')
                    result.mean_rate = ParamSpaceAnalysis2.pooled_mean(model.plot_data.r);
                    for cand = {'br', 'synaptic_output'}
                        if isfield(model.plot_data, cand{1})
                            result.mean_synaptic_output = ...
                                ParamSpaceAnalysis2.pooled_mean(model.plot_data.(cand{1}));
                            break;
                        end
                    end
                end

            catch ME
                result = struct();
                result.success = false;
                result.error_message = ME.message;
                result.config = job.config;
                result.config_idx = job.config_idx;
                result.condition_name = job.condition.name;
                result.network_seed = job.network_seed;
                result.run_duration = toc(run_start);
                result.LLE = NaN;
                result.mean_rate = NaN;
                result.mean_synaptic_output = NaN;

                if verbose_local
                    fprintf('  ERROR config %d: %s\n', job.config_idx, ME.message);
                end
            end
        end
    end

    %% Private Methods
    methods (Access = private)
        function set_default_conditions(obj)
            % SET_DEFAULT_CONDITIONS Initialize the four adaptation conditions

            obj.conditions = { ...
                struct('name', 'no_adaptation', 'n_a_E', 0, 'n_b_E', 0), ...
                struct('name', 'sfa_only',      'n_a_E', 3, 'n_b_E', 0), ...
                struct('name', 'std_only',      'n_a_E', 0, 'n_b_E', 1), ...
                struct('name', 'sfa_and_std',   'n_a_E', 3, 'n_b_E', 1) ...
                };
        end

    end

    %% Grid construction (public: needed to size a run before committing to it)
    methods
        function generate_grid(obj)
            % GENERATE_GRID Create the multi-dimensional parameter grid.
            %
            % run() calls this itself, so a normal sweep never needs to. It is
            % public because materialising the grid WITHOUT running it is how
            % you find out what a run will cost -- num_combinations and
            % numel(shuffled_indices) after this call are the full grid size
            % and the number of points subset_fraction will actually simulate.
            % On a 7-axis grid that difference is the whole planning question.
            %
            % Calling it again rebuilds the grid and redraws the execution
            % order, discarding any previous schedule.

            n_params = length(obj.grid_params);
            obj.param_vectors = cell(1, n_params);

            for i = 1:n_params
                param_name = obj.grid_params{i};
                is_int = ismember(param_name, obj.integer_params);

                if isfield(obj.vector_param_config, param_name)
                    % Vector parameter: use integer indices as grid values
                    obj.param_vectors{i} = 1:obj.n_levels;

                    % Pre-generate vectors for each level
                    vpc = obj.vector_param_config.(param_name);

                    % Generate vary values across levels
                    if strcmp(vpc.level_spacing, 'log')
                        vary_values = logspace(log10(vpc.vary_range(1)), log10(vpc.vary_range(2)), obj.n_levels);
                    else
                        vary_values = linspace(vpc.vary_range(1), vpc.vary_range(2), obj.n_levels);
                    end
                    if is_int
                        vary_values = round(vary_values);
                    end

                    % Generate vectors for each level
                    vectors = cell(obj.n_levels, 1);
                    for lev = 1:obj.n_levels
                        if strcmp(vpc.vary_element, 'last')
                            start_val = vpc.fixed_value;
                            end_val = vary_values(lev);
                        else  % 'first'
                            start_val = vary_values(lev);
                            end_val = vpc.fixed_value;
                        end

                        if strcmp(vpc.spacing, 'log')
                            vectors{lev} = logspace(log10(start_val), log10(end_val), vpc.n_elements);
                        else
                            vectors{lev} = linspace(start_val, end_val, vpc.n_elements);
                        end
                        if is_int
                            vectors{lev} = round(vectors{lev});
                        end
                    end
                    obj.vector_param_lookup.(param_name) = vectors;

                elseif isfield(obj.explicit_vectors, param_name)
                    % Use explicit vector directly
                    obj.param_vectors{i} = obj.explicit_vectors.(param_name);
                else
                    % Generate from range using n_levels
                    param_range = obj.param_ranges.(param_name);
                    if is_int
                        obj.param_vectors{i} = round(linspace(param_range(1), param_range(2), obj.n_levels));
                    else
                        obj.param_vectors{i} = linspace(param_range(1), param_range(2), obj.n_levels);
                    end
                end
            end

            % Generate all combinations using ndgrid
            grid_cells = cell(size(obj.param_vectors));
            [grid_cells{:}] = ndgrid(obj.param_vectors{:});

            obj.num_combinations = numel(grid_cells{1});
            obj.all_configs = cell(obj.num_combinations, 1);

            for i = 1:obj.num_combinations
                config = struct();
                for j = 1:n_params
                    config.(obj.grid_params{j}) = grid_cells{j}(i);
                end
                obj.all_configs{i} = config;
            end

            % Set execution order
            if obj.randomize_order
                rng('shuffle');
                obj.shuffled_indices = randperm(obj.num_combinations);
                fprintf('Generated %d parameter combinations (randomized order)\n', obj.num_combinations);
            else
                obj.shuffled_indices = 1:obj.num_combinations;
                fprintf('Generated %d parameter combinations (sequential order)\n', obj.num_combinations);
            end

            % Thin to a random subset, if asked. shuffled_indices is the single
            % source of truth for what runs -- num_combinations stays the FULL
            % grid size, so results still land at their true grid positions and
            % the unrun points are simply empty. See the subset_fraction
            % property comment for why that is safe.
            obj.apply_subset_fraction();
        end

        function apply_subset_fraction(obj)
            % APPLY_SUBSET_FRACTION Truncate shuffled_indices to a random subset.
            if isempty(obj.subset_fraction) || ~isscalar(obj.subset_fraction) || ...
                    ~isnumeric(obj.subset_fraction) || ~isfinite(obj.subset_fraction) || ...
                    obj.subset_fraction <= 0 || obj.subset_fraction > 1
                error('ParamSpaceAnalysis2:BadSubsetFraction', ...
                    'subset_fraction must be a scalar in (0, 1]; got %s.', ...
                    mat2str(obj.subset_fraction));
            end
            if obj.subset_fraction == 1
                return;   % full grid: shuffled_indices already correct
            end
            if ~obj.randomize_order
                error('ParamSpaceAnalysis2:SubsetNeedsRandomOrder', ...
                    ['subset_fraction = %g requires randomize_order = true.\n' ...
                     'With sequential order a subset is the FIRST %d configs in ' ...
                     'ndgrid order -- a systematic corner of the grid, not a ' ...
                     'random sample of it.'], ...
                    obj.subset_fraction, ceil(obj.subset_fraction * obj.num_combinations));
            end

            % ceil, not round: a fraction small enough to round to zero should
            % still run one point rather than silently running nothing.
            n_run = min(obj.num_combinations, ...
                max(1, ceil(obj.subset_fraction * obj.num_combinations)));
            obj.shuffled_indices = obj.shuffled_indices(1:n_run);
            fprintf(['Subset: running %d of %d combinations (%.1f%% requested, ' ...
                '%.1f%% actual); the rest stay empty.\n'], ...
                n_run, obj.num_combinations, 100*obj.subset_fraction, ...
                100*n_run/obj.num_combinations);
        end
    end

    %% Private Methods (continued)
    methods (Access = private)
        function run_batched_simulation(obj, temp_dir)
            % RUN_BATCHED_SIMULATION Execute simulations in batches with checkpoints

            % Batch over what will actually RUN, not over the full grid.
            % numel(shuffled_indices) == num_combinations unless subset_fraction
            % thinned it, so this is unchanged for a full run.
            n_run = numel(obj.shuffled_indices);
            num_batches = ceil(n_run / obj.batch_size);
            conditions_local = obj.conditions;
            num_conditions = length(conditions_local);

            if n_run < obj.num_combinations
                fprintf('Running %d of %d combinations (subset) in %d batches...\n', ...
                    n_run, obj.num_combinations, num_batches);
            else
                fprintf('Running %d combinations in %d batches...\n', n_run, num_batches);
            end

            for batch_idx = 1:num_batches
                batch_file = fullfile(temp_dir, sprintf('batch_%d.mat', batch_idx));

                % Skip if batch already completed (resume capability)
                if exist(batch_file, 'file')
                    fprintf('Batch %d/%d already completed. Skipping.\n', batch_idx, num_batches);
                    continue;
                end

                start_idx = (batch_idx - 1) * obj.batch_size + 1;
                end_idx = min(batch_idx * obj.batch_size, n_run);
                batch_indices = obj.shuffled_indices(start_idx:end_idx);
                current_batch_size = length(batch_indices);

                fprintf('\n--- Batch %d/%d (configs %d-%d) ---\n', ...
                    batch_idx, num_batches, start_idx, end_idx);

                % Create jobs: each config runs ALL conditions with SAME network seed
                total_jobs = current_batch_size * num_conditions;
                jobs = cell(total_jobs, 1);
                job_idx = 1;

                for k = 1:current_batch_size
                    config_idx = batch_indices(k);
                    config = obj.all_configs{config_idx};

                    % Same network seed for all conditions of this config.
                    % network_seed_offset (0 by default) shifts every seed so
                    % separate runs of the same config use distinct networks.
                    network_seed = config_idx*100 + obj.network_seed_offset;

                    for c_idx = 1:num_conditions
                        job = struct();
                        job.config = config;
                        job.config_idx = config_idx;
                        job.condition = conditions_local{c_idx};
                        job.condition_idx = c_idx;
                        job.network_seed = network_seed;
                        jobs{job_idx} = job;
                        job_idx = job_idx + 1;
                    end
                end

                % Extract values for parfor/for loop
                model_defaults_local = obj.model_defaults;
                grid_params_local = obj.grid_params;
                verbose_local = obj.verbose;
                store_local_lya_local = obj.store_local_lya;
                store_local_lya_dt_local = obj.store_local_lya_dt;
                vector_param_lookup_local = obj.vector_param_lookup;
                model_class_local = obj.model_class;   % char, so parfor-safe

                % Determine execution mode
                run_parallel = obj.use_parallel && canUseParallelPool;

                if obj.use_parallel && ~canUseParallelPool && batch_idx == 1
                    warning('ParamSpaceAnalysis2:NoParallelPool', ...
                        'Parallel pool not available. Falling back to sequential execution.');
                end

                % Run simulation loop
                parallel_results = cell(total_jobs, 1);
                batch_start = tic;

                if run_parallel
                    parfor j = 1:total_jobs
                        parallel_results{j} = ParamSpaceAnalysis2.run_single_job(...
                            jobs{j}, model_defaults_local, grid_params_local, ...
                            verbose_local, store_local_lya_local, store_local_lya_dt_local, ...
                            vector_param_lookup_local, model_class_local);
                    end
                else
                    for j = 1:total_jobs
                        parallel_results{j} = ParamSpaceAnalysis2.run_single_job(...
                            jobs{j}, model_defaults_local, grid_params_local, ...
                            verbose_local, store_local_lya_local, store_local_lya_dt_local, ...
                            vector_param_lookup_local, model_class_local);
                    end
                end

                batch_elapsed = toc(batch_start);

                % Organize results by condition
                batch_results = struct();
                for c_idx = 1:num_conditions
                    cond_name = conditions_local{c_idx}.name;
                    batch_results.(cond_name) = {};
                end

                for j = 1:total_jobs
                    res = parallel_results{j};
                    cond_name = res.condition_name;
                    batch_results.(cond_name){end+1} = res;
                end

                % Save batch checkpoint
                save(batch_file, 'batch_results', 'batch_indices', '-v7.3');

                % Count successes
                n_success = sum(cellfun(@(r) r.success, parallel_results));
                fprintf('Batch %d completed in %.1f min (%d/%d successful)\n', ...
                    batch_idx, batch_elapsed/60, n_success, total_jobs);
            end
        end

        function save_summary(obj)
            % SAVE_SUMMARY Write param_space_summary.mat -- METADATA ONLY.
            %
            % NOT a restore path. It is a readable record of what a run was and
            % how it went (grid, conditions, resolved parameters, per-condition
            % success counts, timings), for inspection and for tooling that wants
            % those facts without loading a whole object.
            %
            % Reconstructing a PSA from it is what caused the two loader bugs
            % fixed on 2026-08-14: two separate methods each restored their own
            % hand-picked subset of these fields, drifted from the twelve written
            % here, and silently dropped model_class and vector_param_lookup.
            % That is why the rule below is a rule. psa_object.mat is the
            % authoritative artifact; load runs with
            % ParamSpaceAnalysis2.from_dir. Do not add restore logic here.

            summary_file = fullfile(obj.output_dir, 'param_space_summary.mat');

            summary_data = struct();
            summary_data.grid_params = obj.grid_params;
            summary_data.param_ranges = obj.param_ranges;
            summary_data.param_vectors = obj.param_vectors;
            summary_data.n_levels = obj.n_levels;
            summary_data.num_combinations = obj.num_combinations;
            summary_data.conditions = obj.conditions;
            summary_data.model_class = obj.model_class;
            summary_data.model_defaults = obj.model_defaults;
            summary_data.resolved_defaults = obj.resolved_defaults;
            summary_data.shuffled_indices = obj.shuffled_indices;
            % So a run directory says outright that it is a partial grid --
            % without this, a thinned histogram is indistinguishable from a
            % full one that mostly failed.
            summary_data.subset_fraction = obj.subset_fraction;
            summary_data.num_run = numel(obj.shuffled_indices);
            summary_data.analysis_start_time = obj.analysis_start_time;
            summary_data.analysis_completed = datestr(now);

            % Compute statistics per condition
            for c_idx = 1:length(obj.conditions)
                cond_name = obj.conditions{c_idx}.name;
                if isfield(obj.results, cond_name)
                    results = obj.results.(cond_name);
                    % n_total is what was ATTEMPTED, not the cell length. The
                    % results cell is always full-grid length so results sit at
                    % their true config_idx, so length(results) on a subsetted
                    % run is the grid size -- which made a flawless 64/64 record
                    % success_rate = 0.0039 and read as a 99.6% failure.
                    n_success = sum(cellfun(@(r) isstruct(r) && isfield(r, 'success') && r.success, results));
                    n_attempted = numel(obj.shuffled_indices);
                    if n_attempted == 0, n_attempted = length(results); end   % legacy/unscheduled
                    summary_data.stats.(cond_name).n_success = n_success;
                    summary_data.stats.(cond_name).n_total = n_attempted;
                    summary_data.stats.(cond_name).n_grid_points = length(results);
                    summary_data.stats.(cond_name).success_rate = n_success / n_attempted;
                end
            end

            save(summary_file, 'summary_data', '-v7.3');
            fprintf('Summary saved to: %s\n', summary_file);
        end
    end
end
