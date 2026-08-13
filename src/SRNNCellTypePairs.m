classdef SRNNCellTypePairs < handle
    %SRNNCELLTYPEPAIRS Stable recurrent nonlinear network with arbitrary cell types.
    %
    % Intrinsic parameters are ordered according to cell_type_names. Synaptic
    % dynamics are configured by named presynaptic-to-postsynaptic routes,
    % for example synapse_config.E.PV.std and synapse_config.E.PV.stf.
    %
    % State ordering is timescale-major within each block:
    %   S = [a{1}(:); ...; a{C}(:); b{1,1}(:); ...; b{C,C}(:); ...
    %        g{1,1}(:); ...; g{C,C}(:); x(:)]
    % Pair arrays use (presynaptic type, postsynaptic type) ordering.
    %
    % Example:
    %   model = SRNNCellTypePairs( ...
    %       'n_cellTypes', 3, 'cell_type_names', {'E','PV','SST'}, ...
    %       'f', [.8 .1 .1], 'mu_tilde', [.08 -.12 -.1], ...
    %       'sigma_tilde', [.02 .02 .02], 'n_a', [3 0 1], ...
    %       'tau_a', {[.25 1 10], [], 2}, ...
    %       'synapse_config', synapse_config);
    %   model.build(); model.run(); model.plot();

    %% Cell types and connectivity
    properties
        n = 300
        indegree = 100
        n_cellTypes
        cell_type_names
        f

        % RMT block statistics in multiples of F = default_val, indexed
        % (POSTsynaptic, PREsynaptic) as in RMTBlocks -- so mu_tilde_relative(1,2)
        % is "type 1 receives from type 2". A 1 x C row is accepted and broadcast
        % down the columns, which is the classical presynaptic-only setting.
        %
        % The absolute mu_tilde / sigma_tilde are Dependent on these; storing the
        % multiplier is what lets them be swept (F depends on n and indegree) and
        % recorded. Dale's law is a COLUMN constraint: signs must stay constant
        % down each column.
        mu_tilde_relative
        sigma_tilde_relative

        % Which network F is computed for. true (default) tracks the current
        % n/indegree, making the bulk radius R independent of n; false pins F to
        % the reference network below. Freezing F does NOT freeze the network --
        % connectivity still tracks n/indegree, only the weight SCALE is pinned.
        F_tracks_network = true
        F_ref_n          = 300
        F_ref_indegree   = 100

        level_of_chaos = 1.0
        rescale_by_abscissa = false
        zrs_mode = 'none'           % 'none' | 'ZRS' | 'SZRS' | 'Partial_SZRS'
    end

    %% Per-type SFA and pair-specific STD/STF
    properties
        n_a
        tau_a
        c
        synapse_config = struct()
        std_zero_floor = false
    end

    %% Shared dynamics and simulation settings
    properties
        tau_d = 0.1

        % The nonlinearity is chosen by NAME and parameterised by S_a / S_c; the
        % activation_function / _derivative handles are Dependent and rebuilt
        % from these, so S_a and S_c always describe the function actually in
        % use. 'logistic' uses S_c only, 'piecewise' uses both, 'tanh' neither.
        % activation_custom = {fn, dfn} overrides the name when non-empty.
        activation = 'logistic'         % 'logistic' | 'piecewise' | 'tanh'
        activation_custom = {}
        S_a = 0.9
        S_c = 0.40                      % matches SRNNModel2's default

        % Per-type setpoints, 1 x C each (a scalar is broadcast to every type).
        % By default every neuron shares the scalar S_c above; give a type a
        % mean and/or a standard deviation here and build() draws a length-n
        % setpoint vector instead:
        %
        %   S_c_i = mu_S_c(q) + sigma_S_c(q) * randn   for each neuron of type q
        %
        % An empty mu_S_c falls back to S_c, so sigma alone adds spread around
        % the shared centre and mu alone gives the types different excitability
        % with no spread. The realized vector is S_c_vec (read-only, empty
        % unless heterogeneity was requested).
        %
        % Only 'logistic' and 'piecewise' take a centre: requesting
        % heterogeneity with 'tanh' or with activation_custom is an error. In
        % heterogeneous mode the activation_function handle is only valid for
        % length-n input.
        mu_S_c = []                     % 1 x C means ([] = use S_c for all types)
        sigma_S_c = []                  % 1 x C std devs ([] = zeros, no spread)

        % Seed for the setpoint draw. Empty derives it from rng_seeds(1), so a
        % sweep over network seeds also varies the setpoints while keeping the
        % draw off the stream that generated W. Set it to pin the draw.
        S_c_seed = []
        fs = 400
        T_range = [0, 50]
        T_plot
        ode_solver = @ode45
        ode_opts
        x0_std = 0.1
        input_config
        u_ex_scale = 1.0
        rng_seeds = [1 2]
        reps = 1
        lya_method = 'benettin'
        lya_T_interval
        filter_local_lya = false
        lya_filter_order = 2
        lya_filter_cutoff = 0.25
        store_full_state = false
        store_decimated_state = true
        plot_deci
        plot_freq = 10
    end

    properties (Dependent)
        alpha
        default_val                     % F = 1/sqrt(N*alpha*(2-alpha)) (Harris 2023)
        mu_tilde                        % C x C absolute block means    = relative * F
        sigma_tilde                     % C x C absolute block std devs = relative * F
        R                               % Bulk spectral radius (generalizes Harris Eq. 18)
        lambda_O                        % Outlier eigenvalues, descending by magnitude

        % Scalar aliases into mu_tilde_relative / sigma_tilde_relative for the
        % two-type case, named (post, pre) -- mu_EI is "E receives from I".
        % These exist so a single block can be a SWEEP AXIS: ParamSpaceAnalysis2
        % passes each grid parameter as a constructor name-value pair, so the
        % thing being swept has to be a settable property in its own right.
        % Reading or writing one on a model with other than 2 types errors.
        mu_EE_relative
        mu_EI_relative
        mu_IE_relative
        mu_II_relative
        sigma_EE_relative
        sigma_EI_relative
        sigma_IE_relative
        sigma_II_relative

        % Two more aliases for the same reason: a sweep axis has to be a
        % settable property, and f / tau_a are a row and a cell array.
        %   f_E     -- the excitatory fraction as a SCALAR; writing it sets
        %              f = [v, 1-v], so the two types stay a partition. Two-type
        %              only, like the block aliases above.
        %   tau_a_E -- tau_a of the FIRST cell type as a plain numeric row, so
        %              add_vector_parameter can sweep the SFA timescales.
        %              Guarded on the first type actually being named E.
        f_E
        tau_a_E
        activation_function             % Built from activation + S_a/S_c
        activation_function_derivative
        n_per_type
        type_indices
        cell_indices
        n_b_pairs
        n_g_pairs
        N_sys_eqs
    end

    properties (SetAccess = protected)
        % Realized per-neuron nonlinearity setpoints (n x 1), drawn by
        % build_setpoints() from mu_S_c / sigma_S_c. EMPTY means every neuron
        % shares the scalar S_c, which is the default and keeps every code path
        % on the scalar branch.
        S_c_vec = []

        W
        is_built = false
        t_ex
        u_ex
        u_interpolant
        S0
        cached_params
        t_out
        S_out
        plot_data
        lya_results
        has_run = false
        dead_state_count = 0
    end

    methods
        function obj = SRNNCellTypePairs(varargin)
            obj.set_defaults();
            if mod(numel(varargin), 2) ~= 0
                error('SRNNCellTypePairs:InvalidInput', ...
                    'Constructor arguments must be name-value pairs.');
            end
            supplied = string(varargin(1:2:end));
            required = ["n_cellTypes", "cell_type_names", "f", ...
                "mu_tilde_relative", "sigma_tilde_relative"];
            missing = required(~ismember(required, supplied));
            if ~isempty(missing)
                error('SRNNCellTypePairs:MissingConfig', ...
                    'Required constructor properties: %s.', strjoin(cellstr(missing), ', '));
            end
            % These became Dependent (computed from the settable names given
            % here), so assigning them would raise MATLAB's generic read-only
            % error. isprop() is true for Dependent properties, so check first.
            renamed = struct( ...
                'mu_tilde',    'mu_tilde_relative', ...
                'sigma_tilde', 'sigma_tilde_relative', ...
                'activation_function',            'activation', ...
                'activation_function_derivative', 'activation');
            for k = 1:2:numel(varargin)
                name = varargin{k};
                if ~(ischar(name) || (isstring(name) && isscalar(name)))
                    error('SRNNCellTypePairs:UnknownProperty', ...
                        'Property names must be character vectors or strings.');
                end
                name = char(name);
                if isfield(renamed, name)
                    error('SRNNCellTypePairs:RenamedProperty', ...
                        ['''%s'' is now a computed (Dependent) property. Set ''%s'' ' ...
                        'instead -- the tildes are held in multiples of ' ...
                        'F = default_val, and the nonlinearity is chosen by name.'], ...
                        name, renamed.(name));
                end
                if ~isprop(obj, name)
                    error('SRNNCellTypePairs:UnknownProperty', ...
                        'Unknown property: %s', name);
                end
                obj.(name) = varargin{k + 1};
            end
            obj.complete_type_defaults();
            obj.validate();
            if isempty(obj.plot_deci)
                obj.plot_deci = max(1, round(obj.fs / obj.plot_freq));
            end
        end

        function val = get.default_val(obj)
            % F = 1/sqrt(N*alpha*(2-alpha)), the factor that yields R = 1 when
            % every tilde equals 1 (Harris 2023). Unchanged by block or
            % multi-type structure: setting all tildes equal reduces the block
            % radius to Harris Eq. (18).
            %
            % F_tracks_network = false pins F to (F_ref_n, F_ref_indegree) so a
            % sweep over n or indegree leaves the weight distribution fixed.
            % This does NOT freeze the network: build_network still passes the
            % real alpha, so connectivity tracks the current n/indegree.
            if obj.F_tracks_network
                n_F = obj.n;         a_F = obj.alpha;
            else
                n_F = obj.F_ref_n;   a_F = obj.F_ref_indegree / obj.F_ref_n;
            end
            val = 1 / sqrt(n_F * a_F * (2 - a_F));
        end

        function val = get.mu_tilde(obj)
            val = SRNNCellTypePairs.scale_tilde_mat( ...
                obj.expand_block(obj.mu_tilde_relative), obj.default_val);
        end

        function val = get.sigma_tilde(obj)
            val = SRNNCellTypePairs.scale_tilde_mat( ...
                obj.expand_block(obj.sigma_tilde_relative), obj.default_val);
        end

        function val = get.R(obj)
            % Perron root of the block variance matrix V(a,b) = n*f_b*sigma_s_sq,
            % generalizing Harris Eq. (18); see RMTBlocks.R.
            V = obj.n * (obj.f .* obj.block_sigma_s_sq());
            val = sqrt(max(0, max(real(eig(V))))) * obj.level_of_chaos;
        end

        function val = get.lambda_O(obj)
            % Outliers of E[W], descending by MAGNITUDE (in an
            % inhibition-dominated network the dominant one is large and
            % negative). Generalizes Harris Eq. (17); see RMTBlocks.lambda_O.
            K = obj.n * (obj.f .* (obj.alpha * obj.mu_tilde));
            lam = eig(K) * obj.level_of_chaos;
            [~, order] = sort(abs(lam), 'descend');
            val = lam(order);
        end

        % --- Two-type scalar block aliases -----------------------------------
        function v = get.mu_EE_relative(obj),    v = obj.get_block('mu', 1, 1); end
        function v = get.mu_EI_relative(obj),    v = obj.get_block('mu', 1, 2); end
        function v = get.mu_IE_relative(obj),    v = obj.get_block('mu', 2, 1); end
        function v = get.mu_II_relative(obj),    v = obj.get_block('mu', 2, 2); end
        function v = get.sigma_EE_relative(obj), v = obj.get_block('sigma', 1, 1); end
        function v = get.sigma_EI_relative(obj), v = obj.get_block('sigma', 1, 2); end
        function v = get.sigma_IE_relative(obj), v = obj.get_block('sigma', 2, 1); end
        function v = get.sigma_II_relative(obj), v = obj.get_block('sigma', 2, 2); end

        function set.mu_EE_relative(obj, v),    obj.set_block('mu', 1, 1, v); end
        function set.mu_EI_relative(obj, v),    obj.set_block('mu', 1, 2, v); end
        function set.mu_IE_relative(obj, v),    obj.set_block('mu', 2, 1, v); end
        function set.mu_II_relative(obj, v),    obj.set_block('mu', 2, 2, v); end
        function set.sigma_EE_relative(obj, v), obj.set_block('sigma', 1, 1, v); end
        function set.sigma_EI_relative(obj, v), obj.set_block('sigma', 1, 2, v); end
        function set.sigma_IE_relative(obj, v), obj.set_block('sigma', 2, 1, v); end
        function set.sigma_II_relative(obj, v), obj.set_block('sigma', 2, 2, v); end

        % --- Scalar / row aliases for the remaining per-type sweep axes -------
        function v = get.f_E(obj)
            obj.assert_two_types();
            if isempty(obj.f)
                v = [];
            else
                v = obj.f(1);
            end
        end

        function set.f_E(obj, v)
            obj.assert_two_types();
            if ~isscalar(v) || ~isfinite(v) || v <= 0 || v >= 1
                error('SRNNCellTypePairs:InvalidParams', ...
                    'f_E must be a scalar strictly between 0 and 1.');
            end
            obj.f = [v, 1 - v];
        end

        function v = get.tau_a_E(obj)
            obj.assert_first_type_is_E();
            if isempty(obj.tau_a)
                v = [];
            else
                v = obj.tau_a{1};
            end
        end

        function set.tau_a_E(obj, v)
            obj.assert_first_type_is_E();
            % tau_a is completed to a 1 x C cell by complete_type_defaults, which
            % runs AFTER the constructor's name-value loop -- so this setter has
            % to be able to create the cell itself.
            if ~iscell(obj.tau_a) || numel(obj.tau_a) ~= obj.n_cellTypes
                obj.tau_a = cell(1, obj.n_cellTypes);
            end
            obj.tau_a{1} = reshape(v, 1, []);
        end

        function val = get.activation_function(obj)
            val = obj.build_activation(1);
        end

        function val = get.activation_function_derivative(obj)
            val = obj.build_activation(2);
        end

        function value = get.alpha(obj)
            value = obj.indegree / obj.n;
        end

        function value = get.n_per_type(obj)
            if isempty(obj.f) || isempty(obj.n)
                value = [];
            else
                value = RMTCellTypes.allocate_counts(obj.n, obj.f);
            end
        end

        function value = get.type_indices(obj)
            if isempty(obj.n_per_type)
                value = {};
            else
                value = RMTCellTypes.make_type_indices(obj.n_per_type);
            end
        end

        function value = get.cell_indices(obj)
            value = struct();
            if isempty(obj.cell_type_names)
                return;
            end
            indices = obj.type_indices;
            for q = 1:obj.n_cellTypes
                value.(obj.cell_type_names{q}) = indices{q};
            end
        end

        function value = get.N_sys_eqs(obj)
            if isempty(obj.n_a) || isempty(obj.n_cellTypes)
                value = obj.n;
            else
                pair = obj.compile_synapse_config();
                pair_counts = sum(pair.n_b + pair.n_g, 2)';
                value = obj.n + sum(obj.n_per_type .* obj.n_a) + ...
                    sum(obj.n_per_type .* pair_counts);
            end
        end

        function value = get.n_b_pairs(obj)
            value = obj.compile_synapse_config().n_b;
        end

        function value = get.n_g_pairs(obj)
            value = obj.compile_synapse_config().n_g;
        end

        function build(obj)
            obj.complete_type_defaults();
            obj.validate();
            obj.build_network();
            obj.build_stimulus();
            obj.cached_params = obj.get_params();
            obj.is_built = true;
            fprintf('SRNNCellTypePairs built successfully. Ready to run.\n');
        end

        function run(obj)
            if ~obj.is_built
                error('SRNNCellTypePairs:NotBuilt', ...
                    'Model must be built before running. Call build() first.');
            end

            params = obj.cached_params;
            params.u_interpolant = obj.u_interpolant;
            dt = 1 / obj.fs;
            if isempty(obj.ode_opts)
                jacobian = @(~, S) SRNNCellTypePairs.compute_Jacobian_fast(S, params);
                obj.ode_opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, ...
                    'MaxStep', dt, 'Jacobian', jacobian);
            end
            rhs = @(t, S) SRNNCellTypePairs.dynamics_fast(t, S, params);

            fprintf('Integrating SRNNCellTypePairs equations\n');
            tic;
            [t_raw, S_raw] = obj.ode_solver(rhs, obj.t_ex, obj.S0, obj.ode_opts);
            fprintf('Integration complete in %.2f seconds.\n', toc);
            if numel(t_raw) ~= numel(obj.t_ex) || ...
                    max(abs(t_raw(:) - obj.t_ex(:))) > 1e-9
                error('SRNNCellTypePairs:TimeMismatch', ...
                    'ODE solver output times do not match requested times.');
            end

            obj.t_out = t_raw;
            obj.S_out = S_raw;
            if ~strcmpi(obj.lya_method, 'none')
                obj.compute_lyapunov();
                if obj.filter_local_lya
                    obj.filter_lyapunov();
                end
            else
                obj.lya_results = struct();
            end
            if obj.store_decimated_state
                obj.decimate_and_unpack();
            end
            if ~obj.store_full_state
                obj.S_out = [];
            end
            obj.has_run = true;
            fprintf('Simulation complete.\n');
        end

        function compute_lyapunov(obj)
            if isempty(obj.S_out)
                error('SRNNCellTypePairs:NoStateData', 'State data are not available.');
            end
            if isempty(obj.lya_T_interval)
                if obj.T_range(2) >= 15
                    obj.lya_T_interval = [obj.T_range(1) + 15, obj.T_range(2)];
                else
                    obj.lya_T_interval = obj.T_range;
                end
            end
            params = obj.cached_params;
            params.u_interpolant = obj.u_interpolant;
            rhs = @(t, S) SRNNCellTypePairs.dynamics_fast(t, S, params);
            obj.lya_results = SRNNCellTypePairs.compute_lyapunov_exponents_internal( ...
                obj.lya_method, obj.S_out, obj.t_out, 1 / obj.fs, obj.fs, ...
                obj.lya_T_interval, params, obj.ode_opts, obj.ode_solver, rhs);
            if isfield(obj.lya_results, 'LLE')
                fprintf('Largest Lyapunov Exponent: %.4f\n', obj.lya_results.LLE);
            end
        end

        function filter_lyapunov(obj)
            if isempty(obj.lya_results) || ~isfield(obj.lya_results, 'lya_fs')
                return;
            end
            Wn = obj.lya_filter_cutoff / (obj.lya_results.lya_fs / 2);
            [bf, af] = butter(obj.lya_filter_order, Wn, 'low');
            if isfield(obj.lya_results, 'local_lya')
                obj.lya_results.local_lya = filtfilt(bf, af, obj.lya_results.local_lya);
            end
            if isfield(obj.lya_results, 'local_LE_spectrum_t')
                for col = 1:size(obj.lya_results.local_LE_spectrum_t, 2)
                    obj.lya_results.local_LE_spectrum_t(:, col) = filtfilt( ...
                        bf, af, obj.lya_results.local_LE_spectrum_t(:, col));
                end
            end
        end

        function params = get_params(obj)
            params = struct();
            params.n = obj.n;
            params.indegree = obj.indegree;
            params.alpha = obj.alpha;
            params.n_cellTypes = obj.n_cellTypes;
            params.cell_type_names = obj.cell_type_names;
            params.f = obj.f;
            params.n_per_type = obj.n_per_type;
            params.type_indices = obj.type_indices;
            params.cell_indices = obj.cell_indices;
            params.mu_tilde = obj.mu_tilde;
            params.sigma_tilde = obj.sigma_tilde;
            params.level_of_chaos = obj.level_of_chaos;
            params.n_a = obj.n_a;
            params.tau_a = obj.tau_a;
            params.c = obj.c;
            pair = obj.compile_synapse_config();
            params.n_b_pairs = pair.n_b;
            params.tau_b_rec = pair.tau_b_rec;
            params.tau_b_rel = pair.tau_b_rel;
            params.std_zero_floor = obj.std_zero_floor;
            params.n_g_pairs = pair.n_g;
            params.tau_g_dec = pair.tau_g_dec;
            params.tau_g_fac = pair.tau_g_fac;
            params.G = pair.G;
            params.tau_d = obj.tau_d;
            params.activation_function = obj.activation_function;
            params.activation_function_derivative = obj.activation_function_derivative;
            % The realized per-neuron setpoints, for diagnostics: the handles
            % above already carry them. Empty in the homogeneous case.
            params.S_c_vec = obj.S_c_vec;
            params.x0_std = obj.x0_std;
            params.N_sys_eqs = obj.N_sys_eqs;
            params.state_layout = SRNNCellTypePairs.make_state_layout( ...
                obj.n, obj.n_per_type, obj.n_a, pair.n_b, pair.n_g);
            params.rng_seeds = obj.rng_seeds;
            if ~isempty(obj.W)
                params.W = obj.W;
            end
        end

        function dS_dt = dynamics(obj, t, S)
            if ~obj.is_built
                error('SRNNCellTypePairs:NotBuilt', 'Call build() before dynamics().');
            end
            params = obj.cached_params;
            params.u_interpolant = obj.u_interpolant;
            dS_dt = SRNNCellTypePairs.dynamics_fast(t, S, params);
        end

        function clear_results(obj)
            obj.t_out = [];
            obj.S_out = [];
            obj.plot_data = [];
            obj.lya_results = [];
            obj.has_run = false;
        end

        function reset(obj)
            obj.is_built = false;
            obj.W = [];
            obj.t_ex = [];
            obj.u_ex = [];
            obj.u_interpolant = [];
            obj.S0 = [];
            obj.cached_params = [];
            obj.dead_state_count = 0;
            obj.clear_results();
        end
    end

    methods (Access = protected)
        function v = get_block(obj, which_stat, post, pre)
            % GET_BLOCK One entry of a _relative block matrix, expanded first so a
            % 1 x C presynaptic row reads correctly.
            obj.assert_two_types();
            M = obj.expand_block(obj.([which_stat '_tilde_relative']));
            v = M(post, pre);
        end

        function set_block(obj, which_stat, post, pre, v)
            % SET_BLOCK Write one entry, expanding a 1 x C row to a full block
            % first so that setting one route does not silently change the other
            % postsynaptic row.
            obj.assert_two_types();
            name = [which_stat '_tilde_relative'];
            M = obj.expand_block(obj.(name));
            if isempty(M)
                error('SRNNCellTypePairs:InvalidParams', ...
                    'Set %s before setting one of its blocks.', name);
            end
            M(post, pre) = v;
            obj.(name) = M;
        end

        function assert_two_types(obj)
            if ~isequal(obj.n_cellTypes, 2)
                error('SRNNCellTypePairs:NotTwoTypes', ...
                    ['The E/I scalar aliases require exactly 2 cell types ' ...
                    '(this model has %s). Use f / mu_tilde_relative / ' ...
                    'sigma_tilde_relative directly.'], ...
                    mat2str(obj.n_cellTypes));
            end
        end

        function assert_first_type_is_E(obj)
            % tau_a_E indexes tau_a{1}, so the name is only honest if type 1 is
            % actually E. Any number of types is fine -- unlike the block
            % aliases, this one does not need a 2 x 2.
            names = obj.cell_type_names;
            % cellstr, because the constructor may not have reached
            % complete_type_defaults yet -- names can still be a char row or a
            % string array at this point rather than the cellstr it ends up as.
            % Anything else malformed falls through to the error below rather
            % than throwing cellstr's; validate() reports it properly.
            usable = ~isempty(names) && (iscellstr(names) || ischar(names) || isstring(names)); %#ok<ISCLSTR>
            if usable, names = cellstr(names); end
            if ~usable || ~strcmp(names{1}, 'E')
                error('SRNNCellTypePairs:NotTypeE', ...
                    ['tau_a_E aliases tau_a{1}, which requires cell_type_names{1} ' ...
                    'to be ''E''. Use tau_a directly.']);
            end
        end

        function M = expand_block(obj, rel)
            % EXPAND_BLOCK Accept a C x C block matrix or a 1 x C presynaptic row.
            % A row is broadcast DOWN the columns, i.e. the value depends only on
            % the presynaptic type -- the classical Harris setting.
            C = obj.n_cellTypes;
            if isempty(rel)
                M = [];
            elseif isvector(rel) && numel(rel) == C
                M = repmat(reshape(rel, 1, []), C, 1);
            else
                M = rel;
            end
        end

        function val = block_sigma_s_sq(obj)
            % Sparse-effective block variances, Harris Eq. (16) blockwise.
            mg = obj.mu_tilde;
            sg = obj.sigma_tilde;
            val = obj.alpha * (1 - obj.alpha) * mg.^2 + obj.alpha * sg.^2;
        end

        function h = build_activation(obj, which_one)
            % BUILD_ACTIVATION Handle for the chosen nonlinearity (1=fn, 2=deriv).
            % Parameters are captured BY VALUE and the handle is rebuilt on every
            % access, so it can never disagree with S_a / S_c.
            if ~isempty(obj.activation_custom)
                SRNNCellTypePairs.check_activation_custom(obj.activation_custom);
                h = obj.activation_custom{which_one};
                return;
            end
            a = obj.S_a;
            % A per-neuron setpoint vector, once build_setpoints() has drawn
            % one, otherwise the shared scalar. In the vector case the handle
            % is only valid for length-n input, since c lines up with x
            % element by element.
            c = SRNNCellTypePairs.pick(obj.S_c_vec, obj.S_c);
            switch obj.activation
                case 'logistic'
                    if which_one == 1
                        h = @(x) SRNNCellTypePairs.logisticSigmoid(x, c);
                    else
                        h = @(x) SRNNCellTypePairs.logisticSigmoidDerivative(x, c);
                    end
                case 'piecewise'
                    if which_one == 1
                        h = @(x) SRNNCellTypePairs.piecewiseSigmoid(x, a, c);
                    else
                        h = @(x) SRNNCellTypePairs.piecewiseSigmoidDerivative(x, a, c);
                    end
                case 'tanh'
                    if which_one == 1
                        h = @SRNNCellTypePairs.tanhActivation;
                    else
                        h = @SRNNCellTypePairs.tanhActivationDerivative;
                    end
                otherwise
                    error('SRNNCellTypePairs:InvalidParams', ...
                        'Unknown activation ''%s''. Valid: %s.', ...
                        char(string(obj.activation)), ...
                        strjoin(SRNNCellTypePairs.activation_names(), ', '));
            end
        end

        function set_defaults(obj)
            % The activation needs no setup: `activation` defaults to 'logistic'
            % as a plain property and the handles are Dependent.
            obj.input_config = struct( ...
                'n_steps', 3, ...
                'step_density', struct(), ...
                'amp', 0.5, ...
                'no_stim_pattern', logical([1 0 1]), ...
                'intrinsic_drive', [], ...
                'positive_only', false);
            obj.T_plot = [];
        end

        function complete_type_defaults(obj)
            C = obj.n_cellTypes;
            if isempty(C) || ~isscalar(C) || C < 1 || C ~= round(C)
                return;
            end
            obj.cell_type_names = cellstr(obj.cell_type_names);
            obj.f = reshape(obj.f, 1, []);
            % The _relative tildes are left in whatever shape they were given:
            % expand_block accepts a 1 x C presynaptic row or a full C x C block.
            if isempty(obj.n_a), obj.n_a = zeros(1, C); end
            if isempty(obj.c), obj.c = repmat(0.15 / 3, 1, C); end
            obj.n_a = reshape(obj.n_a, 1, []);
            obj.c = reshape(obj.c, 1, []);

            % Setpoint statistics: a scalar means "the same for every type".
            % mu_S_c stays empty when unset -- that is what selects the shared
            % scalar S_c -- while sigma_S_c fills in as zeros.
            if isscalar(obj.mu_S_c), obj.mu_S_c = repmat(obj.mu_S_c, 1, C); end
            if isempty(obj.sigma_S_c)
                obj.sigma_S_c = zeros(1, C);
            elseif isscalar(obj.sigma_S_c)
                obj.sigma_S_c = repmat(obj.sigma_S_c, 1, C);
            end
            obj.mu_S_c = reshape(obj.mu_S_c, 1, []);
            obj.sigma_S_c = reshape(obj.sigma_S_c, 1, []);

            if isempty(obj.tau_a), obj.tau_a = cell(1, C); end
            if iscell(obj.tau_a), obj.tau_a = reshape(obj.tau_a, 1, []); end
            if iscell(obj.tau_a) && numel(obj.tau_a) == C
                for q = 1:C
                    if obj.n_a(q) > 0 && isempty(obj.tau_a{q})
                        obj.tau_a{q} = logspace(log10(0.25), log10(10), obj.n_a(q));
                    elseif ~isempty(obj.tau_a{q})
                        obj.tau_a{q} = reshape(obj.tau_a{q}, 1, []);
                    end
                end
            end
            if ~isfield(obj.input_config, 'step_density') || ...
                    isempty(fieldnames(obj.input_config.step_density))
                density = struct();
                for q = 1:C
                    density.(obj.cell_type_names{q}) = 0.2;
                end
                obj.input_config.step_density = density;
            end
        end

        function validate(obj)
            C = obj.n_cellTypes;
            if ~isscalar(C) || ~isfinite(C) || C < 1 || C ~= round(C)
                error('SRNNCellTypePairs:InvalidParams', 'n_cellTypes must be a positive integer.');
            end
            if numel(obj.cell_type_names) ~= C || ...
                    any(~cellfun(@isvarname, obj.cell_type_names)) || ...
                    numel(unique(obj.cell_type_names)) ~= C
                error('SRNNCellTypePairs:InvalidParams', ...
                    'cell_type_names must contain unique valid MATLAB identifiers.');
            end
            if ~isscalar(obj.n) || obj.n < C || obj.n ~= round(obj.n)
                error('SRNNCellTypePairs:InvalidParams', ...
                    'n must be an integer at least as large as n_cellTypes.');
            end
            if ~isscalar(obj.indegree) || obj.indegree <= 0 || obj.indegree > obj.n
                error('SRNNCellTypePairs:InvalidParams', 'indegree must satisfy 0 < indegree <= n.');
            end
            if numel(obj.f) ~= C || any(~isfinite(obj.f)) || any(obj.f <= 0) || ...
                    abs(sum(obj.f) - 1) > 1e-12
                error('SRNNCellTypePairs:InvalidParams', ...
                    'f must contain one positive fraction per type and sum to 1.');
            end
            RMTCellTypes.allocate_counts(obj.n, obj.f);
            mu_b = obj.expand_block(obj.mu_tilde_relative);
            sg_b = obj.expand_block(obj.sigma_tilde_relative);
            if ~isequal(size(mu_b), [C C]) || any(~isfinite(mu_b(:))) || ...
                    ~isequal(size(sg_b), [C C]) || any(~isfinite(sg_b(:))) || ...
                    any(sg_b(:) < 0)
                error('SRNNCellTypePairs:InvalidParams', ...
                    ['mu_tilde_relative and sigma_tilde_relative must each be a ' ...
                    '%d x %d block matrix (post, pre) or a 1 x %d presynaptic ' ...
                    'row, finite, with sigma nonnegative.'], C, C, C);
            end

            % Activation choice, so a typo fails here rather than at the first
            % phi() evaluation inside the solver.
            if ~isempty(obj.activation_custom)
                SRNNCellTypePairs.check_activation_custom(obj.activation_custom);
            elseif ~ismember(obj.activation, SRNNCellTypePairs.activation_names())
                error('SRNNCellTypePairs:InvalidParams', ...
                    'Unknown activation ''%s''. Valid: %s.', ...
                    char(string(obj.activation)), ...
                    strjoin(SRNNCellTypePairs.activation_names(), ', '));
            end

            % Per-type setpoint statistics.
            if ~isempty(obj.mu_S_c) && (numel(obj.mu_S_c) ~= C || ...
                    ~isnumeric(obj.mu_S_c) || any(~isfinite(obj.mu_S_c)))
                error('SRNNCellTypePairs:InvalidParams', ...
                    ['mu_S_c must be empty, a scalar, or a finite 1 x %d row ' ...
                    '(one mean setpoint per cell type).'], C);
            end
            if numel(obj.sigma_S_c) ~= C || ~isnumeric(obj.sigma_S_c) || ...
                    any(~isfinite(obj.sigma_S_c)) || any(obj.sigma_S_c < 0)
                error('SRNNCellTypePairs:InvalidParams', ...
                    ['sigma_S_c must be empty, a scalar, or a finite nonnegative ' ...
                    '1 x %d row (one standard deviation per cell type).'], C);
            end
            if obj.has_setpoint_heterogeneity()
                if ~isempty(obj.activation_custom)
                    error('SRNNCellTypePairs:InvalidParams', ...
                        ['Per-type setpoints (mu_S_c / sigma_S_c) cannot be applied to ' ...
                        'activation_custom: a custom handle takes only x, so the model ' ...
                        'has no way to inject a per-neuron centre. Build the ' ...
                        'heterogeneity into the custom handle itself, or clear it.']);
                end
                if strcmp(obj.activation, 'tanh')
                    error('SRNNCellTypePairs:InvalidParams', ...
                        ['Per-type setpoints (mu_S_c / sigma_S_c) require a nonlinearity ' ...
                        'with a centre; ''tanh'' has none. Use ''logistic'' or ' ...
                        '''piecewise''.']);
                end
            end

            % Reference network for F, when it is pinned.
            if ~obj.F_tracks_network
                if ~isscalar(obj.F_ref_n) || ~isscalar(obj.F_ref_indegree) || ...
                        obj.F_ref_n <= 0 || obj.F_ref_indegree <= 0 || ...
                        obj.F_ref_indegree > obj.F_ref_n
                    error('SRNNCellTypePairs:InvalidParams', ...
                        ['F_ref_n and F_ref_indegree must be positive scalars with ' ...
                        'F_ref_indegree <= F_ref_n.']);
                end
            end

            vector_fields = {'n_a', 'c'};
            for k = 1:numel(vector_fields)
                name = vector_fields{k};
                if ~isnumeric(obj.(name)) || numel(obj.(name)) ~= C || ...
                        any(~isfinite(obj.(name)))
                    error('SRNNCellTypePairs:InvalidParams', ...
                        '%s must have one finite numeric value per cell type.', name);
                end
            end
            if any(obj.n_a < 0 | obj.n_a ~= round(obj.n_a))
                error('SRNNCellTypePairs:InvalidParams', ...
                    'n_a must contain nonnegative integers.');
            end
            if any(obj.c < 0)
                error('SRNNCellTypePairs:InvalidParams', ...
                    'c must be nonnegative.');
            end
            if ~iscell(obj.tau_a) || numel(obj.tau_a) ~= C
                error('SRNNCellTypePairs:InvalidParams', ...
                    'tau_a must be a 1-by-n_cellTypes cell array.');
            end
            for q = 1:C
                if numel(obj.tau_a{q}) ~= obj.n_a(q) || ...
                        any(~isfinite(obj.tau_a{q})) || any(obj.tau_a{q} <= 0)
                    error('SRNNCellTypePairs:InvalidParams', ...
                        'tau_a{%d} must contain n_a(%d) positive values.', q, q);
                end
            end
            obj.compile_synapse_config();
            if obj.T_range(2) <= obj.T_range(1) || obj.fs <= 0 || obj.tau_d <= 0
                error('SRNNCellTypePairs:InvalidParams', ...
                    'T_range, fs, and tau_d must define positive simulation intervals.');
            end
            has_generator = isfield(obj.input_config, 'generator') && ...
                isa(obj.input_config.generator, 'function_handle');
            if ~has_generator
                if ~isfield(obj.input_config, 'step_density') || ...
                        ~isstruct(obj.input_config.step_density)
                    error('SRNNCellTypePairs:InvalidParams', ...
                        'input_config.step_density must be a named struct.');
                end
                density_fields = sort(fieldnames(obj.input_config.step_density));
                if ~isequal(density_fields, sort(obj.cell_type_names(:)))
                    error('SRNNCellTypePairs:InvalidParams', ...
                        'step_density fields must exactly match cell_type_names.');
                end
                for q = 1:C
                    d = obj.input_config.step_density.(obj.cell_type_names{q});
                    if ~isscalar(d) || ~isfinite(d) || d < 0 || d > 1
                        error('SRNNCellTypePairs:InvalidParams', ...
                            'Each step density must satisfy 0 <= density <= 1.');
                    end
                end
            end
        end

        function pair = compile_synapse_config(obj)
            %COMPILE_SYNAPSE_CONFIG Validate named routes and expand scalars.
            C = obj.n_cellTypes;
            pair = struct( ...
                'n_b', zeros(C, C), 'n_g', zeros(C, C), ...
                'tau_b_rec', {cell(C, C)}, 'tau_b_rel', {cell(C, C)}, ...
                'tau_g_dec', {cell(C, C)}, 'tau_g_fac', {cell(C, C)}, ...
                'G', {cell(C, C)});
            config = obj.synapse_config;
            if isempty(config), return; end
            if ~isstruct(config) || ~isscalar(config)
                error('SRNNCellTypePairs:InvalidSynapseConfig', ...
                    'synapse_config must be a scalar named struct.');
            end

            names = obj.cell_type_names;
            unknown_pre = setdiff(fieldnames(config), names);
            if ~isempty(unknown_pre)
                error('SRNNCellTypePairs:InvalidSynapseConfig', ...
                    'Unknown presynaptic type: %s.', unknown_pre{1});
            end
            for pre = 1:C
                pre_name = names{pre};
                if ~isfield(config, pre_name), continue; end
                targets = config.(pre_name);
                if ~isstruct(targets) || ~isscalar(targets)
                    error('SRNNCellTypePairs:InvalidSynapseConfig', ...
                        'synapse_config.%s must be a scalar struct.', pre_name);
                end
                unknown_post = setdiff(fieldnames(targets), names);
                if ~isempty(unknown_post)
                    error('SRNNCellTypePairs:InvalidSynapseConfig', ...
                        'Unknown postsynaptic type on %s route: %s.', ...
                        pre_name, unknown_post{1});
                end
                for post = 1:C
                    post_name = names{post};
                    if ~isfield(targets, post_name), continue; end
                    route = targets.(post_name);
                    if isempty(route), continue; end
                    if ~isstruct(route) || ~isscalar(route)
                        error('SRNNCellTypePairs:InvalidSynapseConfig', ...
                            'Route %s->%s must be a scalar struct.', pre_name, post_name);
                    end
                    unknown_mechanism = setdiff(fieldnames(route), {'std', 'stf'});
                    if ~isempty(unknown_mechanism)
                        error('SRNNCellTypePairs:InvalidSynapseConfig', ...
                            'Unknown field on route %s->%s: %s.', ...
                            pre_name, post_name, unknown_mechanism{1});
                    end

                    if isfield(route, 'std') && ~isempty(route.std) && ...
                            ~(isstruct(route.std) && isempty(fieldnames(route.std)))
                        std_config = route.std;
                        SRNNCellTypePairs.validate_fields(std_config, ...
                            {'tau_rec', 'tau_rel'}, pre_name, post_name, 'std');
                        tau_rec = SRNNCellTypePairs.positive_row(std_config.tau_rec, ...
                            'tau_rec', pre_name, post_name);
                        if isempty(tau_rec)
                            error('SRNNCellTypePairs:InvalidSynapseConfig', ...
                                'STD tau_rec must be nonempty on route %s->%s.', ...
                                pre_name, post_name);
                        end
                        tau_rel = SRNNCellTypePairs.expand_pair_parameter( ...
                            std_config.tau_rel, numel(tau_rec), 'tau_rel', ...
                            pre_name, post_name, false);
                        pair.n_b(pre, post) = numel(tau_rec);
                        pair.tau_b_rec{pre, post} = tau_rec;
                        pair.tau_b_rel{pre, post} = tau_rel;
                    end

                    if isfield(route, 'stf') && ~isempty(route.stf) && ...
                            ~(isstruct(route.stf) && isempty(fieldnames(route.stf)))
                        stf_config = route.stf;
                        SRNNCellTypePairs.validate_fields(stf_config, ...
                            {'tau_dec', 'tau_fac', 'G'}, pre_name, post_name, 'stf');
                        tau_dec = SRNNCellTypePairs.positive_row(stf_config.tau_dec, ...
                            'tau_dec', pre_name, post_name);
                        if isempty(tau_dec)
                            error('SRNNCellTypePairs:InvalidSynapseConfig', ...
                                'STF tau_dec must be nonempty on route %s->%s.', ...
                                pre_name, post_name);
                        end
                        tau_fac = SRNNCellTypePairs.expand_pair_parameter( ...
                            stf_config.tau_fac, numel(tau_dec), 'tau_fac', ...
                            pre_name, post_name, false);
                        G = SRNNCellTypePairs.expand_pair_parameter( ...
                            stf_config.G, numel(tau_dec), 'G', ...
                            pre_name, post_name, true);
                        pair.n_g(pre, post) = numel(tau_dec);
                        pair.tau_g_dec{pre, post} = tau_dec;
                        pair.tau_g_fac{pre, post} = tau_fac;
                        pair.G{pre, post} = G;
                    end
                end
            end
        end

        function build_network(obj)
            rng(obj.rng_seeds(1));
            % RMTBlocks rather than RMTCellTypes: it indexes the statistics by
            % BOTH postsynaptic (row) and presynaptic (column) type, so E->E can
            % differ from E->I, and it supplies R and lambda_O. RMTCellTypes
            % applied the tildes raw with no normalization and had neither.
            % g_mu / g_sigma stay at 1 -- level_of_chaos is applied to the
            % assembled W below, because rescale_by_abscissa divides by the
            % abscissa BEFORE that multiply.
            rmt = RMTBlocks(obj.n);
            rmt.alpha = obj.alpha;
            rmt.f = obj.f;
            rmt.mu_tilde = obj.mu_tilde;
            rmt.sigma_tilde = obj.sigma_tilde;
            rmt.zrs_mode = obj.zrs_mode;
            % RMTBlocks returns a DENSE matrix where RMTCellTypes returned a
            % sparse one. Restore sparse storage: this class is built for sparse
            % connectivity and its state indexing assumes it.
            W_raw = sparse(rmt.W);
            if obj.rescale_by_abscissa
                abscissa = max(real(eig(full(W_raw))));
                if abs(abscissa) <= eps
                    error('SRNNCellTypePairs:InvalidConnectivity', ...
                        'Cannot rescale a matrix with zero spectral abscissa.');
                end
                W_raw = W_raw / abscissa;
            end
            obj.W = obj.level_of_chaos .* W_raw;
            pair = obj.compile_synapse_config();
            obj.dead_state_count = 0;
            for pre = 1:obj.n_cellTypes
                pre_idx = obj.type_indices{pre};
                for post = 1:obj.n_cellTypes
                    n_route_states = pair.n_b(pre, post) + pair.n_g(pre, post);
                    if n_route_states == 0, continue; end
                    post_idx = obj.type_indices{post};
                    has_output = any(obj.W(post_idx, pre_idx) ~= 0, 1);
                    obj.dead_state_count = obj.dead_state_count + ...
                        sum(~has_output) * n_route_states;
                end
            end
            eig_W = eig(full(obj.W));
            fprintf('W created: spectral radius = %.3f, abscissa = %.3f\n', ...
                max(abs(eig_W)), max(real(eig_W)));
            fprintf('Pair-specific dead-end states: %d\n', obj.dead_state_count);

            % Per-neuron nonlinearity setpoints. Drawn here, after W, so the
            % vector exists before build_stimulus() and get_params().
            obj.build_setpoints();
        end

        function tf = has_setpoint_heterogeneity(obj)
            % HAS_SETPOINT_HETEROGENEITY Whether a per-neuron S_c was requested.
            tf = ~isempty(obj.mu_S_c) || any(obj.sigma_S_c > 0);
        end

        function build_setpoints(obj)
            % BUILD_SETPOINTS Draw the per-neuron nonlinearity setpoints S_c_vec
            %
            % S_c_i = mu_S_c(q) + sigma_S_c(q) * randn for each neuron of type
            % q, with an empty mu_S_c falling back to the shared scalar S_c.
            % Leaves S_c_vec EMPTY when no heterogeneity was requested, which
            % keeps every downstream code path on the scalar branch and the
            % results bit-identical to a model without this feature.

            obj.S_c_vec = [];
            if ~obj.has_setpoint_heterogeneity()
                return;
            end

            % Own RNG substream: the state is saved and restored so W, the
            % stimulus and x0 are bit-identical whether or not setpoints are
            % drawn (initialize_state draws x0 from whatever state build left
            % behind). The offset keeps the draw off the stream that produced
            % W's first column, which seeding with rng_seeds(1) would not do.
            seed = obj.S_c_seed;
            if isempty(seed)
                seed = obj.rng_seeds(1) + 104729;
            end
            stream_state = rng;
            rng(seed);
            z = randn(obj.n, 1);
            rng(stream_state);

            C = obj.n_cellTypes;
            idx = obj.type_indices;
            mu = SRNNCellTypePairs.pick(obj.mu_S_c, repmat(obj.S_c, 1, C));
            vals = zeros(obj.n, 1);
            parts = cell(1, C);
            for q = 1:C
                vals(idx{q}) = mu(q) + obj.sigma_S_c(q) * z(idx{q});
                parts{q} = sprintf('%s %.4f +/- %.4f', ...
                    obj.cell_type_names{q}, mu(q), obj.sigma_S_c(q));
            end
            obj.S_c_vec = vals;

            fprintf('Per-neuron S_c drawn (seed %g): %s\n', seed, strjoin(parts, ', '));
        end

        function build_stimulus(obj)
            obj.generate_stimulus();
            obj.u_interpolant = griddedInterpolant( ...
                obj.t_ex, obj.u_ex', 'linear', 'none');
            obj.S0 = SRNNCellTypePairs.initialize_state(obj.get_params());
        end

        function generate_stimulus(obj)
            dt = 1 / obj.fs;
            if isempty(obj.input_config.intrinsic_drive)
                obj.input_config.intrinsic_drive = zeros(obj.n, 1);
            end
            params = obj.get_params();
            [u_stim, t_stim] = SRNNCellTypePairs.generate_external_input( ...
                params, obj.T_range(2), obj.fs, obj.rng_seeds(2), obj.input_config);
            if obj.T_range(1) < 0
                t_pre = (obj.T_range(1):dt:-dt)';
                obj.t_ex = [t_pre; t_stim];
                obj.u_ex = [zeros(obj.n, numel(t_pre)), u_stim];
            else
                keep = t_stim >= obj.T_range(1);
                obj.t_ex = t_stim(keep);
                obj.u_ex = u_stim(:, keep);
            end
            obj.u_ex = obj.u_ex .* obj.u_ex_scale;
        end

        function decimate_and_unpack(obj)
            sample_indices = 1:obj.plot_deci:numel(obj.t_out);
            t_plot = obj.t_out(sample_indices);
            S_plot = obj.S_out(sample_indices, :);
            [x, a, b, r, synaptic_output, g] = ...
                SRNNCellTypePairs.unpack_and_compute_states( ...
                S_plot, obj.cached_params);
            u = SRNNCellTypePairs.split_by_type( ...
                obj.u_ex(:, sample_indices), obj.cached_params, 1);
            range = obj.T_plot;
            if isempty(range), range = obj.T_range; end
            keep = t_plot >= range(1) & t_plot <= range(2);
            obj.plot_data = struct( ...
                't', t_plot(keep), ...
                'u', SRNNCellTypePairs.trim_named(u, 2, keep), ...
                'x', SRNNCellTypePairs.trim_named(x, 2, keep), ...
                'r', SRNNCellTypePairs.trim_named(r, 2, keep), ...
                'a', SRNNCellTypePairs.trim_named(a, 3, keep), ...
                'b', SRNNCellTypePairs.trim_routes(b, 3, keep), ...
                'g', SRNNCellTypePairs.trim_routes(g, 3, keep), ...
                'synaptic_output', ...
                    SRNNCellTypePairs.trim_routes(synaptic_output, 2, keep));
        end
    end

    %% Plotting and spectrum inspection
    methods
        function [fig_handle, ax_handles] = plot(obj, varargin)
            if ~obj.has_run || isempty(obj.plot_data)
                error('SRNNCellTypePairs:NotRun', ...
                    'Run the model with store_decimated_state=true before plotting.');
            end
            range = obj.T_plot;
            for k = 1:2:numel(varargin)
                if strcmpi(varargin{k}, 'T_plot'), range = varargin{k + 1}; end
            end
            if isempty(range), range = obj.T_range; end

            p = obj.plot_data;
            has_a = any(obj.n_a > 0);
            has_b = any(obj.n_b_pairs(:) > 0);
            has_g = any(obj.n_g_pairs(:) > 0);
            has_lya = ~strcmpi(obj.lya_method, 'none') && ~isempty(obj.lya_results);
            n_panels = 4 + has_a + has_b + has_g + has_lya;
            fig_handle = figure('Name', 'SRNNCellTypePairs time series');
            tiledlayout(n_panels, 1);
            ax_handles = gobjects(n_panels, 1);
            panel = 0;

            panel = panel + 1; ax_handles(panel) = nexttile;
            SRNNCellTypePairs.plot_named_series( ...
                p.t, p.u, obj.cell_type_names, 'External input', true);
            panel = panel + 1; ax_handles(panel) = nexttile;
            SRNNCellTypePairs.plot_named_series( ...
                p.t, p.x, obj.cell_type_names, 'Dendritic state', false);
            panel = panel + 1; ax_handles(panel) = nexttile;
            SRNNCellTypePairs.plot_named_series( ...
                p.t, p.r, obj.cell_type_names, 'Firing rate', false);
            panel = panel + 1; ax_handles(panel) = nexttile;
            SRNNCellTypePairs.plot_route_series( ...
                p.t, p.synaptic_output, obj.cell_type_names, ...
                'Synaptic output', true);
            if has_a
                panel = panel + 1; ax_handles(panel) = nexttile;
                a_collapsed = SRNNCellTypePairs.collapse_named(p.a, obj.cell_type_names, 'sum');
                SRNNCellTypePairs.plot_named_series( ...
                    p.t, a_collapsed, obj.cell_type_names, 'SFA (sum)', false);
            end
            if has_b
                panel = panel + 1; ax_handles(panel) = nexttile;
                b_collapsed = SRNNCellTypePairs.collapse_routes(p.b, 'prod');
                SRNNCellTypePairs.plot_route_series( ...
                    p.t, b_collapsed, obj.cell_type_names, 'STD (product)', true);
            end
            if has_g
                panel = panel + 1; ax_handles(panel) = nexttile;
                g_collapsed = SRNNCellTypePairs.collapse_routes(p.g, 'prod');
                SRNNCellTypePairs.plot_route_series( ...
                    p.t, g_collapsed, obj.cell_type_names, 'STF (product)', true);
            end
            if has_lya
                panel = panel + 1; ax_handles(panel) = nexttile;
                if isfield(obj.lya_results, 'local_lya')
                    plot(obj.lya_results.t_lya, obj.lya_results.local_lya, ...
                        'Color', lines(1));
                    hold on;
                    yline(0, '--k');
                    ylabel('\lambda_1');
                    xlim(range);
                    if isfield(obj.lya_results, 'LLE')
                        y_limits = ylim;
                        text(range(2), y_limits(1) + 0.05 * diff(y_limits), ...
                            ['$\lambda_1 = ' sprintf('%.2f', obj.lya_results.LLE) '$'], ...
                            'HorizontalAlignment', 'right', ...
                            'VerticalAlignment', 'bottom', ...
                            'Interpreter', 'latex');
                    end
                    hold off;
                elseif isfield(obj.lya_results, 'local_LE_spectrum_t')
                    plot(obj.lya_results.t_lya, obj.lya_results.local_LE_spectrum_t);
                    ylabel('Local LE');
                end
            end
            linkaxes(ax_handles, 'x');
            xlim(ax_handles(end), range);
            xlabel(ax_handles(end), 'Time (s)');
        end

        function [fig_handle, ax_handles] = plot_celltypes(obj, varargin)
            if ~obj.has_run || isempty(obj.plot_data)
                error('SRNNCellTypePairs:NotRun', ...
                    'Run the model with store_decimated_state=true before plotting.');
            end

            range = obj.T_plot;
            for k = 1:2:numel(varargin)
                if strcmpi(varargin{k}, 'T_plot'), range = varargin{k + 1}; end
            end
            if isempty(range), range = obj.T_range; end

            p = obj.plot_data;
            has_a = any(obj.n_a > 0);
            has_b = any(obj.n_b_pairs(:) > 0);
            has_g = any(obj.n_g_pairs(:) > 0);
            has_lya = ~strcmpi(obj.lya_method, 'none') && ~isempty(obj.lya_results);
            n_rows = 4 + has_a + has_b + has_g + has_lya;
            C = obj.n_cellTypes;
            colors = SRNNCellTypePairs.type_colors(C);
            b_collapsed = SRNNCellTypePairs.collapse_routes(p.b, 'prod');
            g_collapsed = SRNNCellTypePairs.collapse_routes(p.g, 'prod');
            fig_handle = figure('Name', 'SRNNCellTypePairs by presynaptic type');
            layout = tiledlayout(fig_handle, n_rows, C, ...
                'TileSpacing', 'compact', 'Padding', 'compact');
            ax_handles = gobjects(n_rows, C);

            for q = 1:C
                name = obj.cell_type_names{q};
                row = 0;
                basic = {p.u.(name), p.x.(name), p.r.(name)};
                basic_labels = {'External input', 'Dendritic state', 'Firing rate'};
                for k = 1:numel(basic)
                    row = row + 1;
                    tile_index = (row - 1) * C + q;
                    ax = nexttile(layout, tile_index);
                    ax_handles(row, q) = ax;
                    SRNNCellTypePairs.plot_celltype_lines(ax, p.t, basic{k}, colors(q, :));
                    if q == 1, ylabel(ax, basic_labels{k}); end
                    if row == 1, title(ax, name, 'Interpreter', 'none'); end
                end
                row = row + 1;
                ax = nexttile(layout, (row - 1) * C + q);
                ax_handles(row, q) = ax;
                SRNNCellTypePairs.plot_postsynaptic_route_lines( ...
                    ax, p.t, p.synaptic_output.(name), obj.cell_type_names);
                if q == 1, ylabel(ax, 'Synaptic output'); end
                if has_a
                    row = row + 1;
                    ax = nexttile(layout, (row - 1) * C + q);
                    ax_handles(row, q) = ax;
                    if isempty(p.a.(name))
                        SRNNCellTypePairs.show_empty_axis(ax, 'No SFA');
                    else
                        values = reshape(sum(p.a.(name), 2), ...
                            obj.n_per_type(q), []);
                        SRNNCellTypePairs.plot_celltype_lines(ax, p.t, values, colors(q, :));
                    end
                    if q == 1, ylabel(ax, 'SFA (sum)'); end
                end
                if has_b
                    row = row + 1;
                    ax = nexttile(layout, (row - 1) * C + q);
                    ax_handles(row, q) = ax;
                    SRNNCellTypePairs.plot_postsynaptic_route_lines( ...
                        ax, p.t, b_collapsed.(name), obj.cell_type_names);
                    if q == 1, ylabel(ax, 'STD (product)'); end
                end
                if has_g
                    row = row + 1;
                    ax = nexttile(layout, (row - 1) * C + q);
                    ax_handles(row, q) = ax;
                    SRNNCellTypePairs.plot_postsynaptic_route_lines( ...
                        ax, p.t, g_collapsed.(name), obj.cell_type_names);
                    if q == 1, ylabel(ax, 'STF (product)'); end
                end
                if has_lya
                    row = row + 1;
                    ax = nexttile(layout, (row - 1) * C + q);
                    ax_handles(row, q) = ax;
                    if isfield(obj.lya_results, 'local_lya')
                        plot(ax, obj.lya_results.t_lya, obj.lya_results.local_lya, 'k');
                        if q == 1, ylabel(ax, 'Local LLE'); end
                    else
                        plot(ax, obj.lya_results.t_lya, ...
                            obj.lya_results.local_LE_spectrum_t);
                        if q == 1, ylabel(ax, 'Local LE'); end
                    end
                end
            end

            linkaxes(ax_handles(:), 'x');
            xlim(ax_handles(1, 1), range);
            for q = 1:C
                xlabel(ax_handles(end, q), 'Time (s)');
            end
        end

        function [fig_handle, ax_handles] = plot_eigenvalues(obj, times_sec)
            if isempty(obj.S_out)
                error('SRNNCellTypePairs:NoStateData', ...
                    'Set store_full_state=true before running.');
            end
            sample_indices = round((times_sec - obj.t_out(1)) * obj.fs) + 1;
            sample_indices = unique(max(1, min(sample_indices, size(obj.S_out, 1))));
            n_plots = numel(sample_indices);
            fig_handle = figure('Name', 'SRNNCellTypePairs Jacobian eigenvalues');
            ax_handles = gobjects(n_plots, 1);
            n_cols = ceil(sqrt(n_plots));
            n_rows = ceil(n_plots / n_cols);
            for k = 1:n_plots
                ax_handles(k) = subplot(n_rows, n_cols, k);
                J = SRNNCellTypePairs.compute_Jacobian_fast( ...
                    obj.S_out(sample_indices(k), :)', obj.cached_params);
                ev = eig(full(J));
                plot(real(ev), imag(ev), 'ko', 'MarkerSize', 3);
                xline(0, 'r--'); yline(0, 'k-'); axis equal; grid on;
                title(sprintf('t = %.3g s', obj.t_out(sample_indices(k))));
                xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
            end
        end

        function [fig_handle, ax_handles] = plot_W_spectrum(obj)
            if ~obj.is_built
                error('SRNNCellTypePairs:NotBuilt', 'Call build() first.');
            end
            A = -speye(obj.n) + obj.W;
            spectra = {eig(full(A)), eig(full(A / obj.tau_d))};
            titles = {'-I + W', '(-I + W)/\tau_d'};
            fig_handle = figure('Name', 'Empirical connectivity spectra');
            ax_handles = gobjects(2, 1);
            for k = 1:2
                ax_handles(k) = subplot(1, 2, k);
                plot(real(spectra{k}), imag(spectra{k}), 'ko', 'MarkerSize', 3);
                xline(0, 'r--'); yline(0, 'k-'); axis equal; grid on;
                xlabel('Re(\lambda)'); ylabel('Im(\lambda)'); title(titles{k});
            end
        end

        function [fig_handle, ax_handle] = plot_W(obj, ax, clim_val)
            % PLOT_W Image the connectivity matrix with a diverging red/white/blue
            % colormap and a symmetric colour range, so zero is white and the
            % sign of every synapse is readable at a glance.
            %
            % This plots obj.W -- the SCALED matrix, after level_of_chaos and
            % rescale_by_abscissa -- i.e. the network actually simulated, and the
            % one R and lambda_O describe. RMTBlocks.plot_W shows the generator's
            % unscaled output instead, which differs by those factors.
            %
            % Cell-type boundaries are drawn over the image, since the point of
            % this class is that the (post, pre) blocks can differ: each block of
            % the picture is one route.
            %
            %   plot_W()             - new figure, clim from the off-diagonal max
            %   plot_W(ax)           - draw into ax
            %   plot_W(ax, clim_val) - draw into ax with a shared clim, for
            %                          comparing two networks on equal footing
            if ~obj.is_built
                error('SRNNCellTypePairs:NotBuilt', 'Call build() first.');
            end
            if nargin < 2 || isempty(ax)
                fig_handle = figure('Name', 'SRNNCellTypePairs W');
                ax = axes(fig_handle);
            else
                fig_handle = ancestor(ax, 'figure');
            end
            ax_handle = ax;

            W_plot = full(obj.W);

            % Auto range from the OFF-DIAGONAL entries: a diagonal shift, if any,
            % is rendered as-is and simply saturates rather than flattening the
            % rest of the picture.
            if nargin < 3 || isempty(clim_val)
                W_no_diag = W_plot;
                W_no_diag(1:obj.n+1:end) = 0;
                clim_val = ceil(max(abs(W_no_diag(:))) * 10) / 10;
                if clim_val == 0 || ~isfinite(clim_val)
                    clim_val = 0.1;
                end
            end

            imagesc(ax, W_plot);
            colormap(ax, redwhiteblue_colormap(256));
            clim(ax, [-clim_val, clim_val]);
            axis(ax, 'square');
            colorbar(ax);

            % Block boundaries and type labels at the centre of each band.
            counts = obj.n_per_type;
            edges = cumsum(counts);
            hold(ax, 'on');
            for k = 1:numel(edges) - 1
                e = edges(k) + 0.5;
                plot(ax, [e e], [0.5, obj.n + 0.5], 'k-', 'LineWidth', 0.75);
                plot(ax, [0.5, obj.n + 0.5], [e e], 'k-', 'LineWidth', 0.75);
            end
            hold(ax, 'off');
            centres = cumsum(counts) - counts / 2 + 0.5;
            set(ax, 'XTick', centres, 'XTickLabel', obj.cell_type_names, ...
                    'YTick', centres, 'YTickLabel', obj.cell_type_names, ...
                    'TickLength', [0 0]);
            xlabel(ax, 'presynaptic type');
            ylabel(ax, 'postsynaptic type');
            title(ax, sprintf('W  (level\\_of\\_chaos = %g,  R = %.3f)', ...
                obj.level_of_chaos, obj.R));
        end

        function [fig_handle, ax_handle] = plot_weight_histogram(obj, bin_edges)
            if ~obj.is_built
                error('SRNNCellTypePairs:NotBuilt', 'Call build() first.');
            end
            values = nonzeros(obj.W);
            if nargin < 2 || isempty(bin_edges)
                radius = max(abs(values));
                if isempty(radius) || radius == 0, radius = 1; end
                bin_edges = linspace(-radius, radius, 51);
            end
            colors = SRNNCellTypePairs.type_colors(obj.n_cellTypes);
            fig_handle = figure('Name', 'Weights by presynaptic cell type');
            ax_handle = axes(fig_handle); hold(ax_handle, 'on');
            for q = 1:obj.n_cellTypes
                wq = nonzeros(obj.W(:, obj.type_indices{q}));
                histogram(ax_handle, wq, bin_edges, 'DisplayName', ...
                    obj.cell_type_names{q}, 'FaceColor', colors(q, :), ...
                    'EdgeColor', 'none', 'FaceAlpha', 0.45);
            end
            hold(ax_handle, 'off'); legend(ax_handle, 'Location', 'best');
            xlabel(ax_handle, 'Weight'); ylabel(ax_handle, 'Count');
        end
    end

    %% State layout, dynamics, and Jacobian
    methods (Static)
        function names = activation_names()
            % ACTIVATION_NAMES Valid values of the `activation` property.
            names = {'logistic', 'piecewise', 'tanh'};
        end

        function check_activation_custom(c)
            % CHECK_ACTIVATION_CUSTOM Validate the escape-hatch {fn, dfn} pair.
            ok = iscell(c) && numel(c) == 2 && ...
                all(cellfun(@(h) isa(h, 'function_handle'), c));
            if ~ok
                error('SRNNCellTypePairs:InvalidParams', ...
                    ['activation_custom must be a 1x2 cell of function handles, ' ...
                    '{phi, phi_derivative}, or empty to use the named `activation`.']);
            end
        end

        function val = scale_tilde_mat(rel, F)
            % SCALE_TILDE_MAT Elementwise rel*F, keeping exact zeros at zero.
            % A disconnected network has alpha = 0 hence F = Inf, and 0*Inf would
            % be NaN, which would then poison W.
            val = rel * F;
            val(rel == 0) = 0;
        end

        function val = pick(override, fallback)
            % PICK The override when set, else the fallback.
            if isempty(override), val = fallback; else, val = override; end
        end

        function layout = make_state_layout(n, counts, n_a, n_b, n_g)
            C = numel(counts);
            layout = struct('a', {cell(1, C)}, 'b', {cell(1, C)}, ...
                'g', {cell(1, C)}, 'x', []);
            layout.b = cell(C, C);
            layout.g = cell(C, C);
            cursor = 0;
            for q = 1:C
                len = counts(q) * n_a(q);
                layout.a{q} = cursor + (1:len);
                cursor = cursor + len;
            end
            for pre = 1:C
                for post = 1:C
                    len = counts(pre) * n_b(pre, post);
                    layout.b{pre, post} = cursor + (1:len);
                    cursor = cursor + len;
                end
            end
            for pre = 1:C
                for post = 1:C
                    len = counts(pre) * n_g(pre, post);
                    layout.g{pre, post} = cursor + (1:len);
                    cursor = cursor + len;
                end
            end
            layout.x = cursor + (1:n);
            layout.N = cursor + n;
        end

        function S0 = initialize_state(params)
            C = params.n_cellTypes;
            parts = cell(1, C + 2 * C * C + 1);
            cursor = 0;
            for q = 1:params.n_cellTypes
                cursor = cursor + 1;
                parts{cursor} = zeros(params.n_per_type(q) * params.n_a(q), 1);
            end
            for pre = 1:C
                for post = 1:C
                    cursor = cursor + 1;
                    parts{cursor} = ones(params.n_per_type(pre) * ...
                        params.n_b_pairs(pre, post), 1);
                end
            end
            for pre = 1:C
                for post = 1:C
                    cursor = cursor + 1;
                    parts{cursor} = ones(params.n_per_type(pre) * ...
                        params.n_g_pairs(pre, post), 1);
                end
            end
            parts{end} = params.x0_std .* randn(params.n, 1);
            S0 = vertcat(parts{:});
        end

        function validate_fields(config, required, pre_name, post_name, mechanism)
            if ~isstruct(config) || ~isscalar(config)
                error('SRNNCellTypePairs:InvalidSynapseConfig', ...
                    '%s configuration on route %s->%s must be a scalar struct.', ...
                    upper(mechanism), pre_name, post_name);
            end
            fields = fieldnames(config);
            missing = setdiff(required, fields);
            unknown = setdiff(fields, required);
            if ~isempty(missing)
                error('SRNNCellTypePairs:InvalidSynapseConfig', ...
                    'Missing %s field on route %s->%s: %s.', ...
                    upper(mechanism), pre_name, post_name, missing{1});
            end
            if ~isempty(unknown)
                error('SRNNCellTypePairs:InvalidSynapseConfig', ...
                    'Unknown %s field on route %s->%s: %s.', ...
                    upper(mechanism), pre_name, post_name, unknown{1});
            end
        end

        function value = positive_row(value, field_name, pre_name, post_name)
            if ~isnumeric(value) || ~isvector(value) || ...
                    any(~isfinite(value)) || any(value <= 0)
                error('SRNNCellTypePairs:InvalidSynapseConfig', ...
                    '%s on route %s->%s must contain positive finite values.', ...
                    field_name, pre_name, post_name);
            end
            value = reshape(value, 1, []);
        end

        function value = expand_pair_parameter(value, count, field_name, ...
                pre_name, post_name, at_least_one)
            if ~isnumeric(value) || ~isvector(value) || isempty(value) || ...
                    any(~isfinite(value)) || any(value <= 0)
                error('SRNNCellTypePairs:InvalidSynapseConfig', ...
                    '%s on route %s->%s must contain positive finite values.', ...
                    field_name, pre_name, post_name);
            end
            if at_least_one && any(value < 1)
                error('SRNNCellTypePairs:InvalidSynapseConfig', ...
                    'G on route %s->%s must be at least one.', pre_name, post_name);
            end
            if isscalar(value)
                value = repmat(value, 1, count);
            elseif numel(value) ~= count
                error('SRNNCellTypePairs:InvalidSynapseConfig', ...
                    '%s on route %s->%s must be scalar or have %d values.', ...
                    field_name, pre_name, post_name, count);
            else
                value = reshape(value, 1, []);
            end
        end

        function [u_ex, t_ex] = generate_external_input(params, T, fs, rng_seed, config)
            if isfield(config, 'generator') && isa(config.generator, 'function_handle')
                [u_ex, t_ex] = config.generator(params, T, fs, rng_seed, config);
                return;
            end
            rng(rng_seed);
            t_ex = (0:1/fs:T)';
            nt = numel(t_ex);
            n_steps = config.n_steps;
            if numel(config.no_stim_pattern) ~= n_steps
                error('SRNNCellTypePairs:InvalidStimulus', ...
                    'no_stim_pattern must have n_steps elements.');
            end
            if config.positive_only
                amplitudes = config.amp .* abs(randn(params.n, n_steps));
            else
                amplitudes = config.amp .* randn(params.n, n_steps);
            end
            mask = false(params.n, n_steps);
            for q = 1:params.n_cellTypes
                name = params.cell_type_names{q};
                idx = params.type_indices{q};
                mask(idx, :) = rand(numel(idx), n_steps) < config.step_density.(name);
            end
            amplitudes = amplitudes .* mask;
            amplitudes(:, config.no_stim_pattern) = 0;
            u_ex = zeros(params.n, nt);
            edges = round(linspace(1, nt + 1, n_steps + 1));
            for step = 1:n_steps
                idx = edges(step):min(edges(step + 1) - 1, nt);
                u_ex(:, idx) = repmat(amplitudes(:, step), 1, numel(idx));
            end
            drive = config.intrinsic_drive;
            if isscalar(drive), drive = repmat(drive, params.n, 1); end
            if ~isequal(size(drive), [params.n, 1])
                error('SRNNCellTypePairs:InvalidStimulus', ...
                    'intrinsic_drive must be scalar or n-by-1.');
            end
            u_ex = u_ex + drive;
        end

        function dS = dynamics_fast(t, S, params)
            u = params.u_interpolant(t)';
            C = params.n_cellTypes;
            layout = params.state_layout;
            a = cell(1, C);
            b_state = cell(C, C);
            g_state = cell(C, C);
            x = S(layout.x);
            x_eff = x;

            for q = 1:C
                idx = params.type_indices{q};
                nq = params.n_per_type(q);
                if params.n_a(q) > 0
                    a{q} = reshape(S(layout.a{q}), nq, params.n_a(q));
                    x_eff(idx) = x_eff(idx) - params.c(q) .* sum(a{q}, 2);
                else
                    a{q} = zeros(nq, 0);
                end
            end

            rate = params.activation_function(x_eff);
            a_derivatives = cell(1, C);
            for q = 1:C
                idx = params.type_indices{q};
                if params.n_a(q) > 0
                    a_derivatives{q} = ((rate(idx) - a{q}) ./ params.tau_a{q});
                else
                    a_derivatives{q} = [];
                end
            end

            b_derivatives = cell(C, C);
            g_derivatives = cell(C, C);
            recurrent = zeros(params.n, 1);
            for pre = 1:C
                pre_idx = params.type_indices{pre};
                npre = params.n_per_type(pre);
                pre_rate = rate(pre_idx);
                for post = 1:C
                    nb = params.n_b_pairs(pre, post);
                    ng = params.n_g_pairs(pre, post);
                    depression = ones(npre, 1);
                    facilitation = ones(npre, 1);
                    if nb > 0
                        b_state{pre, post} = reshape( ...
                            S(layout.b{pre, post}), npre, nb);
                        depression = prod(b_state{pre, post}, 2);
                        b_derivatives{pre, post} = ...
                            (1 - b_state{pre, post}) ./ params.tau_b_rec{pre, post} ...
                            - (b_state{pre, post} .* pre_rate) ...
                            ./ params.tau_b_rel{pre, post};
                    else
                        b_state{pre, post} = ones(npre, 0);
                        b_derivatives{pre, post} = [];
                    end
                    if ng > 0
                        g_state{pre, post} = reshape( ...
                            S(layout.g{pre, post}), npre, ng);
                        facilitation = prod(g_state{pre, post}, 2);
                        g_derivatives{pre, post} = ...
                            (1 - g_state{pre, post}) ./ params.tau_g_dec{pre, post} ...
                            + ((params.G{pre, post} - g_state{pre, post}) .* pre_rate) ...
                            ./ params.tau_g_fac{pre, post};
                    else
                        g_state{pre, post} = ones(npre, 0);
                        g_derivatives{pre, post} = [];
                    end
                    depression = SRNNCellTypePairs.apply_std_zero_floor_pair( ...
                        depression, params, pre, post);
                    post_idx = params.type_indices{post};
                    recurrent(post_idx) = recurrent(post_idx) + ...
                        params.W(post_idx, pre_idx) * ...
                        (facilitation .* depression .* pre_rate);
                end
            end

            derivatives = [a_derivatives, reshape(b_derivatives', 1, []), ...
                reshape(g_derivatives', 1, []), {(-x + recurrent + u) ./ params.tau_d}];
            for k = 1:numel(derivatives) - 1
                derivatives{k} = derivatives{k}(:);
            end
            dS = vertcat(derivatives{:});
        end

        function b_used = apply_std_zero_floor_pair(b, params, pre, post)
            b_used = b;
            if ~params.std_zero_floor || params.n_b_pairs(pre, post) == 0
                return;
            end
            p_min = prod(params.tau_b_rel{pre, post} ./ ...
                (params.tau_b_rec{pre, post} + params.tau_b_rel{pre, post}));
            b_used = (b - p_min) ./ (1 - p_min);
        end

        function J = compute_Jacobian_fast(S, params)
            C = params.n_cellTypes;
            layout = params.state_layout;
            n = params.n;
            a = cell(1, C);
            b_state = cell(C, C);
            g_state = cell(C, C);
            depression = cell(C, C);
            facilitation = cell(C, C);
            x = S(layout.x);
            x_eff = x;

            for q = 1:C
                idx = params.type_indices{q}; nq = params.n_per_type(q);
                if params.n_a(q) > 0
                    a{q} = reshape(S(layout.a{q}), nq, params.n_a(q));
                    x_eff(idx) = x_eff(idx) - params.c(q) .* sum(a{q}, 2);
                else
                    a{q} = zeros(nq, 0);
                end
            end

            rate = params.activation_function(x_eff);
            rate_prime = params.activation_function_derivative(x_eff);
            J = sparse(layout.N, layout.N);
            row_x = layout.x;
            J(row_x, row_x) = spdiags(-ones(n, 1) ./ params.tau_d, 0, n, n);

            for q = 1:C
                idx = params.type_indices{q}; nq = params.n_per_type(q);
                na = params.n_a(q);
                row_a = layout.a{q};

                if na > 0
                    tau_inv = 1 ./ params.tau_a{q}(:);
                    diag_block = kron(spdiags(-tau_inv, 0, na, na), speye(nq));
                    template = sparse(tau_inv * ones(1, na));
                    coupling = kron(template, spdiags( ...
                        -params.c(q) .* rate_prime(idx), 0, nq, nq));
                    J(row_a, row_a) = diag_block + coupling;
                    vals = kron(tau_inv, rate_prime(idx));
                    J(row_a, row_x) = sparse((1:numel(row_a))', ...
                        repmat(idx(:), na, 1), vals, numel(row_a), n);
                end
            end

            for pre = 1:C
                pre_idx = params.type_indices{pre};
                npre = params.n_per_type(pre);
                pre_rate = rate(pre_idx);
                pre_rate_prime = rate_prime(pre_idx);
                na = params.n_a(pre);
                row_a = layout.a{pre};
                for post = 1:C
                    post_idx = params.type_indices{post};
                    nb = params.n_b_pairs(pre, post);
                    ng = params.n_g_pairs(pre, post);
                    row_b = layout.b{pre, post};
                    row_g = layout.g{pre, post};

                    if nb > 0
                        b_state{pre, post} = reshape( ...
                            S(row_b), npre, nb);
                        depression{pre, post} = prod(b_state{pre, post}, 2);
                    else
                        b_state{pre, post} = ones(npre, 0);
                        depression{pre, post} = ones(npre, 1);
                    end
                    if ng > 0
                        g_state{pre, post} = reshape( ...
                            S(row_g), npre, ng);
                        facilitation{pre, post} = prod(g_state{pre, post}, 2);
                    else
                        g_state{pre, post} = ones(npre, 0);
                        facilitation{pre, post} = ones(npre, 1);
                    end
                    depression_syn = SRNNCellTypePairs.apply_std_zero_floor_pair( ...
                        depression{pre, post}, params, pre, post);

                    if nb > 0
                        if na > 0
                            coeff = b_state{pre, post} .* ...
                                (params.c(pre) .* pre_rate_prime ...
                                ./ params.tau_b_rel{pre, post});
                            stack = sparse(1:numel(row_b), ...
                                repmat((1:npre)', nb, 1), coeff(:), ...
                                numel(row_b), npre);
                            J(row_b, row_a) = kron(sparse(ones(1, na)), stack);
                        end
                        tau_rel_stack = kron(params.tau_b_rel{pre, post}(:), ...
                            ones(npre, 1));
                        diagonal = kron(-1 ./ params.tau_b_rec{pre, post}(:), ...
                            ones(npre, 1)) - repmat(pre_rate, nb, 1) ./ tau_rel_stack;
                        J(row_b, row_b) = spdiags( ...
                            diagonal, 0, numel(row_b), numel(row_b));
                        vals = -(b_state{pre, post}(:) .* ...
                            repmat(pre_rate_prime, nb, 1)) ./ tau_rel_stack;
                        J(row_b, row_x) = sparse(1:numel(row_b), ...
                            repmat(pre_idx(:), nb, 1), vals, numel(row_b), n);
                    end

                    if ng > 0
                        if na > 0
                            coeff = -(params.G{pre, post} - g_state{pre, post}) .* ...
                                (params.c(pre) .* pre_rate_prime ...
                                ./ params.tau_g_fac{pre, post});
                            stack = sparse(1:numel(row_g), ...
                                repmat((1:npre)', ng, 1), coeff(:), ...
                                numel(row_g), npre);
                            J(row_g, row_a) = kron(sparse(ones(1, na)), stack);
                        end
                        tau_fac_stack = kron(params.tau_g_fac{pre, post}(:), ...
                            ones(npre, 1));
                        diagonal = kron(-1 ./ params.tau_g_dec{pre, post}(:), ...
                            ones(npre, 1)) - repmat(pre_rate, ng, 1) ./ tau_fac_stack;
                        J(row_g, row_g) = spdiags( ...
                            diagonal, 0, numel(row_g), numel(row_g));
                        G_stack = kron(params.G{pre, post}(:), ones(npre, 1));
                        vals = ((G_stack - g_state{pre, post}(:)) .* ...
                            repmat(pre_rate_prime, ng, 1)) ./ tau_fac_stack;
                        J(row_g, row_x) = sparse(1:numel(row_g), ...
                            repmat(pre_idx(:), ng, 1), vals, numel(row_g), n);
                    end

                    W_block = params.W(post_idx, pre_idx);
                    synaptic_modulation = facilitation{pre, post} .* depression_syn;
                    rate_block = W_block * spdiags( ...
                        synaptic_modulation .* pre_rate_prime, 0, npre, npre);
                    J(row_x(post_idx), row_x(pre_idx)) = ...
                        J(row_x(post_idx), row_x(pre_idx)) + rate_block ./ params.tau_d;

                    if na > 0
                        replicate = kron(ones(1, na), speye(npre));
                        block = -params.c(pre) .* rate_block;
                        J(row_x(post_idx), row_a) = ...
                            J(row_x(post_idx), row_a) + ...
                            (block * replicate) ./ params.tau_d;
                    end
                    if nb > 0
                        product_derivative = zeros(npre, nb);
                        for m = 1:nb
                            other = [1:m-1, m+1:nb];
                            if isempty(other)
                                product_derivative(:, m) = 1;
                            else
                                product_derivative(:, m) = ...
                                    prod(b_state{pre, post}(:, other), 2);
                            end
                        end
                        gain = 1;
                        if params.std_zero_floor
                            p_min = prod(params.tau_b_rel{pre, post} ./ ...
                                (params.tau_b_rec{pre, post} + ...
                                params.tau_b_rel{pre, post}));
                            gain = 1 / (1 - p_min);
                        end
                        coeff = gain .* facilitation{pre, post} .* ...
                            pre_rate .* product_derivative;
                        D = sparse(repmat((1:npre)', nb, 1), ...
                            1:numel(row_b), coeff(:), npre, numel(row_b));
                        J(row_x(post_idx), row_b) = ...
                            (W_block * D) ./ params.tau_d;
                    end
                    if ng > 0
                        product_derivative = zeros(npre, ng);
                        for k = 1:ng
                            other = [1:k-1, k+1:ng];
                            if isempty(other)
                                product_derivative(:, k) = 1;
                            else
                                product_derivative(:, k) = ...
                                    prod(g_state{pre, post}(:, other), 2);
                            end
                        end
                        coeff = depression_syn .* pre_rate .* product_derivative;
                        D = sparse(repmat((1:npre)', ng, 1), ...
                            1:numel(row_g), coeff(:), npre, numel(row_g));
                        J(row_x(post_idx), row_g) = ...
                            (W_block * D) ./ params.tau_d;
                    end
                end
            end
        end

        function J_array = compute_Jacobian_at_indices(S_out, sample_indices, params)
            J_array = zeros(params.N_sys_eqs, params.N_sys_eqs, numel(sample_indices));
            for k = 1:numel(sample_indices)
                J_array(:, :, k) = full(SRNNCellTypePairs.compute_Jacobian_fast( ...
                    S_out(sample_indices(k), :)', params));
            end
        end
    end

    %% State conversion and plotting helpers
    methods (Static)
        function [x_named, a_named, b_named, r_named, ...
                synaptic_output_named, g_named] = ...
                unpack_and_compute_states(S_out, params)
            nt = size(S_out, 1);
            layout = params.state_layout;
            x_all = S_out(:, layout.x)';
            x_eff = x_all;
            a_named = struct(); b_named = struct(); g_named = struct();

            for q = 1:params.n_cellTypes
                name = params.cell_type_names{q};
                idx = params.type_indices{q};
                nq = params.n_per_type(q);
                if params.n_a(q) > 0
                    aq = reshape(S_out(:, layout.a{q})', nq, params.n_a(q), nt);
                    summed = reshape(sum(aq, 2), nq, nt);
                    x_eff(idx, :) = x_eff(idx, :) - params.c(q) .* summed;
                    a_named.(name) = aq;
                else
                    a_named.(name) = [];
                end
            end
            rate = params.activation_function(x_eff);
            synaptic_output_named = struct();
            for pre = 1:params.n_cellTypes
                pre_name = params.cell_type_names{pre};
                pre_idx = params.type_indices{pre};
                npre = params.n_per_type(pre);
                b_named.(pre_name) = struct();
                g_named.(pre_name) = struct();
                synaptic_output_named.(pre_name) = struct();
                for post = 1:params.n_cellTypes
                    post_name = params.cell_type_names{post};
                    nb = params.n_b_pairs(pre, post);
                    ng = params.n_g_pairs(pre, post);
                    depression = ones(npre, nt);
                    facilitation = ones(npre, nt);
                    if nb > 0
                        bq = reshape(S_out(:, layout.b{pre, post})', ...
                            npre, nb, nt);
                        depression = reshape(prod(bq, 2), npre, nt);
                        b_named.(pre_name).(post_name) = bq;
                    else
                        b_named.(pre_name).(post_name) = [];
                    end
                    if ng > 0
                        gq = reshape(S_out(:, layout.g{pre, post})', ...
                            npre, ng, nt);
                        facilitation = reshape(prod(gq, 2), npre, nt);
                        g_named.(pre_name).(post_name) = gq;
                    else
                        g_named.(pre_name).(post_name) = [];
                    end
                    depression = SRNNCellTypePairs.apply_std_zero_floor_pair( ...
                        depression, params, pre, post);
                    synaptic_output_named.(pre_name).(post_name) = ...
                        facilitation .* depression .* rate(pre_idx, :);
                end
            end
            x_named = SRNNCellTypePairs.split_by_type(x_all, params, 1);
            r_named = SRNNCellTypePairs.split_by_type(rate, params, 1);
        end

        function named = split_by_type(data, params, neuron_dimension)
            named = struct();
            subs = repmat({':'}, 1, ndims(data));
            for q = 1:params.n_cellTypes
                subs{neuron_dimension} = params.type_indices{q};
                named.(params.cell_type_names{q}) = data(subs{:});
            end
        end

        function output = trim_named(input, dimension, mask)
            output = input;
            names = fieldnames(input);
            for q = 1:numel(names)
                value = input.(names{q});
                if isempty(value), continue; end
                subs = repmat({':'}, 1, ndims(value));
                subs{dimension} = mask;
                output.(names{q}) = value(subs{:});
            end
        end

        function output = trim_routes(input, dimension, mask)
            output = input;
            pre_names = fieldnames(input);
            for pre = 1:numel(pre_names)
                post_names = fieldnames(input.(pre_names{pre}));
                for post = 1:numel(post_names)
                    value = input.(pre_names{pre}).(post_names{post});
                    if isempty(value), continue; end
                    subs = repmat({':'}, 1, ndims(value));
                    subs{dimension} = mask;
                    output.(pre_names{pre}).(post_names{post}) = value(subs{:});
                end
            end
        end

        function collapsed = collapse_routes(input, operation)
            collapsed = input;
            pre_names = fieldnames(input);
            for pre = 1:numel(pre_names)
                post_names = fieldnames(input.(pre_names{pre}));
                for post = 1:numel(post_names)
                    value = input.(pre_names{pre}).(post_names{post});
                    if isempty(value), continue; end
                    if strcmp(operation, 'prod')
                        value = prod(value, 2);
                    else
                        value = sum(value, 2);
                    end
                    collapsed.(pre_names{pre}).(post_names{post}) = ...
                        reshape(value, size(value, 1), []);
                end
            end
        end

        function output = collapse_named(input, names, operation)
            output = struct();
            for q = 1:numel(names)
                value = input.(names{q});
                if isempty(value)
                    output.(names{q}) = [];
                elseif strcmp(operation, 'sum')
                    output.(names{q}) = reshape(sum(value, 2), size(value, 1), []);
                else
                    output.(names{q}) = reshape(prod(value, 2), size(value, 1), []);
                end
            end
        end

        function plot_named_series(t, data, names, label, show_legend)
            if nargin < 5, show_legend = true; end
            colors = SRNNCellTypePairs.type_colors(numel(names));
            hold on;
            handles = gobjects(0); labels = {};
            for q = 1:numel(names)
                values = data.(names{q});
                if isempty(values), continue; end
                SRNNCellTypePairs.plot_celltype_lines( ...
                    gca, t, values, colors(q, :));
                handles(end + 1) = plot(nan, nan, '-', 'Color', colors(q, :), ...
                    'LineWidth', 1.5); %#ok<AGROW>
                labels{end + 1} = names{q}; %#ok<AGROW>
            end
            hold off;
            ylabel(label);
            if show_legend && ~isempty(handles)
                legend(handles, labels, 'Location', 'best');
            end
        end

        function plot_route_series(t, data, names, label, show_legend)
            if nargin < 5, show_legend = true; end
            colors = lines(max(1, numel(names)^2));
            hold on;
            handles = gobjects(0); labels = {};
            color_index = 0;
            for pre = 1:numel(names)
                pre_name = names{pre};
                for post = 1:numel(names)
                    post_name = names{post};
                    values = data.(pre_name).(post_name);
                    if isempty(values), continue; end
                    color_index = color_index + 1;
                    handles(end + 1) = plot(t, mean(values, 1), ...
                        'Color', colors(color_index, :)); %#ok<AGROW>
                    labels{end + 1} = sprintf('%s->%s', ...
                        pre_name, post_name); %#ok<AGROW>
                end
            end
            hold off;
            ylabel(label);
            if show_legend && ~isempty(handles)
                legend(handles, labels, 'Interpreter', 'none', 'Location', 'best');
            end
        end

        function plot_postsynaptic_route_lines(ax, t, routes, post_names)
            % Plot every presynaptic-neuron trace, colored by target type.
            colors = SRNNCellTypePairs.type_colors(numel(post_names));
            handles = gobjects(0); labels = {};
            hold(ax, 'on');
            for post = 1:numel(post_names)
                name = post_names{post};
                values = routes.(name);
                if isempty(values), continue; end
                route_lines = SRNNCellTypePairs.plot_celltype_lines( ...
                    ax, t, values, colors(post, :));
                set(route_lines, 'Tag', ['PostType_' name]);
                handles(end + 1) = plot(ax, nan, nan, '-', ...
                    'Color', colors(post, :), 'LineWidth', 1.5); %#ok<AGROW>
                labels{end + 1} = name; %#ok<AGROW>
            end
            hold(ax, 'off');
            if isempty(handles)
                SRNNCellTypePairs.show_empty_axis(ax, 'No route state');
            else
                legend(ax, handles, labels, 'Interpreter', 'none', 'Location', 'best');
            end
        end

        function show_empty_axis(ax, label)
            text(ax, 0.5, 0.5, label, 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'Color', [0.4 0.4 0.4]);
        end

        function C = type_colors(n)
            % TYPE_COLORS Per-cell-type colours: lines(), with the first two
            % swapped.
            %
            % MATLAB's lines() starts blue then orange-red. With cell types
            % conventionally ordered {E, I}, that paints excitation blue and
            % inhibition red -- the opposite of the usual reading. Swapping the
            % first two rows gives E a warm colour and I a cool one, and leaves
            % every later type untouched.
            %
            % Used only where the index IS a cell type. Per-neuron accents and
            % per-route colours keep plain lines(), since their index means
            % something else.
            C = lines(n);
            if n >= 2
                C([1 2], :) = C([2 1], :);
            end
        end

        function line_handles = plot_celltype_lines(ax, t, values, type_color)
            %PLOT_CELLTYPE_LINES Blend cell-type and neuron-specific colors.
            %
            % The cell type supplies the dominant color while lines() adds a
            % repeatable accent for each neuron. This keeps populations easy
            % to identify while making overlapping neurons distinguishable.
            n_neurons = size(values, 1);
            neuron_accents = lines(n_neurons);
            neuron_colors = 0.5 .* repmat(type_color, n_neurons, 1) ...
                + 0.5 .* neuron_accents;

            line_handles = plot(ax, t, values');
            for neuron = 1:numel(line_handles)
                set(line_handles(neuron), ...
                    'Color', neuron_colors(neuron, :), ...
                    'HandleVisibility', 'off');
            end
        end

        function y = piecewiseSigmoid(x, a, c)
            % c may be a scalar or an array that broadcasts against x (an n x 1
            % per-neuron S_c_vec against an n x 1 state or an n x nt
            % trajectory). The curve is a pure TRANSLATION in c -- every
            % threshold is c + const and the linear branch is (x - c) + 0.5 --
            % so a per-neuron centre is exactly the c = 0 curve evaluated at
            % x - c. That keeps one implementation of the branch structure and
            % leaves the scalar path bit-identical.
            if ~isscalar(c)
                y = SRNNCellTypePairs.piecewiseSigmoid(x - c, a, 0);
                return;
            end
            if a < 0 || a > 1
                error('SRNNCellTypePairs:InvalidActivation', 'a must be between 0 and 1.');
            end
            a = a / 2;
            if a == 0.5
                y = min(max((x - c) + 0.5, 0), 1);
                return;
            end
            y = zeros(size(x));
            k = 0.5 / (1 - 2*a);
            x1 = c + a - 1; x2 = c - a; x3 = c + a; x4 = c + 1 - a;
            left = x >= x1 & x < x2;
            linear = x >= x2 & x <= x3;
            right = x > x3 & x <= x4;
            y(left) = k .* (x(left) - x1).^2;
            y(linear) = (x(linear) - c) + 0.5;
            y(right) = 1 - k .* (x(right) - x4).^2;
            y(x > x4) = 1;
        end

        function dy = piecewiseSigmoidDerivative(x, a, c)
            % As with piecewiseSigmoid, c may be a scalar or an array that
            % broadcasts against x, and the derivative is likewise a pure
            % translation in c.
            if ~isscalar(c)
                dy = SRNNCellTypePairs.piecewiseSigmoidDerivative(x - c, a, 0);
                return;
            end
            if a < 0 || a > 1
                error('SRNNCellTypePairs:InvalidActivation', 'a must be between 0 and 1.');
            end
            a = a / 2;
            dy = zeros(size(x), 'like', x);
            if a == 0.5
                dy(x >= c - 0.5 & x <= c + 0.5) = 1;
                return;
            end
            k = 0.5 / (1 - 2*a);
            x1 = c + a - 1; x2 = c - a; x3 = c + a; x4 = c + 1 - a;
            left = x >= x1 & x < x2;
            linear = x >= x2 & x <= x3;
            right = x > x3 & x <= x4;
            dy(left) = 2 .* k .* (x(left) - x1);
            dy(linear) = 1;
            dy(right) = -2 .* k .* (x(right) - x4);
        end

        function y = logisticSigmoid(x, c)
            if nargin < 2, c = 0; end
            y = 1 ./ (1 + exp(-4 .* (x - c)));
        end

        function dy = logisticSigmoidDerivative(x, c)
            if nargin < 2, c = 0; end
            value = SRNNCellTypePairs.logisticSigmoid(x, c);
            dy = 4 .* value .* (1 - value);
        end

        function y = tanhActivation(x)
            y = tanh(x);
        end

        function dy = tanhActivationDerivative(x)
            dy = 1 - tanh(x).^2;
        end
    end

    %% Lyapunov algorithms
    methods (Static, Access = protected)
        function results = compute_lyapunov_exponents_internal(method, S_out, t_out, ...
                dt, fs, interval, params, opts, ode_solver, rhs)
            results = struct();
            switch lower(method)
                case 'benettin'
                    lya_dt = 0.02;
                case 'qr'
                    lya_dt = 0.1;
                case 'none'
                    return;
                otherwise
                    error('SRNNCellTypePairs:UnknownLyapunovMethod', ...
                        'Unknown Lyapunov method: %s', method);
            end
            factor = lya_dt / dt;
            if abs(factor - round(factor)) > 1e-11 || factor < 3
                error('SRNNCellTypePairs:InvalidLyapunovStep', ...
                    'Lyapunov interval must be an integer multiple of at least 3*dt.');
            end

            if strcmpi(method, 'benettin')
                [LLE, local, finite, times] = SRNNCellTypePairs.benettin_algorithm_internal( ...
                    S_out, t_out, dt, fs, 1e-3, interval, lya_dt, ...
                    params, opts, rhs, ode_solver);
                results.LLE = LLE;
                results.local_lya = local;
                results.finite_lya = finite;
                results.t_lya = times;
            else
                [spectrum, local, finite, times] = ...
                    SRNNCellTypePairs.lyapunov_spectrum_qr_internal( ...
                    S_out, t_out, lya_dt, params, ode_solver, opts, fs);
                [spectrum, order] = sort(real(spectrum), 'descend');
                results.LE_spectrum = spectrum;
                results.local_LE_spectrum_t = local(:, order);
                results.finite_LE_spectrum_t = finite(:, order);
                results.t_lya = times;
                results.sort_idx = order;
                results.params.N_sys_eqs = params.N_sys_eqs;
                fprintf('Lyapunov Dimension: %.2f\n', ...
                    SRNNCellTypePairs.compute_kaplan_yorke_dimension_internal(spectrum));
            end
            results.lya_dt = lya_dt;
            results.lya_fs = 1 / lya_dt;
        end

        function [LLE, local_lya, finite_lya, t_lya] = ...
                benettin_algorithm_internal(X, t, dt, fs, d0, interval, ...
                lya_dt, ~, ode_options, dynamics_func, ode_solver)
            decimation = round(lya_dt * fs);
            tau = dt * decimation;
            t_lya = t(1:decimation:end);
            if ~isempty(t_lya) && t_lya(end) + tau > t(end)
                t_lya(end) = [];
            end
            n_intervals = numel(t_lya);
            local_lya = zeros(n_intervals, 1);
            finite_lya = nan(n_intervals, 1);
            perturbation = randn(size(X, 2), 1);
            perturbation = d0 .* perturbation ./ norm(perturbation);
            accumulated = 0;

            for k = 1:n_intervals
                start_index = (k - 1) * decimation + 1;
                end_index = start_index + decimation;
                perturbed_start = X(start_index, :)' + perturbation;
                detailed_times = t_lya(k) + (0:dt:tau);
                [~, perturbed] = ode_solver( ...
                    dynamics_func, detailed_times, perturbed_start, ode_options);
                delta = perturbed(end, :)' - X(end_index, :)';
                distance = norm(delta);
                local_lya(k) = log(distance / d0) / tau;
                if ~isfinite(local_lya(k)) || distance <= eps
                    warning('SRNNCellTypePairs:LyapunovDiverged', ...
                        'Lyapunov trajectory diverged at t=%g.', t_lya(k));
                    local_lya(k:end) = [];
                    finite_lya(k:end) = [];
                    t_lya(k:end) = [];
                    break;
                end
                perturbation = d0 .* delta ./ distance;
                if t_lya(k) >= max(0, interval(1)) && t_lya(k) < interval(2)
                    accumulated = accumulated + log(distance / d0);
                    elapsed = t_lya(k) + tau - max(0, interval(1));
                    finite_lya(k) = accumulated / max(elapsed, eps);
                end
            end
            valid = finite_lya(isfinite(finite_lya));
            if isempty(valid), LLE = 0; else, LLE = valid(end); end
        end

        function [spectrum, local_spectrum, finite_spectrum, t_lya] = ...
                lyapunov_spectrum_qr_internal(X, t, lya_dt, params, ...
                ode_solver, ode_options, fs)
            N = params.N_sys_eqs;
            interpolants = cell(N, 1);
            for state = 1:N
                interpolants{state} = griddedInterpolant(t, X(:, state), 'pchip');
            end
            dt = 1 / fs;
            decimation = round(lya_dt / dt);
            tau = decimation * dt;
            sample_indices = 1:decimation:numel(t);
            t_lya = t(sample_indices);
            if ~isempty(t_lya) && t_lya(end) + tau > t(end) + eps(t(end))
                t_lya(end) = [];
            end
            n_intervals = numel(t_lya);
            local_spectrum = zeros(n_intervals, N);
            finite_spectrum = nan(n_intervals, N);
            Q = eye(N);
            accumulated = zeros(N, 1);
            elapsed = 0;
            variational_options = odeset(ode_options, 'Jacobian', []);

            for k = 1:n_intervals
                t_start = t_lya(k);
                t_end = min(t_start + tau, t(end));
                psi0 = Q(:);
                variational = @(tt, psi) SRNNCellTypePairs.variational_equations( ...
                    tt, psi, interpolants, N, params);
                [~, solution] = ode_solver(variational, [t_start, t_end], ...
                    psi0, variational_options);
                evolved = reshape(solution(end, :)', N, N);
                if any(~isfinite(evolved(:)))
                    warning('SRNNCellTypePairs:LyapunovDiverged', ...
                        'QR trajectory diverged at t=%g.', t_start);
                    t_lya(k:end) = [];
                    local_spectrum(k:end, :) = [];
                    finite_spectrum(k:end, :) = [];
                    break;
                end
                [Q, R] = qr(evolved);
                log_diag = log(max(abs(diag(R)), realmin));
                local_spectrum(k, :) = (log_diag ./ (t_end - t_start))';
                if t_start >= 0
                    accumulated = accumulated + log_diag;
                    elapsed = elapsed + (t_end - t_start);
                    finite_spectrum(k, :) = (accumulated ./ elapsed)';
                end
            end
            if elapsed > 0
                spectrum = accumulated ./ elapsed;
            else
                spectrum = nan(N, 1);
            end
        end

        function derivative = variational_equations(t, psi, interpolants, N, params)
            state = zeros(N, 1);
            for k = 1:N
                state(k) = interpolants{k}(t);
            end
            J = SRNNCellTypePairs.compute_Jacobian_fast(state, params);
            derivative = reshape(J * reshape(psi, N, N), [], 1);
        end

        function dimension = compute_kaplan_yorke_dimension_internal(lambda)
            lambda = sort(lambda, 'descend');
            sums = cumsum(lambda);
            j = find(sums >= 0, 1, 'last');
            if isempty(j)
                dimension = 0;
            elseif j == numel(lambda)
                dimension = numel(lambda);
            else
                dimension = j + sums(j) / abs(lambda(j + 1));
            end
        end
    end
end
