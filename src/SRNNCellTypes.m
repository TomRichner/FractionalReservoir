classdef SRNNCellTypes < handle
    %SRNNCELLTYPES Stable recurrent nonlinear network with arbitrary cell types.
    %
    % Cell-type-specific parameters are ordered according to cell_type_names.
    % Numeric vectors hold one scalar per type; tau_a, tau_b_rec, and
    % tau_g_dec are cell arrays because their per-type vectors may have
    % different lengths.
    %
    % State ordering is timescale-major within each type:
    %   S = [a{1}(:); ...; a{C}(:); b{1}(:); ...; b{C}(:); ...
    %        g{1}(:); ...; g{C}(:); x(:)]
    %
    % Example:
    %   model = SRNNCellTypes( ...
    %       'n_cellTypes', 3, 'cell_type_names', {'E','PV','SST'}, ...
    %       'f', [.8 .1 .1], 'mu_tilde', [.08 -.12 -.1], ...
    %       'sigma_tilde', [.02 .02 .02], 'n_a', [3 0 1], ...
    %       'tau_a', {[.25 1 10], [], 2}, 'n_b', [1 0 2], ...
    %       'tau_b_rec', {1, [], [.5 2]});
    %   model.build(); model.run(); model.plot();

    %% Cell types and connectivity
    properties
        n = 300
        indegree = 100
        n_cellTypes
        cell_type_names
        f
        mu_tilde
        sigma_tilde
        level_of_chaos = 1.0
        rescale_by_abscissa = false
    end

    %% Per-type SFA, STD, and STF
    properties
        n_a
        tau_a
        c
        n_b
        tau_b_rec
        tau_b_rel
        std_zero_floor = false
        n_g
        tau_g_dec
        tau_g_fac
        G
    end

    %% Shared dynamics and simulation settings
    properties
        tau_d = 0.1
        activation_function
        activation_function_derivative
        S_a = 0.9
        S_c = 0.40                      % matches SRNNModel2's default
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
        n_per_type
        type_indices
        cell_indices
        N_sys_eqs
    end

    properties (SetAccess = protected)
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
    end

    methods
        function obj = SRNNCellTypes(varargin)
            obj.set_defaults();
            if mod(numel(varargin), 2) ~= 0
                error('SRNNCellTypes:InvalidInput', ...
                    'Constructor arguments must be name-value pairs.');
            end
            supplied = string(varargin(1:2:end));
            required = ["n_cellTypes", "cell_type_names", "f", ...
                "mu_tilde", "sigma_tilde"];
            missing = required(~ismember(required, supplied));
            if ~isempty(missing)
                error('SRNNCellTypes:MissingConfig', ...
                    'Required constructor properties: %s.', strjoin(cellstr(missing), ', '));
            end
            for k = 1:2:numel(varargin)
                name = varargin{k};
                if ~(ischar(name) || (isstring(name) && isscalar(name))) || ...
                        ~isprop(obj, char(name))
                    error('SRNNCellTypes:UnknownProperty', ...
                        'Unknown property: %s', string(name));
                end
                obj.(char(name)) = varargin{k + 1};
            end
            obj.complete_type_defaults();
            obj.validate();
            if isempty(obj.plot_deci)
                obj.plot_deci = max(1, round(obj.fs / obj.plot_freq));
            end
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
            if isempty(obj.n_a) || isempty(obj.n_b) || isempty(obj.n_g)
                value = obj.n;
            else
                value = obj.n + sum(obj.n_per_type .* ...
                    (obj.n_a + obj.n_b + obj.n_g));
            end
        end

        function build(obj)
            obj.complete_type_defaults();
            obj.validate();
            obj.build_network();
            obj.build_stimulus();
            obj.cached_params = obj.get_params();
            obj.is_built = true;
            fprintf('SRNNCellTypes built successfully. Ready to run.\n');
        end

        function run(obj)
            if ~obj.is_built
                error('SRNNCellTypes:NotBuilt', ...
                    'Model must be built before running. Call build() first.');
            end

            params = obj.cached_params;
            params.u_interpolant = obj.u_interpolant;
            dt = 1 / obj.fs;
            if isempty(obj.ode_opts)
                jacobian = @(~, S) SRNNCellTypes.compute_Jacobian_fast(S, params);
                obj.ode_opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, ...
                    'MaxStep', dt, 'Jacobian', jacobian);
            end
            rhs = @(t, S) SRNNCellTypes.dynamics_fast(t, S, params);

            fprintf('Integrating SRNNCellTypes equations\n');
            tic;
            [t_raw, S_raw] = obj.ode_solver(rhs, obj.t_ex, obj.S0, obj.ode_opts);
            fprintf('Integration complete in %.2f seconds.\n', toc);
            if numel(t_raw) ~= numel(obj.t_ex) || ...
                    max(abs(t_raw(:) - obj.t_ex(:))) > 1e-9
                error('SRNNCellTypes:TimeMismatch', ...
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
                error('SRNNCellTypes:NoStateData', 'State data are not available.');
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
            rhs = @(t, S) SRNNCellTypes.dynamics_fast(t, S, params);
            obj.lya_results = SRNNCellTypes.compute_lyapunov_exponents_internal( ...
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
            params.n_b = obj.n_b;
            params.tau_b_rec = obj.tau_b_rec;
            params.tau_b_rel = obj.tau_b_rel;
            params.std_zero_floor = obj.std_zero_floor;
            params.n_g = obj.n_g;
            params.tau_g_dec = obj.tau_g_dec;
            params.tau_g_fac = obj.tau_g_fac;
            params.G = obj.G;
            params.tau_d = obj.tau_d;
            params.activation_function = obj.activation_function;
            params.activation_function_derivative = obj.activation_function_derivative;
            params.x0_std = obj.x0_std;
            params.N_sys_eqs = obj.N_sys_eqs;
            params.state_layout = SRNNCellTypes.make_state_layout( ...
                obj.n, obj.n_per_type, obj.n_a, obj.n_b, obj.n_g);
            params.rng_seeds = obj.rng_seeds;
            if ~isempty(obj.W)
                params.W = obj.W;
            end
        end

        function dS_dt = dynamics(obj, t, S)
            if ~obj.is_built
                error('SRNNCellTypes:NotBuilt', 'Call build() before dynamics().');
            end
            params = obj.cached_params;
            params.u_interpolant = obj.u_interpolant;
            dS_dt = SRNNCellTypes.dynamics_fast(t, S, params);
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
            obj.clear_results();
        end
    end

    methods (Access = protected)
        function set_defaults(obj)
            obj.activation_function = @(x) SRNNCellTypes.logisticSigmoid(x, obj.S_c);
            obj.activation_function_derivative = ...
                @(x) SRNNCellTypes.logisticSigmoidDerivative(x, obj.S_c);
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
            obj.mu_tilde = reshape(obj.mu_tilde, 1, []);
            obj.sigma_tilde = reshape(obj.sigma_tilde, 1, []);
            if isempty(obj.n_a), obj.n_a = zeros(1, C); end
            if isempty(obj.n_b), obj.n_b = zeros(1, C); end
            if isempty(obj.n_g), obj.n_g = zeros(1, C); end
            if isempty(obj.c), obj.c = repmat(0.15 / 3, 1, C); end
            if isempty(obj.tau_b_rel), obj.tau_b_rel = repmat(0.25, 1, C); end
            if isempty(obj.tau_g_fac), obj.tau_g_fac = repmat(0.25, 1, C); end
            if isempty(obj.G), obj.G = repmat(2, 1, C); end
            obj.n_a = reshape(obj.n_a, 1, []);
            obj.n_b = reshape(obj.n_b, 1, []);
            obj.n_g = reshape(obj.n_g, 1, []);
            obj.c = reshape(obj.c, 1, []);
            obj.tau_b_rel = reshape(obj.tau_b_rel, 1, []);
            obj.tau_g_fac = reshape(obj.tau_g_fac, 1, []);
            obj.G = reshape(obj.G, 1, []);

            if isempty(obj.tau_a), obj.tau_a = cell(1, C); end
            if isempty(obj.tau_b_rec), obj.tau_b_rec = cell(1, C); end
            if isempty(obj.tau_g_dec), obj.tau_g_dec = cell(1, C); end
            if iscell(obj.tau_a), obj.tau_a = reshape(obj.tau_a, 1, []); end
            if iscell(obj.tau_b_rec), obj.tau_b_rec = reshape(obj.tau_b_rec, 1, []); end
            if iscell(obj.tau_g_dec), obj.tau_g_dec = reshape(obj.tau_g_dec, 1, []); end
            if iscell(obj.tau_a) && numel(obj.tau_a) == C
                for q = 1:C
                    if obj.n_a(q) > 0 && isempty(obj.tau_a{q})
                        obj.tau_a{q} = logspace(log10(0.25), log10(10), obj.n_a(q));
                    elseif ~isempty(obj.tau_a{q})
                        obj.tau_a{q} = reshape(obj.tau_a{q}, 1, []);
                    end
                end
            end
            if iscell(obj.tau_b_rec) && numel(obj.tau_b_rec) == C
                for q = 1:C
                    if obj.n_b(q) == 1 && isempty(obj.tau_b_rec{q})
                        obj.tau_b_rec{q} = 1;
                    elseif ~isempty(obj.tau_b_rec{q})
                        obj.tau_b_rec{q} = reshape(obj.tau_b_rec{q}, 1, []);
                    end
                end
            end
            if iscell(obj.tau_g_dec) && numel(obj.tau_g_dec) == C
                for q = 1:C
                    if obj.n_g(q) == 1 && isempty(obj.tau_g_dec{q})
                        obj.tau_g_dec{q} = 1;
                    elseif ~isempty(obj.tau_g_dec{q})
                        obj.tau_g_dec{q} = reshape(obj.tau_g_dec{q}, 1, []);
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
                error('SRNNCellTypes:InvalidParams', 'n_cellTypes must be a positive integer.');
            end
            if numel(obj.cell_type_names) ~= C || ...
                    any(~cellfun(@isvarname, obj.cell_type_names)) || ...
                    numel(unique(obj.cell_type_names)) ~= C
                error('SRNNCellTypes:InvalidParams', ...
                    'cell_type_names must contain unique valid MATLAB identifiers.');
            end
            if ~isscalar(obj.n) || obj.n < C || obj.n ~= round(obj.n)
                error('SRNNCellTypes:InvalidParams', ...
                    'n must be an integer at least as large as n_cellTypes.');
            end
            if ~isscalar(obj.indegree) || obj.indegree <= 0 || obj.indegree > obj.n
                error('SRNNCellTypes:InvalidParams', 'indegree must satisfy 0 < indegree <= n.');
            end
            if numel(obj.f) ~= C || any(~isfinite(obj.f)) || any(obj.f <= 0) || ...
                    abs(sum(obj.f) - 1) > 1e-12
                error('SRNNCellTypes:InvalidParams', ...
                    'f must contain one positive fraction per type and sum to 1.');
            end
            RMTCellTypes.allocate_counts(obj.n, obj.f);
            if numel(obj.mu_tilde) ~= C || any(~isfinite(obj.mu_tilde)) || ...
                    numel(obj.sigma_tilde) ~= C || any(~isfinite(obj.sigma_tilde)) || ...
                    any(obj.sigma_tilde < 0)
                error('SRNNCellTypes:InvalidParams', ...
                    'mu_tilde and sigma_tilde must have one valid value per type.');
            end

            vector_fields = {'n_a', 'n_b', 'n_g', 'c', 'tau_b_rel', ...
                'tau_g_fac', 'G'};
            for k = 1:numel(vector_fields)
                name = vector_fields{k};
                if ~isnumeric(obj.(name)) || numel(obj.(name)) ~= C || ...
                        any(~isfinite(obj.(name)))
                    error('SRNNCellTypes:InvalidParams', ...
                        '%s must have one finite numeric value per cell type.', name);
                end
            end
            if any(obj.n_a < 0 | obj.n_a ~= round(obj.n_a)) || ...
                    any(obj.n_b < 0 | obj.n_b ~= round(obj.n_b)) || ...
                    any(obj.n_g < 0 | obj.n_g ~= round(obj.n_g))
                error('SRNNCellTypes:InvalidParams', ...
                    'n_a, n_b, and n_g must contain nonnegative integers.');
            end
            if any(obj.c < 0) || any(obj.tau_b_rel <= 0) || ...
                    any(obj.tau_g_fac <= 0) || any(obj.G < 1)
                error('SRNNCellTypes:InvalidParams', ...
                    ['c must be nonnegative, tau_b_rel and tau_g_fac must ' ...
                    'be positive, and G must be at least one.']);
            end
            if ~iscell(obj.tau_a) || numel(obj.tau_a) ~= C || ...
                    ~iscell(obj.tau_b_rec) || numel(obj.tau_b_rec) ~= C || ...
                    ~iscell(obj.tau_g_dec) || numel(obj.tau_g_dec) ~= C
                error('SRNNCellTypes:InvalidParams', ...
                    ['tau_a, tau_b_rec, and tau_g_dec must be ' ...
                    '1-by-n_cellTypes cell arrays.']);
            end
            for q = 1:C
                if numel(obj.tau_a{q}) ~= obj.n_a(q) || ...
                        any(~isfinite(obj.tau_a{q})) || any(obj.tau_a{q} <= 0)
                    error('SRNNCellTypes:InvalidParams', ...
                        'tau_a{%d} must contain n_a(%d) positive values.', q, q);
                end
                if numel(obj.tau_b_rec{q}) ~= obj.n_b(q) || ...
                        any(~isfinite(obj.tau_b_rec{q})) || any(obj.tau_b_rec{q} <= 0)
                    error('SRNNCellTypes:InvalidParams', ...
                        'tau_b_rec{%d} must contain n_b(%d) positive values.', q, q);
                end
                if numel(obj.tau_g_dec{q}) ~= obj.n_g(q) || ...
                        any(~isfinite(obj.tau_g_dec{q})) || any(obj.tau_g_dec{q} <= 0)
                    error('SRNNCellTypes:InvalidParams', ...
                        'tau_g_dec{%d} must contain n_g(%d) positive values.', q, q);
                end
            end
            if obj.T_range(2) <= obj.T_range(1) || obj.fs <= 0 || obj.tau_d <= 0
                error('SRNNCellTypes:InvalidParams', ...
                    'T_range, fs, and tau_d must define positive simulation intervals.');
            end
            has_generator = isfield(obj.input_config, 'generator') && ...
                isa(obj.input_config.generator, 'function_handle');
            if ~has_generator
                if ~isfield(obj.input_config, 'step_density') || ...
                        ~isstruct(obj.input_config.step_density)
                    error('SRNNCellTypes:InvalidParams', ...
                        'input_config.step_density must be a named struct.');
                end
                density_fields = sort(fieldnames(obj.input_config.step_density));
                if ~isequal(density_fields, sort(obj.cell_type_names(:)))
                    error('SRNNCellTypes:InvalidParams', ...
                        'step_density fields must exactly match cell_type_names.');
                end
                for q = 1:C
                    d = obj.input_config.step_density.(obj.cell_type_names{q});
                    if ~isscalar(d) || ~isfinite(d) || d < 0 || d > 1
                        error('SRNNCellTypes:InvalidParams', ...
                            'Each step density must satisfy 0 <= density <= 1.');
                    end
                end
            end
        end

        function build_network(obj)
            rng(obj.rng_seeds(1));
            rmt = RMTCellTypes(obj.n, obj.alpha, obj.f, ...
                obj.mu_tilde, obj.sigma_tilde);
            W_raw = rmt.W;
            if obj.rescale_by_abscissa
                abscissa = max(real(eig(full(W_raw))));
                if abs(abscissa) <= eps
                    error('SRNNCellTypes:InvalidConnectivity', ...
                        'Cannot rescale a matrix with zero spectral abscissa.');
                end
                W_raw = W_raw / abscissa;
            end
            obj.W = obj.level_of_chaos .* W_raw;
            eig_W = eig(full(obj.W));
            fprintf('W created: spectral radius = %.3f, abscissa = %.3f\n', ...
                max(abs(eig_W)), max(real(eig_W)));
        end

        function build_stimulus(obj)
            obj.generate_stimulus();
            obj.u_interpolant = griddedInterpolant( ...
                obj.t_ex, obj.u_ex', 'linear', 'none');
            obj.S0 = SRNNCellTypes.initialize_state(obj.get_params());
        end

        function generate_stimulus(obj)
            dt = 1 / obj.fs;
            if isempty(obj.input_config.intrinsic_drive)
                obj.input_config.intrinsic_drive = zeros(obj.n, 1);
            end
            params = obj.get_params();
            [u_stim, t_stim] = SRNNCellTypes.generate_external_input( ...
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
            [x, a, b, r, br, g] = SRNNCellTypes.unpack_and_compute_states( ...
                S_plot, obj.cached_params);
            u = SRNNCellTypes.split_by_type( ...
                obj.u_ex(:, sample_indices), obj.cached_params, 1);
            range = obj.T_plot;
            if isempty(range), range = obj.T_range; end
            keep = t_plot >= range(1) & t_plot <= range(2);
            obj.plot_data = struct( ...
                't', t_plot(keep), ...
                'u', SRNNCellTypes.trim_named(u, 2, keep), ...
                'x', SRNNCellTypes.trim_named(x, 2, keep), ...
                'r', SRNNCellTypes.trim_named(r, 2, keep), ...
                'a', SRNNCellTypes.trim_named(a, 3, keep), ...
                'b', SRNNCellTypes.trim_named(b, 3, keep), ...
                'g', SRNNCellTypes.trim_named(g, 3, keep), ...
                'br', SRNNCellTypes.trim_named(br, 2, keep));
        end
    end

    %% Plotting and spectrum inspection
    methods
        function [fig_handle, ax_handles] = plot(obj, varargin)
            if ~obj.has_run || isempty(obj.plot_data)
                error('SRNNCellTypes:NotRun', ...
                    'Run the model with store_decimated_state=true before plotting.');
            end
            range = obj.T_plot;
            for k = 1:2:numel(varargin)
                if strcmpi(varargin{k}, 'T_plot'), range = varargin{k + 1}; end
            end
            if isempty(range), range = obj.T_range; end

            p = obj.plot_data;
            has_a = any(obj.n_a > 0);
            has_b = any(obj.n_b > 0);
            has_g = any(obj.n_g > 0);
            has_lya = ~strcmpi(obj.lya_method, 'none') && ~isempty(obj.lya_results);
            n_panels = 4 + has_a + has_b + has_g + has_lya;
            fig_handle = figure('Name', 'SRNNCellTypes time series');
            tiledlayout(n_panels, 1);
            ax_handles = gobjects(n_panels, 1);
            panel = 0;

            panel = panel + 1; ax_handles(panel) = nexttile;
            SRNNCellTypes.plot_named_series( ...
                p.t, p.u, obj.cell_type_names, 'External input', true);
            panel = panel + 1; ax_handles(panel) = nexttile;
            SRNNCellTypes.plot_named_series( ...
                p.t, p.x, obj.cell_type_names, 'Dendritic state', false);
            panel = panel + 1; ax_handles(panel) = nexttile;
            SRNNCellTypes.plot_named_series( ...
                p.t, p.r, obj.cell_type_names, 'Firing rate', false);
            panel = panel + 1; ax_handles(panel) = nexttile;
            SRNNCellTypes.plot_named_series( ...
                p.t, p.br, obj.cell_type_names, 'Synaptic output', false);
            if has_a
                panel = panel + 1; ax_handles(panel) = nexttile;
                a_collapsed = SRNNCellTypes.collapse_named(p.a, obj.cell_type_names, 'sum');
                SRNNCellTypes.plot_named_series( ...
                    p.t, a_collapsed, obj.cell_type_names, 'SFA (sum)', false);
            end
            if has_b
                panel = panel + 1; ax_handles(panel) = nexttile;
                b_collapsed = SRNNCellTypes.collapse_named(p.b, obj.cell_type_names, 'prod');
                SRNNCellTypes.plot_named_series( ...
                    p.t, b_collapsed, obj.cell_type_names, 'STD (product)', false);
            end
            if has_g
                panel = panel + 1; ax_handles(panel) = nexttile;
                g_collapsed = SRNNCellTypes.collapse_named(p.g, obj.cell_type_names, 'prod');
                SRNNCellTypes.plot_named_series( ...
                    p.t, g_collapsed, obj.cell_type_names, 'STF (product)', false);
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
            %PLOT_CELLTYPES Plot each cell type in a separate subplot column.
            %
            % Rows match plot(): external input, dendritic state, firing rate,
            % synaptic output, optional SFA, STD, and STF, and optional Lyapunov
            % output. The network-level Lyapunov trace is repeated in every
            % cell-type column to preserve the aligned grid.
            %
            % Usage:
            %   model.plot_celltypes()
            %   model.plot_celltypes('T_plot', [10 40])

            if ~obj.has_run || isempty(obj.plot_data)
                error('SRNNCellTypes:NotRun', ...
                    'Run the model with store_decimated_state=true before plotting.');
            end

            range = obj.T_plot;
            for k = 1:2:numel(varargin)
                if strcmpi(varargin{k}, 'T_plot'), range = varargin{k + 1}; end
            end
            if isempty(range), range = obj.T_range; end

            p = obj.plot_data;
            has_a = any(obj.n_a > 0);
            has_b = any(obj.n_b > 0);
            has_g = any(obj.n_g > 0);
            has_lya = ~strcmpi(obj.lya_method, 'none') && ~isempty(obj.lya_results);

            series = {p.u, p.x, p.r, p.br};
            row_labels = {'External input', 'Dendritic state', ...
                'Firing rate', 'Synaptic output'};
            empty_labels = {'', '', '', ''};
            b_row = [];
            g_row = [];
            if has_a
                series{end + 1} = SRNNCellTypes.collapse_named( ...
                    p.a, obj.cell_type_names, 'sum');
                row_labels{end + 1} = 'SFA (sum)';
                empty_labels{end + 1} = 'No SFA';
            end
            if has_b
                series{end + 1} = SRNNCellTypes.collapse_named( ...
                    p.b, obj.cell_type_names, 'prod');
                b_row = numel(series);
                row_labels{end + 1} = 'STD (product)';
                empty_labels{end + 1} = 'No STD';
            end
            if has_g
                series{end + 1} = SRNNCellTypes.collapse_named( ...
                    p.g, obj.cell_type_names, 'prod');
                g_row = numel(series);
                row_labels{end + 1} = 'STF (product)';
                empty_labels{end + 1} = 'No STF';
            end

            n_series_rows = numel(series);
            n_rows = n_series_rows + has_lya;
            C = obj.n_cellTypes;
            colors = lines(C);
            fig_handle = figure('Name', 'SRNNCellTypes by cell type');
            layout = tiledlayout(fig_handle, n_rows, C, ...
                'TileSpacing', 'compact', 'Padding', 'compact');
            ax_handles = gobjects(n_rows, C);

            for row = 1:n_series_rows
                for q = 1:C
                    tile_index = (row - 1) * C + q;
                    ax = nexttile(layout, tile_index);
                    ax_handles(row, q) = ax;
                    name = obj.cell_type_names{q};
                    values = series{row}.(name);

                    % Inactive synaptic mechanisms are stored as constant-one
                    % readouts. Show them as disabled instead of plotting them.
                    inactive_b = ~isempty(b_row) && row == b_row && obj.n_b(q) == 0;
                    inactive_g = ~isempty(g_row) && row == g_row && obj.n_g(q) == 0;
                    if inactive_b || inactive_g
                        values = [];
                    end

                    if isempty(values)
                        text(ax, 0.5, 0.5, empty_labels{row}, ...
                            'Units', 'normalized', ...
                            'HorizontalAlignment', 'center', ...
                            'VerticalAlignment', 'middle', ...
                            'Color', [0.4 0.4 0.4]);
                    else
                        SRNNCellTypes.plot_celltype_lines( ...
                            ax, p.t, values, colors(q, :));
                    end
                    if row == 1
                        title(ax, name, 'Interpreter', 'none');
                    end
                    if q == 1
                        ylabel(ax, row_labels{row});
                    end
                    if row < n_rows
                        set(ax, 'XTickLabel', []);
                    end
                end
            end

            if has_lya
                row = n_rows;
                for q = 1:C
                    tile_index = (row - 1) * C + q;
                    ax = nexttile(layout, tile_index);
                    ax_handles(row, q) = ax;
                    if isfield(obj.lya_results, 'local_lya')
                        plot(ax, obj.lya_results.t_lya, ...
                            obj.lya_results.local_lya, 'k');
                        if q == 1, ylabel(ax, 'Local LLE'); end
                    elseif isfield(obj.lya_results, 'local_LE_spectrum_t')
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
                error('SRNNCellTypes:NoStateData', ...
                    'Set store_full_state=true before running.');
            end
            sample_indices = round((times_sec - obj.t_out(1)) * obj.fs) + 1;
            sample_indices = unique(max(1, min(sample_indices, size(obj.S_out, 1))));
            n_plots = numel(sample_indices);
            fig_handle = figure('Name', 'SRNNCellTypes Jacobian eigenvalues');
            ax_handles = gobjects(n_plots, 1);
            n_cols = ceil(sqrt(n_plots));
            n_rows = ceil(n_plots / n_cols);
            for k = 1:n_plots
                ax_handles(k) = subplot(n_rows, n_cols, k);
                J = SRNNCellTypes.compute_Jacobian_fast( ...
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
                error('SRNNCellTypes:NotBuilt', 'Call build() first.');
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

        function [fig_handle, ax_handle] = plot_weight_histogram(obj, bin_edges)
            if ~obj.is_built
                error('SRNNCellTypes:NotBuilt', 'Call build() first.');
            end
            values = nonzeros(obj.W);
            if nargin < 2 || isempty(bin_edges)
                radius = max(abs(values));
                if isempty(radius) || radius == 0, radius = 1; end
                bin_edges = linspace(-radius, radius, 51);
            end
            colors = lines(obj.n_cellTypes);
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
        function layout = make_state_layout(n, counts, n_a, n_b, n_g)
            if nargin < 5, n_g = zeros(size(counts)); end
            C = numel(counts);
            layout = struct('a', {cell(1, C)}, 'b', {cell(1, C)}, ...
                'g', {cell(1, C)}, 'x', []);
            cursor = 0;
            for q = 1:C
                len = counts(q) * n_a(q);
                layout.a{q} = cursor + (1:len);
                cursor = cursor + len;
            end
            for q = 1:C
                len = counts(q) * n_b(q);
                layout.b{q} = cursor + (1:len);
                cursor = cursor + len;
            end
            for q = 1:C
                len = counts(q) * n_g(q);
                layout.g{q} = cursor + (1:len);
                cursor = cursor + len;
            end
            layout.x = cursor + (1:n);
            layout.N = cursor + n;
        end

        function S0 = initialize_state(params)
            parts = cell(1, 3 * params.n_cellTypes + 1);
            for q = 1:params.n_cellTypes
                parts{q} = zeros(params.n_per_type(q) * params.n_a(q), 1);
                parts{params.n_cellTypes + q} = ...
                    ones(params.n_per_type(q) * params.n_b(q), 1);
                parts{2 * params.n_cellTypes + q} = ...
                    ones(params.n_per_type(q) * params.n_g(q), 1);
            end
            parts{end} = params.x0_std .* randn(params.n, 1);
            S0 = vertcat(parts{:});
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
                error('SRNNCellTypes:InvalidStimulus', ...
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
                error('SRNNCellTypes:InvalidStimulus', ...
                    'intrinsic_drive must be scalar or n-by-1.');
            end
            u_ex = u_ex + drive;
        end

        function dS = dynamics_fast(t, S, params)
            u = params.u_interpolant(t)';
            C = params.n_cellTypes;
            layout = params.state_layout;
            a = cell(1, C); b_state = cell(1, C); g_state = cell(1, C);
            x = S(layout.x);
            x_eff = x;
            depression = ones(params.n, 1);
            facilitation = ones(params.n, 1);

            for q = 1:C
                idx = params.type_indices{q};
                nq = params.n_per_type(q);
                if params.n_a(q) > 0
                    a{q} = reshape(S(layout.a{q}), nq, params.n_a(q));
                    x_eff(idx) = x_eff(idx) - params.c(q) .* sum(a{q}, 2);
                else
                    a{q} = zeros(nq, 0);
                end
                if params.n_b(q) > 0
                    b_state{q} = reshape(S(layout.b{q}), nq, params.n_b(q));
                    depression(idx) = prod(b_state{q}, 2);
                else
                    b_state{q} = ones(nq, 0);
                end
                if params.n_g(q) > 0
                    g_state{q} = reshape(S(layout.g{q}), nq, params.n_g(q));
                    facilitation(idx) = prod(g_state{q}, 2);
                else
                    g_state{q} = ones(nq, 0);
                end
            end

            rate = params.activation_function(x_eff);
            depression_syn = SRNNCellTypes.apply_std_zero_floor(depression, params);
            derivatives = cell(1, 3 * C + 1);
            for q = 1:C
                idx = params.type_indices{q};
                if params.n_a(q) > 0
                    derivatives{q} = ((rate(idx) - a{q}) ./ params.tau_a{q});
                else
                    derivatives{q} = [];
                end
                if params.n_b(q) > 0
                    derivatives{C + q} = (1 - b_state{q}) ./ params.tau_b_rec{q} ...
                        - (b_state{q} .* rate(idx)) ./ params.tau_b_rel(q);
                else
                    derivatives{C + q} = [];
                end
                if params.n_g(q) > 0
                    derivatives{2 * C + q} = ...
                        (1 - g_state{q}) ./ params.tau_g_dec{q} ...
                        + ((params.G(q) - g_state{q}) .* rate(idx)) ...
                        ./ params.tau_g_fac(q);
                else
                    derivatives{2 * C + q} = [];
                end
            end
            derivatives{end} = (-x + params.W * ...
                (facilitation .* depression_syn .* rate) + u) ...
                ./ params.tau_d;
            for k = 1:3*C
                derivatives{k} = derivatives{k}(:);
            end
            dS = vertcat(derivatives{:});
        end

        function b_used = apply_std_zero_floor(b, params)
            b_used = b;
            if ~params.std_zero_floor, return; end
            for q = 1:params.n_cellTypes
                if params.n_b(q) == 0, continue; end
                p_min = prod(params.tau_b_rel(q) ./ ...
                    (params.tau_b_rec{q} + params.tau_b_rel(q)));
                idx = params.type_indices{q};
                b_used(idx, :) = (b(idx, :) - p_min) ./ (1 - p_min);
            end
        end

        function J = compute_Jacobian_fast(S, params)
            C = params.n_cellTypes;
            layout = params.state_layout;
            n = params.n;
            a = cell(1, C); b_state = cell(1, C); g_state = cell(1, C);
            P = cell(1, C); Q = cell(1, C);
            x = S(layout.x);
            x_eff = x;
            depression = ones(n, 1);
            facilitation = ones(n, 1);

            for q = 1:C
                idx = params.type_indices{q}; nq = params.n_per_type(q);
                if params.n_a(q) > 0
                    a{q} = reshape(S(layout.a{q}), nq, params.n_a(q));
                    x_eff(idx) = x_eff(idx) - params.c(q) .* sum(a{q}, 2);
                else
                    a{q} = zeros(nq, 0);
                end
                if params.n_b(q) > 0
                    b_state{q} = reshape(S(layout.b{q}), nq, params.n_b(q));
                    P{q} = prod(b_state{q}, 2);
                    depression(idx) = P{q};
                else
                    b_state{q} = ones(nq, 0); P{q} = ones(nq, 1);
                end
                if params.n_g(q) > 0
                    g_state{q} = reshape(S(layout.g{q}), nq, params.n_g(q));
                    Q{q} = prod(g_state{q}, 2);
                    facilitation(idx) = Q{q};
                else
                    g_state{q} = ones(nq, 0); Q{q} = ones(nq, 1);
                end
            end

            rate = params.activation_function(x_eff);
            rate_prime = params.activation_function_derivative(x_eff);
            depression_syn = SRNNCellTypes.apply_std_zero_floor(depression, params);
            synaptic_modulation = facilitation .* depression_syn;
            J = sparse(layout.N, layout.N);
            row_x = layout.x;

            for q = 1:C
                idx = params.type_indices{q}; nq = params.n_per_type(q);
                na = params.n_a(q); nb = params.n_b(q); ng = params.n_g(q);
                row_a = layout.a{q}; row_b = layout.b{q}; row_g = layout.g{q};

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

                if nb > 0
                    if na > 0
                        coeff = b_state{q} .* (params.c(q) .* rate_prime(idx) ...
                            ./ params.tau_b_rel(q));
                        stack = sparse(1:numel(row_b), repmat((1:nq)', nb, 1), ...
                            coeff(:), numel(row_b), nq);
                        J(row_b, row_a) = kron(sparse(ones(1, na)), stack);
                    end
                    diagonal = kron(-1 ./ params.tau_b_rec{q}(:), ones(nq, 1)) ...
                        + repmat(-rate(idx) ./ params.tau_b_rel(q), nb, 1);
                    J(row_b, row_b) = spdiags(diagonal, 0, numel(row_b), numel(row_b));
                    vals = -(b_state{q}(:) .* repmat(rate_prime(idx), nb, 1)) ...
                        ./ params.tau_b_rel(q);
                    J(row_b, row_x) = sparse(1:numel(row_b), ...
                        repmat(idx(:), nb, 1), vals, numel(row_b), n);
                end

                if ng > 0
                    if na > 0
                        coeff = -(params.G(q) - g_state{q}) .* ...
                            (params.c(q) .* rate_prime(idx) ./ params.tau_g_fac(q));
                        stack = sparse(1:numel(row_g), repmat((1:nq)', ng, 1), ...
                            coeff(:), numel(row_g), nq);
                        J(row_g, row_a) = kron(sparse(ones(1, na)), stack);
                    end
                    diagonal = kron(-1 ./ params.tau_g_dec{q}(:), ones(nq, 1)) ...
                        + repmat(-rate(idx) ./ params.tau_g_fac(q), ng, 1);
                    J(row_g, row_g) = spdiags(diagonal, 0, numel(row_g), numel(row_g));
                    vals = ((params.G(q) - g_state{q}(:)) .* ...
                        repmat(rate_prime(idx), ng, 1)) ./ params.tau_g_fac(q);
                    J(row_g, row_x) = sparse(1:numel(row_g), ...
                        repmat(idx(:), ng, 1), vals, numel(row_g), n);
                end

                if na > 0
                    replicate = kron(ones(1, na), speye(nq));
                    block = -params.c(q) .* params.W(:, idx) * spdiags( ...
                        synaptic_modulation(idx) .* rate_prime(idx), 0, nq, nq);
                    J(row_x, row_a) = (block * replicate) ./ params.tau_d;
                end
                if nb > 0
                    product_derivative = zeros(nq, nb);
                    for m = 1:nb
                        other = [1:m-1, m+1:nb];
                        if isempty(other)
                            product_derivative(:, m) = 1;
                        else
                            product_derivative(:, m) = prod(b_state{q}(:, other), 2);
                        end
                    end
                    gain = 1;
                    if params.std_zero_floor
                        p_min = prod(params.tau_b_rel(q) ./ ...
                            (params.tau_b_rec{q} + params.tau_b_rel(q)));
                        gain = 1 / (1 - p_min);
                    end
                    coeff = gain .* facilitation(idx) .* rate(idx) .* product_derivative;
                    D = sparse(repmat((1:nq)', nb, 1), 1:numel(row_b), ...
                        coeff(:), nq, numel(row_b));
                    J(row_x, row_b) = (params.W(:, idx) * D) ./ params.tau_d;
                end
                if ng > 0
                    product_derivative = zeros(nq, ng);
                    for k = 1:ng
                        other = [1:k-1, k+1:ng];
                        if isempty(other)
                            product_derivative(:, k) = 1;
                        else
                            product_derivative(:, k) = prod(g_state{q}(:, other), 2);
                        end
                    end
                    coeff = depression_syn(idx) .* rate(idx) .* product_derivative;
                    D = sparse(repmat((1:nq)', ng, 1), 1:numel(row_g), ...
                        coeff(:), nq, numel(row_g));
                    J(row_x, row_g) = (params.W(:, idx) * D) ./ params.tau_d;
                end
            end

            J(row_x, row_x) = spdiags(-ones(n, 1) ./ params.tau_d, 0, n, n) ...
                + params.W * spdiags(synaptic_modulation .* rate_prime, 0, n, n) ...
                ./ params.tau_d;
        end

        function J_array = compute_Jacobian_at_indices(S_out, sample_indices, params)
            J_array = zeros(params.N_sys_eqs, params.N_sys_eqs, numel(sample_indices));
            for k = 1:numel(sample_indices)
                J_array(:, :, k) = full(SRNNCellTypes.compute_Jacobian_fast( ...
                    S_out(sample_indices(k), :)', params));
            end
        end
    end

    %% State conversion and plotting helpers
    methods (Static)
        function [x_named, a_named, b_named, r_named, br_named, g_named] = ...
                unpack_and_compute_states(S_out, params)
            nt = size(S_out, 1);
            layout = params.state_layout;
            x_all = S_out(:, layout.x)';
            x_eff = x_all;
            depression = ones(params.n, nt);
            facilitation = ones(params.n, nt);
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
                if params.n_b(q) > 0
                    bq = reshape(S_out(:, layout.b{q})', nq, params.n_b(q), nt);
                    depression(idx, :) = reshape(prod(bq, 2), nq, nt);
                    b_named.(name) = bq;
                else
                    b_named.(name) = ones(nq, 1, nt);
                end
                if params.n_g(q) > 0
                    gq = reshape(S_out(:, layout.g{q})', nq, params.n_g(q), nt);
                    facilitation(idx, :) = reshape(prod(gq, 2), nq, nt);
                    g_named.(name) = gq;
                else
                    g_named.(name) = ones(nq, 1, nt);
                end
            end
            rate = params.activation_function(x_eff);
            br = facilitation .* ...
                SRNNCellTypes.apply_std_zero_floor(depression, params) .* rate;
            x_named = SRNNCellTypes.split_by_type(x_all, params, 1);
            r_named = SRNNCellTypes.split_by_type(rate, params, 1);
            br_named = SRNNCellTypes.split_by_type(br, params, 1);
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
            colors = lines(numel(names));
            hold on;
            handles = gobjects(0); labels = {};
            for q = 1:numel(names)
                values = data.(names{q});
                if isempty(values), continue; end
                SRNNCellTypes.plot_celltype_lines( ...
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
            if a < 0 || a > 1
                error('SRNNCellTypes:InvalidActivation', 'a must be between 0 and 1.');
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
            if a < 0 || a > 1
                error('SRNNCellTypes:InvalidActivation', 'a must be between 0 and 1.');
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
            value = SRNNCellTypes.logisticSigmoid(x, c);
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
                    error('SRNNCellTypes:UnknownLyapunovMethod', ...
                        'Unknown Lyapunov method: %s', method);
            end
            factor = lya_dt / dt;
            if abs(factor - round(factor)) > 1e-11 || factor < 3
                error('SRNNCellTypes:InvalidLyapunovStep', ...
                    'Lyapunov interval must be an integer multiple of at least 3*dt.');
            end

            if strcmpi(method, 'benettin')
                [LLE, local, finite, times] = SRNNCellTypes.benettin_algorithm_internal( ...
                    S_out, t_out, dt, fs, 1e-3, interval, lya_dt, ...
                    params, opts, rhs, ode_solver);
                results.LLE = LLE;
                results.local_lya = local;
                results.finite_lya = finite;
                results.t_lya = times;
            else
                [spectrum, local, finite, times] = ...
                    SRNNCellTypes.lyapunov_spectrum_qr_internal( ...
                    S_out, t_out, lya_dt, params, ode_solver, opts, fs);
                [spectrum, order] = sort(real(spectrum), 'descend');
                results.LE_spectrum = spectrum;
                results.local_LE_spectrum_t = local(:, order);
                results.finite_LE_spectrum_t = finite(:, order);
                results.t_lya = times;
                results.sort_idx = order;
                results.params.N_sys_eqs = params.N_sys_eqs;
                fprintf('Lyapunov Dimension: %.2f\n', ...
                    SRNNCellTypes.compute_kaplan_yorke_dimension_internal(spectrum));
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
                    warning('SRNNCellTypes:LyapunovDiverged', ...
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
                variational = @(tt, psi) SRNNCellTypes.variational_equations( ...
                    tt, psi, interpolants, N, params);
                [~, solution] = ode_solver(variational, [t_start, t_end], ...
                    psi0, variational_options);
                evolved = reshape(solution(end, :)', N, N);
                if any(~isfinite(evolved(:)))
                    warning('SRNNCellTypes:LyapunovDiverged', ...
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
            J = SRNNCellTypes.compute_Jacobian_fast(state, params);
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
