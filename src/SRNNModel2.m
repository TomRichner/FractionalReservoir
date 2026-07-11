classdef SRNNModel2 < SRNNModelBase
    % SRNNMODEL2 Excitatory/Inhibitory rate reservoir (concrete SRNNModelBase subclass).
    %
    % Implements the E/I model: RMT connectivity, the S = [a_E;a_I;b_E;b_I;x]
    % state layout, and the dynamics_fast / compute_Jacobian_fast numeric kernels
    % (plugged into the base via the eval_dynamics / eval_jacobian seam). The
    % type-agnostic machinery (build/run/lyapunov/decimation/plotting, activation
    % library) is inherited from SRNNModelBase.
    %
    % E/I-SPECIFIC MEMBERS kept here: RMT connectivity (build_network) and stimulus
    % (build_stimulus/generate_stimulus/generate_external_input), state pack/unpack
    % (initialize_state/unpack_and_compute_states), the RHS/Jacobian kernels
    % (dynamics_fast/compute_Jacobian_fast/apply_std_zero_floor), get_params/validate,
    % and the E/I time-series plotters. The Lyapunov algorithms, decimation, activation
    % functions, and eigenvalue plotting now live in SRNNModelBase.
    %
    % Usage:
    %   model = SRNNModel2('n_a_E', 3, 'n_b_E', 1);
    %   model.build();
    %   model.run();
    %   model.plot();
    %
    % See also: SRNN_ESN_reservoir


    %% Network Architecture Properties
    properties
        f = 0.5                     % Fraction of excitatory neurons

        % RMT tilde-notation parameters (Harris 2023)
        mu_E_tilde                  % Normalized excitatory mean (default: 1/(alpha*sqrt(n)))
        mu_I_tilde                  % Normalized inhibitory mean (default: -1/(alpha*sqrt(n)))
        sigma_E_tilde               % Normalized excitatory std dev (default: 1/(alpha*sqrt(n)))
        sigma_I_tilde               % Normalized inhibitory std dev (default: 1/(alpha*sqrt(n)))
        E_W = 0                     % Mean offset: added to both mu_E_tilde and mu_I_tilde
        zrs_mode = 'none'           % ZRS mode: 'none', 'ZRS', 'SZRS', 'Partial_SZRS'

    end

    %% Spike-Frequency Adaptation (SFA) Properties
    properties
        n_a_E = 0                   % Number of adaptation timescales for E neurons
        n_a_I = 0                   % Number of adaptation timescales for I neurons
        tau_a_E                     % Adaptation time constants for E neurons (1 x n_a_E)
        tau_a_I                     % Adaptation time constants for I neurons (1 x n_a_I)
        c_E = 0.15/3                 % Adaptation scaling for E neurons
        c_I = 0.15/3                  % Adaptation scaling for I neurons
    end

    %% Short-Term Depression (STD) Properties
    properties
        n_b_E = 0                   % Number of STD timescales for E neurons
        n_b_I = 0                   % Number of STD timescales for I neurons
        tau_b_E_rec = 1             % STD recovery time constants for E neurons (1 x n_b_E vector)
        tau_b_E_rel = 0.25          % STD release time constant for E neurons (scalar, shared across timescales)
        tau_b_I_rec = 1             % STD recovery time constants for I neurons (1 x n_b_I vector)
        tau_b_I_rel = 0.25          % STD release time constant for I neurons (scalar, shared across timescales)
        std_zero_floor = false      % If true, rescale synaptic prod_m(b_m) -> (P-P_min)/(1-P_min) so (prod_m b_m)*r reaches 0 at r=1
    end

    %% RMT Dependent Properties (E/I-specific; type-agnostic props live in SRNNModelBase)
    properties (Dependent)
        default_val         % Normalization factor F = 1/sqrt(N*alpha*(2-alpha)) (Harris 2023)
        mu_se               % Sparse excitatory mean
        mu_si               % Sparse inhibitory mean
        sigma_se            % Sparse excitatory std dev
        sigma_si            % Sparse inhibitory std dev
        R                   % Theoretical spectral radius (Harris 2023 Eq 18)
        n_E                 % Number of excitatory neurons = round(f*n)
        n_I                 % Number of inhibitory neurons = n - n_E
        E_indices           % Indices of E neurons = 1:n_E
        I_indices           % Indices of I neurons = (n_E+1):n
        N_sys_eqs           % Total state dimension
    end

    %% Constructor
    methods
        function obj = SRNNModel2(varargin)
            % SRNNMODEL2 Constructor. Forwards to the shared base constructor,
            % which applies defaults (via the set_defaults override below), parses
            % name-value pairs, and derives plot_deci.
            %
            % Usage:
            %   model = SRNNModel2()  % All defaults (n=300, indegree=100)
            %   model = SRNNModel2('n', 200, 'level_of_chaos', 2.0)
            %   model = SRNNModel2('n_a_E', 3, 'n_b_E', 1)
            obj@SRNNModelBase(varargin{:});
        end
    end

    %% Subclass hooks (implement the SRNNModelBase abstract contract)
    methods (Access = protected)
        function set_defaults(obj)
            % SET_DEFAULTS Extend the base (type-agnostic) defaults with the
            % E/I-specific input_config fields. Net result is identical to the
            % original single set_defaults.
            set_defaults@SRNNModelBase(obj);
            obj.input_config.step_density_E = 0.15;   % Fraction of E neurons receiving input
            obj.input_config.step_density_I = 0;      % Fraction of I neurons receiving input (0 = E only)
        end

        function dS_dt = eval_dynamics(obj, t, S, params)
            % EVAL_DYNAMICS RHS seam: plug the E/I kernel into the base machinery.
            dS_dt = SRNNModel2.dynamics_fast(t, S, params);
        end

        function J = eval_jacobian(obj, S, params)
            % EVAL_JACOBIAN Jacobian seam: plug the E/I kernel into the base machinery.
            J = SRNNModel2.compute_Jacobian_fast(S, params);
        end
    end

    %% Dependent Property Getters
    methods
        function val = get.default_val(obj)
            % DEFAULT_VAL Normalization factor F = 1/sqrt(N*alpha*(2-alpha))
            % Scaling factor which yields R=1 when all tilde parameters are equal.
            % See parameter_table.md for derivation (Harris 2023).
            val = 1 / sqrt(obj.n * obj.alpha * (2 - obj.alpha));
        end

        function val = get.mu_se(obj)
            if isempty(obj.mu_E_tilde)
                val = NaN;
            else
                val = obj.alpha * (obj.mu_E_tilde + obj.E_W);
            end
        end

        function val = get.mu_si(obj)
            if isempty(obj.mu_I_tilde)
                val = NaN;
            else
                val = obj.alpha * (obj.mu_I_tilde + obj.E_W);
            end
        end

        function val = get.sigma_se(obj)
            if isempty(obj.sigma_E_tilde) || isempty(obj.mu_E_tilde)
                val = NaN;
            else
                mu_eff = obj.mu_E_tilde + obj.E_W;
                val = sqrt(obj.alpha * (1 - obj.alpha) * mu_eff^2 + obj.alpha * obj.sigma_E_tilde^2);
            end
        end

        function val = get.sigma_si(obj)
            if isempty(obj.sigma_I_tilde) || isempty(obj.mu_I_tilde)
                val = NaN;
            else
                mu_eff = obj.mu_I_tilde + obj.E_W;
                val = sqrt(obj.alpha * (1 - obj.alpha) * mu_eff^2 + obj.alpha * obj.sigma_I_tilde^2);
            end
        end

        function val = get.R(obj)
            if isnan(obj.sigma_se) || isnan(obj.sigma_si)
                val = NaN;
            else
                val = sqrt(obj.n * (obj.f * obj.sigma_se^2 + (1 - obj.f) * obj.sigma_si^2)) * obj.level_of_chaos;
            end
        end

        function val = get.n_E(obj)
            val = round(obj.f * obj.n);
        end

        function val = get.n_I(obj)
            val = obj.n - obj.n_E;
        end

        function val = get.E_indices(obj)
            val = 1:obj.n_E;
        end

        function val = get.I_indices(obj)
            val = (obj.n_E + 1):obj.n;
        end

        function val = get.N_sys_eqs(obj)
            nE = obj.n_E;
            nI = obj.n_I;
            val = nE * obj.n_a_E + nI * obj.n_a_I + nE * obj.n_b_E + nI * obj.n_b_I + obj.n;
        end
    end

    %% Public Methods
    methods
        function [fig_handle, ax_handles] = plot(obj, varargin)
            % PLOT Generate time series plots for SRNN simulation
            %
            % Usage:
            %   model.plot()
            %   model.plot('T_plot', [10, 40])
            %   [fig, axes] = model.plot()

            if ~obj.has_run
                error('SRNNModel:NotRun', 'Model must be run before plotting. Call run() first.');
            end

            if isempty(obj.plot_data)
                error('SRNNModel:NoPlotData', 'Plot data not available. Set store_decimated_state=true.');
            end

            % Parse optional arguments
            T_plot_arg = obj.T_plot;
            for i = 1:2:length(varargin)
                if strcmpi(varargin{i}, 'T_plot')
                    T_plot_arg = varargin{i+1};
                end
            end

            if isempty(T_plot_arg)
                T_plot_arg = obj.T_range;
            end

            params = obj.cached_params;

            [fig_handle, ax_handles] = SRNNModel2.plot_SRNN_tseries(obj.plot_data.t, obj.plot_data.u, obj.plot_data.x, obj.plot_data.r, obj.plot_data.a, obj.plot_data.b, obj.plot_data.br, params, obj.lya_results, obj.lya_method, T_plot_arg);
        end

        function [fig_handle, ax_handles] = plot_W_spectrum(obj)
            % PLOT_W_SPECTRUM Plot eigenvalue spectra of -I+W and the LTI Jacobian
            %
            % This method creates a 2-panel figure showing:
            %   Left:  Eigenvalues of (-I + W) (unscaled Jacobian)
            %   Right: Eigenvalues of (-I + W)/tau_d (the LTI Jacobian)
            %
            % The theoretical spectral radius R (from Harris 2023 Eq 18) is shown
            % as a circle centered at -1. For the Jacobian, the circle is
            % centered at -1/tau_d and scaled by 1/tau_d.
            %
            % Usage:
            %   model.build();
            %   model.plot_W_spectrum();

            if ~obj.is_built
                error('SRNNModel:NotBuilt', 'Model must be built before plotting spectrum. Call build() first.');
            end

            % Compute eigenvalues
            J_unscaled = -eye(obj.n) + obj.W;
            eigs_unscaled = eig(J_unscaled);
            J_lti = J_unscaled / obj.tau_d;
            eigs_J = eig(J_lti);

            % Theoretical predictions
            R_W = obj.R;  % Already scaled by level_of_chaos
            outlier_threshold = 1.04;

            % Create figure
            fig_handle = figure('Position', [100, 300, 900, 400]);
            ax_handles = gobjects(2, 1);

            %% Left panel: (-I + W) eigenvalues (unscaled Jacobian)
            ax_handles(1) = subplot(1, 2, 1);
            center_unscaled = -1;  % Shifted by -I
            obj.plot_spectrum_helper(ax_handles(1), eigs_unscaled, R_W, center_unscaled, outlier_threshold);
            title(ax_handles(1), sprintf('-I + W eigenvalues (R = %.2f)', R_W), 'FontWeight', 'bold');
            xlabel(ax_handles(1), 'Re(\lambda)');
            ylabel(ax_handles(1), 'Im(\lambda)');

            % Add stability line at Re = 0
            hold(ax_handles(1), 'on');
            yl = ylim(ax_handles(1));
            plot(ax_handles(1), [0, 0], yl, 'r--', 'LineWidth', 1.5);
            hold(ax_handles(1), 'off');

            %% Right panel: LTI Jacobian eigenvalues (-I + W)/tau_d
            ax_handles(2) = subplot(1, 2, 2);
            % For the Jacobian: center shifts to -1/tau_d, radius scales by 1/tau_d
            R_J = R_W / obj.tau_d;
            center_J = -1 / obj.tau_d;
            obj.plot_spectrum_helper(ax_handles(2), eigs_J, R_J, center_J, outlier_threshold);
            title(ax_handles(2), sprintf('(-I + W)/\\tau_d eigenvalues (R_J = %.2f)', R_J), 'FontWeight', 'bold');
            xlabel(ax_handles(2), 'Re(\lambda)');
            ylabel(ax_handles(2), 'Im(\lambda)');

            % Add stability line at Re = 0
            hold(ax_handles(2), 'on');
            yl = ylim(ax_handles(2));
            plot(ax_handles(2), [0, 0], yl, 'r--', 'LineWidth', 1.5);
            hold(ax_handles(2), 'off');
        end

        function [fig_handle, ax_handle] = plot_weight_histogram(obj, bin_edges, show_legend)
            % PLOT_WEIGHT_HISTOGRAM Histogram of W entries split by E/I presynaptic population
            %
            % Overlays the excitatory (columns 1:n_E, red) and inhibitory
            % (columns n_E+1:end, blue) nonzero-weight distributions, with
            % population-mean markers below the x-axis. The means are measured
            % from the plotted (scaled) entries, so they stay correct under
            % level_of_chaos / E_W / rescale_by_abscissa. Values outside the
            % bin range are clipped to the first/last bin so large entries
            % surface as edge spikes rather than vanishing.
            %
            % Adapted from RMT.plot_weight_histogram (ConnectivityAdaptation repo).
            %
            % Usage:
            %   model.build();
            %   model.plot_weight_histogram();
            %   model.plot_weight_histogram(linspace(-0.2, 0.2, 51));
            %   model.plot_weight_histogram([], false);   % omit legend

            if ~obj.is_built
                error('SRNNModel:NotBuilt', 'Model must be built before plotting weights. Call build() first.');
            end
            if nargin < 3 || isempty(show_legend), show_legend = true; end

            W_full = full(obj.W);
            if nargin < 2 || isempty(bin_edges)
                r = max(abs(W_full(:)));
                if r == 0, r = 1; end
                bin_edges = linspace(-r, r, 51);
            end
            edge_lo = bin_edges(1);
            edge_hi = bin_edges(end);

            % Split by presynaptic population (columns), nonzero entries only
            W_E = W_full(:, 1:obj.n_E);       W_E = W_E(W_E ~= 0);
            W_I = W_full(:, obj.n_E+1:end);   W_I = W_I(W_I ~= 0);

            fig_handle = figure('Position', [100, 300, 560, 420], 'Color', 'white');
            ax_handle = gca;
            hold(ax_handle, 'on');
            histogram(ax_handle, min(max(W_E, edge_lo), edge_hi), bin_edges, ...
                'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
            histogram(ax_handle, min(max(W_I, edge_lo), edge_hi), bin_edges, ...
                'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.6);

            if show_legend
                legend(ax_handle, 'E', 'I', 'Location', 'northeast');
                legend(ax_handle, 'boxoff');
            end

            % Population-mean markers just below the x-axis
            y_bottom = ax_handle.YLim(1);
            if ~isempty(W_E)
                text(ax_handle, mean(W_E), y_bottom, '$\mu_E$', 'Interpreter', 'latex', ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                    'Color', [0.8 0.2 0.2], 'FontSize', 11, 'Clipping', 'on');
            end
            if ~isempty(W_I)
                text(ax_handle, mean(W_I), y_bottom, '$\mu_I$', 'Interpreter', 'latex', ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                    'Color', [0.2 0.4 0.8], 'FontSize', 11, 'Clipping', 'on');
            end
            hold(ax_handle, 'off');

            xlabel(ax_handle, 'Weight');
            ylabel(ax_handle, 'Count');
            title(ax_handle, 'Synaptic weight distribution (nonzero entries)', 'FontWeight', 'bold');
            box(ax_handle, 'off');
            set(ax_handle, 'Color', 'none');
        end

        function params = get_params(obj)
            % GET_PARAMS Return params struct for compatibility with existing functions
            %
            % This method creates a struct containing all parameters needed
            % by functions like SRNN_reservoir, dynamics_fast, etc.

            params = struct();

            % Network architecture
            params.n = obj.n;
            params.f = obj.f;
            params.indegree = obj.indegree;
            params.alpha = obj.alpha;
            params.level_of_chaos = obj.level_of_chaos;

            % RMT parameters
            params.mu_E_tilde = obj.mu_E_tilde;
            params.mu_I_tilde = obj.mu_I_tilde;
            params.sigma_E_tilde = obj.sigma_E_tilde;
            params.sigma_I_tilde = obj.sigma_I_tilde;
            params.E_W = obj.E_W;
            params.R = obj.R;

            % Computed E/I params
            params.n_E = obj.n_E;
            params.n_I = obj.n_I;
            params.E_indices = obj.E_indices;
            params.I_indices = obj.I_indices;
            params.N_sys_eqs = obj.N_sys_eqs;

            % Adaptation params
            params.n_a_E = obj.n_a_E;
            params.n_a_I = obj.n_a_I;
            params.tau_a_E = obj.tau_a_E;
            params.tau_a_I = obj.tau_a_I;
            params.c_E = obj.c_E;
            params.c_I = obj.c_I;

            % STD params
            params.n_b_E = obj.n_b_E;
            params.n_b_I = obj.n_b_I;
            params.tau_b_E_rec = obj.tau_b_E_rec;
            params.tau_b_E_rel = obj.tau_b_E_rel;
            params.tau_b_I_rec = obj.tau_b_I_rec;
            params.tau_b_I_rel = obj.tau_b_I_rel;
            params.std_zero_floor = obj.std_zero_floor;

            % Dynamics
            params.tau_d = obj.tau_d;
            params.activation_function = obj.activation_function;
            params.activation_function_derivative = obj.activation_function_derivative;
            params.x0_std = obj.x0_std;

            % Connection matrix (if built)
            if ~isempty(obj.W)
                params.W = obj.W;
            end

            % RNG seeds
            params.rng_seeds = obj.rng_seeds;
        end

    end

    %% Protected Build Sub-Methods (overridable by subclasses)
    methods (Access = protected)
        function build_network(obj)
            % BUILD_NETWORK Create the connectivity matrix W
            %
            % Sets RNG seed, fills in default tilde/tau parameters if empty,
            % then creates and scales the W matrix via RMTMatrix.
            % Subclasses can override for custom connectivity.

            % Set RNG seed for network generation
            rng(obj.rng_seeds(1));

            % Compute RMT tilde defaults if not set
            % F = 1/sqrt(N*alpha*(2-alpha)), see parameter_table.md
            F = obj.default_val;

            if isempty(obj.mu_E_tilde),    obj.mu_E_tilde = 3*F;     end % 3*F
            if isempty(obj.mu_I_tilde),    obj.mu_I_tilde = -4*F;    end % -4*F
            if isempty(obj.sigma_E_tilde), obj.sigma_E_tilde = 1*F;    end
            if isempty(obj.sigma_I_tilde), obj.sigma_I_tilde = 1*F;    end

            % Compute tau_a arrays if n_a > 0 but tau_a not set
            if obj.n_a_E > 0 && isempty(obj.tau_a_E)
                obj.tau_a_E = logspace(log10(0.25), log10(10), obj.n_a_E);
            end
            if obj.n_a_I > 0 && isempty(obj.tau_a_I)
                obj.tau_a_I = logspace(log10(0.25), log10(10), obj.n_a_I);
            end

            % Create W matrix using RMTMatrix
            rmt = RMTMatrix(obj.n);
            rmt.alpha = obj.alpha;
            rmt.f = obj.f;
            rmt.mu_tilde_e = obj.mu_E_tilde + obj.E_W;
            rmt.mu_tilde_i = obj.mu_I_tilde + obj.E_W;
            rmt.sigma_tilde_e = obj.sigma_E_tilde;
            rmt.sigma_tilde_i = obj.sigma_I_tilde;
            rmt.zrs_mode = obj.zrs_mode;

            W_raw = rmt.W;

            % Scale W
            if obj.rescale_by_abscissa
                W_eigs = eig(W_raw);
                abscissa_0 = max(real(W_eigs));
                gamma = 1 / abscissa_0;
                obj.W = obj.level_of_chaos * gamma * W_raw;
            else
                obj.W = obj.level_of_chaos * W_raw;
            end

            % Report info
            W_eigs_scaled = eig(obj.W);
            fprintf('W matrix created: spectral radius = %.3f, abscissa = %.3f, theoretical R = %.3f\n', ...
                max(abs(W_eigs_scaled)), max(real(W_eigs_scaled)), obj.R);
        end

        function build_stimulus(obj)
            % BUILD_STIMULUS Generate external stimulus, interpolant, and initial state
            %
            % Default implementation for SRNNModel2: generates sparse step
            % stimulus using generate_external_input, builds a linear
            % griddedInterpolant, and initializes the state vector S0.
            % SRNN_ESN_reservoir overrides this to generate ESN-specific stimulus.

            % Generate external stimulus
            obj.generate_stimulus();

            % Build interpolant for external input (avoids persistent variables)
            obj.u_interpolant = griddedInterpolant(obj.t_ex, obj.u_ex', 'linear', 'none');

            % Initialize state vector
            params_init = obj.get_params();
            obj.S0 = obj.initialize_state(params_init);
        end

    end

    %% Private Methods
    methods (Access = protected)
        function validate(obj)
            % VALIDATE Check parameter consistency and constraints

            % Check n_E and n_I
            if obj.n_E < 1
                error('SRNNModel:InvalidParams', 'n_E must be >= 1. Current: %d (n=%d, f=%.2f)', obj.n_E, obj.n, obj.f);
            end

            if obj.n_I < 1
                warning('SRNNModel:NoInhibitory', 'No inhibitory neurons (n_I=%d). Network may be unstable.', obj.n_I);
            end

            % Check adaptation consistency
            if obj.n_a_E > 0 && isempty(obj.tau_a_E)
                error('SRNNModel:InvalidParams', 'tau_a_E must be set when n_a_E > 0');
            end
            if obj.n_a_I > 0 && isempty(obj.tau_a_I)
                error('SRNNModel:InvalidParams', 'tau_a_I must be set when n_a_I > 0');
            end

            % Check STD consistency: tau_b_*_rec must be a 1 x n_b_* vector,
            % tau_b_*_rel a shared scalar. Coerce recovery constants to row vectors.
            if obj.n_b_E > 0
                if numel(obj.tau_b_E_rec) ~= obj.n_b_E
                    error('SRNNModel:InvalidParams', ...
                        'tau_b_E_rec must be a 1 x n_b_E vector (n_b_E = %d, numel(tau_b_E_rec) = %d).', ...
                        obj.n_b_E, numel(obj.tau_b_E_rec));
                end
                if ~isscalar(obj.tau_b_E_rel)
                    error('SRNNModel:InvalidParams', 'tau_b_E_rel must be a scalar (shared release across timescales).');
                end
                obj.tau_b_E_rec = reshape(obj.tau_b_E_rec, 1, []);
            end
            if obj.n_b_I > 0
                if numel(obj.tau_b_I_rec) ~= obj.n_b_I
                    error('SRNNModel:InvalidParams', ...
                        'tau_b_I_rec must be a 1 x n_b_I vector (n_b_I = %d, numel(tau_b_I_rec) = %d).', ...
                        obj.n_b_I, numel(obj.tau_b_I_rec));
                end
                if ~isscalar(obj.tau_b_I_rel)
                    error('SRNNModel:InvalidParams', 'tau_b_I_rel must be a scalar (shared release across timescales).');
                end
                obj.tau_b_I_rec = reshape(obj.tau_b_I_rec, 1, []);
            end

            % Check T_range
            if obj.T_range(2) <= obj.T_range(1)
                error('SRNNModel:InvalidParams', 'T_range(2) must be > T_range(1)');
            end

            % Check level_of_chaos
            if obj.level_of_chaos <= 0
                warning('SRNNModel:InvalidParams', 'level_of_chaos should be > 0. Current: %.2f', obj.level_of_chaos);
            end
        end

        function generate_stimulus(obj)
            % GENERATE_STIMULUS Create external input u_ex using generate_external_input

            dt = 1 / obj.fs;
            T_stim = obj.T_range(2);

            % Set intrinsic drive if not specified
            if isempty(obj.input_config.intrinsic_drive)
                obj.input_config.intrinsic_drive = zeros(obj.n, 1);
            end

            % Create params struct for generate_external_input
            params_stim = struct('n', obj.n, 'f', obj.f, ...
                'E_indices', obj.E_indices, 'I_indices', obj.I_indices);

            % Generate stimulus
            [u_stim, t_stim] = SRNNModel2.generate_external_input(params_stim, T_stim, obj.fs, obj.rng_seeds(2), obj.input_config);

            % Handle negative start time (prepend zeros for settling)
            if obj.T_range(1) < 0
                t_pre = (obj.T_range(1):dt:-dt)';
                u_pre = zeros(obj.n, length(t_pre));
                obj.t_ex = [t_pre; t_stim];
                obj.u_ex = [u_pre, u_stim];
            else
                % Slice if start time is positive
                indices = t_stim >= obj.T_range(1);
                obj.t_ex = t_stim(indices);
                obj.u_ex = u_stim(:, indices);
            end

            % Apply scaling
            obj.u_ex = obj.u_ex .* obj.u_ex_scale;

            fprintf('External stimulus generated: %d time points, %d neurons\n', length(obj.t_ex), obj.n);
        end

        function decimate_and_unpack(obj)
            % DECIMATE_AND_UNPACK Decimate state data and unpack for plotting

            params = obj.cached_params;

            % Decimate
            [t_plot, S_plot, plot_indices] = obj.decimate_states(obj.t_out, obj.S_out, obj.plot_deci);

            % Unpack state vector and compute firing rates
            [x_plot, a_plot, b_plot, r_plot, br_plot] = obj.unpack_and_compute_states(S_plot, params);

            % Split external input into E and I
            u_ex_plot = obj.u_ex(:, plot_indices);
            u_plot.E = u_ex_plot(obj.E_indices, :);
            u_plot.I = u_ex_plot(obj.I_indices, :);

            % Trim to T_plot if specified
            T_plot_range = obj.T_plot;
            if isempty(T_plot_range)
                T_plot_range = obj.T_range;
            end

            keep_mask = t_plot >= T_plot_range(1) & t_plot <= T_plot_range(2);

            t_plot = t_plot(keep_mask);
            u_plot.E = u_plot.E(:, keep_mask);
            u_plot.I = u_plot.I(:, keep_mask);
            x_plot = obj.trim_struct_data(x_plot, 2, keep_mask);
            r_plot = obj.trim_struct_data(r_plot, 2, keep_mask);
            b_plot = obj.trim_struct_data(b_plot, 3, keep_mask);  % b.E/b.I are n x n_b x nt (time = dim 3)
            br_plot = obj.trim_struct_data(br_plot, 2, keep_mask);
            a_plot = obj.trim_struct_data(a_plot, 3, keep_mask);

            % Trim Lyapunov results if present
            if ~isempty(obj.lya_results) && isfield(obj.lya_results, 't_lya')
                keep_mask_lya = obj.lya_results.t_lya >= T_plot_range(1) & obj.lya_results.t_lya <= T_plot_range(2);
                obj.lya_results.t_lya = obj.lya_results.t_lya(keep_mask_lya);

                if isfield(obj.lya_results, 'local_lya')
                    obj.lya_results.local_lya = obj.lya_results.local_lya(keep_mask_lya);
                end
                if isfield(obj.lya_results, 'finite_lya')
                    obj.lya_results.finite_lya = obj.lya_results.finite_lya(keep_mask_lya);
                end
                if isfield(obj.lya_results, 'local_LE_spectrum_t')
                    obj.lya_results.local_LE_spectrum_t = obj.lya_results.local_LE_spectrum_t(keep_mask_lya, :);
                end
                if isfield(obj.lya_results, 'finite_LE_spectrum_t')
                    obj.lya_results.finite_LE_spectrum_t = obj.lya_results.finite_LE_spectrum_t(keep_mask_lya, :);
                end
            end

            % Store plot data
            obj.plot_data = struct();
            obj.plot_data.t = t_plot;
            obj.plot_data.u = u_plot;
            obj.plot_data.x = x_plot;
            obj.plot_data.r = r_plot;
            obj.plot_data.a = a_plot;
            obj.plot_data.b = b_plot;
            obj.plot_data.br = br_plot;
        end

        function S0 = initialize_state(~, params)
            % INITIALIZE_STATE Initialize state vector for SRNN
            %
            % Initializes the complete state vector for an SRNN with adaptation and
            % short-term depression. The state is organized as:
            %   S0 = [a_E(:); a_I(:); b_E(:); b_I(:); x(:)]
            %
            %   - a_E, a_I: Adaptation variables (initialized to 0)
            %   - b_E, b_I: STD variables (initialized to 1, no depression)
            %   - x: Dendritic states (initialized to small random values)

            % Initialize adaptation states for E neurons
            a0_E = [];
            if params.n_a_E > 0
                a0_E = zeros(params.n_E * params.n_a_E, 1);
            end

            % Initialize adaptation states for I neurons
            a0_I = [];
            if params.n_a_I > 0
                a0_I = zeros(params.n_I * params.n_a_I, 1);
            end

            % Initialize STD states for E neurons (b variables start at 1.0, no depression)
            b0_E = [];
            if params.n_b_E > 0
                b0_E = ones(params.n_E * params.n_b_E, 1);
            end

            % Initialize STD states for I neurons
            b0_I = [];
            if params.n_b_I > 0
                b0_I = ones(params.n_I * params.n_b_I, 1);
            end

            % Initialize dendritic states (small random values, or 0 if x0_std=0)
            x0_std = SRNNModel2.safe_get_param(params, 'x0_std', 0.1);
            x0 = x0_std .* randn(params.n, 1);

            % Pack state vector: [a_E; a_I; b_E; b_I; x]
            S0 = [a0_E; a0_I; b0_E; b0_I; x0];
        end

        function [x, a, b, r, br] = unpack_and_compute_states(~, S_out, params, a_zeros_b_ones)
            % UNPACK_AND_COMPUTE_STATES Unpack state vector and compute dependent variables
            %
            % Unpacks the state trajectory S_out into individual state variables,
            % splits them into excitatory and inhibitory components, and computes
            % the firing rate r and synaptic output br.
            %
            % Inputs:
            %   S_out          - State trajectory (nt x N_sys_eqs) or column vector
            %   params         - Struct containing network parameters
            %   a_zeros_b_ones - (Optional) If true, returns a as zeros and b as ones

            % Handle optional parameter
            if nargin < 4
                a_zeros_b_ones = false;
            end

            nt = size(S_out, 1);
            current_idx = 0;

            %% Unpack adaptation states

            % Unpack adaptation states for E neurons (a_E)
            len_a_E = params.n_E * params.n_a_E;
            if len_a_E > 0
                a_E_ts = reshape(S_out(:, current_idx + (1:len_a_E))', params.n_E, params.n_a_E, nt);
            else
                a_E_ts = [];
            end
            current_idx = current_idx + len_a_E;

            % Unpack adaptation states for I neurons (a_I)
            len_a_I = params.n_I * params.n_a_I;
            if len_a_I > 0
                a_I_ts = reshape(S_out(:, current_idx + (1:len_a_I))', params.n_I, params.n_a_I, nt);
            else
                a_I_ts = [];
            end
            current_idx = current_idx + len_a_I;

            %% Unpack STD variables (b)

            % Unpack b states for E neurons (b_E), n_E x n_b_E x nt (timescale-major)
            len_b_E = params.n_E * params.n_b_E;
            if len_b_E > 0
                b_E_ts = reshape(S_out(:, current_idx + (1:len_b_E))', params.n_E, params.n_b_E, nt);
            else
                b_E_ts = [];
            end
            current_idx = current_idx + len_b_E;

            % Unpack b states for I neurons (b_I), n_I x n_b_I x nt (timescale-major)
            len_b_I = params.n_I * params.n_b_I;
            if len_b_I > 0
                b_I_ts = reshape(S_out(:, current_idx + (1:len_b_I))', params.n_I, params.n_b_I, nt);
            else
                b_I_ts = [];
            end
            current_idx = current_idx + len_b_I;

            %% Unpack dendritic states (x)
            x_ts = S_out(:, current_idx + (1:params.n))';  % n x nt

            %% Compute firing rates with adaptation and STD

            % Compute effective dendritic state (x_eff = x - c * sum(a))
            x_eff_ts = x_ts;  % n x nt

            % Apply adaptation effect to E neurons (scaled by c_E)
            if params.n_E > 0 && params.n_a_E > 0 && ~isempty(a_E_ts)
                % sum(a_E_ts, 2) is n_E x 1 x nt, need to squeeze to n_E x nt
                sum_a_E = squeeze(sum(a_E_ts, 2));  % n_E x nt
                if size(sum_a_E, 1) ~= params.n_E  % Handle case where nt=1
                    sum_a_E = sum_a_E';
                end
                x_eff_ts(params.E_indices, :) = x_eff_ts(params.E_indices, :) - params.c_E * sum_a_E;
            end

            % Apply adaptation effect to I neurons (scaled by c_I)
            if params.n_I > 0 && params.n_a_I > 0 && ~isempty(a_I_ts)
                sum_a_I = squeeze(sum(a_I_ts, 2));  % n_I x nt
                if size(sum_a_I, 1) ~= params.n_I
                    sum_a_I = sum_a_I';
                end
                x_eff_ts(params.I_indices, :) = x_eff_ts(params.I_indices, :) - params.c_I * sum_a_I;
            end

            % Apply STD effect: the synaptic depression factor is the product
            % over timescales, prod_m b_{i,m} (collapse the timescale dim like
            % SFA sums it above, but multiplicatively).
            b_ts = ones(params.n, nt);  % Initialize b = 1 for all neurons (no depression)
            if params.n_b_E > 0 && ~isempty(b_E_ts)
                prod_b_E = squeeze(prod(b_E_ts, 2));  % n_E x nt
                if size(prod_b_E, 1) ~= params.n_E   % Handle case where nt=1
                    prod_b_E = prod_b_E';
                end
                b_ts(params.E_indices, :) = prod_b_E;
            end
            if params.n_b_I > 0 && ~isempty(b_I_ts)
                prod_b_I = squeeze(prod(b_I_ts, 2));  % n_I x nt
                if size(prod_b_I, 1) ~= params.n_I
                    prod_b_I = prod_b_I';
                end
                b_ts(params.I_indices, :) = prod_b_I;
            end

            % Compute firing rates: r = phi(x_eff) (raw rate)
            r_ts = params.activation_function(x_eff_ts);  % n x nt

            % Compute synaptic output: br = b .* r (presynaptically depressed).
            % Apply the optional STD zero-floor here so br matches the drive
            % used in dynamics_fast; b_ts (the plotted state) stays raw.
            b_syn_ts = SRNNModel2.apply_std_zero_floor(b_ts, params);
            br_ts = b_syn_ts .* r_ts; % n x nt

            %% Split into E and I components and package into structs

            % x: dendritic states
            x.E = x_ts(params.E_indices, :);  % n_E x nt
            x.I = x_ts(params.I_indices, :);  % n_I x nt

            % a: adaptation variables
            a.E = a_E_ts;  % n_E x n_a_E x nt (or empty)
            a.I = a_I_ts;  % n_I x n_a_I x nt (or empty)

            % b: STD variables (full per-timescale array, mirroring a)
            if isempty(b_E_ts)
                b.E = ones(params.n_E, max(1, params.n_b_E), nt);  % Default to no depression
            else
                b.E = b_E_ts;  % n_E x n_b_E x nt
            end

            if isempty(b_I_ts)
                b.I = ones(params.n_I, max(1, params.n_b_I), nt);  % Default to no depression
            else
                b.I = b_I_ts;  % n_I x n_b_I x nt
            end

            % r: firing rates
            r.E = r_ts(params.E_indices, :);  % n_E x nt
            r.I = r_ts(params.I_indices, :);  % n_I x nt

            % br: synaptic output
            br.E = br_ts(params.E_indices, :);
            br.I = br_ts(params.I_indices, :);

            %% Override with zeros/ones if requested (for Jacobian computation)
            if a_zeros_b_ones
                % Return x as simple array instead of struct (n x nt)
                x = x_ts;

                % Replace a with zeros (n x nt)
                a = zeros(params.n, nt);

                % Replace b with ones (n x nt)
                b = ones(params.n, nt);

                % Return r as simple array instead of struct (n x nt)
                r = r_ts;

                % Return br as simple array (equal to r since b=1)
                br = r_ts;
            end
        end
    end

    methods (Static, Access = protected)
        function b_used = apply_std_zero_floor(b, params)
            % APPLY_STD_ZERO_FLOOR Rescale synaptic availability so full
            % depression drives the synaptic output b*r -> 0.
            %
            % Operates on the collapsed per-neuron product P = prod_m b_{i,m}.
            % Maps the dynamic range [P_min, 1] onto [0, 1], where
            % P_min = prod_m tau_rel/(tau_rec_m + tau_rel) is the product of the
            % per-timescale b-ODE asymptotes at r = 1 (a per-population scalar).
            % This is a readout-only transform: the b state and its ODE are
            % unchanged. At rest (all b = 1) P = 1 and b_used = 1, so baseline
            % behavior is preserved. Accepts b as an n x 1 vector or n x nt
            % matrix (P_min broadcasts over columns).
            b_used = b;
            if ~isfield(params, 'std_zero_floor') || ~params.std_zero_floor
                return;
            end
            if params.n_b_E > 0
                Pmin_E = prod(params.tau_b_E_rel ./ (params.tau_b_E_rec + params.tau_b_E_rel));
                b_used(params.E_indices, :) = (b(params.E_indices, :) - Pmin_E) / (1 - Pmin_E);
            end
            if params.n_b_I > 0
                Pmin_I = prod(params.tau_b_I_rel ./ (params.tau_b_I_rec + params.tau_b_I_rel));
                b_used(params.I_indices, :) = (b(params.I_indices, :) - Pmin_I) / (1 - Pmin_I);
            end
        end

        function dS_dt = dynamics_fast(t, S, params)
            % DYNAMICS_FAST Static method for fast ODE evaluation
            %
            % This static method avoids OOP overhead by taking all parameters
            % as a struct. The interpolant must be in params.u_interpolant.
            %
            % Implements:
            %   dx_i/dt = (-x_i + sum_j(w_ij * r_j) + u_i) / tau_d
            %   r_i = b_i * phi(x_i - c * sum_k(a_i,k))
            %   da_i,k/dt = (-a_i,k + r_i) / tau_k
            %   db_i/dt = (1 - b_i) / tau_rec - (b_i * r_i) / tau_rel
            %
            % State organization: S = [a_E(:); a_I(:); b_E(:); b_I(:); x(:)]

            %% Interpolate external input
            u = params.u_interpolant(t)';  % n x 1

            %% Load parameters
            n = params.n;
            n_E = params.n_E;
            n_I = params.n_I;
            E_indices = params.E_indices;
            I_indices = params.I_indices;

            n_a_E = params.n_a_E;
            n_a_I = params.n_a_I;
            n_b_E = params.n_b_E;
            n_b_I = params.n_b_I;

            W = params.W;
            tau_d = params.tau_d;
            tau_a_E = params.tau_a_E;
            tau_a_I = params.tau_a_I;
            tau_b_E_rec = params.tau_b_E_rec;
            tau_b_E_rel = params.tau_b_E_rel;
            tau_b_I_rec = params.tau_b_I_rec;
            tau_b_I_rel = params.tau_b_I_rel;

            c_E = params.c_E;
            c_I = params.c_I;
            activation_fn = params.activation_function;

            %% Unpack state variables
            current_idx = 0;

            % Adaptation states for E neurons (a_E)
            len_a_E = n_E * n_a_E;
            if len_a_E > 0
                a_E = reshape(S(current_idx + (1:len_a_E)), n_E, n_a_E);
            else
                a_E = [];
            end
            current_idx = current_idx + len_a_E;

            % Adaptation states for I neurons (a_I)
            len_a_I = n_I * n_a_I;
            if len_a_I > 0
                a_I = reshape(S(current_idx + (1:len_a_I)), n_I, n_a_I);
            else
                a_I = [];
            end
            current_idx = current_idx + len_a_I;

            % STD states for E neurons (b_E), reshaped to n_E x n_b_E (timescale-major)
            len_b_E = n_E * n_b_E;
            if len_b_E > 0
                b_E = reshape(S(current_idx + (1:len_b_E)), n_E, n_b_E);
            else
                b_E = [];
            end
            current_idx = current_idx + len_b_E;

            % STD states for I neurons (b_I), reshaped to n_I x n_b_I (timescale-major)
            len_b_I = n_I * n_b_I;
            if len_b_I > 0
                b_I = reshape(S(current_idx + (1:len_b_I)), n_I, n_b_I);
            else
                b_I = [];
            end
            current_idx = current_idx + len_b_I;

            % Dendritic states (x)
            x = S(current_idx + (1:n));

            %% Compute firing rates
            x_eff = x;

            % Apply adaptation effect to E neurons
            if n_E > 0 && n_a_E > 0 && ~isempty(a_E)
                x_eff(E_indices) = x_eff(E_indices) - c_E * sum(a_E, 2);
            end

            % Apply adaptation effect to I neurons
            if n_I > 0 && n_a_I > 0 && ~isempty(a_I)
                x_eff(I_indices) = x_eff(I_indices) - c_I * sum(a_I, 2);
            end

            % Apply STD effect: the synaptic depression factor is the product
            % of the per-timescale resources, prod_m b_{i,m} (n_E x 1 / n_I x 1).
            b = ones(n, 1);
            if n_b_E > 0 && ~isempty(b_E)
                b(E_indices) = prod(b_E, 2);
            end
            if n_b_I > 0 && ~isempty(b_I)
                b(I_indices) = prod(b_I, 2);
            end

            r = activation_fn(x_eff);

            % Optional STD zero-floor: rescale b in the synaptic readout only
            % (the b ODE below still uses the raw b state).
            b_syn = SRNNModel2.apply_std_zero_floor(b, params);

            %% Compute derivatives
            dx_dt = (-x + W * (b_syn .* r) + u) / tau_d;

            da_E_dt = [];
            if n_E > 0 && n_a_E > 0 && ~isempty(a_E)
                da_E_dt = (r(E_indices) - a_E) ./ tau_a_E;
            end

            da_I_dt = [];
            if n_I > 0 && n_a_I > 0 && ~isempty(a_I)
                da_I_dt = (r(I_indices) - a_I) ./ tau_a_I;
            end

            % db_{i,m}/dt: per-timescale recovery (1 x n_b vector) with a shared
            % release scalar. b_* is n x n_b, tau_b_*_rec is 1 x n_b, r(idx) is
            % n x 1 and broadcasts across timescale columns (as with da/dt).
            db_E_dt = [];
            if n_E > 0 && n_b_E > 0 && ~isempty(b_E)
                db_E_dt = (1 - b_E) ./ tau_b_E_rec - (b_E .* r(E_indices)) / tau_b_E_rel;
            end

            db_I_dt = [];
            if n_I > 0 && n_b_I > 0 && ~isempty(b_I)
                db_I_dt = (1 - b_I) ./ tau_b_I_rec - (b_I .* r(I_indices)) / tau_b_I_rel;
            end

            %% Pack derivatives
            dS_dt = [da_E_dt(:); da_I_dt(:); db_E_dt(:); db_I_dt(:); dx_dt];
        end

        %% ====================================================================
        %                    INTERNAL LYAPUNOV FUNCTIONS
        % =====================================================================
        % Internalized from ConnectivityAdaptation to avoid path conflicts.

    end

    %% ====================================================================
    %              INTERNALIZED ACTIVATION FUNCTIONS
    % =====================================================================
    % Internalized from src/nonlinearities/ to make SRNNModel2 standalone.


    %% ====================================================================
    %              INTERNALIZED STIMULUS GENERATION
    % =====================================================================
    % Internalized from src/generate_stimulus/ to make SRNNModel2 standalone.

    methods (Static, Access = protected)
        function [u_ex, t_ex] = generate_external_input(params, T, fs, rng_seed, input_config)
            % GENERATE_EXTERNAL_INPUT Generate external input for SRNN simulation.
            % Internalized from src/generate_stimulus/generate_external_input.m

            % Check for custom generator function
            if isfield(input_config, 'generator') && isa(input_config.generator, 'function_handle')
                [u_ex, t_ex] = input_config.generator(params, T, fs, rng_seed, input_config);
                return;
            end

            % Standard sparse step stimulus generation
            rng(rng_seed);

            dt = 1 / fs;
            t_ex = (0:dt:T)';
            nt = length(t_ex);

            n_steps = input_config.n_steps;
            step_density = input_config.step_density;
            amp = input_config.amp;
            no_stim_pattern = input_config.no_stim_pattern;
            intrinsic_drive = input_config.intrinsic_drive;

            step_period = fix(T / n_steps);
            step_length = round(step_period * fs);

            if isfield(input_config, 'positive_only') && input_config.positive_only
                random_sparse_step = amp * abs(randn(params.n, n_steps));
            else
                random_sparse_step = amp * randn(params.n, n_steps);
            end

            sparse_mask = false(params.n, n_steps);

            if isfield(input_config, 'step_density_E')
                density_E = input_config.step_density_E;
            else
                density_E = step_density;
            end

            if isfield(input_config, 'step_density_I')
                density_I = input_config.step_density_I;
            else
                density_I = step_density;
            end

            if isfield(params, 'E_indices') && isfield(params, 'I_indices')
                E_idx = params.E_indices;
                I_idx = params.I_indices;
            else
                f_val = 0.5;
                if isfield(params, 'f')
                    f_val = params.f;
                end
                n_E = round(params.n * f_val);
                E_idx = 1:n_E;
                I_idx = (n_E + 1):params.n;
            end

            sparse_mask(E_idx, :) = rand(length(E_idx), n_steps) < density_E;
            sparse_mask(I_idx, :) = rand(length(I_idx), n_steps) < density_I;

            random_sparse_step = random_sparse_step .* sparse_mask;
            random_sparse_step(:, no_stim_pattern) = 0;

            u_ex = zeros(params.n, nt);
            for step_idx = 1:n_steps
                start_idx = (step_idx - 1) * step_length + 1;
                end_idx = min(step_idx * step_length, nt);
                if start_idx > nt
                    break;
                end
                u_ex(:, start_idx:end_idx) = repmat(random_sparse_step(:, step_idx), 1, end_idx - start_idx + 1);
            end

            u_ex = u_ex + intrinsic_drive;
        end
    end

    %% ====================================================================
    %              INTERNALIZED JACOBIAN COMPUTATION
    % =====================================================================
    % Internalized from src/algorithms/Jacobian/ to make SRNNModel2 standalone.

    methods (Static)
        function J = compute_Jacobian_fast(S, params)
            % COMPUTE_JACOBIAN_FAST Sparse/vectorized Jacobian assembly for the SRNN system.
            % Internalized from src/algorithms/Jacobian/compute_Jacobian_fast.m

            n = params.n;
            n_E = params.n_E;
            n_I = params.n_I;
            E_indices = params.E_indices;
            I_indices = params.I_indices;

            n_a_E = params.n_a_E;
            n_a_I = params.n_a_I;
            n_b_E = params.n_b_E;
            n_b_I = params.n_b_I;

            W = params.W;
            tau_d = params.tau_d;
            tau_a_E = params.tau_a_E;
            tau_a_I = params.tau_a_I;
            tau_b_E_rec = params.tau_b_E_rec;
            tau_b_E_rel = params.tau_b_E_rel;
            tau_b_I_rec = params.tau_b_I_rec;
            tau_b_I_rel = params.tau_b_I_rel;

            c_E = SRNNModel2.safe_get_param(params, 'c_E', 1.0);
            c_I = SRNNModel2.safe_get_param(params, 'c_I', 1.0);

            if ~isfield(params, 'activation_function_derivative') || ...
               ~isa(params.activation_function_derivative, 'function_handle')
                error('compute_Jacobian_fast:MissingActivationFunctionDerivative', ...
                      'params.activation_function_derivative must be provided as a function handle');
            end
            phi_prime = params.activation_function_derivative;

            if ~isfield(params, 'activation_function') || ...
               ~isa(params.activation_function, 'function_handle')
                error('compute_Jacobian_fast:MissingActivationFunction', ...
                      'params.activation_function must be provided as a function handle');
            end
            phi_fun = params.activation_function;

            %% Unpack state variables
            current_idx = 0;
            len_a_E = n_E * n_a_E;
            len_a_I = n_I * n_a_I;
            len_b_E = n_E * n_b_E;
            len_b_I = n_I * n_b_I;

            if len_a_E > 0
                a_E = reshape(S(current_idx + (1:len_a_E)), n_E, n_a_E);
            else
                a_E = [];
            end
            current_idx = current_idx + len_a_E;

            if len_a_I > 0
                a_I = reshape(S(current_idx + (1:len_a_I)), n_I, n_a_I);
            else
                a_I = [];
            end
            current_idx = current_idx + len_a_I;

            if len_b_E > 0
                b_E = reshape(S(current_idx + (1:len_b_E)), n_E, n_b_E);
            else
                b_E = [];
            end
            current_idx = current_idx + len_b_E;

            if len_b_I > 0
                b_I = reshape(S(current_idx + (1:len_b_I)), n_I, n_b_I);
            else
                b_I = [];
            end
            current_idx = current_idx + len_b_I;

            x = S(current_idx + (1:n));

            %% Effective potentials and rates
            x_eff = x;
            if len_a_E > 0
                x_eff(E_indices) = x_eff(E_indices) - c_E * sum(a_E, 2);
            end
            if len_a_I > 0
                x_eff(I_indices) = x_eff(I_indices) - c_I * sum(a_I, 2);
            end

            % Collapsed per-neuron depression factor P = prod_m b_{i,m}.
            b = ones(n, 1);
            P_E = [];
            P_I = [];
            if len_b_E > 0
                P_E = prod(b_E, 2);        % n_E x 1
                b(E_indices) = P_E;
            end
            if len_b_I > 0
                P_I = prod(b_I, 2);        % n_I x 1
                b(I_indices) = P_I;
            end

            phi_x_eff = phi_fun(x_eff);
            phi_prime_x_eff = phi_prime(x_eff);

            % STD zero-floor affects only the synaptic drive W*(b_syn.*r) in
            % dx/dt, so only the dx/dt (row_x) Jacobian blocks below use b_syn /
            % the per-population gain g = d(P_syn)/dP = 1/(1-P_min), where
            % P_min = prod_m tau_rel/(tau_rec_m+tau_rel) (a scalar). The b ODE
            % and SFA blocks use the raw b state and are unchanged.
            b_syn = SRNNModel2.apply_std_zero_floor(b, params);
            g_b_E = 1;
            g_b_I = 1;
            if isfield(params, 'std_zero_floor') && params.std_zero_floor
                if len_b_E > 0
                    Pmin_E = prod(tau_b_E_rel ./ (tau_b_E_rec + tau_b_E_rel));
                    g_b_E = 1 / (1 - Pmin_E);
                end
                if len_b_I > 0
                    Pmin_I = prod(tau_b_I_rel ./ (tau_b_I_rec + tau_b_I_rel));
                    g_b_I = 1 / (1 - Pmin_I);
                end
            end

            %% Dimensions and indexing
            N_sys_eqs = len_a_E + len_a_I + len_b_E + len_b_I + n;

            row_a_E = 1:len_a_E;
            row_a_I = len_a_E + (1:len_a_I);
            row_b_E = len_a_E + len_a_I + (1:len_b_E);
            row_b_I = len_a_E + len_a_I + len_b_E + (1:len_b_I);
            row_x   = len_a_E + len_a_I + len_b_E + len_b_I + (1:n);

            col_a_E = row_a_E;
            col_a_I = row_a_I;
            col_b_E = row_b_E;
            col_b_I = row_b_I;
            col_x   = row_x;

            W_sparse = sparse(W);
            J = sparse(N_sys_eqs, N_sys_eqs);

            %% SFA blocks (E) — Jacobian of da_{ik}/dt = (r_i - a_{ik})/tau_k,
            % r_i = phi(x_eff_i), x_eff_i = x_i - c_E*sum_k a_{ik}. Adaptation is
            % driven by the raw rate r (no b). Timescale-major state ordering:
            % a_{i,k} lives at index (k-1)*n_E + i.
            if len_a_E > 0
                tau_inv_E = 1 ./ tau_a_E(:);
                gamma_E = c_E * phi_prime_x_eff(E_indices);   % n_E x 1, no b
                % a -> a: diagonal -1/tau_k plus within-neuron cross-timescale coupling
                diag_block_E = kron(spdiags(-tau_inv_E, 0, n_a_E, n_a_E), speye(n_E));
                row_template_E = sparse(tau_inv_E * ones(1, n_a_E));   % (k,k') -> tau_inv_k
                coupling_block_E = kron(row_template_E, spdiags(-gamma_E, 0, n_E, n_E));
                J(row_a_E, col_a_E) = diag_block_E + coupling_block_E;

                % a -> x: (1/tau_k)*phi'_i
                beta_E = phi_prime_x_eff(E_indices);          % n_E x 1, no b
                vals = kron(tau_inv_E, beta_E);               % (k-1)*n_E+i -> tau_inv_k*phi'_i
                rows = (1:len_a_E)';
                cols = repmat(E_indices(:), n_a_E, 1);        % (k-1)*n_E+i -> E_indices(i)
                J(row_a_E, col_x) = sparse(rows, cols, vals, len_a_E, n);
                % a -> b block is zero: da/dt does not depend on b.
            end

            %% SFA blocks (I) — same structure as E (raw-rate drive, ts-major)
            if len_a_I > 0
                tau_inv_I = 1 ./ tau_a_I(:);
                gamma_I = c_I * phi_prime_x_eff(I_indices);   % n_I x 1, no b
                diag_block_I = kron(spdiags(-tau_inv_I, 0, n_a_I, n_a_I), speye(n_I));
                row_template_I = sparse(tau_inv_I * ones(1, n_a_I));
                coupling_block_I = kron(row_template_I, spdiags(-gamma_I, 0, n_I, n_I));
                J(row_a_I, col_a_I) = diag_block_I + coupling_block_I;

                beta_I = phi_prime_x_eff(I_indices);          % n_I x 1, no b
                vals = kron(tau_inv_I, beta_I);
                rows = (1:len_a_I)';
                cols = repmat(I_indices(:), n_a_I, 1);
                J(row_a_I, col_x) = sparse(rows, cols, vals, len_a_I, n);
                % a -> b block is zero.
            end

            %% STD blocks (E) — Jacobian of db_{i,m}/dt = (1-b_{i,m})/tau_rec_m
            % - b_{i,m}*r_i/tau_rel, r_i = phi(x_eff_i) (raw rate). Each timescale
            % is an independent ODE; only the shared rate r_i couples them via a
            % and x. Timescale-major: b_{i,m} lives at index (m-1)*n_E + i.
            if len_b_E > 0
                phi_prime_E = phi_prime_x_eff(E_indices);   % n_E x 1
                % b -> a: -b_{i,m}/tau_rel * dr_i/da_{ik} = b_{i,m}*c_E*phi'_i/tau_rel.
                % Same for every SFA timescale k, so replicate across the n_a_E col-blocks.
                if len_a_E > 0
                    coeff_a_E = b_E .* (c_E .* phi_prime_E / tau_b_E_rel);   % n_E x n_b_E
                    stack_a_E = sparse(1:len_b_E, repmat((1:n_E)', n_b_E, 1), coeff_a_E(:), len_b_E, n_E);
                    J(row_b_E, col_a_E) = kron(sparse(ones(1, n_a_E)), stack_a_E);
                end
                % b -> b: diagonal, -1/tau_rec_m - r_i/tau_rel
                diag_vals_b_E = kron(-1 ./ tau_b_E_rec(:), ones(n_E, 1)) ...
                              + repmat(-phi_x_eff(E_indices) / tau_b_E_rel, n_b_E, 1);
                J(row_b_E, col_b_E) = spdiags(diag_vals_b_E, 0, len_b_E, len_b_E);
                % b -> x: -b_{i,m}/tau_rel * phi'_i
                vals_bx_E = -(b_E(:) .* repmat(phi_prime_E, n_b_E, 1)) / tau_b_E_rel;
                J(row_b_E, col_x) = sparse(1:len_b_E, repmat(E_indices(:), n_b_E, 1), vals_bx_E, len_b_E, n);
            end

            %% STD blocks (I) — same structure as E (raw-rate drive, ts-major)
            if len_b_I > 0
                phi_prime_I = phi_prime_x_eff(I_indices);
                if len_a_I > 0
                    coeff_a_I = b_I .* (c_I .* phi_prime_I / tau_b_I_rel);
                    stack_a_I = sparse(1:len_b_I, repmat((1:n_I)', n_b_I, 1), coeff_a_I(:), len_b_I, n_I);
                    J(row_b_I, col_a_I) = kron(sparse(ones(1, n_a_I)), stack_a_I);
                end
                diag_vals_b_I = kron(-1 ./ tau_b_I_rec(:), ones(n_I, 1)) ...
                              + repmat(-phi_x_eff(I_indices) / tau_b_I_rel, n_b_I, 1);
                J(row_b_I, col_b_I) = spdiags(diag_vals_b_I, 0, len_b_I, len_b_I);
                vals_bx_I = -(b_I(:) .* repmat(phi_prime_I, n_b_I, 1)) / tau_b_I_rel;
                J(row_b_I, col_x) = sparse(1:len_b_I, repmat(I_indices(:), n_b_I, 1), vals_bx_I, len_b_I, n);
            end

            %% dx/dt blocks
            if len_a_E > 0
                replicate_a_E = kron(ones(1, n_a_E), speye(n_E));   % ts-major cols
                block = -c_E * W_sparse(:, E_indices) * spdiags(b_syn(E_indices) .* phi_prime_x_eff(E_indices), 0, n_E, n_E);
                J(row_x, col_a_E) = (block * replicate_a_E) / tau_d;
            end

            if len_a_I > 0
                replicate_a_I = kron(ones(1, n_a_I), speye(n_I));   % ts-major cols
                block = -c_I * W_sparse(:, I_indices) * spdiags(b_syn(I_indices) .* phi_prime_x_eff(I_indices), 0, n_I, n_I);
                J(row_x, col_a_I) = (block * replicate_a_I) / tau_d;
            end

            % x -> b: d(dx_j/dt)/db_{i,m} = (1/tau_d) W_ji r_i g dP_i/db_{i,m},
            % with P_i = prod_m b_{i,m} so dP_i/db_{i,m} = P_i/b_{i,m}. Column
            % (m,i) scales W(:,E_i) by r_i*g*(P_i/b_{i,m}) (ts-major cols).
            if len_b_E > 0
                coeff_xb_E = (g_b_E * (phi_x_eff(E_indices) .* P_E)) ./ b_E;   % n_E x n_b_E
                D_E = sparse(repmat((1:n_E)', n_b_E, 1), 1:len_b_E, coeff_xb_E(:), n_E, len_b_E);
                J(row_x, col_b_E) = (W_sparse(:, E_indices) * D_E) / tau_d;
            end

            if len_b_I > 0
                coeff_xb_I = (g_b_I * (phi_x_eff(I_indices) .* P_I)) ./ b_I;   % n_I x n_b_I
                D_I = sparse(repmat((1:n_I)', n_b_I, 1), 1:len_b_I, coeff_xb_I(:), n_I, len_b_I);
                J(row_x, col_b_I) = (W_sparse(:, I_indices) * D_I) / tau_d;
            end

            diag_term = spdiags(-ones(n,1)/tau_d, 0, n, n);
            gain_diag = spdiags(b_syn .* phi_prime_x_eff, 0, n, n);
            J(row_x, col_x) = diag_term + (W_sparse * gain_diag) / tau_d;
        end

    end

    methods (Static, Access = private)
        function value = safe_get_param(params, field, default_value)
            % SAFE_GET_PARAM Helper to get a field from params with a default.
            if isfield(params, field)
                value = params.(field);
            else
                value = default_value;
            end
        end
    end

    %% ====================================================================
    %              INTERNALIZED PLOTTING FUNCTIONS
    % =====================================================================
    % Internalized from src/plotting/ to make SRNNModel2 standalone.

    methods (Static, Access = protected)
        function [fig_handle, ax_handles] = plot_SRNN_tseries(t_out, u, x, r, a, b, br, params, lya_results, Lya_method, T_plot)
            % PLOT_SRNN_TSERIES Create comprehensive time series plots for SRNN simulation.
            % Internalized from src/plotting/plot_SRNN_tseries.m

            if nargin < 11
                T_plot = [];
            end

            % Determine which subplots are needed
            has_adaptation = params.n_a_E > 0 || params.n_a_I > 0;
            has_std_vars = (params.n_b_E > 0 || params.n_b_I > 0) && ~isempty(b);
            has_synaptic_output = ~isempty(br);
            has_lyapunov = ~strcmpi(Lya_method, 'none');
            has_firing_rate = ~isempty(r);

            n_plots = 2;  % Always: External input, Dendritic states
            if has_firing_rate, n_plots = n_plots + 1; end
            if has_synaptic_output, n_plots = n_plots + 1; end
            if has_adaptation, n_plots = n_plots + 1; end
            if has_std_vars, n_plots = n_plots + 1; end
            if has_lyapunov, n_plots = n_plots + 1; end

            fig_handle = figure();
            tiledlayout(n_plots, 1);
            ax_handles = [];

            % Always: External input
            ax_handles(end+1) = nexttile;
            SRNNModel2.plot_external_input(t_out, u);
            set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');

            % Always: Dendritic states
            ax_handles(end+1) = nexttile;
            plot_mean = false;
            if isfield(params, 'plot_mean_dendrite')
                plot_mean = params.plot_mean_dendrite;
            end
            SRNNModel2.plot_dendritic_state(t_out, x, plot_mean);
            set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');

            % Conditionally: Firing rates
            if has_firing_rate
                ax_handles(end+1) = nexttile;
                SRNNModel2.plot_firing_rate(t_out, r);
                set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');
            end

            % Conditionally: Synaptic output (br)
            if has_synaptic_output
                ax_handles(end+1) = nexttile;
                SRNNModel2.plot_synaptic_output(t_out, br);
                set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');
            end

            % Conditionally: Adaptation variables
            if has_adaptation
                ax_handles(end+1) = nexttile;
                SRNNModel2.plot_adaptation(t_out, a, params);
                set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');
            end

            % Conditionally: STD variables (b)
            if has_std_vars
                ax_handles(end+1) = nexttile;
                SRNNModel2.plot_std_variable(t_out, b, params);
                set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');
            end

            % Conditionally: Lyapunov exponent(s)
            if has_lyapunov
                ax_handles(end+1) = nexttile;
                if strcmpi(Lya_method, 'benettin')
                    SRNNModel2.plot_lyapunov(lya_results, Lya_method, {'local', 'EOC', 'value'});
                else
                    SRNNModel2.plot_lyapunov(lya_results, Lya_method);
                end
                set(gca, 'XTick', [], 'XTickLabel', [], 'XColor', 'white');
            end

            linkaxes(ax_handles,'x');

            if ~isempty(T_plot)
                xlim(ax_handles(end), T_plot);
            end

            % Add time scale bar overlay in lower right of last subplot
            axes(ax_handles(end));
            hold on;
            xlims = xlim;
            ylims = ylim;

            scale_bar_length = round(0.1 * (xlims(2) - xlims(1)));
            if scale_bar_length < 1
                scale_bar_length = 0.1 * (xlims(2) - xlims(1));
            end

            x_end = xlims(1) + 0.95 * (xlims(2) - xlims(1));
            x_start = x_end - scale_bar_length;
            y_pos = ylims(1) + 0.10 * (ylims(2) - ylims(1));
            plot([x_start, x_end], [y_pos, y_pos], 'k-', 'LineWidth', 4);

            text_x = (x_start + x_end) / 2;
            text_y = ylims(1) + 0.05 * (ylims(2) - ylims(1));
            text(text_x, text_y, sprintf('%g seconds', scale_bar_length), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
            hold off;
        end

        function plot_external_input(t, u)
            % PLOT_EXTERNAL_INPUT Plot external input for E and I neurons.
            % Internalized from src/plotting/plot_external_input.m
            cmap_I = SRNNModel2.inhibitory_colormap(8);
            cmap_E = SRNNModel2.excitatory_colormap(8);
            SRNNModel2.plot_lines_with_colormap(t, u.I, cmap_I);
            hold on;
            SRNNModel2.plot_lines_with_colormap(t, u.E, cmap_E);
            hold off;
            ylabel('stim');
            yl = ylim;
            ylim(yl+[-0.05 0])
            yticks([-1 0 1]);
        end

        function plot_dendritic_state(t, x, plot_mean)
            % PLOT_DENDRITIC_STATE Plot dendritic states for E and I neurons.
            % Internalized from src/plotting/plot_dendritic_state.m
            if nargin < 3, plot_mean = false; end
            cmap_I = SRNNModel2.inhibitory_colormap(8);
            cmap_E = SRNNModel2.excitatory_colormap(8);
            SRNNModel2.plot_lines_with_colormap(t, x.I, cmap_I);
            hold on;
            SRNNModel2.plot_lines_with_colormap(t, x.E, cmap_E);
            if plot_mean
                mean_x = mean(x.E, 1);
                plot(t, mean_x, 'k', 'LineWidth', 3);
            end
            hold off;
            ylabel('dendrite');
        end

        function plot_firing_rate(t, r)
            % PLOT_FIRING_RATE Plot firing rates for E and I neurons.
            % Internalized from src/plotting/plot_firing_rate.m
            cmap_I = SRNNModel2.inhibitory_colormap(8);
            cmap_E = SRNNModel2.excitatory_colormap(8);
            SRNNModel2.plot_lines_with_colormap(t, r.I, cmap_I);
            hold on;
            SRNNModel2.plot_lines_with_colormap(t, r.E, cmap_E);
            hold off;
            ylabel('firing rate');
            yticks([0, 1]);
            ylim([0, 1]);
        end

        function plot_synaptic_output(t, br)
            % PLOT_SYNAPTIC_OUTPUT Plot synaptic output (br = b .* r) for E and I neurons.
            % Internalized from src/plotting/plot_synaptic_output.m
            cmap_I = SRNNModel2.inhibitory_colormap(8);
            cmap_E = SRNNModel2.excitatory_colormap(8);
            SRNNModel2.plot_lines_with_colormap(t, br.I, cmap_I);
            hold on;
            SRNNModel2.plot_lines_with_colormap(t, br.E, cmap_E);
            hold off;
            ylabel('synaptic output');
            yticks([0, 1]);
            ylim([0, 1]);
        end

        function plot_adaptation(t, a, params)
            % PLOT_ADAPTATION Plot adaptation variables for E and I neurons.
            % Internalized from src/plotting/plot_adaptation.m
            cmap_I = SRNNModel2.inhibitory_colormap(8);
            cmap_E = SRNNModel2.excitatory_colormap(8);
            has_adaptation = false;

            if ~isempty(a.I) && params.n_a_I > 0
                a_I_sum = sum(a.I, 2);
                a_I_summed = reshape(a_I_sum, params.n_I, []);
                SRNNModel2.plot_lines_with_colormap(t, a_I_summed, cmap_I);
                has_adaptation = true;
            end

            if ~isempty(a.E) && params.n_a_E > 0
                if has_adaptation, hold on; end
                a_E_sum = sum(a.E, 2);
                a_E_summed = reshape(a_E_sum, params.n_E, []);
                SRNNModel2.plot_lines_with_colormap(t, a_E_summed, cmap_E);
                has_adaptation = true;
            end

            if has_adaptation
                hold off;
                ylabel('adaptation');
            else
                text(0.5, 0.5, 'No adaptation variables', 'HorizontalAlignment', 'center');
                axis off;
            end
        end

        function plot_std_variable(t, b, params)
            % PLOT_STD_VARIABLE Plot short-term depression variables for E and I neurons.
            % Internalized from src/plotting/plot_std_variable.m
            cmap_I = SRNNModel2.inhibitory_colormap(8);
            cmap_E = SRNNModel2.excitatory_colormap(8);
            has_std = false;

            if params.n_b_I > 0 && ~isempty(b.I) && size(b.I, 1) > 0
                if ~all(b.I(:) == 1)
                    b_I_prod = reshape(prod(b.I, 2), params.n_I, []);   % collapse timescales
                    SRNNModel2.plot_lines_with_colormap(t, b_I_prod, cmap_I);
                    has_std = true;
                end
            end

            if params.n_b_E > 0 && ~isempty(b.E) && size(b.E, 1) > 0
                if ~all(b.E(:) == 1)
                    if has_std, hold on; end
                    b_E_prod = reshape(prod(b.E, 2), params.n_E, []);   % collapse timescales
                    SRNNModel2.plot_lines_with_colormap(t, b_E_prod, cmap_E);
                    has_std = true;
                end
            end

            if has_std
                hold off;
                ylabel('depression');
                ylim([0, 1]);
                yticks([0, 1]);
            else
                text(0.5, 0.5, 'No STD variables', 'HorizontalAlignment', 'center');
                axis off;
            end
        end

        function cmap = inhibitory_colormap(n_colors)
            % INHIBITORY_COLORMAP Custom colormap for inhibitory neurons.
            % Internalized from src/plotting/inhibitory_colormap.m
            if nargin < 1, n_colors = 8; end
            base_palette = [
                0.00, 0.45, 0.74;
                0.00, 0.75, 1.00;
                0.20, 0.47, 0.62;
                0.00, 0.50, 0.50;
                0.30, 0.75, 0.93;
                0.25, 0.62, 0.75;
                0.00, 0.80, 0.80;
                0.15, 0.55, 0.65;
                ];
            n_base = size(base_palette, 1);
            if n_colors == n_base
                cmap = base_palette;
            elseif n_colors < n_base
                indices = round(linspace(1, n_base, n_colors));
                cmap = base_palette(indices, :);
            else
                x_base = linspace(1, n_colors, n_base);
                x_new = 1:n_colors;
                cmap = zeros(n_colors, 3);
                for i = 1:3
                    cmap(:, i) = interp1(x_base, base_palette(:, i), x_new, 'pchip');
                end
                cmap = max(0, min(1, cmap));
            end
        end

        function cmap = excitatory_colormap(n_colors)
            % EXCITATORY_COLORMAP Custom colormap for excitatory neurons.
            % Internalized from src/plotting/excitatory_colormap.m
            if nargin < 1, n_colors = 8; end
            base_palette = [
                1.00, 0.00, 0.00;
                1.00, 0.75, 0.00;
                0.85, 0.20, 0.45;
                0.90, 0.10, 0.60;
                0.90, 0.55, 0.00;
                0.55, 0.27, 0.27;
                0.86, 0.08, 0.24;
                0.60, 0.15, 0.45;
                ];
            n_base = size(base_palette, 1);
            if n_colors == n_base
                cmap = base_palette;
            elseif n_colors < n_base
                indices = round(linspace(1, n_base, n_colors));
                cmap = base_palette(indices, :);
            else
                x_base = linspace(1, n_colors, n_base);
                x_new = 1:n_colors;
                cmap = zeros(n_colors, 3);
                for i = 1:3
                    cmap(:, i) = interp1(x_base, base_palette(:, i), x_new, 'pchip');
                end
                cmap = max(0, min(1, cmap));
            end
        end
    end
end

