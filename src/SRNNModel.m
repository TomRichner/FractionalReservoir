classdef SRNNModel < handle
    % SRNNMODEL Spiking Rate Neural Network Model class
    %
    % This class encapsulates all parameters, simulation logic, Lyapunov
    % analysis, and plotting for SRNN simulations with spike-frequency
    % adaptation (SFA) and short-term depression (STD).
    %
    % Usage:
    %   model = SRNNModel('n', 100, 'level_of_chaos', 1.8);
    %   model.build();
    %   model.run();
    %   model.plot();
    %
    % See also: SRNN_reservoir, create_W_matrix, compute_lyapunov_exponents

    %% Network Architecture Properties
    properties
        n = 100                     % Total number of neurons
        f = 0.5                     % Fraction of excitatory neurons
        indegree = 20               % Expected in-degree

        % RMT tilde-notation parameters (Harris 2023)
        mu_E_tilde                  % Normalized excitatory mean (default: 1/(alpha*sqrt(n)))
        mu_I_tilde                  % Normalized inhibitory mean (default: -1/(alpha*sqrt(n)))
        sigma_E_tilde               % Normalized excitatory std dev (default: 1/(alpha*sqrt(n)))
        sigma_I_tilde               % Normalized inhibitory std dev (default: 1/(alpha*sqrt(n)))
        E_W = 0                     % Mean offset: added to both mu_E_tilde and mu_I_tilde
        zrs_mode = 'none'           % ZRS mode: 'none', 'ZRS', 'SZRS', 'Partial_SZRS'

        level_of_chaos = 1.0        % Scaling factor for W
        rescale_by_abscissa = false % Whether to apply 1/abscissa_0 scaling
        row_zero_W = true           % Whether to apply row-mean centering to W (legacy, not used with RMTMatrix)
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
        n_b_E = 0                   % Number of STD timescales for E neurons (0 or 1)
        n_b_I = 0                   % Number of STD timescales for I neurons (0 or 1)
        tau_b_E_rec = 1             % STD recovery time constant for E neurons (s)
        tau_b_E_rel = 0.25          % STD release time constant for E neurons (s)
        tau_b_I_rec = 1             % STD recovery time constant for I neurons (s)
        tau_b_I_rel = 0.25          % STD release time constant for I neurons (s)
    end

    %% Dynamics Properties
    properties
        tau_d = 0.1                 % Dendritic time constant (s)
        activation_function         % Activation function handle
        activation_function_derivative  % Derivative of activation function
        S_a = 0.9                   % Activation function parameter a
        S_c = 0.35                  % Activation function parameter c (center)
    end

    %% Simulation Settings Properties
    properties
        fs = 200                    % Sampling frequency (Hz)
        T_range = [0, 50]           % Simulation time interval [start, end]
        T_plot                      % Plotting time interval (defaults to T_range)
        ode_solver = @ode45         % ODE solver function handle
        ode_opts                    % ODE solver options struct
    end

    %% Input Configuration Properties
    properties
        input_config                % Struct with stimulus parameters
        u_ex_scale = 1.0            % Scaling factor for external input
        rng_seeds = [1 2]    % RNG seeds [network, stimulus, etc]
    end

    %% Lyapunov Settings Properties
    properties
        lya_method = 'benettin'     % Lyapunov method: 'benettin', 'qr', or 'none'
        lya_T_interval              % Time interval for Lyapunov computation
        filter_local_lya = false    % Whether to filter local Lyapunov exponent
        lya_filter_order = 2        % Butterworth filter order
        lya_filter_cutoff = 0.25    % Normalized cutoff frequency (fraction of Nyquist)
    end

    %% Storage Options Properties
    properties
        store_full_state = false    % Whether to keep full S_out in memory
        store_decimated_state = true % Whether to keep decimated plot data
        plot_deci = 20              % Decimation factor for plotting
    end

    %% RMT Dependent Properties (computed from tilde parameters)
    properties (Dependent)
        alpha               % Sparsity = indegree/n
        mu_se               % Sparse excitatory mean
        mu_si               % Sparse inhibitory mean
        sigma_se            % Sparse excitatory std dev
        sigma_si            % Sparse inhibitory std dev
        R                   % Theoretical spectral radius (Harris 2023 Eq 18)
    end

    %% Computed Properties (SetAccess = protected for subclass access)
    properties (SetAccess = protected)
        W                           % Connection matrix (n x n)
        n_E                         % Number of excitatory neurons
        n_I                         % Number of inhibitory neurons
        E_indices                   % Indices of E neurons
        I_indices                   % Indices of I neurons
        N_sys_eqs                   % Total state dimension
        is_built = false            % Flag indicating if network is initialized

        % External input (generated by build)
        t_ex                        % Time vector for external input
        u_ex                        % External input matrix (n x nt)
        u_interpolant               % griddedInterpolant for external input (avoids persistent vars)
        S0                          % Initial state vector
        cached_params               % Cached params struct (set by build)
    end

    %% Results Properties (conditionally stored)
    properties (SetAccess = protected)
        t_out                       % Time vector from ODE solver
        S_out                       % State trajectory (nt x N_sys_eqs)
        plot_data                   % Struct with decimated data for plotting
        lya_results                 % Lyapunov analysis results struct
        has_run = false             % Flag indicating if simulation has run
    end

    %% Constructor
    methods
        function obj = SRNNModel(varargin)
            % SRNNMODEL Constructor with name-value pairs
            %
            % Usage:
            %   model = SRNNModel()  % All defaults
            %   model = SRNNModel('n', 200, 'level_of_chaos', 2.0)
            %   model = SRNNModel('n_a_E', 3, 'n_b_E', 1)

            % Set default values
            obj.set_defaults();

            % Parse name-value pairs
            for i = 1:2:length(varargin)
                if isprop(obj, varargin{i})
                    obj.(varargin{i}) = varargin{i+1};
                else
                    warning('SRNNModel:UnknownProperty', 'Unknown property: %s', varargin{i});
                end
            end
        end
    end

    %% Dependent Property Getters
    methods
        function val = get.alpha(obj)
            val = obj.indegree / obj.n;
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
                val = sqrt(obj.n * (obj.f * obj.sigma_se^2 + (1 - obj.f) * obj.sigma_si^2));
            end
        end
    end

    %% Public Methods
    methods
        function build(obj)
            % BUILD Initialize the network: compute derived params, create W, generate stimulus
            %
            % This method must be called before run(). It:
            %   1. Computes derived parameters (n_E, n_I, indices, etc.)
            %   2. Creates the weight matrix W using RMTMatrix
            %   3. Optionally scales W based on level_of_chaos
            %   4. Generates external stimulus u_ex
            %   5. Initializes the state vector S0

            % Set RNG seed for network generation
            rng(obj.rng_seeds(1));

            % Compute derived parameters
            obj.compute_derived_params();

            % Compute RMT tilde defaults if not set
            alpha_val = obj.indegree / obj.n;
            % Default val D derived to give R=1: 1 = D * sqrt(n * alpha * (2 - alpha))
            default_val = 1 / sqrt(obj.n * alpha_val * (2 - alpha_val));

            if isempty(obj.mu_E_tilde),    obj.mu_E_tilde = default_val;     end
            if isempty(obj.mu_I_tilde),    obj.mu_I_tilde = -default_val;    end
            if isempty(obj.sigma_E_tilde), obj.sigma_E_tilde = default_val;  end
            if isempty(obj.sigma_I_tilde), obj.sigma_I_tilde = default_val;  end

            % Compute tau_a arrays if n_a > 0 but tau_a not set
            if obj.n_a_E > 0 && isempty(obj.tau_a_E)
                obj.tau_a_E = logspace(log10(0.25), log10(10), obj.n_a_E);
            end
            if obj.n_a_I > 0 && isempty(obj.tau_a_I)
                obj.tau_a_I = logspace(log10(0.25), log10(10), obj.n_a_I);
            end

            % Create W matrix using RMTMatrix
            rmt = RMTMatrix(obj.n);
            rmt.alpha = alpha_val;
            rmt.f = obj.f;
            rmt.mu_tilde_e = obj.mu_E_tilde + obj.E_W;   % Apply mean offset
            rmt.mu_tilde_i = obj.mu_I_tilde + obj.E_W;   % Apply mean offset
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
                max(abs(W_eigs_scaled)), max(real(W_eigs_scaled)), rmt.R * obj.level_of_chaos);

            % Generate external stimulus
            obj.generate_stimulus();

            % Build interpolant for external input (avoids persistent variables)
            obj.u_interpolant = griddedInterpolant(obj.t_ex, obj.u_ex', 'linear', 'none');

            % Initialize state vector
            params_init = obj.get_params();
            obj.S0 = initialize_state(params_init);

            % Validate configuration
            obj.validate();

            % Cache params struct for fast access in run/plot methods
            obj.cached_params = obj.get_params();

            obj.is_built = true;
            fprintf('Model built successfully. Ready to run.\n');
        end

        function run(obj)
            % RUN Execute the ODE simulation
            %
            % This method integrates the SRNN equations and optionally computes
            % Lyapunov exponents. Results are stored based on storage options.

            if ~obj.is_built
                error('SRNNModel:NotBuilt', 'Model must be built before running. Call build() first.');
            end

            % Use cached params struct
            params = obj.cached_params;

            % Set up ODE options
            dt = 1 / obj.fs;
            if isempty(obj.ode_opts)
                jac_wrapper = @(t, S) compute_Jacobian_fast(S, params);
                obj.ode_opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-7, 'MaxStep', dt, 'Jacobian', jac_wrapper);
            end

            % Define RHS function using closure (avoids OOP overhead)
            % Cache interpolant and params to avoid property access on every call
            u_interp = obj.u_interpolant;
            params.u_interpolant = u_interp;  % Add to params for dynamics_fast
            rhs = @(t, S) SRNNModel.dynamics_fast(t, S, params);

            % Integrate
            fprintf('Integrating equations\n');
            tic
            [t_raw, S_raw] = obj.ode_solver(rhs, obj.t_ex, obj.S0, obj.ode_opts);
            integration_time = toc;
            fprintf('Integration complete in %.2f seconds.\n', integration_time);

            % Verify output times match input times
            if length(t_raw) ~= length(obj.t_ex) || max(abs(t_raw - obj.t_ex)) > 1e-9
                error('SRNNModel:TimeMismatch', 'ODE solver output times do not match input times. Max diff: %.2e', max(abs(t_raw(:) - obj.t_ex(:))));
            end

            % Store temporarily for Lyapunov and decimation
            obj.t_out = t_raw;
            obj.S_out = S_raw;

            % Compute Lyapunov exponents
            if ~strcmpi(obj.lya_method, 'none')
                obj.compute_lyapunov();
                % Filter local Lyapunov exponent (before decimation to avoid edge effects)
                if obj.filter_local_lya
                    obj.filter_lyapunov();
                end
            end

            % Decimate and unpack for plotting
            if obj.store_decimated_state
                obj.decimate_and_unpack();
            end

            % Clear full state if not storing (free memory)
            if ~obj.store_full_state
                obj.S_out = [];
            end

            obj.has_run = true;
            fprintf('Simulation complete.\n');
        end

        function compute_lyapunov(obj)
            % COMPUTE_LYAPUNOV Compute Lyapunov exponents based on lya_method

            if isempty(obj.S_out)
                error('SRNNModel:NoStateData', 'State data not available. Set store_full_state=true or call before clearing.');
            end

            dt = 1 / obj.fs;
            params = obj.cached_params;

            % Set Lyapunov time interval
            if isempty(obj.lya_T_interval)
                obj.lya_T_interval = [obj.T_range(1) + 2, obj.T_range(2)]; % by default skip the first two seconds
            end

            % Define RHS function using closure (avoids OOP overhead)
            params.u_interpolant = obj.u_interpolant;
            rhs = @(t, S) SRNNModel.dynamics_fast(t, S, params);

            fprintf('Computing Lyapunov exponents using %s method\n', obj.lya_method);
            obj.lya_results = compute_lyapunov_exponents(obj.lya_method, obj.S_out, obj.t_out, dt, obj.fs, obj.lya_T_interval, params, obj.ode_opts, obj.ode_solver, rhs, obj.t_ex, obj.u_ex);

            if isfield(obj.lya_results, 'LLE')
                fprintf('Largest Lyapunov Exponent: %.4f\n', obj.lya_results.LLE);
            end
        end

        function filter_lyapunov(obj)
            % FILTER_LYAPUNOV Apply lowpass filter to local Lyapunov exponent
            %
            % This filters the local_lya signal BEFORE decimation to avoid
            % edge effects that would occur if filtering after trimming.
            % Uses a Butterworth filter with parameters from obj properties.

            if isempty(obj.lya_results)
                return;
            end

            % Design filter (cutoff in Hz, normalized by Nyquist = lya_fs/2)
            Wn = obj.lya_filter_cutoff / (obj.lya_results.lya_fs / 2);
            [b, a] = butter(obj.lya_filter_order, Wn, 'low');

            % Filter local_lya for Benettin method
            if isfield(obj.lya_results, 'local_lya') && ~isempty(obj.lya_results.local_lya)
                obj.lya_results.local_lya = filtfilt(b, a, obj.lya_results.local_lya);
            end

            % Filter local_LE_spectrum_t for QR method (each column)
            if isfield(obj.lya_results, 'local_LE_spectrum_t') && ~isempty(obj.lya_results.local_LE_spectrum_t)
                for col = 1:size(obj.lya_results.local_LE_spectrum_t, 2)
                    obj.lya_results.local_LE_spectrum_t(:, col) = filtfilt(b, a, obj.lya_results.local_LE_spectrum_t(:, col));
                end
            end
        end

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

            [fig_handle, ax_handles] = plot_SRNN_tseries(obj.plot_data.t, obj.plot_data.u, obj.plot_data.x, obj.plot_data.r, obj.plot_data.a, obj.plot_data.b, obj.plot_data.br, params, obj.lya_results, obj.lya_method, T_plot_arg);
        end

        function [fig_handle, ax_handles] = plot_eigenvalues(obj, J_times_sec)
            % PLOT_EIGENVALUES Plot Jacobian eigenvalues at specified times
            %
            % Usage:
            %   model.plot_eigenvalues([5, 10, 15])  % Times in seconds

            if isempty(obj.S_out)
                error('SRNNModel:NoStateData', 'State data required. Set store_full_state=true.');
            end

            params = obj.cached_params;

            % Convert times to indices
            J_times = round((J_times_sec - obj.t_out(1)) * obj.fs) + 1;
            J_times = unique(max(1, min(J_times, size(obj.S_out, 1))));

            fprintf('Computing Jacobian at %d time points\n', length(J_times));
            J_array = compute_Jacobian_at_indices(obj.S_out, J_times, params);

            % Compute eigenvalues
            n_plots = length(J_times);
            eigenvalues_all = cell(n_plots, 1);
            for i = 1:n_plots
                eigenvalues_all{i} = eig(J_array(:,:,i));
            end

            % Determine subplot layout
            if n_plots <= 4
                n_rows = 1;
                n_cols = n_plots;
            else
                n_cols = ceil(sqrt(n_plots));
                n_rows = ceil(n_plots / n_cols);
            end

            % Compute global axis limits
            all_real = [];
            all_imag = [];
            for i = 1:n_plots
                evals = eigenvalues_all{i};
                all_real = [all_real; real(evals)];
                all_imag = [all_imag; imag(evals)];
            end
            global_xlim = [min(all_real), max(all_real)];
            global_ylim = [min(all_imag), max(all_imag)];
            x_range = diff(global_xlim);
            y_range = diff(global_ylim);
            global_xlim = global_xlim + [-0.1, 0.1] * x_range;
            global_ylim = global_ylim + [-0.1, 0.1] * y_range;

            % Create figure
            fig_handle = figure('Position', [1312, 526, 600, 360]);
            ax_handles = zeros(n_plots, 1);

            for i = 1:n_plots
                ax = subplot(n_rows, n_cols, i);
                evals = eigenvalues_all{i};
                time_val = obj.t_out(J_times(i));
                ax_handles(i) = plot_eigenvalues(evals, ax, time_val, global_xlim, global_ylim);
                set(ax_handles(i), 'Color', 'none');
            end

            linkaxes(ax_handles, 'xy');
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
            R_W = obj.R * obj.level_of_chaos;  % Scale by level_of_chaos since W is scaled
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

            % Dynamics
            params.tau_d = obj.tau_d;
            params.activation_function = obj.activation_function;
            params.activation_function_derivative = obj.activation_function_derivative;

            % Connection matrix (if built)
            if ~isempty(obj.W)
                params.W = obj.W;
            end

            % RNG seeds
            params.rng_seeds = obj.rng_seeds;
        end

        function clear_results(obj)
            % CLEAR_RESULTS Free memory by clearing stored state data

            obj.t_out = [];
            obj.S_out = [];
            obj.plot_data = [];
            obj.lya_results = [];
            obj.has_run = false;
            fprintf('Results cleared.\n');
        end

        function reset(obj)
            % RESET Clear built state to allow rebuilding with new parameters
            %
            % Usage:
            %   model.reset();
            %   model.n = 200;  % Change parameters
            %   model.build();  % Rebuild with new settings

            obj.is_built = false;
            obj.W = [];
            obj.u_interpolant = [];
            obj.t_ex = [];
            obj.u_ex = [];
            obj.S0 = [];
            obj.cached_params = [];
            obj.clear_results();
            fprintf('Model reset. Modify parameters and call build() to reinitialize.\n');
        end

        function dS_dt = dynamics(obj, t, S)
            % DYNAMICS Compute the right-hand side of the SRNN ODE system
            %
            % This is a convenience method that wraps dynamics_fast.
            % For performance-critical code (e.g., ODE integration), use
            % dynamics_fast directly with a params struct.

            params = obj.cached_params;
            params.u_interpolant = obj.u_interpolant;
            dS_dt = SRNNModel.dynamics_fast(t, S, params);
        end
    end

    %% Private/Protected Methods
    methods (Access = private)
        function set_defaults(obj)
            % SET_DEFAULTS Initialize all properties to default values

            % Set default activation function (piecewiseSigmoid)
            obj.activation_function = @(x) piecewiseSigmoid(x, obj.S_a, obj.S_c);
            obj.activation_function_derivative = @(x) piecewiseSigmoidDerivative(x, obj.S_a, obj.S_c);

            % Set default input configuration
            obj.input_config = struct();
            obj.input_config.n_steps = 3;
            obj.input_config.step_density = 0.2;
            obj.input_config.amp = 0.5;
            obj.input_config.no_stim_pattern = false(1, 3);
            obj.input_config.no_stim_pattern(1:2:end) = true;
            obj.input_config.intrinsic_drive = [];  % Will be set in build
            obj.input_config.positive_only = false;  % Default: allow positive and negative amplitudes
            % step_density_E/I not set by default (uses step_density uniformly)

            % T_plot defaults to T_range (set in build if not specified)
            obj.T_plot = [];
        end

        function compute_derived_params(obj)
            % COMPUTE_DERIVED_PARAMS Calculate n_E, n_I, indices, N_sys_eqs

            obj.n_E = round(obj.f * obj.n);
            obj.n_I = obj.n - obj.n_E;
            obj.E_indices = 1:obj.n_E;
            obj.I_indices = obj.n_E + 1:obj.n;

            % Calculate total state dimension
            obj.N_sys_eqs = obj.n_E * obj.n_a_E + obj.n_I * obj.n_a_I + obj.n_E * obj.n_b_E + obj.n_I * obj.n_b_I + obj.n;
        end

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
            [u_stim, t_stim] = generate_external_input(params_stim, T_stim, obj.fs, obj.rng_seeds(2), obj.input_config);

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
            [t_plot, S_plot, plot_indices] = decimate_states(obj.t_out, obj.S_out, obj.plot_deci);

            % Unpack state vector and compute firing rates
            [x_plot, a_plot, b_plot, r_plot, br_plot] = unpack_and_compute_states(S_plot, params);

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
            b_plot = obj.trim_struct_data(b_plot, 2, keep_mask);
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

        function plot_spectrum_helper(~, ax, eigs, R, center, outlier_threshold)
            % PLOT_SPECTRUM_HELPER Helper function for plotting eigenvalue spectra
            %
            % Inputs:
            %   ax               - Axes handle to plot on
            %   eigs             - Vector of eigenvalues
            %   R                - Theoretical spectral radius
            %   center           - Center of the spectral disc (real part)
            %   outlier_threshold - Multiplier for R to classify far outliers

            % Compute distances from center for all eigenvalues
            distances = abs(eigs - center);

            % Plot interior eigenvalues (within R) as black circles
            mSize = 4;
            interior_mask = distances <= R;
            interior_eigs = eigs(interior_mask);
            plot(ax, real(interior_eigs), imag(interior_eigs), 'ko', 'MarkerSize', mSize, 'MarkerFaceColor', 'none', 'LineWidth', 0.5);
            hold(ax, 'on');

            % Plot theoretical radius circle (Eq 18)
            theta = linspace(0, 2*pi, 100);
            plot(ax, center + R*cos(theta), R*sin(theta), 'k-', 'LineWidth', 2);

            % Plot near outlier eigenvalues (between R and outlier_threshold*R) as black Xs
            near_outlier_mask = (distances > R) & (distances <= outlier_threshold * R);
            near_outlier_eigs = eigs(near_outlier_mask);
            if ~isempty(near_outlier_eigs)
                plot(ax, real(near_outlier_eigs), imag(near_outlier_eigs), 'kx', 'MarkerSize', mSize, 'LineWidth', 0.5);
            end

            % Plot far outlier eigenvalues (beyond outlier_threshold*R) as green filled circles
            far_outlier_mask = distances > outlier_threshold * R;
            far_outlier_eigs = eigs(far_outlier_mask);
            if ~isempty(far_outlier_eigs)
                plot(ax, real(far_outlier_eigs), imag(far_outlier_eigs), 'o', 'MarkerSize', mSize, 'MarkerFaceColor', [0 .7 0], 'MarkerEdgeColor', [0 .7 0]);
            end

            % Add axis lines through origin
            xl = xlim(ax);
            yl = ylim(ax);
            plot(ax, xl, [0 0], 'k-', 'LineWidth', 0.5);
            plot(ax, [0 0], yl, 'k-', 'LineWidth', 0.5);

            grid(ax, 'on');
            axis(ax, 'equal');
            hold(ax, 'off');
        end
    end

    methods (Static, Access = protected)
        function s_out = trim_struct_data(s_in, dim, mask)
            % TRIM_STRUCT_DATA Helper to trim fields of a struct along a dimension
            s_out = s_in;
            fields = fieldnames(s_in);
            for i = 1:length(fields)
                val = s_in.(fields{i});
                if ~isempty(val)
                    if dim == 2
                        s_out.(fields{i}) = val(:, mask);
                    elseif dim == 3
                        s_out.(fields{i}) = val(:, :, mask);
                    end
                end
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

            % STD states for E neurons (b_E)
            len_b_E = n_E * n_b_E;
            if len_b_E > 0
                b_E = S(current_idx + (1:len_b_E));
            else
                b_E = [];
            end
            current_idx = current_idx + len_b_E;

            % STD states for I neurons (b_I)
            len_b_I = n_I * n_b_I;
            if len_b_I > 0
                b_I = S(current_idx + (1:len_b_I));
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

            % Apply STD effect
            b = ones(n, 1);
            if n_b_E > 0 && ~isempty(b_E)
                b(E_indices) = b_E;
            end
            if n_b_I > 0 && ~isempty(b_I)
                b(I_indices) = b_I;
            end

            r = activation_fn(x_eff);

            %% Compute derivatives
            dx_dt = (-x + W * (b .* r) + u) / tau_d;

            da_E_dt = [];
            if n_E > 0 && n_a_E > 0 && ~isempty(a_E)
                da_E_dt = (r(E_indices) - a_E) ./ tau_a_E;
            end

            da_I_dt = [];
            if n_I > 0 && n_a_I > 0 && ~isempty(a_I)
                da_I_dt = (r(I_indices) - a_I) ./ tau_a_I;
            end

            db_E_dt = [];
            if n_E > 0 && n_b_E > 0 && ~isempty(b_E)
                db_E_dt = (1 - b_E) / tau_b_E_rec - (b_E .* r(E_indices)) / tau_b_E_rel;
            end

            db_I_dt = [];
            if n_I > 0 && n_b_I > 0 && ~isempty(b_I)
                db_I_dt = (1 - b_I) / tau_b_I_rec - (b_I .* r(I_indices)) / tau_b_I_rel;
            end

            %% Pack derivatives
            dS_dt = [da_E_dt(:); da_I_dt(:); db_E_dt(:); db_I_dt(:); dx_dt];
        end
    end
end

