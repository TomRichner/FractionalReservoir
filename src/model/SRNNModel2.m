classdef SRNNModel2 < handle
    % SRNNMODEL2 Spiking Rate Neural Network Model class (FractionalResevoir version)
    %
    % Copied from ConnectivityAdaptation/StabilityAnalysis/src/SRNNModel.m
    % Renamed to SRNNModel2 for use in FractionalResevoir project.
    %
    % DIVERGENT INTERNAL FUNCTIONS (to avoid path conflicts):
    %   - plot_eigenvalues_helper (static) - eigenvalue scatter plot
    %   - decimate_states (static) - state decimation for plotting
    %   - initialize_state (static) - state vector initialization
    %   - filter_lyapunov (instance) - local Lyapunov filtering
    %   - compute_lyapunov_exponents_internal (static) - main Lyapunov dispatcher
    %   - benettin_algorithm_internal (static) - Benettin's LLE algorithm
    %   - lyapunov_spectrum_qr_internal (static) - QR full spectrum method
    %   - lyapunov_sample_grid (static) - shared Benettin/QR segment grid
    %   - compute_kaplan_yorke_dimension_internal (static) - Kaplan-Yorke dimension
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
        n = 300                     % Total number of neurons
        f = 0.5                     % Fraction of excitatory neurons
        indegree = 100              % Expected in-degree
        
        % RMT tilde-notation parameters (Harris 2023), expressed as multiples of
        % the normalization factor F = default_val. The absolute mu_*_tilde /
        % sigma_*_tilde / E_W are Dependent properties computed from these; see
        % the Dependent block below. Storing the multiplier rather than the
        % product is what lets these be swept (F depends on n and indegree) and
        % frozen into a run's recorded parameter set.
        mu_E_tilde_relative    =  3     % Normalized excitatory mean, in multiples of F
        mu_I_tilde_relative    = -4     % Normalized inhibitory mean, in multiples of F
        sigma_E_tilde_relative =  1     % Normalized excitatory std dev, in multiples of F
        sigma_I_tilde_relative =  1     % Normalized inhibitory std dev, in multiples of F
        E_W_relative           =  0     % Mean offset (multiples of F): added to both mu tildes

        % BLOCK overrides, indexed (POSTsynaptic, PREsynaptic) as in RMTBlocks.
        % Empty means "use the column shorthand above", i.e. the value depends
        % only on the presynaptic type, which is the classical Harris (2023)
        % setting. Setting one of these lets, say, E->E differ from E->I:
        %
        %   mu_EE = E receives from E      mu_EI = E receives from I
        %   mu_IE = I receives from E      mu_II = I receives from I
        %
        % Dale's law is a COLUMN constraint, so signs must stay constant down
        % each column (mu_EE and mu_IE share the sign of mu_E_tilde_relative).
        % Like the shorthands, these are multiples of F.
        mu_EE_relative    = []
        mu_EI_relative    = []
        mu_IE_relative    = []
        mu_II_relative    = []
        sigma_EE_relative = []
        sigma_EI_relative = []
        sigma_IE_relative = []
        sigma_II_relative = []
        zrs_mode = 'none'           % ZRS mode: 'none', 'ZRS', 'SZRS', 'Partial_SZRS'

        % Which network F is computed for. When true (default), F follows the
        % CURRENT n and indegree, which makes the theoretical spectral radius R
        % exactly independent of n. When false, F is pinned to the reference
        % network below, so sweeping n or indegree leaves the weight
        % distribution fixed and lets R vary with network size. See default_val.
        F_tracks_network = true     % F follows the current n/indegree
        F_ref_n          = 300      % Reference network size for F when F_tracks_network is false
        F_ref_indegree   = 100      % Reference in-degree for F when F_tracks_network is false
        
        level_of_chaos = 1.0        % Scaling factor for W
        rescale_by_abscissa = false % Whether to apply 1/abscissa_0 scaling
        check_connectivity = false  % If true, build_network reports the digraph connectivity class of W
    end
    
    %% Spike-Frequency Adaptation (SFA) Properties
    properties
        n_a_E = 0                   % Number of adaptation timescales for E neurons
        n_a_I = 0                   % Number of adaptation timescales for I neurons
        tau_a_E                     % Adaptation time constants for E neurons (1 x n_a_E)
        tau_a_I                     % Adaptation time constants for I neurons (1 x n_a_I)
        c_E = 0.15    % TOTAL adaptation budget for E neurons (divided by n_a_E)
        c_I = 0.15    % TOTAL adaptation budget for I neurons (divided by n_a_I)
        % Was 0.15/3, hand-divided on the assumption n_a = 3. The model now
        % divides by the number of timescales actually in use (effective_c), so
        % keeping the hand-division would apply it twice.
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
    
    %% Dynamics Properties
    properties
        tau_d = 0.1                 % Dendritic time constant (s)

        % The nonlinearity is chosen by NAME and parameterised by S_a / S_c; the
        % activation_function / activation_function_derivative handles are
        % Dependent and rebuilt from these (see the Dependent block below).
        % Storing the choice as data rather than as a handle is what keeps S_a
        % and S_c from silently disagreeing with the function actually in use,
        % and lets a parameter preset express the nonlinearity without handles.
        %
        % 'logistic'  - logisticSigmoid(x, S_c); smooth, unit slope at centre.
        %               The default: more robust near-edge-of-chaos stability
        %               than the piecewise variant. Ignores S_a.
        % 'piecewise' - piecewiseSigmoid(x, S_a, S_c); piecewise linear/quadratic.
        % 'tanh'      - tanhActivation(x); ignores both S_a and S_c.
        activation = 'logistic'     % 'logistic' | 'piecewise' | 'tanh'
        S_a = 0.9                   % Activation function parameter a (piecewise only)
        S_c = 0.4                   % Activation function parameter c, the centre (piecewise, logistic)

        % Per-neuron setpoints. By default every neuron shares the scalar S_c
        % above. Give a population a mean and/or a standard deviation here and
        % build() draws a length-n S_c vector instead:
        %
        %   S_c_i = mu_S_c_<pop> + sigma_S_c_<pop> * randn,   pop = E or I
        %
        % An empty mu falls back to S_c, so setting sigma alone adds spread
        % around the shared centre, and setting mu alone gives E and I
        % different excitability with no spread. The realized vector is
        % S_c_vec (read-only, empty unless heterogeneity was requested).
        %
        % Only 'logistic' and 'piecewise' take a centre: requesting
        % heterogeneity with 'tanh' or with activation_custom is an error, not
        % a silent no-op. Note that in heterogeneous mode the
        % activation_function handle is only valid for length-n input.
        mu_S_c_E = []               % Mean setpoint for E neurons ([] = use S_c)
        mu_S_c_I = []               % Mean setpoint for I neurons ([] = use S_c)
        sigma_S_c_E = 0             % Std dev of the E setpoint (0 = no spread)
        sigma_S_c_I = 0             % Std dev of the I setpoint (0 = no spread)

        % Seed for the setpoint draw. Empty derives it from rng_seeds(1) (so a
        % reps sweep, which varies the network seed, also varies the setpoints)
        % while keeping it well away from the stream that generated W. Set it
        % to pin the draw across networks.
        S_c_seed = []

        % Escape hatch for a nonlinearity that is not one of the named three:
        % a 1x2 cell {fn, dfn} of function handles. When non-empty it overrides
        % `activation`. Leave empty for the normal case, which keeps the
        % recorded parameter set free of function handles.
        activation_custom = {}
    end
    
    %% Simulation Settings Properties
    properties
        fs = 400                    % Sampling frequency (Hz)
        T_range = [0, 50]           % Simulation time interval [start, end]
        T_plot                      % Plotting time interval (defaults to T_range)
        % The integrator is chosen BY NAME, not by passing a handle, for the
        % same reason `activation` is: the name is data, so it survives into
        % resolved_defaults, compares cleanly in same_config, and a sweep can
        % carry it without smuggling a function handle through a struct.
        %
        %   'ode45'  - adaptive Runge-Kutta 4(5) (the default)
        %   'ode15s' - adaptive stiff solver
        %   'rk4'    - fixed-step classic RK4 (ode_rk4); much faster than
        %              ode45 with noisy forcing, since there is no step-size
        %              control to thrash
        %
        % Passing a function handle raises SRNNModel:RenamedProperty naming the
        % string that replaces it.
        %
        % NOTE the QR Lyapunov method always integrates its variational
        % equation with ode45 regardless of this setting -- see compute_lyapunov.
        ode_solver = 'ode45'        % 'ode45' | 'ode15s' | 'rk4'
        ode_opts                    % ODE solver options struct
        x0_std = 0.1                % Std dev of the random initial dendritic state x0 (0 = deterministic x0=0)

        % Additive Wiener noise on the dendritic state:
        %
        %   dx_i = (-x_i + sum_j w_ij b_j r_j + u_i)/tau_d dt + sigma_u_noise/tau_d dW_i
        %
        % sigma_u_noise is INPUT-REFERRED: it is in the same units as u, so it
        % is directly comparable to intrinsic_drive and the stimulus amplitude
        % ("the noise is 20% of the DC drive"), and migrating the older
        % noise-through-the-stimulus route is the identity rather than a
        % conversion. The raw diffusion coefficient the integrator actually
        % multiplies dW by is the Dependent sigma_x_raw = sigma_u_noise/tau_d;
        % x_noise_std reports the same thing as a stationary std of x.
        %
        % Noise enters ONLY x. That keeps the diffusion constant, which is what
        % makes Ito and Stratonovich coincide, kills the Milstein term, leaves
        % the QR variational equation untouched, and -- because the two
        % trajectories share one noise path -- makes the additive noise cancel
        % in Benettin's difference, so the LLE stays measurable at any noise
        % level rather than only at small ones.
        %
        % sigma_u_noise > 0 REQUIRES a stochastic integrator ('euler', 'heun'
        % or 'sra1'); asking for noise with ode45/ode15s/rk4 is an error, not a
        % silent drop. sigma_u_noise = 0 skips noise generation entirely and is
        % bit-identical to the deterministic model.
        sigma_u_noise = 0

        % Seed for the noise draw. Empty derives it from rng_seeds(1) with an
        % offset well away from the streams that build W and the S_c setpoints,
        % so a reps sweep varies the noise realisation too. Under
        % ParamSpaceAnalysis2 that gives the same noise path across all
        % conditions at a grid point (they share a network seed) while varying
        % it across reps and across grid points. Set it to pin the realisation.
        noise_seed = []
    end
    
    %% Input Configuration Properties
    properties
        input_config                % Struct with stimulus parameters
        u_ex_scale = 1.0            % Scaling factor for external input
        rng_seeds = [1 2]           % RNG seeds [network, stimulus, etc]
        reps = 1                    % Repetition index (reserved for future use; typically unused)
    end
    
    %% Lyapunov Settings Properties
    properties
        lya_method = 'benettin'     % Lyapunov method: 'benettin', 'qr', or 'none'
        lya_T_interval              % Time interval for Lyapunov computation

        % Rescaling interval (s): how far the perturbation is allowed to grow
        % (Benettin) or the variational flow to run (QR) before renormalising.
        % EMPTY takes the per-method default -- 0.02 for 'benettin', 0.1 for
        % 'qr' -- which is what the hardcoded values used to be. The two differ
        % because QR renormalises a whole N x N basis per step and Benettin only
        % one vector, so QR can afford far fewer, larger steps.
        %
        % Must be an integer multiple of 1/fs and at least 3/fs; both are
        % checked, with the failing ratio reported. Larger values are cheaper
        % but coarser: the local exponent is a mean over the interval, and an
        % interval long enough for the perturbation to leave the linear regime
        % biases the estimate. Note lya_warmup is in seconds, so changing
        % lya_dt does not change how much alignment time the warmup buys.
        lya_dt = []

        % Seconds of Lyapunov iteration BEFORE lya_T_interval(1) used to let the
        % Benettin perturbation (or the QR basis Q) align with the leading
        % direction. Nothing is accumulated during it, so it costs compute and
        % buys alignment. Iteration starts at
        %
        %   max(T_range(1), lya_T_interval(1) - lya_warmup)
        %
        % i.e. the request is CLAMPED to the data that exists rather than
        % abandoned: if the simulation does not start early enough, the warmup
        % is shortened (with a warning) instead of dropped, since the samples
        % before the accumulation window are free alignment either way.
        %
        % TOO SHORT A WARMUP BIASES THE EXPONENT, and the two methods do not
        % need the same amount. Measured on the test_benettin_vs_qr network
        % (n=40, indegree=20, n_a_E=3, n_b_E=1, T_range=[-10,15]):
        %
        %   warmup:   0.5s     1s      2s      5s     10s     25s
        %   Benettin -0.1318 -0.1177 -0.1062 -0.1015 -0.1013 -0.1013
        %   QR            -       -  -0.2197 -0.1121 -0.1004 -0.1003
        %
        % Convergence tracks PHYSICAL time (the spectral gap), not the number
        % of renormalisations: QR reaches its plateau at the same ~10 s as
        % Benettin despite taking 5x fewer steps to get there. The default of 5
        % is comfortable for Benettin (0.2% from plateau) but NOT converged for
        % QR (~12% out) -- raise it for QR work, and expect systems with slow
        % adaptation (the tau_a sweeps reach 60 s) to need more still.
        lya_warmup = 5

        filter_local_lya = false    % Whether to filter local Lyapunov exponent
        lya_filter_order = 2        % Butterworth filter order
        lya_filter_cutoff = 0.25    % Normalized cutoff frequency (fraction of Nyquist)
    end
    
    %% Storage Options Properties
    properties
        store_full_state = false    % Whether to keep full S_out in memory
        store_decimated_state = true % Whether to keep decimated plot data
        plot_deci                   % Decimation factor for plotting (computed from fs/plot_freq if not set)
        plot_freq = 10              % Target plotting frequency (Hz)
    end
    
    %% RMT Dependent Properties (computed from tilde parameters)
    properties (Dependent)
        alpha               % Sparsity = indegree/n
        default_val         % Normalization factor F = 1/sqrt(N*alpha*(2-alpha)) (Harris 2023)
        mu_E_tilde          % Normalized excitatory mean    = mu_E_tilde_relative    * F
        mu_I_tilde          % Normalized inhibitory mean    = mu_I_tilde_relative    * F
        sigma_E_tilde       % Normalized excitatory std dev = sigma_E_tilde_relative * F
        sigma_I_tilde       % Normalized inhibitory std dev = sigma_I_tilde_relative * F
        E_W                 % Mean offset added to both mu tildes = E_W_relative * F
        mu_tilde_block      % 2x2 absolute block means (post, pre), including E_W
        sigma_tilde_block   % 2x2 absolute block std devs (post, pre)
        lambda_O            % Outlier eigenvalues of E[W], by descending magnitude
        activation_function            % Nonlinearity handle, built from activation + S_a/S_c
        activation_function_derivative % Its derivative, built the same way
        sigma_x_raw         % Raw diffusion coefficient = sigma_u_noise / tau_d
        x_noise_std         % Nominal stationary std of x from noise alone (W=0, u=0)
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
    
    %% Computed Properties (SetAccess = protected for subclass access)
    properties (SetAccess = protected)
        W                           % Connection matrix (n x n)
        is_built = false            % Flag indicating if network is initialized
        
        % External input (generated by build)
        t_ex                        % Time vector for external input
        u_ex                        % External input matrix (n x nt)
        u_interpolant               % griddedInterpolant for external input (avoids persistent vars)
        S0                          % Initial state vector
        cached_params               % Cached params struct (set by build)

        % Pre-generated Wiener increments for the current run, or empty. Built
        % at the top of run() rather than in build(): it is large (n x nt x 2
        % doubles, ~96 MB at n=300, T=50 s, fs=400), it is only needed while
        % integrating, and keeping it out of build() keeps it out of the
        % build-comparison in verify_shared_build. It stays alive through
        % compute_lyapunov -- Benettin re-integrates segments and must see the
        % same increments -- and is cleared alongside S_out. Never saved: it is
        % regenerable from noise_seed.
        noise_increments = []

        % Realized per-neuron nonlinearity setpoints (n x 1), drawn by
        % build_setpoints() from mu_S_c_* / sigma_S_c_*. EMPTY means every
        % neuron shares the scalar S_c, which is the default and keeps every
        % code path on the scalar branch.
        S_c_vec = []
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
        function obj = SRNNModel2(varargin)
            % SRNNMODEL2 Constructor with name-value pairs
            %
            % Usage:
            %   model = SRNNModel2()  % All defaults (n=300, indegree=100)
            %   model = SRNNModel2('n', 200, 'level_of_chaos', 2.0)
            %   model = SRNNModel2('n_a_E', 3, 'n_b_E', 1)
            
            % Set default values
            obj.set_defaults();
            
            % These properties became Dependent (computed from the settable ones
            % named here), so assigning them would raise MATLAB's generic
            % read-only error. Name the replacement and how to use it instead --
            % isprop() is true for Dependent properties, so check this first.
            renamed = struct( ...
                'mu_E_tilde', struct('set', 'mu_E_tilde_relative', ...
                    'hint', 'which holds the value in multiples of F = default_val (''mu_E_tilde_relative'', 3 replaces ''mu_E_tilde'', 3*F)'), ...
                'mu_I_tilde', struct('set', 'mu_I_tilde_relative', ...
                    'hint', 'which holds the value in multiples of F = default_val'), ...
                'sigma_E_tilde', struct('set', 'sigma_E_tilde_relative', ...
                    'hint', 'which holds the value in multiples of F = default_val'), ...
                'sigma_I_tilde', struct('set', 'sigma_I_tilde_relative', ...
                    'hint', 'which holds the value in multiples of F = default_val'), ...
                'E_W', struct('set', 'E_W_relative', ...
                    'hint', 'which holds the value in multiples of F = default_val'), ...
                'activation_function', struct('set', 'activation', ...
                    'hint', 'one of ''logistic'', ''piecewise'', ''tanh'', parameterised by S_a/S_c (or set activation_custom for a bespoke nonlinearity)'), ...
                'activation_function_derivative', struct('set', 'activation', ...
                    'hint', 'the derivative follows the same choice; there is no separate setting'), ...
                'sigma_x_raw', struct('set', 'sigma_u_noise', ...
                    'hint', 'which holds the noise INPUT-REFERRED, in the units of u (sigma_x_raw = sigma_u_noise/tau_d)'), ...
                'x_noise_std', struct('set', 'sigma_u_noise', ...
                    'hint', 'which holds the noise input-referred; x_noise_std = sigma_u_noise/sqrt(2*tau_d) is the stationary std it implies'));

            % Parse name-value pairs
            for i = 1:2:length(varargin)
                name = varargin{i};
                if ischar(name) && isfield(renamed, name)
                    error('SRNNModel:RenamedProperty', ...
                        '''%s'' is now a computed (Dependent) property. Set ''%s'' instead, %s.', ...
                        name, renamed.(name).set, renamed.(name).hint);
                elseif ischar(name) && strcmp(name, 'ode_solver')
                    % Fail here rather than at the first solver call, and give
                    % the handle-to-name mapping while the caller can see it.
                    check_ode_solver(varargin{i+1}, 'SRNNModel');
                    obj.ode_solver = varargin{i+1};
                elseif isprop(obj, name)
                    obj.(name) = varargin{i+1};
                else
                    warning('SRNNModel:UnknownProperty', 'Unknown property: %s', name);
                end
            end
            
            % Compute plot_deci from fs and plot_freq if not explicitly set
            if isempty(obj.plot_deci)
                obj.plot_deci = round(obj.fs / obj.plot_freq);
            end
        end
    end
    
    %% Dependent Property Getters
    methods
        function val = get.alpha(obj)
            val = obj.indegree / obj.n;
        end
        
        function val = get.default_val(obj)
            % DEFAULT_VAL Normalization factor F = 1/sqrt(N*alpha*(2-alpha))
            % Scaling factor which yields R=1 when all tilde parameters are equal.
            % See docs/EquationsParametersDocs/cell_type_pair_equations.md for the derivation (Harris 2023).
            %
            % F_tracks_network = true (default) computes F from the CURRENT n and
            % indegree. Because the n*alpha in R cancels against this F, the
            % theoretical spectral radius R is then exactly independent of n --
            % that is the point of the normalization.
            %
            % F_tracks_network = false pins F to the reference network
            % (F_ref_n, F_ref_indegree), so a sweep over n or indegree leaves the
            % weight distribution fixed and lets R vary with network size.
            %
            % NOTE: freezing F does NOT freeze the network. build() passes the
            % REAL obj.alpha to RMTMatrix, so sparsity and connectivity still
            % track the current n/indegree; only the scale of the weight
            % distribution is pinned.
            if obj.F_tracks_network
                n_F = obj.n;
                alpha_F = obj.alpha;
            else
                n_F = obj.F_ref_n;
                alpha_F = obj.F_ref_indegree / obj.F_ref_n;
            end
            val = 1 / sqrt(n_F * alpha_F * (2 - alpha_F));
        end

        function val = get.mu_E_tilde(obj)
            val = SRNNModel2.scale_tilde(obj.mu_E_tilde_relative, obj.default_val);
        end

        function val = get.mu_I_tilde(obj)
            val = SRNNModel2.scale_tilde(obj.mu_I_tilde_relative, obj.default_val);
        end

        function val = get.sigma_E_tilde(obj)
            val = SRNNModel2.scale_tilde(obj.sigma_E_tilde_relative, obj.default_val);
        end

        function val = get.sigma_I_tilde(obj)
            val = SRNNModel2.scale_tilde(obj.sigma_I_tilde_relative, obj.default_val);
        end

        function val = get.E_W(obj)
            val = SRNNModel2.scale_tilde(obj.E_W_relative, obj.default_val);
        end

        function val = f_vector(obj)
            % Population fractions as the 1 x 2 row [f_E, f_I] that the block
            % formulas expect (they scale COLUMNS, i.e. presynaptic types).
            val = [obj.f, 1 - obj.f];
        end

        function val = block_mu_s(obj)
            % Sparse-effective block means, Harris Eq. (15) blockwise.
            val = obj.alpha * obj.mu_tilde_block;
        end

        function val = block_sigma_s_sq(obj)
            % Sparse-effective block variances, Harris Eq. (16) blockwise. The
            % first term is why sparsity couples the MEAN into the bulk radius.
            mg = obj.mu_tilde_block;
            sg = obj.sigma_tilde_block;
            val = obj.alpha * (1 - obj.alpha) * mg.^2 + obj.alpha * sg.^2;
        end

        function val = get.mu_tilde_block(obj)
            % 2x2 absolute block means, (post, pre). A block override wins over
            % the column shorthand; E_W is added to every entry, as it was when
            % this was two scalars.
            rel = [SRNNModel2.pick(obj.mu_EE_relative, obj.mu_E_tilde_relative), ...
                   SRNNModel2.pick(obj.mu_EI_relative, obj.mu_I_tilde_relative); ...
                   SRNNModel2.pick(obj.mu_IE_relative, obj.mu_E_tilde_relative), ...
                   SRNNModel2.pick(obj.mu_II_relative, obj.mu_I_tilde_relative)];
            val = SRNNModel2.scale_tilde_mat(rel, obj.default_val) + obj.E_W;
        end

        function val = get.sigma_tilde_block(obj)
            rel = [SRNNModel2.pick(obj.sigma_EE_relative, obj.sigma_E_tilde_relative), ...
                   SRNNModel2.pick(obj.sigma_EI_relative, obj.sigma_I_tilde_relative); ...
                   SRNNModel2.pick(obj.sigma_IE_relative, obj.sigma_E_tilde_relative), ...
                   SRNNModel2.pick(obj.sigma_II_relative, obj.sigma_I_tilde_relative)];
            val = SRNNModel2.scale_tilde_mat(rel, obj.default_val);
        end

        function val = get.lambda_O(obj)
            % Outlier eigenvalues of E[W], descending by MAGNITUDE (an outlier is
            % defined by lying outside the bulk disk, and in an
            % inhibition-dominated network the dominant one is large and
            % negative). Generalizes Harris Eq. (17); see RMTBlocks.lambda_O.
            % Scaled by level_of_chaos because build_network multiplies W by it.
            K = obj.n * (obj.f_vector .* obj.block_mu_s);
            lam = eig(K) * obj.level_of_chaos;
            [~, order] = sort(abs(lam), 'descend');
            val = lam(order);
        end

        function val = get.activation_function(obj)
            val = obj.build_activation(1);
        end

        function val = get.activation_function_derivative(obj)
            val = obj.build_activation(2);
        end
        
        % The tilde properties are now always defined (they are computed from the
        % _relative multipliers), so these no longer need an isempty -> NaN guard.
        function val = get.mu_se(obj)
            val = obj.alpha * (obj.mu_E_tilde + obj.E_W);
        end

        function val = get.mu_si(obj)
            val = obj.alpha * (obj.mu_I_tilde + obj.E_W);
        end

        function val = get.sigma_se(obj)
            mu_eff = obj.mu_E_tilde + obj.E_W;
            val = sqrt(obj.alpha * (1 - obj.alpha) * mu_eff^2 + obj.alpha * obj.sigma_E_tilde^2);
        end

        function val = get.sigma_si(obj)
            mu_eff = obj.mu_I_tilde + obj.E_W;
            val = sqrt(obj.alpha * (1 - obj.alpha) * mu_eff^2 + obj.alpha * obj.sigma_I_tilde^2);
        end
        
        function val = get.R(obj)
            % Bulk spectral radius. Generalizes Harris Eq. (18) to block
            % statistics: the transition is at the Perron root of the block
            % variance matrix V(a,b) = n * f_b * sigma_s_sq(a,b), and R is its
            % square root (see RMTBlocks.R). Reduces exactly to the old
            % column-uniform formula, since every row of V is then identical and
            % the dominant eigenvalue is the trace.
            V = obj.n * (obj.f_vector .* obj.block_sigma_s_sq);
            val = sqrt(max(0, max(real(eig(V))))) * obj.level_of_chaos;
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
        
        function val = get.sigma_x_raw(obj)
            % The coefficient the integrator multiplies dW by. sigma_u_noise is
            % input-referred (units of u), and u enters dx/dt divided by tau_d.
            val = obj.sigma_u_noise / obj.tau_d;
        end

        function val = get.x_noise_std(obj)
            % Stationary std of x driven by noise alone: for the OU part
            % dx = -x/tau_d dt + sigma_x_raw dW, var = sigma_x_raw^2*tau_d/2.
            % NOMINAL only -- it ignores W and u, exactly as x0_std does -- but
            % it is the number to compare against S_c and the sigmoid width
            % when judging whether the noise is large enough to matter.
            val = obj.sigma_u_noise / sqrt(2 * obj.tau_d);
        end

        function val = get.N_sys_eqs(obj)
            nE = obj.n_E;
            nI = obj.n_I;
            val = nE * obj.n_a_E + nI * obj.n_a_I + nE * obj.n_b_E + nI * obj.n_b_I + obj.n;
        end
    end
    
    %% Public Methods
    methods
        function build(obj)
            % BUILD Initialize the network: create W, generate stimulus, cache params
            %
            % This method must be called before run(). It delegates to three
            % protected sub-methods that subclasses can override:
            %   1. build_network()   - Create W matrix
            %   2. build_stimulus()  - Generate external stimulus, interpolant, S0
            %   3. finalize_build()  - Validate and cache params
            
            obj.build_network();
            obj.build_stimulus();
            obj.finalize_build();
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
                jac_wrapper = @(t, S) SRNNModel2.compute_Jacobian_fast(S, params);
                obj.ode_opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'MaxStep', dt, 'Jacobian', jac_wrapper);
            end
            
            % Define RHS function using closure (avoids OOP overhead)
            % Cache interpolant and params to avoid property access on every call
            u_interp = obj.u_interpolant;
            params.u_interpolant = u_interp;  % Add to params for dynamics_fast
            rhs = @(t, S) SRNNModel2.dynamics_fast(t, S, params);

            % Pre-generate the Wiener increments (a no-op at sigma_u_noise = 0).
            % Must happen before both the trajectory and compute_lyapunov, since
            % Benettin re-integrates segments against the same increments.
            obj.build_noise();

            % Integrate
            fprintf('Integrating equations\n');
            tic
            [t_raw, S_raw] = obj.integrate(rhs, obj.t_ex, obj.S0);
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

            % The noise tensor has done its job (trajectory + Benettin) and is
            % the largest thing on the object; drop it. It is regenerable from
            % noise_seed, so nothing is lost.
            obj.noise_increments = [];

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
                if obj.T_range(2) >= 15
                    obj.lya_T_interval = [obj.T_range(1) + 15, obj.T_range(2)]; % skip first 15 seconds
                else
                    obj.lya_T_interval = [obj.T_range(1), obj.T_range(2)]; % fallback for short simulations
                end
            end
            
            % Define RHS function using closure (avoids OOP overhead)
            params.u_interpolant = obj.u_interpolant;
            rhs = @(t, S) SRNNModel2.dynamics_fast(t, S, params);
            
            fprintf('Computing Lyapunov exponents using %s method\n', obj.lya_method);
            % Benettin re-integrates trajectory segments, so it must use the
            % same integrator as the trajectory itself -- that is what makes the
            % discretisation error common to both trajectories and cancel in the
            % difference. (The QR path integrates a variational equation on a
            % 2-point span instead and always uses ode45; see below.)
            solver = resolve_solver(obj.ode_solver, obj.noise_increments, 'SRNNModel');
            obj.lya_results = SRNNModel2.compute_lyapunov_exponents_internal(obj.lya_method, obj.S_out, obj.t_out, dt, obj.fs, obj.lya_T_interval, obj.lya_warmup, obj.lya_dt, params, obj.ode_opts, solver, rhs);
            
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
            
            [fig_handle, ax_handles] = SRNNModel2.plot_SRNN_tseries(obj.plot_data.t, obj.plot_data.u, obj.plot_data.x, obj.plot_data.r, obj.plot_data.a, obj.plot_data.b, obj.plot_data.br, params, obj.lya_results, obj.lya_method, T_plot_arg);
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
            J_array = SRNNModel2.compute_Jacobian_at_indices(obj.S_out, J_times, params);
            
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
                ax_handles(i) = SRNNModel2.plot_eigenvalues_helper(evals, ax, time_val, global_xlim, global_ylim, [], 1/obj.tau_d);
                set(ax_handles(i), 'Color', 'none');
            end
            
            linkaxes(ax_handles, 'xy');
        end

        function [evals_all, evals_by_time] = eigenvalue_time_series(obj, J_times_sec, varargin)
            % EIGENVALUE_TIME_SERIES Jacobian eigenvalues sampled through a run
            %
            % Samples the instantaneous Jacobian at the requested times and
            % returns every eigenvalue seen. Because the network is nonlinear,
            % the Jacobian (and thus its eigenvalues) changes as the state
            % evolves; pooling the eigenvalues over many sampled times reveals
            % how much time they spend in each region of the complex plane.
            %
            % Requires store_full_state = true (so obj.S_out is available).
            %
            % The Jacobian build + eig at each sampled time are independent, so
            % they run under parfor by default (Parallel Computing Toolbox). The
            % state rows are sliced out first, so workers receive only the
            % n_times x N samples, not the full S_out trajectory.
            %
            % Usage:
            %   evals = model.eigenvalue_time_series(linspace(0, 20, 150));
            %   evals = model.eigenvalue_time_series(t, 'use_parallel', false);
            %
            % Name-value options:
            %   'use_parallel' - use parfor when a pool is available (default true)
            %
            % Outputs:
            %   evals_all     - column vector of ALL eigenvalues (all times pooled)
            %   evals_by_time - cell array, one entry per sampled time index

            p = inputParser;
            p.addParameter('use_parallel', true);
            p.parse(varargin{:});
            use_parallel = p.Results.use_parallel;

            if isempty(obj.S_out)
                error('SRNNModel:NoStateData', 'State data required. Set store_full_state=true.');
            end

            params = obj.cached_params;

            % Convert times to state-row indices (same mapping as plot_eigenvalues)
            J_times = round((J_times_sec - obj.t_out(1)) * obj.fs) + 1;
            J_times = unique(max(1, min(J_times, size(obj.S_out, 1))));

            % Slice the needed state rows up front so parfor ships only these
            % (n_times x N) rows to the workers, not the full S_out trajectory.
            S_sel = obj.S_out(J_times, :);
            n_times = size(S_sel, 1);
            evals_by_time = cell(n_times, 1);

            fprintf('Computing Jacobian eigenvalues at %d time points\n', n_times);

            % Fuse Jacobian build + eig per time point (avoids the large
            % N x N x n_times J_array). Independent across i -> parfor.
            run_parallel = use_parallel && canUseParallelPool;
            if use_parallel && ~canUseParallelPool
                warning('SRNNModel:NoParallelPool', ...
                    'Parallel pool not available. Falling back to sequential execution.');
            end

            if run_parallel
                parfor i = 1:n_times
                    J = full(SRNNModel2.compute_Jacobian_fast(S_sel(i, :)', params));
                    evals_by_time{i} = eig(J);
                end
            else
                for i = 1:n_times
                    J = full(SRNNModel2.compute_Jacobian_fast(S_sel(i, :)', params));
                    evals_by_time{i} = eig(J);
                end
            end
            evals_all = vertcat(evals_by_time{:});
        end

        function [fig_handle, ax_handle] = plot_eigenvalue_heatmap(obj, J_times_sec, varargin)
            % PLOT_EIGENVALUE_HEATMAP Occupancy heatmap of Jacobian eigenvalues
            %
            % Samples the Jacobian at J_times_sec (seconds), pools the
            % eigenvalues, and renders a Gaussian-smoothed 2-D density over the
            % complex plane. The stability line at Re = 0 is overlaid so the
            % time spent in the locally-unstable half-plane (Re > 0) is visible.
            %
            % Usage:
            %   model.plot_eigenvalue_heatmap(linspace(0, 20, 150));
            %   model.plot_eigenvalue_heatmap(t, 'grid_res', 300, 'use_log', false);
            %
            % Name-value options:
            %   'grid_res'   - number of bins per axis (default 200)
            %   'sigma_bins' - Gaussian smoothing width in bins (default 1.5)
            %   're_lim'     - [min max] real-axis limits (default from data)
            %   'im_lim'     - [min max] imag-axis limits (default from data)
            %   'use_log'    - log10(1+density) color scale (default true)
            %   'use_parallel' - parallelize the eig loop (default true)

            p = inputParser;
            p.addParameter('grid_res', 200);
            p.addParameter('sigma_bins', 1.5);
            p.addParameter('re_lim', []);
            p.addParameter('im_lim', []);
            p.addParameter('use_log', true);
            p.addParameter('use_parallel', true);
            p.parse(varargin{:});
            opt = p.Results;

            evals = obj.eigenvalue_time_series(J_times_sec, 'use_parallel', opt.use_parallel);

            % Axis limits: from data (with ~10% pad) unless overridden
            re = real(evals); im = imag(evals);
            re_lim = opt.re_lim; im_lim = opt.im_lim;
            if isempty(re_lim)
                re_lim = [min(re), max(re)] + [-0.1, 0.1] * max(range(re), eps);
            end
            if isempty(im_lim)
                im_lim = [min(im), max(im)] + [-0.1, 0.1] * max(range(im), eps);
            end

            re_edges = linspace(re_lim(1), re_lim(2), opt.grid_res + 1);
            im_edges = linspace(im_lim(1), im_lim(2), opt.grid_res + 1);

            D = SRNNModel2.compute_eigenvalue_density(evals, re_edges, im_edges, opt.sigma_bins);

            if opt.use_log
                clim = [0, max(log10(1 + D(:)))];
            else
                clim = [0, max(D(:))];
            end
            if clim(2) <= clim(1), clim(2) = clim(1) + 1; end

            fig_handle = figure('Position', [1312, 526, 560, 460]);
            ax_handle = axes(fig_handle);
            SRNNModel2.plot_eigenvalue_heatmap_helper(ax_handle, D, re_edges, im_edges, clim, opt.use_log);
            cb = colorbar(ax_handle);
            if opt.use_log
                cb.Label.String = 'log_{10}(1 + density)';
            else
                cb.Label.String = 'density';
            end
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
        
        function [fig_handle, ax_handle] = plot_W(obj, ax, clim_val)
            % PLOT_W Image the connectivity matrix with a diverging red/white/blue
            % colormap and a symmetric colour range, so zero is white and the
            % sign of every synapse is readable at a glance.
            %
            % Plots obj.W -- the SCALED matrix, after level_of_chaos and
            % rescale_by_abscissa -- i.e. the network actually simulated, and the
            % one R and lambda_O describe. RMTBlocks.plot_W shows the generator's
            % unscaled output instead, which differs by those factors.
            %
            % E/I boundaries are drawn over the image, so each quadrant is one
            % route: a raised mu_EE_relative reddens the top-left block only.
            %
            %   plot_W()             - new figure, clim from the off-diagonal max
            %   plot_W(ax)           - draw into ax
            %   plot_W(ax, clim_val) - draw into ax with a shared clim, for
            %                          comparing two networks on equal footing
            if ~obj.is_built
                error('SRNNModel:NotBuilt', 'Call build() first.');
            end
            if nargin < 2 || isempty(ax)
                fig_handle = figure('Name', 'SRNNModel2 W');
                ax = axes(fig_handle);
            else
                fig_handle = ancestor(ax, 'figure');
            end
            ax_handle = ax;

            W_plot = full(obj.W);

            % Auto range from the OFF-DIAGONAL entries, so a diagonal shift (if
            % any) saturates rather than flattening the rest of the picture.
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

            % E/I boundary and labels at the centre of each band.
            counts = [obj.n_E, obj.n_I];
            edge = obj.n_E + 0.5;
            hold(ax, 'on');
            plot(ax, [edge edge], [0.5, obj.n + 0.5], 'k-', 'LineWidth', 0.75);
            plot(ax, [0.5, obj.n + 0.5], [edge edge], 'k-', 'LineWidth', 0.75);
            hold(ax, 'off');
            centres = cumsum(counts) - counts / 2 + 0.5;
            set(ax, 'XTick', centres, 'XTickLabel', {'E', 'I'}, ...
                    'YTick', centres, 'YTickLabel', {'E', 'I'}, ...
                    'TickLength', [0 0]);
            xlabel(ax, 'presynaptic type');
            ylabel(ax, 'postsynaptic type');
            title(ax, sprintf('W  (level\\_of\\_chaos = %g,  R = %.3f)', ...
                obj.level_of_chaos, obj.R));
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
        
        function [cls, info] = checkConnectivityClass(obj, varargin)
            % CHECKCONNECTIVITYCLASS Digraph connectivity class of the network W
            %
            % Classifies the directed graph defined by the nonzero entries of W
            % as one of 'strong', 'unilateral', 'weak', or 'disconnected'. Call
            % after build() (or after build_network has assigned W).
            %
            % Usage:
            %   model.checkConnectivityClass()          % prints a summary
            %   cls = model.checkConnectivityClass()
            %   [cls, info] = model.checkConnectivityClass('Tolerance', 1e-12)
            %
            % Name-value options are passed through to the static
            % SRNNModel2.connectivity_class; see that method for the full list
            % of info fields.
            %
            % In addition to the fields set by connectivity_class, info carries
            % the *nominal* (configured) quantities, so the realized draw can be
            % compared against what the parameters predict:
            %   indegree_nominal, alpha_nominal, p_strong_est_nominal,
            %   expected_sources_sinks_nominal
            %
            % This is a purely structural check: it ignores weight magnitude and
            % sign. A strongly connected W whose only path between two halves
            % runs through a 1e-6 weight is functionally disconnected. Use
            % 'Tolerance' to screen out negligible weights.
            %
            % See also: connectivity_class, strong_connectivity_probability

            if isempty(obj.W)
                error('SRNNModel2:NotBuilt', ...
                    ['W is empty. Call build() (or build_network) before ' ...
                     'checkConnectivityClass. To classify an arbitrary matrix, ' ...
                     'use the static SRNNModel2.connectivity_class(A).']);
            end

            [cls, info] = SRNNModel2.connectivity_class(obj.W, varargin{:});

            % Nominal (configured) counterparts to the realized values above.
            info.indegree_nominal = obj.indegree;
            info.alpha_nominal    = obj.alpha;
            [info.p_strong_est_nominal, info.expected_sources_sinks_nominal] = ...
                SRNNModel2.strong_connectivity_probability(obj.n, obj.indegree);

            if nargout == 0
                fprintf(['Connectivity class: %s  (n = %d, indegree = %g, ' ...
                         'alpha = %.4f)\n'], cls, obj.n, obj.indegree, obj.alpha);
                fprintf(['  %d SCC(s), largest spans %d/%d nodes (%.1f%%); ' ...
                         '%d weak component(s)\n'], ...
                    info.n_scc, info.scc_sizes(1), obj.n, ...
                    100*info.largest_scc_frac, info.n_weak_comp);
                fprintf(['  %d node(s) with zero in-degree, %d with zero ' ...
                         'out-degree; realized mean in-degree = %.2f\n'], ...
                    info.n_sources, info.n_sinks, info.mean_indegree);
                fprintf(['  P(strongly connected) ~ %.3f for these settings ' ...
                         '(expected %.2f source/sink nodes)\n'], ...
                    info.p_strong_est_nominal, info.expected_sources_sinks_nominal);
                clear cls;   % suppress ans display
            end
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
            % The realized per-neuron setpoints, for diagnostics: the handles
            % above already carry them. Empty in the homogeneous case.
            params.S_c_vec = obj.S_c_vec;
            params.x0_std = obj.x0_std;
            
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
            dS_dt = SRNNModel2.dynamics_fast(t, S, params);
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
            
            % The RMT tildes need no default-filling here: they are Dependent on
            % the _relative multipliers and F, so they are always defined and are
            % recomputed whenever n or indegree changes. (The old version filled
            % them in place only when empty, which meant a rebuild after changing
            % n silently kept the F from the previous n.)

            % Validate the activation choice up front, so a typo fails here
            % rather than at the first phi() evaluation deep inside the solver.
            if ~isempty(obj.activation_custom)
                SRNNModel2.check_activation_custom(obj.activation_custom);
            elseif ~ismember(obj.activation, SRNNModel2.activation_names())
                error('SRNNModel:InvalidParams', ...
                    'Unknown activation ''%s''. Valid: %s.', ...
                    char(string(obj.activation)), ...
                    strjoin(SRNNModel2.activation_names(), ', '));
            end

            % Validate the reference network used when F_tracks_network is false.
            % This is a pairwise constraint, so it cannot live in a set method.
            if ~obj.F_tracks_network
                if ~isscalar(obj.F_ref_n) || ~isscalar(obj.F_ref_indegree) || ...
                        obj.F_ref_n <= 0 || obj.F_ref_indegree <= 0
                    error('SRNNModel:InvalidParams', ...
                        'F_ref_n and F_ref_indegree must be positive scalars.');
                end
                if obj.F_ref_indegree > obj.F_ref_n
                    error('SRNNModel:InvalidParams', ...
                        ['F_ref_indegree (%g) must not exceed F_ref_n (%g); the ' ...
                        'reference in-degree cannot be larger than the reference network.'], ...
                        obj.F_ref_indegree, obj.F_ref_n);
                end
            end

            % Compute tau_a arrays if n_a > 0 but tau_a not set
            if obj.n_a_E > 0 && isempty(obj.tau_a_E)
                obj.tau_a_E = logspace(log10(0.25), log10(10), obj.n_a_E);
            end
            if obj.n_a_I > 0 && isempty(obj.tau_a_I)
                obj.tau_a_I = logspace(log10(0.25), log10(10), obj.n_a_I);
            end
            
            % Create W matrix using RMTBlocks, which indexes the statistics by
            % BOTH postsynaptic (row) and presynaptic (column) type, so E->E can
            % differ from E->I. With no block overrides set the blocks are
            % column-uniform and this reproduces the old RMTMatrix result
            % exactly -- including bit-for-bit, since RMTBlocks consumes the RNG
            % in the same order (randn(N) then rand(N) via update_sparsity).
            %
            % g_mu / g_sigma are left at 1: level_of_chaos is still applied to
            % the assembled W below, because rescale_by_abscissa divides by the
            % abscissa BEFORE that multiply and folding the scale into the
            % generator would not be equivalent.
            rmt = RMTBlocks(obj.n);
            rmt.alpha = obj.alpha;
            rmt.f = obj.f;
            rmt.mu_tilde = obj.mu_tilde_block;
            rmt.sigma_tilde = obj.sigma_tilde_block;
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

            % Optional structural check. Relevant at low indegree, where the
            % network can fail to be strongly connected -- usually because a
            % neuron or two ends up a pure source or sink.
            if obj.check_connectivity
                [cls, cinfo] = obj.checkConnectivityClass();
                fprintf(['W connectivity: %s (%d SCC, largest %.0f%% of nodes; ' ...
                         '%d source / %d sink; P(strong) ~ %.2f for n = %d, ' ...
                         'indegree = %g)\n'], ...
                    cls, cinfo.n_scc, 100*cinfo.largest_scc_frac, ...
                    cinfo.n_sources, cinfo.n_sinks, cinfo.p_strong_est_nominal, ...
                    obj.n, obj.indegree);

                % SCC sizes, descending. Truncated so a large network with many
                % singleton components does not flood the command window.
                n_show  = min(20, numel(cinfo.scc_sizes));
                scc_str = mat2str(cinfo.scc_sizes(1:n_show));
                if numel(cinfo.scc_sizes) > n_show
                    scc_str = sprintf('%s ... (%d more)', scc_str, ...
                        numel(cinfo.scc_sizes) - n_show);
                end
                fprintf('  SCC sizes: %s\n', scc_str);
            end

            % Per-neuron nonlinearity setpoints. Drawn here, after W, so the
            % vector exists before build_stimulus() caches params.
            obj.build_setpoints();
        end

        function build_setpoints(obj)
            % BUILD_SETPOINTS Draw the per-neuron nonlinearity setpoints S_c_vec
            %
            % S_c_i = mu_S_c_<pop> + sigma_S_c_<pop> * randn, with an empty mu
            % falling back to the shared scalar S_c. Leaves S_c_vec EMPTY when
            % no heterogeneity was requested, which keeps every downstream code
            % path on the scalar branch and the results bit-identical to a
            % model without this feature.

            obj.S_c_vec = [];

            % Validated unconditionally, so a nonsensical value (a negative
            % standard deviation, say) is caught even though it would not by
            % itself switch heterogeneity on.
            SRNNModel2.check_setpoint_stat(obj.mu_S_c_E,    'mu_S_c_E',    true);
            SRNNModel2.check_setpoint_stat(obj.mu_S_c_I,    'mu_S_c_I',    true);
            SRNNModel2.check_setpoint_stat(obj.sigma_S_c_E, 'sigma_S_c_E', false);
            SRNNModel2.check_setpoint_stat(obj.sigma_S_c_I, 'sigma_S_c_I', false);

            has_mu    = ~isempty(obj.mu_S_c_E) || ~isempty(obj.mu_S_c_I);
            has_sigma = any([obj.sigma_S_c_E, obj.sigma_S_c_I] > 0);
            if ~(has_mu || has_sigma)
                return;
            end

            % Only 'logistic' and 'piecewise' take a centre. Failing loudly
            % here beats silently ignoring the request.
            if ~isempty(obj.activation_custom)
                error('SRNNModel:InvalidParams', ...
                    ['Per-neuron setpoints (mu_S_c_* / sigma_S_c_*) cannot be applied ' ...
                    'to activation_custom: a custom handle takes only x, so the model ' ...
                    'has no way to inject a per-neuron centre. Build the heterogeneity ' ...
                    'into the custom handle itself, or clear activation_custom.']);
            end
            if strcmp(obj.activation, 'tanh')
                error('SRNNModel:InvalidParams', ...
                    ['Per-neuron setpoints (mu_S_c_* / sigma_S_c_*) require a ' ...
                    'nonlinearity with a centre; ''tanh'' has none. Use ''logistic'' ' ...
                    'or ''piecewise''.']);
            end

            % Own RNG substream: the state is saved and restored so W, the
            % stimulus and x0 are bit-identical whether or not setpoints are
            % drawn (initialize_state draws x0 from whatever state build left
            % behind). The offset keeps the draw off the stream that produced
            % W's first column, which seeding with rng_seeds(1) directly would
            % not do.
            seed = obj.S_c_seed;
            if isempty(seed)
                seed = obj.rng_seeds(1) + 104729;
            end
            stream_state = rng;
            rng(seed);
            z = randn(obj.n, 1);
            rng(stream_state);

            mu_E    = SRNNModel2.pick(obj.mu_S_c_E, obj.S_c);
            mu_I    = SRNNModel2.pick(obj.mu_S_c_I, obj.S_c);
            sigma_E = obj.sigma_S_c_E;
            sigma_I = obj.sigma_S_c_I;

            vals = zeros(obj.n, 1);
            vals(obj.E_indices) = mu_E + sigma_E * z(obj.E_indices);
            vals(obj.I_indices) = mu_I + sigma_I * z(obj.I_indices);
            obj.S_c_vec = vals;

            fprintf(['Per-neuron S_c drawn (seed %g): E %.4f +/- %.4f, ' ...
                'I %.4f +/- %.4f\n'], seed, mu_E, sigma_E, mu_I, sigma_I);
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
        
        function finalize_build(obj)
            % FINALIZE_BUILD Validate configuration and cache params
            %
            % Called after build_network() and build_stimulus(). Validates
            % the configuration and caches the params struct for fast access
            % during ODE integration. Subclasses typically do not override.
            
            % Validate configuration
            obj.validate();
            
            % Cache params struct for fast access in run/plot methods
            obj.cached_params = obj.get_params();
            
            obj.is_built = true;
            fprintf('Model built successfully. Ready to run.\n');
        end
    end
    
    %% Private Methods
    methods (Access = protected)
        function [t_out, S_out] = integrate(obj, rhs, tspan, S0)
            % INTEGRATE Run the trajectory integrator named by ode_solver.
            %
            % The single place the main trajectory is stepped, so that
            % SRNN_ESN_reservoir -- which does not go through run(), it has its
            % own run_reservoir_esn -- gets any integrator work for free rather
            % than silently missing it.
            solver = resolve_solver(obj.ode_solver, obj.noise_increments, 'SRNNModel');
            [t_out, S_out] = solver(rhs, tspan, S0, obj.ode_opts);
        end

        function build_noise(obj)
            % BUILD_NOISE Pre-generate the Wiener increments for this run.
            %
            % Unit-variance normals, not scaled increments: sde_fixed_step
            % forms dW = sqrt(h)*xi at use time, which keeps the stored numbers
            % independent of the step size.
            obj.noise_increments = [];
            if obj.sigma_u_noise == 0
                return;     % no allocation, and bit-identical to the ODE path
            end

            seed = obj.noise_seed;
            if isempty(seed)
                % A different offset from the S_c draw's 104729, so the W,
                % setpoint and noise streams cannot overlap.
                seed = obj.rng_seeds(1) + 224737;
            end

            % Own RNG substream, saved and restored, so W, the stimulus, the
            % setpoints and x0 are bit-identical whether or not noise is drawn.
            n_steps = numel(obj.t_ex) - 1;
            stream_state = rng;
            rng(seed);
            xi1 = randn(obj.n, n_steps);
            xi2 = randn(obj.n, n_steps);
            rng(stream_state);

            % Noise reaches x only: the trailing n entries of the state vector.
            N = obj.N_sys_eqs;
            obj.noise_increments = struct( ...
                'xi1', xi1, 'xi2', xi2, ...
                't0', obj.t_ex(1), 'fs', obj.fs, ...
                'sigma', obj.sigma_x_raw, ...
                'idx', (N - obj.n + 1):N);
        end

        function h = build_activation(obj, which_one)
            % BUILD_ACTIVATION Handle for the chosen nonlinearity (1=fn, 2=derivative)
            %
            % Parameters are captured BY VALUE, not by capturing obj: the handle
            % is rebuilt on every property access so it can never be stale, and
            % capturing the model would drag the whole object into anything that
            % saves the handle (a saved obj-capturing handle costs ~3x as much,
            % and would embed an SRNNModel2 in every recorded parameter set).
            %
            % Each nonlinearity binds only the parameters it actually takes, so
            % S_a is genuinely unused by 'logistic' and both are unused by
            % 'tanh' -- which was already true, just implicit, when these were
            % hand-built handles.
            if ~isempty(obj.activation_custom)
                SRNNModel2.check_activation_custom(obj.activation_custom);
                h = obj.activation_custom{which_one};
                return;
            end

            a = obj.S_a;
            % A per-neuron setpoint vector, once build_setpoints() has drawn
            % one, otherwise the shared scalar. In the vector case the handle
            % is only valid for length-n input, since c lines up with x
            % element by element.
            c = SRNNModel2.pick(obj.S_c_vec, obj.S_c);
            switch obj.activation
                case 'logistic'
                    if which_one == 1
                        h = @(x) SRNNModel2.logisticSigmoid(x, c);
                    else
                        h = @(x) SRNNModel2.logisticSigmoidDerivative(x, c);
                    end
                case 'piecewise'
                    if which_one == 1
                        h = @(x) SRNNModel2.piecewiseSigmoid(x, a, c);
                    else
                        h = @(x) SRNNModel2.piecewiseSigmoidDerivative(x, a, c);
                    end
                case 'tanh'
                    if which_one == 1
                        h = @SRNNModel2.tanhActivation;
                    else
                        h = @SRNNModel2.tanhActivationDerivative;
                    end
                otherwise
                    error('SRNNModel:InvalidParams', ...
                        ['Unknown activation ''%s''. Valid: %s. (Set activation_custom ' ...
                        'for a nonlinearity outside this set.)'], ...
                        char(string(obj.activation)), ...
                        strjoin(SRNNModel2.activation_names(), ', '));
            end
        end

        function set_defaults(obj)
            % SET_DEFAULTS Initialize all properties to default values
            
            % The activation needs no setup here: `activation` defaults to
            % 'logistic' as a plain property, and the handles are Dependent.

            % Set default input configuration
            obj.input_config = struct();
            obj.input_config.n_steps = 3;
            obj.input_config.step_density = 0.2;
            obj.input_config.amp = 0.5;
            obj.input_config.no_stim_pattern = false(1, 3);
            obj.input_config.no_stim_pattern(1:2:end) = true;
            obj.input_config.intrinsic_drive = [];  % Will be set in build
            obj.input_config.positive_only = false;  % Default: allow positive and negative amplitudes
            obj.input_config.step_density_E = 0.15;   % Fraction of E neurons receiving input
            obj.input_config.step_density_I = 0;      % Fraction of I neurons receiving input (0 = E only)
            
            % T_plot defaults to T_range (set in build if not specified)
            obj.T_plot = [];
        end
        
        function validate(obj)
            % VALIDATE Check parameter consistency and constraints

            % Catch an ode_solver assigned after construction, so a typo fails
            % here rather than at the first solver call inside run().
            check_ode_solver(obj.ode_solver, 'SRNNModel');
            check_noise_settings(obj.sigma_u_noise, obj.ode_solver, 'SRNNModel');

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
        
        function [t_plot, S_plot, indices] = decimate_states(~, t_out, S_out, deci)
            % DECIMATE_STATES Decimates state trajectory for plotting
            %
            % Inputs:
            %   t_out - Time vector
            %   S_out - State matrix (nt x N)
            %   deci  - Decimation factor (integer)
            %
            % Outputs:
            %   t_plot  - Decimated time vector
            %   S_plot  - Decimated state matrix
            %   indices - Indices used for decimation
            
            indices = 1:deci:length(t_out);
            t_plot = t_out(indices);
            S_plot = S_out(indices, :);
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
    end

    methods (Static)
        function [x, a, b, r, br] = unpack_and_compute_states(S_out, params, a_zeros_b_ones)
            % UNPACK_AND_COMPUTE_STATES Unpack state vector and compute dependent variables
            %
            % Unpacks the state trajectory S_out into individual state variables,
            % splits them into excitatory and inhibitory components, and computes
            % the firing rate r and synaptic output br.
            %
            % Public static: this is a pure function of (S_out, params) and touches
            % no object state, so analysis scripts can call it directly rather than
            % re-deriving the state layout themselves. Matches the convention already
            % used by SRNNCellTypePairs. Because MATLAB does not pass
            % the object when a static method is dot-called on an instance, existing
            % obj.unpack_and_compute_states(S, params) call sites work unchanged.
            %
            % Inputs:
            %   S_out          - State trajectory (nt x N_sys_eqs)
            %   params         - Struct containing network parameters (see get_params)
            %   a_zeros_b_ones - (Optional, INTERNAL) If true, returns a as zeros and
            %                    b as ones, and returns x/r/br as plain n x nt arrays
            %                    instead of E/I structs. Used by the Jacobian path.

            % Handle optional parameter
            if nargin < 3
                a_zeros_b_ones = false;
            end

            % Guard against a params struct that does not describe this S_out.
            % Without this, a mismatched pair unpacks silently into wrong numbers.
            if size(S_out, 2) ~= params.N_sys_eqs
                error('SRNNModel2:unpack_and_compute_states:StateSizeMismatch', ...
                    ['S_out has %d columns but params.N_sys_eqs = %d. ' ...
                     'S_out must be nt x N_sys_eqs (time in rows).'], ...
                    size(S_out, 2), params.N_sys_eqs);
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
                x_eff_ts(params.E_indices, :) = x_eff_ts(params.E_indices, :) - ...
                    SRNNModel2.effective_c(params, 'c_E', 'n_a_E') * sum_a_E;
            end
            
            % Apply adaptation effect to I neurons (scaled by c_I)
            if params.n_I > 0 && params.n_a_I > 0 && ~isempty(a_I_ts)
                sum_a_I = squeeze(sum(a_I_ts, 2));  % n_I x nt
                if size(sum_a_I, 1) ~= params.n_I
                    sum_a_I = sum_a_I';
                end
                x_eff_ts(params.I_indices, :) = x_eff_ts(params.I_indices, :) - ...
                    SRNNModel2.effective_c(params, 'c_I', 'n_a_I') * sum_a_I;
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
        
        function ax = plot_eigenvalues_helper(eigenvalues, ax, time_value, x_lim, y_lim, circle_params, scalebar_length)
            % PLOT_EIGENVALUES_HELPER Plot eigenvalue distribution on complex plane
            % Internalized from src/plotting/plot_eigenvalues.m (Option B: 3-tier outlier coloring)
            
            if nargin < 7, scalebar_length = 1; end
            if nargin < 6, circle_params = []; end
            if nargin < 5, y_lim = []; end
            if nargin < 4, x_lim = []; end
            
            axes(ax);
            mSize = 4;
            hold on;
            
            has_circle = ~isempty(circle_params) && isfield(circle_params, 'center') && isfield(circle_params, 'radius');
            
            if has_circle
                R = circle_params.radius;
                xc = real(circle_params.center);
                yc = imag(circle_params.center);
                
                if isfield(circle_params, 'outlier_threshold')
                    outlier_threshold = circle_params.outlier_threshold;
                else
                    outlier_threshold = 1.04;
                end
                
                distances = abs(eigenvalues - circle_params.center);
                
                % Interior eigenvalues (within R): black unfilled circles
                interior_mask = distances <= R;
                interior_eigs = eigenvalues(interior_mask);
                plot(real(interior_eigs), imag(interior_eigs), 'ko', 'MarkerSize', mSize, 'MarkerFaceColor', 'none', 'LineWidth', 0.5);
                
                % Near outlier eigenvalues (between R and outlier_threshold*R): black Xs
                near_outlier_mask = (distances > R) & (distances <= outlier_threshold * R);
                near_outlier_eigs = eigenvalues(near_outlier_mask);
                if ~isempty(near_outlier_eigs)
                    plot(real(near_outlier_eigs), imag(near_outlier_eigs), 'kx', 'MarkerSize', mSize, 'LineWidth', 0.5);
                end
                
                % Far outlier eigenvalues (beyond outlier_threshold*R): green filled circles
                far_outlier_mask = distances > outlier_threshold * R;
                far_outlier_eigs = eigenvalues(far_outlier_mask);
                if ~isempty(far_outlier_eigs)
                    plot(real(far_outlier_eigs), imag(far_outlier_eigs), 'o', 'MarkerSize', mSize, 'MarkerFaceColor', [0 .7 0], 'MarkerEdgeColor', [0 .7 0]);
                end
                
                % Draw theoretical radius as solid black circle
                theta = linspace(0, 2*pi, 100);
                plot(xc + R*cos(theta), yc + R*sin(theta), 'k-', 'LineWidth', 2);
            else
                % No circle params: plot all eigenvalues as black unfilled circles
                plot(real(eigenvalues), imag(eigenvalues), 'ko', 'MarkerSize', mSize, 'MarkerFaceColor', 'none', 'LineWidth', 0.5);
            end
            
            if isempty(x_lim), x_lim = xlim; end
            if x_lim(2) < 0, x_lim(2) = 0.05; end
            if isempty(y_lim), y_lim = ylim; end
            
            axis off;
            hold on;
            h_x = plot(x_lim, [0, 0], 'k', 'LineWidth', 1.25);
            h_y = plot([0, 0], y_lim, 'k', 'LineWidth', 1.25);
            uistack([h_x, h_y], 'bottom');
            
            text(1.02*x_lim(2), 0, ' Re($\lambda$)', 'Interpreter', 'latex', 'VerticalAlignment', 'middle');
            text(0, y_lim(2), 'Im($\lambda$)', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
            
            % Scale bar in the lower left
            sb_x = x_lim(1) + 0.05 * diff(x_lim);
            sb_y = y_lim(1) + 0.05 * diff(y_lim);
            plot([sb_x, sb_x + scalebar_length], [sb_y, sb_y], 'k-', 'LineWidth', 2);
            text(sb_x + scalebar_length/2, sb_y, num2str(scalebar_length, '%g'), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            
            xlim(x_lim);
            ylim(y_lim);
            hold off;
            axis equal;
        end

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
            %   dx_i/dt = (-x_i + sum_j(w_ij * b_j * r_j) + u_i) / tau_d
            %   r_i = phi(x_i - c * sum_k(a_i,k))
            %   da_i,k/dt = (-a_i,k + r_i) / tau_k
            %   db_i/dt = (1 - b_i) / tau_rec - (b_i * r_i) / tau_rel
            %
            % Note the placement of b: r_i is the PRE-depression rate, and
            % depression enters as the product b_j*r_j in the recurrent sum
            % (presynaptic, multiplicative). Both SFA and STD are therefore
            % driven by the raw rate r_i, not by b_i*r_i.
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
            
            c_E = SRNNModel2.effective_c(params, 'c_E', 'n_a_E');
            c_I = SRNNModel2.effective_c(params, 'c_I', 'n_a_I');
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
        
        function lya_results = compute_lyapunov_exponents_internal(Lya_method, S_out, t_out, dt, fs, T_interval, lya_warmup, lya_dt, params, opts, ode_solver, rhs_func)
            % Compute Lyapunov exponents using Benettin or QR method.
            % Internalized from ConnectivityAdaptation/src/algorithms/Lyapunov/compute_lyapunov_exponents.m

            lya_results = struct();

            if strcmpi(Lya_method, 'none')
                return;
            end

            lya_dt = SRNNModel2.resolve_lya_dt(lya_dt, Lya_method, dt, 'SRNNModel');

            lya_fs = 1 / lya_dt;
            
            switch lower(Lya_method)
                case 'benettin'
                    fprintf('Computing largest Lyapunov exponent using Benettin''s algorithm...\n');
                    d0 = 1e-3;
                    tic
                    [LLE, local_lya, finite_lya, t_lya] = SRNNModel2.benettin_algorithm_internal(S_out, t_out, dt, fs, d0, T_interval, lya_dt, lya_warmup, opts, rhs_func, ode_solver);
                    toc
                    lya_results.LLE = LLE;
                    lya_results.local_lya = local_lya;
                    lya_results.finite_lya = finite_lya;
                    lya_results.t_lya = t_lya;
                    lya_results.lya_dt = lya_dt;
                    lya_results.lya_fs = lya_fs;
                    
                case 'qr'
                    fprintf('Computing full Lyapunov spectrum using QR decomposition method...\n');
                    tic
                    jacobian_wrapper = @(tt, S, p) SRNNModel2.compute_Jacobian_fast(S, p);
                    [LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya] = SRNNModel2.lyapunov_spectrum_qr_internal(S_out, t_out, lya_dt, params, ode_solver, opts, jacobian_wrapper, T_interval, lya_warmup, params.N_sys_eqs, fs);
                    toc
                    fprintf('Lyapunov Dimension: %.2f\n', SRNNModel2.compute_kaplan_yorke_dimension_internal(LE_spectrum));
                    lya_results.LE_spectrum = LE_spectrum;
                    lya_results.local_LE_spectrum_t = local_LE_spectrum_t;
                    lya_results.finite_LE_spectrum_t = finite_LE_spectrum_t;
                    lya_results.t_lya = t_lya;
                    lya_results.params.N_sys_eqs = params.N_sys_eqs;
                    
                    [sorted_LE, sort_idx] = sort(real(lya_results.LE_spectrum), 'descend');
                    lya_results.LE_spectrum = sorted_LE;
                    lya_results.local_LE_spectrum_t = lya_results.local_LE_spectrum_t(:, sort_idx);
                    lya_results.finite_LE_spectrum_t = lya_results.finite_LE_spectrum_t(:, sort_idx);
                    lya_results.sort_idx = sort_idx;
                    lya_results.lya_dt = lya_dt;
                    lya_results.lya_fs = lya_fs;
                    fprintf('Largest Lyapunov Exponent (sorted): %.4f\n', lya_results.LE_spectrum(1));
                    
                otherwise
                    error('Unknown Lyapunov method: %s', Lya_method);
            end
        end
        
        function [LLE, local_lya, finite_lya, t_lya] = benettin_algorithm_internal(X, t, dt, fs, d0, T, lya_dt, lya_warmup, ode_options, dynamics_func, ode_solver)
            % Benettin's algorithm to compute the largest Lyapunov exponent.
            % Internalized from ConnectivityAdaptation/src/algorithms/Lyapunov/benettin_algorithm.m
            %
            % T is lya_T_interval. Accumulation runs over [max(0,T(1)), T(2)];
            % iteration starts lya_warmup seconds earlier so the perturbation
            % can align with the leading direction first.

            if ~isscalar(lya_dt) || ~isnumeric(lya_dt) || lya_dt <= 0
                error('lya_dt must be a positive scalar.');
            end

            deci_lya = round(lya_dt * fs);
            if deci_lya < 1
                error('lya_dt * fs must result in at least 1 sample per Lyapunov interval.');
            end

            tau_lya = dt * deci_lya;
            tol = dt * 1e-6;

            [idx_all, t_lya, acc_start] = SRNNModel2.lyapunov_sample_grid( ...
                t, dt, deci_lya, tau_lya, T, lya_warmup, 'SRNNModel');
            nt_lya = numel(t_lya);

            local_lya = zeros(nt_lya, 1);
            finite_lya = nan(nt_lya, 1);
            sum_log_stretching_factors = 0;

            n_state = size(X, 2);
            rnd_IC = randn(n_state, 1);
            pert = (rnd_IC ./ norm(rnd_IC)) .* d0;

            for k = 1:nt_lya
                idx_start = idx_all(k);
                idx_end = idx_start + deci_lya;

                X_k_pert = X(idx_start, :).' + pert;

                t_seg_detailed = t_lya(k) + (0:dt:tau_lya);

                [~, X_pert_output_all_steps] = ode_solver(dynamics_func, t_seg_detailed, X_k_pert, ode_options);

                X_pert_end = X_pert_output_all_steps(end, :).';
                X_end = X(idx_end, :).';

                delta = X_pert_end - X_end;
                d_k = norm(delta);
                local_lya(k) = log(d_k / d0) / tau_lya;

                if ~isfinite(local_lya(k))
                    warning('System diverged at t=%f. Truncating results.', t_lya(k));
                    if k > 1
                        last_valid = finite_lya(1:k-1);
                        last_valid = last_valid(~isnan(last_valid));
                        if ~isempty(last_valid)
                            LLE = last_valid(end);
                        else
                            LLE = 0;
                        end
                    else
                        LLE = 0;
                    end
                    local_lya(k:end) = [];
                    finite_lya(k:end) = [];
                    t_lya(k:end) = [];
                    return;
                end

                pert = (delta ./ d_k) .* d0;

                % Accumulate only inside lya_T_interval; the warmup samples
                % before acc_start align the perturbation but do not count.
                if t_lya(k) >= acc_start - tol && t_lya(k) < T(2)
                    sum_log_stretching_factors = sum_log_stretching_factors + log(d_k / d0);
                    elapsed = t_lya(k) + tau_lya - acc_start;
                    finite_lya(k, 1) = sum_log_stretching_factors / max(elapsed, eps);
                end
            end

            last_valid = finite_lya(~isnan(finite_lya));
            if ~isempty(last_valid)
                LLE = last_valid(end);
            else
                % Segments existed but none fell inside the window. Reporting a
                % bare 0 here would read as "edge of chaos" rather than "not
                % measured", so say what happened.
                warning('SRNNModel:LyapunovWindowEmpty', ...
                    ['No Lyapunov segment fell inside lya_T_interval = [%g, %g] ' ...
                     '(segment length %g s); LLE is reported as 0. Widen the ' ...
                     'window or reduce lya_dt.'], T(1), T(2), tau_lya);
                LLE = 0;
            end
        end

        function [idx_all, t_lya, acc_start] = lyapunov_sample_grid(t, dt, deci_lya, tau_lya, T, lya_warmup, err_id_prefix)
            % Shared Benettin/QR sample grid: which trajectory rows start a
            % Lyapunov segment, and where accumulation begins.
            %
            % Accumulation begins at acc_start = max(0, T(1)) -- lya_T_interval
            % is honoured, and a negative T_range pre-window is never counted.
            % Iteration begins lya_warmup seconds earlier, CLAMPED to t(1): a
            % warmup that would reach before the simulation start is shortened
            % (with a warning) rather than dropped, because those samples are
            % free alignment. Segments that would run past the data or past
            % T(2) are trimmed.

            tol = dt * 1e-6;
            acc_start = max(0, T(1));

            requested_start = acc_start - lya_warmup;
            if requested_start < t(1) - tol
                available = max(0, acc_start - t(1));
                if lya_warmup > 0
                    warning([err_id_prefix ':LyapunovWarmupClamped'], ...
                        ['lya_warmup = %g s would start the Lyapunov iteration at ' ...
                         't = %g s, before the simulation start t = %g s. Using %g s ' ...
                         'of warmup instead.'], ...
                        lya_warmup, requested_start, t(1), available);
                end
                i_start = 1;
            else
                i_start = round((requested_start - t(1)) / dt) + 1;
                i_start = max(1, min(i_start, numel(t)));
            end

            idx_all = i_start:deci_lya:numel(t);
            t_lya = t(idx_all);
            t_lya = t_lya(:);

            % A segment needs deci_lya further samples, and must not extend
            % past the end of the requested Lyapunov window.
            last_allowed = min(t(end), T(2));
            keep = (idx_all + deci_lya <= numel(t)) & ...
                (t_lya(:).' + tau_lya <= last_allowed + tol);
            idx_all = idx_all(keep);
            t_lya = t_lya(keep);

            if isempty(t_lya)
                error([err_id_prefix ':NoLyapunovIntervals'], ...
                    ['No Lyapunov intervals fit in lya_T_interval = [%g, %g] with ' ...
                     'lya_dt giving tau = %g s over a trajectory spanning [%g, %g].'], ...
                    T(1), T(2), tau_lya, t(1), t(end));
            end
        end
        
        function lya_dt = resolve_lya_dt(lya_dt, method, dt, err_id_prefix)
            % Resolve the lya_dt property: empty means the per-method default.
            % The two differ because QR renormalises an N x N basis per step
            % while Benettin renormalises one vector, so QR can afford far
            % fewer, larger steps.
            if isempty(lya_dt)
                if strcmpi(method, 'benettin')
                    lya_dt = 0.02;
                else
                    lya_dt = 0.1;
                end
            end

            if ~isscalar(lya_dt) || ~isnumeric(lya_dt) || ~isfinite(lya_dt) || lya_dt <= 0
                error([err_id_prefix ':InvalidLyapunovStep'], ...
                    'lya_dt must be a positive finite scalar (or empty for the %s default).', ...
                    lower(method));
            end

            factor = lya_dt / dt;
            if abs(round(factor) - factor) > 1e-11
                error([err_id_prefix ':InvalidLyapunovStep'], ...
                    ['lya_dt must be an integer multiple of dt = 1/fs. ' ...
                     'lya_dt = %g and dt = %g give a ratio of %g.'], lya_dt, dt, factor);
            end
            if factor < 3
                error([err_id_prefix ':InvalidLyapunovStep'], ...
                    ['lya_dt must be at least 3*dt. lya_dt = %g and dt = %g give a ' ...
                     'ratio of %g; increase fs or increase lya_dt.'], lya_dt, dt, factor);
            end
        end

        function [LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya_vec] = lyapunov_spectrum_qr_internal(X_fid_traj, t_fid_traj, lya_dt_interval, params, ~, ode_options_main, jacobian_func_handle, T_full_interval, lya_warmup, N_states_sys, fs_fid)
            % QR method for full Lyapunov spectrum.
            % Internalized from ConnectivityAdaptation/src/algorithms/Lyapunov/lyapunov_spectrum_qr.m
            %
            % Honours T_full_interval (= lya_T_interval) the same way Benettin
            % does: accumulate over [max(0,T(1)), T(2)], iterate from lya_warmup
            % seconds earlier so Q can align with the leading subspace first.

            fiducial_interpolants = cell(N_states_sys, 1);
            for i = 1:N_states_sys
                fiducial_interpolants{i} = griddedInterpolant(t_fid_traj, X_fid_traj(:, i), 'pchip');
            end
            
            dt_fid = 1 / fs_fid;
            
            fid_dt = diff(t_fid_traj);
            if ~isempty(fid_dt)
                nominal_dt = median(fid_dt);
                max_dt_dev = max(abs(fid_dt - nominal_dt));
                tol_dt_dev = max(1e-4 * max(nominal_dt, eps(nominal_dt)), eps(dt_fid));
                if max_dt_dev > tol_dt_dev
                    error('Fiducial trajectory timestamps are not uniformly spaced. Max deviation %.3g s exceeds tolerance %.3g s.', max_dt_dev, tol_dt_dev);
                end
            end
            
            deci_lya = round(lya_dt_interval / dt_fid);
            if deci_lya == 0
                error('lya_dt_interval is too small.');
            end
            tau_lya = dt_fid * deci_lya;
            
            [~, t_lya_vec, acc_start] = SRNNModel2.lyapunov_sample_grid( ...
                t_fid_traj, dt_fid, deci_lya, tau_lya, T_full_interval, lya_warmup, 'SRNNModel');
            tol_acc = dt_fid * 1e-6;
            nt_lya = numel(t_lya_vec);

            Q_current = eye(N_states_sys);
            sum_log_R_diag = zeros(N_states_sys, 1);
            
            local_LE_spectrum_t = zeros(nt_lya, N_states_sys);
            finite_LE_spectrum_t = nan(nt_lya, N_states_sys);
            
            total_positive_time_accumulated = 0;
            ode_options_var = ode_options_main;
            
            for k = 1:nt_lya
                t_start_segment = t_lya_vec(k);
                t_end_segment = t_start_segment + tau_lya;
                t_end_segment = min(t_end_segment, t_fid_traj(end));
                
                current_segment_duration = t_end_segment - t_start_segment;
                if current_segment_duration <= eps
                    if k > 1
                        local_LE_spectrum_t(k, :) = local_LE_spectrum_t(k-1, :);
                        finite_LE_spectrum_t(k, :) = finite_LE_spectrum_t(k-1, :);
                    else
                        local_LE_spectrum_t(k, :) = NaN;
                        finite_LE_spectrum_t(k, :) = NaN;
                    end
                    continue;
                end
                
                t_span_ode = [t_start_segment, t_end_segment];
                Psi0_vec = reshape(Q_current, [], 1);
                
                variational_eqs = @(tt, current_Psi_vec) SRNNModel2.variational_eqs_ode_internal(tt, current_Psi_vec, fiducial_interpolants, N_states_sys, jacobian_func_handle, params);
                % ode45 rather than the model's ode_solver: this is a 2-POINT
                % span, which the fixed-step solvers reject outright (see
                % ode_rk4.m), so QR has only ever worked with an adaptive
                % solver. The variational equation is deterministic and
                % independent of how the fiducial trajectory was produced, so
                % pinning it here is not a loss of generality.
                [~, Psi_solution_vec] = ode45(variational_eqs, t_span_ode, Psi0_vec, ode_options_var);
                
                Psi_evolved_matrix = reshape(Psi_solution_vec(end, :)', [N_states_sys, N_states_sys]);
                
                if any(~isfinite(Psi_evolved_matrix(:)))
                    warning('System diverged at t=%f. Truncating results.', t_start_segment);
                    if total_positive_time_accumulated > eps
                        LE_spectrum = sum_log_R_diag / total_positive_time_accumulated;
                    else
                        LE_spectrum = nan(N_states_sys, 1);
                    end
                    t_lya_vec(k:end) = [];
                    local_LE_spectrum_t(k:end, :) = [];
                    finite_LE_spectrum_t(k:end, :) = [];
                    return;
                end
                
                [Q_new, R_segment] = qr(Psi_evolved_matrix);
                diag_R = diag(R_segment);
                
                log_abs_diag_R = log(abs(diag_R));
                valid_diag_R = abs(diag_R) > eps;
                
                current_local_LEs = zeros(N_states_sys, 1);
                current_local_LEs(valid_diag_R) = log_abs_diag_R(valid_diag_R) / current_segment_duration;
                current_local_LEs(~valid_diag_R) = -Inf;
                local_LE_spectrum_t(k, :) = current_local_LEs';
                
                % Accumulate only inside lya_T_interval; the warmup segments
                % before acc_start align Q but do not count.
                in_window = t_start_segment >= acc_start - tol_acc && ...
                    t_start_segment < T_full_interval(2);
                if in_window
                    sum_log_R_diag(valid_diag_R) = sum_log_R_diag(valid_diag_R) + log_abs_diag_R(valid_diag_R);
                    total_positive_time_accumulated = total_positive_time_accumulated + current_segment_duration;
                end

                if total_positive_time_accumulated > eps
                    finite_LE_spectrum_t(k, :) = (sum_log_R_diag / total_positive_time_accumulated)';
                elseif k > 1 && in_window
                    finite_LE_spectrum_t(k, :) = finite_LE_spectrum_t(k-1, :);
                else
                    finite_LE_spectrum_t(k, :) = NaN;
                end
                
                Q_current = Q_new;
            end
            
            if total_positive_time_accumulated > eps
                LE_spectrum = sum_log_R_diag / total_positive_time_accumulated;
            else
                warning('No accumulation over positive time for global LEs.');
                LE_spectrum = nan(N_states_sys, 1);
            end
        end
        
        function dPsi_vec_dt = variational_eqs_ode_internal(tt, current_Psi_vec, fiducial_interpolants, N_states_sys, jacobian_func_handle, params)
            % Variational ODE for QR method.
            X_fid_at_tt = zeros(N_states_sys, 1);
            for state_idx = 1:N_states_sys
                X_fid_at_tt(state_idx) = fiducial_interpolants{state_idx}(tt);
            end
            J_matrix = jacobian_func_handle(tt, X_fid_at_tt, params);
            Psi_matrix = reshape(current_Psi_vec, [N_states_sys, N_states_sys]);
            dPsi_matrix_dt = J_matrix * Psi_matrix;
            dPsi_vec_dt = reshape(dPsi_matrix_dt, [], 1);
        end
        
        function D_KY = compute_kaplan_yorke_dimension_internal(lambda)
            % Compute Kaplan-Yorke (Lyapunov) dimension from spectrum.
            lambda = sort(lambda, 'descend');
            cumsum_lambda = cumsum(lambda);
            j = find(cumsum_lambda >= 0, 1, 'last');
            
            if isempty(j)
                D_KY = 0;
            elseif j == length(lambda)
                D_KY = length(lambda);
            else
                D_KY = j + cumsum_lambda(j) / abs(lambda(j + 1));
            end
        end
    end
    
    %% ====================================================================
    %              CONNECTIVITY STRUCTURE (graph diagnostics)
    % =====================================================================
    % Structural checks on W. Pure functions of a matrix; the instance-level
    % entry point is checkConnectivityClass.

    methods (Static)
        function [cls, info] = connectivity_class(A, varargin)
            % CONNECTIVITY_CLASS Classify the digraph defined by a matrix
            %
            % [cls, info] = SRNNModel2.connectivity_class(A, ...)
            %
            % Returns one of, in order of decreasing strength:
            %   'strong'       - every node reaches every other node
            %   'unilateral'   - for every pair (i,j) there is a path i->j OR j->i
            %   'weak'         - the underlying undirected graph is connected
            %   'disconnected' - none of the above
            % The classes nest (strong subset unilateral subset weak), so the
            % strongest applicable label is returned.
            %
            % Name-value options:
            %   'Tolerance' (0)    - entries with abs(A) <= Tolerance are treated
            %                        as absent edges. Useful for screening out
            %                        weights too small to matter dynamically.
            %   'Transpose' (true) - interpret A(i,j) as an edge j->i, which is
            %                        this model's convention (dx_i/dt includes
            %                        sum_j w_ij r_j). MATLAB's digraph(A) reads
            %                        A(i,j) as i->j, so the adjacency is
            %                        transposed before the graph is built. All
            %                        four classes are invariant under edge
            %                        reversal, so cls is unaffected, but the
            %                        per-node source/sink diagnostics are not.
            %
            % Self-loops are removed before analysis: they never affect the
            % connectivity class, and RMTMatrix's sparsity mask includes the
            % diagonal, so they would otherwise inflate the degree counts.
            %
            % info fields (realized structure):
            %   class, is_strong, is_unilateral, is_weakly_connected
            %   n_scc, scc_sizes (descending), largest_scc_frac, scc_bins
            %   n_weak_comp, weak_bins
            %   in_degree, out_degree      (1 x n, self-loops excluded)
            %   zero_indegree, zero_outdegree (index vectors)
            %   n_sources, n_sinks, n_selfloops, n_edges
            %   mean_indegree, alpha_hat   (realized density over n*(n-1) slots)
            %   p_strong_est, expected_sources_sinks
            %       Erdos-Renyi estimate evaluated at the REALIZED density; see
            %       strong_connectivity_probability.
            %
            % largest_scc_frac is usually the field of interest: it separates
            % "one big SCC plus a couple of dangling source/sink neurons" from a
            % genuinely fractured reservoir.
            %
            % See also: checkConnectivityClass, strong_connectivity_probability

            p = inputParser;
            p.FunctionName = 'SRNNModel2.connectivity_class';
            addRequired(p, 'A', @(x) isnumeric(x) || islogical(x));
            addParameter(p, 'Tolerance', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
            addParameter(p, 'Transpose', true, @(x) islogical(x) || isnumeric(x));
            parse(p, A, varargin{:});
            tol       = p.Results.Tolerance;
            do_transp = logical(p.Results.Transpose);

            if ndims(A) ~= 2 || size(A, 1) ~= size(A, 2) %#ok<ISMAT>
                error('SRNNModel2:connectivity_class:NotSquare', ...
                    'A must be a square matrix; got %s.', mat2str(size(A)));
            end
            n = size(A, 1);
            if n == 0
                error('SRNNModel2:connectivity_class:Empty', 'A must be non-empty.');
            end

            % Structural adjacency, oriented so that Adj(i,j) means i -> j.
            Adj = abs(full(A)) > tol;
            if do_transp
                Adj = Adj.';
            end
            n_selfloops = nnz(diag(Adj));
            Adj(1:n+1:end) = false;      % self-loops do not affect the class

            G = digraph(Adj);

            scc_bins  = conncomp(G, 'Type', 'strong');
            weak_bins = conncomp(G, 'Type', 'weak');
            n_scc       = max(scc_bins);
            n_weak_comp = max(weak_bins);

            is_strong = (n_scc == 1);
            is_weak   = (n_weak_comp == 1);

            % Unilateral <=> the condensation DAG has a Hamiltonian path
            % <=> its topological order is unique <=> consecutive nodes in any
            % topological order are joined by an edge. O(V+E). A disconnected
            % condensation fails this automatically, so no special-casing.
            if is_strong
                is_unilateral = true;
            else
                Gc  = condensation(G);
                ord = toposort(Gc);
                is_unilateral = true;
                for k = 1:(numnodes(Gc) - 1)
                    if findedge(Gc, ord(k), ord(k+1)) == 0
                        is_unilateral = false;
                        break;
                    end
                end
            end

            if is_strong
                cls = 'strong';
            elseif is_unilateral
                cls = 'unilateral';
            elseif is_weak
                cls = 'weak';
            else
                cls = 'disconnected';
            end

            if nargout < 2
                return;
            end

            in_deg  = full(sum(Adj, 1));    % Adj(i,j) = i->j, so column sums are in-degrees
            out_deg = full(sum(Adj, 2)).';

            scc_sizes = sort(accumarray(scc_bins(:), 1), 'descend').';

            info = struct();
            info.class               = cls;
            info.is_strong           = is_strong;
            info.is_unilateral       = is_unilateral;
            info.is_weakly_connected = is_weak;

            info.n_scc            = n_scc;
            info.scc_sizes        = scc_sizes;
            info.largest_scc_frac = scc_sizes(1) / n;
            info.scc_bins         = scc_bins;
            info.n_weak_comp      = n_weak_comp;
            info.weak_bins        = weak_bins;

            info.in_degree      = in_deg;
            info.out_degree     = out_deg;
            info.zero_indegree  = find(in_deg  == 0);
            info.zero_outdegree = find(out_deg == 0);
            info.n_sources      = numel(info.zero_indegree);
            info.n_sinks        = numel(info.zero_outdegree);
            info.n_selfloops    = n_selfloops;
            info.n_edges        = nnz(Adj);

            info.mean_indegree = info.n_edges / n;
            if n > 1
                info.alpha_hat = info.n_edges / (n * (n - 1));
            else
                info.alpha_hat = 0;
            end

            [info.p_strong_est, info.expected_sources_sinks] = ...
                SRNNModel2.strong_connectivity_probability(n, info.alpha_hat * n);
        end

        function [p_strong, expected_sources_sinks] = strong_connectivity_probability(n, indegree)
            % STRONG_CONNECTIVITY_PROBABILITY P(strongly connected) for D(n,p)
            %
            % [p, lambda] = SRNNModel2.strong_connectivity_probability(n, indegree)
            %
            % Estimates the probability that an Erdos-Renyi directed random
            % graph on n nodes with expected in-degree `indegree` (edge
            % probability alpha = indegree/n, as built by RMTMatrix) is strongly
            % connected, and returns lambda = the expected number of nodes with
            % zero in-degree or zero out-degree.
            %
            % At the densities used here the dominant obstruction to strong
            % connectivity is a single node of in- or out-degree zero, not a
            % large structural split (Palasti: P(strong) -> exp(-2c) when
            % n*(1-p)^(n-1) -> c). This estimate accounts only for that
            % obstruction and treats the 2n zero-degree events as independent,
            % so it is an APPROXIMATION, not a bound: ignoring larger splits
            % biases it high, while the independence assumption (those events
            % are in fact negatively correlated) biases it low. Empirically the
            % two roughly cancel -- at n = 50, indegree = 4 it returns 0.184
            % against a Monte-Carlo value of ~0.21 (see
            % scripts/tests/test_connectivity_class.m). Treat it as good to a
            % few percent in the alpha ~ log(n)/n regime, and as a rough guide
            % below it.
            %
            % Example (the fig_stim_engages_adaptation bursting network):
            %   [p, lam] = SRNNModel2.strong_connectivity_probability(50, 4)
            %   -> p ~ 0.18, lam ~ 1.68
            % i.e. only about one draw in five is strongly connected, and one to
            % two neurons are expected to be pure sources or sinks. Mean degree
            % log(n) ~ 3.9 is the connectivity threshold, so indegree = 4 sits
            % right on it.
            %
            % See also: connectivity_class, checkConnectivityClass

            alpha = indegree / n;
            alpha = min(max(alpha, 0), 1);

            % P(a given node has in-degree 0), excluding the self-loop slot.
            p_zero = (1 - alpha)^(n - 1);

            expected_sources_sinks = 2 * n * p_zero;
            p_strong = (1 - p_zero)^(2 * n);
        end
    end

    %% ====================================================================
    %              INTERNALIZED ACTIVATION FUNCTIONS
    % =====================================================================
    % Internalized from the standalone nonlinearity files  to make SRNNModel2 standalone.

    methods (Static)
        % compute_eigenvalue_density and plot_eigenvalue_heatmap_helper were
        % MOVED to src/plotting/ on 2026-09-02. Both are pure functions over
        % arrays -- no model state, no class coupling -- and their only caller,
        % fig_eig_heatmap, applies them to SRNNCellTypePairs data, so living
        % here made the paper pipeline depend on this class for no reason.

        function y = piecewiseSigmoid(x, a, c)
            % PIECEWISESIGMOID A piecewise linear/quadratic sigmoid activation function.
            % Internalized from piecewiseSigmoid.m
            %
            % c may be a scalar (all neurons share a centre) or an array that
            % broadcasts against x -- typically an n x 1 per-neuron S_c_vec
            % against an n x 1 state or an n x nt trajectory.
            %
            % The curve is a pure TRANSLATION in c: every branch threshold is
            % c + const and the linear branch is (x - c) + 0.5, so it depends
            % on x and c only through x - c. A per-neuron centre is therefore
            % exactly the c = 0 curve evaluated at x - c, which keeps one
            % implementation of the branch structure and leaves the scalar path
            % below bit-identical.
            if ~isscalar(c)
                y = SRNNModel2.piecewiseSigmoid(x - c, a, 0);
                return;
            end

            if a < 0 || a > 1
                error('Parameter "a" must be between 0 and 1.');
            end
            a = a / 2;

            if a == 0.5
                y_linear = (x - c) + 0.5;
                y = min(max(y_linear, 0), 1);
            else
                y = zeros(size(x));
                k = 0.5 / (1 - 2*a);
                x1 = c + a - 1;
                x2 = c - a;
                x3 = c + a;
                x4 = c + 1 - a;
                
                mask_left_quad = (x >= x1) & (x < x2);
                mask_linear = (x >= x2) & (x <= x3);
                mask_right_quad = (x > x3) & (x <= x4);
                mask_right_sat = (x > x4);
                
                if any(mask_left_quad, 'all')
                    y(mask_left_quad) = k * (x(mask_left_quad) - x1).^2;
                end
                if any(mask_linear, 'all')
                    y(mask_linear) = (x(mask_linear) - c) + 0.5;
                end
                if any(mask_right_quad, 'all')
                    y(mask_right_quad) = 1 - k * (x(mask_right_quad) - x4).^2;
                end
                if any(mask_right_sat, 'all')
                    y(mask_right_sat) = 1;
                end
            end
        end
        
        function dy = piecewiseSigmoidDerivative(x, a, c)
            % PIECEWISESIGMOIDDERIVATIVE First derivative of the piecewise sigmoid.
            % Internalized from piecewiseSigmoidDerivative.m
            %
            % As with piecewiseSigmoid, c may be a scalar or an array that
            % broadcasts against x, and the derivative is likewise a pure
            % translation in c.
            if ~isscalar(c)
                dy = SRNNModel2.piecewiseSigmoidDerivative(x - c, a, 0);
                return;
            end

            if a < 0 || a > 1
                error('Parameter "a" must be between 0 and 1.');
            end
            a = a / 2;
            dy = zeros(size(x), 'like', x);

            if a == 0.5
                breakpoint1 = c - 0.5;
                breakpoint2 = c + 0.5;
                mask_linear = (x >= breakpoint1) & (x <= breakpoint2);
                dy(mask_linear) = 1;
            else
                k = 0.5 / (1 - 2*a);
                x1 = c + a - 1;
                x2 = c - a;
                x3 = c + a;
                x4 = c + 1 - a;
                
                mask_left_quad = (x >= x1) & (x < x2);
                mask_linear = (x >= x2) & (x <= x3);
                mask_right_quad = (x > x3) & (x <= x4);
                
                if any(mask_left_quad)
                    dy(mask_left_quad) = 2 * k * (x(mask_left_quad) - x1);
                end
                if any(mask_linear)
                    dy(mask_linear) = 1;
                end
                if any(mask_right_quad)
                    dy(mask_right_quad) = -2 * k * (x(mask_right_quad) - x4);
                end
            end
        end
        
        function y = logisticSigmoid(x, c)
            % LOGISTICSIGMOID Logistic sigmoid activation with unit slope at its center.
            %   y = logisticSigmoid(x, c) = 1/(1+exp(-4*(x-c)))
            %   Range (0,1); y(c)=0.5; slope 1 at x=c (the factor 4 sets that).
            %   c is an optional center/bias (default 0): c>0 shifts the
            %   inflection to positive x, so at the resting point (x~0) the
            %   network sits on the LOWER part of the curve. It may be a scalar
            %   or any array that broadcasts against x (an n x 1 per-neuron
            %   S_c_vec against an n x 1 state or an n x nt trajectory); the
            %   expression is elementwise, so no special case is needed.
            if nargin < 2, c = 0; end
            y = 1 ./ (1 + exp(-4 * (x - c)));
        end
        
        function dy = logisticSigmoidDerivative(x, c)
            % LOGISTICSIGMOIDDERIVATIVE First derivative of logisticSigmoid.
            %   dy = 4*sigma*(1-sigma) with sigma = 1/(1+exp(-4*(x-c))).
            %   c is an optional center/bias (default 0); dy(c)=1.
            if nargin < 2, c = 0; end
            sigmoid_val = 1 ./ (1 + exp(-4 * (x - c)));
            dy = 4 * sigmoid_val .* (1 - sigmoid_val);
        end
        
        function y = tanhActivation(x)
            % TANHACTIVATION Hyperbolic tangent activation function.
            % Internalized from tanhActivation.m
            y = tanh(x);
        end
        
        function dy = tanhActivationDerivative(x)
            % TANHACTIVATIONDERIVATIVE First derivative of the hyperbolic tangent.
            % Internalized from tanhActivationDerivative.m
            dy = 1 - tanh(x).^2;
        end
    end
    
    %% ====================================================================
    %              INTERNALIZED STIMULUS GENERATION
    % =====================================================================
    % Internalized from the standalone stimulus files  to make SRNNModel2 standalone.
    
    methods (Static, Access = protected)
        function [u_ex, t_ex] = generate_external_input(params, T, fs, rng_seed, input_config)
            % GENERATE_EXTERNAL_INPUT Generate external input for SRNN simulation.
            % Internalized from generate_external_input.m
            
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
    % Internalized from the standalone Jacobian files  to make SRNNModel2 standalone.
    
    methods (Static)
        function J = compute_Jacobian_fast(S, params)
            % COMPUTE_JACOBIAN_FAST Sparse/vectorized Jacobian assembly for the SRNN system.
            % Internalized from compute_Jacobian_fast.m
            
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
            
            c_E = SRNNModel2.effective_c(params, 'c_E', 'n_a_E');
            c_I = SRNNModel2.effective_c(params, 'c_I', 'n_a_I');
            
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
        
        function J_array = compute_Jacobian_at_indices(S_out, J_times, params)
            % COMPUTE_JACOBIAN_AT_INDICES Computes Jacobian matrices at multiple time indices.
            % Internalized from compute_Jacobian_at_indices.m
            
            N_sys_eqs = size(S_out, 2);
            n_times = length(J_times);
            
            nt = size(S_out, 1);
            if any(J_times < 1) || any(J_times > nt)
                error('J_times contains invalid indices. Must be between 1 and %d', nt);
            end
            
            J_array = zeros(N_sys_eqs, N_sys_eqs, n_times);
            
            for i = 1:n_times
                S = S_out(J_times(i), :)';
                J_array(:,:,i) = full(SRNNModel2.compute_Jacobian_fast(S, params));
            end
        end
    end
    
    % The solver registry -- solver_names, deterministic_solver_names,
    % stochastic_solver_names, resolve_solver, check_ode_solver and
    % check_noise_settings -- used to be statics here, duplicated verbatim on
    % SRNNCellTypePairs. They now live as shared functions in
    % src/model/integrators/ and are called by bare name from both classes.
    %
    % check_noise_settings is the one that mattered: SRNNCellTypePairs.validate
    % called it, and validate() runs from that class's constructor AND build(),
    % so every Pairs model in the repo depended on this class existing.
    %
    % This class passes 'SRNNModel' as the err_id_prefix, NOT 'SRNNModel2' --
    % a pre-existing inconsistency that test_ode_solver_name asserts, kept
    % deliberately.

    methods (Static, Access = private)
        function names = activation_names()
            % ACTIVATION_NAMES The valid values of the `activation` property.
            names = {'logistic', 'piecewise', 'tanh'};
        end

        function check_activation_custom(c)
            % CHECK_ACTIVATION_CUSTOM Validate the escape-hatch {fn, dfn} pair.
            ok = iscell(c) && numel(c) == 2 && ...
                all(cellfun(@(h) isa(h, 'function_handle'), c));
            if ~ok
                error('SRNNModel:InvalidParams', ...
                    ['activation_custom must be a 1x2 cell of function handles, ' ...
                    '{phi, phi_derivative}, or empty to use the named `activation`.']);
            end
        end

        function val = pick(override, fallback)
            % PICK The block override when set, else the column shorthand.
            if isempty(override), val = fallback; else, val = override; end
        end

        function check_setpoint_stat(val, name, allow_empty)
            % CHECK_SETPOINT_STAT Validate one mu_S_c_* / sigma_S_c_* entry.
            % A mu may be empty (meaning "fall back to S_c"); a sigma may not,
            % and must be nonnegative.
            if isempty(val)
                if allow_empty, return; end
                error('SRNNModel:InvalidParams', '%s must not be empty.', name);
            end
            if ~isnumeric(val) || ~isscalar(val) || ~isfinite(val)
                error('SRNNModel:InvalidParams', ...
                    '%s must be a finite numeric scalar.', name);
            end
            if ~allow_empty && val < 0
                error('SRNNModel:InvalidParams', ...
                    '%s must be nonnegative (it is a standard deviation).', name);
            end
        end

        function val = scale_tilde_mat(rel, F)
            % SCALE_TILDE_MAT Elementwise rel*F, keeping exact zeros at zero.
            % Same rationale as scale_tilde: alpha = 0 gives F = Inf, and 0*Inf
            % would be NaN.
            val = rel * F;
            val(rel == 0) = 0;
        end

        function val = scale_tilde(relative, F)
            % SCALE_TILDE relative * F, keeping an exact zero at zero.
            % A disconnected model (indegree = 0, e.g. the single-neuron scripts)
            % has alpha = 0, so F = 1/sqrt(0) = Inf and a plain 0*F would be NaN,
            % which would then poison W and the dynamics. A zero multiplier means
            % "this term is absent" regardless of the normalization.
            if relative == 0
                val = 0;
            else
                val = relative * F;
            end
        end

        function value = safe_get_param(params, field, default_value)
            % SAFE_GET_PARAM Helper to get a field from params with a default.
            if isfield(params, field)
                value = params.(field);
            else
                value = default_value;
            end
        end

        function c = effective_c(params, c_field, n_a_field)
            % EFFECTIVE_C The adaptation scaling actually applied: c / K.
            %
            % c is the TOTAL adaptation budget and K the number of timescales
            % sharing it. Every a_k relaxes to the rate, so sum_k a_k -> K*r at
            % steady state whatever the taus are; without the division the total
            % adaptation would scale linearly with K, and changing the number of
            % timescales in a condition would silently change adaptation
            % STRENGTH as well as timescale STRUCTURE.
            %
            % One helper rather than the expression inline, because four call
            % sites need it (the trajectory unpack, the dynamics, and two
            % Jacobian extractions) and a partial edit would leave the analytic
            % Jacobian disagreeing with the RHS.
            %
            % max(1, .) guards K = 0, where there are no a-states and the term is
            % absent anyway. (min(1, .) would not: min(1,0) is 0.)
            K = SRNNModel2.safe_get_param(params, n_a_field, 1);
            c = SRNNModel2.safe_get_param(params, c_field, 1.0) / max(1, K);
        end
    end
    
    %% ====================================================================
    %              INTERNALIZED PLOTTING FUNCTIONS
    % =====================================================================
    % Internalized from src/plotting/ to make SRNNModel2 standalone.
    
    methods (Static, Access = protected)
        function [fig_handle, ax_handles] = plot_SRNN_tseries(t_out, u, x, r, a, b, br, params, lya_results, Lya_method, T_plot)
            % PLOT_SRNN_TSERIES Create comprehensive time series plots for SRNN simulation.
            % Internalized from plot_SRNN_tseries.m (since deleted)
            
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
            % Internalized from plot_adaptation.m (since deleted)
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
            % Internalized from plot_std_variable.m (since deleted)
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
        
        function plot_lyapunov(lya_results, Lya_method, plot_options)
            % PLOT_LYAPUNOV Plot Lyapunov exponent(s) on current axes.
            % Internalized from src/plotting/plot_lyapunov.m
            
            if nargin < 3
                plot_options = {'local', 'asym', 'EOC', 'value'};
            end
            
            if strcmpi(Lya_method, 'benettin')
                valid_options = {'local', 'asym', 'EOC', 'value'};
                if ~iscell(plot_options)
                    error('plot_options must be a cell array of strings');
                end
                for i = 1:length(plot_options)
                    if strcmpi(plot_options{i}, 'filtered')
                        error(['The ''filtered'' option has been removed. ', ...
                            'Filtering now occurs in SRNNModel before decimation. ', ...
                            'Set model.filter_local_lya = true and use ''local'' instead.']);
                    end
                    if ~any(strcmpi(plot_options{i}, valid_options))
                        error('Invalid plot_option: %s. Valid options are: %s', ...
                            plot_options{i}, strjoin(valid_options, ', '));
                    end
                end
            end
            
            if strcmpi(Lya_method, 'benettin')
                plot_local = any(strcmpi('local', plot_options));
                plot_asym = any(strcmpi('asym', plot_options));
                plot_EOC = any(strcmpi('EOC', plot_options));
                plot_value = any(strcmpi('value', plot_options));
                
                legend_entries = {};
                plot_started = false;
                
                if plot_local
                    colors = lines(1);
                    plot(lya_results.t_lya, lya_results.local_lya, 'Color', colors(1,:))
                    hold on
                    plot_started = true;
                    legend_entries{end+1} = 'Local LLE';
                end
                
                if plot_asym
                    if ~plot_started, hold on; plot_started = true; end
                    plot([lya_results.t_lya(1), lya_results.t_lya(end)], ...
                        [lya_results.LLE, lya_results.LLE], '--r', 'LineWidth', 1.5);
                end
                
                if plot_EOC
                    if ~plot_started, hold on; plot_started = true; end
                    plot([lya_results.t_lya(1), lya_results.t_lya(end)], [0, 0], '--k');
                end
                
                ylabel('\lambda_1')
                
                if plot_value
                    ylims = ylim;
                    xlims = xlim;
                    text_y = 0.05 * (ylims(2) - ylims(1));
                    text_x = xlims(2);
                    text(text_x, text_y, ['$\lambda_1 = ' sprintf('%.2f', lya_results.LLE) '$'], ...
                        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
                        'Interpreter', 'latex');
                end
                
                hold off
                box off
                
            elseif strcmpi(Lya_method, 'qr')
                plot_data = lya_results.local_LE_spectrum_t(:, end:-1:1);
                line_handles = plot(lya_results.t_lya, plot_data);
                line_handles = line_handles(end:-1:1);
                
                hold on
                yline(0, '--k')
                ylabel('\lambda_1')
                
                legend_count = min(5, lya_results.params.N_sys_eqs);
                legend_entries = cell(1, legend_count);
                for i = 1:legend_count
                    legend_entries{i} = sprintf('\\lambda_{%d} = %.3f', i, lya_results.LE_spectrum(i));
                end
                legend(line_handles(1:legend_count), legend_entries, 'Location', 'best')
                hold off
                box off
            end
        end
        
        function plot_lines_with_colormap(t, data, cmap, varargin)
            % PLOT_LINES_WITH_COLORMAP Plot multiple lines with explicit colormap assignment.
            % Internalized from src/plotting/plot_lines_with_colormap.m
            if isempty(data), return; end
            n_lines = size(data, 1);
            n_colors = size(cmap, 1);
            hold on;
            for i = 1:n_lines
                color_idx = mod(i - 1, n_colors) + 1;
                plot(t, data(i, :), 'Color', cmap(color_idx, :), varargin{:});
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

