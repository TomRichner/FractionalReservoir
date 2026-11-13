classdef SRNN_ESN_reservoir < SRNNCellTypePairs
    % SRNN_ESN_RESERVOIR Echo State Network reservoir with memory capacity measurement
    %
    % This class extends SRNNCellTypePairs to provide Echo State Network (ESN)
    % functionality with the ability to measure memory capacity using the
    % protocol described in Memory_capacity_protocol.md.
    %
    % Usage:
    %   esn = SRNN_ESN_reservoir('n', 100, 'level_of_chaos', 1.8, ...
    %       'n_cellTypes', 2, 'cell_type_names', {'E','I'}, 'f', [0.6 0.4], ...
    %       'mu_tilde_relative', mu, 'sigma_tilde_relative', sg);
    %   esn.build();
    %   MC = esn.run_memory_capacity();
    %
    % Memory Capacity Protocol:
    %   - Drive reservoir with scalar random input u(t) ~ U(0,1)
    %   - For each delay d, train a linear readout to reconstruct u(t-d)
    %   - Compute R^2_d between true delayed input and readout output
    %   - Memory capacity = sum of R^2_d for d = 1 to d_max
    %
    % REPARENTED 2026-09-02, from SRNNModel2 to SRNNCellTypePairs. Memory
    % capacity was the last part of the paper on a different model class, which
    % meant two figures showed a network built by a different class from the
    % other fifteen and the methods section had to say so. Nothing about the
    % protocol needed SRNNModel2: the ESN overrides build_stimulus (protected on
    % both), writes S_out/t_out/S0/t_ex/u_interpolant/noise_increments (protected
    % SetAccess on both), and calls integrate/build_noise/get_params/
    % dynamics_fast (present on both). The only genuine coupling was eight lines
    % that named the E and I fields directly.
    %
    % CONSEQUENCE FOR CALLERS: the constructor now requires the SRNNCellTypePairs
    % arguments that have no sensible default -- n_cellTypes, cell_type_names, f,
    % mu_tilde_relative, sigma_tilde_relative -- and adaptation is configured with
    % tau_a (a 1 x C cell) and synapse_config rather than n_a_E / n_b_E. Note the
    % inversion: on SRNNModel2 n_a_E was set and tau_a_E auto-filled; here n_a is
    % Dependent on tau_a and is read-only.
    %
    % See also: SRNNCellTypePairs, run_memory_capacity, Memory_capacity_protocol.md

    %% ESN Input Properties
    properties
        f_in = 0.1              % Fraction of neurons receiving input
        sigma_in = 0.5          % Input weight scaling parameter
        rng_seed_input = 3      % RNG seed for input weight generation
        input_type = 'one_over_f'    % Input type: 'white', 'bandlimited', 'one_over_f', or 'sample_hold'
        u_f_cutoff = []         % Cutoff frequency for bandlimited input (Hz)
        % If empty, defaults to 1/(2*pi*tau_d)
        u_alpha = 1             % Spectral exponent for 1/f^alpha noise (default=1 for pink noise)
        T_hold = []             % 'sample_hold' hold duration (s); if empty, defaults to 3*tau_d
        u_scale = 1             % Stimulus amplitude scaling
        u_offset = 0            % Stimulus DC offset
    end

    %% Memory Capacity Protocol Properties
    properties
        T_wash = 1000           % Washout samples (discard transients)
        T_train = 5000          % Training samples
        T_test = 5000           % Test samples
        d_max = 70              % Maximum delay for memory capacity
        eta = 1e-7              % Ridge regression regularization
    end

    %% Dependent Properties
    properties (Dependent)
        dt_sample               % Sampling time step (= 1/fs)
    end

    %% Build-output Properties (set during build, not user-modifiable)
    properties (SetAccess = protected)
        W_in                    % Input weight vector (n x 1), set during build
        u_scalar                % Scalar input sequence (T_total x 1), set during build
        hold_len = 1            % Samples per held input value (1 unless input_type='sample_hold')
    end

    %% Memory Capacity Results
    properties (SetAccess = private)
        mc_results              % Struct with memory capacity results
    end

    %% Constructor
    methods
        function obj = SRNN_ESN_reservoir(varargin)
            % SRNN_ESN_RESERVOIR Constructor with name-value pairs
            %
            % Usage:
            %   esn = SRNN_ESN_reservoir()  % All defaults
            %   esn = SRNN_ESN_reservoir('n', 200, 'd_max', 100)

            % Call superclass constructor first (MATLAB requirement).
            % SRNNCellTypePairs REQUIRES n_cellTypes, cell_type_names, f,
            % mu_tilde_relative and sigma_tilde_relative -- the class is general
            % over cell types, so there is no sensible default for them and it
            % errors rather than guessing.
            obj = obj@SRNNCellTypePairs(varargin{:});

            % Define ESN-specific property names (not in SRNNCellTypePairs)
            esn_props = {'f_in', 'sigma_in', 'rng_seed_input', ...
                'T_wash', 'T_train', 'T_test', 'd_max', 'eta', ...
                'input_type', 'u_f_cutoff', 'u_alpha', 'u_scale', 'u_offset', 'T_hold'};

            % Parse ESN-specific name-value pairs
            for i = 1:2:length(varargin)
                if ismember(varargin{i}, esn_props)
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
    end

    %% Dependent Property Getters
    methods
        function val = get.dt_sample(obj)
            val = 1 / obj.fs;
        end
    end

    %% Public Methods
    methods
        function generate_input_weights(obj)
            % GENERATE_INPUT_WEIGHTS Create the input weight vector W_in
            %
            % Creates a sparse input weight vector where f_in fraction of
            % neurons receive input, with weights drawn uniformly from
            % [-sigma_in/2, sigma_in/2].

            rng(obj.rng_seed_input);

            % Initialize input weights to zero
            obj.W_in = zeros(obj.n, 1);

            % Select neurons to receive input
            n_input = round(obj.f_in * obj.n);
            input_neurons = randperm(obj.n, n_input);

            % Assign random weights centered at 0 with spread sigma_in
            obj.W_in(input_neurons) = obj.sigma_in * (rand(n_input, 1) - 0.5);

            vprintf(obj.verbose, 'verbose', 'Input weights generated: %d neurons receive input (%.1f%%)\n', ...
                n_input, 100 * n_input / obj.n);
        end

        function [MC, R2_d, mc_results] = run_memory_capacity(obj, varargin)
            % RUN_MEMORY_CAPACITY Measure memory capacity of the reservoir
            %
            % [MC, R2_d, results] = run_memory_capacity()
            % [MC, R2_d, results] = run_memory_capacity('verbose', true)
            % [MC, R2_d, results] = run_memory_capacity('store_timeseries', false)
            %
            % Name-value options:
            %   'verbose'          - print progress. Default: the inherited
            %                        verbose property at 'verbose' (see verbose_level);
            %                        a logical or a level name overrides it
            %   'store_timeseries' - default true. Set FALSE for a lean MC-only run
            %                        that does NOT keep the full-resolution state
            %                        time series (x/a/b/r/br, u_ex) or S_out and
            %                        skips the LLE -- much lower memory for batch /
            %                        parallel runs. MC and R2_d are identical either
            %                        way. Do not use with plot_esn_timeseries, which
            %                        needs the stored states.
            %
            % Outputs:
            %   MC        - Total memory capacity (sum of R^2_d)
            %   R2_d      - R^2 for each delay d (1 x d_max)
            %   results   - Struct with detailed results
            %
            % Protocol:
            %   1. Generate scalar random input sequence u ~ U(0,1)
            %   2. Run reservoir with piecewise-constant input
            %   3. Collect reservoir states (firing rates)
            %   4. Discard washout period
            %   5. Train linear readouts for each delay via ridge regression
            %   6. Compute R^2_d on test set
            %   7. Sum R^2_d to get memory capacity

            if ~obj.is_built
                error('SRNN_ESN_reservoir:NotBuilt', ...
                    'Reservoir must be built first. Call build().');
            end

            % Parse optional arguments
            verbose = verbose_level(obj.verbose) >= 2;
            store_timeseries = true;   % false => lean MC-only path (frees the
                                       % full-res states/trajectory, skips the LLE)
            readout_signal = 'rate';   % 'rate' => read out r = phi(x_eff); 'synaptic'
                                       % => read out the synaptic output, r times
                                       % the route's output factor
            for i = 1:2:length(varargin)
                switch lower(varargin{i})
                    case 'verbose',          verbose = verbose_level(varargin{i+1}) >= 2;
                    case 'store_timeseries', store_timeseries = varargin{i+1};
                    case 'readout_signal',   readout_signal = varargin{i+1};
                end
            end

            if verbose
                fprintf('Running memory capacity measurement...\n');
                fprintf('  Washout: %d, Train: %d, Test: %d samples\n', ...
                    obj.T_wash, obj.T_train, obj.T_test);
                fprintf('  Max delay: %d\n', obj.d_max);
            end

            %% Step 1: Run reservoir and collect states
            if verbose
                fprintf('  Running reservoir simulation...\n');
            end

            obj.run_reservoir_esn();
            t_all = obj.t_out;

            % Unpack states using standard utility. SRNNCellTypePairs returns
            % SIX outputs, all keyed by cell-type name, and its fifth is called
            % synaptic_output rather than br -- with STF in play the quantity is
            % not b*r.
            [x_all, a_all, b_all, r_all, syn_all, g_all] = ...
                SRNNCellTypePairs.unpack_and_compute_states(obj.S_out, obj.cached_params);

            % The signal the readout sees, every cell type stacked into one
            % n x T matrix in TYPE ORDER: 'rate' or 'synaptic'. One definition,
            % readout_from_states, shared with run_memory_capacity_lean.
            params = obj.cached_params;
            R_all = SRNN_ESN_reservoir.readout_from_states(r_all, syn_all, params, readout_signal);

            %% Step 2b: Sample-and-hold -> decimate to one state per hold so MC is
            % measured in hold units (canonical i.i.d.-input MC). Other input
            % types keep sample resolution (hold_len = 1, all locals = obj.*).
            % The decimation, the split and the delay grid are
            % mc_readout_plan's, shared with run_memory_capacity_lean.
            plan = obj.mc_readout_plan();
            if obj.hold_len > 1
                R_all = R_all(:, plan.idx);         % last sample of each hold
            end

            % Lean path: the full-resolution states and the raw trajectory are no
            % longer needed once R_all (decimated for sample_hold) is built. Free
            % them now to keep per-worker memory low in batch/parallel MC runs.
            if ~store_timeseries
                x_all = []; a_all = []; b_all = []; r_all = [];
                syn_all = []; g_all = [];
                obj.S_out = [];
            end

            %% Steps 3-6: washout, split, one ridge readout per delay, R^2, MC.
            % score_readout is the ONE implementation of the scoring; the lean
            % path (run_memory_capacity_lean) calls the same method.
            [MC, R2_d, weights_all, predictions, u_test] = ...
                obj.score_readout(R_all, plan, verbose);
            % Real times of the (possibly decimated) test samples. predictions
            % index into 1:n_test, so the delay panels plot against t_pred, not
            % the full-resolution t_test (which would mis-scale under sample_hold).
            t_pred = t_all(plan.idx(plan.te_idx));                % n_test x 1

            %% Store results including time series for plotting
            mc_results = obj.pack_mc_results(MC, R2_d, plan, weights_all, ...
                predictions, t_pred, u_test);

            if store_timeseries
                % Full-resolution state time series + LLE, for plot_esn_timeseries.
                test_start_idx = obj.T_wash + obj.T_train + 1;
                mc_results.t_test = t_all(test_start_idx:end);
                mc_results.u_ex = obj.u_ex;   % actual neural input (n x T)
                mc_results.x = x_all;
                mc_results.a = a_all;
                mc_results.b = b_all;
                mc_results.r = r_all;
                % Named synaptic_output, not br: it nests by ROUTE (.E.E) and
                % with STF in play the quantity is not b*r. g is kept alongside
                % it so a stored run can show facilitation separately.
                mc_results.synaptic_output = syn_all;
                mc_results.g = g_all;

                %% Step 7: Compute Lyapunov exponents using parent class method
                if ~strcmpi(obj.lya_method, 'none')
                    obj.lya_T_interval = [t_all(obj.T_wash + 1), t_all(end)];  % after washout
                    obj.compute_lyapunov();
                end
                mc_results.lya_results = obj.lya_results;
            else
                % Lean path: heavy states/trajectory were freed above; the LLE
                % needs S_out and is not part of MC, so it is skipped.
                mc_results.lya_results = [];
            end

            obj.mc_results = mc_results;
        end

        function [MC, R2_d, mc_results] = run_memory_capacity_lean(obj, varargin)
            % RUN_MEMORY_CAPACITY_LEAN Memory capacity without the state trajectory.
            %
            % [MC, R2_d, results] = run_memory_capacity_lean()
            % [...] = run_memory_capacity_lean('readout_signal', 'synaptic', ...
            %             'segment_len', 600, 'verbose', false)
            %
            % The same task, input, integrator and scoring as
            % run_memory_capacity, which integrates the whole sequence in one
            % call and holds every state at every step (at n = 500 and 760 s,
            % 152,000 rows by about 4,000 states, about 5 GB, then the
            % per-variable copies) only to keep one readout column per hold.
            % This path integrates SEGMENT BY SEGMENT (segment_len samples,
            % default 600), keeps the readout signal at the samples the scoring
            % uses, accumulates each neuron's time-mean rate, and carries the
            % last state into the next segment. What it holds is the readout
            % matrix (n by the number of holds), one segment of states, and the
            % input the reservoir built (u_ex and its interpolant, n by the
            % number of samples each: the largest thing left).
            %
            % EXACT, NOT APPROXIMATE, with a fixed-step integrator ('rk4',
            % 'euler', 'heun', 'sra1'). Those take each step from the time
            % grid they are handed (h = t(k+1) - t(k), stages at t(k) + c h)
            % and from the state alone, so cutting the grid into segments that
            % share their end points changes no step: the states are the same
            % bits, and so are the readout matrix, the R^2 curve and MC
            % (test_mc_lean_path asserts isequal). With noise on it is exact
            % too: the Wiener increments are generated once for the whole
            % sequence (build_noise) and sde_fixed_step addresses them by
            % ABSOLUTE time, so a segment consumes the increments the single
            % call would have. (The noise tensor is 2 n by the number of
            % steps, so with noise on this path is exact but less lean.) With
            % an adaptive integrator ('ode45', 'ode15s') a segment boundary
            % restarts the step-size control, and the result differs within
            % the tolerance.
            %
            % THE SCORING IS SHARED: mc_readout_plan (decimation, split, delay
            % grid), readout_from_states (which signal is read out),
            % score_readout (ridge readout per delay, R^2, MC) and
            % pack_mc_results are the methods run_memory_capacity calls.
            %
            % Name-value options:
            %   'readout_signal'  'rate' (default) | 'synaptic', as run_memory_capacity
            %   'segment_len'     samples integrated per call (default 600; >= 2)
            %   'verbose'         as run_memory_capacity
            %   'trace_window'    [t_lo t_hi] in seconds, default [] (nothing).
            %                     Keeps, for the samples in that window only,
            %                     what a figure of the time series draws:
            %                     results.trace with t, the input u_ex, and
            %                     x, a, b, r, synaptic_output, g keyed as
            %                     run_memory_capacity's stored series are.
            %
            % results has run_memory_capacity's fields (lya_results empty, no
            % state series), and:
            %   rate_mean    n x 1, each neuron's rate averaged over the
            %                training and test windows, at full sample
            %                resolution, neurons in type order
            %   rate_window  [first last] sample of that average
            %   lean         true
            %
            % See also: run_memory_capacity, score_readout, run_mc_gain_sweep

            if ~obj.is_built
                error('SRNN_ESN_reservoir:NotBuilt', ...
                    'Reservoir must be built first. Call build().');
            end
            verbose = verbose_level(obj.verbose) >= 2;
            readout_signal = 'rate';
            segment_len = 600;
            trace_window = [];
            for i = 1:2:length(varargin)
                switch lower(varargin{i})
                    case 'verbose',        verbose = verbose_level(varargin{i+1}) >= 2;
                    case 'readout_signal', readout_signal = varargin{i+1};
                    case 'segment_len',    segment_len = varargin{i+1};
                    case 'trace_window',   trace_window = varargin{i+1};
                    otherwise
                        error('SRNN_ESN_reservoir:BadOption', ...
                            'Unknown option ''%s''.', varargin{i});
                end
            end
            if ~(isscalar(segment_len) && segment_len >= 2 && segment_len == round(segment_len))
                error('SRNN_ESN_reservoir:BadOption', ...
                    'segment_len must be an integer of at least 2 (a fixed-step integrator needs three grid points).');
            end

            plan = obj.mc_readout_plan();
            params = obj.cached_params;
            t_all = obj.t_ex;
            T_total = numel(t_all);

            % The right-hand side and the noise, exactly as run_reservoir_esn
            % sets them up.
            if isempty(obj.ode_opts)
                obj.ode_opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-5, 'MaxStep', 1 / obj.fs);
            end
            rhs_params = params;
            rhs_params.u_interpolant = obj.u_interpolant;
            rhs = @(t, S) SRNNCellTypePairs.dynamics_fast(t, S, rhs_params);
            obj.build_noise();

            % Which samples are readout columns, and which count towards the
            % time-mean rate (the training and test windows, full resolution).
            n_cols = numel(plan.idx);
            col_of = zeros(T_total, 1);
            col_of(plan.idx) = 1:n_cols;
            L = obj.hold_len;
            rate_lo = plan.n_wash * L + 1;
            rate_hi = min((plan.n_wash + plan.n_train + plan.n_test) * L, T_total);
            R_dec = zeros(obj.n, n_cols);
            rate_sum = zeros(obj.n, 1);
            n_rate = 0;

            % Optional trace: the samples of one window, for a figure.
            trace_rows = [];
            if ~isempty(trace_window)
                trace_rows = find(t_all >= trace_window(1) & t_all <= trace_window(2));
            end
            S_trace = zeros(numel(trace_rows), numel(obj.S0));
            n_traced = 0;

            if verbose
                fprintf('Running memory capacity (lean): %d samples in segments of %d\n', ...
                    T_total, segment_len);
            end
            S_cur = obj.S0;
            s0 = 1;
            while s0 < T_total
                s1 = min(s0 + segment_len, T_total);
                if T_total - s1 < 2
                    s1 = T_total;       % never leave a stub of fewer than three points
                end
                [~, S_seg] = obj.integrate(rhs, t_all(s0:s1), S_cur);
                S_cur = S_seg(end, :).';
                if s0 == 1
                    rows = (s0:s1).';
                else
                    rows = (s0 + 1:s1).';       % row s0 was the previous segment's last
                    S_seg = S_seg(2:end, :);
                end

                [~, ~, ~, r_seg, syn_seg] = ...
                    SRNNCellTypePairs.unpack_and_compute_states(S_seg, params);
                is_col = col_of(rows) > 0;
                if any(is_col)
                    R_seg = SRNN_ESN_reservoir.readout_from_states(r_seg, syn_seg, params, readout_signal);
                    R_dec(:, col_of(rows(is_col))) = R_seg(:, is_col);
                end
                in_rate = rows >= rate_lo & rows <= rate_hi;
                if any(in_rate)
                    rate_seg = SRNN_ESN_reservoir.stack_by_type(r_seg, params);
                    rate_sum = rate_sum + sum(rate_seg(:, in_rate), 2);
                    n_rate = n_rate + nnz(in_rate);
                end
                if ~isempty(trace_rows)
                    in_trace = rows >= trace_rows(1) & rows <= trace_rows(end);
                    k = nnz(in_trace);
                    S_trace(n_traced + (1:k), :) = S_seg(in_trace, :);
                    n_traced = n_traced + k;
                end
                s0 = s1;
            end
            obj.noise_increments = [];
            obj.t_out = t_all;
            obj.S_out = [];

            [MC, R2_d, weights_all, predictions, u_test] = ...
                obj.score_readout(R_dec, plan, verbose);
            t_pred = t_all(plan.idx(plan.te_idx));
            mc_results = obj.pack_mc_results(MC, R2_d, plan, weights_all, ...
                predictions, t_pred, u_test);
            mc_results.lya_results = [];
            mc_results.rate_mean = rate_sum / max(1, n_rate);
            mc_results.rate_window = [rate_lo, rate_hi];
            mc_results.lean = true;
            if ~isempty(trace_rows)
                tr = struct();
                tr.t = t_all(trace_rows);
                tr.rows = trace_rows;
                tr.u_ex = obj.u_ex(:, trace_rows);
                [tr.x, tr.a, tr.b, tr.r, tr.synaptic_output, tr.g] = ...
                    SRNNCellTypePairs.unpack_and_compute_states(S_trace, params);
                mc_results.trace = tr;
            end
            obj.mc_results = mc_results;
        end

        function plan = mc_readout_plan(obj)
            % MC_READOUT_PLAN Which samples are scored, and how they are split.
            %
            % One definition for run_memory_capacity and
            % run_memory_capacity_lean. Sample-and-hold input is scored in HOLD
            % units: one readout column per hold, at its last sample (canonical
            % i.i.d.-input MC). Other input types keep sample resolution.
            %
            %   idx          full-resolution sample of each readout column
            %   u            the input at those samples
            %   n_wash, n_train, n_test   the split, in readout columns
            %   d_max_eff    delays scored;  delay_scale  seconds per delay
            %   te_idx       readout columns of the test window
            T_total = numel(obj.u_scalar);
            plan = struct();
            if obj.hold_len > 1
                plan.idx = obj.hold_len:obj.hold_len:T_total;  % last sample of each hold
                plan.u = obj.u_scalar(plan.idx);
                plan.n_wash  = round(obj.T_wash  / obj.hold_len);
                plan.n_train = round(obj.T_train / obj.hold_len);
                plan.n_test  = round(obj.T_test  / obj.hold_len);
                plan.d_max_eff   = max(1, floor(obj.d_max / obj.hold_len));
                plan.delay_scale = obj.hold_len / obj.fs;   % seconds per hold-delay
            else
                plan.idx = 1:T_total;               % identity (sample resolution)
                plan.u = obj.u_scalar;
                plan.n_wash = obj.T_wash;  plan.n_train = obj.T_train;  plan.n_test = obj.T_test;
                plan.d_max_eff = obj.d_max;  plan.delay_scale = 1 / obj.fs;
            end
            % Clamp the test window so the requested split fits the sequence.
            n_avail = numel(plan.u) - plan.n_wash;
            if plan.n_train + plan.n_test > n_avail
                plan.n_test = max(1, n_avail - plan.n_train);
            end
            plan.te_idx = (plan.n_wash + plan.n_train + 1):(plan.n_wash + plan.n_train + plan.n_test);
        end

        function [MC, R2_d, weights_all, predictions, u_test] = score_readout(obj, R_all, plan, verbose)
            % SCORE_READOUT Memory capacity from the readout matrix.
            %
            % R_all is n by the number of readout columns (plan.idx), however
            % it was collected. Discards the washout, splits into training and
            % test, trains one ridge readout per delay on the training window
            % and scores its R^2 on the test window; MC is the sum.

            n_wash = plan.n_wash;  n_train = plan.n_train;  n_test = plan.n_test;
            d_max_eff = plan.d_max_eff;

            %% Step 3: Discard washout and split into train/test
            % Reservoir states after washout
            R_eff = R_all(:, (n_wash + 1):end);  % n x (>= T_eff)
            u_eff = plan.u((n_wash + 1):end);    % (>= T_eff) x 1

            % Split into training and test
            R_train = R_eff(:, 1:n_train)';                       % n_train x n
            R_test  = R_eff(:, (n_train + 1):(n_train + n_test))';% n_test x n
            u_train = u_eff(1:n_train);                           % n_train x 1
            u_test  = u_eff((n_train + 1):(n_train + n_test));    % n_test x 1

            %% Step 4-5: Train readouts and compute R^2 for each delay
            if verbose
                fprintf('  Training readouts for %d delays...\n', d_max_eff);
            end

            R2_d = zeros(1, d_max_eff);
            weights_all = zeros(obj.n, d_max_eff);

            % Store predictions for selected delays (for plotting)
            predictions = struct();

            for d = 1:d_max_eff
                % Build delayed targets for training
                % Target: u(t-d) where t is current time
                % For training indices (d+1):n_train, target is u(1:(n_train-d))
                train_indices = (d + 1):n_train;
                target_indices = 1:(n_train - d);

                if length(train_indices) < 10
                    warning('SRNN_ESN_reservoir:InsufficientData', ...
                        'Delay %d has insufficient training data.', d);
                    continue;
                end

                X_train = R_train(train_indices, :);
                y_train = u_train(target_indices);

                % Train linear readout with ridge regression
                w_d = obj.train_linear_readout(X_train, y_train, obj.eta);
                weights_all(:, d) = w_d;

                % Build delayed targets for test
                % For test indices (d+1):n_test, target is u_test(1:(n_test-d))
                test_indices = (d + 1):n_test;
                test_target_indices = 1:(n_test - d);

                if length(test_indices) < 10
                    continue;
                end

                X_test = R_test(test_indices, :);
                y_test_true = u_test(test_target_indices);

                % Compute predictions
                y_test_pred = X_test * w_d;

                % Store predictions for this delay
                predictions(d).y_true = y_test_true;
                predictions(d).y_pred = y_test_pred;
                predictions(d).t_indices = test_indices;

                % Compute R^2
                R2_d(d) = obj.compute_R2(y_test_true, y_test_pred);

                if verbose && mod(d, 10) == 0
                    fprintf('    Delay %d: R^2 = %.4f\n', d, R2_d(d));
                end
            end

            %% Step 6: Compute total memory capacity
            MC = sum(R2_d);

            if verbose
                fprintf('  Total Memory Capacity: %.4f\n', MC);
            end
        end

        function mc_results = pack_mc_results(obj, MC, R2_d, plan, weights_all, predictions, t_pred, u_test)
            % PACK_MC_RESULTS The result fields both memory-capacity paths share.
            mc_results = struct();
            mc_results.MC = MC;
            mc_results.R2_d = R2_d;
            mc_results.d = 1:plan.d_max_eff;            % delay index (holds if sample_hold, else samples)
            mc_results.hold_len = obj.hold_len;         % samples per held value (1 = sample resolution)
            mc_results.delay_s = (1:plan.d_max_eff) * plan.delay_scale;  % delay axis in seconds
            mc_results.weights = weights_all;
            mc_results.T_wash = obj.T_wash;
            mc_results.T_train = obj.T_train;
            mc_results.T_test = obj.T_test;
            mc_results.eta = obj.eta;
            mc_results.u_scalar = obj.u_scalar;   % full input (staircase for sample_hold)
            mc_results.predictions = predictions;
            mc_results.t_pred = t_pred;   % time axis matching predictions (hold-aware)
            mc_results.u_test = u_test;
        end

        function run_reservoir_esn(obj)
            % RUN_RESERVOIR_ESN Run reservoir with single ODE integration
            %
            % Uses pre-built stimulus (t_ex, u_interpolant, S0) from build().
            %
            % Results are stored in inherited properties:
            %   obj.S_out - State trajectory (T x N_sys_eqs)
            %   obj.t_out - Time vector (T x 1)

            params = obj.cached_params;
            dt = 1 / obj.fs;

            % Set up ODE options (matching SRNNModel pattern)
            if isempty(obj.ode_opts)
                obj.ode_opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-5, 'MaxStep', dt);
            end

            % Define RHS function using static method (avoids OOP overhead)
            params.u_interpolant = obj.u_interpolant;
            rhs = @(t, S) SRNNCellTypePairs.dynamics_fast(t, S, params);

            % Integrate entire trajectory at once
            vprintf(obj.verbose, 'verbose', '  Integrating ESN dynamics...\n');
            tic
            % Pre-generate Wiener increments (a no-op at sigma_u_noise = 0).
            % run_reservoir_esn does not go through SRNNCellTypePairs.run, so it has
            % to do this itself; integrate() then picks the tensor up.
            obj.build_noise();
            [obj.t_out, obj.S_out] = obj.integrate(rhs, obj.t_ex, obj.S0);
            obj.noise_increments = [];
            integration_time = toc;
            vprintf(obj.verbose, 'verbose', '  Integration complete in %.2f seconds.\n', integration_time);
        end

        function [fig_handle, ax_handles] = plot_memory_capacity(obj, varargin)
            % PLOT_MEMORY_CAPACITY Plot memory capacity results
            %
            % [fig, axes] = plot_memory_capacity()
            %
            % Creates a figure showing:
            %   - R^2_d vs delay d
            %   - Total MC value

            if isempty(obj.mc_results)
                error('SRNN_ESN_reservoir:NoResults', ...
                    'No memory capacity results. Run run_memory_capacity() first.');
            end

            fig_handle = figure();

            % Plot R^2 vs delay
            ax1 = subplot(1, 2, 1);
            bar(obj.mc_results.d, obj.mc_results.R2_d, 'FaceColor', [0.3, 0.6, 0.9]);
            xlabel('Delay d (samples)');
            ylabel('R^2_d');
            title('Memory Capacity by Delay');
            grid on;
            ax_handles(1) = ax1;

            % Plot cumulative MC
            ax2 = subplot(1, 2, 2);
            cumMC = cumsum(obj.mc_results.R2_d);
            plot(obj.mc_results.d, cumMC, 'b-', 'LineWidth', 2);
            xlabel('Delay d (samples)');
            ylabel('Cumulative MC');
            title(sprintf('Cumulative Memory Capacity (Total: %.2f)', obj.mc_results.MC));
            grid on;
            ax_handles(2) = ax2;

            sgtitle(sprintf('Memory Capacity Analysis (MC = %.2f)', obj.mc_results.MC));
        end

        function [fig_handle, ax_handles] = plot_esn_timeseries(obj, delays_to_plot, varargin)
            % PLOT_ESN_TIMESERIES Plot comprehensive time series for ESN reservoir
            %
            % [fig, axes] = plot_esn_timeseries(delays_to_plot)
            % [fig, axes] = plot_esn_timeseries(delays_to_plot, 'title', 'My Title')
            %
            % Inputs:
            %   delays_to_plot - Vector of delays to show (e.g., [1, 10, 30, 50])
            %
            % Creates a figure with:
            %   - u(t): Scalar input
            %   - x(t): Dendritic states (line plots)
            %   - r(t): Firing rates (line plots)
            %   - br(t): Synaptic output (if STD enabled)
            %   - a(t): Adaptation states (if SFA enabled)
            %   - b(t): STD states (if STD enabled)
            %   - Lyapunov exponent (if computed)
            %   - For each delay d: u(t-d) vs y_d(t) overlay

            if isempty(obj.mc_results) || ~isfield(obj.mc_results, 'x')
                error('SRNN_ESN_reservoir:NoResults', ...
                    'No time series data. Run run_memory_capacity() first.');
            end

            if nargin < 2 || isempty(delays_to_plot)
                delays_to_plot = [1, 10, 30, 50];
            end

            % Filter delays that exist in results
            delays_to_plot = delays_to_plot(delays_to_plot <= obj.d_max);
            n_delays = length(delays_to_plot);

            % Parse optional arguments
            fig_title = '';
            for i = 1:2:length(varargin)
                if strcmpi(varargin{i}, 'title')
                    fig_title = varargin{i+1};
                end
            end

            % Stored data, all keyed by cell-type name.
            t = obj.mc_results.t_test;
            u_ex = obj.mc_results.u_ex;   % Actual neural input (n x T)
            x = obj.mc_results.x;
            r = obj.mc_results.r;
            syn = obj.mc_results.synaptic_output;   % nests by route: .E.E
            a = obj.mc_results.a;
            b = obj.mc_results.b;                   % nests by route: .E.E
            params = obj.cached_params;
            names = params.cell_type_names;
            C = params.n_cellTypes;

            % Slice to test period only
            test_start_idx = obj.T_wash + obj.T_train + 1;
            test_end_idx = size(x.(names{1}), 2);
            test_indices = test_start_idx:test_end_idx;

            % Determine which subplots are needed
            has_adaptation = any(params.n_a > 0);
            has_std = any(params.n_b_pairs(:) > 0);

            % Check for Lyapunov results
            has_lyapunov = isfield(obj.mc_results, 'lya_results') && ...
                ~isempty(obj.mc_results.lya_results) && ...
                isfield(obj.mc_results.lya_results, 'LLE');

            % Calculate number of base plots
            n_base_plots = 3;  % u(t), x(t), r(t) always present
            if has_std
                n_base_plots = n_base_plots + 1;  % synaptic output
            end
            if has_adaptation
                n_base_plots = n_base_plots + 1;  % a(t)
            end
            if has_std
                n_base_plots = n_base_plots + 1;  % b(t)
            end
            if has_lyapunov
                n_base_plots = n_base_plots + 1;  % Lyapunov exponent
            end

            n_total_plots = n_base_plots + n_delays;

            % Create figure
            fig_height = min(100 + 80 * n_total_plots, 1000);
            fig_handle = figure();
            tiledlayout(n_total_plots, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

            ax_handles = [];

            % Per-type colours and layering, from the model class's own helpers
            % rather than a local E/I pair -- type 1 warm, type 2 cool, further
            % types from lines(). draw_order is n:-1:1, so type 1 is drawn LAST
            % and lands on top; colours stay indexed by type, so the reversal is
            % layering only.
            tcolors = SRNNCellTypePairs.type_colors(C);
            order   = SRNNCellTypePairs.draw_order(C);

            %% Plot 1: Actual neural input u_ex(t)
            ax_handles(end+1) = nexttile;
            % Plot actual neural input for neurons that receive input
            u_ex_test = u_ex(:, test_indices);
            % Only plot neurons with non-zero input weights
            input_mask = any(u_ex_test ~= 0, 2);
            u_ex_active = u_ex_test(input_mask, :);
            if ~isempty(u_ex_active)
                % Use a neutral colormap for input
                n_active = size(u_ex_active, 1);
                % The standalone primitives in src/plotting/, not a model class's
                % static copy: this class no longer inherits SRNNModel2's, and
                % SRNNCellTypePairs has none. Calling the shared functions is
                % what keeps the ESN independent of either class's plotting.
                cmap_input = parula(max(n_active, 8));
                plot_lines_with_colormap(t, u_ex_active, cmap_input);
            end
            ylabel('u_{ex}(t)');
            set(gca, 'XTickLabel', []);
            grid on;

            %% Plot 2: Dendritic states x(t) as line plots
            ax_handles(end+1) = nexttile;
            obj.plot_types(t, x, test_indices, names, tcolors, order);
            ylabel('dendrite');
            set(gca, 'XTickLabel', []);

            %% Plot 3: Firing rates r(t) as line plots
            ax_handles(end+1) = nexttile;
            obj.plot_types(t, r, test_indices, names, tcolors, order);
            ylabel('firing rate');
            ylim([0, 1]);
            yticks([0, 1]);
            set(gca, 'XTickLabel', []);

            %% Plot 4 (conditional): Synaptic output
            % Per route on this class, so take the FIRST outgoing route of each
            % presynaptic type. Unlike the readout this does not assert
            % redundancy: a diagnostic plot showing one representative route is
            % informative even when the routes differ, and refusing to draw
            % would be worse than drawing a labelled subset.
            if has_std
                ax_handles(end+1) = nexttile;
                syn_first = struct();
                for q = 1:C
                    posts = fieldnames(syn.(names{q}));
                    syn_first.(names{q}) = syn.(names{q}).(posts{1});
                end
                obj.plot_types(t, syn_first, test_indices, names, tcolors, order);
                ylabel('synaptic output');
                ylim([0, 1]);
                yticks([0, 1]);
                set(gca, 'XTickLabel', []);
            end

            %% Plot 5 (conditional): Adaptation states a(t)
            if has_adaptation
                ax_handles(end+1) = nexttile;
                a_sum = struct();
                for q = 1:C
                    aq = a.(names{q});
                    if isempty(aq) || params.n_a(q) == 0
                        a_sum.(names{q}) = zeros(0, numel(test_indices));
                        continue
                    end
                    % Sum across timescales: (nq x n_a x T) -> (nq x T)
                    s = reshape(sum(aq(:, :, test_indices), 2), ...
                        params.n_per_type(q), numel(test_indices));
                    a_sum.(names{q}) = s;
                end
                obj.plot_types(t, a_sum, [], names, tcolors, order);
                ylabel('adaptation');
                set(gca, 'XTickLabel', []);
                grid on;
            end

            %% Plot 6 (conditional): depression states b(t)
            if has_std
                ax_handles(end+1) = nexttile;
                b_prod = struct();
                for q = 1:C
                    % Each route's depression FACTOR, then the product across
                    % the type's routes. The factor is the route model's
                    % (SRNNCellTypePairs.std_route_factor), not prod(b, 2)
                    % formed here: on a plain product route the two are the
                    % same array, on a chain route the factor is b1, the
                    % releasable pool, and NOT b1 b2. (A product route with
                    % combine = 'geomean' or the zero floor also draws its
                    % factor now, where it drew the bare product before.)
                    % posts is in type order (unpack_and_compute_states
                    % fills it so), so p is the postsynaptic type's index.
                    posts = fieldnames(b.(names{q}));
                    acc = [];
                    for p = 1:numel(posts)
                        bq = b.(names{q}).(posts{p});
                        if isempty(bq); continue; end
                        pr = reshape(SRNNCellTypePairs.std_route_factor(bq, params, q, p), ...
                            params.n_per_type(q), []);
                        if isempty(acc); acc = pr; else; acc = acc .* pr; end
                    end
                    if isempty(acc) || all(acc(:) == 1)
                        b_prod.(names{q}) = zeros(0, numel(test_indices));
                    else
                        b_prod.(names{q}) = acc(:, test_indices);
                    end
                end
                obj.plot_types(t, b_prod, [], names, tcolors, order);
                ylabel('depression');
                ylim([0, 1]);
                yticks([0, 1]);
                set(gca, 'XTickLabel', []);
                grid on;
            end

            %% Plot (conditional): Lyapunov exponent
            if has_lyapunov
                ax_handles(end+1) = nexttile;
                plot_lyapunov(obj.mc_results.lya_results, 'benettin', {'local', 'EOC'});
                % Add LLE value to subplot title
                title(sprintf('\\lambda_1 = %.2f', obj.mc_results.lya_results.LLE), 'FontWeight', 'normal');
                set(gca, 'XTickLabel', []);
            end

            %% Delay reconstruction plots
            colors_pred = lines(n_delays);

            for i = 1:n_delays
                d = delays_to_plot(i);
                ax_handles(end+1) = nexttile;

                if d <= length(obj.mc_results.predictions) && ...
                        isfield(obj.mc_results.predictions(d), 'y_true') && ...
                        ~isempty(obj.mc_results.predictions(d).y_true)

                    y_true = obj.mc_results.predictions(d).y_true;
                    y_pred = obj.mc_results.predictions(d).y_pred;
                    t_indices = obj.mc_results.predictions(d).t_indices;
                    % Use the prediction-resolution time axis (hold-aware); falls
                    % back to full-res t_test for older results without t_pred.
                    if isfield(obj.mc_results, 't_pred')
                        t_delay = obj.mc_results.t_pred(t_indices);
                    else
                        t_delay = t(t_indices);
                    end

                    % Plot true delayed input
                    plot(t_delay, y_true, 'k-', 'LineWidth', 0.8, 'DisplayName', sprintf('u(t-%d)', d));
                    hold on;
                    % Plot prediction
                    plot(t_delay, y_pred, '-', 'Color', colors_pred(i,:), 'LineWidth', 1.2, ...
                        'DisplayName', sprintf('\\hat{y}_{%d}(t)', d));
                    hold off;

                    R2 = obj.mc_results.R2_d(d);
                    title(sprintf('Delay d=%d: R^2=%.3f', d, R2), 'FontWeight', 'normal');
                    ylabel('Value');
                    legend('Location', 'best', 'FontSize', 7);
                else
                    text(0.5, 0.5, sprintf('No data for delay %d', d), ...
                        'HorizontalAlignment', 'center');
                    title(sprintf('Delay d=%d', d), 'FontWeight', 'normal');
                end

                if i < n_delays
                    set(gca, 'XTickLabel', []);
                else
                    xlabel('Time (s)');
                end
                grid on;
            end

            % Link all axes
            linkaxes(ax_handles, 'x');

            % Limit x-axis to test period
            xlim([t(1), t(end)]);

            % Add figure title
            if ~isempty(fig_title)
                sgtitle(fig_title);
            else
                sgtitle(sprintf('ESN Time Series (MC = %.2f)', obj.mc_results.MC));
            end
        end

        function reset(obj)
            % RESET Clear built state and memory capacity results

            reset@SRNNCellTypePairs(obj);
            obj.W_in = [];
            obj.u_scalar = [];
            obj.mc_results = [];
            vprintf(obj.verbose, 'verbose', 'ESN reservoir reset.\n');
        end
    end

    %% Protected Build Sub-Methods (override parent)
    methods (Access = protected)
        function build_stimulus(obj)
            % BUILD_STIMULUS Generate ESN-specific stimulus at build time
            %
            % Overrides SRNNCellTypePairs.build_stimulus() to generate:
            %   1. Input weight vector W_in
            %   2. Scalar input sequence u_scalar (white/bandlimited/1_over_f)
            %   3. Neural input matrix u_ex = W_in * u_scalar'
            %   4. Piecewise-constant griddedInterpolant for ODE solver
            %   5. Initial state vector S0

            % 1. Generate input weights W_in
            obj.generate_input_weights();

            % 2. Generate scalar input sequence
            T_total = obj.T_wash + obj.T_train + obj.T_test;
            rng(obj.rng_seeds(2));  % Use stimulus seed for reproducibility
            obj.hold_len = 1;       % reset; only 'sample_hold' overrides this

            if strcmpi(obj.input_type, 'bandlimited')
                % Generate zero-mean white noise for filtering
                u_raw = rand(T_total, 1) - 0.5;

                % Determine cutoff frequency
                if isempty(obj.u_f_cutoff)
                    f_cut = 1 / (2 * pi * obj.tau_d);
                else
                    f_cut = obj.u_f_cutoff;
                end

                % Design 3rd-order Butterworth low-pass filter
                f_nyq = obj.fs / 2;
                [b_filt, a_filt] = butter(3, f_cut / f_nyq, 'low');

                % Apply zero-phase filtering
                u_filtered = filtfilt(b_filt, a_filt, u_raw);

                % Normalize to [-0.5, 0.5] (zero-mean), then apply scaling and offset
                u_normalized = (u_filtered - min(u_filtered)) / (max(u_filtered) - min(u_filtered)) - 0.5;
                obj.u_scalar = obj.u_offset + obj.u_scale * u_normalized;

                vprintf(obj.verbose, 'verbose', 'ESN stimulus: bandlimited input (f_cutoff = %.2f Hz)\n', f_cut);

            elseif strcmpi(obj.input_type, 'one_over_f')
                % Generate 1/f^alpha noise using Fourier filtering method
                alpha_val = obj.u_alpha;

                N = T_total;
                X = randn(N, 1) + 1i * randn(N, 1);

                df = obj.fs / N;
                f = (0:(N-1))' * df;
                f(f > obj.fs/2) = f(f > obj.fs/2) - obj.fs;
                f = abs(f);
                f(f < df) = df;

                scale_factor = f.^(-alpha_val/2);
                X_shaped = X .* scale_factor;

                X_shaped(1) = real(X_shaped(1));
                if mod(N, 2) == 0
                    X_shaped(N/2 + 1) = real(X_shaped(N/2 + 1));
                end

                u_raw = real(ifft(X_shaped));
                u_normalized = (u_raw - min(u_raw)) / (max(u_raw) - min(u_raw)) - 0.5;
                obj.u_scalar = obj.u_offset + obj.u_scale * u_normalized;

                vprintf(obj.verbose, 'verbose', 'ESN stimulus: 1/f^%.2f noise input\n', alpha_val);

            elseif strcmpi(obj.input_type, 'white')
                obj.u_scalar = obj.u_offset + obj.u_scale * (rand(T_total, 1) - 0.5);
                obj.hold_len = 1;
                vprintf(obj.verbose, 'verbose', 'ESN stimulus: white noise input\n');

            elseif strcmpi(obj.input_type, 'sample_hold')
                % i.i.d. values held for hold_len samples (staircase). Puts input
                % energy at frequencies the network can follow and keeps values
                % i.i.d. at the hold rate (canonical MC in hold units).
                if isempty(obj.T_hold)
                    T_hold_eff = 3 * obj.tau_d;
                else
                    T_hold_eff = obj.T_hold;
                end
                L = max(1, round(T_hold_eff * obj.fs));
                obj.hold_len = L;
                n_holds = ceil(T_total / L);
                held_vals = obj.u_offset + obj.u_scale * (rand(n_holds, 1) - 0.5);
                u_full = repelem(held_vals, L);
                obj.u_scalar = u_full(1:T_total);
                vprintf(obj.verbose, 'verbose', 'ESN stimulus: sample-and-hold input (T_hold = %.3f s, %d samples/hold)\n', ...
                    T_hold_eff, L);
            else
                error('SRNN_ESN_reservoir:InvalidInputType', ...
                    'Unknown input_type ''%s''. Valid options: ''white'', ''bandlimited'', ''one_over_f'', ''sample_hold''', ...
                    obj.input_type);
            end

            % 3. Map scalar input to neural input
            dt = 1 / obj.fs;
            obj.t_ex = (0:(T_total-1))' * dt;
            obj.u_ex = obj.W_in * obj.u_scalar';  % n x T

            % 4. Create linear interpolant for ODE solver
            obj.u_interpolant = griddedInterpolant(obj.t_ex, obj.u_ex', 'linear', 'none');

            % 5. Initialize state vector.
            % Called through the class, not the instance: initialize_state is a
            % public STATIC on SRNNCellTypePairs where it was a protected
            % INSTANCE method on SRNNModel2. obj.initialize_state(params) happens
            % to resolve either way in MATLAB, which is exactly why it is worth
            % writing explicitly rather than leaving to that.
            params_init = obj.get_params();
            obj.S0 = SRNNCellTypePairs.initialize_state(params_init);

            vprintf(obj.verbose, 'verbose', 'ESN stimulus built: %d samples, %d neurons receive input\n', ...
                T_total, sum(obj.W_in ~= 0));
        end
    end

    %% Static Methods for Ridge Regression and R^2
    methods (Access = protected)
        function plot_types(~, t, named, test_indices, names, tcolors, order)
            % PLOT_TYPES Draw one per-cell-type signal, correctly layered.
            %
            % `named` is a struct keyed by cell-type name, each field nq x T (or
            % nq x numel(test_indices) when the caller has already sliced, in
            % which case pass test_indices = []). An empty field is skipped, so
            % a type with no adaptation or no depression simply does not draw.
            %
            % Types are drawn in `order` (n:-1:1), so type 1 lands on top --
            % the same back-to-front convention SRNNCellTypePairs.plot_celltypes
            % uses wherever types share an axes. Colours stay indexed BY TYPE,
            % not by draw position, so the reversal is layering only.
            %
            % Per-neuron accents come from neuron_colors, which varies lightness
            % within the type's hue. The previous E/I version used two fixed
            % colormaps, which cannot generalize past two types.
            drawn = false;
            for k = order
                data = named.(names{k});
                if isempty(data); continue; end
                if ~isempty(test_indices)
                    data = data(:, test_indices);
                end
                if isempty(data); continue; end
                if drawn; hold on; end
                cmap = SRNNCellTypePairs.neuron_colors(tcolors(k, :), size(data, 1));
                plot_lines_with_colormap(t, data, cmap);
                drawn = true;
            end
            if drawn; hold off; end
        end
    end

    methods (Static)
        function R = readout_from_states(r_all, syn_all, params, readout_signal)
            % READOUT_FROM_STATES The n x T signal the linear readout sees.
            %
            % r_all and syn_all are unpack_and_compute_states' rate and
            % synaptic_output for any set of time samples. Every cell type's
            % signal is stacked in TYPE ORDER, which is the order type_indices
            % uses -- so row i is neuron i, and the readout weights stay
            % interpretable. One definition for both memory-capacity paths.
            %
            % 'rate'     => firing rate r = phi(x_eff) (adaptation baked into
            %               x_eff, but depression NOT included).
            % 'synaptic' => r times the route's output factor, the
            %               depressed/facilitated presynaptic output, which
            %               exposes the STD state to the linear readout. The
            %               factor is the route MODEL's (route_model): on a
            %               product route prod(b) * prod(g), on a chain
            %               route b1, the releasable pool. It is not formed
            %               here: unpack_and_compute_states takes it from
            %               the model's `output`, so syn_all is right for
            %               every route model. Without STD it equals r, so
            %               those conditions are unchanged.
            switch lower(readout_signal)
                case 'rate'
                    R = SRNN_ESN_reservoir.stack_by_type(r_all, params);
                case 'synaptic'
                    % Per-route on this class: a presynaptic neuron has one
                    % synaptic output PER TARGET TYPE. They coincide only when
                    % all of that neuron's outgoing routes carry identical
                    % synaptic dynamics -- which this refuses to assume.
                    SRNN_ESN_reservoir.assert_route_redundancy(params);
                    first_route = struct();
                    for q = 1:params.n_cellTypes
                        name = params.cell_type_names{q};
                        posts = fieldnames(syn_all.(name));
                        first_route.(name) = syn_all.(name).(posts{1});
                    end
                    R = SRNN_ESN_reservoir.stack_by_type(first_route, params);
                otherwise
                    error('SRNN_ESN_reservoir:BadReadoutSignal', ...
                        'readout_signal must be ''rate'' or ''synaptic'' (got ''%s'').', readout_signal);
            end
        end

        function R = stack_by_type(named, params)
            % STACK_BY_TYPE Concatenate a per-cell-type signal into one n x T.
            %
            % Rows come out in TYPE ORDER, which is the order type_indices
            % assigns neurons, so row i of the result is neuron i and a trained
            % readout weight can still be attributed to a neuron.
            parts = cell(1, params.n_cellTypes);
            for q = 1:params.n_cellTypes
                parts{q} = named.(params.cell_type_names{q});
            end
            R = vertcat(parts{:});
        end

        function assert_route_redundancy(params)
            % ASSERT_ROUTE_REDUNDANCY The 'synaptic' readout needs each neuron to
            % have ONE synaptic output. On this class it has one per target type.
            %
            % theta_j^{q->s} = r_j * F^{q->s}, with F the route's output factor
            % (prod(b) * prod(g) on a product route, b1 on a chain route), so a
            % presynaptic neuron of type q transmits a different signal down each
            % outgoing route unless those routes share their synaptic dynamics.
            % When they do share them the coincidence is EXACT, not approximate:
            % a route's state obeys an ODE that depends only on the route's
            % parameters and the presynaptic neuron's own rate, from a common
            % initial state -- identical ODE, identical driver, identical
            % initial condition, so bit-identical trajectories. Reading route 1
            % is then not a shortcut but the same array.
            %
            % PER PRESYNAPTIC TYPE, not globally. Neuron j of type q only has
            % routes q->1 ... q->C; what other types do cannot affect it. A
            % global check would reject a perfectly well-defined configuration
            % such as "E->* depressed, I->* not", or "E->* a chain, I->* a
            % product".
            %
            % WHAT IS COMPARED is SRNNCellTypePairs.route_fields, the one
            % definition of "the same route" (routes_identical uses it too):
            % the legacy STD/STF rows, the route scale, and the route MODEL's
            % own signature -- which model, and every parameter of it. Until
            % 2026-11-12 this class kept a private list (n_b, tau_rec, tau_rel,
            % tau_dec, tau_fac, G, scale, std_geomean), written when depression
            % was always a product. It could not see a chain's n2 at all, so
            % two chain routes differing only in n2 PASSED and the readout
            % silently read one of two different signals; and it told a chain
            % from a two-component product only by the accident that the chain
            % stores one tau_rel where the product stores two.
            %
            % Checked from the CONFIG rather than by comparing trajectories: the
            % config determines the answer exactly, and comparing npre x nt
            % arrays per route on every trial would be wasted work. The numerical
            % claim is asserted once, in test_esn_route_redundancy.
            for pre = 1:params.n_cellTypes
                pre_name = params.cell_type_names{pre};
                ref = SRNNCellTypePairs.route_fields(params, pre, 1);
                for post = 2:params.n_cellTypes
                    here = SRNNCellTypePairs.route_fields(params, pre, post);
                    if ~isequal(ref, here)
                        error('SRNN_ESN_reservoir:AmbiguousSynapticOutput', ...
                            ['Presynaptic type %s has non-identical synaptic ' ...
                             'dynamics across its outgoing routes:\n' ...
                             '  %s->%s : %s\n  %s->%s : %s\n' ...
                             'They differ in:\n%s\n' ...
                             'A neuron of type %s therefore transmits a ' ...
                             'different signal to each target, so the ' ...
                             '''synaptic'' readout has no single value to ' ...
                             'read.\nUse readout_signal = ''rate'', or give ' ...
                             'every %s->* route the same std/stf settings.'], ...
                            pre_name, ...
                            pre_name, params.cell_type_names{1}, ...
                            SRNN_ESN_reservoir.describe_route(ref), ...
                            pre_name, params.cell_type_names{post}, ...
                            SRNN_ESN_reservoir.describe_route(here), ...
                            SRNN_ESN_reservoir.describe_route_difference(ref, here), ...
                            pre_name, pre_name);
                    end
                end
            end
        end

        function sig = route_signature(params, pre, post)
            % ROUTE_SIGNATURE The synaptic dynamics of one route, as comparable
            % data: SRNNCellTypePairs.route_fields, under the name this class
            % has always had for it. There is ONE definition and it is the
            % model class's; this is kept only so the name still resolves. It
            % covers STD and STF, the route scale (a weight, not part of theta,
            % but two routes of one presynaptic type scaled differently would
            % make "the synaptic output of neuron j" ambiguous again) and the
            % route model's signature (the model, and what only it has: a
            % chain's n2).
            sig = SRNNCellTypePairs.route_fields(params, pre, post);
        end

        function s = describe_route(sig)
            % DESCRIBE_ROUTE One-line human description, for the error above.
            %
            % sig is a route_fields struct. Numbers are printed with as many
            % digits as it takes to read them back exactly (number_text), so
            % two routes that differ at all never print alike: the rounded
            % three digits this used before showed two identical lines under
            % "non-identical" whenever the difference was small.
            model = '';
            if isfield(sig, 'signature') && isfield(sig.signature, 'model')
                model = sig.signature.model;
            end
            num = @SRNN_ESN_reservoir.number_text;
            parts = {};
            coupled = false;
            if strcmp(model, 'pool')
                % A vesicle-pool route (route_pool). Its depression is the
                % chain, a single pool or absent, and its facilitation, if
                % any, release-coupled; the model's signature says which
                % (std_model, stf_model). Until 2026-11-12 every pool route
                % was described as a chain, reading tau_rec(2) and n2,
                % which a single pool and a route without depression do
                % not have: two such routes that differed made the
                % redundancy check throw MATLAB:badsubscript from here
                % where it should have raised AmbiguousSynapticOutput.
                std_model = 'none';
                if isfield(sig.signature, 'std_model'), std_model = sig.signature.std_model; end
                coupled = isfield(sig.signature, 'stf_model') && strcmp(sig.signature.stf_model, 'chain');
                if strcmp(std_model, 'chain')
                    % tau_rec holds [releasable-pool refill, reserve
                    % refill] and tau_rel ONE value, the releasable pool's.
                    parts{end+1} = sprintf(['chained STD(releasable-pool refill tau_rec1=%s, ' ...
                        'reserve refill tau_rec2=%s, release tau_rel=%s, n2=%s)'], ...
                        num(sig.tau_rec(1)), num(sig.tau_rec(2)), num(sig.tau_rel), ...
                        num(sig.signature.n2));
                elseif sig.n_b > 0
                    parts{end+1} = sprintf('single-pool STD(refill tau_rec=%s, release tau_rel=%s)', ...
                        num(sig.tau_rec), num(sig.tau_rel));
                end
            elseif sig.n_b > 0
                extra = '';
                if sig.std_geomean
                    extra = [extra ', combine=geomean'];
                end
                if isfield(sig, 'signature') && isfield(sig.signature, 'zero_floor') ...
                        && sig.signature.zero_floor
                    extra = [extra ', zero floor'];
                end
                parts{end+1} = sprintf('product STD(tau_rec=%s, tau_rel=%s%s)', ...
                    num(sig.tau_rec), num(sig.tau_rel), extra);
            end
            if coupled
                % Release-coupled: one g, and a calcium trace where the
                % route has one (n_g is then 2 and tau_h is stated).
                trace = '';
                if isfield(sig.signature, 'tau_h') && ~isempty(sig.signature.tau_h)
                    trace = sprintf(', calcium trace tau_h=%s', num(sig.signature.tau_h));
                end
                parts{end+1} = sprintf('release-coupled STF(tau_dec=%s, tau_fac=%s, G=%s%s)', ...
                    num(sig.tau_dec), num(sig.tau_fac), num(sig.G), trace);
            elseif sig.n_g > 0
                parts{end+1} = sprintf('STF(tau_dec=%s, tau_fac=%s, G=%s)', ...
                    num(sig.tau_dec), num(sig.tau_fac), num(sig.G));
            end
            if isempty(parts)
                s = 'no short-term plasticity';
            else
                s = strjoin(parts, ' + ');
            end
            if isfield(sig, 'scale') && sig.scale ~= 1
                s = sprintf('%s, scale=%s', s, num(sig.scale));
            end
        end

        function s = describe_route_difference(ref, here)
            % DESCRIBE_ROUTE_DIFFERENCE Name every field in which two routes differ.
            %
            % ref and here are route_fields structs. One line per differing
            % field, both values in full, and the relative difference where
            % both are numbers of one size -- so a difference of 1e-6 reads
            % as one. The fields are the route model's signature (model,
            % time constants, a chain's n2, a product's combine and floor)
            % and whatever route_fields carries beside it (the counts, the
            % scale); a field one route does not have prints as (absent).
            % Written generically over the two structs' fields, so a field
            % added to a model's signature is reported with no change here.
            [names, a, b] = SRNN_ESN_reservoir.comparable_fields(ref, here);
            lines = {};
            for k = 1:numel(names)
                if isequal(a{k}, b{k})
                    continue
                end
                line = sprintf('  %s: %s vs %s', names{k}, ...
                    SRNN_ESN_reservoir.value_text(a{k}), ...
                    SRNN_ESN_reservoir.value_text(b{k}));
                if isnumeric(a{k}) && isnumeric(b{k}) && ~isempty(a{k}) ...
                        && isequal(size(a{k}), size(b{k}))
                    va = double(a{k}(:)); vb = double(b{k}(:));
                    rel = max(abs(va - vb) ./ max(abs(va), abs(vb)));
                    line = sprintf('%s (relative difference %.3g)', line, rel);
                end
                lines{end+1} = line; %#ok<AGROW>
            end
            if isempty(lines)
                s = '  (no single field; the two routes'' data are not isequal)';
            else
                s = strjoin(lines, '\n');
            end
        end

        function [names, a, b] = comparable_fields(ref, here)
            % The fields of two route_fields structs, side by side. The route
            % model's signature is listed field by field; a top-level field
            % the signature repeats under the same name (tau_rec, ...) is
            % listed once, from the signature. A missing field is the char
            % '(absent)'.
            absent = '(absent)';
            sig_a = struct(); sig_b = struct();
            if isfield(ref, 'signature'),  sig_a = ref.signature;  end
            if isfield(here, 'signature'), sig_b = here.signature; end
            sig_names = union(fieldnames(sig_a), fieldnames(sig_b), 'stable');
            top_names = setdiff(union(fieldnames(ref), fieldnames(here), 'stable'), ...
                [{'signature'}; sig_names(:)], 'stable');
            names = [sig_names(:); top_names(:)]';
            a = cell(size(names)); b = cell(size(names));
            for k = 1:numel(names)
                if k <= numel(sig_names)
                    sa = sig_a; sb = sig_b;
                else
                    sa = ref; sb = here;
                end
                if isfield(sa, names{k}), a{k} = sa.(names{k}); else, a{k} = absent; end
                if isfield(sb, names{k}), b{k} = sb.(names{k}); else, b{k} = absent; end
            end
        end

        function s = value_text(v)
            % One field's value, for an error message.
            if ischar(v)
                if strcmp(v, '(absent)'), s = v; else, s = ['''' v '''']; end
            elseif islogical(v) && isscalar(v)
                if v, s = 'true'; else, s = 'false'; end
            elseif isnumeric(v) || islogical(v)
                s = SRNN_ESN_reservoir.number_text(double(v));
            else
                s = strtrim(evalc('disp(v)'));
            end
        end

        function s = number_text(v)
            % NUMBER_TEXT A numeric array as text that reads back EXACTLY.
            %
            % Each element is printed with the fewest significant digits (up
            % to 17) whose text converts back to the same double, so 0.25 is
            % '0.25' and 2*(1 + 1e-6) is '2.000002', never '2'. Two values
            % that differ therefore never print alike.
            if isempty(v)
                s = '[]';
                return
            end
            parts = cell(1, numel(v));
            for k = 1:numel(v)
                for digits = 1:17
                    parts{k} = sprintf('%.*g', digits, v(k));
                    if str2double(parts{k}) == v(k) || isnan(v(k))
                        break
                    end
                end
            end
            s = strjoin(parts, ' ');
            if numel(v) > 1
                s = ['[' s ']'];
            end
        end

        function w = train_linear_readout(X, y, eta)
            % TRAIN_LINEAR_READOUT Train a linear readout using ridge regression
            %
            % w = train_linear_readout(X, y, eta)
            %
            % Inputs:
            %   X   - Design matrix (T x n_features)
            %   y   - Target vector (T x 1)
            %   eta - Regularization parameter
            %
            % Outputs:
            %   w   - Weight vector (n_features x 1)
            %
            % Closed-form solution:
            %   w = (X'X + eta*I)^(-1) * X'y

            n_features = size(X, 2);
            XTX = X' * X;
            XTy = X' * y;
            w = (XTX + eta * eye(n_features)) \ XTy;
        end

        function R2 = compute_R2(y_true, y_pred)
            % COMPUTE_R2 Compute coefficient of determination
            %
            % R2 = compute_R2(y_true, y_pred)
            %
            % Computes R^2 as per the memory capacity protocol:
            %   R^2 = cov(y_true, y_pred)^2 / (var(y_true) * var(y_pred))
            %
            % This is the squared correlation coefficient.

            % Ensure column vectors
            y_true = y_true(:);
            y_pred = y_pred(:);

            % Handle degenerate cases
            var_true = var(y_true);
            var_pred = var(y_pred);

            if var_true < 1e-12 || var_pred < 1e-12
                R2 = 0;
                return;
            end

            % Compute squared correlation
            cov_matrix = cov(y_true, y_pred);
            cov_val = cov_matrix(1, 2);

            R2 = (cov_val^2) / (var_true * var_pred);

            % Clamp to [0, 1]
            R2 = max(0, min(1, R2));
        end

        function verify_shared_build(esn_array, expected_to_differ, also_check_protected)
            % VERIFY_SHARED_BUILD Verify that built ESN objects share identical configuration
            % Internalized from verify_shared_build
            %
            % SRNN_ESN_reservoir.verify_shared_build(esn_array, expected_to_differ, also_check_protected)
            %
            % After building multiple ESN objects for a comparison experiment, this
            % function verifies that all properties that SHOULD be shared are identical,
            % and that properties expected to differ DO actually differ.

            if numel(esn_array) < 2
                vprintf(esn_array{1}.verbose, 'verbose', 'verify_shared_build: only 1 object, nothing to compare.\n');
                return;
            end

            ref = esn_array{1};
            mc = metaclass(ref);
            n_obj = numel(esn_array);

            always_skip = {'S0', 'cached_params', 'mc_results', 'u_interpolant', ...
                           'ode_opts', 't_out', 'S_out', 'plot_data', 'lya_results', ...
                           'noise_increments', 'verbose', ...
                           'tau_a_matrix', 'tau_a_clipped'};   % derived from tau_a (2026-09-14)

            n_checked = 0;
            n_matched = 0;
            n_skipped = 0;
            checked_names = {};
            skipped_names = {};

            for p = 1:numel(mc.PropertyList)
                prop = mc.PropertyList(p);
                name = prop.Name;

                if prop.Dependent
                    skipped_names{end+1} = name; %#ok<AGROW>
                    n_skipped = n_skipped + 1;
                    continue;
                end

                if ismember(name, always_skip)
                    skipped_names{end+1} = name; %#ok<AGROW>
                    n_skipped = n_skipped + 1;
                    continue;
                end

                if ismember(name, expected_to_differ)
                    skipped_names{end+1} = name; %#ok<AGROW>
                    n_skipped = n_skipped + 1;
                    continue;
                end

                is_public_get = strcmp(prop.GetAccess, 'public');
                is_in_also_check = ismember(name, also_check_protected);

                if ~is_public_get && ~is_in_also_check
                    skipped_names{end+1} = name; %#ok<AGROW>
                    n_skipped = n_skipped + 1;
                    continue;
                end

                for i = 2:n_obj
                    obj = esn_array{i};
                    val_ref = ref.(name);
                    val_obj = obj.(name);

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

            for k = 1:numel(expected_to_differ)
                name = expected_to_differ{k};
                prop_meta = findobj(mc.PropertyList, 'Name', name);
                if isempty(prop_meta)
                    warning('verify_shared_build:UnknownProperty', ...
                        'Property ''%s'' in expected_to_differ does not exist on this class.', name);
                    continue;
                end

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

            vprintf(ref.verbose, 'verbose', 'verify_shared_build: %d properties checked, all matched across %d conditions.\n', ...
                n_checked, n_obj);
            vprintf(ref.verbose, 'verbose', '  Checked: %s\n', strjoin(checked_names, ', '));
            vprintf(ref.verbose, 'verbose', '  Expected to differ: %s\n', strjoin(expected_to_differ, ', '));
            if ~isempty(also_check_protected)
                vprintf(ref.verbose, 'verbose', '  Also verified (protected): %s\n', strjoin(also_check_protected, ', '));
            end
        end
    end
end
