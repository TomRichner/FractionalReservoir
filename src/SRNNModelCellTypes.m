classdef SRNNModelCellTypes < SRNNModelBase
    % SRNNMODELCELLTYPES  Cell-type-resolved rate reservoir (Pyr/Pvalb/Sst/Vip).
    %
    % A concrete SRNNModelBase subclass (sibling of SRNNModel2) that models K cell
    % types with a PER-NEURON type-index representation and PER-EDGE (pre-type x
    % post-type) short-term plasticity, parameterized from the Allen/Campagnola 2022
    % matrices (via load_campagnola_matrices).
    %
    % Mechanisms:
    %   - SFA  (a): per neuron, intrinsic firing adaptation. State a is n x n_a.
    %               x_eff = x - c .* sum(a,2);  da/dt = (r - a)./tau_a.
    %   - STD  (b): per (presynaptic neuron j, post-type T). State b is n x K (n_b=1).
    %               db/dt = (1-b)./tau_rec - (p.*b.*r)./tau_rel.  (driven by r_j; depletion coupled to p)
    %   - STF  (p): per (presynaptic neuron j, post-type T). State p is n x K (n_u=1).
    %               p is a dynamic release probability (Tsodyks-Markram u, renamed to p
    %               because u denotes the external input here).
    %               dp/dt = (p0-p)./tau_f + kappa.*(1-p).*r.  (NEW; driven by r_j)
    % Synaptic efficacy of j -> post-type T is eff(j,T) = (p(j,T)/p0(j,T)) * b(j,T),
    % so STF-off (p=p0) reduces to pure STD and rest (p=p0,b=1) gives eff=1.
    %
    % The recurrent drive is post-type structured: neuron i (type T_i) receives
    %   drive_i = sum_j W_ij * eff(j,T_i) * r_j.
    % Using the full (pre,post) data-matrix element param(type(j),T) automatically
    % realizes the paper's alignment rules (excitatory dynamics vary with post-type,
    % inhibitory mostly with pre-type).
    %
    % State layout  S = [a(:); b(:); p(:); x(:)]   (x last), with
    %   a -> (n, n_a),  b -> (n, K),  p -> (n, K),  x -> (n,1).
    %   N_sys_eqs = n*n_a + n*K*n_b + n*K*n_u + n.
    %
    % Scope of the scaffold: SFA supports multiple timescales (n_a>=0); STD/STF are
    % single-timescale (n_b,n_u in {0,1}). The optional STD zero-floor is not applied.
    %
    % Lifecycle:  model = SRNNModelCellTypes(...); model.build(); model.run(); model.plot();
    %
    % See also: SRNNModelBase, SRNNModel2, load_campagnola_matrices

    %% Cell-type configuration
    properties
        type_names = {'pyr', 'pvalb', 'sst', 'vip'}   % K types; type 1 is excitatory (Pyr)
        type_fractions = [0.80, 0.08, 0.07, 0.05]     % fraction of neurons per type (renormalized)
        use_campagnola_data = true                    % parameterize from load_campagnola_matrices
    end

    %% Mechanism configuration
    properties
        n_a = 1                  % # SFA timescales for ADAPTING types (0 disables SFA everywhere)
        n_b = 1                  % # STD timescales (0 or 1)
        n_u = 1                  % # STF timescales (0 or 1)
        tau_a = 1.0              % fallback SFA time constant(s), 1 x n_a (s), if no per-type fit
        c_gain = 0.7             % maps per-type adaptation_index -> SFA strength c
        w_cv = 1.0               % per-edge weight heterogeneity (std / |mean|)
        tau_b_rel_ref = 0.5      % reference STD release tau (s); scaled down by depression amount
        sfa_min_index = 0.01     % a type with adaptation_index below this is non-adapting (n_a=0)
    end

    %% Per-type parameter tables (K x K unless noted; filled at build if empty)
    properties
        conn_prob      % K x K  connection probability (pre -> post)
        psp_amp        % K x K  signed synaptic weight (pre -> post)
        dep_tau        % K x K  STD recovery tau (s)
        dep_amount     % K x K  STD depression amount [0,1)
        rel_prob       % K x K  baseline release probability p0 (STF rest value)
        fac_tau        % K x K  STF facilitation tau (s)
        kappa          % K x K  STF facilitation rate kappa (from ml_facilitation_amount)
        adapt_index    % K x 1  per-type adaptation index (SFA strength source)
        tau_a_type     % K x 1  per-type fitted SFA tau (s); NaN entries fall back to tau_a
        tau_d_type     % K x 1  per-type membrane tau (s) from Campagnola C.tau; settable override
    end

    %% Computed (set at build)
    properties (SetAccess = protected)
        type_of        % n x 1  integer type label per neuron
        is_exc         % n x 1  logical, true for excitatory (Pyr) neurons
        n_types        % number of cell types (K)
        adapting       % n x 1  logical, true for neurons of adapting types (carry SFA state)
        ad_idx         % n_ad x 1  global indices of adapting neurons
        n_ad           % number of adapting neurons (SFA state is n_ad x n_a)
    end

    %% Constructor
    methods
        function obj = SRNNModelCellTypes(varargin)
            % Forwards to the shared base constructor (defaults + name-value parse).
            obj@SRNNModelBase(varargin{:});
        end
    end

    %% Subclass hooks
    methods (Access = protected)
        function set_defaults(obj)
            % SET_DEFAULTS Activation, input_config, and plotting defaults.
            obj.activation_function            = @(x) SRNNModelBase.logisticSigmoid(x, obj.S_c);
            obj.activation_function_derivative = @(x) SRNNModelBase.logisticSigmoidDerivative(x, obj.S_c);

            % Cell-typed W uses signed psp_amplitude weights (tiny, O(1e-4)), so
            % normalize the spectral abscissa by default: level_of_chaos then sets
            % the edge-of-chaos scale directly (overridable via name-value).
            obj.rescale_by_abscissa = true;

            obj.input_config = struct();
            obj.input_config.n_steps        = 3;
            obj.input_config.amp            = 0.5;
            obj.input_config.step_density   = 0.15;   % fraction of neurons driven per step
            obj.input_config.no_stim_pattern = false(1, 3);
            obj.input_config.no_stim_pattern(1:2:end) = true;
            obj.input_config.intrinsic_drive = [];
            obj.input_config.positive_only   = false;

            obj.T_plot = [];
        end

        function build_network(obj)
            % BUILD_NETWORK Assign types, gather per-type parameters, build block W.
            rng(obj.rng_seeds(1));
            obj.n_types = numel(obj.type_names);
            obj.assign_types();
            obj.load_parameter_tables();          % fills the K x K tables (data or defaults)

            % Per-type SFA: a type is adapting iff its adaptation_index >= sfa_min_index (and
            % SFA is enabled globally, n_a>0). Only adapting neurons carry SFA state (n_a=0 for
            % the rest, e.g. fast-spiking Pvalb) -> the "adapting-neuron subset" ragged layout.
            type_adapts = (obj.adapt_index(:) >= obj.sfa_min_index);
            obj.adapting = (obj.n_a > 0) & type_adapts(obj.type_of);
            obj.ad_idx = find(obj.adapting);
            obj.n_ad = numel(obj.ad_idx);

            n = obj.n; K = obj.n_types; t = obj.type_of;

            % Block-structured signed connectivity. Column j (presynaptic) carries the
            % sign of its pre-type; each (pre,post) block has its own p and weight.
            W = zeros(n, n);
            for tp = 1:K                          % presynaptic type (columns)
                cols = find(t == tp);
                for ts = 1:K                      % postsynaptic type (rows)
                    rows = find(t == ts);
                    if isempty(rows) || isempty(cols), continue; end
                    p  = min(max(obj.conn_prob(tp, ts), 0), 1);
                    mu = obj.psp_amp(tp, ts);
                    sig = obj.w_cv * abs(mu);
                    m = numel(rows); q = numel(cols);
                    blk = (mu + sig * randn(m, q)) .* (rand(m, q) < p);
                    W(rows, cols) = blk;
                end
            end

            % Scale W (same edge-of-chaos convention as SRNNModel2).
            if obj.rescale_by_abscissa
                abscissa_0 = max(real(eig(W)));
                if abscissa_0 <= 0, abscissa_0 = max(abs(eig(W))); end
                if abscissa_0 == 0, abscissa_0 = 1; end
                obj.W = obj.level_of_chaos * (1 / abscissa_0) * W;
            else
                obj.W = obj.level_of_chaos * W;
            end

            W_eigs = eig(obj.W);
            fprintf('CellTypes W: n=%d, K=%d, spectral radius=%.3f, abscissa=%.3f\n', ...
                n, K, max(abs(W_eigs)), max(real(W_eigs)));
        end

        function build_stimulus(obj)
            % BUILD_STIMULUS Sparse step external input, interpolant, and S0.
            dt = 1 / obj.fs;
            T  = obj.T_range(2);
            t_stim = (0:dt:T)';
            nt = numel(t_stim);

            rng(obj.rng_seeds(2));
            ic = obj.input_config;
            n = obj.n;
            step_period = fix(T / ic.n_steps);
            step_len = round(step_period * obj.fs);
            if ic.positive_only
                steps = ic.amp * abs(randn(n, ic.n_steps));
            else
                steps = ic.amp * randn(n, ic.n_steps);
            end
            mask = rand(n, ic.n_steps) < ic.step_density;
            steps = steps .* mask;
            steps(:, ic.no_stim_pattern) = 0;

            u_stim = zeros(n, nt);
            for s = 1:ic.n_steps
                a0 = (s - 1) * step_len + 1;
                b0 = min(s * step_len, nt);
                if a0 > nt, break; end
                u_stim(:, a0:b0) = repmat(steps(:, s), 1, b0 - a0 + 1);
            end
            if ~isempty(ic.intrinsic_drive)
                u_stim = u_stim + ic.intrinsic_drive;
            end

            % Prepend zeros for a negative settling window (as SRNNModel2 does).
            if obj.T_range(1) < 0
                t_pre = (obj.T_range(1):dt:-dt)';
                obj.t_ex = [t_pre; t_stim];
                obj.u_ex = [zeros(n, numel(t_pre)), u_stim];
            else
                keep = t_stim >= obj.T_range(1);
                obj.t_ex = t_stim(keep);
                obj.u_ex = u_stim(:, keep);
            end
            obj.u_ex = obj.u_ex * obj.u_ex_scale;

            obj.u_interpolant = griddedInterpolant(obj.t_ex, obj.u_ex', 'linear', 'none');
            obj.S0 = obj.initialize_state(obj.get_params());
            fprintf('CellTypes stimulus: %d time points, %d neurons\n', numel(obj.t_ex), n);
        end

        function validate(obj)
            % VALIDATE Parameter-consistency checks.
            if obj.n < obj.n_types
                error('SRNNModelCellTypes:TooFewNeurons', 'n (%d) must be >= K (%d).', obj.n, obj.n_types);
            end
            if ~ismember(obj.n_b, [0 1]) || ~ismember(obj.n_u, [0 1])
                error('SRNNModelCellTypes:BadTimescales', 'n_b and n_u must be 0 or 1 (single-timescale STD/STF).');
            end
            if obj.n_a < 0
                error('SRNNModelCellTypes:BadTimescales', 'n_a must be >= 0.');
            end
            if obj.T_range(2) <= obj.T_range(1)
                error('SRNNModelCellTypes:BadT', 'T_range(2) must be > T_range(1).');
            end
        end

        function dS_dt = eval_dynamics(~, t, S, params)
            dS_dt = SRNNModelCellTypes.dynamics_fast_ct(t, S, params);
        end

        function J = eval_jacobian(~, S, params)
            J = SRNNModelCellTypes.compute_Jacobian_fast_ct(S, params);
        end

        function decimate_and_unpack(obj)
            % DECIMATE_AND_UNPACK Decimate the trajectory and store plot_data.
            params = obj.cached_params;
            [t_plot, S_plot, plot_idx] = obj.decimate_states(obj.t_out, obj.S_out, obj.plot_deci);
            st = SRNNModelCellTypes.unpack_states_ct(S_plot, params);   % fields x,a,b,u,r,br (n x nt)

            Tp = obj.T_plot; if isempty(Tp), Tp = obj.T_range; end
            keep = t_plot >= Tp(1) & t_plot <= Tp(2);

            pd = struct();
            pd.t       = t_plot(keep);
            pd.u_ext   = obj.u_ex(:, plot_idx); pd.u_ext = pd.u_ext(:, keep);
            pd.x       = st.x(:, keep);
            pd.r       = st.r(:, keep);
            pd.br      = st.br(:, keep);
            pd.b       = st.b(:, :, keep);          % n x K x nt (per post-type)
            pd.p       = st.p(:, :, keep);          % n x K x nt  (STF release-prob state)
            if ~isempty(st.a), pd.a = st.a(:, :, keep); else, pd.a = []; end
            pd.type_of = obj.type_of;
            pd.type_names = obj.type_names;
            obj.plot_data = pd;
        end

        function S0 = initialize_state(~, params)
            % INITIALIZE_STATE  a=0, b=1 (no depression), p=p0 (facilitation at rest), x random.
            n = params.n; K = params.K;
            a0 = zeros(params.n_ad * params.n_a, 1);       % SFA state only for adapting neurons
            if params.n_b > 0, b0 = ones(n * K, 1); else, b0 = []; end
            if params.n_u > 0, u0 = params.p0_mat(:); else, u0 = []; end
            x0 = params.x0_std .* randn(n, 1);
            S0 = [a0; b0; u0; x0];
        end
    end

    %% Public
    methods
        function params = get_params(obj)
            % GET_PARAMS Pack everything the kernels / base machinery read.
            params = struct();
            params.n = obj.n;
            params.K = obj.n_types;
            params.n_a = obj.n_a; params.n_b = obj.n_b; params.n_u = obj.n_u;
            % SFA state is ragged: only the n_ad adapting neurons carry a (n_ad x n_a).
            params.n_ad = obj.n_ad; params.ad_idx = obj.ad_idx; params.adapting = obj.adapting;
            params.N_sys_eqs = obj.n_ad * obj.n_a + obj.n * obj.n_types * obj.n_b + ...
                               obj.n * obj.n_types * obj.n_u + obj.n;
            % Dendritic time constant, PER NEURON (n x 1): each neuron uses its type's membrane tau.
            if ~isempty(obj.tau_d_type) && ~isempty(obj.type_of)
                params.tau_d = obj.tau_d_type(obj.type_of);      % n x 1
            else
                params.tau_d = repmat(obj.tau_d, obj.n, 1);      % scalar fallback (data disabled)
            end
            params.activation_function = obj.activation_function;
            params.activation_function_derivative = obj.activation_function_derivative;
            params.x0_std = obj.x0_std;
            params.rng_seeds = obj.rng_seeds;
            params.type_of = obj.type_of;
            params.is_exc = obj.is_exc;

            % SFA time constant, PER NEURON (n x n_a): each neuron uses its type's fitted tau.
            if obj.n_a > 0 && ~isempty(obj.type_of) && ~isempty(obj.tau_a_type)
                params.tau_a = repmat(obj.tau_a_type(obj.type_of), 1, obj.n_a);   % n x n_a
            else
                ta = obj.tau_a; if isscalar(ta), ta = repmat(ta, 1, max(obj.n_a, 1)); end
                params.tau_a = repmat(reshape(ta, 1, []), obj.n, 1);              % n x n_a fallback
            end
            if ~isempty(obj.adapt_index) && ~isempty(obj.type_of)
                params.c_vec = obj.c_gain * obj.adapt_index(obj.type_of);   % n x 1
            else
                params.c_vec = zeros(obj.n, 1);
            end

            % STD / STF per (neuron j, post-type T): gather K x K -> n x K by pre-type row.
            if ~isempty(obj.type_of)
                t = obj.type_of;
                params.tau_b_rec_mat = obj.dep_tau(t, :);
                params.tau_b_rel_mat = max(obj.tau_b_rel_ref * (1 - min(max(obj.dep_amount(t, :), 0), 0.95)), 0.05);
                params.p0_mat        = min(max(obj.rel_prob(t, :), 0.05), 0.95);
                params.tau_f_mat     = obj.fac_tau(t, :);
                params.kappa_mat     = max(obj.kappa(t, :), 0);
            end

            if ~isempty(obj.W), params.W = obj.W; end
        end

        function [fig, ax] = plot(obj)
            % PLOT Stacked per-type time series (all colored by cell type): external
            % stimulus, dendritic x, firing rate r, synaptic output, SFA, STD, STF, and
            % the local Lyapunov exponent. Each panel shows every neuron's trace (faint,
            % in its type's color) plus a bold per-type mean; panels for a mechanism are
            % included only when it is enabled.
            if ~obj.has_run || isempty(obj.plot_data)
                error('SRNNModelCellTypes:NotRun', 'Run the model before plotting.');
            end
            pd = obj.plot_data; t = pd.t; ntype = pd.type_of;
            n = size(pd.x, 1);
            cmap = SRNNModelCellTypes.type_colors(obj.n_types);

            % Panel list: {ylabel, data (n x nt)}; mechanism panels reduced to per-neuron.
            panels = {'stim', pd.u_ext; 'dendrite', pd.x; 'firing rate', pd.r; ...
                      'synaptic output', pd.br};
            if obj.n_a > 0 && isfield(pd, 'a') && ~isempty(pd.a)
                panels(end+1, :) = {'SFA', reshape(sum(pd.a, 2), n, [])};          % sum over timescales
            end
            if obj.n_b > 0
                panels(end+1, :) = {'STD', reshape(mean(pd.b, 2), n, [])};         % mean over post-types
            end
            if obj.n_u > 0
                panels(end+1, :) = {'STF', reshape(mean(pd.p, 2), n, [])};         % mean over post-types
            end
            has_lya = ~strcmpi(obj.lya_method, 'none') && ~isempty(obj.lya_results);
            n_plots = size(panels, 1) + double(has_lya);

            fig = figure('Color', 'w', 'Name', 'SRNNModelCellTypes');
            tl = tiledlayout(fig, n_plots, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
            ax = gobjects(n_plots, 1);

            for k = 1:size(panels, 1)
                ax(k) = nexttile(tl); hold(ax(k), 'on');
                h = SRNNModelCellTypes.plot_type_traces(ax(k), t, panels{k, 2}, ntype, cmap);
                ylabel(ax(k), panels{k, 1}); box(ax(k), 'off');
                if k == 1
                    idx = find(isgraphics(h));
                    lg = legend(h(idx), obj.type_names(idx)); lg.Layout.Tile = 'east';
                end
                if k < n_plots, set(ax(k), 'XTickLabel', []); end
            end

            if has_lya
                ax(end) = nexttile(tl); hold(ax(end), 'on');
                if strcmpi(obj.lya_method, 'benettin')
                    SRNNModelBase.plot_lyapunov(obj.lya_results, obj.lya_method, {'local', 'EOC', 'value'});
                else
                    SRNNModelBase.plot_lyapunov(obj.lya_results, obj.lya_method);
                end
                box(ax(end), 'off');
            end

            linkaxes(ax, 'x');
            Tp = obj.T_plot; if isempty(Tp), Tp = [t(1), t(end)]; end
            xlim(ax(end), Tp);
            xlabel(ax(end), 'time (s)');
        end
    end

    %% Private helpers
    methods (Access = protected)
        function assign_types(obj)
            % ASSIGN_TYPES Contiguous per-type blocks from type_fractions (each >= 1 neuron).
            K = obj.n_types; n = obj.n;
            fr = obj.type_fractions(:)'; fr = fr / sum(fr);
            cnt = max(round(fr * n), 1);
            % fix rounding so sum == n
            while sum(cnt) > n, [~, i] = max(cnt); cnt(i) = cnt(i) - 1; end
            while sum(cnt) < n, [~, i] = max(fr);  cnt(i) = cnt(i) + 1; fr(i) = -inf; end
            t = zeros(n, 1); pos = 0;
            for k = 1:K
                t(pos + (1:cnt(k))) = k; pos = pos + cnt(k);
            end
            obj.type_of = t;
            obj.is_exc = (t == 1);
        end

        function load_parameter_tables(obj)
            % LOAD_PARAMETER_TABLES Fill the K x K tables from Campagnola data or defaults.
            K = obj.n_types;
            haveData = false;
            if obj.use_campagnola_data && exist('load_campagnola_matrices', 'file') == 2
                try
                    C = load_campagnola_matrices();
                    haveData = true;
                catch ME
                    warning('SRNNModelCellTypes:NoData', 'load_campagnola_matrices failed (%s); using defaults.', ME.message);
                end
            end

            if haveData
                cp   = C.conn_prob_adj;
                psp  = C.psp_amplitude;
                dtau = C.ml_depression_tau;
                damt = C.ml_depression_amount;
                relp = C.ml_release_prob;
                ftau = C.ml_facilitation_tau;
                famt = C.ml_facilitation_amount;
                adix = C.sfa_adaptation_index(:);
                % NaN fallbacks (under-sampled elements, e.g. Vip->Pyr)
                cp   = obj.fillnan(cp, 0.05);
                psp  = obj.fill_psp_nan(psp);
                dtau = min(max(obj.fillnan(dtau, 1.0), 0.05), 5);
                damt = min(max(obj.fillnan(damt, 0.3), 0), 0.95);
                relp = min(max(obj.fillnan(relp, 0.2), 0.05), 0.95);
                ftau = min(max(obj.fillnan(ftau, 1.0), 0.05), 5);
                famt = max(obj.fillnan(famt, 0.2), 0);
                adix = obj.fillnan(adix, 0.03);
                if isfield(C, 'sfa_tau'), stau = C.sfa_tau(:); else, stau = nan(K, 1); end
                % Per-type membrane tau (dendritic tau_d). Data is REQUIRED here.
                if ~isfield(C, 'tau') || any(~isfinite(C.tau(:))) || any(C.tau(:) <= 0)
                    error('SRNNModelCellTypes:NoMembraneTau', ...
                        ['Campagnola per-type membrane tau (C.tau) is missing or invalid; ', ...
                         'cannot set per-type tau_d. Set use_campagnola_data=false to run without data.']);
                end
                mtau = C.tau(:);
            else
                % Synthetic defaults so the class runs stand-alone. Pyr (row 1) excitatory (+).
                cp   = 0.1 * ones(K) + diag(-0.05 * ones(K, 1));
                psp  = -3e-4 * ones(K); psp(1, :) = 3e-4;
                dtau = 1.0 * ones(K);
                damt = 0.3 * ones(K);
                relp = 0.2 * ones(K);
                ftau = 1.0 * ones(K);
                famt = 0.2 * ones(K);
                adix = 0.03 * ones(K, 1);
                stau = nan(K, 1);
                mtau = [];   % no data -> tau_d_type stays empty; get_params uses scalar tau_d
            end

            % Per-type SFA tau: fitted value where available, else the scalar tau_a fallback.
            ta_scalar = obj.tau_a; if ~isscalar(ta_scalar), ta_scalar = ta_scalar(1); end
            stau(isnan(stau)) = ta_scalar;

            if isempty(obj.conn_prob),   obj.conn_prob   = cp;   end
            if isempty(obj.psp_amp),     obj.psp_amp     = psp;  end
            if isempty(obj.dep_tau),     obj.dep_tau     = dtau; end
            if isempty(obj.dep_amount),  obj.dep_amount  = damt; end
            if isempty(obj.rel_prob),    obj.rel_prob    = relp; end
            if isempty(obj.fac_tau),     obj.fac_tau     = ftau; end
            if isempty(obj.kappa),       obj.kappa       = famt; end
            if isempty(obj.adapt_index), obj.adapt_index = adix(:); end
            if isempty(obj.tau_a_type),  obj.tau_a_type  = stau(:); end
            if isempty(obj.tau_d_type),  obj.tau_d_type  = mtau(:); end
        end

        function M = fill_psp_nan(~, M)
            % Fill NaN synaptic weights with a signed default (pre-type = row sign).
            mag = median(abs(M(~isnan(M))), 'omitnan'); if isempty(mag) || mag == 0, mag = 3e-4; end
            for r = 1:size(M, 1)
                sgn = 1; if r > 1, sgn = -1; end   % row 1 (pyr) excitatory (+), others (-)
                col = isnan(M(r, :));
                M(r, col) = sgn * mag;
            end
        end
    end

    methods (Static)
        function M = fillnan(M, val)
            M(isnan(M)) = val;
        end

        function c = type_colors(K)
            % TYPE_COLORS Fixed cell-type palette (Pyr/Pvalb/Sst/Vip), matching the
            % project's other figures; falls back to lines() for K > 4.
            base = [0.1216 0.4667 0.7059;   % pyr   #1f77b4
                    0.8392 0.1529 0.1569;   % pvalb #d62728
                    0.1725 0.6275 0.1725;   % sst   #2ca02c
                    0.5804 0.4039 0.7412];  % vip   #9467bd
            if nargin < 1 || isempty(K), K = 4; end
            if K <= 4, c = base(1:K, :); else, c = lines(K); end
        end

        function h = plot_type_traces(ax, t, Y, type_of, cmap)
            % PLOT_TYPE_TRACES Plot every neuron's trace (Y is n x nt) faint in its cell
            % type's color, plus a bold per-type mean. Returns the K mean-line handles
            % (gobjects placeholder for absent types) for building a legend.
            K = size(cmap, 1);
            h = gobjects(K, 1);
            for T = 1:K
                sel = (type_of == T);
                if ~any(sel), continue; end
                col = cmap(T, :);
                plot(ax, t, Y(sel, :)', 'Color', [col 0.25], 'LineWidth', 0.5);   % all traces, faint
                h(T) = plot(ax, t, mean(Y(sel, :), 1), 'Color', col, 'LineWidth', 2);  % per-type mean
            end
        end

        function dS_dt = dynamics_fast_ct(t, S, params)
            % DYNAMICS_FAST_CT Per-neuron RHS with post-type-structured synaptic drive.
            u_ext = params.u_interpolant(t)';           % n x 1
            n = params.n; K = params.K;
            n_a = params.n_a; n_b = params.n_b; n_u = params.n_u;
            W = params.W; tau_d = params.tau_d; phi = params.activation_function;
            type_of = params.type_of;
            n_ad = params.n_ad; ad_idx = params.ad_idx;

            % --- unpack ---  (SFA state is ragged: only adapting neurons carry a)
            idx = 0;
            len_a = n_ad * n_a;
            if len_a > 0
                a_ad = reshape(S(idx + (1:len_a)), n_ad, n_a);
                a = zeros(n, n_a); a(ad_idx, :) = a_ad;    % scatter to full for x_eff
            else
                a_ad = []; a = [];
            end
            idx = idx + len_a;
            len_b = n * K * n_b;
            if len_b > 0, b = reshape(S(idx + (1:len_b)), n, K); else, b = []; end
            idx = idx + len_b;
            len_u = n * K * n_u;
            if len_u > 0, p = reshape(S(idx + (1:len_u)), n, K); else, p = []; end
            idx = idx + len_u;
            x = S(idx + (1:n));

            % --- rates ---
            x_eff = x;
            if len_a > 0, x_eff = x_eff - params.c_vec .* sum(a, 2); end
            r = phi(x_eff);                              % n x 1

            % --- per-edge synaptic gains (n x K) ---
            if len_b > 0, b_mat = b; else, b_mat = ones(n, K); end
            if len_u > 0, p_gain = p ./ params.p0_mat; else, p_gain = ones(n, K); end
            eff = p_gain .* b_mat;                        % n x K
            % effective release prob for STD depletion: state p if STF on, else frozen p0
            if len_u > 0, p_eff = p; else, p_eff = params.p0_mat; end   % n x K

            % --- post-type-structured recurrent drive ---
            drive = zeros(n, 1);
            for T = 1:K
                dr = W * (eff(:, T) .* r);
                sel = (type_of == T);
                drive(sel) = dr(sel);
            end
            dx = (-x + drive + u_ext) ./ tau_d;

            % --- state ODEs ---  (da only for adapting neurons; per-neuron tau)
            if len_a > 0, da = (r(ad_idx) - a_ad) ./ params.tau_a(ad_idx, :); else, da = []; end
            if len_b > 0
                db = (1 - b_mat) ./ params.tau_b_rec_mat - (p_eff .* b_mat .* r) ./ params.tau_b_rel_mat;
            else
                db = [];
            end
            if len_u > 0
                dp = (params.p0_mat - p) ./ params.tau_f_mat + params.kappa_mat .* (1 - p) .* r;
            else
                dp = [];
            end

            dS_dt = [da(:); db(:); dp(:); dx];
        end

        function J = compute_Jacobian_fast_ct(S, params)
            % COMPUTE_JACOBIAN_FAST_CT Analytic Jacobian for dynamics_fast_ct.
            n = params.n; K = params.K;
            n_a = params.n_a; n_b = params.n_b; n_u = params.n_u;
            W = sparse(params.W); tau_d = params.tau_d;
            type_of = params.type_of;
            phi = params.activation_function; phip = params.activation_function_derivative;

            n_ad = params.n_ad; ad_idx = params.ad_idx;
            len_a = n * n_a; len_b = n * K * n_b; len_u = n * K * n_u;
            N = len_a + len_b + len_u + n;            % FULL assembly size (compressed at the end)
            sfa_state_len = n_ad * n_a;               % SFA states actually in S (ragged)

            % --- unpack ---  (SFA state is ragged: scatter adapting-only a into full n x n_a)
            idx = 0;
            if sfa_state_len > 0
                a_ad = reshape(S(idx + (1:sfa_state_len)), n_ad, n_a);
                a = zeros(n, n_a); a(ad_idx, :) = a_ad;
            else
                a = [];
            end
            idx = idx + sfa_state_len;
            if len_b > 0, b_mat = reshape(S(idx + (1:len_b)), n, K); else, b_mat = ones(n, K); end
            idx = idx + len_b;
            if len_u > 0, p = reshape(S(idx + (1:len_u)), n, K); else, p = []; end
            idx = idx + len_u;
            x = S(idx + (1:n));

            x_eff = x;
            if len_a > 0, x_eff = x_eff - params.c_vec .* sum(a, 2); end
            r = phi(x_eff);                       % n x 1
            pr = phip(x_eff);                     % n x 1  (phi')
            c = params.c_vec;                     % n x 1

            if len_u > 0, p_gain = p ./ params.p0_mat; else, p_gain = ones(n, K); end
            eff = p_gain .* b_mat;                % n x K
            if len_u > 0, p_eff = p; else, p_eff = params.p0_mat; end   % STD depletion release prob

            % --- row/col index blocks ---
            row_a = 1:len_a;
            row_b = len_a + (1:len_b);
            row_u = len_a + len_b + (1:len_u);
            row_x = len_a + len_b + len_u + (1:n);

            J = sparse(N, N);

            % ---------- SFA rows (per-neuron tau; full n x n_a, compressed to adapting later) ----------
            if len_a > 0
                TA = params.tau_a;                                 % n x n_a (per-neuron SFA tau)
                gamma = c .* pr;                                   % n x 1
                % da/da: block(k,l) diag = -(1/tau_{:,k}).*c.*pr  (+ -1/tau_{:,k} on k==l)
                rda = zeros(n * n_a * n_a, 1); cda = rda; vda = rda; e = 0;
                for k = 1:n_a
                    invk = 1 ./ TA(:, k);                          % n x 1
                    for l = 1:n_a
                        v = -invk .* gamma;
                        if k == l, v = v - invk; end
                        rng_e = e + (1:n); e = e + n;
                        rda(rng_e) = (k - 1) * n + (1:n)';
                        cda(rng_e) = (l - 1) * n + (1:n)';
                        vda(rng_e) = v;
                    end
                end
                J(row_a, row_a) = sparse(rda, cda, vda, len_a, len_a);
                % da/dx: entry (k,j) -> x_j = (1/tau_{j,k}).*pr_j
                vals_ax = zeros(len_a, 1);
                for k = 1:n_a, vals_ax((k - 1) * n + (1:n)) = pr ./ TA(:, k); end
                J(row_a, row_x) = sparse((1:len_a)', repmat((1:n)', n_a, 1), vals_ax, len_a, n);
            end

            % ---------- STD rows (per (i,T)) ----------
            if len_b > 0
                trec = params.tau_b_rec_mat; trel = params.tau_b_rel_mat;   % n x K
                rrep = repmat(r, K, 1); prrep = repmat(pr, K, 1); crep = repmat(c .* pr, K, 1);
                diag_b = -1 ./ trec(:) - (p_eff(:) .* rrep) ./ trel(:);
                J(row_b, row_b) = spdiags(diag_b, 0, len_b, len_b);         % db/db
                vals_bx = -(p_eff(:) .* b_mat(:) ./ trel(:)) .* prrep;      % db/dx
                J(row_b, row_x) = sparse((1:len_b)', repmat((1:n)', K, 1), vals_bx, len_b, n);
                if len_a > 0                                                % db/da
                    coeff_ba = (p_eff(:) .* b_mat(:) ./ trel(:)) .* crep;
                    stack = sparse((1:len_b)', repmat((1:n)', K, 1), coeff_ba, len_b, n);
                    J(row_b, row_a) = kron(sparse(ones(1, n_a)), stack);
                end
                if len_u > 0                                                % db/dp (depletion coupling)
                    diag_bp = -(b_mat(:) .* rrep) ./ trel(:);
                    J(row_b, row_u) = spdiags(diag_bp, 0, len_b, len_u);
                end
            end

            % ---------- STF rows (per (i,T)) ----------
            if len_u > 0
                tf = params.tau_f_mat; kappa = params.kappa_mat;           % n x K
                rrep = repmat(r, K, 1); prrep = repmat(pr, K, 1); crep = repmat(c .* pr, K, 1);
                diag_u = -1 ./ tf(:) - kappa(:) .* rrep;
                J(row_u, row_u) = spdiags(diag_u, 0, len_u, len_u);        % dp/dp
                vals_ux = (kappa(:) .* (1 - p(:))) .* prrep;              % dp/dx
                J(row_u, row_x) = sparse((1:len_u)', repmat((1:n)', K, 1), vals_ux, len_u, n);
                if len_a > 0                                                % dp/da
                    coeff_ua = -(kappa(:) .* (1 - p(:))) .* crep;
                    stack = sparse((1:len_u)', repmat((1:n)', K, 1), coeff_ua, len_u, n);
                    J(row_u, row_a) = kron(sparse(ones(1, n_a)), stack);
                end
            end

            % ---------- dx rows (post-type partitioned) ----------
            Jxx = -spdiags(1 ./ tau_d(:), 0, n, n);   % tau_d is n x 1 (per-neuron membrane tau)
            Jxa = sparse(n, len_a);
            Jxb = sparse(n, len_b);
            Jxu = sparse(n, len_u);
            for T = 1:K
                rows = find(type_of == T);
                if isempty(rows), continue; end
                gcol = eff(:, T) .* pr;                          % n x 1  (dx/dx column scale)
                Jxx(rows, :) = Jxx(rows, :) + (W(rows, :) .* gcol') ./ tau_d(rows);
                if len_a > 0
                    acol = (-c .* pr .* eff(:, T));              % n x 1
                    blk = (W(rows, :) .* acol') ./ tau_d(rows);  % m x n
                    Jxa(rows, :) = repmat(blk, 1, n_a);
                end
                if len_b > 0
                    bcol = r .* p_gain(:, T);                    % dx/db_{:,T}
                    Jxb(rows, (T - 1) * n + (1:n)) = (W(rows, :) .* bcol') ./ tau_d(rows);
                end
                if len_u > 0
                    pcol = r .* (b_mat(:, T) ./ params.p0_mat(:, T));   % dx/dp_{:,T}
                    Jxu(rows, (T - 1) * n + (1:n)) = (W(rows, :) .* pcol') ./ tau_d(rows);
                end
            end
            J(row_x, row_x) = Jxx;
            if len_a > 0, J(row_x, row_a) = Jxa; end
            if len_b > 0, J(row_x, row_b) = Jxb; end
            if len_u > 0, J(row_x, row_u) = Jxu; end

            % Compress the full n x n_a SFA block to the adapting-neuron subset (ragged state):
            % keep only adapting a rows/cols; deleting the non-adapting a rows+cols yields the
            % exact ragged Jacobian. keep_a ordered by timescale to match a_ad(:).
            if n_a > 0
                keep_a = reshape((0:n_a - 1) * n + ad_idx, [], 1);   % n_ad*n_a
                keep = [keep_a; ((len_a + 1):N)'];
                J = J(keep, keep);
            end
        end

        function st = unpack_states_ct(S_out, params)
            % UNPACK_STATES_CT Reconstruct x,a,b,u,r,br as n x (K x) nt arrays.
            nt = size(S_out, 1);
            n = params.n; K = params.K;
            n_a = params.n_a; n_b = params.n_b; n_u = params.n_u;
            idx = 0;
            n_ad = params.n_ad; ad_idx = params.ad_idx;
            len_a = n_ad * n_a;                        % ragged SFA state length (guard below)
            if len_a > 0
                a_ad = reshape(S_out(:, idx + (1:len_a))', n_ad, n_a, nt);
                a = zeros(n, n_a, nt); a(ad_idx, :, :) = a_ad;   % scatter to full for plotting
            else
                a = [];
            end
            idx = idx + len_a;
            len_b = n * K * n_b;
            if len_b > 0, b = reshape(S_out(:, idx + (1:len_b))', n, K, nt); else, b = ones(n, K, nt); end
            idx = idx + len_b;
            len_u = n * K * n_u;
            if len_u > 0, p = reshape(S_out(:, idx + (1:len_u))', n, K, nt); else, p = ones(n, K, nt); end
            idx = idx + len_u;
            x = S_out(:, idx + (1:n))';                          % n x nt

            x_eff = x;
            if len_a > 0
                sa = reshape(sum(a, 2), n, nt);
                x_eff = x_eff - params.c_vec .* sa;
            end
            r = params.activation_function(x_eff);               % n x nt
            % synaptic output magnitude proxy: efficacy to each post-type times r
            if len_u > 0
                ug = p ./ params.p0_mat;                        % n x K x nt
            else
                ug = ones(n, K, nt);
            end
            eff = ug .* b;                                       % n x K x nt
            br = reshape(mean(eff, 2), n, nt) .* r;              % mean efficacy over post-types * r

            st = struct('x', x, 'a', a, 'b', b, 'p', p, 'r', r, 'br', br);
        end
    end
end
