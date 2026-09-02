classdef RMTBlocks < handle
    % RMTBlocks - Sparse random connectivity with block-structured mean AND variance.
    %
    % Generalizes Harris et al. (2023) from presynaptic-only (column) statistics to
    % full D x D block statistics indexed by BOTH the postsynaptic (row) and the
    % presynaptic (column) cell type:
    %
    %   W = S .* ( A .* (g_sigma * Sigma) + g_mu * M ) + shift * eye(N)
    %
    %   Sigma_ij = sigma_tilde( type(i), type(j) )      N x N, block-constant
    %   M_ij     = mu_tilde(    type(i), type(j) )      N x N, block-constant
    %   A_ij     ~ N(0,1) i.i.d.,   S_ij ~ Bernoulli(alpha)
    %
    % INDEX CONVENTION: (a <- b) = row a POSTsynaptic, column b PREsynaptic, matching
    % A in dx/dt = A*x. So for D = 2 with type 1 = E and type 2 = I:
    %
    %   mu_tilde(1,1) = mu_EE   E receives from E
    %   mu_tilde(1,2) = mu_EI   E receives from I
    %   mu_tilde(2,1) = mu_IE   I receives from E
    %   mu_tilde(2,2) = mu_II   I receives from I
    %
    % Dale's law is a COLUMN constraint: signs must be constant down each column of
    % mu_tilde (all weights from an E neuron are positive, from an I neuron negative).
    % Magnitudes may differ freely between rows.
    %
    % WHY ELEMENTWISE. Harris Eq. (6) writes W = S o (A*D + M) with D diagonal. Right-
    % multiplication by a diagonal is pure COLUMN scaling and is structurally incapable
    % of row dependence, so block sigma requires the elementwise A .* Sigma above. Note
    % also that Ahmadian et al. (2015)'s M + LJR ensemble only covers SEPARABLE variance
    % (sigma_ij = l_i * r_j), which cannot express four independent sigma blocks; the
    % relevant ensemble is Aljadeff, Renfrew & Stern (2015, J. Math. Phys.), which allows
    % an arbitrary D x D block variance matrix.
    %
    % THEORETICAL PREDICTIONS (see lambda_O and R below) compose Harris's sparse moment
    % matching with Aljadeff's block spectral radius. No published result covers sparse
    % + block mean + block variance simultaneously, so these are predictions to be
    % VERIFIED NUMERICALLY, not citable theorems. They reduce exactly to Harris
    % Eqs. (17)-(18) when the blocks are column-uniform, which is the correctness check
    % (see tests/test_RMTBlocks_equivalence.m).
    %
    % Backward compatibility: with D = 2 the scalar aliases mu_tilde_e, mu_tilde_i,
    % sigma_tilde_e, sigma_tilde_i (and E, I, mu_se, mu_si, sigma_se_sq, sigma_si_sq)
    % behave exactly as in the original RMT class, so existing scripts run unchanged.
    %
    % Usage:
    %   rmt = RMTBlocks(300);
    %   rmt.alpha = 1/3;
    %   rmt.f = 0.5;                                  % scalar -> [0.5 0.5]
    %   rmt.mu_tilde    = [ 0.23 -0.31 ;  0.23 -0.31];   % (a <- b)
    %   rmt.sigma_tilde = [ 0.08  0.08 ;  0.08  0.08];
    %   W = rmt.W;
    %   rmt.display_parameters();
    %
    % References:
    %   Harris, Meffin, Burkitt & Peterson (2023) Phys. Rev. Res. 5, 043132.
    %   Aljadeff, Stern & Sharpee (2015) Phys. Rev. Lett. 114, 088101.
    %   Aljadeff, Renfrew & Stern (2015) J. Math. Phys. 56, 103502.
    %   Ahmadian, Fumarola & Miller (2015) Phys. Rev. E 91, 012820.
    %
    % See also: RMT, RMTMatrix

    %% Stored Properties
    properties
        N               % System size (number of neurons)
        alpha           % Sparsity / connection probability (0 < alpha <= 1)

        f               % 1 x D vector of population fractions, sums to 1

        % Block statistics in Harris tilde notation, indexed (postsynaptic, presynaptic)
        mu_tilde        % D x D normalized block means
        sigma_tilde     % D x D normalized block standard deviations (nonnegative)

        % Global scale factors. g_mu = g_sigma reproduces a single "level_of_chaos"
        % multiplier on all of W. Splitting them lets the mean structure and the
        % disorder be swept independently.
        g_mu            % Scale applied to mu_tilde (default 1)
        g_sigma         % Scale applied to sigma_tilde (default 1)

        % Internal matrices
        A               % Base random matrix (Gaussian, mean 0, var 1)
        S               % Sparsity mask (logical)

        % Control flags
        zrs_mode        % 'none', 'ZRS', 'SZRS', 'Partial_SZRS'
        shift           % Scalar diagonal shift added to W (default 0)

        % Visualization
        description
        outlier_threshold  % Multiplier on R separating near/far outliers (default 1.04)
    end

    %% Dependent Properties
    properties (Dependent)
        % Type bookkeeping
        n_types         % D = numel(f)
        n_per_type      % 1 x D neuron counts (may contain zeros)
        type_indices    % 1 x D cell of contiguous neuron index ranges
        type_id         % 1 x N type label per neuron

        % Expanded N x N block-constant matrices
        Sigma           % Sigma_ij = sigma_tilde(type(i), type(j))
        M               % M_ij     = mu_tilde(type(i), type(j))     (Harris Eq. 6)

        % Weight matrix (computed on access)
        W

        % Sparse-effective block statistics (Harris Eqs. 15, 16, applied blockwise,
        % and including g_mu / g_sigma)
        mu_s            % D x D  alpha * g_mu * mu_tilde
        sigma_s_sq      % D x D  alpha(1-alpha)(g_mu*mu_tilde)^2 + alpha(g_sigma*sigma_tilde)^2

        % Theoretical predictions
        lambda_O        % D x 1 outlier eigenvalues, dominant first (generalizes Eq. 17)
        R               % Bulk spectral radius (generalizes Eq. 18)

        % ---- Backward-compatible D = 2 aliases ----
        E               % Logical index of type 1 (excitatory)
        I               % Logical index of type 2 (inhibitory)
        mu_tilde_e      % Column 1 of mu_tilde    (weights FROM E)
        mu_tilde_i      % Column 2 of mu_tilde    (weights FROM I)
        sigma_tilde_e   % Column 1 of sigma_tilde
        sigma_tilde_i   % Column 2 of sigma_tilde
        mu_se           % Column 1 of mu_s        (Harris Eq. 15)
        mu_si           % Column 2 of mu_s
        sigma_se_sq     % Column 1 of sigma_s_sq  (Harris Eq. 16)
        sigma_si_sq     % Column 2 of sigma_s_sq
    end

    %% Private Properties
    properties (Access = private)
        eigenvalues_cache   % Cached eigenvalue computation
        eigenvalues_valid   % Flag indicating if cache is valid
        warned_zrs          % One-shot flag for the ZRS-with-blocks warning
    end

    methods
        function obj = RMTBlocks(N, varargin)
            % RMTBLOCKS Constructor.
            %
            %   rmt = RMTBlocks(N)
            %   rmt = RMTBlocks(N, 'f', [...], 'mu_tilde', [...], 'sigma_tilde', [...])
            %
            % With no name-value pairs the defaults reproduce RMT(N) exactly: D = 2,
            % alpha = 1, f = [0.5 0.5], zero means, sigma_tilde = 1/sqrt(N) everywhere.
            % The RNG is consumed in the same order as RMT(N) (randn(N) then rand(N)),
            % so scripts that rely on a seeded stream reproduce identically.

            obj.N = N;

            % --- Defaults (no RNG consumed) ---
            obj.alpha = 1.0;                        % A is empty here, so no resample
            % mu_tilde before f: set.f derives D from mu_tilde when given a scalar.
            obj.mu_tilde = zeros(2, 2);
            obj.sigma_tilde = ones(2, 2) / sqrt(N);
            obj.f = [0.5, 0.5];
            obj.g_mu = 1;
            obj.g_sigma = 1;
            obj.zrs_mode = 'none';
            obj.shift = 0;
            obj.description = '';
            obj.outlier_threshold = 1.04;

            obj.eigenvalues_cache = [];
            obj.eigenvalues_valid = false;
            obj.warned_zrs = false;

            % --- Random matrices (same order as RMT.m) ---
            obj.A = randn(N, N);        % Mean 0, Var 1
            obj.update_sparsity();

            % --- Optional name-value overrides ---
            if ~isempty(varargin)
                if mod(numel(varargin), 2) ~= 0
                    error('RMTBlocks:InvalidInput', ...
                        'Optional arguments must be name-value pairs.');
                end
                for k = 1:2:numel(varargin)
                    name = varargin{k};
                    if ~ischar(name) && ~isstring(name)
                        error('RMTBlocks:InvalidInput', ...
                            'Property names must be character vectors or strings.');
                    end
                    name = char(name);
                    if ~isprop(obj, name)
                        error('RMTBlocks:UnknownProperty', ...
                            'Unknown property "%s".', name);
                    end
                    obj.(name) = varargin{k+1};
                end
            end
        end

        %% ---------------- Type bookkeeping ----------------
        function val = get.n_types(obj)
            val = numel(obj.f);
        end

        function val = get.n_per_type(obj)
            % Cumulative rounding. At D = 2 this yields [round(f1*N), N-round(f1*N)],
            % identical to RMT's round(f*N) convention, which is what makes bit-for-bit
            % equivalence with the original class possible. Always sums to N exactly.
            % Zero counts are permitted (an empty population is legal).
            boundaries = round(cumsum(obj.f) * obj.N);
            boundaries(end) = obj.N;                % guard against rounding drift
            val = diff([0, boundaries]);
        end

        function val = get.type_indices(obj)
            counts = obj.n_per_type;
            val = cell(1, numel(counts));
            first = 1;
            for q = 1:numel(counts)
                last = first + counts(q) - 1;
                val{q} = first:last;                % empty range when counts(q) == 0
                first = last + 1;
            end
        end

        function val = get.type_id(obj)
            val = repelem(1:obj.n_types, obj.n_per_type);
        end

        %% ---------------- Expanded matrices ----------------
        function val = get.Sigma(obj)
            obj.validate_consistency();
            tid = obj.type_id;
            val = obj.sigma_tilde(tid, tid);
        end

        function val = get.M(obj)
            obj.validate_consistency();
            tid = obj.type_id;
            val = obj.mu_tilde(tid, tid);
        end

        function val = get.W(obj)
            obj.validate_consistency();
            obj.warn_if_zrs_with_blocks();

            % Random component (pre-mask) and mean component, both scaled.
            % A .* Sigma is the elementwise generalization of Harris's A*D.
            AS = obj.A .* (obj.g_sigma * obj.Sigma);
            Mg = obj.g_mu * obj.M;

            switch obj.zrs_mode
                case 'none'
                    % Harris Eq. (6): W = S o (A*D + M)
                    val = obj.S .* (AS + Mg);

                case 'ZRS'
                    % Dense ZRS via projection operator (Harris Eqs. 24, 25).
                    % P acts on the random component only, so the mean structure -
                    % and hence the outliers - survive.
                    if obj.alpha < 1
                        warning('RMTBlocks:SparsityWarning', ...
                            ['Using ''ZRS'' (projection) with a sparse matrix. ', ...
                             'This destroys sparsity. Consider ''SZRS''.']);
                    end

                    u = ones(obj.N, 1);
                    P = eye(obj.N) - (u * u') / obj.N;       % Eq. (24)

                    val = (AS * P) + Mg;                     % Eq. (25)

                    if obj.alpha < 1
                        val = obj.S .* val;
                    end

                case 'SZRS'
                    % Sparse zero row sum (Harris Eqs. 30, 31): subtract each row's
                    % mean over its nonzero entries. Applied to BOTH components, so it
                    % also removes the imposed imbalance (Harris explicitly notes this
                    % at Eq. 32's motivation) - use 'Partial_SZRS' to keep the means.
                    W_base = obj.S .* (AS + Mg);

                    row_sums = sum(W_base, 2);
                    row_counts = sum(obj.S, 2);
                    row_counts(row_counts == 0) = 1;         % avoid divide-by-zero
                    W_bar_i = row_sums ./ row_counts;        % Eq. (31)

                    B = obj.S .* W_bar_i;
                    val = W_base - B;                        % Eq. (30)

                case 'Partial_SZRS'
                    % Harris Eq. (32): correct the random component only, preserving
                    % the imposed imbalance. This is the mode that composes correctly
                    % with block-structured means.
                    J_base = obj.S .* AS;
                    M_base = obj.S .* Mg;

                    J_row_sums = sum(J_base, 2);
                    row_counts = sum(obj.S, 2);
                    row_counts(row_counts == 0) = 1;
                    J_bar_i = J_row_sums ./ row_counts;

                    B_partial = obj.S .* J_bar_i;
                    val = (J_base - B_partial) + M_base;

                otherwise
                    error('RMTBlocks:InvalidZRSMode', ...
                        'Unknown zrs_mode ''%s''.', obj.zrs_mode);
            end

            % Apply diagonal shift
            val = val + obj.shift * eye(obj.N);
        end

        %% ---------------- Sparse-effective block statistics ----------------
        function val = get.mu_s(obj)
            % Harris Eq. (15) applied blockwise: mu_s = alpha * mu_tilde
            val = obj.alpha * (obj.g_mu * obj.mu_tilde);
        end

        function val = get.sigma_s_sq(obj)
            % Harris Eq. (16) applied blockwise:
            %   sigma_s^2 = alpha(1-alpha) mu_tilde^2 + alpha sigma_tilde^2
            % The first term is why sparsity couples the MEAN into the bulk radius.
            mg = obj.g_mu * obj.mu_tilde;
            sg = obj.g_sigma * obj.sigma_tilde;
            val = obj.alpha * (1 - obj.alpha) * mg.^2 + obj.alpha * sg.^2;
        end

        %% ---------------- Theoretical predictions ----------------
        function val = get.lambda_O(obj)
            % Outlier eigenvalues. E[W] is block-constant with entries mu_s(a,b), so
            % its action on type-constant vectors is the D x D matrix
            %
            %   K(a,b) = n_b * mu_s(a,b) = N * f_b * mu_s(a,b)
            %
            % and the nonzero eigenvalues of E[W] are the eigenvalues of K.
            %
            % Returns a D x 1 vector sorted by descending MAGNITUDE. With D cell types
            % there are up to D outliers, so there is no single "the" outlier. Sorting
            % by |lambda| rather than by real part is deliberate: an outlier is defined
            % by lying outside the bulk disk of radius R, and in an inhibition-dominated
            % network the dominant outlier is large and NEGATIVE, so a real-part sort
            % would rank a meaningless zero eigenvalue above it.
            %
            % Column-uniform reduction: mu_s(a,b) = mu_s(b) makes every row of K
            % identical, hence rank 1, so the eigenvalues are
            %   { N[f*mu_se + (1-f)*mu_si], 0 }  =  { Harris Eq. (17), 0 }.
            %
            % Empty populations: f_b = 0 zeroes column b of K, so the surviving
            % eigenvalues are exactly the lower-D answer.
            %
            % NOTE: this prediction assumes zrs_mode = 'none'. Dense 'ZRS' and
            % 'Partial_SZRS' preserve the mean structure so it still applies, but
            % 'SZRS' zeroes all row sums and invalidates it (see warn_if_zrs_with_blocks).
            obj.validate_consistency();
            K = obj.N * (obj.f .* obj.mu_s);        % f is 1 x D -> scales columns
            lam = eig(K);
            [~, order] = sort(abs(lam), 'descend');
            val = lam(order);
        end

        function val = get.R(obj)
            % Bulk spectral radius. Aljadeff et al. (2015) give the transition at
            % Lambda_1 = 1 where Lambda_1 is the largest eigenvalue (Perron root) of
            % the block variance matrix; in Harris's normalization that matrix is
            %
            %   V(a,b) = N * f_b * sigma_s_sq(a,b)
            %
            % and the radius is R = sqrt(Lambda_1). V is entrywise nonnegative, so
            % Perron-Frobenius guarantees a real nonnegative dominant eigenvalue.
            %
            % Column-uniform reduction: every row of V is identical, hence rank 1, so
            % Lambda_1 = trace = N[f*sigma_se^2 + (1-f)*sigma_si^2] and
            % R = sqrt(that) = Harris Eq. (18).
            %
            % R remains the meaningful radius under every zrs_mode: the ZRS conditions
            % exist precisely to pull stray eigenvalues back INSIDE this radius.
            obj.validate_consistency();
            V = obj.N * (obj.f .* obj.sigma_s_sq);
            lam = eig(V);
            val = sqrt(max(0, max(real(lam))));
        end

        %% ---------------- Backward-compatible D = 2 aliases ----------------
        function val = get.E(obj)
            obj.assert_two_types('E');
            val = false(obj.N, 1);
            counts = obj.n_per_type;
            val(1:counts(1)) = true;
        end

        function val = get.I(obj)
            obj.assert_two_types('I');
            val = ~obj.E;
        end

        function val = get.mu_tilde_e(obj)
            val = obj.uniform_column(obj.mu_tilde, 1, 'mu_tilde_e');
        end

        function set.mu_tilde_e(obj, val)
            obj.assert_two_types('mu_tilde_e');
            obj.mu_tilde(:, 1) = val;
        end

        function val = get.mu_tilde_i(obj)
            val = obj.uniform_column(obj.mu_tilde, 2, 'mu_tilde_i');
        end

        function set.mu_tilde_i(obj, val)
            obj.assert_two_types('mu_tilde_i');
            obj.mu_tilde(:, 2) = val;
        end

        function val = get.sigma_tilde_e(obj)
            val = obj.uniform_column(obj.sigma_tilde, 1, 'sigma_tilde_e');
        end

        function set.sigma_tilde_e(obj, val)
            obj.assert_two_types('sigma_tilde_e');
            obj.sigma_tilde(:, 1) = val;
        end

        function val = get.sigma_tilde_i(obj)
            val = obj.uniform_column(obj.sigma_tilde, 2, 'sigma_tilde_i');
        end

        function set.sigma_tilde_i(obj, val)
            obj.assert_two_types('sigma_tilde_i');
            obj.sigma_tilde(:, 2) = val;
        end

        function val = get.mu_se(obj)
            val = obj.uniform_column(obj.mu_s, 1, 'mu_se');
        end

        function val = get.mu_si(obj)
            val = obj.uniform_column(obj.mu_s, 2, 'mu_si');
        end

        function val = get.sigma_se_sq(obj)
            val = obj.uniform_column(obj.sigma_s_sq, 1, 'sigma_se_sq');
        end

        function val = get.sigma_si_sq(obj)
            val = obj.uniform_column(obj.sigma_s_sq, 2, 'sigma_si_sq');
        end

        %% ---------------- Property setters (with cache invalidation) ----------------
        function set.alpha(obj, val)
            if ~isscalar(val) || ~isnumeric(val) || ~isfinite(val) || val <= 0 || val > 1
                error('RMTBlocks:InvalidParams', 'alpha must satisfy 0 < alpha <= 1.');
            end
            obj.alpha = val;
            if ~isempty(obj.A)
                obj.update_sparsity();
            end
            obj.invalidate_eigenvalues();
        end

        function set.f(obj, val)
            val = reshape(val, 1, []);

            if isscalar(val)
                % Scalar convenience: only meaningful for D <= 2. D is derived from
                % mu_tilde rather than from the current f, so that set_types can assign
                % the block matrices first and then pass a scalar f without the old D
                % leaking in (set_types(1, 1x1, 1x1) must give f = 1, not [1 0]).
                if ~isempty(obj.mu_tilde)
                    D = size(obj.mu_tilde, 1);
                elseif ~isempty(obj.f)
                    D = numel(obj.f);
                else
                    D = 2;                          % construction-time default
                end
                if D == 2
                    val = [val, 1 - val];
                elseif D == 1
                    if abs(val - 1) > 1e-12
                        error('RMTBlocks:InvalidParams', ...
                            'With a single cell type f must equal 1 (got %g).', val);
                    end
                else
                    error('RMTBlocks:InvalidParams', ...
                        ['Scalar f is only allowed when D <= 2 (current D = %d). ', ...
                         'Use set_types(f, mu_tilde, sigma_tilde) to change D.'], D);
                end
            end

            if any(~isfinite(val)) || any(val < 0)
                error('RMTBlocks:InvalidParams', ...
                    'f must contain finite nonnegative fractions.');
            end
            if abs(sum(val) - 1) > 1e-12
                error('RMTBlocks:InvalidParams', ...
                    'f must sum to 1 (got %.12g).', sum(val));
            end

            obj.f = val;
            obj.invalidate_eigenvalues();
        end

        function set.mu_tilde(obj, val)
            if ~isnumeric(val) || ~ismatrix(val) || size(val, 1) ~= size(val, 2)
                error('RMTBlocks:InvalidParams', 'mu_tilde must be a square matrix.');
            end
            if any(~isfinite(val(:)))
                error('RMTBlocks:InvalidParams', 'mu_tilde must be finite.');
            end
            obj.mu_tilde = val;
            obj.invalidate_eigenvalues();
        end

        function set.sigma_tilde(obj, val)
            if ~isnumeric(val) || ~ismatrix(val) || size(val, 1) ~= size(val, 2)
                error('RMTBlocks:InvalidParams', 'sigma_tilde must be a square matrix.');
            end
            if any(~isfinite(val(:))) || any(val(:) < 0)
                error('RMTBlocks:InvalidParams', ...
                    'sigma_tilde must be finite and nonnegative.');
            end
            obj.sigma_tilde = val;
            obj.invalidate_eigenvalues();
        end

        function set.g_mu(obj, val)
            if ~isscalar(val) || ~isnumeric(val) || ~isfinite(val)
                error('RMTBlocks:InvalidParams', 'g_mu must be a finite scalar.');
            end
            obj.g_mu = val;
            obj.invalidate_eigenvalues();
        end

        function set.g_sigma(obj, val)
            if ~isscalar(val) || ~isnumeric(val) || ~isfinite(val)
                error('RMTBlocks:InvalidParams', 'g_sigma must be a finite scalar.');
            end
            obj.g_sigma = val;
            obj.invalidate_eigenvalues();
        end

        function set.zrs_mode(obj, val)
            valid_modes = {'none', 'ZRS', 'SZRS', 'Partial_SZRS'};
            if ~ismember(val, valid_modes)
                error('RMTBlocks:InvalidZRSMode', ...
                    'Invalid ZRS mode. Valid choices: %s', strjoin(valid_modes, ', '));
            end
            obj.zrs_mode = val;
            obj.warned_zrs = false;
            obj.invalidate_eigenvalues();
        end

        function set.shift(obj, val)
            obj.shift = val;
            obj.invalidate_eigenvalues();
        end

        function set.A(obj, val)
            obj.A = val;
            obj.invalidate_eigenvalues();
        end

        function set.S(obj, val)
            obj.S = val;
            obj.invalidate_eigenvalues();
        end

        %% ---------------- Parameter setters ----------------
        function set_types(obj, f, mu_tilde, sigma_tilde)
            % SET_TYPES Set the number of cell types and all block statistics atomically.
            %
            %   rmt.set_types(f, mu_tilde, sigma_tilde)
            %
            % This is the only way to CHANGE D. Setting f, mu_tilde and sigma_tilde one
            % at a time works fine when D is unchanged, but changing D piecemeal would
            % leave the object transiently inconsistent, so use this instead.
            D = numel(f);
            if ~isequal(size(mu_tilde), [D, D]) || ~isequal(size(sigma_tilde), [D, D])
                error('RMTBlocks:InconsistentTypes', ...
                    ['numel(f) = %d requires mu_tilde and sigma_tilde to be %dx%d ', ...
                     '(got %s and %s).'], D, D, D, ...
                    mat2str(size(mu_tilde)), mat2str(size(sigma_tilde)));
            end

            % Assign the matrices first so a scalar f never sees a stale D.
            obj.mu_tilde = mu_tilde;
            obj.sigma_tilde = sigma_tilde;
            obj.f = reshape(f, 1, []);
        end

        function set_params(obj, mu_tilde_e, mu_tilde_i, sigma_tilde_e, sigma_tilde_i, f, alpha)
            % SET_PARAMS Backward-compatible D = 2 setter in Harris notation.
            %
            % Signature preserved from RMT/RMTMatrix. Each scalar sets a whole COLUMN
            % of the corresponding block matrix.
            if nargin > 1, obj.mu_tilde_e = mu_tilde_e; end
            if nargin > 2, obj.mu_tilde_i = mu_tilde_i; end
            if nargin > 3, obj.sigma_tilde_e = sigma_tilde_e; end
            if nargin > 4, obj.sigma_tilde_i = sigma_tilde_i; end
            if nargin > 5, obj.f = f; end
            if nargin > 6, obj.alpha = alpha; end
        end

        function set_alpha(obj, alpha)
            % SET_ALPHA Set alpha (sparsity) independently.
            obj.alpha = alpha;
        end

        function set_zrs_mode(obj, mode)
            % SET_ZRS_MODE Set the zero row-sum mode.
            obj.zrs_mode = mode;
        end

        function sigma_tilde_i = compute_sigma_tilde_i_for_target_variance(obj, target_variance)
            % COMPUTE_SIGMA_TILDE_I_FOR_TARGET_VARIANCE Solve for sigma_tilde_i so that
            % Var(W) hits a target (Harris Eq. 14 solved for sigma_si^2).
            %
            % Dense (alpha = 1), D = 2, column-uniform only. A block-general version
            % (holding the Perron root of V fixed) is future work.
            obj.assert_two_types('compute_sigma_tilde_i_for_target_variance');

            if obj.alpha < 1
                error('RMTBlocks:SparseCaseNotSupported', ...
                    ['compute_sigma_tilde_i_for_target_variance only supports dense ', ...
                     'matrices (alpha=1). For alpha<1, a more complex solver is needed.']);
            end

            % For the dense case, sigma_se^2 = sigma_tilde_e^2
            sigma_se_sq_local = obj.sigma_tilde_e^2;
            f_E = obj.f(1);

            sigma_si_sq_local = (target_variance - f_E * sigma_se_sq_local) / (1 - f_E);

            if sigma_si_sq_local < 0
                error('RMTBlocks:InvalidTargetVariance', ...
                    ['Target variance %.4f is too small. Minimum achievable is %.4f ', ...
                     '(when sigma_tilde_i=0).'], target_variance, f_E * sigma_se_sq_local);
            end

            sigma_tilde_i = sqrt(sigma_si_sq_local);
        end

        %% ---------------- Internal updates ----------------
        function update_sparsity(obj)
            % UPDATE_SPARSITY Redraw the Bernoulli connection mask.
            obj.S = rand(obj.N, obj.N) < obj.alpha;
            obj.invalidate_eigenvalues();
        end

        %% ---------------- Display and diagnostics ----------------
        function display_parameters(obj)
            % DISPLAY_PARAMETERS Print set / sparse-effective / measured statistics
            % blockwise, alongside the theoretical predictions.
            obj.validate_consistency();

            W_mat = full(obj.W);

            % Remove the diagonal, which may carry the shift
            W_no_diag = W_mat;
            W_no_diag(1:obj.N+1:end) = NaN;

            D = obj.n_types;
            counts = obj.n_per_type;
            idx = obj.type_indices;
            mu_s_mat = obj.mu_s;
            sig_s_mat = obj.sigma_s_sq;

            fprintf('\n========== RMTBlocks Parameter Summary ==========\n');
            if ~isempty(obj.description)
                fprintf('Description:          %s\n', obj.description);
            end
            fprintf('Mode:                 %s\n', obj.zrs_mode);
            fprintf('N:                    %d\n', obj.N);
            fprintf('alpha:                %.4f\n', obj.alpha);
            fprintf('D (cell types):       %d\n', D);
            fprintf('f:                    %s\n', mat2str(obj.f, 4));
            fprintf('n_per_type:           %s\n', mat2str(counts));
            if obj.g_mu ~= 1 || obj.g_sigma ~= 1
                fprintf('g_mu / g_sigma:       %.4f / %.4f\n', obj.g_mu, obj.g_sigma);
            end
            if obj.shift ~= 0
                fprintf('shift:                %.4f\n', obj.shift);
            end

            fprintf('\n--- Blocks (row = postsynaptic, col = presynaptic) ---\n');
            % The measured columns are computed over NONZERO entries only, so at
            % alpha < 1 they track mu_tilde rather than mu_s = alpha*mu_tilde. The
            % (NZ) tag matters: without it the two columns look contradictory.
            fprintf('%-14s %10s %10s %12s %12s %12s %12s\n', ...
                'block(a<-b)', 'mu_tilde', 'sig_tilde', 'mu_s', 'sqrt(sig_s2)', ...
                'meas.mu(NZ)', 'meas.sd(NZ)');
            fprintf('%s\n', repmat('-', 1, 88));
            for a = 1:D
                for b = 1:D
                    if counts(a) == 0 || counts(b) == 0
                        meas_mu = NaN;
                        meas_sd = NaN;
                    else
                        blk = W_no_diag(idx{a}, idx{b});
                        vals = blk(~isnan(blk) & (blk ~= 0));
                        if isempty(vals)
                            meas_mu = NaN;
                            meas_sd = NaN;
                        else
                            meas_mu = mean(vals);
                            meas_sd = std(vals, 1);
                        end
                    end
                    fprintf('%-14s %10.4f %10.4f %12.4f %12.4f %12.4f %12.4f\n', ...
                        sprintf('(%d<-%d)', a, b), ...
                        obj.mu_tilde(a, b), obj.sigma_tilde(a, b), ...
                        mu_s_mat(a, b), sqrt(sig_s_mat(a, b)), ...
                        meas_mu, meas_sd);
                end
            end

            fprintf('\n--- Theoretical predictions ---\n');
            lam = obj.lambda_O;
            for k = 1:numel(lam)
                if abs(imag(lam(k))) > 1e-12
                    fprintf('lambda_O(%d):          %.4f %+.4fi\n', k, real(lam(k)), imag(lam(k)));
                else
                    fprintf('lambda_O(%d):          %.4f\n', k, real(lam(k)));
                end
            end
            fprintf('R (radius):           %.4f\n', obj.R);
            if ~strcmp(obj.zrs_mode, 'none')
                fprintf('(predictions assume zrs_mode = ''none''; see warning notes)\n');
            end
            fprintf('==================================================\n\n');
        end

        %% ---------------- Eigenvalue computation (with caching) ----------------
        function eigs_out = get_eigenvalues(obj)
            % GET_EIGENVALUES Return eigenvalues of W, computing only if the cache is stale.
            if ~obj.eigenvalues_valid || isempty(obj.eigenvalues_cache)
                obj.eigenvalues_cache = eig(obj.W);
                obj.eigenvalues_valid = true;
            end
            eigs_out = obj.eigenvalues_cache;
        end

        function compute_eigenvalues(obj)
            % COMPUTE_EIGENVALUES Force recomputation and cache update.
            obj.eigenvalues_cache = eig(obj.W);
            obj.eigenvalues_valid = true;
        end

        function eigs_out = eigenvalues(obj)
            % EIGENVALUES Property-like access to the cached eigenvalues.
            eigs_out = obj.get_eigenvalues();
        end

        function val = bulk_edge(obj, quantile_level)
            % BULK_EDGE Robust empirical estimate of the bulk radius.
            %
            %   val = rmt.bulk_edge()            % 99th percentile of |lambda|
            %   val = rmt.bulk_edge(q)           % q-th percentile, q in [0, 100]
            %
            % Uses a quantile rather than max(|lambda|) deliberately. Harris documents
            % a small number of "local eigenvalue outliers" escaping the circular
            % support (his Sec. III D 2), which is exactly why his Eq. (18) radius is
            % introduced as approximate. Comparing a prediction against max(|lambda|)
            % would reject a formula that is in fact correct.
            if nargin < 2 || isempty(quantile_level)
                quantile_level = 99;
            end
            lam = obj.get_eigenvalues();
            val = RMTBlocks.percentile(abs(lam), quantile_level);
        end

        %% ---------------- Plotting ----------------
        function plot_spectrum(obj, ax)
            % PLOT_SPECTRUM Eigenvalues of W with the predicted radius circle overlaid.
            % Interior eigenvalues are black circles, near-outliers (within
            % outlier_threshold*R) black crosses, and far outliers green discs.
            if nargin < 2 || isempty(ax)
                figure; ax = gca;
            end

            eigs_val = obj.get_eigenvalues();

            R_val = obj.R;
            xc = obj.shift;
            yc = 0;

            distances = abs(eigs_val - xc - 1i*yc);

            mSize = 4;
            interior_mask = distances <= R_val;
            interior_eigs = eigs_val(interior_mask);
            plot(ax, real(interior_eigs), imag(interior_eigs), 'ko', ...
                'MarkerSize', mSize, 'MarkerFaceColor', 'none', 'LineWidth', 0.5);
            hold(ax, 'on');

            theta = linspace(0, 2*pi, 100);
            plot(ax, xc + R_val*cos(theta), yc + R_val*sin(theta), 'k-', 'LineWidth', 2);

            near_outlier_mask = (distances > R_val) & (distances <= obj.outlier_threshold * R_val);
            near_outlier_eigs = eigs_val(near_outlier_mask);
            if ~isempty(near_outlier_eigs)
                plot(ax, real(near_outlier_eigs), imag(near_outlier_eigs), 'kx', ...
                    'MarkerSize', mSize, 'LineWidth', 0.5);
            end

            far_outlier_mask = distances > obj.outlier_threshold * R_val;
            far_outlier_eigs = eigs_val(far_outlier_mask);
            if ~isempty(far_outlier_eigs)
                plot(ax, real(far_outlier_eigs), imag(far_outlier_eigs), 'o', ...
                    'MarkerSize', mSize, 'MarkerFaceColor', [0 .7 0], ...
                    'MarkerEdgeColor', [0 .7 0]);
            end

            xlabel(ax, 'Re(\lambda)');
            ylabel(ax, 'Im(\lambda)');
            grid(ax, 'on');
            axis(ax, 'equal');
            hold(ax, 'off');
        end

        function plot_W(obj, ax, clim_val)
            % PLOT_W Display W with a diverging red/white/blue colormap and a
            % symmetric, optionally shared color range. The diagonal shift (if any) is
            % rendered as-is and saturates at the clim extremes.
            %
            %   plot_W()             - new figure, auto clim from off-diagonal max
            %   plot_W(ax)           - draw into ax, auto clim
            %   plot_W(ax, clim_val) - draw into ax with shared clim
            if nargin < 2 || isempty(ax)
                figure; ax = gca;
            end

            W_plot = full(obj.W);

            if nargin < 3 || isempty(clim_val)
                W_no_diag = W_plot;
                W_no_diag(1:obj.N+1:end) = 0;
                clim_val = ceil(max(abs(W_no_diag(:))) * 10) / 10;
                if clim_val == 0
                    clim_val = 0.1;
                end
            end

            imagesc(ax, W_plot);
            colormap(ax, redwhiteblue_colormap(256));
            clim(ax, [-clim_val, clim_val]);
            axis(ax, 'square');
            set(ax, 'XTick', [], 'YTick', []);
            box(ax, 'off');
            set(ax, 'Color', 'none', ...
                    'XColor', 'white', 'YColor', 'white', 'Layer', 'bottom');
        end

        function plot_weight_histogram(obj, ax, bin_edges, show_legend, show_yaxis)
            % PLOT_WEIGHT_HISTOGRAM Histogram of the nonzero entries of W, split by
            % PRESYNAPTIC type (i.e. by column block). With a single non-empty
            % population a plain gray histogram is drawn; with two, the classic E/I
            % red/blue overlay with tilde-mu markers; with more, one color per type.
            % Values outside the bin range are clipped to the first/last bin so that
            % diagonal-shift entries surface as edge spikes rather than vanishing.
            %
            %   plot_weight_histogram()
            %   plot_weight_histogram(ax)
            %   plot_weight_histogram(ax, bin_edges)
            %   plot_weight_histogram(ax, bin_edges, show_legend)
            %   plot_weight_histogram(ax, bin_edges, show_legend, show_yaxis)
            if nargin < 2 || isempty(ax)
                figure; ax = gca;
            end

            W_full = full(obj.W);

            if nargin < 3 || isempty(bin_edges)
                r = max(abs(W_full(:)));
                if r == 0, r = 1; end
                bin_edges = linspace(-r, r, 51);
            end
            if nargin < 4 || isempty(show_legend), show_legend = true; end
            if nargin < 5 || isempty(show_yaxis),  show_yaxis  = true; end

            edge_lo = bin_edges(1);
            edge_hi = bin_edges(end);

            counts = obj.n_per_type;
            idx = obj.type_indices;
            present = find(counts > 0);

            if numel(present) <= 1
                % Single population: one gray histogram (matches the old f == 1 path).
                w = W_full(:);
                w = w(w ~= 0);
                w = min(max(w, edge_lo), edge_hi);
                histogram(ax, w, bin_edges, ...
                    'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
            else
                if numel(present) == 2
                    colors = [0.8 0.2 0.2; 0.2 0.4 0.8];    % E red, I blue
                    labels = {'E', 'I'};
                    mu_labels = {'$\tilde{\mu}_E$', '$\tilde{\mu}_I$'};
                else
                    colors = lines(numel(present));
                    labels = arrayfun(@(b) sprintf('type %d', b), present, ...
                        'UniformOutput', false);
                    mu_labels = arrayfun(@(b) sprintf('$\\tilde{\\mu}_{%d}$', b), ...
                        present, 'UniformOutput', false);
                end

                hold(ax, 'on');
                for k = 1:numel(present)
                    b = present(k);
                    w_b = W_full(:, idx{b});
                    w_b = w_b(w_b ~= 0);
                    w_b = min(max(w_b, edge_lo), edge_hi);
                    histogram(ax, w_b, bin_edges, ...
                        'FaceColor', colors(k, :), 'EdgeColor', 'none', 'FaceAlpha', 0.6);
                end
                hold(ax, 'off');

                if show_legend
                    legend(ax, labels, 'Location', 'northeast');
                    legend(ax, 'boxoff');
                end

                % Markers point at the off-diagonal bulk, which is unaffected by
                % obj.shift (shift only displaces diagonal entries, and those are
                % absorbed into the first/last bin by the clipping above). A marker is
                % only meaningful when the column is uniform, i.e. when all
                % postsynaptic types see the same mean from that presynaptic type.
                y_bottom = ax.YLim(1);
                for k = 1:numel(present)
                    b = present(k);
                    col = obj.mu_tilde(:, b);
                    if all(col == col(1))
                        text(ax, obj.g_mu * col(1), y_bottom, mu_labels{k}, ...
                            'Interpreter', 'latex', ...
                            'HorizontalAlignment', 'center', ...
                            'VerticalAlignment', 'top', 'FontSize', 10, ...
                            'Clipping', 'on');
                    end
                end
            end

            xlabel(ax, 'Weight');
            ylabel(ax, 'Count');
            box(ax, 'off');
            set(ax, 'Color', 'none');

            if ~show_yaxis
                ax.YAxis.Visible = 'off';
            end
        end

        function plot_row_sums(obj, ax, xlim_val)
            % PLOT_ROW_SUMS Plot row sums of W as a vertical black trace.
            %   x = row sum (includes diagonal shift), y = row index (1 at top)
            if nargin < 2 || isempty(ax)
                figure; ax = gca;
            end

            row_sums = full(sum(obj.W, 2));
            row_idx = (1:obj.N)';

            if nargin < 3 || isempty(xlim_val)
                xlim_val = max(abs(row_sums)) * 1.1;
                if xlim_val == 0
                    xlim_val = 1;
                end
            end

            plot(ax, row_sums, row_idx, 'k-', 'LineWidth', 1);
            set(ax, 'YDir', 'reverse');
            xlim(ax, [-xlim_val, xlim_val]);
            ylim(ax, [1, obj.N]);
            axis(ax, 'off');
            set(ax, 'Color', 'none');
        end

        %% ---------------- Deep copy ----------------
        function new_obj = copy(obj)
            % COPY Deep copy. The property assignment order mirrors RMT.copy so that
            % the RNG stream advances identically (the constructor draws randn(N) and
            % rand(N), and set.alpha redraws the mask), keeping seeded scripts that
            % chain copy() calls reproducible against the original class.
            new_obj = RMTBlocks(obj.N);

            new_obj.alpha = obj.alpha;
            new_obj.f = obj.f;
            new_obj.mu_tilde = obj.mu_tilde;
            new_obj.sigma_tilde = obj.sigma_tilde;
            new_obj.g_mu = obj.g_mu;
            new_obj.g_sigma = obj.g_sigma;
            new_obj.A = obj.A;
            new_obj.S = obj.S;
            new_obj.zrs_mode = obj.zrs_mode;
            new_obj.shift = obj.shift;
            new_obj.description = obj.description;
            new_obj.outlier_threshold = obj.outlier_threshold;
            new_obj.eigenvalues_cache = obj.eigenvalues_cache;
            new_obj.eigenvalues_valid = obj.eigenvalues_valid;
        end

        %% ---------------- Introspection ----------------
        function tf = is_column_uniform(obj)
            % IS_COLUMN_UNIFORM True when every column of mu_tilde and sigma_tilde is
            % constant down its rows, i.e. when the statistics are presynaptic-only and
            % the object is equivalent to the original Harris/RMT parameterization.
            tf = all(obj.mu_tilde == obj.mu_tilde(1, :), 'all') && ...
                 all(obj.sigma_tilde == obj.sigma_tilde(1, :), 'all');
        end
    end

    methods (Static)
        function val = percentile(x, p)
            % PERCENTILE Linear-interpolation percentile, implemented locally so the
            % class carries no Statistics Toolbox dependency.
            %
            %   val = RMTBlocks.percentile(x, p)   p in [0, 100]
            x = sort(x(:));
            n = numel(x);
            if n == 0
                val = NaN;
                return;
            end
            if n == 1
                val = x;
                return;
            end
            % Midpoint convention: the k-th of n sorted values sits at 100*(k-0.5)/n.
            pos = p / 100 * n + 0.5;
            pos = min(max(pos, 1), n);
            lo = floor(pos);
            hi = ceil(pos);
            if lo == hi
                val = x(lo);
            else
                val = x(lo) + (pos - lo) * (x(hi) - x(lo));
            end
        end
    end

    methods (Access = private)
        function invalidate_eigenvalues(obj)
            obj.eigenvalues_valid = false;
        end

        function validate_consistency(obj)
            % Cross-property validation is done lazily, at use time, so that f and the
            % block matrices can be assigned in any order (Fig_1_RMT_examples.m sets
            % mu_tilde_e before f, for example).
            D = numel(obj.f);
            if ~isequal(size(obj.mu_tilde), [D, D])
                error('RMTBlocks:InconsistentTypes', ...
                    ['numel(f) = %d but mu_tilde is %s. Use ', ...
                     'set_types(f, mu_tilde, sigma_tilde) to change the number of ', ...
                     'cell types.'], D, mat2str(size(obj.mu_tilde)));
            end
            if ~isequal(size(obj.sigma_tilde), [D, D])
                error('RMTBlocks:InconsistentTypes', ...
                    ['numel(f) = %d but sigma_tilde is %s. Use ', ...
                     'set_types(f, mu_tilde, sigma_tilde) to change the number of ', ...
                     'cell types.'], D, mat2str(size(obj.sigma_tilde)));
            end
        end

        function assert_two_types(obj, name)
            if numel(obj.f) ~= 2
                error('RMTBlocks:NotTwoTypes', ...
                    ['''%s'' is a D = 2 convenience alias, but this object has ', ...
                     'D = %d cell types. Use the D x D properties directly.'], ...
                    name, numel(obj.f));
            end
        end

        function val = uniform_column(obj, mat, col, name)
            % Read a column of a D x D block matrix as a scalar, which is only
            % well-defined when every postsynaptic type sees the same value.
            obj.assert_two_types(name);
            c = mat(:, col);
            if ~all(c == c(1))
                error('RMTBlocks:AmbiguousAlias', ...
                    ['''%s'' is ambiguous: column %d is non-uniform (%s). The scalar ', ...
                     'aliases assume presynaptic-only statistics; read the D x D ', ...
                     'property directly instead.'], name, col, mat2str(c(:).', 4));
            end
            val = c(1);
        end

        function warn_if_zrs_with_blocks(obj)
            % One-shot warning per object.
            if obj.warned_zrs || strcmp(obj.zrs_mode, 'none') || obj.is_column_uniform()
                return;
            end
            obj.warned_zrs = true;
            warning('RMTBlocks:UntestedZRSWithBlocks', ...
                ['zrs_mode ''%s'' has not been validated with non-uniform block ', ...
                 'means/variances. SZRS zeroes all row sums, which no longer removes ', ...
                 'all mean-induced outliers when M has rank > 1 (a second outlier can ', ...
                 'survive); Partial_SZRS preserves block means by construction and is ', ...
                 'the safer choice. lambda_O assumes zrs_mode = ''none''.'], obj.zrs_mode);
        end
    end
end
