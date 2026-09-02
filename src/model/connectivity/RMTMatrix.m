classdef RMTMatrix < handle
    % RMTMatrix - Random Matrix Theory class for W matrix construction
    % Based on Harris et al. (2023). Simplified version of RMT4 for SRNNModel.

    properties
        N               % System size
        alpha           % Sparsity/connection probability (0 < alpha <= 1)
        f               % Fraction of excitatory neurons

        % Normalized population statistics (tilde notation from Harris 2023)
        mu_tilde_e      % Normalized mean of excitatory population
        mu_tilde_i      % Normalized mean of inhibitory population
        sigma_tilde_e   % Normalized std dev of excitatory population
        sigma_tilde_i   % Normalized std dev of inhibitory population

        % Internal matrices
        A               % Base random matrix (Gaussian, mean 0, var 1)
        S               % Sparsity mask (logical)

        % Control flags
        zrs_mode        % 'none', 'ZRS', 'SZRS', 'Partial_SZRS'
    end

    properties (Dependent)
        % Population indices
        E               % Logical index for Excitatory neurons
        I               % Logical index for Inhibitory neurons

        % Low-rank structure M = u * v' (Eq 12)
        u               % Left vector: ones(N,1)
        v               % Right vector: population means

        % Variance structure (Eq 11)
        D               % Diagonal variance matrix

        % Weight/Jacobian matrix
        W               % Jacobian matrix (computed on access)

        % Sparse statistics (Eq 15, 16)
        mu_se           % Sparse excitatory mean
        mu_si           % Sparse inhibitory mean
        sigma_se_sq     % Sparse excitatory variance
        sigma_si_sq     % Sparse inhibitory variance

        % Theoretical predictions (Eq 17, 18)
        lambda_O        % Outlier eigenvalue
        R               % Spectral radius
    end

    methods
        function obj = RMTMatrix(N)
            % RMTMatrix Constructor
            obj.N = N;

            % Defaults
            obj.alpha = 1.0;
            obj.f = 0.5;
            obj.mu_tilde_e = 0;
            obj.mu_tilde_i = 0;
            obj.sigma_tilde_e = 1/sqrt(N);
            obj.sigma_tilde_i = 1/sqrt(N);
            obj.zrs_mode = 'none';

            % Initialize random matrices
            obj.A = randn(N, N);
            obj.update_sparsity();
        end

        %% Dependent Property Getters
        function val = get.E(obj)
            val = false(obj.N, 1);
            val(1:round(obj.f * obj.N)) = true;
        end

        function val = get.I(obj)
            val = ~obj.E;
        end

        function val = get.u(obj)
            val = ones(obj.N, 1);
        end

        function val = get.v(obj)
            val = zeros(obj.N, 1);
            E_idx = obj.E;
            val(E_idx) = obj.mu_tilde_e;
            val(~E_idx) = obj.mu_tilde_i;
        end

        function val = get.D(obj)
            D_vec = zeros(obj.N, 1);
            E_idx = obj.E;
            D_vec(E_idx) = obj.sigma_tilde_e;
            D_vec(~E_idx) = obj.sigma_tilde_i;
            val = diag(D_vec);
        end

        function val = get.W(obj)
            % Construct weight matrix W based on Harris 2023 equations
            D = obj.D;
            M = obj.u * obj.v';

            switch obj.zrs_mode
                case 'none'
                    W_dense = (obj.A * D) + M;
                    val = obj.S .* W_dense;

                case 'ZRS'
                    if obj.alpha < 1
                        warning('RMTMatrix:SparsityWarning', ...
                            'Using ZRS with sparse matrix destroys sparsity. Consider SZRS.');
                    end
                    P = eye(obj.N) - (obj.u * obj.u') / obj.N;
                    val = (obj.A * D * P) + M;
                    if obj.alpha < 1
                        val = obj.S .* val;
                    end

                case 'SZRS'
                    W_base = obj.S .* ((obj.A * D) + M);
                    row_sums = sum(W_base, 2);
                    row_counts = sum(obj.S, 2);
                    row_counts(row_counts == 0) = 1;
                    W_bar_i = row_sums ./ row_counts;
                    B = obj.S .* W_bar_i;
                    val = W_base - B;

                case 'Partial_SZRS'
                    J_base = obj.S .* (obj.A * D);
                    M_base = obj.S .* M;
                    J_row_sums = sum(J_base, 2);
                    row_counts = sum(obj.S, 2);
                    row_counts(row_counts == 0) = 1;
                    J_bar_i = J_row_sums ./ row_counts;
                    B_partial = obj.S .* J_bar_i;
                    val = (J_base - B_partial) + M_base;
            end
        end

        function val = get.mu_se(obj)
            val = obj.alpha * obj.mu_tilde_e;
        end

        function val = get.mu_si(obj)
            val = obj.alpha * obj.mu_tilde_i;
        end

        function val = get.sigma_se_sq(obj)
            val = obj.alpha * (1 - obj.alpha) * obj.mu_tilde_e^2 + obj.alpha * obj.sigma_tilde_e^2;
        end

        function val = get.sigma_si_sq(obj)
            val = obj.alpha * (1 - obj.alpha) * obj.mu_tilde_i^2 + obj.alpha * obj.sigma_tilde_i^2;
        end

        function val = get.lambda_O(obj)
            val = obj.N * (obj.f * obj.mu_se + (1 - obj.f) * obj.mu_si);
        end

        function val = get.R(obj)
            val = sqrt(obj.N * (obj.f * obj.sigma_se_sq + (1 - obj.f) * obj.sigma_si_sq));
        end

        %% Property Setters
        function set.alpha(obj, val)
            obj.alpha = val;
            if ~isempty(obj.A)
                obj.update_sparsity();
            end
        end

        %% Internal Updates
        function update_sparsity(obj)
            obj.S = rand(obj.N, obj.N) < obj.alpha;
        end

        %% Convenience setters
        function set_params(obj, mu_tilde_e, mu_tilde_i, sigma_tilde_e, sigma_tilde_i, f, alpha)
            if nargin > 1, obj.mu_tilde_e = mu_tilde_e; end
            if nargin > 2, obj.mu_tilde_i = mu_tilde_i; end
            if nargin > 3, obj.sigma_tilde_e = sigma_tilde_e; end
            if nargin > 4, obj.sigma_tilde_i = sigma_tilde_i; end
            if nargin > 5, obj.f = f; end
            if nargin > 6, obj.alpha = alpha; end
        end

        function set_zrs_mode(obj, mode)
            valid_modes = {'none', 'ZRS', 'SZRS', 'Partial_SZRS'};
            if ~ismember(mode, valid_modes)
                error('Invalid ZRS mode. Valid: %s', strjoin(valid_modes, ', '));
            end
            obj.zrs_mode = mode;
        end
    end
end
