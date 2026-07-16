classdef RMTCellTypes < handle
    %RMTCELLTYPES Sparse random connectivity for arbitrary cell types.
    %
    % Columns of W are presynaptic neurons. For a neuron belonging to cell
    % type q, each potential connection is present with probability alpha and
    % its weight is drawn from
    %
    %   mu_tilde(q) + sigma_tilde(q) * randn().
    %
    % Usage:
    %   rmt = RMTCellTypes(N, alpha, f, mu_tilde, sigma_tilde);
    %   W = rmt.W;

    properties (SetAccess = private)
        N
        alpha
        f
        mu_tilde
        sigma_tilde
        n_cellTypes
        n_per_type
        type_indices
        W
    end

    methods
        function obj = RMTCellTypes(N, alpha, f, mu_tilde, sigma_tilde)
            if nargin ~= 5
                error('RMTCellTypes:InvalidInput', ...
                    'Expected RMTCellTypes(N, alpha, f, mu_tilde, sigma_tilde).');
            end

            obj.N = N;
            obj.alpha = alpha;
            obj.f = reshape(f, 1, []);
            obj.mu_tilde = reshape(mu_tilde, 1, []);
            obj.sigma_tilde = reshape(sigma_tilde, 1, []);
            obj.n_cellTypes = numel(obj.f);
            obj.validate();

            obj.n_per_type = RMTCellTypes.allocate_counts(obj.N, obj.f);
            obj.type_indices = RMTCellTypes.make_type_indices(obj.n_per_type);
            obj.W = obj.generate_W();
        end
    end

    methods (Access = private)
        function validate(obj)
            if ~isscalar(obj.N) || ~isnumeric(obj.N) || ~isfinite(obj.N) || ...
                    obj.N < 1 || obj.N ~= round(obj.N)
                error('RMTCellTypes:InvalidParams', 'N must be a positive integer.');
            end
            if ~isscalar(obj.alpha) || ~isnumeric(obj.alpha) || ...
                    ~isfinite(obj.alpha) || obj.alpha <= 0 || obj.alpha > 1
                error('RMTCellTypes:InvalidParams', 'alpha must satisfy 0 < alpha <= 1.');
            end
            if obj.n_cellTypes < 1 || any(~isfinite(obj.f)) || any(obj.f <= 0)
                error('RMTCellTypes:InvalidParams', 'f must contain positive finite fractions.');
            end
            if abs(sum(obj.f) - 1) > 1e-12
                error('RMTCellTypes:InvalidParams', 'f must sum to 1.');
            end
            if obj.N < obj.n_cellTypes
                error('RMTCellTypes:InvalidParams', ...
                    'N must be at least the number of cell types.');
            end
            if numel(obj.mu_tilde) ~= obj.n_cellTypes || ...
                    any(~isfinite(obj.mu_tilde))
                error('RMTCellTypes:InvalidParams', ...
                    'mu_tilde must have one finite value per cell type.');
            end
            if numel(obj.sigma_tilde) ~= obj.n_cellTypes || ...
                    any(~isfinite(obj.sigma_tilde)) || any(obj.sigma_tilde < 0)
                error('RMTCellTypes:InvalidParams', ...
                    'sigma_tilde must have one nonnegative finite value per cell type.');
            end
        end

        function W = generate_W(obj)
            type_id = repelem(1:obj.n_cellTypes, obj.n_per_type);
            mu_by_column = obj.mu_tilde(type_id);
            sigma_by_column = obj.sigma_tilde(type_id);

            mask = rand(obj.N, obj.N) < obj.alpha;
            weights = randn(obj.N, obj.N) .* sigma_by_column + mu_by_column;
            W = sparse(weights .* mask);
        end
    end

    methods (Static)
        function counts = allocate_counts(N, f)
            %ALLOCATE_COUNTS Convert fractions to counts by largest remainder.
            f = reshape(f, 1, []);
            raw_counts = N .* f;
            counts = floor(raw_counts);
            remainder = N - sum(counts);

            if remainder > 0
                fractional = raw_counts - counts;
                ordering = sortrows([-fractional(:), (1:numel(f))'], [1, 2]);
                recipients = ordering(1:remainder, 2);
                counts(recipients) = counts(recipients) + 1;
            end

            if any(counts < 1)
                error('RMTCellTypes:InvalidParams', ...
                    'Every cell type must receive at least one neuron.');
            end
        end

        function indices = make_type_indices(counts)
            %MAKE_TYPE_INDICES Return contiguous neuron indices by type.
            indices = cell(1, numel(counts));
            first_index = 1;
            for q = 1:numel(counts)
                last_index = first_index + counts(q) - 1;
                indices{q} = first_index:last_index;
                first_index = last_index + 1;
            end
        end
    end
end
