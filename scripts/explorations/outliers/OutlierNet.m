classdef OutlierNet < SRNNCellTypePairs
    %OUTLIERNET SRNNCellTypePairs with low-rank structure added to the drawn weight matrix (outlier explorations only).
    %
    %   model = OutlierNet(<the usual name-value pairs>, 'added_structure', s)
    %
    % WHY A SUBCLASS. W and cached_params are protected on the model class,
    % and the outlier exploration needs networks whose weight matrix has
    % structure the generator does not make: an eigenvalue outside the
    % bulk carried by a random pattern of neurons instead of by a whole
    % population. Overriding build_network is the one place where W can be
    % changed before get_params() folds it into everything the dynamics,
    % the Jacobian and the Lyapunov code read, so nothing in the model class
    % is touched. With added_structure empty this is the parent class
    % exactly.
    %
    % THE ADDITION (added_structure, a struct)
    %   a, w    real and imaginary part of the intended eigenvalue, in units
    %           of the structural bulk radius R of the drawn matrix
    %           (obj.R). w = 0 gives rank one, a real eigenvalue.
    %   seed    seed of the private stream the patterns are drawn from
    %           (default rng_seeds(1) + 90001). The global stream is not
    %           touched, so W, the setpoints, the time constants, the
    %           initial state and the noise are those of the network
    %           without the addition.
    % P = U G U' with U an n x 2 orthonormal pair of random Gaussian
    % patterns over ALL neurons (not tied to cell type) and
    % G = R [a -w; w a], so the nonzero eigenvalues of P are R (a +- i w)
    % (for w = 0, P = R a u u'). P is a normal matrix: the left and right
    % eigenvectors of the pair are the same patterns. Added to the bulk,
    % the pair stays close to R (a +- i w) when that is outside the bulk.
    %
    % WHAT IT IS NOT. P is DENSE and IGNORES DALE'S LAW: every neuron gets
    % a small weight of either sign onto every other (about a R / n per
    % entry, a few percent of the standard deviation of an entry of W), so
    % an excitatory neuron's column is no longer one-signed and the
    % connectivity is no longer sparse. W is kept in sparse STORAGE because
    % the class's code assumes that type.
    %
    % See also: outlier_network, outlier_run, SRNNCellTypePairs

    properties
        % struct with fields a, w (and optionally seed), or empty for none
        added_structure = []
    end

    properties (SetAccess = protected)
        % What was added by the last build: struct with U (n x 2 or n x 1),
        % G, target (the intended eigenvalue, units of W), R (the radius
        % used), seed. Empty when nothing was added.
        added = []
    end

    methods
        function obj = OutlierNet(varargin)
            obj = obj@SRNNCellTypePairs(varargin{:});
        end
    end

    methods (Access = protected)
        function build_network(obj)
            build_network@SRNNCellTypePairs(obj);
            obj.added = [];
            s = obj.added_structure;
            if isempty(s); return; end
            seed = obj.rng_seeds(1) + 90001;
            if isfield(s, 'seed') && ~isempty(s.seed); seed = s.seed; end
            stream = RandStream('mt19937ar', 'Seed', seed);
            R = obj.R;
            if s.w == 0
                U = randn(stream, obj.n, 1);
                U = U ./ norm(U);
                G = R * s.a;
            else
                [U, ~] = qr(randn(stream, obj.n, 2), 0);
                G = R * [s.a, -s.w; s.w, s.a];
            end
            obj.W = sparse(full(obj.W) + U * G * U');
            obj.added = struct('U', U, 'G', G, 'target', R * complex(s.a, abs(s.w)), ...
                'R', R, 'seed', seed);
        end
    end
end
