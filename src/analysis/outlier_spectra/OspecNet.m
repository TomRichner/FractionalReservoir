classdef OspecNet < SRNNCellTypePairs
    %OSPECNET SRNNCellTypePairs for the outlier_spectra stage: optional low-rank structure added to W, and a run that keeps only the sampled states.
    %
    %   model = OspecNet(<the usual name-value pairs>, 'added_structure', s)
    %   model.build();
    %   [S_keep, lya] = model.run_sampled(keep_idx, chunk_s)
    %
    % Ported from scripts/explorations/outliers/OutlierNet.m (FD-040), which
    % stays where it is; the addition to W is the same code. run_sampled is
    % new. With added_structure empty the model is the parent class's,
    % weights and all.
    %
    % WHY A SUBCLASS. W, the noise tensor and the trajectory are protected
    % on the model class. Overriding build_network is the one place where W
    % can be changed before get_params() folds it into everything the
    % dynamics, the Jacobian and the Lyapunov code read; and a method of a
    % subclass can integrate the model's own equations with the model's
    % own noise without holding the whole trajectory. Nothing in the model
    % class is touched.
    %
    % RUN_SAMPLED: THE SAME TRAJECTORY AS run(), WITHOUT HOLDING IT.
    % run() keeps every step of every state variable until it returns: at
    % 260 s, 500 Hz and 4000 state variables 4.2 GB, per pool worker.
    % run_sampled integrates the same time grid (t_ex) in pieces of
    % chunk_s seconds through the same integrator (obj.integrate) against
    % the same pre-generated noise (obj.build_noise), each piece started
    % from the last row of the one before, and keeps only the rows asked
    % for. The fixed-step integrators address the noise by absolute time
    % and carry no state between steps but y, so the rows are those of
    % run() BIT FOR BIT (test_outlier_spectra_stage asserts it). It needs a
    % fixed-step integrator (rk4, or the stochastic schemes): an adaptive
    % one would choose other steps in each piece.
    %
    % The Lyapunov exponents, when lya_method is not 'none', are the class's
    % own (compute_lyapunov, then lya_summary) on the stretch of trajectory
    % they need and no more: from lya_warmup before lya_T_interval(1) to
    % lya_T_interval(2), which must be set. That stretch IS held at full
    % rate while they are computed (25 s at 500 Hz and 4000 variables:
    % 0.4 GB). The top-K method reads the stored states on a grid laid out
    % from where the iteration starts, so the result is the one run() gives
    % for the same interval (asserted by the same test).
    %
    % WHAT IS STILL HELD for the whole run: the noise (two n x steps
    % arrays, 1.0 GB at 260 s, 500 Hz, 500 neurons) and the external input
    % with its interpolant (n x steps each, another 1.0 GB), both made by
    % the parent class.
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
    % (for w = 0, P = R a u u'). P is DENSE and IGNORES DALE'S LAW; W is
    % kept in sparse STORAGE because the class's code assumes that type.
    %
    % See also: ospec_network, ospec_run, ospec_rebuild_model, SRNNCellTypePairs

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
        function obj = OspecNet(varargin)
            obj = obj@SRNNCellTypePairs(varargin{:});
        end

        function [S_keep, lya] = run_sampled(obj, keep_idx, chunk_s)
            %RUN_SAMPLED Integrate the whole time grid; return the states at rows keep_idx of it (N x numel(keep_idx)) and the Lyapunov summary.
            %
            %   keep_idx  ascending row numbers of obj.t_ex
            %   chunk_s   seconds integrated per piece (default 10)
            %   S_keep    N x numel(keep_idx): the state at t_ex(keep_idx)
            %   lya       lya_summary() of the class, or an empty struct
            %             when lya_method is 'none'
            if nargin < 3 || isempty(chunk_s); chunk_s = 10; end
            if ~obj.is_built
                error('OspecNet:NotBuilt', 'Model must be built before running. Call build() first.');
            end
            if ~ismember(obj.ode_solver, [{'rk4'}, stochastic_solver_names()])
                error('OspecNet:Solver', ...
                    'run_sampled integrates in pieces and needs a fixed-step integrator; ode_solver is ''%s''.', ...
                    obj.ode_solver);
            end
            nt = numel(obj.t_ex);
            keep_idx = keep_idx(:)';
            if isempty(keep_idx) || any(diff(keep_idx) <= 0) || keep_idx(1) < 1 || keep_idx(end) > nt ...
                    || any(keep_idx ~= round(keep_idx))
                error('OspecNet:KeepIdx', 'keep_idx must be ascending row numbers of t_ex (1 to %d).', nt);
            end

            params = obj.cached_params;
            params.u_interpolant = obj.u_interpolant;
            N = params.N_sys_eqs;
            rhs = @(t, S) SRNNCellTypePairs.dynamics_fast(t, S, params);
            obj.build_noise();

            % The stretch the Lyapunov estimate reads, at full rate.
            want_lya = ~strcmpi(obj.lya_method, 'none');
            l1 = 0; l2 = -1; S_lya = [];
            if want_lya
                if isempty(obj.lya_T_interval)
                    error('OspecNet:LyaInterval', 'run_sampled needs lya_T_interval set when lya_method is not ''none''.');
                end
                t_from = max(obj.t_ex(1), obj.lya_T_interval(1) - obj.lya_warmup);
                l1 = round((t_from - obj.t_ex(1)) * obj.fs) + 1;
                l2 = min(nt, round((obj.lya_T_interval(2) - obj.t_ex(1)) * obj.fs) + 1);
                S_lya = zeros(l2 - l1 + 1, N);
            end

            S_keep = zeros(N, numel(keep_idx));
            per = max(2, round(chunk_s * obj.fs));
            S = obj.S0(:);
            i0 = 1;
            if keep_idx(1) == 1; S_keep(:, 1) = S; end
            if want_lya && l1 == 1; S_lya(1, :) = S.'; end
            while i0 < nt
                i1 = min(nt, i0 + per);
                if nt - i1 < 2; i1 = nt; end         % the integrator needs three points: no short last piece
                [~, Y] = obj.integrate(rhs, obj.t_ex(i0:i1), S);
                % Rows i0+1 .. i1 are new in this piece (row i0 was the last of the one before).
                sel = keep_idx > i0 & keep_idx <= i1;
                S_keep(:, sel) = Y(keep_idx(sel) - i0 + 1, :).';
                a = max(l1, i0 + 1); b = min(l2, i1);
                if a <= b
                    S_lya(a - l1 + 1:b - l1 + 1, :) = Y(a - i0 + 1:b - i0 + 1, :);
                end
                S = Y(end, :).';
                i0 = i1;
            end
            clear Y

            lya = struct();
            if want_lya
                obj.t_out = obj.t_ex(l1:l2);
                obj.S_out = S_lya;
                clear S_lya
                obj.compute_lyapunov();
                lya = obj.lya_summary();
                obj.S_out = [];
                obj.t_out = [];
            end
            obj.noise_increments = [];
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
