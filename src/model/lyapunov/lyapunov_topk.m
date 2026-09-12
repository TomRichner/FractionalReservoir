function R = lyapunov_topk(X, t, fs, lya_dt, interval, lya_warmup, K, jac_fn, params, opts)
%LYAPUNOV_TOPK The K largest Lyapunov exponents by the discrete QR method.
%   R = lyapunov_topk(X, t, fs, lya_dt, interval, lya_warmup, K, jac_fn, params, opts)
%
%   Propagates an orthonormal N x K basis through the variational equation
%   dY/dt = J(t) Y along a STORED fiducial trajectory, re-orthonormalises it
%   every lya_dt seconds with a thin QR, and accumulates log diag(R) /
%   elapsed time. This is the "diagonal" algorithm of Shimada & Nagashima
%   (1979) and Benettin et al. (1980) -- Skokos (2010) Table 3, Geist et al.
%   (1990) eq. 18, Dieci & Van Vleck (1995) eqs. 2.13-2.16, Engelken, Wolf &
%   Abbott (2023) App. B -- and the one all four sources agree on for K < N:
%   continuous orthonormalisation is numerically delicate once Q Q' ~= I
%   (Dieci sec. 2.3, Skokos eq. 92, Geist Fig. 5).
%
%   Shared by SRNNModel2 and SRNNCellTypePairs, which are duck-typed siblings;
%   each passes its own Jacobian and its own sample-grid helper, so the
%   classes cannot drift apart here. K = N gives the full spectrum by this
%   propagator; the classes' 'qr' method (ode45 on N^2 equations) remains the
%   independent reference it was verified against.
%
%   INPUTS
%     X, t        nt x N fiducial states on the fs grid, and their times
%     fs          sampling rate of X (1/dt)
%     lya_dt      re-orthonormalisation interval (s); an integer multiple of dt
%     interval    [t_acc_start, t_end]: accumulation window (see the classes'
%                 lya_T_interval); iteration starts lya_warmup before it
%     lya_warmup  seconds of alignment before accumulation
%     K           number of vectors, 1 <= K <= N
%     jac_fn      @(S, params) -> N x N sparse Jacobian at state S
%     params      handed to jac_fn unchanged
%     opts        struct: .seed (initial basis), .err_id_prefix (char),
%                 .grid_fn (@(t, dt, decimation, tau, interval, warmup) ->
%                 [idx_all, t_lya, acc_start]; the class's lyapunov_sample_grid),
%                 .jac_times (optional, @(S, Y, params) -> J(S) * Y without
%                 assembling J; when given, jac_fn is never called),
%                 .verbose (logical, default false)
%
%   OUTPUT struct R
%     LE_spectrum            K x 1, sorted descending
%     local_LE_spectrum_t    n_int x K, per-segment rates (sorted columns)
%     finite_LE_spectrum_t   n_int x K, running finite-time exponents; NaN
%                            before accumulation starts
%     t_lya                  n_int x 1 segment start times
%     sort_idx, LLE (= LE_spectrum(1)), K, N
%     h_KS, h_KS_bits        sum of the positive exponents (nats/s, bits/s):
%                            Pesin's identity, valid only if K reaches past
%                            the last positive exponent (see n_positive)
%     n_positive             how many of the K are > 0; if n_positive == K
%                            the positive part may be truncated
%     D_KY, D_KY_resolved    Kaplan-Yorke dimension, and whether the
%                            cumulative sum crossed zero within K (else NaN)
%     cond_t, cond_max       R(1,1)/R(K,K) per segment and its maximum: the
%                            conditioning diagnostic of Engelken et al. eq. B4;
%                            keep it well below 1e8 by shortening lya_dt
%     orth_defect_max        max ||Q'Q - I|| after re-orthonormalisation
%                            (computed when K <= 500)
%     Q_final                N x K basis at the end (a warm start)
%     lya_dt, lya_fs, seconds
%
%   PROPAGATION. Heun (RK2) on the fiducial grid with the Jacobian taken at
%   the STORED states at t_i and t_{i+1}: no interpolants, no adaptive
%   solver. Second order, the same order as the trajectory's own scheme
%   (SRA1's drift is Ralston RK2). Two ways to apply J, same arithmetic:
%     * opts.jac_times given (SRNNCellTypePairs.jacobian_times): two
%       matrix-free products J(S) Y per step, O(nnz(W) K + N K), no matrix.
%     * otherwise (SRNNModel2): one sparse Jacobian assembled per step,
%       reused as the next step's start, plus two sparse x dense products.
%       At N = 4000 the ASSEMBLY (~17 ms) dominates the product (~1 ms), which
%       is why the matrix-free path exists.
%   Per segment: one N x K QR, O(N K^2).
%
%   SIGN CONVENTION. MATLAB's qr does not make diag(R) positive; the columns
%   of Q and rows of R are flipped so it is (Geist A.6, Dieci 2.15). Without
%   that, log(diag(R)) is undefined and Q's orientation would drift.
%
%   INITIAL BASIS. Householder QR of a Gaussian N x K draw, seeded from
%   opts.seed with the twister generator under a saved/restored RNG state.
%   Random, not [e_1 .. e_K]: the coordinate basis can pin a trivial exponent
%   to the wrong slot (Geist sec. 2.2). randn fills column-major and QR is
%   column-sequential, so the first K columns of a K' > K run with the same
%   seed start identical, evolve identically (the propagation is linear and
%   column-wise) and are re-orthonormalised identically: spectra are NESTED,
%   and test_lyapunov_topk asserts it to round-off.
%
%   See also SRNNCellTypePairs, SRNNModel2, test_lyapunov_topk,
%   docs/notes/Lyapunov_estimation_methods.md

    t0_wall = tic;
    pfx = 'lyapunov_topk';
    if isfield(opts, 'err_id_prefix') && ~isempty(opts.err_id_prefix)
        pfx = opts.err_id_prefix;
    end
    verbose = isfield(opts, 'verbose') && opts.verbose;
    matrix_free = isfield(opts, 'jac_times') && ~isempty(opts.jac_times);
    if matrix_free
        jac_times = opts.jac_times;
    end

    N  = size(X, 2);
    nt = size(X, 1);
    if ~isscalar(K) || ~isnumeric(K) || ~isfinite(K) || K ~= round(K) || K < 1 || K > N
        error([pfx ':InvalidLyapunovK'], ...
            'lya_K must be an integer in [1, %d] (0 selects N); got %s.', N, mat2str(K));
    end
    if ~isfield(opts, 'grid_fn') || isempty(opts.grid_fn)
        error([pfx ':InvalidParams'], 'opts.grid_fn (the sample-grid helper) is required.');
    end

    dt         = 1 / fs;
    decimation = round(lya_dt * fs);
    tau        = decimation * dt;
    tol        = dt * 1e-6;
    [idx_all, t_lya, acc_start] = opts.grid_fn(t, dt, decimation, tau, interval, lya_warmup);
    n_int = numel(t_lya);

    % Initial orthonormal basis, seeded and generator-pinned; RNG restored.
    seed = 0;
    if isfield(opts, 'seed') && ~isempty(opts.seed); seed = opts.seed; end
    stream_state = rng;
    rng(seed, 'twister');
    G = randn(N, K);             % column-major fill: the first K columns of the N x N draw
    rng(stream_state);
    [Q, R0] = qr(G, 0);          % Householder QR is column-sequential, so Q(:,1:K) is
    s0 = sign(diag(R0)); s0(s0 == 0) = 1;   % the same for every K' >= K (nesting)
    Q = Q .* s0';

    local_spec  = zeros(n_int, K);
    finite_spec = nan(n_int, K);
    cond_t      = nan(n_int, 1);
    accumulated = zeros(K, 1);
    elapsed     = 0;
    orth_defect_max = 0;
    check_orth  = K <= 500;

    J1 = [];                     % Jacobian at the segment's first grid point
    i_prev_end = -1;
    for k = 1:n_int
        i0 = idx_all(k);
        i1 = i0 + decimation;
        if i1 > nt; break; end
        Y = Q;
        if matrix_free
            for i = i0:i1 - 1
                h  = t(i + 1) - t(i);
                F0 = jac_times(X(i, :)', Y, params);
                Yp = Y + h * F0;
                Y  = Y + (h / 2) * (F0 + jac_times(X(i + 1, :)', Yp, params));
            end
        else
            if i0 ~= i_prev_end || isempty(J1)
                J1 = jac_fn(X(i0, :)', params);
            end
            for i = i0:i1 - 1
                J0 = J1;
                J1 = jac_fn(X(i + 1, :)', params);
                h  = t(i + 1) - t(i);
                F0 = J0 * Y;
                Yp = Y + h * F0;
                Y  = Y + (h / 2) * (F0 + J1 * Yp);
            end
            i_prev_end = i1;
        end

        if any(~isfinite(Y(:)))
            warning([pfx ':LyapunovDiverged'], ...
                'Top-K trajectory diverged at t = %g.', t_lya(k));
            t_lya(k:end) = []; local_spec(k:end, :) = []; finite_spec(k:end, :) = []; cond_t(k:end) = [];
            break;
        end

        [Q, Rk] = qr(Y, 0);
        s = sign(diag(Rk)); s(s == 0) = 1;
        Q  = Q .* s';
        d  = abs(diag(Rk));
        d  = max(d, realmin);
        if check_orth
            orth_defect_max = max(orth_defect_max, norm(Q' * Q - eye(K)));
        end
        cond_t(k) = d(1) / d(end);

        seg = t(i1) - t(i0);
        log_d = log(d);
        local_spec(k, :) = (log_d / seg)';
        if t_lya(k) >= acc_start - tol && t_lya(k) < interval(2)
            accumulated = accumulated + log_d;
            elapsed     = elapsed + seg;
            finite_spec(k, :) = (accumulated / elapsed)';
        end
        if verbose && (k == 1 || mod(k, 50) == 0)
            fprintf('  [topk] segment %d/%d t = %.2f  lambda_1 = %+.4f  cond %.2e\n', ...
                k, n_int, t_lya(k), accumulated(1) / max(elapsed, eps), cond_t(k));
        end
    end

    if elapsed > 0
        spectrum = accumulated / elapsed;
    else
        warning([pfx ':LyapunovWindowEmpty'], ...
            'No top-K segment fell inside lya_T_interval = [%g, %g]; the spectrum is all NaN.', ...
            interval(1), interval(2));
        spectrum = nan(K, 1);
    end

    [spectrum, order] = sort(spectrum, 'descend');
    R = struct();
    R.LE_spectrum          = spectrum;
    R.local_LE_spectrum_t  = local_spec(:, order);
    R.finite_LE_spectrum_t = finite_spec(:, order);
    R.t_lya    = t_lya;
    R.sort_idx = order;
    R.LLE      = spectrum(1);
    R.K = K;  R.N = N;
    pos = spectrum(spectrum > 0);
    R.n_positive = numel(pos);
    R.h_KS       = sum(pos);
    R.h_KS_bits  = R.h_KS / log(2);
    [R.D_KY, R.D_KY_resolved] = kaplan_yorke_partial(spectrum, N);
    R.cond_t   = cond_t;
    R.cond_max = max(cond_t, [], 'omitnan');
    R.orth_defect_max = orth_defect_max;
    R.Q_final  = Q(:, order);
    R.lya_dt   = lya_dt;
    R.lya_fs   = 1 / lya_dt;
    R.seconds  = toc(t0_wall);
end

function [D, resolved] = kaplan_yorke_partial(lambda, N)
% Kaplan-Yorke dimension from a (possibly partial) descending spectrum.
% Resolved only if the cumulative sum crosses zero within the K given, or
% if K = N and it never does (then D = N).
    sums = cumsum(lambda);
    K = numel(lambda);
    j = find(sums >= 0, 1, 'last');
    if isempty(j)
        D = 0; resolved = true;                     % lambda_1 < 0: a point attractor
    elseif j == K
        if K == N
            D = N; resolved = true;                 % volume never contracts: not dissipative
        else
            D = NaN; resolved = false;              % crossing lies beyond the K tracked
        end
    else
        D = j + sums(j) / abs(lambda(j + 1)); resolved = true;
    end
end
