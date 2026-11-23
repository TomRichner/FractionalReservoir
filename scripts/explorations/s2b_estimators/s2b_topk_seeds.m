function C = s2b_topk_seeds(t, S, fs, params, interval, warmup, seeds, K, lya_dt)
%S2B_TOPK_SEEDS The estimator of the paper runs (lyapunov_topk, as
% SRNNCellTypePairs calls it: K vectors, QR every lya_dt, the Jacobian along
% the STORED trajectory) repeated with different seeds of the initial basis,
% on one trajectory. interval = accumulation window, warmup = seconds of
% alignment before it. Returns per seed the leading local rate the figures
% use (local_LE_spectrum_t(:, 1) inside the window: the column with the
% largest window average), which raw column that is, the spectrum, and the
% leading vector at the START and the END of the window (start: a second
% pass over the warm-up alone, time shifted so that the same basis is
% iterated).
opts = struct('err_id_prefix', 's2b', 'verbose', false, 'grid_fn', @grid, ...
    'jac_times', @(St, Y, p) SRNNCellTypePairs.jacobian_times(St, Y, p));
jac = @(St, p) SRNNCellTypePairs.compute_Jacobian_fast(St, p);
ts = interval(1) - warmup;
rows = t >= ts - 1e-9 & t <= interval(1) + 1e-9;
C = struct('seed', {}, 'LLE', {}, 'spec', {}, 't', {}, 'loc', {}, 'raw_col', {}, 'q_end', {}, 'q_start', {}, 'seconds', {});
for k = 1:numel(seeds)
    opts.seed = seeds(k);
    t0 = tic;
    R = lyapunov_topk(S, t, fs, lya_dt, interval, warmup, K, jac, params, opts);
    in = isfinite(R.finite_LE_spectrum_t(:, 1));
    c = struct('seed', seeds(k), 'LLE', R.LLE, 'spec', R.LE_spectrum(:)', 't', R.t_lya(in), ...
        'loc', R.local_LE_spectrum_t(in, 1), 'raw_col', R.sort_idx(1), 'q_end', R.Q_final(:, 1), ...
        'q_start', [], 'seconds', 0);
    if warmup > 0
        W = lyapunov_topk(S(rows, :), t(rows) - ts, fs, lya_dt, [0, warmup], 0, K, jac, params, opts);
        Qraw = zeros(size(W.Q_final));  Qraw(:, W.sort_idx) = W.Q_final;
        c.q_start = Qraw(:, c.raw_col);
    end
    c.seconds = toc(t0);
    C(k) = c;
end
end

function [idx_all, t_lya, acc_start] = grid(t, dt, decimation, tau, interval, lya_warmup)
% SRNNCellTypePairs.lyapunov_sample_grid (protected there), without its warning.
tol = dt * 1e-6;
acc_start = max(0, interval(1));
i_start = max(1, min(round((acc_start - lya_warmup - t(1)) / dt) + 1, numel(t)));
idx_all = i_start:decimation:numel(t);
t_lya = t(idx_all);  t_lya = t_lya(:);
keep = (idx_all + decimation <= numel(t)) & (t_lya(:).' + tau <= min(t(end), interval(2)) + tol);
idx_all = idx_all(keep);  t_lya = t_lya(keep);
end
