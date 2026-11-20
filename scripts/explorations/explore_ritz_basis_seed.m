function out = explore_ritz_basis_seed(cfg)
% EXPLORE_RITZ_BASIS_SEED How much do the Ritz values depend on the initial basis? (FD-033)
%
%   out = EXPLORE_RITZ_BASIS_SEED()
%   out = EXPLORE_RITZ_BASIS_SEED('K', 50, 'seeds', [424243 7])
%
% The Ritz overlay of Figure 3 draws the eigenvalues of the Jacobian
% projected on the top-K Lyapunov subspace. Ritz values depend only on the
% subspace, but in a stable adapting network the top-K subspace has its
% edge inside a nearly continuous band of slow modes, where the gap between
% neighboring exponents is far too small for any affordable alignment time
% to settle which directions are in. This runs the eig-heatmap stage's
% reference trajectory at its 'essential' / 'fast' size (T 40 s, 200 Hz,
% seeds [1 2], states sampled as the stage samples them) for each condition
% and makes the Ritz pass from two initial bases, in the raw and the
% type_sd norm, and reports per condition:
%   * the largest principal angle between the two tracked subspaces, and how
%     many of the K principal angles are below 1 degree (the part of the
%     subspace both bases agree on),
%   * per state, the largest distance between a Ritz value and its partner
%     in the other basis's set (greedy pairing): median and max over states,
%   * the same summary numbers the figure's table gives, for both bases
%     (share with Re > 0, median Re, share complex), so one can see whether
%     the picture, as opposed to single points, depends on the basis.
% No eigenvalues of the full Jacobian are computed, so it takes about a
% minute per condition. Nothing is saved.
%
% See also: run_eig_heatmap, fig_eig_heatmap, SRNNCellTypePairs.projected_jacobian

arguments
    cfg.preset_name (1,:) char = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain'
    cfg.conditions = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'}
    cfg.K (1,1) double = 50
    cfg.seeds (1,2) double = [424243, 7]      % the first is the stage's (rng_seeds(1) + 424242)
    cfg.n_samples (1,1) double = 40
    cfg.T (1,1) double = 40
    cfg.fs (1,1) double = 200
end
setup_paths();
window = [cfg.T / 2, cfg.T];
scalings = {'raw', 'type_sd'};
out = struct('name', cfg.conditions, 'rows', []);
for c = 1:numel(cfg.conditions)
    m = build_from_preset(cfg.preset_name, cfg.conditions{c}, 'T_range', [0 cfg.T], 'fs', cfg.fs, ...
        'rng_seeds', [1 2], 'lya_method', 'none', 'store_full_state', true, 'verbose', 'near-none');
    m.run();
    params = m.get_params();
    J_times = linspace(window(1), cfg.T, cfg.n_samples);
    idx = arrayfun(@(tt) find(m.t_out >= tt, 1, 'first'), J_times);
    idx = unique(idx);
    scales = SRNNCellTypePairs.state_type_scales(m.S_out, min(idx):numel(m.t_out), params);
    fprintf('\n%s (state dimension %d, K = %d, %d states)\n', cfg.conditions{c}, params.state_layout.N, ...
        min(cfg.K, params.state_layout.N), numel(idx));
    rows = struct('scaling', scalings, 'angle_max_deg', [], 'n_below_1deg', [], 'moved_median', [], ...
        'moved_max', [], 'share_pos', [], 'median_re', [], 'share_complex', [], 'seconds', []);
    for s = 1:numel(scalings)
        Z = cell(1, 2); Qend = cell(1, 2); secs = zeros(1, 2);
        for b = 1:2
            % The pass, with the basis at the last sampled state kept for the angles.
            o = struct('seed', cfg.seeds(b), 'err_id_prefix', 'explore_ritz_basis_seed', ...
                'grid_fn', @plain_grid, 'H_idx', idx(:), 'H_return_Q', true);
            if strcmp(scalings{s}, 'type_sd')
                d = scales.d;
                o.jac_times = @(S, Y, p) SRNNCellTypePairs.jacobian_times(S, d .* Y, p) ./ d;
            else
                o.jac_times = @(S, Y, p) SRNNCellTypePairs.jacobian_times(S, Y, p);
            end
            K = min(cfg.K, params.state_layout.N);
            r = lyapunov_topk(m.S_out, m.t_out, m.fs, 0.05, [m.t_out(idx(1)), m.t_out(end)], ...
                m.t_out(idx(1)) - m.t_out(1), K, @(S, p) error('unused'), params, o);
            secs(b) = r.seconds;
            Z{b} = complex(zeros(K, numel(idx)));
            for k = 1:numel(idx); Z{b}(:, k) = eig(r.H_t(:, :, k)); end
            Qend{b} = r.Q_t(:, :, end);
        end
        ang = acosd(min(1, svd(Qend{1}' * Qend{2})));           % principal angles at the last state
        moved = arrayfun(@(k) matched_distance(Z{1}(:, k), Z{2}(:, k)), 1:numel(idx));
        cx = @(z) mean(abs(imag(z(:))) > 1e-9 * max(1, abs(z(:))));
        rows(s).angle_max_deg = max(ang); rows(s).n_below_1deg = nnz(ang < 1);
        rows(s).moved_median = median(moved); rows(s).moved_max = max(moved);
        rows(s).share_pos = [mean(real(Z{1}(:)) > 0), mean(real(Z{2}(:)) > 0)];
        rows(s).median_re = [median(real(Z{1}(:))), median(real(Z{2}(:)))];
        rows(s).share_complex = [cx(Z{1}), cx(Z{2})];
        rows(s).seconds = secs;
        fprintf(['  %-7s largest principal angle %5.1f deg, %2d of %d below 1 deg | Ritz values moved: median %.3g, ' ...
            'max %.3g 1/s | Re > 0: %.1f%% / %.1f%% | median Re %+.3f / %+.3f | complex %.0f%% / %.0f%% | %.0f + %.0f s\n'], ...
            scalings{s}, max(ang), nnz(ang < 1), K, median(moved), max(moved), 100 * rows(s).share_pos, ...
            rows(s).median_re, 100 * rows(s).share_complex, secs);
    end
    out(c).rows = rows;
end
end

function [idx_all, t_lya, acc_start] = plain_grid(t, dt, decimation, tau, interval, lya_warmup)
% The classes' segment grid (lyapunov_sample_grid is protected): start
% lya_warmup before the window, whole segments only.
acc_start = max(0, interval(1));
i_start = max(1, round((acc_start - lya_warmup - t(1)) / dt) + 1);
idx_all = i_start:decimation:numel(t);
t_lya = t(idx_all); t_lya = t_lya(:);
keep = (idx_all + decimation <= numel(t)) & (t_lya(:).' + tau <= min(t(end), interval(2)) + dt * 1e-6);
idx_all = idx_all(keep); t_lya = t_lya(keep);
end

function d = matched_distance(a, b)
% Largest distance after pairing two sets of complex numbers greedily.
a = a(:); b = b(:); d = 0;
for k = 1:numel(a)
    [dmin, j] = min(abs(b - a(k)));
    d = max(d, dmin); b(j) = [];
end
end
