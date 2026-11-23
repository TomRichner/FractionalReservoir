function R = soied_run_network(seed, S)
% SOIED_RUN_NETWORK One IED network (Fig 8's model), a baseline block and a stimulation block: discharges, frozen-gain spectra, full-Jacobian network mode.
%
%   R = soied_run_network(seed, S)
%
% seed   network seed of the ied_stim bundle (42, 73, 101, ...)
% S      settings (soied_run): block_len, settle, fs, amplitude,
%        coarse_dt, dense_dt, dense_win, full_dt, n_top, log
%
% THE NETWORK AND BLOCKS are the ied_stim bundle's (ied_stim_blocks,
% run_block): ied_base_config with fs = S.fs, ied_build_network, a
% 10-s model warmup, S.settle s of settling, then S.block_len s analysed.
% Block 1: uniform input 0; block 2: uniform input S.amplitude, constant
% through the block (the 'essential' mode's numbering, n_no = 1). Noise
% seed base.noise_seed + (seed - 42) + 1000*block, as the bundle. The
% detector is ied_detect_events on plot_data over the analysed window,
% exactly as run_block calls it.
%
% FROZEN-GAIN SPECTRA. M = W diag(d) = I + tau_d J_xx (ospec_effective_matrix)
% at every S.coarse_dt s of the analysed window and every S.dense_dt s in
% [t_ev + dense_win(1), t_ev + dense_win(2)] around each detected event.
% Per sample, all eigenvalues (single) and for the S.n_top rightmost
% (Im >= 0): coherence coh10 = norm^2 fraction of the unit eigenvector in
% the span of the ten population indicators (E and I of each subnet; a
% random vector has 10/350 = 0.029), coh2 (the E/I span), and the share of
% |v|^2 on each subnet. Bulk radius predicted from the gains on the 2-type
% blocks (ospec R_eff) and on the 10-population blocks (R10), the
% eigenvalues of the 10 x 10 gain-weighted mean matrix K10, mean gains
% (E/I x E/I) and mean rates.
%
% FULL JACOBIAN every S.full_dt s: dense eig with vectors of the
% 1190 x 1190 matrix (ospec_full_spectrum goes to eigs above N = 1000,
% which does not converge here); network modes have >= 10% of the
% eigenvector's norm^2 in the dendritic rows.
%
% See also: soied_run, soied_measures, soied_main_figure, ied_stim_blocks

base = ied_base_config(); base.fs = S.fs;
R = struct('seed', seed, 'settings', S);
for b = 1:2
    t0 = tic;
    A = S.amplitude * (b == 2);
    cfg = base; cfg.seed = seed;
    cfg.noise_seed = base.noise_seed + (seed - base.seed) + 1000 * b;
    cfg.duration = S.settle + S.block_len;
    ic = struct('dc_levels', A, 'hold_dur', cfg.duration, 'ramp_dur', 0, ...
        'noise_intensity', 0, 'intrinsic_drive', [], 'generator', @dc_staircase_stimulus);
    [m, groups] = ied_build_network(cfg, ic, 'lya_method', 'none', ...
        'T_range', [-cfg.warmup cfg.duration], 'T_plot', [0 cfg.duration], ...
        'store_full_state', true);
    m.run();
    sim_s = toc(t0);
    params = m.get_params();
    ws = ospec_weight_stats(params);
    P = pops(params, groups);

    % Detector, as run_block.
    p = m.plot_data; t = p.t(:)';
    x = [p.x.E; p.x.I]; r = [p.r.E; p.r.I];
    sel = t >= S.settle & t < cfg.duration;
    det = ied_detect_events(t(sel), x(:, sel), r(:, sel), groups, struct());
    G = max(groups); module_x = zeros(G, nnz(sel));
    for g = 1:G, module_x(g, :) = mean(x(groups == g, sel), 1); end
    ev = det.events;
    B = struct('block', b, 'amplitude', A, 'noise_seed', cfg.noise_seed, 'sim_s', sim_s);
    B.signals = struct('t', single(t(sel)), 'module_x', single(module_x), ...
        'mean_x_E', single(mean(x(m.type_indices{1}, sel), 1)), ...
        'mean_x_I', single(mean(x(m.type_indices{2}, sel), 1)), ...
        'mean_r_E', single(mean(r(m.type_indices{1}, sel), 1)));
    B.events = struct('time', [ev.time], 'width', [ev.width], ...
        'recruitment', [ev.recruitment], 'groups', {{ev.groups}}, 'peak_x', [ev.peak_x]);
    B.det_metrics = det.metrics;

    % Sample grid on the full-rate state.
    tf = m.t_out(:);
    tc = (S.settle:S.coarse_dt:cfg.duration - 1e-9)';
    td = [];
    for e = 1:numel(ev)
        td = [td; (ev(e).time + S.dense_win(1):S.dense_dt:ev(e).time + S.dense_win(2))']; %#ok<AGROW>
    end
    ts = unique(round([tc; td] / S.dense_dt) * S.dense_dt);
    ts = ts(ts >= S.settle & ts < cfg.duration);
    idx = round((ts - tf(1)) * S.fs) + 1;
    ts = tf(idx);
    is_coarse = ismember(round(ts / S.dense_dt), round(tc / S.dense_dt));
    if b == 1
        kc = find(ts >= S.settle + 5, 1);
        if isempty(kc), kc = ceil(numel(ts) / 2); end
        R.check = m_check(m.S_out(idx(kc), :)', ts(kc), m, params, ws);
        R.groups = groups; R.tau_d = params.tau_d; R.P = P;
        R.struct = struct('ev_W', eig(ws.Wf), 'R_struct', ws.R_struct, ...
            'K10_eig', eig(P.mu10 .* P.npop(:)'));
        R.ospec_full_probe = ospec_full_probe(m.S_out(idx(1), :)', params, ws);
    end
    t1 = tic;
    B.spec = spectral_series(m.S_out(idx, :), ts, params, ws, P, S.n_top);
    B.spec.is_coarse = is_coarse;
    B.spec_s = toc(t1);
    logline(S, sprintf('seed %d block %d: sim %.0f s, %d events, %d spectra in %.0f s', ...
        seed, b, sim_s, numel(ev), numel(ts), B.spec_s));
    t1 = tic;
    tj = (S.settle:S.full_dt:cfg.duration - 1e-9)';
    jdx = round((tj - tf(1)) * S.fs) + 1;
    B.full = full_series(m.S_out(jdx, :), tf(jdx), params, ws, P);
    B.full_s = toc(t1);
    logline(S, sprintf('seed %d block %d: %d full Jacobians in %.0f s', seed, b, numel(tj), B.full_s));
    R.block(b) = B;
    clear m
end
end

%% ------------------------------------------------------------------------
function P = pops(params, groups)
% The ten populations: subnet g E is pop 2g-1, subnet g I is pop 2g.
n = params.n; G = max(groups);
type_of = zeros(n, 1);
for q = 1:params.n_cellTypes, type_of(params.type_indices{q}) = q; end
pop = 2 * (groups(:) - 1) + type_of;
U10 = zeros(n, 2 * G);
for k = 1:2 * G, U10(pop == k, k) = 1 / sqrt(nnz(pop == k)); end
Wf = full(params.W);
mu10 = zeros(2 * G); s210 = zeros(2 * G);
for a = 1:2 * G
    for c = 1:2 * G
        blk = Wf(pop == a, pop == c);
        mu10(a, c) = mean(blk(:)); s210(a, c) = var(blk(:), 1);
    end
end
P = struct('pop', pop, 'type_of', type_of, 'groups', groups(:), 'U10', U10, ...
    'npop', accumarray(pop, 1)', 'mu10', mu10, 's210', s210, 'G', G);
end

function Q = spectral_series(Sx, ts, params, ws, P, n_top)
ns = numel(ts); n = params.n; G = P.G;
Q = struct('t', ts(:)', 'lam', complex(nan(ns, n_top, 'single')), ...
    'coh10', nan(ns, n_top, 'single'), 'coh2', nan(ns, n_top, 'single'), ...
    'share', nan(ns, n_top, G, 'single'), 'R_eff', nan(ns, 1), 'R10', nan(ns, 1), ...
    'K10_eig', complex(nan(ns, 2 * G)), 'dbar', nan(ns, 2, 2), 'rate', nan(ns, 2), ...
    'rate_sub', nan(ns, G), 'ev', single(complex(nan(ns, n))));
for q = 1:ns
    [Mf, d, ~, pred, pop] = ospec_effective_matrix(Sx(q, :)', params, ws);
    [V, D] = eig(Mf, 'vector');
    Q.ev(q, :) = single(D).';
    up = find(imag(D) >= 0);
    [~, o] = sort(real(D(up)), 'descend');
    pick = up(o(1:min(n_top, numel(o))));
    Vt = V(:, pick); Vt = Vt ./ vecnorm(Vt);
    Q.lam(q, 1:numel(pick)) = single(D(pick)).';
    Q.coh10(q, 1:numel(pick)) = single(sum(abs(P.U10' * Vt).^2, 1));
    Q.coh2(q, 1:numel(pick)) = single(sum(abs(ws.U' * Vt).^2, 1));
    e = abs(Vt).^2;
    for g = 1:G
        Q.share(q, 1:numel(pick), g) = single(sum(e(P.groups == g, :), 1));
    end
    % 10-population predictions from the gains (d is C x n: gain on the
    % route to the postsynaptic type).
    K10 = zeros(2 * G); V10 = zeros(2 * G);
    for a = 1:2 * G
        ta = mod(a - 1, 2) + 1;
        for c = 1:2 * G
            gj = d(ta, P.pop == c);
            K10(a, c) = P.mu10(a, c) * sum(gj);
            V10(a, c) = P.s210(a, c) * sum(gj.^2);
        end
    end
    k10 = eig(K10);
    [~, o] = sort(real(k10), 'descend');
    Q.K10_eig(q, :) = k10(o).';
    Q.R10(q) = sqrt(max(0, max(real(eig(V10)))));
    Q.R_eff(q) = pred.R_eff;
    Q.dbar(q, :, :) = pop.dbar;
    Q.rate(q, :) = pop.rate;
    Q.rate_sub(q, :) = accumarray(P.groups, pop.r, [G 1], @mean)';
end
end

function F = full_series(Sx, ts, params, ws, P)
ns = numel(ts); nk = 12;
F = struct('t', ts(:)', 'abscissa', nan(ns, 1), 'kin_edge', nan(ns, 1), ...
    'net', complex(nan(ns, 1)), 'net_xfrac', nan(ns, 1), 'net_coh10', nan(ns, 1), ...
    'net_share', nan(ns, P.G), 'netc', complex(nan(ns, 1)), 'netc_coh10', nan(ns, 1), ...
    'top_lam', complex(nan(ns, nk)), 'top_xfrac', nan(ns, nk), 'top_coh10', nan(ns, nk), ...
    'n_net_unstable', nan(ns, 1));
for q = 1:ns
    J = SRNNCellTypePairs.compute_Jacobian_fast(Sx(q, :)', params);
    [V, D] = eig(full(J), 'vector');
    V = V ./ vecnorm(V);
    e = abs(V).^2;
    xfrac = sum(e(ws.is_x, :), 1)';
    vx = V(ws.is_x, :);
    coh10 = (sum(abs(P.U10' * vx).^2, 1) ./ max(sum(abs(vx).^2, 1), realmin))';
    [~, o] = sort(real(D), 'descend');
    F.abscissa(q) = real(D(o(1)));
    F.kin_edge(q) = ospec_kinetic_edge(J, ws);
    F.top_lam(q, :) = D(o(1:nk)).'; F.top_xfrac(q, :) = xfrac(o(1:nk))'; F.top_coh10(q, :) = coh10(o(1:nk))';
    net = o(xfrac(o) >= 0.1);
    F.n_net_unstable(q) = nnz(real(D(net)) > 0);
    k = net(1);
    F.net(q) = complex(real(D(k)), abs(imag(D(k))));
    F.net_xfrac(q) = xfrac(k); F.net_coh10(q) = coh10(k);
    F.net_share(q, :) = accumarray(P.groups, e(ws.is_x, k), [P.G 1])' / sum(e(ws.is_x, k));
    kc = net(find(coh10(net) >= 0.2, 1));
    if ~isempty(kc)
        F.netc(q) = complex(real(D(kc)), abs(imag(D(kc)))); F.netc_coh10(q) = coh10(kc);
    end
end
end

function c = ospec_full_probe(S, params, ws)
% Does ospec_full_spectrum accept this model class? (It takes eigs above N = 1000.)
try
    f = ospec_full_spectrum(S, params, ws, 12);
    c = struct('ok', true, 'n_converged', f.n_converged, 'net', f.net, ...
        'kin_edge', f.kin_edge, 'abscissa', f.abscissa, 'seconds', f.seconds, 'msg', '');
catch err
    c = struct('ok', false, 'n_converged', NaN, 'net', NaN, 'kin_edge', NaN, ...
        'abscissa', NaN, 'seconds', NaN, 'msg', err.message);
end
end

function c = m_check(S, tt, m, params, ws)
% M against W diag(F phi') formed from the route outputs, and against
% I + tau_d times a central finite difference of dynamics_fast (as so_run_network).
[Mf, ~, J, ~, pop] = ospec_effective_matrix(S, params, ws);
p = params; p.u_interpolant = m.u_interpolant;
L = p.state_layout; n = p.n;
route = SRNNCellTypePairs.active_routes(p);
Mi = zeros(n);
for pre = 1:p.n_cellTypes
    for post = 1:p.n_cellTypes
        Rt = route{pre, post};
        pidx = p.type_indices{pre}; qidx = p.type_indices{post};
        [~, Fo] = Rt.model.rhs(reshape(S(Rt.idx), p.n_per_type(pre), Rt.P.ns), pop.r(pidx), Rt.P);
        Mi(qidx, pidx) = full(p.W(qidx, pidx)) .* (Fo(:) .* pop.rp(pidx))';
    end
end
h = 1e-7; Jfd = zeros(n);
for j = 1:n
    Sp = S; Sp(L.x(j)) = Sp(L.x(j)) + h;
    Sm = S; Sm(L.x(j)) = Sm(L.x(j)) - h;
    fp = SRNNCellTypePairs.dynamics_fast(tt, Sp, p);
    fm = SRNNCellTypePairs.dynamics_fast(tt, Sm, p);
    Jfd(:, j) = (fp(L.x) - fm(L.x)) / (2 * h);
end
c = struct('t', tt, 'max_abs_M', max(abs(Mf), [], 'all'), ...
    'M_vs_I_plus_tau_Jxx', max(abs(Mf - (eye(n) + p.tau_d * full(J(L.x, L.x)))), [], 'all'), ...
    'M_vs_W_diag_F_rp', max(abs(Mf - Mi), [], 'all'), ...
    'M_vs_I_plus_tau_Jxx_FD', max(abs(Mf - (eye(n) + p.tau_d * Jfd)), [], 'all'));
end

function logline(S, msg)
fprintf('%s %s\n', char(datetime('now', 'Format', 'HH:mm:ss')), msg);
if isfield(S, 'log') && ~isempty(S.log)
    fid = fopen(S.log, 'a'); fprintf(fid, '%s %s\n', char(datetime('now', 'Format', 'HH:mm:ss')), msg); fclose(fid);
end
end
