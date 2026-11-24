% S2B_05_FIGS  Three figures from the saved results of s2b_01 .. s2b_04.
[ddir, fdir] = s2b_out_dirs();
W = load(fullfile(ddir, 's2b_01_window.mat'));
L = load(fullfile(ddir, 's2b_02_long.mat'));
S = load(fullfile(ddir, 's2b_03_seeds.mat'));
U = load(fullfile(ddir, 's2b_04_warmup.mat'));
cols = lines(6);  ttl = struct('sfa1_std1', 'Single-timescale', 'sfa3_stdc2', 'Multiple-timescale');
names = {'sfa1_std1', 'sfa3_stdc2'};

%% Figure 1: the stage's window; the trace depends on the vector, not on the method
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1300 900]);
tl = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for c = 1:2
    o = W.(names{c});
    ax = nexttile(tl, c);  hold(ax, 'on');
    plot(ax, o.saved.benettin.t_lya, o.saved.benettin.local_lya, 'Color', [0 0.45 0.74], 'LineWidth', 1.6);
    plot(ax, o.saved.qr.t_lya, o.saved.qr.local_LE_spectrum_t(:, 1), 'Color', [0.25 0.25 0.25], 'LineWidth', 0.8);
    ylim(ax, [-4 4]);  xlim(ax, [-8 4]);
    title(ax, sprintf('%s, 30 neurons, seeds [1 2]: panel B as saved by the stage', ttl.(names{c})));
    if c == 1; ylabel(ax, 'local rate (s^{-1})'); legend(ax, {'Benettin (random vector, 20 ms)', 'QR \lambda_1 (first vector = e_1, 100 ms)'}, 'Location', 'southwest'); end
    ax = nexttile(tl, 2 + c);  hold(ax, 'on');
    for j = 1:6; plot(ax, o.R100.t, o.R100.loc_v(:, j), 'Color', cols(j, :), 'LineWidth', 0.9); end
    plot(ax, o.R100.t, o.R100.loc_qr(:, 1), 'k', 'LineWidth', 1.6);
    ylim(ax, [-4 4]);  xlim(ax, [-8 4]);
    title(ax, 'six random vectors, each evolved alone (colours), and the QR first vector (black)');
    if c == 1; ylabel(ax, 'local rate (s^{-1})'); end
    ax = nexttile(tl, 4 + c);  hold(ax, 'on');
    plot(ax, o.R100.t, o.R100.loc_qr(:, 1), 'k', 'LineWidth', 2.4);
    plot(ax, o.fd_t, o.fd_e1, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.0);
    plot(ax, o.ref_t, o.ref_e1, ':', 'Color', [0.2 0.7 0.2], 'LineWidth', 1.6);
    plot(ax, o.Rseg{1}.t, o.Rseg{1}.loc_v(:, 1), 'Color', [0.6 0.75 0.95], 'LineWidth', 2.2);
    plot(ax, o.fd_t, o.fd_rand, 'Color', [0 0.2 0.6], 'LineWidth', 0.8);
    plot(ax, o.Rseg{5}.t, o.Rseg{5}.loc_v(:, 1), '--', 'Color', [0.5 0 0.5], 'LineWidth', 1.0);
    ylim(ax, [-4 4]);  xlim(ax, [-8 4]);  xlabel(ax, 'time (s)');
    title(ax, 'same vector, different method / segment / precision');
    if c == 1; ylabel(ax, 'local rate (s^{-1})'); legend(ax, {'QR first vector, 100 ms (RK4 variational)', ...
        'Benettin of the class (d_0 = 10^{-3}, ode45, 20 ms) started from e_1', 'ode45 10^{-10} variational reference, e_1', ...
        'random vector 1, linear, 5 ms', 'random vector 1, Benettin of the class, 20 ms', 'random vector 1, linear, 100 ms'}, ...
        'Location', 'southwest', 'FontSize', 7); end
end
exportgraphics(fig, fullfile(fdir, 'S2B_1_window_vectors.png'), 'Resolution', 150);  close(fig);

%% Figure 2: the long run
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1300 950]);
tl = tiledlayout(fig, 4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for c = 1:2
    o = L.(names{c});
    ax = nexttile(tl, c);  semilogy(ax, o.t, o.speed, 'k');  xline(ax, [-8 4], ':');
    title(ax, sprintf('%s, 30 neurons, seeds [1 2], no noise: speed |dS/dt|', ttl.(names{c})));  xlim(ax, [-8 200]);
    ax = nexttile(tl, 2 + c);  hold(ax, 'on');
    for j = 1:6; plot(ax, o.t, o.loc_v(:, j), 'Color', cols(j, :)); end
    plot(ax, o.t, o.loc_qr(:, 1), 'k');  xlim(ax, [-8 60]);  ylim(ax, [-4 4]);
    title(ax, 'local rate: six random vectors (colours), QR first vector (black); first 60 s');
    if c == 1; ylabel(ax, 's^{-1}'); end
    ax = nexttile(tl, 4 + c);  plot(ax, o.t, o.proj_q1);  ylim(ax, [0 1]);  xlim(ax, [-8 200]);
    title(ax, '|cos| between each random vector and the QR first vector');  xlabel(ax, 'time (s)');
    ax = nexttile(tl, 6 + c);  sp = sort(mean(o.loc_qr(o.t >= 100, :)), 'descend');
    stem(ax, 1:20, sp(1:20), 'filled', 'Color', 'k');  hold(ax, 'on');
    if o.speed(end) < 1e-6; plot(ax, 1:20, real(o.eig_end(1:20)), 'ro'); legend(ax, {'QR exponents over [100, 200] s', 'Re eig J at the fixed point'}, 'Location', 'southwest'); end
    title(ax, sprintf('20 leading exponents; \\lambda_1 - \\lambda_2 = %.3f s^{-1}, %d within 0.05 s^{-1} of \\lambda_1', sp(1) - sp(2), nnz(sp >= sp(1) - 0.05)));
    xlabel(ax, 'index');  ylabel(ax, 's^{-1}');
end
exportgraphics(fig, fullfile(fdir, 'S2B_2_long_run.png'), 'Resolution', 150);  close(fig);

%% Figure 3: other networks, and the warm-up needed
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1300 520]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl);  hold(ax, 'on');
r = S.rows;  mk = struct('sfa1_std1', [0.8 0.6 0], 'sfa3_stdc2', [0 0.45 0.74]);
for k = 1:numel(r)
    sus = r(k).speed60 > 1e-2;                       % activity sustained at 60 s
    h = plot(ax, median(r(k).cos_pair), median(r(k).rms_v_q1), 'o', 'Color', mk.(r(k).cname), 'MarkerSize', 6 + 6 * (r(k).n == 60), 'LineWidth', 1.3);
    if sus; set(h, 'MarkerFaceColor', mk.(r(k).cname)); end
    text(ax, median(r(k).cos_pair) + 0.015, median(r(k).rms_v_q1), sprintf('[%d %d]', r(k).seeds), 'FontSize', 7);
end
xlabel(ax, 'median |cos| between random vectors at t = 4 s (1 = one direction found)');
ylabel(ax, 'rms (random vector - QR first vector) local rate, t in [-4, 4] (s^{-1})');
title(ax, {'8 networks x 2 conditions, the stage''s window', 'yellow single, blue multiple; filled = activity sustained at 60 s; large = 60 neurons'});
ax = nexttile(tl);  hold(ax, 'on');
for k = 1:numel(U.out)
    o = U.out(k);
    semilogy(ax, o.w0 + 4, max(median(o.rms_w, 2), 1e-4), '-o', 'LineWidth', 1.3, 'DisplayName', ...
        sprintf('%s n = %d [%d %d]: \\lambda_1 - \\lambda_2 = %.3f', ttl.(o.cname), o.n, o.seeds, max(o.lam) - max(o.lam(o.lam < max(o.lam)))));
end
set(ax, 'YScale', 'log');  legend(ax, 'Location', 'southwest', 'FontSize', 7);
xlabel(ax, 'centre of the 8-s window (s); the vectors start at -8 s');  ylabel(ax, 'median rms (random - QR first vector) local rate (s^{-1})');
title(ax, 'sustained activity: how long until the two estimators agree');
exportgraphics(fig, fullfile(fdir, 'S2B_3_seeds_warmup.png'), 'Resolution', 150);  close(fig);
close all force
fprintf('figures in %s\n', fdir);
