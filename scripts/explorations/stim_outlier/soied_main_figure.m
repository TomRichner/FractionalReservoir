function info = soied_main_figure(R, Mz, out_dir)
% SOIED_MAIN_FIGURE Discharges, the coherent frozen-gain mode and the full Jacobian's network mode in Fig 8's IED network, baseline beside stimulation.
%
%   info = soied_main_figure(R, Mz, out_dir)
%
% R from soied_run, Mz = soied_measures(R). Writes soied_main.png (300 dpi)
% and .fig to out_dir; info holds the chosen windows and states.
%
% Columns: a 20-s stretch of the baseline block and of the stimulation
% block (the window holding the most discharges). Rows:
%  (A) the subnet-mean dendritic potentials (E mean, I mean, five subnets,
%      stacked as Fig 8 D) with the detector's marks
%  (B) Re of the rightmost coherent mode of M (coh10 >= Mz.thr; black),
%      Re lambda_1(M) (any mode; gray), predicted bulk radius R_eff
%      (dashed), boundary Re = 1 (red)
%  (C) the full Jacobian's rightmost network mode (1/s), boundary 0
%  (D) zoom on one discharge: its subnet's mean potential and rc, dense
%      grid; event-triggered medians of rc (frozen, all discharges) and of
%      the full network mode (6 discharges, 6 control times)
%  (E) eigenvalues of M as rates (lambda - 1)/tau_d in the complex plane
%      (complex_plane_axes; the imaginary axis is the boundary) at a quiet
%      baseline state, just before the zoomed discharge, and during
%      stimulation; the coherent mode ringed; one scale
%
% See also: soied_measures, soied_run, complex_plane_axes

if ~isfolder(out_dir), mkdir(out_dir); end
fs = 12; win = 20; tau = R.tau_d; thr = Mz.thr;
f = figure('Visible', 'off', 'Units', 'pixels', 'Position', [50 50 1500 1500], 'Color', 'w');
col = {[0 0 0], [0.2 0.4 0.8]};
red = [0.8 0.1 0.1];
info = struct();
L = 0.07; Wc = 0.41; gap = 0.06;
for b = 1:2
    B = R.block(b); Q = B.spec; te = B.events.time;
    t = double(B.signals.t);
    % window with the most discharges
    starts = t(1):1:t(end) - win;
    cnt = arrayfun(@(a) sum(te >= a & te < a + win), starts);
    [~, k] = max(cnt); t0 = starts(k); w = [t0 t0 + win];
    info.window(b, :) = w;
    x0 = L + (b - 1) * (Wc + gap);
    % (A) traces
    ax = axes(f, 'Position', [x0 0.80 Wc 0.17]); hold(ax, 'on');
    rows = [double(B.signals.mean_x_E); double(B.signals.mean_x_I); double(B.signals.module_x)];
    sel = t >= w(1) & t <= w(2); sp = 1.0; nR = size(rows, 1);
    for r = 1:nR
        tr = rows(r, sel);
        c = [0.15 0.15 0.15]; if r == 2, c = [0.25 0.4 0.75]; end
        if r > 2, c = [0.55 0.1 0.1] + (r - 3) * [0.08 0.08 0.05]; end
        plot(ax, t(sel), tr - median(tr) + (nR - r) * sp, '-', 'Color', c, 'LineWidth', 0.8);
    end
    in = te(te >= w(1) & te <= w(2));
    plot(ax, in, ((nR - 1) * sp + 0.6) * ones(size(in)), 'v', 'Color', [0.15 0.15 0.15], ...
        'MarkerFaceColor', [0.15 0.15 0.15], 'MarkerSize', 6);
    xlim(ax, w); ylim(ax, [-0.6 (nR - 1) * sp + 0.9]);
    set(ax, 'YTick', (0:nR - 1) * sp, 'YTickLabel', fliplr({'E', 'I', 'S1', 'S2', 'S3', 'S4', 'S5'}), ...
        'XTick', [], 'FontSize', fs, 'Box', 'off', 'XColor', 'none');
    if b == 1
        title(ax, 'Baseline (input 0)', 'FontSize', fs + 2, 'FontWeight', 'normal');
        ylabel(ax, 'mean x', 'FontSize', fs);
    else
        title(ax, sprintf('Stimulation (uniform input %.1f)', B.amplitude), 'FontSize', fs + 2, 'FontWeight', 'normal');
    end
    % (B) frozen-gain modes
    ax2 = axes(f, 'Position', [x0 0.58 Wc 0.19]); hold(ax2, 'on');
    ts = Q.t(:); c = Q.is_coarse(:); s2 = c & ts >= w(1) & ts <= w(2);
    rc = coherent_re(Q, thr);
    for e = in(:)', plot(ax2, [e e], [0.4 2], '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.8); end
    plot(ax2, ts(s2), real(Q.lam(s2, 1)), '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.7);
    plot(ax2, ts(s2), rc(s2), '-', 'Color', [0 0 0], 'LineWidth', 1.0);
    plot(ax2, ts(s2), Q.R_eff(s2), '--', 'Color', [0.2 0.5 0.2], 'LineWidth', 1.2);
    plot(ax2, w, [1 1], '-', 'Color', red, 'LineWidth', 1.2);
    xlim(ax2, w); ylim(ax2, [0.4 2]);
    set(ax2, 'FontSize', fs, 'Box', 'off', 'XTick', [], 'XColor', 'none');
    if b == 1
        ylabel(ax2, 'Re \lambda(M)', 'FontSize', fs);
    end
    % (C) full Jacobian network mode
    ax3 = axes(f, 'Position', [x0 0.44 Wc 0.11]); hold(ax3, 'on');
    F = B.full; sf = F.t >= w(1) & F.t <= w(2);
    for e = in(:)', plot(ax3, [e e], [-3 5], '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.8); end
    plot(ax3, w, [0 0], '-', 'Color', red, 'LineWidth', 1.2);
    plot(ax3, F.t(sf), real(F.net(sf)), 'o-', 'Color', [0.3 0.1 0.5], 'MarkerFaceColor', [0.3 0.1 0.5], ...
        'MarkerSize', 4, 'LineWidth', 0.8);
    xlim(ax3, w); ylim(ax3, [-3 5]);
    set(ax3, 'FontSize', fs, 'Box', 'off', 'TickDir', 'out');
    xlabel(ax3, 'time in block (s)', 'FontSize', fs);
    if b == 1, ylabel(ax3, {'full J network', 'mode Re (s^{-1})'}, 'FontSize', fs); end
end
% legend for (B), drawn as text-free line keys in the right panel's corner
axk = axes(f, 'Position', [0.89 0.70 0.10 0.07], 'Visible', 'off'); hold(axk, 'on');
hk(1) = plot(axk, nan, nan, '-', 'Color', [0 0 0], 'LineWidth', 1);
hk(2) = plot(axk, nan, nan, '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.7);
hk(3) = plot(axk, nan, nan, '--', 'Color', [0.2 0.5 0.2], 'LineWidth', 1.2);
hk(4) = plot(axk, nan, nan, '-', 'Color', red, 'LineWidth', 1.2);
legend(axk, hk, {sprintf('coherent mode (coh_{10} \\geq %.1f)', thr), '\lambda_1, any mode', 'bulk radius R_{eff}', 'boundary'}, ...
    'FontSize', fs - 2, 'Box', 'off', 'Location', 'north');

% (D) zoom on one discharge of the baseline window
B = R.block(1); Q = B.spec; te = B.events.time; ts = Q.t(:); rc = coherent_re(Q, thr);
w = info.window(1, :);
cand = find(te >= w(1) + 1 & te <= w(2) - 1);
[~, kk] = max(B.events.peak_x(cand)); ez = cand(kk); tz = te(ez);
g = B.events.groups{ez}(1);
info.zoom_event = tz; info.zoom_group = B.events.groups{ez};
zw = [tz - 0.8 tz + 0.4];
t = double(B.signals.t); sel = t >= zw(1) & t <= zw(2);
ax = axes(f, 'Position', [L 0.24 0.25 0.14]); hold(ax, 'on');
yyaxis(ax, 'left');
plot(ax, t(sel) - tz, double(B.signals.module_x(g, sel)), '-', 'Color', [0.55 0.1 0.1], 'LineWidth', 1.2);
ylabel(ax, sprintf('mean x, S%d', g), 'FontSize', fs); ax.YColor = [0.55 0.1 0.1];
yyaxis(ax, 'right'); hold(ax, 'on');
sz = ts >= zw(1) & ts <= zw(2);
plot(ax, ts(sz) - tz, rc(sz), '.-', 'Color', [0 0 0], 'MarkerSize', 6, 'LineWidth', 0.8);
plot(ax, zw - tz, [1 1], '-', 'Color', red, 'LineWidth', 1);
ylabel(ax, 'Re, coherent mode of M', 'FontSize', fs); ax.YColor = [0 0 0];
xlim(ax, zw - tz); set(ax, 'FontSize', fs, 'Box', 'off');
xlabel(ax, 'time from discharge peak (s)', 'FontSize', fs);
% pre-discharge state: the max of rc in [-0.3, 0]
sp = find(ts >= tz - 0.3 & ts <= tz);
[~, i] = max(rc(sp)); i_pre = sp(i); info.t_pre = ts(i_pre);
plot(ax, ts(i_pre) - tz, rc(i_pre), 'o', 'Color', [0.9 0.5 0], 'MarkerSize', 8, 'LineWidth', 1.5);

% event-triggered medians
ax = axes(f, 'Position', [0.40 0.24 0.22 0.14]); hold(ax, 'on');
wb = Mz.block(1);
plot(ax, wb.trig_lags, wb.trig_rc, '-', 'Color', [0 0 0], 'LineWidth', 1.4);
plot(ax, [-0.6 0.2], [1 1], '-', 'Color', red, 'LineWidth', 1);
ylabel(ax, 'Re, coherent mode of M', 'FontSize', fs); ylim(ax, [0.8 1.8]);
xlabel(ax, 'time from discharge peak (s)', 'FontSize', fs);
set(ax, 'FontSize', fs, 'Box', 'off');
ax = axes(f, 'Position', [0.71 0.24 0.22 0.14]); hold(ax, 'on');
plot(ax, Mz.ef.lags, Mz.ef.ev_net_med, 'o-', 'Color', [0.3 0.1 0.5], 'MarkerFaceColor', [0.3 0.1 0.5], 'MarkerSize', 4, 'LineWidth', 1.2);
plot(ax, Mz.ef.lags, Mz.ef.ctrl_net_med, 's--', 'Color', [0.6 0.6 0.6], 'MarkerSize', 4, 'LineWidth', 1.0);
plot(ax, [-0.6 0.2], [0 0], '-', 'Color', red, 'LineWidth', 1);
ylabel(ax, 'full J network mode Re (s^{-1})', 'FontSize', fs);
xlabel(ax, 'time from discharge peak (s)', 'FontSize', fs);
legend(ax, {sprintf('%d discharges', Mz.ef.n), sprintf('%d control times', Mz.ef.n)}, 'Box', 'off', 'FontSize', fs - 2, 'Location', 'northwest');
set(ax, 'FontSize', fs, 'Box', 'off');

% (E) complex planes
Qb = R.block(1).spec; Qs = R.block(2).spec;
cb = Qb.is_coarse(:); tsb = Qb.t(:); rcb = coherent_re(Qb, thr);
quiet = find(cb & arrayfun(@(x) all(abs(te - x) >= 1), tsb));
[~, i] = min(abs(rcb(quiet) - median(rcb(quiet), 'omitnan'))); i_q = quiet(i);
rcs = coherent_re(Qs, thr); cs = find(Qs.is_coarse(:));
[~, i] = min(abs(rcs(cs) - median(rcs(cs), 'omitnan'))); i_s = cs(i);
info.t_quiet = tsb(i_q); info.t_stim = Qs.t(i_s);
states = {Qb, i_q; Qb, i_pre; Qs, i_s};
opts = {'re_lim', [-26 9], 'im_lim', [-17 17], 'bulk_center', -10, 'bulk_radius', 12, 'font_size', fs};
for p = 1:3
    Qp = states{p, 1}; ii = states{p, 2};
    ax = axes(f, 'Position', [L + (p - 1) * 0.31 0.005 0.27 0.19]); hold(ax, 'on');
    cp = complex_plane_axes(ax, opts{:});
    z = (double(Qp.ev(ii, :)) - 1) / tau;
    plot(ax, cp.x(real(z)), imag(z), '.', 'Color', [0.3 0.3 0.3], 'MarkerSize', 5);
    re = real(double(Qp.lam(ii, :))); re(Qp.coh10(ii, :) < thr) = -Inf; [~, k] = max(re);
    zc = (double(Qp.lam(ii, k)) - 1) / tau;
    plot(ax, cp.x(real([zc conj(zc)])), imag([zc conj(zc)]), 'o', 'Color', [0.9 0.5 0], 'MarkerSize', 9, 'LineWidth', 1.5);
    th = linspace(0, 2 * pi, 200); rr = Qp.R_eff(ii) / tau;
    plot(ax, cp.x(-1 / tau + rr * cos(th)), rr * sin(th), '--', 'Color', [0.2 0.5 0.2], 'LineWidth', 1);
end
% panel letters
lt = {'(A)', 0.80 + 0.17; '(B)', 0.58 + 0.19; '(C)', 0.44 + 0.11; '(D)', 0.24 + 0.14; '(E)', 0.02 + 0.17};
for k = 1:size(lt, 1)
    annotation(f, 'textbox', [0.005 lt{k, 2} - 0.01 0.04 0.03], 'String', lt{k, 1}, ...
        'EdgeColor', 'none', 'FontSize', fs + 2);
end
info.png = fullfile(out_dir, 'soied_main.png');
exportgraphics(f, info.png, 'Resolution', 300);
savefig(f, fullfile(out_dir, 'soied_main.fig'));
close(f);
end

function rc = coherent_re(Q, thr)
re = real(double(Q.lam)); re(Q.coh10 < thr | isnan(re)) = -Inf;
rc = max(re, [], 2); rc(~isfinite(rc)) = NaN;
end
