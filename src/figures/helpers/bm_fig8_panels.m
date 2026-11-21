function P = bm_fig8_panels(fig, rect, R, opts)
%BM_FIG8_PANELS Draw Figure 8 A to C (the bursting network's mode) into a rectangle of a figure.
%
%   P = bm_fig8_panels(fig, rect, R)
%   P = bm_fig8_panels(fig, rect, R, 'font_size', 9, ...)
%
% fig   the figure (Units are set to inches)
% rect  [x y w h] in inches, the area the three panels fill
% R     one network of the burst_mode stage (bm_load_runs(...).R{1}: the
%       stage's first seed)
%
% A  the whole run: the input (mean over neurons); mean E and I rates, the
%    frozen-adaptation arm dashed after the step; the real part of the
%    rightmost mode as a rate (1/s), gains frozen (gray) and the rightmost
%    network mode of the full Jacobian (the bursting condition's color),
%    the frozen arm's frozen-gain mode dashed. A light band marks the burst
%    of B; the four sampled states of C are marked on A and B.
% B  that burst, from 1 s before its onset to 0.5 s after: rates, and the
%    two modes on the full Jacobian's dense grid (every dt_sample there).
% C  the eigenvalues of M as rates z = (lambda - 1)/tau_d in the complex
%    plane (complex_plane_axes: the imaginary axis is the boundary, the
%    bulk is centered at -1/tau_d), one scale across the five panels:
%    the structure (W alone: every gain 1), a quiet moment, just before the
%    burst, during it, and under stimulation. The bulk circle has the
%    radius predicted from that state's gains (R_eff/tau_d; R_struct for
%    W); the rightmost eigenvalue is drawn larger.
%
% THE FIVE STATES (rules, not choices made by eye): the burst of B is the
% one before the step whose onset is closest to the middle of the baseline
% window; quiet = 0.9 s before its onset; before = the maximum of the
% frozen-gain mode within tol (0.6 s) before onset; during = the burst's
% peak (the mean E rate's maximum); stimulation = 5 s before the end of
% the run.
%
% Options: font_size (9), annot_size (8), line_scale (1, multiplies every
% line width), plane_gap (0.12 in), colors (struct, see below).
%
% P: handles ax_A (3 axes, top to bottom), ax_B (2), ax_C (5), and info:
% the chosen times and the window of the complex plane (for the table).
%
% See also: fig_burst_mode, ied_stim_fig_combined, complex_plane_axes

arguments
    fig
    rect (1,4) double
    R (1,1) struct
    opts.font_size (1,1) double = 9
    opts.annot_size (1,1) double = 8
    opts.line_scale (1,1) double = 1
    opts.plane_gap (1,1) double = 0.12
end
st = manuscript_style();
fs = opts.font_size; lw = opts.line_scale;
tau = R.tau_d; H = R.t_step;
I = R.intact; F = R.full; Z = R.frozen;
ec = SRNNCellTypePairs.type_colors(2);
cond_col = st.condition_color(R.condition);
mode_col = 0.6 * cond_col;                        % the condition's hue, darkened to read as a line
fg_col = [0.62 0.62 0.62];                        % gains frozen
purple = [0.494 0.247 0.596]; green = [0 0.62 0.45];   % no stimulation / stimulation (Figure 8's)
state_col = {[0 0 0], 0.45 * purple + 0.55, purple, 0.6 * purple, green};
state_name = {'Structure', 'Quiet', 'Before burst', 'Burst', 'Stimulation'};

%% The burst of B and the five states
on = R.bursts.t_on; off = R.bursts.t_off;
cand = find(on < H & on - 1 > 0);
if isempty(cand)
    error('bm_fig8_panels:NoBurst', 'The network has no burst before the step to zoom on.');
end
[~, i] = min(abs(on(cand) - mean(R.windows.before)));
b = cand(i);
zoom = [on(b) - 1, on(b) + 0.5];
re1 = real(I.z1);
k = find(I.t >= on(b) - 0.6 & I.t <= on(b));
[~, j] = max(re1(k));
t_state = [NaN, on(b) - 0.9, I.t(k(j)), R.bursts.peak_t(b), I.t(end) - 5];
idx = arrayfun(@(tt) nearest_index(I.t, tt), t_state(2:end));

%% Geometry (inches)
set(fig, 'Units', 'inches');
x0 = rect(1); y0 = rect(2); W = rect(3); Ht = rect(4);
gap = opts.plane_gap; left = 0.52;
pw = (W - left - 4 * gap) / 5;

% complex plane window, shared
zs = cell(1, 5); cen = -1 / tau; rad = zeros(1, 5);
zs{1} = (R.structure.ev_W(:) - 1) / tau; rad(1) = R.structure.R_struct / tau;
for q = 2:5
    zs{q} = (double(I.ev(idx(q - 1), :)).' - 1) / tau;
    rad(q) = I.R_eff(idx(q - 1)) / tau;
end
zall = vertcat(zs{:}); rmax = max(rad);
re_lim = [cen - 1.15 * rmax, max(cen + 1.15 * rmax, 0.2 * rmax)];
re_lim(2) = max(re_lim(2), max(real(zall)) + 0.08 * rmax);
im_lim = [-1.3 1.3] * max(rmax, max(abs(imag(zall))));      % room above the bulk for the Im name
step = nice_step(rmax / 2);
% Few ticks on small panels: on the real axis the multiples of step well
% clear of the largest bulk (a quarter radius), on the imaginary axis the
% first multiple outside it.
mult = @(lo, hi) step * (ceil(lo / step):floor(hi / step));
re_ticks = mult(cen + 1.25 * rmax, re_lim(2));
if isempty(re_ticks); re_ticks = mult(cen + 1.02 * rmax, re_lim(2)); end
im_t = mult(1.0001 * rmax, im_lim(2));
if isempty(im_t)                                  % widen the window to the next multiple
    im_t = step * ceil(1.0001 * rmax / step);
    im_lim = [-1.08 1.08] * im_t;
end
im_ticks = [-im_t(1), im_t(1)];
axes_args = {'re_lim', re_lim, 'im_lim', im_lim, 'bulk_center', cen, 'bulk_radius', rmax, ...
    're_ticks', re_ticks, 'im_ticks', im_ticks, 'font_size', opts.annot_size, 'line_width', 0.6 * lw};
g = complex_plane_axes([], axes_args{:});
ph = pw * diff(g.y_lim) / diff(g.x_lim);
c_title = 0.2;
yC = y0 + 0.05;
yAB = yC + ph + c_title + 0.42;                   % room for B's x label and C's titles
hA = [0.32 0.8 0.95]; gA = 0.1;
top = y0 + Ht;
avail = top - yAB;
scale = (avail - 2 * gA) / sum(hA);
hA = hA * scale;
ya = [yAB + hA(3) + hA(2) + 2 * gA, yAB + hA(3) + gA, yAB];
wA = 0.60 * W - left; xB = x0 + 0.71 * W; wB = W - 0.71 * W - 0.05;

%% A
P.ax_A = gobjects(1, 3);
for r = 1:3
    P.ax_A(r) = axes(fig, 'Units', 'inches', 'Position', [x0 + left, ya(r), wA, hA(r)]);
    hold(P.ax_A(r), 'on');
end
ax = P.ax_A(1);
plot(ax, I.t, I.u, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.8 * lw);
ylim(ax, [min(I.u) - 0.05, max(I.u) + 0.1]);
yticks(ax, unique(round([0 R.dc_level], 2)));
ylabel(ax, 'Input');
ax = P.ax_A(2);
band(ax, zoom, [0 1]);
h = gobjects(1, 3);
for q = 2:-1:1
    h(q) = plot(ax, I.t, I.rate(:, q), '-', 'Color', ec(q, :), 'LineWidth', 0.8 * lw);
    after = Z.t >= H;
    plot(ax, Z.t(after), Z.rate(after, q), '--', 'Color', ec(q, :), 'LineWidth', 0.8 * lw);
end
h(3) = plot(ax, NaN, NaN, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.8 * lw);
ylim(ax, [0 1]); yticks(ax, [0 0.5 1]);
ylabel(ax, 'Mean rate');
lg = legend(ax, h, {'E', 'I', 'Adaptation frozen'}, 'Location', 'northwest', 'Box', 'off', ...
    'FontSize', opts.annot_size, 'NumColumns', 3, 'AutoUpdate', 'off'); lg.ItemTokenSize = [12 8];
ax = P.ax_A(3);
zl = mode_limits(re1, real(F.net), real(Z.z1));
band(ax, zoom, zl);
yline(ax, 0, ':', 'Color', [0 0 0], 'LineWidth', 0.6 * lw);
hm = gobjects(1, 2);
hm(1) = plot(ax, I.t, re1, '-', 'Color', fg_col, 'LineWidth', 0.6 * lw);
after = Z.t >= H;
plot(ax, Z.t(after), real(Z.z1(after)), '--', 'Color', fg_col, 'LineWidth', 0.8 * lw);
c = F.is_coarse;
hm(2) = plot(ax, F.t(c), real(F.net(c)), '-', 'Color', mode_col, 'LineWidth', 0.8 * lw);
ylim(ax, zl);
ylabel(ax, 'Re \lambda (s^{-1})');
xlabel(ax, 'Time (s)');
lg = legend(ax, hm, {'Gains frozen', 'Full Jacobian'}, 'Location', 'northeast', 'Box', 'off', ...
    'FontSize', opts.annot_size, 'AutoUpdate', 'off'); lg.ItemTokenSize = [12 8];
for r = 1:3
    xlim(P.ax_A(r), [I.t(1) I.t(end)]);
    if r < 3; P.ax_A(r).XTickLabel = {}; end
end
for q = 2:5
    plot(P.ax_A(3), t_state(q), zl(2), 'v', 'MarkerSize', 4.5, 'MarkerFaceColor', state_col{q}, ...
        'MarkerEdgeColor', 'none', 'Clipping', 'off');
end

%% B
P.ax_B = gobjects(1, 2);
for r = 1:2
    P.ax_B(r) = axes(fig, 'Units', 'inches', 'Position', [xB, ya(r + 1), wB, hA(r + 1)]);
    hold(P.ax_B(r), 'on');
end
ax = P.ax_B(1);
band(ax, [on(b) off(b)], [0 1]);
kz = I.t >= zoom(1) & I.t <= zoom(2);
for q = 2:-1:1
    plot(ax, I.t(kz), I.rate(kz, q), '-', 'Color', ec(q, :), 'LineWidth', 1.0 * lw);
end
ylim(ax, [0 1]); yticks(ax, [0 0.5 1]); ax.XTickLabel = {};
ax = P.ax_B(2);
band(ax, [on(b) off(b)], zl);
yline(ax, 0, ':', 'Color', [0 0 0], 'LineWidth', 0.6 * lw);
plot(ax, I.t(kz), re1(kz), '-', 'Color', fg_col, 'LineWidth', 0.8 * lw);
kf = F.t >= zoom(1) & F.t <= zoom(2);
plot(ax, F.t(kf), real(F.net(kf)), '-', 'Color', mode_col, 'LineWidth', 1.0 * lw);
for q = 2:4
    plot(ax, t_state(q), zl(2), 'v', 'MarkerSize', 4.5, 'MarkerFaceColor', state_col{q}, ...
        'MarkerEdgeColor', 'none', 'Clipping', 'off');
end
ylim(ax, zl); xlabel(ax, 'Time (s)');
for r = 1:2; xlim(P.ax_B(r), zoom); end

%% C
P.ax_C = gobjects(1, 5);
th = linspace(0, 2 * pi, 241);
for q = 1:5
    axq = axes(fig, 'Units', 'inches', 'Position', [x0 + left + (q - 1) * (pw + gap), yC, pw, ph]);
    names = {'re_label', '', 'im_label', ''};
    if q == 5; names{2} = 'Re (s^{-1})'; end
    if q == 1; names{4} = 'Im (s^{-1})'; end
    cp = complex_plane_axes(axq, axes_args{:}, names{:});
    hold(axq, 'on');
    plot(axq, cp.x(cen + rad(q) * cos(th)), rad(q) * sin(th), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6 * lw);
    z = zs{q};
    scatter(axq, cp.x(real(z)), imag(z), 5 * lw, [0.25 0.25 0.25], 'filled', 'MarkerEdgeColor', 'none');
    [~, kr] = max(real(z));
    zr = z(abs(z - z(kr)) < 1e-9 | abs(z - conj(z(kr))) < 1e-9);   % the rightmost eigenvalue and its conjugate
    scatter(axq, cp.x(real(zr)), imag(zr), 22 * lw, state_col{q}, 'filled', 'MarkerEdgeColor', 'none');
    set(axq, 'XLim', cp.x_lim, 'YLim', cp.y_lim);
    text(axq, mean(cp.x_lim), cp.y_lim(2) + 0.04 * diff(cp.y_lim), state_name{q}, 'Color', state_col{q}, ...
        'FontSize', fs, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    P.ax_C(q) = axq;
end

%% Style
for a = [P.ax_A P.ax_B]
    set(a, 'FontSize', fs, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 0.6 * lw, ...
        'TitleFontSizeMultiplier', 1, 'LabelFontSizeMultiplier', 1);
end
P.info = struct('burst', b, 'zoom', zoom, 't_state', t_state, 'state_name', {state_name}, ...
    'z_right', cellfun(@(z) max(real(z)), zs), 'rad', rad, 're_lim', re_lim, 'im_lim', im_lim, ...
    'mode_lim', zl);
end

%% ------------------------------------------------------------------------
function k = nearest_index(t, tt)
[~, k] = min(abs(t - tt));
end

function band(ax, x, y)
patch(ax, x([1 2 2 1]), y([1 1 2 2]), [0.92 0.92 0.92], 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function zl = mode_limits(varargin)
v = cat(1, varargin{:});
v = v(isfinite(v));
lo = floor((min(v) - 1) / 5) * 5; hi = ceil((max(v) + 1) / 5) * 5;
zl = [lo, hi];
end

function s = nice_step(raw)
e = 10^floor(log10(raw)); m = raw / e;
if m < 1.5, s = e; elseif m < 3.5, s = 2 * e; elseif m < 7.5, s = 5 * e; else, s = 10 * e; end
end
