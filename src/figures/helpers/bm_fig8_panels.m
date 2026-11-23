function P = bm_fig8_panels(fig, rect, R, opts)
%BM_FIG8_PANELS Draw Figure 8 G and H (the bursting network's mode) into a rectangle of a figure.
%
%   P = bm_fig8_panels(fig, rect, R)
%   P = bm_fig8_panels(fig, rect, R, 'font_size', 9, ...)
%
% fig   the figure (Units are set to inches)
% rect  [x y w h] in inches, the area the panels fill
% R     one network of the burst_mode stage (bm_load_runs(...).R{1}: the
%       stage's first seed)
%
% The letters are Figure 8's of Plan 3 (Tom, 2026-11-22: patients A to C,
% modular model D to F, then this example). Until then they were A to C.
%
% G, left: the whole run: the input (mean over neurons); mean E and I
%    rates; the real part of the rightmost mode as a rate (1/s), gains
%    frozen (a wide pale line of the bursting condition's hue) and the
%    rightmost network mode of the full Jacobian (a thin black line on top
%    of it). After the step the same network WITHOUT the step continues
%    dashed: rates and both modes as they would have been with no input
%    (the stage's no-stimulation continuation, R.nostim: same network, same
%    noise; Tom, 2026-11-22, note 3). The adaptation-frozen arm is no longer
%    drawn here (note 4; S10 B shows it). A light band marks the burst of
%    the zoom; the four sampled moments of H are marked along the top by
%    triangles in their colors.
% G, right (the zoom): that burst, from 1 s before its onset to 0.5 s
%    after: rates, and the two modes on the full Jacobian's dense grid
%    (every dt_sample there); the three moments of H in this window as dots
%    of their colors on the frozen-gain mode.
% H  the eigenvalues of M as rates z = (lambda - 1)/tau_d in the complex
%    plane (complex_plane_axes: the imaginary axis is the boundary, the
%    bulk is centered at -1/tau_d), one scale across the five panels:
%    the structure (W alone: every gain 1), a quiet moment, just before the
%    burst, during it, and under stimulation, each name in its color. The
%    bulk circle has the radius predicted from that state's gains
%    (R_eff/tau_d; R_struct for W); the rightmost eigenvalue is drawn
%    larger in its moment's color. Ticks only outside the bulk, with room
%    for their labels (the imaginary-axis tick at least 1.25 bulk radii up:
%    at the first multiple above the radius the Structure plane's circle ran
%    through its label).
%
% THE COLORS of the moments (Tom, 2026-11-22, note 5: as the mockup, more
% colorful than black, purples and green): structure black, quiet blue,
% before burst orange, burst dark red, stimulation the figure's green.
%
% THE FIVE STATES (rules, not choices made by eye): the burst of the zoom
% is the one before the step whose onset is closest to the middle of the
% baseline window; quiet = 0.9 s before its onset; before = the maximum of
% the frozen-gain mode within tol (0.6 s) before onset; during = the
% burst's peak (the mean E rate's maximum); stimulation = 5 s before the
% end of the run.
%
% THE KEY: one row above G (E, I, gains frozen, full Jacobian, without
% stimulation), placed by the caller through P.legend (in G's axes there
% is no stretch free of traces once the continuation runs after the step).
%
% Options: font_size (9), annot_size (8), left (0.52), c_title (0.2), ab_room
% (0.42; inches, the room for labels, so that larger type can be given room
% on the common canvas, stim_fig8_compose), line_scale (1, multiplies every
% line width), plane_gap (0.12 in), key_room (0.3 in above G for the key).
%
% P: handles ax_A (3 axes of G's time course, top to bottom), ax_B (2, the
% zoom), ax_C (5 planes of H), legend, and info: the chosen times and the
% window of the complex plane (for the table).
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
    opts.left (1,1) double = 0.52      % inches left of G and H for their y labels
    opts.c_title (1,1) double = 0.2    % inches above H for its state names
    opts.ab_room (1,1) double = 0.42   % inches between H's names and G (its time axis)
    opts.key_room (1,1) double = 0.3   % inches above G for the key
end
if ~isfield(R, 'nostim')
    error('bm_fig8_panels:NoContinuation', ['The burst_mode file of seed %s has no no-stimulation continuation ' ...
        '(R.nostim); run the burst_mode stage again (bm_run_network writes it since 2026-11-22).'], mat2str(R.seed));
end
st = manuscript_style();
fs = opts.font_size; lw = opts.line_scale;
tau = R.tau_d; H = R.t_step;
I = R.intact; F = R.full; N = R.nostim;
ec = [st.E_color; st.I_color];                    % the manuscript's one E and one I color
cond_col = st.condition_color(R.condition);
% The two modes must be told apart where they run together (polish list,
% 2026-11-22): the frozen-gain mode is a WIDE PALE line of the condition's
% hue, the full Jacobian's a THIN BLACK line on top of it.
mode_col = [0 0 0];
fg_col = cond_col + 0.35 * (1 - cond_col);
fg_w = 2.2; mode_w = 0.7;
green = [0 0.62 0.45];                            % stimulation (Figure 8's)
state_col = {[0 0 0], [0.17 0.45 0.75], [0.95 0.55 0.1], [0.72 0.1 0.12], green};
state_name = {'Structure', 'Quiet', 'Before burst', 'Burst', 'Stimulation'};
gray_dash = [0.35 0.35 0.35];                     % the key's entry for the continuation

%% The burst of the zoom and the five states
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
gap = opts.plane_gap; left = opts.left;
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
im_lim = [-1.3 1.3] * max(rmax, max(abs(imag(zall))));
step = nice_step(rmax / 2);
% Few ticks on small panels, all well outside the largest bulk: on the real
% axis the multiples of step beyond a quarter radius past its edge, on the
% imaginary axis the first multiple above 1.25 radii (the label must clear
% the circle, not only the tick).
mult = @(lo, hi) step * (ceil(lo / step):floor(hi / step));
re_ticks = mult(cen + 1.25 * rmax, re_lim(2));
if isempty(re_ticks); re_ticks = mult(cen + 1.02 * rmax, re_lim(2)); end
im_t = mult(1.25 * rmax, im_lim(2));
if isempty(im_t)                                  % widen the window to the next multiple
    im_t = step * ceil(1.25 * rmax / step);
    im_lim = [-1.08 1.08] * im_t;
end
im_ticks = im_t(1);   % the upper one only: the lower one's label met the Re name of the last plane
axes_args = {'re_lim', re_lim, 'im_lim', im_lim, 'bulk_center', cen, 'bulk_radius', rmax, ...
    're_ticks', re_ticks, 'im_ticks', im_ticks, 'font_size', opts.annot_size, 'line_width', 0.6 * lw};
g = complex_plane_axes([], axes_args{:});
ph = pw * diff(g.y_lim) / diff(g.x_lim);
c_title = opts.c_title;
yC = y0 + 0.27;                                   % the Re name under the last plane
yAB = yC + ph + c_title + opts.ab_room;           % room for G's time axis and H's names
hA = [0.45 0.8 0.95]; gA = 0.1;
top = y0 + Ht - opts.key_room;
avail = top - yAB;
scale = (avail - 2 * gA) / sum(hA);
hA = hA * scale;
ya = [yAB + hA(3) + hA(2) + 2 * gA, yAB + hA(3) + gA, yAB];
wA = 0.60 * W - left; xB = x0 + 0.71 * W; wB = W - 0.71 * W - 0.05;

%% G, the time course
P.ax_A = gobjects(1, 3);
for r = 1:3
    P.ax_A(r) = axes(fig, 'Units', 'inches', 'Position', [x0 + left, ya(r), wA, hA(r)]);
    hold(P.ax_A(r), 'on');
end
ax = P.ax_A(1);
plot(ax, I.t, I.u, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.8 * lw);
ylim(ax, [min(I.u) - 0.05, max(I.u) + 0.1]);
yticks(ax, unique(round([0 R.dc_level], 2)));
ylabel(ax, 'Input', 'Rotation', 0, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');   % across: turned, it ran into the rate's name below
ax = P.ax_A(2);
band(ax, zoom, [0 1]);
after = N.t >= H;
for q = 2:-1:1
    plot(ax, N.t(after), N.rate(after, q), '--', 'Color', ec(q, :), 'LineWidth', 0.6 * lw);
    plot(ax, I.t, I.rate(:, q), '-', 'Color', ec(q, :), 'LineWidth', 0.8 * lw);
end
ylim(ax, [0 1]); yticks(ax, [0 1]);
ylabel(ax, 'Rate');   % 'Mean rate' ran into the mode's name below at 9 in
ax = P.ax_A(3);
naf = N.full.t >= H;
zl = mode_limits(re1, real(F.net), real(N.z1(after)), real(N.full.net(naf)));
band(ax, zoom, zl);
yline(ax, 0, ':', 'Color', [0 0 0], 'LineWidth', 0.6 * lw);
plot(ax, N.t(after), real(N.z1(after)), '--', 'Color', fg_col, 'LineWidth', 1.2 * lw);
plot(ax, I.t, re1, '-', 'Color', fg_col, 'LineWidth', fg_w * lw);
plot(ax, N.full.t(naf), real(N.full.net(naf)), '--', 'Color', mode_col, 'LineWidth', mode_w * lw);
c = F.is_coarse;
plot(ax, F.t(c), real(F.net(c)), '-', 'Color', mode_col, 'LineWidth', mode_w * lw);
ylim(ax, zl);
ylabel(ax, 'Re \lambda (s^{-1})');
xlabel(ax, 'Time (s)');
for r = 1:3
    xlim(P.ax_A(r), [I.t(1) I.t(end)]);
    if r < 3; P.ax_A(r).XTickLabel = {}; end
end
for q = 2:5
    plot(P.ax_A(3), t_state(q), zl(2), 'v', 'MarkerSize', 5.5, 'MarkerFaceColor', state_col{q}, ...
        'MarkerEdgeColor', 'none', 'Clipping', 'off');
end
% The key: proxies in the rate axes (a legend takes the objects of its own axes).
ax = P.ax_A(2);
hk = [plot(ax, NaN, NaN, '-', 'Color', ec(1, :), 'LineWidth', 0.8 * lw), ...
    plot(ax, NaN, NaN, '-', 'Color', ec(2, :), 'LineWidth', 0.8 * lw), ...
    plot(ax, NaN, NaN, '-', 'Color', fg_col, 'LineWidth', fg_w * lw), ...
    plot(ax, NaN, NaN, '-', 'Color', mode_col, 'LineWidth', mode_w * lw), ...
    plot(ax, NaN, NaN, '--', 'Color', gray_dash, 'LineWidth', 0.8 * lw)];
P.legend = legend(ax, hk, {'E', 'I', 'Gains frozen', 'Full Jacobian', 'Without stimulation'}, ...
    'Orientation', 'horizontal', 'Box', 'off', 'FontSize', opts.annot_size, 'AutoUpdate', 'off');
P.legend.ItemTokenSize = [18 8];

%% G, the zoom
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
ylim(ax, [0 1]); yticks(ax, [0 1]); ax.XTickLabel = {};
ax = P.ax_B(2);
band(ax, [on(b) off(b)], zl);
yline(ax, 0, ':', 'Color', [0 0 0], 'LineWidth', 0.6 * lw);
plot(ax, I.t(kz), re1(kz), '-', 'Color', fg_col, 'LineWidth', fg_w * lw);
kf = F.t >= zoom(1) & F.t <= zoom(2);
plot(ax, F.t(kf), real(F.net(kf)), '-', 'Color', mode_col, 'LineWidth', mode_w * lw);
for q = 2:4
    plot(ax, t_state(q), re1(idx(q - 1)), 'o', 'MarkerSize', 6.5, 'MarkerFaceColor', state_col{q}, ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
end
ylim(ax, zl); xlabel(ax, 'Time (s)');
for r = 1:2; xlim(P.ax_B(r), zoom); end

%% H
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
    if q == 5
        % The Re name under the plane, at its right end: under the real
        % axis' tick labels (complex_plane_axes' place) it is wider than the
        % stretch right of the imaginary axis and the axis ran through it.
        tr = findall(axq, 'Type', 'text');
        tr = tr(arrayfun(@(t) any(strcmp(cellstr(t.String), 'Re (s^{-1})')), tr));
        set(tr, 'String', 'Re (s^{-1})', 'Position', [cp.x_lim(2), cp.y_lim(1), 0], ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    end
    P.ax_C(q) = axq;
end

%% Style
for a = [P.ax_A P.ax_B]
    set(a, 'FontSize', fs, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 0.6 * lw, ...
        'TitleFontSizeMultiplier', 1, 'LabelFontSizeMultiplier', 1);
end
% The key in one row above G, at the left edge of its time course.
drawnow;
set([P.legend P.ax_A(1)], 'Units', 'inches');
q = P.ax_A(1).Position;
P.legend.Location = 'none';
P.legend.Position(1:2) = [q(1) - 0.1, q(2) + q(4) + 0.03];
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
