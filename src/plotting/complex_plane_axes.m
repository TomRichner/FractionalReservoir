function cp = complex_plane_axes(ax, opts)
% COMPLEX_PLANE_AXES Axes for eigenvalues in the complex plane: through the origin, 1:1, with an optional break.
%
%   cp = COMPLEX_PLANE_AXES(ax, 're_lim', [-32 12], 'im_lim', [-22 22])
%   cp = COMPLEX_PLANE_AXES(ax, ..., 'far_lim', [28 42])                 a break in the real axis
%   cp = COMPLEX_PLANE_AXES(ax, ..., 'bulk_center', -10, 'bulk_radius', 17.5)   no ticks inside the bulk
%   cp = COMPLEX_PLANE_AXES([], ...)                                     geometry and mapping only, nothing drawn
%
% Tom's preferences for complex-plane eigenvalue plots (2026-11-20):
%
%   - The real axis is drawn along Im = 0 and the imaginary axis along
%     Re = 0, with no box and none of MATLAB's edge axes (the axes object is
%     switched off and the two lines are drawn by hand, as plot_eigenvalues
%     and the discs of fig_main1_concepts already do). The imaginary axis is
%     then also the stability boundary: no separate zero line is needed.
%   - One unit of Re is as long as one unit of Im, in the main window and in
%     the segment beyond a break.
%   - Tick marks and labels only OUTSIDE the bulk, never at the origin.
%   - The axis names sit at the ends of the axes, not centered on a box.
%
% THE BREAK. With far_lim = [f1 f2], the stretch of the real axis between
% re_lim(2) and f1 is left out: the main window [re_lim(1), re_lim(2)] is
% followed by a gap and then by the far segment [f1, f2], at the same scale.
% The real axis line stops at the gap, with a pair of slanted strokes on
% each side of it. Because the plotted x is then not the data x, EVERYTHING
% the caller draws must go through the mapping this function returns:
%
%   plot(ax, cp.x(real(z)), imag(z), '.')
%
% cp.x ERRORS (complex_plane_axes:InBreak) when any value lies strictly
% inside the removed stretch, so a point there cannot vanish silently, and
% cannot be drawn at a wrong place either. NaN passes through. A value left
% of re_lim(1) or right of the last limit is mapped as usual and then
% clipped by the axes like on any plot; cp.check(re) counts those. A LINE
% whose vertices lie on both sides of the break would be drawn straight
% across the gap: split it (cp.x cannot see that).
%
% THE SAME SCALE IN SEVERAL PANELS. The geometry depends only on the
% options, so calling this with the same options on axes of the same size
% gives panels with identical limits, identical ticks and the same length
% per unit. The limits are set explicitly and the plot box is given the
% aspect ratio of the limits, so a panel does not change size with its
% data. The far segment is part of that geometry: it takes its room in
% every panel, and by default its axis segment and ticks are drawn in every
% panel too (an empty segment says "nothing is out here"); 'draw_far',
% false leaves the room empty in a panel, without changing its size.
%
% OPTIONS
%   re_lim       [lo hi]  the main window on the real axis (required)
%   im_lim       [lo hi]  the imaginary axis (required)
%   far_lim      [f1 f2]  the segment beyond the break, f1 > re_lim(2);
%                         [] (default): no break
%   gap          width of the gap, in axis units (default: 5% of the
%                plotted width)
%   draw_far     draw the far segment's axis, break marks and ticks (true)
%   bulk_center, bulk_radius   keep ticks out of the disk: the real axis is
%                kept clear on real(center) +/- radius, the imaginary axis
%                on imag(center) +/- radius (the whole height of the disk,
%                not only where the disk crosses the axis)
%   re_clear, im_clear   [lo hi] intervals to keep clear, instead of or in
%                addition to the bulk (the union is kept clear)
%   tick_step    spacing of the ticks, multiples of it (default: a round
%                number near a sixth of the larger span)
%   re_ticks, im_ticks   explicit tick values, instead of the multiples. A
%                value inside a clear interval, in the removed stretch or
%                outside the limits is an error, not a silent drop.
%   tick_min_sep the smallest plotted distance between neighboring ticks on
%                the real axis (default 0.6 of tick_step); a far tick closer
%                than this to the last main tick is not drawn
%   re_label, im_label   the names at the ends ('Re (s^{-1})', 'Im (s^{-1})';
%                '' for none)
%   re_label_side  'below' (default): the Re name under the tick labels at
%                the right end of the real axis, clear of anything drawn ON
%                the real axis, which is where the outliers of these plots
%                are. 'above': just above the axis line at its right end
%                (it then collides with a mark near the end of the axis).
%   font_size (12), line_width (1), color ([0 0 0])
%
% RETURNS cp, a struct:
%   x        @(re) plotted x of a real part (errors inside the break)
%   x_inv    @(xp) the real part of a plotted x (NaN inside the gap)
%   in_break @(re) logical, true strictly inside the removed stretch
%   check    @(re) struct of counts: n_in_break, n_left, n_right
%   re_ticks, im_ticks   the tick values drawn (data units)
%   x_lim, y_lim         the plotted limits (y_lim is im_lim)
%   re_lim, far_lim, gap, re_clear, im_clear, tick_step
%   handles  the graphics objects drawn (empty when ax is [])
%
% Where to put labels of your own: the tick labels are below the real axis
% and left of the imaginary axis; the Re name is under the tick labels at the
% right end of the real axis and the Im name right of the top of the
% imaginary axis. Calling
% this again on the same axes replaces what it drew before.
%
% See also: plot_eigenvalues, fig_main1_concepts, plot_box_position

arguments
    ax
    opts.re_lim (1,2) double
    opts.im_lim (1,2) double
    opts.far_lim double = []
    opts.gap double = []
    opts.draw_far (1,1) logical = true
    opts.bulk_center double = []
    opts.bulk_radius double = []
    opts.re_clear double = []
    opts.im_clear double = []
    opts.tick_step double = []
    opts.re_ticks double = []
    opts.im_ticks double = []
    opts.tick_min_sep double = []
    opts.re_label (1,:) char = 'Re (s^{-1})'
    opts.im_label (1,:) char = 'Im (s^{-1})'
    opts.re_label_side (1,:) char {mustBeMember(opts.re_label_side, {'below', 'above'})} = 'below'
    opts.font_size (1,1) double = 12
    opts.line_width (1,1) double = 1
    opts.color (1,3) double = [0 0 0]
end

%% Geometry
if ~isfield(opts, 're_lim') || ~isfield(opts, 'im_lim')
    error('complex_plane_axes:MissingLimits', 're_lim and im_lim are required.');
end
G.re_lim = opts.re_lim; G.im_lim = opts.im_lim;
if ~(G.re_lim(1) < G.re_lim(2)) || ~(G.im_lim(1) < G.im_lim(2))
    error('complex_plane_axes:BadLimits', 're_lim and im_lim must be [lo hi] with lo < hi.');
end
if ~(G.re_lim(1) <= 0 && G.re_lim(2) >= 0 && G.im_lim(1) <= 0 && G.im_lim(2) >= 0)
    error('complex_plane_axes:OriginOutside', 'The origin must lie inside re_lim and im_lim: the axes are drawn through it.');
end
G.far_lim = opts.far_lim(:).';
has_far = ~isempty(G.far_lim);
if has_far
    if numel(G.far_lim) ~= 2 || ~(G.far_lim(1) < G.far_lim(2)) || ~(G.far_lim(1) > G.re_lim(2))
        error('complex_plane_axes:BadFarLimits', 'far_lim must be [f1 f2] with re_lim(2) < f1 < f2.');
    end
    far_w = diff(G.far_lim);
    if isempty(opts.gap)
        G.gap = 0.05 * (diff(G.re_lim) + far_w) / 0.95;
    else
        G.gap = opts.gap;
    end
    if ~(isscalar(G.gap) && G.gap > 0)
        error('complex_plane_axes:BadGap', 'gap must be a positive scalar.');
    end
    G.x_lim = [G.re_lim(1), G.re_lim(2) + G.gap + far_w];
else
    G.gap = 0;
    G.x_lim = G.re_lim;
end
G.y_lim = G.im_lim;

%% Intervals to keep clear of ticks
G.re_clear = zeros(0, 2); G.im_clear = zeros(0, 2);
if ~isempty(opts.bulk_center) || ~isempty(opts.bulk_radius)
    if ~isscalar(opts.bulk_center) || ~isscalar(opts.bulk_radius) || ~(opts.bulk_radius >= 0)
        error('complex_plane_axes:BadBulk', 'bulk_center and bulk_radius are given together, as scalars.');
    end
    G.re_clear(end + 1, :) = real(opts.bulk_center) + [-1 1] * opts.bulk_radius;
    G.im_clear(end + 1, :) = imag(opts.bulk_center) + [-1 1] * opts.bulk_radius;
end
if ~isempty(opts.re_clear); G.re_clear(end + 1, :) = opts.re_clear(:).'; end
if ~isempty(opts.im_clear); G.im_clear(end + 1, :) = opts.im_clear(:).'; end

%% Ticks
if isempty(opts.tick_step)
    G.tick_step = round_step(max(diff(G.x_lim), diff(G.y_lim)) / 6);
else
    G.tick_step = opts.tick_step;
end
if isempty(opts.tick_min_sep); min_sep = 0.6 * G.tick_step; else; min_sep = opts.tick_min_sep; end
if isempty(opts.re_ticks)
    re_main = multiples(G.tick_step, G.re_lim);
    re_main = re_main(~in_clear(re_main, G.re_clear) & re_main ~= 0);
    re_far = [];
    if has_far && opts.draw_far
        re_far = multiples(G.tick_step, G.far_lim);
        re_far = re_far(~in_clear(re_far, G.re_clear));
        if ~isempty(re_main)                                % not on top of the last tick before the break
            re_far = re_far(map_re(G, re_far) - map_re(G, re_main(end)) >= min_sep);
        end
    end
    G.re_ticks = [re_main, re_far];
else
    t = opts.re_ticks(:).';
    bad = in_clear(t, G.re_clear) | is_in_break(G, t) | t < G.re_lim(1) | t > last_re(G);
    if any(bad)
        error('complex_plane_axes:BadTick', 're_ticks %s: inside a clear interval, in the removed stretch or outside the limits.', mat2str(t(bad)));
    end
    G.re_ticks = t;
end
if isempty(opts.im_ticks)
    t = multiples(G.tick_step, G.im_lim);
    G.im_ticks = t(~in_clear(t, G.im_clear) & t ~= 0);
else
    t = opts.im_ticks(:).';
    bad = in_clear(t, G.im_clear) | t < G.im_lim(1) | t > G.im_lim(2);
    if any(bad)
        error('complex_plane_axes:BadTick', 'im_ticks %s: inside a clear interval or outside the limits.', mat2str(t(bad)));
    end
    G.im_ticks = t;
end

%% What the caller gets
cp = G;
cp.x = @(re) map_re(G, re);
cp.x_inv = @(xp) unmap_re(G, xp);
cp.in_break = @(re) is_in_break(G, re);
cp.check = @(re) check_re(G, re);
cp.handles = gobjects(1, 0);
if isempty(ax); return; end

%% Draw
delete(findall(ax, 'Tag', 'complex_plane_axes'));
hold(ax, 'on');
set(ax, 'XLim', G.x_lim, 'YLim', G.y_lim, 'DataAspectRatio', [1 1 1], ...
    'PlotBoxAspectRatio', [diff(G.x_lim), diff(G.y_lim), 1], 'Visible', 'off', 'Clipping', 'on');
set(ax.Title, 'Visible', 'on');                       % a column title survives the hidden axes
tl = 0.018 * diff(G.y_lim);                           % tick length, and the unit of every offset below
line_args = {'Color', opts.color, 'LineWidth', opts.line_width, 'Tag', 'complex_plane_axes', ...
    'HandleVisibility', 'off', 'Clipping', 'off'};
text_args = {'FontSize', opts.font_size, 'Color', opts.color, 'Tag', 'complex_plane_axes', ...
    'HandleVisibility', 'off', 'Interpreter', 'tex'};
h = gobjects(1, 0);

% The two axes. The real axis stops at the gap.
h(end + 1) = plot(ax, G.re_lim, [0 0], '-', line_args{:});
h(end + 1) = plot(ax, [0 0], G.im_lim, '-', line_args{:});
x_end = G.re_lim(2);
if has_far && opts.draw_far
    xa = G.re_lim(2); xb = G.re_lim(2) + G.gap; x_end = G.x_lim(2);
    h(end + 1) = plot(ax, [xb, x_end], [0 0], '-', line_args{:});
    s = 1.3 * tl;                                     % the break marks: a slanted stroke on each side of the gap
    for x0 = [xa, xb]
        h(end + 1) = plot(ax, x0 + [-0.4 0.4] * s, [-1 1] * s, '-', line_args{:}); %#ok<AGROW>
    end
end

% Ticks: down from the real axis, left from the imaginary axis.
for t = G.re_ticks
    xp = map_re(G, t);
    h(end + 1) = plot(ax, [xp xp], [0 -tl], '-', line_args{:}); %#ok<AGROW>
    h(end + 1) = text(ax, xp, -1.5 * tl, tick_label(t), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', text_args{:}); %#ok<AGROW>
end
for t = G.im_ticks
    h(end + 1) = plot(ax, [0 -tl], [t t], '-', line_args{:}); %#ok<AGROW>
    h(end + 1) = text(ax, -1.6 * tl, t, tick_label(t), 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', text_args{:}); %#ok<AGROW>
end

% The names, at the ends.
if ~isempty(opts.re_label)
    if strcmp(opts.re_label_side, 'below')              % a blank first line: under the tick labels, whatever the scale
        h(end + 1) = text(ax, x_end, -1.5 * tl, {' ', opts.re_label}, 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', text_args{:});
    else
        h(end + 1) = text(ax, x_end, 1.2 * tl, opts.re_label, 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'bottom', text_args{:});
    end
end
if ~isempty(opts.im_label)
    h(end + 1) = text(ax, 1.2 * tl, G.im_lim(2), opts.im_label, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', text_args{:});
end
cp.handles = h;
end

%% ------------------------------------------------------------------------
function xp = map_re(G, re)
% Data real part to plotted x. Errors inside the removed stretch.
bad = is_in_break(G, re);
if any(bad(:))
    v = re(bad);
    error('complex_plane_axes:InBreak', ...
        '%d value(s) lie inside the removed stretch of the real axis (%g, %g), for example %g. Move the break.', ...
        nnz(bad), G.re_lim(2), G.far_lim(1), v(1));
end
xp = re;
if ~isempty(G.far_lim)
    far = re >= G.far_lim(1);
    xp(far) = re(far) - G.far_lim(1) + G.re_lim(2) + G.gap;
end
end

function re = unmap_re(G, xp)
% Plotted x back to the real part; NaN inside the gap.
re = xp;
if ~isempty(G.far_lim)
    x_far = G.re_lim(2) + G.gap;
    far = xp >= x_far;
    re(far) = xp(far) - x_far + G.far_lim(1);
    re(xp > G.re_lim(2) & xp < x_far) = NaN;
end
end

function tf = is_in_break(G, re)
if isempty(G.far_lim)
    tf = false(size(re));
else
    tf = re > G.re_lim(2) & re < G.far_lim(1);
end
end

function hi = last_re(G)
if isempty(G.far_lim); hi = G.re_lim(2); else; hi = G.far_lim(2); end
end

function c = check_re(G, re)
re = re(~isnan(re));
c.n_in_break = nnz(is_in_break(G, re));
c.n_left = nnz(re < G.re_lim(1));
c.n_right = nnz(re > last_re(G));
end

function t = multiples(step, lim)
t = (ceil(lim(1) / step - 1e-9):floor(lim(2) / step + 1e-9)) * step;
end

function tf = in_clear(t, intervals)
% Strictly inside any of the intervals (a tick on the edge of one is outside the bulk).
tf = false(size(t));
for k = 1:size(intervals, 1)
    tf = tf | (t > intervals(k, 1) & t < intervals(k, 2));
end
end

function step = round_step(raw)
% The round number (1, 2 or 5 times a power of ten) nearest to raw from below.
p = 10 ^ floor(log10(raw));
m = raw / p;
if m >= 5; step = 5 * p; elseif m >= 2; step = 2 * p; else; step = p; end
end

function s = tick_label(t)
s = sprintf('%g', abs(t));
if t < 0; s = [char(8722) s]; end                      % a minus sign, not a hyphen
end
