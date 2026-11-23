function h = column_dividers(target, prefix, n_col, yspan, varargin)
% COLUMN_DIVIDERS The gray lines between the condition columns of a grouped figure, placed clear of the tick labels.
%
%   h = COLUMN_DIVIDERS(sep_ax, prefix, n_col, yspan, 'Color', c, 'LineWidth', w)
%   COLUMN_DIVIDERS(fig)
%
% The first form draws n_col - 1 vertical lines into sep_ax, an invisible
% axes over the whole figure with XLim and YLim [0 1] (so its data units
% are normalized figure units), from yspan(1) to yspan(2), and places them.
% The axes of column c are those whose Tag is <prefix>_r<row>_c<c>
% (fig_main2_dynamics' 'dynamics', fig_main3_stability's 'stability').
% Name-value pairs go to plot (Color, LineWidth).
%
% The second form places again every divider of fig. Call it whenever the
% figure's size, the axes' positions or the type have changed since the
% lines were drawn: fig_grouped_main copies the panels onto a canvas of
% another size, where type in points takes a different share of the
% width. Each line keeps what it needs in its UserData, so it survives
% copyobj.
%
% WHERE A DIVIDER GOES (Tom, 2026-11-22, on Figures 2 and 3: the dividers
% must clear the y tick labels of every row). Between columns c and c + 1,
% halfway between the rightmost drawn extent of column c (over every row:
% the plot box plus its tight inset, its title, and the box of any legend
% or colour bar of those axes) and the leftmost drawn extent of column
% c + 1 (the plot box less its tight inset, which holds its y tick labels
% and y label). Only what lies beside the line's own height counts, so a
% column title above the line's top does not move it. Both are measured
% on the figure as drawn, in pixels, so the place follows the
% canvas width, the type sizes and the tick labels; no offset is typed in.
% When the columns leave no gap the line goes to the left edge of column
% c + 1's tick labels and a warning (column_dividers:NoGap) says so.
%
% See also: fig_main2_dynamics, fig_main3_stability, fig_grouped_main

if nargin == 1
    fig = ancestor(target, 'figure');
    lines = findall(fig, 'Type', 'line', 'Tag', 'column_divider');
    if isempty(lines); h = lines; return; end
    drawnow;
    for k = 1:numel(lines)
        place(lines(k));
    end
    h = lines;
    return
end
sep_ax = target;
h = gobjects(n_col - 1, 1);
for c = 1:n_col - 1
    h(c) = plot(sep_ax, [0 0], yspan, varargin{:}, 'Tag', 'column_divider', ...
        'UserData', struct('prefix', prefix, 'left_column', c));
end
drawnow;
for c = 1:n_col - 1
    place(h(c));
end
end

%% ------------------------------------------------------------------------
function place(line)
u = line.UserData;
fig = ancestor(line, 'figure');
ax = line.Parent;
p = getpixelposition(ax, true);
% The line's own height, in figure pixels: only what lies beside it counts.
ys = p(2) + (sort(line.YData) - ax.YLim(1)) / diff(ax.YLim) * p(4);
left = column_boxes(fig, u.prefix, u.left_column, ys);
right = column_boxes(fig, u.prefix, u.left_column + 1, ys);
if isempty(left) || isempty(right)
    warning('column_dividers:NoColumn', 'No axes tagged %s_r*_c%d or _c%d beside the divider: it is left where it is.', ...
        u.prefix, u.left_column, u.left_column + 1);
    return
end
r = max(left(:, 3));       % right edge of column c, tick labels, labels and keys included
l = min(right(:, 1));      % left edge of column c + 1's tick labels and y label
if r < l
    x = (r + l) / 2;
else
    x = l;
    warning('column_dividers:NoGap', ['Columns %d and %d of %s overlap by %.0f px: the divider is put at ' ...
        'the left edge of column %d''s tick labels.'], u.left_column, u.left_column + 1, u.prefix, r - l, u.left_column + 1);
end
line.XData = ax.XLim(1) + (x - p(1)) / p(3) * diff(ax.XLim) * [1 1];
end

function B = column_boxes(fig, prefix, c, ys)
% [x_left y_bottom x_right y_top] in figure pixels of what column c draws
% between the heights ys: per visible axes its plot box with tick labels and
% axis labels (the tight inset with the title hidden) and, apart, its title
% (a title may be wider than its column and sit above the divider's top);
% per legend or colour bar of those axes its box.
pat = ['^' regexptranslate('escape', prefix) '_r\d+_c' num2str(c) '$'];
ax = findall(fig, '-isa', 'matlab.graphics.axis.Axes', 'Visible', 'on');
ax = ax(~cellfun(@isempty, regexp(get(ax, {'Tag'}), pat, 'once')));
B = zeros(0, 4);
for k = 1:numel(ax)
    a = ax(k);
    p = getpixelposition(a, true);
    t = a.Title; shown = t.Visible;
    t.Visible = 'off';
    units = a.Units; a.Units = 'pixels'; ti = a.TightInset; a.Units = units;
    t.Visible = shown;
    B(end + 1, :) = [p(1) - ti(1), p(2) - ti(2), p(1) + p(3) + ti(3), p(2) + p(4) + ti(4)]; %#ok<AGROW>
    if strcmp(shown, 'on') && ~isempty(t.String)
        tu = t.Units; t.Units = 'pixels'; e = t.Extent; t.Units = tu;
        B(end + 1, :) = [p(1) + e(1), p(2) + e(2), p(1) + e(1) + e(3), p(2) + e(2) + e(4)]; %#ok<AGROW>
    end
end
for o = findall(fig, 'Type', 'legend', '-or', 'Type', 'colorbar')'
    if strcmp(o.Visible, 'off') || isempty(o.Axes) || ~any(o.Axes == ax); continue; end
    q = getpixelposition(o, true);
    B(end + 1, :) = [q(1), q(2), q(1) + q(3), q(2) + q(4)]; %#ok<AGROW>
end
B = B(B(:, 4) > ys(1) & B(:, 2) < ys(2), :);   % beside the line
end
