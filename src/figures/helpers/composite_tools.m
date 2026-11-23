function T = composite_tools()
% COMPOSITE_TOOLS The shared steps of a composite figure, as a struct of function handles.
%
%   T = COMPOSITE_TOOLS();
%   fig = T.canvas(1100, 900, cfg.visible);
%   C   = T.copy(src_fig, fig);            % axes, legends, colour bars, in pixels
%   T.place(C.axes, [x y w h], C.legends); % the block's plot boxes fill the rectangle
%   T.letter(C.axes(1), 'A');
%   T.style(fig);
%   out = composite_save(fig, cfg, out_dir, tag, heading, caption, parts, sources);
%
% A composite figure of package F of Plan 2 (fig_main7_memory,
% fig_supp1_grouped, _2, _4, _5, _9) calls the functions that draw its
% panels, invisibly and unsaved, copies their axes onto one canvas and
% places them in pixels. So a panel is drawn by ONE function, the one that
% owns it, and a change of colour or title there reaches every composite.
% Nothing here simulates or reads run data.
%
% BESIDE composite_figure (package D). composite_figure stacks WHOLE source
% figures in rows at their own size. The composites here take SOME of a
% source figure's axes (the reconstruction axes of fig_memory_capacity
% without its other two, the trace axes of the memory example) and put them
% at chosen places beside axes of other figures, which composite_figure does
% not do; hence these steps. The figure is saved by the same
% composite_save (save_figure_stable, with its export of figures larger
% than the screen, and the same <tag>_table.md), and the letters look the
% same as composite_letter's.
%
%   canvas(W, H, visible)      a white figure canvas_w wide (manuscript_style)
%                              and H pixels tall, with the size recorded
%                              for composite_save ('composite_size'). W is
%                              the width the caller lays its panels out
%                              on, in pixels; fit stretches that layout
%                              across the canvas. Keep it under 0.9 of the
%                              screen height where the panels allow: a
%                              taller figure is exported by printing, on an
%                              uncropped canvas.
%   copy(src, dst, attached)   copies the axes, legends and colour bars of
%                              the figure src (or the handles src, of one
%                              figure) to dst; with handles, the legends and
%                              colour bars that belong to those axes come
%                              too unless attached is false. Each copy is in
%                              pixel units at its pixel position in the
%                              source. Returns a struct: all, axes, legends,
%                              colorbars (the copies, axes in the source's
%                              order of findall).
%   place(axs, rect, extras)   moves the axes so that the bounding box of
%                              their plot boxes is rect (pixels), keeping
%                              their arrangement; extras (colour bars,
%                              legends) follow the same map, a legend
%                              keeping its size.
%   grid(axs)                  the axes as a rows x columns array, read from
%                              their positions (top row first).
%   letter(ax, ch, dx, dy)     the panel letter '(A)' at (dx, dy) pixels from
%                              the top-left corner of the plot box of ax
%                              (default -70, 30), letter_fs, normal weight.
%   style(fig)                 the house type, by role (manuscript_style):
%                              tick labels, legends and free text tick_fs,
%                              axis and colour-bar labels label_fs, titles
%                              title_fs in normal weight; axes lines
%                              axis_lw; no boxes around legends.
%   fit(fig)                   the common canvas (Plan 2, 2026-11-22): the
%                              figure becomes canvas_w wide and what is laid
%                              out in pixels on the canvas's own width is
%                              stretched across it, horizontally only. Axes
%                              take the new left edge and width; a legend
%                              or colour bar keeps its size while its anchor
%                              follows (a legend its centre, a colour bar
%                              its left edge); type is never scaled. Call
%                              it once the layout is done and before the
%                              letters, so that a letter keeps its offset
%                              in pixels from its panel.
%   axes_table(ax, what)       the numbers drawn in an axes as markdown, for
%                              a section of composite_save's table: a line
%                              saying what it is, then one row per line,
%                              scatter and histogram.
%   grid_table(what, G, row_levels, col_levels, row_name, col_name)
%                              a matrix on two axes as a markdown table.
%   condition_title(name, two_line)
%                              the condition's title in sentence case ('No
%                              adaptation'), on two lines when asked.
%   retitle(axs, names, two_line)
%                              titles each axes with its condition, in
%                              sentence case and the condition's colour
%                              (manuscript_style), in normal weight.
%   close(figs)                closes the ones that still exist.
%   sentence(strings)          legend entries or titles with their first
%                              word capitalized ('noise-free' becomes
%                              'Noise-free'); text that opens with a
%                              symbol ('n_a = 3', 'lambda_1') is left.
%   margin(fig, px)            keeps the canvas down to px pixels above its
%                              bottom edge in the cropped export.
%                              exportgraphics crops to the drawn pixels,
%                              so a bottom row's x labels end 2 px from
%                              the image edge (S4, 2026-11-22). A short
%                              line at mid-width, px above the bottom, in
%                              a white that is not the background's
%                              (0.99), is counted by the crop and is not
%                              seen. A line in the background's own white
%                              is not counted (tried).
%
% See also: composite_save, composite_figure, fig_main7_memory, fig_grouped_main

T = struct('canvas', @canvas, 'copy', @copy_graphics, 'place', @place_block, 'grid', @axes_grid, ...
    'letter', @panel_letter, 'style', @house_style, ...
    'axes_table', @axes_table, 'grid_table', @grid_table, ...
    'condition_title', @condition_title, 'retitle', @retitle, 'close', @close_valid, ...
    'sentence', @sentence_case, 'margin', @crop_margin, 'fit', @fit_width);
end

%% ------------------------------------------------------------------------
function fig = canvas(W, H, visible)
% The window is made canvas_w wide at once; W is the width the caller lays
% out on, which fit stretches across canvas_w. (Widening the window later
% was tried, 2026-11-22: an export in a run of several cut the PNG at the
% old width, the right-hand legends cut through.)
vis = 'off';
if visible; vis = 'on'; end
st = manuscript_style();
W1 = round(st.canvas_w * get(groot, 'ScreenPixelsPerInch'));
fig = figure('Visible', vis, 'Color', 'w', 'Units', 'pixels', 'Position', [40 40 W1 H]);
% composite_save puts this size back before it saves (a window larger than
% the screen can be shrunk by MATLAB at a moment of its choosing).
setappdata(fig, 'composite_size', [W1 H]);
setappdata(fig, 'composite_layout_width', W);
drawnow;
end

function C = copy_graphics(src, dst, attached)
if nargin < 3; attached = true; end
if isscalar(src) && isgraphics(src, 'figure')
    h = findall(src, '-isa', 'matlab.graphics.axis.Axes', '-or', 'Type', 'legend', '-or', 'Type', 'colorbar');
else
    h = src(:);
    if attached
        owner = ancestor(h(1), 'figure');
        extra = findall(owner, 'Type', 'legend', '-or', 'Type', 'colorbar');
        keep = false(size(extra));
        for k = 1:numel(extra)
            keep(k) = ~isempty(extra(k).Axes) && any(h == extra(k).Axes);
        end
        h = [h; extra(keep)];
    end
end
drawnow;
rect = zeros(numel(h), 4);
for k = 1:numel(h); rect(k, :) = getpixelposition(h(k), true); end
new = copyobj(h, dst);
is_ax = false(size(new)); is_lg = is_ax; is_cb = is_ax;
for k = 1:numel(new)
    is_ax(k) = isa(new(k), 'matlab.graphics.axis.Axes');
    is_lg(k) = isa(new(k), 'matlab.graphics.illustration.Legend');
    is_cb(k) = isa(new(k), 'matlab.graphics.illustration.ColorBar');
    new(k).Units = 'pixels';
    if is_ax(k); new(k).PositionConstraint = 'innerposition'; end
    new(k).Position = rect(k, :);
end
C = struct('all', new, 'axes', new(is_ax), 'legends', new(is_lg), 'colorbars', new(is_cb));
end

function place_block(axs, rect, extras)
if nargin < 3; extras = gobjects(0); end
axs = axs(:); extras = extras(:);
set([axs; extras], 'Units', 'pixels');
P = zeros(numel(axs), 4);
for k = 1:numel(axs); P(k, :) = axs(k).Position; end
b = [min(P(:, 1)), min(P(:, 2)), max(P(:, 1) + P(:, 3)), max(P(:, 2) + P(:, 4))];
sx = rect(3) / (b(3) - b(1)); sy = rect(4) / (b(4) - b(2));
map = @(p) [rect(1) + (p(1) - b(1)) * sx, rect(2) + (p(2) - b(2)) * sy, p(3) * sx, p(4) * sy];
for k = 1:numel(axs); axs(k).Position = map(P(k, :)); end
for k = 1:numel(extras)
    p = extras(k).Position;
    q = map(p);
    if isa(extras(k), 'matlab.graphics.illustration.Legend')
        % A legend keeps its size; its centre follows the map.
        q = [q(1) + q(3) / 2 - p(3) / 2, q(2) + q(4) / 2 - p(4) / 2, p(3), p(4)];
    end
    extras(k).Position = q;
end
end

function G = axes_grid(axs)
% Rows by the bottom edge (top row first), columns by the left edge. Axes
% whose edges agree within 5 px share a row or a column.
axs = axs(:);
P = zeros(numel(axs), 4);
for k = 1:numel(axs); P(k, :) = getpixelposition(axs(k), true); end
rows = cluster_1d(-P(:, 2));
cols = cluster_1d(P(:, 1));
G = gobjects(max(rows), max(cols));
for k = 1:numel(axs)
    assert(~isgraphics(G(rows(k), cols(k))), 'composite_tools:NotAGrid', ...
        'Two axes sit in row %d, column %d: the axes are not laid out as a grid.', rows(k), cols(k));
    G(rows(k), cols(k)) = axs(k);
end
end

function id = cluster_1d(v)
[s, order] = sort(v);
step = [1; diff(s) > 5];
id = zeros(size(v));
id(order) = cumsum(step);
end

function h = panel_letter(ax, ch, dx, dy)
if nargin < 3 || isempty(dx); dx = -70; end
if nargin < 4 || isempty(dy); dy = 30; end
p = getpixelposition(ax);
st = manuscript_style();
h = text(ax, dx / p(3), 1 + dy / p(4), sprintf('(%s)', ch), 'Units', 'normalized', 'Clipping', 'off', ...
    'FontSize', st.letter_fs, 'FontWeight', 'normal', 'Interpreter', 'none', 'Tag', 'composite_panel_letter', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end

function house_style(fig)
% Type by role (Plan2_canvas_standard.md): everything at tick_fs first, then
% each piece of type its role's size by the shared canvas_type (axis names
% and colour-bar labels label_fs, titles title_fs, panel letters
% letter_fs; the roles are canvas_audit's, the ones test_canvas_standard
% checks), each set on the label itself (the multipliers are 1).
st = manuscript_style();
ax = findall(fig, '-isa', 'matlab.graphics.axis.Axes');
set(ax, 'FontSize', st.tick_fs, 'LineWidth', st.axis_lw, 'TitleFontSizeMultiplier', 1, 'LabelFontSizeMultiplier', 1);
set(findall(fig, '-property', 'FontSize'), 'FontSize', st.tick_fs);
for k = 1:numel(ax)
    ax(k).Title.FontWeight = 'normal';
    if isprop(ax(k), 'Subtitle'); ax(k).Subtitle.FontWeight = 'normal'; end
end
canvas_type(fig);
set(findall(fig, 'Type', 'legend'), 'Box', 'off');
set(findall(fig, 'Type', 'colorbar'), 'Box', 'off', 'LineWidth', st.axis_lw);
end

function fit_width(fig)
% Stretch a layout made on the canvas's own width across canvas_w (see the header).
st = manuscript_style();
W1 = round(st.canvas_w * get(groot, 'ScreenPixelsPerInch'));
fig.Units = 'pixels';
W0 = getappdata(fig, 'composite_layout_width');
if isempty(W0); W0 = fig.Position(3); end
s = W1 / W0;
setappdata(fig, 'composite_layout_width', W1);      % a second call does nothing
if abs(s - 1) < 1e-9; return; end
h = findall(fig, '-isa', 'matlab.graphics.axis.Axes', '-or', 'Type', 'legend', '-or', 'Type', 'colorbar');
% A legend or colour bar that MATLAB places beside its axes follows them.
placed = arrayfun(@(o) ~isa(o, 'matlab.graphics.axis.Axes') && ~any(strcmp(o.Location, {'none', 'manual'})), h);
h = h(~placed);
units = get(h, {'Units'});
P = zeros(numel(h), 4);
for k = 1:numel(h)
    h(k).Units = 'pixels';
    P(k, :) = h(k).Position;
end
sz = getappdata(fig, 'composite_size');
if ~isempty(sz) && any(abs(fig.Position(3:4) - sz) > 1e-6); fig.Position(3:4) = sz; end
drawnow;
for k = 1:numel(h)
    p = P(k, :);
    if isa(h(k), 'matlab.graphics.axis.Axes')
        q = [p(1) * s, p(2), p(3) * s, p(4)];
    elseif isa(h(k), 'matlab.graphics.illustration.Legend')
        q = [(p(1) + p(3) / 2) * s - p(3) / 2, p(2), p(3), p(4)];
    else
        q = [p(1) * s, p(2), p(3), p(4)];
    end
    h(k).Position = q;
    if ~strcmp(units{k}, 'pixels'); h(k).Units = units{k}; end
end
drawnow;
end

function txt = axes_table(ax, what)
% One markdown row per drawn series. Up to 24 points are listed; a longer
% series is summarised (its count, its x range, and the median and range
% of y).
lines = {sprintf('%s. Read back from the axes drawn; x: %s; y: %s.', what, ...
    plain(ax.XLabel.String), plain(ax.YLabel.String)), ''};
lines{end+1} = '| series | n | x | y |';
lines{end+1} = '|---|---|---|---|';
h = flipud(findall(ax, 'Type', 'line', '-or', 'Type', 'scatter', '-or', 'Type', 'histogram'));
n_rows = 0;
for k = 1:numel(h)
    if isa(h(k), 'matlab.graphics.chart.primitive.Histogram')
        x = h(k).BinEdges(1:end-1) + diff(h(k).BinEdges) / 2;
        y = h(k).Values;
    else
        x = double(h(k).XData); y = double(h(k).YData);
    end
    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);
    if isempty(x); continue; end
    name = plain(h(k).DisplayName);
    if isempty(name); name = plain(h(k).Tag); end
    if isempty(name); name = sprintf('%s %d', h(k).Type, k); end
    if numel(x) <= 24
        lines{end+1} = sprintf('| %s | %d | %s | %s |', name, numel(x), nums(x), nums(y)); %#ok<AGROW>
    else
        lines{end+1} = sprintf('| %s | %d | %.4g to %.4g | median %.4g, range %.4g to %.4g |', ...
            name, numel(x), min(x), max(x), median(y), min(y), max(y)); %#ok<AGROW>
    end
    n_rows = n_rows + 1;
end
if n_rows == 0; lines{end+1} = '| (no line, scatter or histogram) | | | |'; end
lines{end+1} = '';
txt = strjoin(lines, newline);
end

function txt = grid_table(what, G, row_levels, col_levels, row_name, col_name)
lines = {sprintf('%s. Rows: %s; columns: %s.', what, row_name, col_name), '', ...
    ['| | ' strjoin(arrayfun(@(v) sprintf('%.3g', v), col_levels, 'UniformOutput', false), ' | ') ' |'], ...
    ['|---|' repmat('---|', 1, numel(col_levels))]};
for r = numel(row_levels):-1:1
    lines{end+1} = sprintf('| %.3g | %s |', row_levels(r), ...
        strjoin(arrayfun(@(v) sprintf('%.3g', v), G(r, :), 'UniformOutput', false), ' | ')); %#ok<AGROW>
end
lines{end+1} = '';
txt = strjoin(lines, newline);
end

function s = nums(v)
s = strjoin(arrayfun(@(x) sprintf('%.4g', x), v(:)', 'UniformOutput', false), ', ');
end

function s = plain(x)
% A label as one line of plain text, safe in a markdown table cell.
if iscell(x) || isstring(x); x = strjoin(cellstr(x), ' '); end
if size(x, 1) > 1; x = strjoin(cellstr(x), ' '); end
s = strtrim(strrep(char(x), '|', '/'));
end

function t = condition_title(name, two_line)
if nargin < 2; two_line = false; end
t = grouped_condition_label(name, 'sentence');
t = t{1};
if two_line
    k = strfind(t, ' adaptation');
    if ~isempty(k) && ~startsWith(t, 'No ')
        t = {t(1:k(1) - 1), t(k(1) + 1:end)};
    end
end
end

function close_valid(h)
h = h(isgraphics(h));
if ~isempty(h); close(h); end
end

function retitle(axs, names, two_line)
if nargin < 3; two_line = false; end
st = manuscript_style();
names = cellstr(names);
assert(numel(axs) == numel(names), 'composite_tools:Titles', ...
    '%d axes to title with %d condition names.', numel(axs), numel(names));
for k = 1:numel(axs)
    color = [0 0 0];
    if isKey(st.condition_color, names{k}); color = st.condition_color(names{k}); end
    title(axs(k), condition_title(names{k}, two_line), 'FontWeight', 'normal', 'FontSize', 14, ...
        'Color', color, 'Interpreter', 'tex');
end
end

function s = sentence_case(s)
% The first word capitalized when the text opens with two lowercase letters
% (a word; 'n_a' and 'lambda' are symbols and stay).
was_char = ischar(s);
s = cellstr(s);
for k = 1:numel(s)
    if ~isempty(regexp(s{k}, '^[a-z][a-z]', 'once')); s{k}(1) = upper(s{k}(1)); end
end
if was_char; s = s{1}; end
end

function h = crop_margin(fig, px)
% An almost-white line at mid-width, px above the canvas bottom (see the header).
if nargin < 2; px = 2; end
p = getpixelposition(fig);
y = px / p(4);
h = annotation(fig, 'line', 0.5 + [0, 3 / p(3)], [y y], 'Color', [0.99 0.99 0.99], 'LineWidth', 0.5, ...
    'Tag', 'composite_crop_margin');
end
