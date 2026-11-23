function [fig, B] = composite_figure(rows, opts)
% COMPOSITE_FIGURE One figure from the panels of several source figures, stacked in rows, at their own size.
%
%   [fig, B] = COMPOSITE_FIGURE(rows)
%   [fig, B] = COMPOSITE_FIGURE(rows, 'visible', false, 'row_gap', 0.3)
%
% The composite main and supplementary figures of the spectral analyses
% (fig_main5_mechanism, fig_main6_amplification, fig_supp6_grouped to
% fig_supp8_grouped) are built from figures that already exist: each source
% function is called with save = false, and its axes and legends are
% copied here (copyobj, as fig_grouped_main does), so a panel is drawn by
% one piece of code wherever it appears.
%
%   rows   cell array, one element per row, top to bottom. Each element is
%          a struct array of BLOCKS, left to right, with fields
%            fig     the source figure (closed here once it is copied)
%            width   inches the block is given; [] or absent: the source's
%                    own width. The block is the source scaled uniformly to
%                    that width, so 1:1 panels stay 1:1. Text keeps its
%                    point size (only positions scale).
%            strip   (optional) inches kept free above this row for its
%                    letter, in place of opts.strip; [] or absent: opts.strip.
%                    0 for a row whose source keeps a title strip of its own
%                    where the letter goes (S7 panel C). The row takes the
%                    largest of its blocks' values.
%
% Options
%   visible   false
%   margin    [left right top bottom] inches around everything ([0.15 0.15 0.1 0.15])
%   row_gap   inches between rows (0.25)
%   col_gap   inches between the blocks of a row (0.35)
%   strip     inches kept free above every row for its panel letters (0.32)
%   align     'left' (default) or 'center': rows narrower than the canvas
%             centered in it
%   width     inches: the figure is this wide, whatever its rows need (the
%             common canvas, manuscript_style canvas_w; a warning when the
%             rows need more). [] (default): as wide as the widest row.
%             The caller sizes the blocks (their width field) to fill it.
%
% Returns the figure (units inches) and B, a cell like rows, of structs
%   rect      [x y w h] of the block on the figure, inches
%   ax        the copied axes, sorted top row first, left to right
%   legends   the copied legends
%   top       the y of the block's top edge (the letter strip is above it)
%
% Letters are added by the caller with COMPOSITE_LETTER, which puts "(A)" in
% the strip above a block, at the block's left edge or at the left edge of
% one of its axes.
%
% A composite is usually taller than the screen; save_figure_stable exports
% such a figure whole (save_some_figs_to_folder_2).
%
% See also: composite_letter, composite_tables, fig_grouped_main

arguments
    rows (1,:) cell
    opts.visible (1,1) logical = false
    opts.margin (1,4) double = [0.15 0.15 0.1 0.15]
    opts.row_gap (1,1) double = 0.25
    opts.col_gap (1,1) double = 0.35
    opts.strip (1,1) double = 0.32
    opts.align (1,:) char {mustBeMember(opts.align, {'left', 'center'})} = 'left'
    opts.width double = []
end

ppi = get(groot, 'ScreenPixelsPerInch');
n_rows = numel(rows);
row_w = zeros(1, n_rows); row_h = zeros(1, n_rows); row_strip = opts.strip * ones(1, n_rows);
S = cell(1, n_rows);
for r = 1:n_rows
    blocks = rows{r};
    if isfield(blocks, 'strip')
        v = [blocks.strip];
        if ~isempty(v); row_strip(r) = max(v); end
    end
    S{r} = struct('scale', {}, 'w', {}, 'h', {});
    for k = 1:numel(blocks)
        drawnow;
        px = getpixelposition(blocks(k).fig);
        w0 = px(3) / ppi; h0 = px(4) / ppi;
        s = 1;
        if isfield(blocks, 'width') && ~isempty(blocks(k).width); s = blocks(k).width / w0; end
        S{r}(k) = struct('scale', s, 'w', w0 * s, 'h', h0 * s);
    end
    row_w(r) = sum([S{r}.w]) + opts.col_gap * (numel(blocks) - 1);
    row_h(r) = max([S{r}.h]);
end
m = opts.margin;
fig_w = m(1) + max(row_w) + m(2);
if ~isempty(opts.width)
    if fig_w > opts.width + 1e-6
        warning('composite_figure:TooWide', 'The rows need %.2f in; the canvas is %.2f in wide.', fig_w, opts.width);
    end
    fig_w = opts.width;
end
span = fig_w - m(1) - m(2);                          % the width the rows are aligned in
fig_h = m(3) + sum(row_h) + sum(row_strip) + (n_rows - 1) * opts.row_gap + m(4);
vis = 'off'; if opts.visible; vis = 'on'; end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 fig_w fig_h], 'Visible', vis);
% The size the layout is made for. A window larger than the screen can be
% shrunk by MATLAB at a moment of its choosing (it happened once while the
% letters were being placed, 2026-11-21); composite_letter and
% composite_save put this size back before they measure or save.
setappdata(fig, 'composite_size', [fig_w fig_h]);

B = cell(1, n_rows);
y_top = fig_h - m(3);
for r = 1:n_rows
    blocks = rows{r};
    y_top = y_top - row_strip(r);                    % the letters' strip
    x = m(1);
    if strcmp(opts.align, 'center'); x = m(1) + (span - row_w(r)) / 2; end
    for k = 1:numel(blocks)
        s = S{r}(k);
        rect = [x, y_top - s.h, s.w, s.h];
        [ax, lg] = copy_block(blocks(k).fig, fig, rect, s.scale, ppi);
        B{r}(k) = struct('rect', rect, 'ax', ax, 'legends', lg, 'top', y_top);
        if isgraphics(blocks(k).fig); close(blocks(k).fig); end
        x = x + s.w + opts.col_gap;
    end
    y_top = y_top - row_h(r) - opts.row_gap;
end
end

%% ------------------------------------------------------------------------
function [ax, lg] = copy_block(src, dst, rect, scale, ppi)
% Axes, legends and colorbars of src into rect of dst (inches), together so
% that a legend keeps its axes.
drawnow;
h = findall(src, '-isa', 'matlab.graphics.axis.Axes', '-or', 'Type', 'legend', '-or', 'Type', 'colorbar');
ax = gobjects(0); lg = gobjects(0);
if isempty(h); return; end
pos = zeros(numel(h), 4);
for k = 1:numel(h); pos(k, :) = getpixelposition(h(k), true) / ppi; end
new = copyobj(h, dst);
for k = 1:numel(new)
    new(k).Units = 'inches';
    if isa(new(k), 'matlab.graphics.axis.Axes'); new(k).PositionConstraint = 'innerposition'; end
    new(k).Position = [rect(1:2) + pos(k, 1:2) * scale, pos(k, 3:4) * scale];
    if isa(new(k), 'matlab.graphics.axis.Axes')
        % copyobj returns the title and the axis names to the axes' own
        % font size (seen 2026-11-21: 14-point condition titles came out
        % at the hidden axes' 10 points). Give them back what the source set.
        for name = {'Title', 'XLabel', 'YLabel'}
            set(new(k).(name{1}), 'FontSize', h(k).(name{1}).FontSize, 'Color', h(k).(name{1}).Color, ...
                'FontWeight', h(k).(name{1}).FontWeight, 'Visible', h(k).(name{1}).Visible);
        end
    end
end
is_ax = arrayfun(@(o) isa(o, 'matlab.graphics.axis.Axes'), new);
ax = new(is_ax); lg = new(strcmp(get(new, {'Type'}), 'legend'));
p = cell2mat(arrayfun(@(a) a.Position, ax(:), 'UniformOutput', false));
% Top row first, then left to right (rows told apart by a tenth of an inch).
[~, order] = sortrows([-round(10 * (p(:, 2) + p(:, 4))), p(:, 1)]);
ax = ax(order);
end
