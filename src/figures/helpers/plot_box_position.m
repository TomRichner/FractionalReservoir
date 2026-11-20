function pos = plot_box_position(ax)
% PLOT_BOX_POSITION The plot box of an axes, in normalized units of its figure.
%
%   pos = PLOT_BOX_POSITION(ax)
%
% With a constrained aspect ratio (axis image, axis equal) the plot box is
% smaller than the axes' Position: it is the largest box of that aspect
% ratio centered in the inner position. This returns [x y w h] of the box
% that is actually drawn (tightPosition), for placing something (an inset)
% relative to it. It works for an axes in a tiled layout too, whose
% Position the layout manages.
%
% See also: fig_eig_heatmap, fig_main3_stability, tightPosition
drawnow;
pos = tightPosition(ax);
if ~strcmp(ax.Units, 'normalized')
    fig = ancestor(ax, 'figure');
    px = hgconvertunits(fig, pos, ax.Units, 'pixels', ax.Parent);
    fig_px = getpixelposition(fig);
    pos = px ./ [fig_px(3:4), fig_px(3:4)];
end
end
