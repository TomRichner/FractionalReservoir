function extend_canvas_below(fig, ax_all, margin)
% EXTEND_CANVAS_BELOW Grow a figure downward until its x labels are on the canvas.
%
%   EXTEND_CANVAS_BELOW(fig, ax_all)          % margin 0.01 of the figure height
%   EXTEND_CANVAS_BELOW(fig, ax_all, margin)
%
% Finds the lowest x-label bottom over ax_all (figure-normalized). When it is
% below `margin`, the figure is made taller by the overrun, added at the
% bottom, and every axes keeps its size and place in PIXELS (their normalized
% Positions are rescaled). So a figure whose labels already fit is left
% exactly as it was, and one that did not gains only the strip it needs.
% Annotations drawn before the call are not moved, so call it before drawing
% any. Written 2026-11-12 (FD-006) for the parameter-space sheets, whose
% 300-px canvas left the 20-pt x label ('lambda_1 (1/s)') printed cut in half,
% on the chained bundle and on the product reference alike. Assumes normalized
% axes Positions.
%
% A label that ends flush with the PNG's edge is not necessarily cut:
% exportgraphics crops to the content, and a label without descenders
% ('I -> I') then touches the edge and looks clipped in a downscaled preview.
% Check at full resolution before reaching for this.
%
% See also: fig_param_space_allStd
if nargin < 3; margin = 0.01; end
low = inf;
for a = 1:numel(ax_all)
    ax = ax_all(a);
    xl = get(ax, 'XLabel');
    if isempty(get(xl, 'String')); continue; end
    p = get(ax, 'Position');
    old_units = get(xl, 'Units');
    set(xl, 'Units', 'normalized');
    e = get(xl, 'Extent');
    set(xl, 'Units', old_units);
    low = min(low, p(2) + e(2) * p(4));
end
if ~isfinite(low) || low >= margin
    return
end
old_units = get(fig, 'Units');
set(fig, 'Units', 'pixels');
fp = get(fig, 'Position');
H  = fp(4);
d  = (margin - low) * H;                 % pixels to add below
H2 = H + d;
set(fig, 'Position', [fp(1), fp(2) - d, fp(3), H2]);
for a = 1:numel(ax_all)
    p = get(ax_all(a), 'Position');
    set(ax_all(a), 'Position', [p(1), (p(2) * H + d) / H2, p(3), p(4) * H / H2]);
end
set(fig, 'Units', old_units);
end
