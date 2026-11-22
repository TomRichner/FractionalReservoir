function h = draw_offscale_markers(ax, x, medians, cond_spec, y_window, marker_size, step_share)
% DRAW_OFFSCALE_MARKERS Small triangles on the top and bottom edge of an axes where a curve has left the fixed y window.
%
%   h = DRAW_OFFSCALE_MARKERS(ax, x, medians, cond_spec, y_window, marker_size, step_share)
%
%   ax           the axes (hold on), its y limits y_window
%   x            1 x n sweep values
%   medians      1 x C cell, medians{c} the n values of condition c's curve
%   cond_spec    1 x C struct with fields name and color
%   y_window     [lo hi], the window the axes shows
%   marker_size  points
%   step_share   share of the x range between two conditions' markers at
%                one sweep value and edge
%   h            the markers: lines tagged 'offscale_marker', UserData
%                .condition (name), .x (the sweep value), .side (1 top,
%                2 bottom)
%
% A triangle pointing up on the top edge at every x(k) where a condition's
% median lies above y_window, pointing down on the bottom edge where it lies
% below, in the condition's color. Conditions off scale at the same x and
% edge are set side by side, centered on x(k), and a group at an end of the
% sweep moves inward as one so that it stays inside the x limits. The
% markers sit just inside the edge (2.5% of the window), unclipped, so one
% is whole at the first and last sweep value and none reaches a tick label.
% They are kept out of the legend. NaN medians are not marked. The axis
% scale is not touched.
%
% Used by fig_sensitivity_medians (its 'offscale_markers', default true),
% and through it by main Figure 4, panel (A), and the supplementary sheet.
%
% See also: fig_sensitivity_medians, mark_default_value
h = gobjects(0);
xl = [min(x), max(x)];
dx = step_share * diff(xl);
inset = 0.025 * diff(y_window);
edge_y = [y_window(2) - inset, y_window(1) + inset];
marks = {'^', 'v'};
for k = 1:numel(x)
    for side = 1:2
        if side == 1
            off = find(cellfun(@(y) y(k) > y_window(2), medians));
        else
            off = find(cellfun(@(y) y(k) < y_window(1), medians));
        end
        n_off = numel(off);
        if n_off == 0; continue; end
        xs = x(k) + dx * ((1:n_off) - (n_off + 1) / 2);
        xs = xs - max(0, max(xs) - xl(2)) + max(0, xl(1) - min(xs));   % a group at an end moves inward as one
        for j = 1:n_off
            c = off(j);
            h(end + 1) = plot(ax, xs(j), edge_y(side), marks{side}, 'Color', cond_spec(c).color, ...
                'MarkerFaceColor', cond_spec(c).color, 'MarkerSize', marker_size, 'LineStyle', 'none', ...
                'Clipping', 'off', 'HandleVisibility', 'off', 'Tag', 'offscale_marker', ...
                'UserData', struct('condition', cond_spec(c).name, 'x', x(k), 'side', side)); %#ok<AGROW>
        end
    end
end
end
