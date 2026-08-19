function mark_default_value(ax, xd, color, lw, frac)
% MARK_DEFAULT_VALUE Short vertical tick rising off the x-axis at xd.
%
%   MARK_DEFAULT_VALUE(ax, xd, color, lw, frac)
%
% Marks the preset's default for this panel -- the network the sweep departs
% from. Spans ylim(1) to ylim(1) + frac*diff(ylim), i.e. it rises from the
% bottom axis (these axes are YDir normal, so ylim(1) IS the bottom edge).
%
% Skipped when the default lies outside the swept range: a marker clamped to
% the edge would assert the default sits at the end of the sweep when it
% actually sits beyond it.
%
% Extracted verbatim from the local subfunction of
% Fig_sensitivity_analysis_allStd.m so fig_sensitivity_medians can share it.
%
% See also: preset_default_values, apply_percent_axis
    xl = xlim(ax);
    if xd < min(xl) || xd > max(xl)
        return;
    end
    yl = ylim(ax);
    y_top = yl(1) + frac * diff(yl);

    was_held = ishold(ax);
    hold(ax, 'on');
    plot(ax, [xd, xd], [yl(1), y_top], '-', ...
        'Color', color, 'LineWidth', lw, ...
        'Clipping', 'on', 'HandleVisibility', 'off');
    if ~was_held
        hold(ax, 'off');
    end

    % Adding a line can widen the limits, which would shift every other element
    % on the panel; pin them back to what the data set.
    xlim(ax, xl);
    ylim(ax, yl);
end
