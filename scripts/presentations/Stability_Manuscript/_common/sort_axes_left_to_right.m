function ax_sorted = sort_axes_left_to_right(fig)
% SORT_AXES_LEFT_TO_RIGHT A figure's axes ordered by x-position.
%
%   ax = SORT_AXES_LEFT_TO_RIGHT(fig)
%
% findobj returns axes in creation/stacking order, which is not the order they
% appear on screen. Every figure that copies panels out of a regenerated plot
% needs them left-to-right, because that is the order the conditions were drawn
% in and therefore the order the column titles assume.
%
% Was a local subfunction in two figure scripts.

ax = findobj(fig, 'Type', 'axes');
if isempty(ax)
    ax_sorted = ax;
    return
end
p = cell2mat(get(ax, 'Position'));
if size(p, 1) == 1
    % get() on a single handle returns a 1x4 row, not a cell; cell2mat above
    % still yields 1x4, so this is only a guard against sorting a scalar.
    ax_sorted = ax;
    return
end
[~, order] = sort(p(:, 1));
ax_sorted = ax(order);
end
