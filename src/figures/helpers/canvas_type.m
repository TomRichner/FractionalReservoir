function A = canvas_type(fig, exceptions)
% CANVAS_TYPE Set every piece of type on a composite figure to its role's size on the common canvas.
%
%   A = CANVAS_TYPE(fig)
%   A = CANVAS_TYPE(fig, exceptions)
%
% The roles and sizes are canvas_audit's (tick labels, legends and free
% text tick_fs; axis and colorbar labels label_fs; titles title_fs; panel
% letters letter_fs; manuscript_style). Items matching an exception
% ({pattern, reason}, as in canvas_audit) keep the size they have. An
% axes' labels and title are given their size after the axes', so the
% axes' FontSize does not carry over into them (LabelFontSizeMultiplier
% and TitleFontSizeMultiplier are set to 1 as well).
%
% Called by the composites of group A (Plan 2, canvas standard) once their
% panels and letters are in place, before composite_save. Returns the audit
% after the change (every row ok).
%
% See also: canvas_audit, composite_figure, manuscript_style

arguments
    fig (1,1) matlab.ui.Figure
    exceptions (:,2) cell = cell(0, 2)
end
A = canvas_audit(fig, exceptions);
fix = ~A.ok | ~cellfun(@isempty, A.reason);
fix = fix & cellfun(@isempty, A.reason);            % exceptions keep their size
is_ax = strcmp(A.kind, 'axes');
for k = find(fix & is_ax)'
    set(A.h(k), 'LabelFontSizeMultiplier', 1, 'TitleFontSizeMultiplier', 1);
    % Setting an axes' FontSize resets its labels' and title's; keep them
    % and set them again below.
    keep = {A.h(k).XLabel.FontSize, A.h(k).YLabel.FontSize, A.h(k).Title.FontSize};
    A.h(k).FontSize = A.want(k);
    A.h(k).XLabel.FontSize = keep{1}; A.h(k).YLabel.FontSize = keep{2}; A.h(k).Title.FontSize = keep{3};
end
for k = find(fix & ~is_ax)'
    A.h(k).FontSize = A.want(k);
end
A = canvas_audit(fig, exceptions);
end
