function rgb = condition_rgb_by_column(ax, pos, fallback, st)
% CONDITION_RGB_BY_COLUMN The condition colour of each axes of a sheet whose columns are conditions.
%
%   rgb = CONDITION_RGB_BY_COLUMN(ax, pos, fallback)
%   rgb = CONDITION_RGB_BY_COLUMN(ax, pos, fallback, st)
%
%   ax        axes handles of one sheet (any order)
%   pos       their positions, numel(ax) x 4 (normalized; [left bottom w h])
%   fallback  1 x 3 colour for a column whose condition is not known
%   st        manuscript_style() (optional; read when absent)
%   rgb       numel(ax) x 3, the colour of each axes' column
%
% The sheets assembled from ParamSpaceAnalysis2's plotters (the sensitivity
% sheets, the joint-sample histograms) title each column with its condition
% title (srnn_condition_titles) and draw every condition's data in one shared
% colour. This reads the condition of each COLUMN from the title of any axes
% in it (grouped_condition_label 'name_of_title'), so the data can be drawn in
% the condition's manuscript_style colour (Plan 2, package A, 2026-11-21: no
% condition in a shared, black or gray colour). Columns are grouped by left
% edge. A column with no title, or a title no condition carries, gets the
% fallback, so an unknown condition degrades to the old look, not an error.
%
% See also: grouped_condition_label, manuscript_style, srnn_condition_titles
if nargin < 4 || isempty(st); st = manuscript_style(); end
n = numel(ax);
rgb = repmat(reshape(fallback, 1, 3), n, 1);
if n == 0; return; end
[~, ~, col_of] = uniquetol(pos(:, 1), 0.01);
for c = reshape(unique(col_of), 1, [])
    members = find(col_of == c);
    name = '';
    for m = reshape(members, 1, [])
        t = get(get(ax(m), 'Title'), 'String');
        if iscell(t); t = strjoin(t, ' '); end
        if isempty(strtrim(char(t))); continue; end
        try
            hit = grouped_condition_label({char(t)}, 'name_of_title');
            name = hit{1};
            break
        catch
        end
    end
    if ~isempty(name) && isKey(st.condition_color, name)
        rgb(members, :) = repmat(st.condition_color(name), numel(members), 1);
    end
end
end
