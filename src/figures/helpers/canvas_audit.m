function A = canvas_audit(fig, exceptions)
% CANVAS_AUDIT Every piece of type on a composite figure, by role, against the common canvas of manuscript_style.
%
%   A = CANVAS_AUDIT(fig)
%   A = CANVAS_AUDIT(fig, exceptions)
%
% The standard (the manuscript repository's
% docs/Advancing_the_narrative/Plan2_canvas_standard.md, Tom 2026-11-22):
% every main and supplementary composite is canvas_w wide and at most
% canvas_max_h tall, and its type has one size per role:
%
%   role      what                                           size
%   tick      axes (their tick labels), legends, colorbar     tick_fs
%             ticks, free text in axes, annotations
%   label     axis names (also those complex_plane_axes       label_fs
%             writes as text, 'Re (s^{-1})'), colorbar labels
%   title     axes titles and subtitles                       title_fs
%   letter    panel letters, "(A)"                             letter_fs
%
% exceptions  n x 2 cell {pattern, reason}: an item whose name matches the
%             regular expression pattern is allowed another size; the
%             reason is kept in A (column reason) so that a test names it.
%             The name of a text in an axes is '<text> | axes <tag> "<title>"'.
%
% A is a table, one row per visible item with text: kind (axes, xlabel,
% ylabel, title, subtitle, legend, colorbar, cblabel, text, annotation),
% role, name, size (points), want (points), ok (size == want, or an
% exception), reason, h (the graphics object). Its UserData holds the
% figure's size in inches (fig_in) and the thinnest visible axes line and
% data line (min_axis_lw, min_line_lw).
%
% Nothing is changed in the figure; canvas_type sets the sizes.
%
% See also: canvas_type, manuscript_style, test_canvas_standard_a

arguments
    fig (1,1) matlab.ui.Figure
    exceptions (:,2) cell = cell(0, 2)
end
st = manuscript_style();
want = struct('tick', st.tick_fs, 'label', st.label_fs, 'title', st.title_fs, 'letter', st.letter_fs);
rows = cell(0, 5);                                   % kind, role, name, size, handle

for ax = findall(fig, '-isa', 'matlab.graphics.axis.Axes')'
    if strcmp(ax.Tag, 'print_whole_check'); continue; end
    ax_name = sprintf('axes %s "%s"', ax.Tag, flat(ax.Title.String));
    if strcmp(ax.Visible, 'on') && has_ticks(ax)
        rows(end + 1, :) = {'axes', 'tick', ax_name, ax.FontSize, ax}; %#ok<AGROW>
    end
    for nm = {'XLabel', 'YLabel'}
        h = ax.(nm{1});
        if shows(h) && (strcmp(ax.Visible, 'on') || strcmp(h.Visible, 'on'))
            rows(end + 1, :) = {lower(nm{1}), 'label', [flat(h.String) ' | ' ax_name], h.FontSize, h}; %#ok<AGROW>
        end
    end
    if shows(ax.Title)
        rows(end + 1, :) = {'title', 'title', flat(ax.Title.String), ax.Title.FontSize, ax.Title}; %#ok<AGROW>
    end
    if isprop(ax, 'Subtitle') && shows(ax.Subtitle)
        rows(end + 1, :) = {'subtitle', 'title', flat(ax.Subtitle.String), ax.Subtitle.FontSize, ax.Subtitle}; %#ok<AGROW>
    end
    for t = findall(ax, 'Type', 'text')'
        if is_label_of(t, ax) || ~shows(t); continue; end
        s = flat(t.String);
        role = 'tick';
        if is_letter(t, s)
            role = 'letter';
        elseif strcmp(t.Tag, 'complex_plane_axes') && ~isempty(regexp(s, '^(Re|Im)\>', 'once'))
            role = 'label';                          % complex_plane_axes writes the axis names as text
        end
        rows(end + 1, :) = {'text', role, [s ' | ' ax_name], t.FontSize, t}; %#ok<AGROW>
    end
end
for lg = findall(fig, 'Type', 'legend')'
    if strcmp(lg.Visible, 'off'); continue; end
    rows(end + 1, :) = {'legend', 'tick', strjoin(cellstr(lg.String), ', '), lg.FontSize, lg}; %#ok<AGROW>
end
for cb = findall(fig, 'Type', 'colorbar')'
    if strcmp(cb.Visible, 'off'); continue; end
    rows(end + 1, :) = {'colorbar', 'tick', 'colorbar', cb.FontSize, cb}; %#ok<AGROW>
    if shows(cb.Label)
        rows(end + 1, :) = {'cblabel', 'label', flat(cb.Label.String), cb.Label.FontSize, cb.Label}; %#ok<AGROW>
    end
end
for an = findall(fig, 'Type', 'textboxshape')'
    if ~shows(an); continue; end
    s = flat(an.String);
    role = 'tick';
    if is_letter(an, s); role = 'letter'; end
    rows(end + 1, :) = {'annotation', role, s, an.FontSize, an}; %#ok<AGROW>
end

n = size(rows, 1);
kind = rows(:, 1); role = rows(:, 2); name = rows(:, 3);
sz = zeros(n, 1); w = zeros(n, 1);
for k = 1:n; sz(k) = rows{k, 4}; w(k) = want.(role{k}); end
h = gobjects(n, 1);
for k = 1:n; h(k) = rows{k, 5}; end
ok = abs(sz - w) < 0.01;
reason = repmat({''}, n, 1);
for k = 1:n
    for e = 1:size(exceptions, 1)
        if ~isempty(regexp(name{k}, exceptions{e, 1}, 'once'))
            reason{k} = exceptions{e, 2};
            ok(k) = true;
            break
        end
    end
end
A = table(kind, role, name, sz, w, ok, reason, h, 'VariableNames', {'kind', 'role', 'name', 'size', 'want', 'ok', 'reason', 'h'});

u = fig.Units; fig.Units = 'inches'; fsz = fig.Position(3:4); fig.Units = u;
vis_ax = findall(fig, '-isa', 'matlab.graphics.axis.Axes', 'Visible', 'on');
vis_ax = vis_ax(~strcmp(get(vis_ax, {'Tag'}), 'print_whole_check'));
alw = arrayfun(@(a) a.LineWidth, vis_ax);
lns = findall(fig, 'Type', 'line', 'Visible', 'on');
llw = arrayfun(@(l) l.LineWidth, lns);
A.Properties.UserData = struct('fig_in', fsz, 'min_axis_lw', min([alw(:); Inf]), 'min_line_lw', min([llw(:); Inf]));
end

%% ------------------------------------------------------------------------
function tf = has_ticks(ax)
tf = (~isempty(ax.XTick) && ~isempty(ax.XTickLabel) && strcmp(ax.XAxis(1).Visible, 'on')) || ...
     (~isempty(ax.YTick) && ~isempty(ax.YTickLabel) && strcmp(ax.YAxis(1).Visible, 'on'));
end

function tf = shows(h)
tf = strcmp(h.Visible, 'on') && ~isempty(strtrim(flat(h.String)));
end

function tf = is_label_of(t, ax)
tf = any(t == [ax.XLabel ax.YLabel ax.ZLabel ax.Title]) || (isprop(ax, 'Subtitle') && t == ax.Subtitle);
end

function tf = is_letter(h, s)
tf = startsWith(h.Tag, 'panel_letter') || any(strcmp(h.Tag, {'letter', 'composite_panel_letter'})) || ...
    ~isempty(regexp(s, '^\([A-Z]\)$', 'once'));
end

function s = flat(x)
if isempty(x); s = ''; return; end
if iscell(x) || isstring(x); x = strjoin(cellstr(x), ' '); end
if size(x, 1) > 1; x = strjoin(cellstr(x), ' '); end
s = strtrim(char(x));
end
