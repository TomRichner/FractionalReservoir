function files = stim_fig8_compose(r, B, out_dir, tag, clinical_note, model_note, visible)
%STIM_FIG8_COMPOSE Figure 8 in the new order, on the common canvas: the bursting network's panels above the clinical and modular panels.
%
%   files = stim_fig8_compose(r, B, out_dir, tag, clinical_note, model_note, visible)
%
% r     fig_stim_combined's output, drawn with 'save' false (r.fig open,
%       the mirror layout: clinical traces, rate and PSD; model traces,
%       rate and PSD)
% B     bm_load_runs' struct (the burst_mode stage; B.R{1} is drawn)
%
% THE COMMON CANVAS (Plan 2, canvas standard, Tom 2026-11-22: "the figures
% should have similar canvas widths so that the letters appear
% proportionally the same"). The figure is manuscript_style's canvas_w
% (13 in) wide and at most canvas_max_h tall, its type by role (canvas_type:
% tick labels, legends and text tick_fs, axis names label_fs, titles
% title_fs, letters letter_fs). Until then it was 7 in wide with 9-pt type
% and 12-pt letters, its panels stacked in the full width one under the
% other, and the canvas shrunk to the screen with every size scaled by as
% much. Now the same nine panels in the same order, laid out again:
%
%   A, B, C  the bursting network (bm_fig8_panels, drawn into the top band
%            with the room for 14-pt type): A the whole run, B one burst
%            beside it, C the five complex planes under them
%   D | E/F  the clinical example traces, no stimulation beside
%            stimulation (left), the IED rate (E) over the SOZ PSD (F) at
%            the right
%   G | H/I  the same for the modular model
%
% The clinical and modular axes are fig_stim_combined's, copied with their
% legends and placed in inches here; nothing in the clinical repository is
% changed. fig_stim_combined thinned their lines and marks for its 7-in
% canvas (times 0.6, scatter 0.35); here that is undone, so that no hairline
% stands beside 14-pt type. The titles, tick labels and legend entries are
% put in sentence case and lose their counts (restyle_titles); the
% provenance file lists the originals.
%
% A figure taller than the screen is exported whole by save_figure_stable
% (save_some_figs_to_folder_2 prints it at its own size), as every composite
% is; the canvas is no longer scaled.
%
% Called by ied_stim_fig_combined; the output holds clinical recordings
% (never git add -f).
%
% See also: ied_stim_fig_combined, fig_stim_combined, bm_fig8_panels, canvas_type

st = manuscript_style();
old = r.fig;
guard = onCleanup(@() close_if(old));

%% Geometry (inches), bottom up
W = st.canvas_w;
xT = [1.0, 0.304 * W]; wT = 0.204 * W;                    % trace columns of D and G (names left of the first)
xR = 0.673 * W; wR = W - 0.15 - xR;                  % right column: E, F, H, I
wLg = 1.75;                                      % E's and H's keys, right of their axes
hR = 1.6;                                       % E, F, H, I axes height
yI = 0.75;                 yH = yI + hR + 0.85;
topG = yH + hR;            yG = 1.05;            % G's traces down to its scale bar
yF = topG + 0.35 + 0.75;   yE = yF + hR + 0.85;
topD = yE + hR;            yD = yF + 0.3;
band = [0.1, topD + 0.45, W - 0.2, 6.25];        % A to C
Ht = band(2) + band(4) + 0.3;
assert(Ht <= st.canvas_max_h, 'stim_fig8_compose:Tall', 'Figure 8 would be %.2f in tall (canvas_max_h %.2f).', Ht, st.canvas_max_h);

vis = 'off'; if visible; vis = 'on'; end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [0.2 0.2 W Ht], 'Visible', vis);
setappdata(fig, 'composite_size', [W Ht]);

%% A to C
P = bm_fig8_panels(fig, band, B.R{1}, 'font_size', st.tick_fs, 'annot_size', st.tick_fs, 'line_scale', 1.8, ...
    'left', 1.05, 'c_title', 0.35, 'ab_room', 1.05, 'plane_gap', 0.3);

%% D to I: fig_stim_combined's axes, placed
ax = findall(old, 'Type', 'axes');
set(ax, 'Units', 'normalized');
p = cell2mat(get(ax, 'Position'));
[~, order] = sortrows([-round(20 * (p(:, 2) + p(:, 4))), p(:, 1)]);   % top row first, then left to right
ax = ax(order);
assert(numel(ax) == 8, 'stim_fig8_compose:Axes', 'Expected the 8 axes of the mirror layout, found %d.', numel(ax));
lgs = findall(old, 'Type', 'legend');
place = {ax(1), [xT(1), yD, wT, topD - yD]; ax(2), [xT(2), yD, wT, topD - yD]; ...       % D
         ax(3), [xR, yE, wR - wLg, hR]; ax(4), [xR, yF, wR, hR]; ...                     % E, F
         ax(5), [xT(1), yG, wT, topG - yG]; ax(6), [xT(2), yG, wT, topG - yG]; ...       % G
         ax(7), [xR, yH, wR - wLg, hR]; ax(8), [xR, yI, wR, hR]};                        % H, I
new_ax = gobjects(1, 8); new_lg = cell(1, 8);
for k = 1:8
    a = place{k, 1};
    lg = lgs(arrayfun(@(l) isequal(l.Axes, a), lgs));
    if isempty(lg)
        na = copyobj(a, fig);
    else
        h = copyobj([lg(1) a], fig); na = h(2); new_lg{k} = h(1);
    end
    na.Units = 'inches'; na.PositionConstraint = 'innerposition'; na.Position = place{k, 2};
    new_ax(k) = na;
end
unthin(new_ax);
changes = restyle_titles(new_ax, [new_lg{:}]);

%% Type by role, then the keys, the rate axes and the letters
canvas_type(fig, fig8_type_exceptions());
drawnow;
for k = [3 7]                                    % E and H: rates as % of no stimulation
    a = new_ax(k);
    % Every other tick: at 14 pt the nine labels 0 to 200 % ran into each other.
    keep = ismember(a.YTickLabel, {'0', '10%', '100%'});
    if any(keep); set(a, 'YTick', a.YTick(keep), 'YTickLabel', a.YTickLabel(keep)); end
    % The condition names under the ticks, the long one on two lines: on one
    % line MATLAB turned them by 30 degrees into the panel below, and a tick
    % label cannot hold two lines (a newline makes it two labels), so they
    % are text in the axes, at tick-label size.
    names = a.XTickLabel; xt = a.XTick;
    a.XTickLabel = {};
    for j = 1:numel(xt)
        words = strsplit(names{j}, ' ');
        if numel(words) > 1; s = {strjoin(words(1:end - 1), ' '), words{end}}; else; s = names{j}; end
        text(a, (xt(j) - a.XLim(1)) / diff(a.XLim), -0.05, s, 'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontSize', manuscript_style().tick_fs, 'Clipping', 'off', 'Tag', 'fig8_condition_label');
    end
end
for k = [2 6]                                    % D and G: the amplitude bar's label right of the bar
    a = new_ax(k);
    for t = findall(a, 'Type', 'text')'
        if strcmp(t.HorizontalAlignment, 'right') && ~isempty(regexp(char(t.String), '^\d.*V$', 'once'))
            t.Position(1) = a.XLim(2) + 0.015 * diff(a.XLim);
            t.HorizontalAlignment = 'left'; t.Clipping = 'off';
        end
    end
end
for k = [3 7]                                    % E's and H's keys right of their axes
    lg = new_lg{k};
    if isempty(lg); continue; end
    lg.Units = 'inches'; lg.Location = 'none';
    q = new_ax(k).Position;
    lg.Position(1:2) = [q(1) + q(3) + 0.15, q(2) + (q(4) - lg.Position(4)) / 2];
end
for k = [4 8]                                    % the PSD keys inside, bottom left
    if ~isempty(new_lg{k}); new_lg{k}.Location = 'southwest'; end
end
lx = 0.05;                                       % letters of the left column at the canvas edge
lx_r = xR - 1.3;                                % of the right column, left of the y labels
letter(fig, 'A', lx, P.ax_A(1)); letter(fig, 'B', P.ax_B(1).Position(1) - 0.75, P.ax_B(1));
letter(fig, 'C', lx, P.ax_C(1), 0.33);
letter(fig, 'D', lx, new_ax(1), 0.3); letter(fig, 'E', lx_r, new_ax(3)); letter(fig, 'F', lx_r, new_ax(4));
letter(fig, 'G', lx, new_ax(5), 0.3); letter(fig, 'H', lx_r, new_ax(7)); letter(fig, 'I', lx_r, new_ax(8));

%% Save
if ~isfolder(out_dir); mkdir(out_dir); end
fig.Position(3:4) = [W Ht];
save_figure_stable(out_dir, tag, fig);
prov = fullfile(out_dir, [tag '_provenance.md']);
write_provenance(prov, r, B, clinical_note, model_note, W, Ht, changes);
tbl = fullfile(out_dir, [tag '_table.md']);
bm_write_table(tbl, B, P.info, 'Figure 8 A to C');
files = [fullfile(out_dir, strcat(tag, {'.png', '.fig', '.svg'})), {prov, tbl}];
if ~visible; close(fig); end
end

%% ------------------------------------------------------------------------
function ex = fig8_type_exceptions()
% The type that may differ from its role's size (canvas_audit), with why.
% None: the channel and row names, scale bars and event marks of the
% traces are tick-sized text like every other.
ex = cell(0, 2);
end

function letter(fig, ch, x_in, ax, raise)
% The panel letter, its bottom at the top of ax plus raise (inches; the
% height of a title above the axes), its left at x_in.
if nargin < 5; raise = 0.05; end
st = manuscript_style();
u = ax.Units; ax.Units = 'inches'; p = ax.Position; ax.Units = u;
annotation(fig, 'textbox', 'Units', 'inches', 'Position', [x_in, p(2) + p(4) + raise, 0.5, 0.3], ...
    'String', ['(' ch ')'], 'FontSize', st.letter_fs, 'FontName', 'Helvetica', 'LineStyle', 'none', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Margin', 0, 'Tag', 'letter');
end

function unthin(axs)
% fig_stim_combined drew its lines at 0.6 times their width (at least 0.5)
% and its scatter marks at 0.35 times their area, for a 7-in canvas; on the
% common canvas they get their width back.
for h = findall(axs, 'Type', 'line')'
    h.LineWidth = h.LineWidth / 0.6; h.MarkerSize = h.MarkerSize / 0.6;
end
for h = findall(axs, 'Type', 'constantline')'
    h.LineWidth = h.LineWidth / 0.6;
end
for h = findall(axs, 'Type', 'scatter')'
    h.SizeData = h.SizeData / 0.35;
end
st = manuscript_style();
set(axs, 'LineWidth', st.axis_lw);
end

function close_if(f)
if ~isempty(f) && isvalid(f); close(f); end
end

function write_provenance(file, r, B, clinical_note, model_note, W, Ht, changes)
fid = fopen(file, 'w'); guard = onCleanup(@() fclose(fid));
fprintf(fid, ['# Figure 8 (new order) provenance\n\nWritten by `scripts/paper/stim_fig8_compose.m` (through ' ...
    '`ied_stim_fig_combined`) on %s. Panels D to I are the axes of `fig_stim_combined` (mirror layout), copied, ' ...
    'placed on the common canvas and given its type by role; A to C are drawn from the burst_mode stage by ' ...
    '`bm_fig8_panels`.\n\nContains clinical recordings: do not add to git (`git add -f`) or publish outside ' ...
    'the manuscript.\n\n## Sources\n\n'], char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
src = r.sources;
fprintf(fid, '- Bursting network: `%s` (first network, seed %s)\n', B.index_file, mat2str(B.R{1}.seed));
fprintf(fid, '- Clinical traces: `%s`\n- Clinical rate: `%s`\n- Clinical PSD: `%s`\n- Model: `%s`\n', ...
    src.traces, src.rate, src.psd, src.model);
fprintf(fid, '- Clinical note: %s\n- Model note: %s\n', clinical_note, model_note);
fprintf(fid, '\n## Panel letters\n\n| Letter | Shows | Source |\n|---|---|---|\n');
fprintf(fid, ['| (A) | bursting network: input; mean E and I rates, adaptation frozen dashed; real part of the ' ...
    'rightmost mode, gains frozen and full Jacobian (1/s) | burst_mode, bm_fig8_panels |\n' ...
    '| (B) | one burst: the mode crosses the boundary before the rates take off | same |\n' ...
    '| (C) | eigenvalues of M as rates: structure, quiet, before burst, burst, stimulation | same |\n']);
m = r.letters; new = 'DEFGHI';
for k = 1:size(m, 1)
    fprintf(fid, '| (%s) | %s | %s (was %s) |\n', new(k), m{k, 2}, m{k, 3}, m{k, 1});
end
if ~isempty(changes)
    fprintf(fid, ['\n## Titles restyled on the composite\n\nThe saved panel figures carry these titles, tick labels ' ...
        'and legend entries; the composite shows them in sentence case, without counts (an empty new title: the ' ...
        'title left the panel, and what it said belongs in the caption).\n\n| In the panel figure | On Figure 8 |\n|---|---|\n']);
    for k = 1:size(changes, 1)
        fprintf(fid, '| %s | %s |\n', changes{k, 1}, changes{k, 2});
    end
end
fprintf(fid, ['\n## Canvas\n\n%.2f x %.2f in, the common canvas (manuscript_style canvas_w), exported by ' ...
    'save_figure_stable.\n'], W, Ht);
end

function changes = restyle_titles(axs, lgs)
% The titles, the condition tick labels and the legend entries of the
% clinical and modular panels, as the paper's style has them: sentence
% case, no counts and no parentheses (polish list, 2026-11-22). They are
% written in the saved panel figures that fig_stim_combined loads (the
% clinical ones in another repository), so they are restyled here, on the
% composite. A title that only states a statistic and a count ("Median, IQR
% over subjects (n = 12)") leaves the panel; the provenance file keeps every
% original. Returns the changes, n x 2 {was, is}.
changes = cell(0, 2);
for ax = axs(:)'
    was = ax.Title.String;
    if ischar(was) && ~isempty(was)
        is = restyled(was);
        if ~strcmp(is, was); ax.Title.String = is; changes(end + 1, :) = {was, is}; end %#ok<AGROW>
    end
    lab = ax.XTickLabel;
    if iscell(lab) && ~isempty(lab) && all(cellfun(@(s) ~isempty(regexp(s, '^[A-Za-z]', 'once')), lab))
        new = cellfun(@restyled, lab, 'UniformOutput', false);
        if ~isequal(new, lab)
            ax.XTickLabel = new;
            changes = [changes; [lab(:), new(:)]]; %#ok<AGROW>
        end
    end
end
for lg = lgs(:)'                              % D to I only: "Full Jacobian" of A keeps its capital
    lab = lg.String;
    new = cellfun(@restyled, lab, 'UniformOutput', false);
    if ~isequal(new, lab)
        lg.String = new;
        changes = [changes; [lab(:), new(:)]]; %#ok<AGROW>
    end
end
if isempty(changes); return; end
[~, k] = unique(strcat(changes(:, 1), '|', changes(:, 2)), 'stable');
changes = changes(k, :);
end

function s = restyled(s)
if startsWith(s, 'Median, IQR'); s = ''; return; end
s = strtrim(regexprep(s, '\s*\([^)]*\)', ''));          % "(n = 12)", "(Uniform Input 0.2)"
words = strsplit(s, ' ');
for k = 2:numel(words)                                  % all but the first word; acronyms keep their capitals
    w = words{k};
    if ~isempty(regexp(w, '^[A-Z][a-z]+[,.]?$', 'once')); words{k} = [lower(w(1)) w(2:end)]; end
end
s = strjoin(words, ' ');
end
