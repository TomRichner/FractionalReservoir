function files = stim_fig8_compose(r, B, out_dir, tag, clinical_note, model_note, visible)
%STIM_FIG8_COMPOSE Figure 8 on the common canvas: the patients, the modular model, then the bursting example.
%
%   files = stim_fig8_compose(r, B, out_dir, tag, clinical_note, model_note, visible)
%
% r     fig_stim_combined's output, drawn with 'save' false (r.fig open,
%       the mirror layout: clinical traces, rate and PSD; model traces,
%       rate and PSD)
% B     bm_load_runs' struct (the burst_mode stage; B.R{1} is drawn)
%
% THE ORDER (Tom, 2026-11-22, Plan 2 notes on Figure 8, note 1: "I do want
% the real data first, the human data first, followed by model data"):
%
%   A | B/C  the patients: SEEG traces, no stimulation beside stimulation
%            (left); the IED rate (B) over the SOZ spectrum (C) at the right
%   D | E/F  the same three for the modular model
%   G        the bursting example: the time course and, at its right, the
%            zoom on one burst (bm_fig8_panels), with one key above
%   H        its five complex planes
%
% No panels I or J. Until 2026-11-22 the bursting example came first (A to
% C) and the patients and the model followed (D to I).
%
% THE CANVAS: manuscript_style's canvas_w (9 in) wide and at most
% canvas_max_h (11.8 in) tall, type by role (canvas_type). Laid out for
% 9 in, in inches, bottom up. The rate and spectrum panels sit in a narrow
% right column (the rate panels have no legends now, notes 9 and 14), so
% the traces have most of the width.
%
% The clinical and modular axes are fig_stim_combined's, copied and placed
% here; nothing in the clinical repository is changed. The patient panel
% (cscs_stim_stability, 2026-11-22) comes with its titles, 5 s of signal,
% marks, circles and scale bar as Tom accepted them: none of these is
% restyled here beyond the type sizes by role. fig_stim_combined thinned
% lines and marks for its 7-in canvas (times 0.6, scatter 0.35); that is
% undone here. Titles and tick labels are put in sentence case without
% counts (restyle_titles); the provenance file lists the originals. Each
% spectrum (C, F) gets a labelled tick at every decade it spans, at least
% two (psd_decades; Tom saw a single labelled tick).
%
% Called by ied_stim_fig_combined; the output holds clinical recordings
% (never git add -f; images stay under figs/).
%
% See also: ied_stim_fig_combined, fig_stim_combined, bm_fig8_panels, canvas_type

st = manuscript_style();
old = r.fig;
guard = onCleanup(@() close_if(old));

%% Geometry (inches), bottom up
% The right column sets the height of the upper four rows: rate over
% spectrum, the rate's condition names on one line under it, the patient
% spectrum's tick labels and x label (each spectrum names its frequency
% axis since the Plan 3 review). G and H take what remains of canvas_max_h.
W = st.canvas_w;
hR = 1.0;          % rate and spectrum axes
gR = 0.45;         % between them: the rate panel's condition names, one line
xlabC = 0.58;      % under C: tick labels and the x label (Plan 3 review: C had none)
xlabF = 0.58;      % under F: tick labels and the x label
wR = 1.6; xR = W - 0.22 - wR;                    % right column: B, C, E, F (room right for the "100" under C and F)
xN = 0.95;                                       % left edge of the traces (channel names left of it)
gT = 0.2; wT = (xR - 1.5 - xN - gT) / 2;         % two trace strips; the stimulation strip's bars and the
xT = [xN, xN + wT + gT];                         % rate panels' y labels right of them
rest = 0.1 + 0.15 + xlabF + 2 * hR + gR + 0.35 + xlabC + 2 * hR + gR + 0.42;
hBand = st.canvas_max_h - 0.05 - rest;           % G and H
band = [0.1, 0.1, W - 0.2, hBand];
yF = band(2) + hBand + 0.15 + xlabF;  yE = yF + hR + gR;  topM = yE + hR;
yD = yF + 0.1;
yC = topM + 0.35 + xlabC;             yB = yC + hR + gR;  topP = yB + hR;
yA = yC + 0.1;
Ht = topP + 0.42;
assert(Ht <= st.canvas_max_h + 1e-9, 'stim_fig8_compose:Tall', 'Figure 8 would be %.2f in tall (canvas_max_h %.2f).', Ht, st.canvas_max_h);

vis = 'off'; if visible; vis = 'on'; end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [0.2 0.2 W Ht], 'Visible', vis);
setappdata(fig, 'composite_size', [W Ht]);

%% G and H
P = bm_fig8_panels(fig, band, B.R{1}, 'font_size', st.tick_fs, 'annot_size', st.tick_fs, 'line_scale', 1.8, ...
    'left', 1.05, 'c_title', 0.35, 'ab_room', 0.62, 'plane_gap', 0.25, 'key_room', 0.34);

%% A to F: fig_stim_combined's axes, placed
ax = findall(old, 'Type', 'axes');
set(ax, 'Units', 'normalized');
p = cell2mat(get(ax, 'Position'));
[~, order] = sortrows([-round(20 * (p(:, 2) + p(:, 4))), p(:, 1)]);   % top row first, then left to right
ax = ax(order);
assert(numel(ax) == 8, 'stim_fig8_compose:Axes', 'Expected the 8 axes of the mirror layout, found %d.', numel(ax));
lgs = findall(old, 'Type', 'legend');
place = {ax(1), [xT(1), yA, wT, topP - yA]; ax(2), [xT(2), yA, wT, topP - yA]; ...       % A
         ax(3), [xR, yB, wR, hR]; ax(4), [xR, yC, wR, hR]; ...                            % B, C
         ax(5), [xT(1), yD, wT, topM - yD]; ax(6), [xT(2), yD, wT, topM - yD]; ...       % D
         ax(7), [xR, yE, wR, hR]; ax(8), [xR, yF, wR, hR]};                               % E, F
new_ax = gobjects(1, 8); new_lg = cell(1, 8);
for k = 1:8
    a = place{k, 1};
    lg = lgs(arrayfun(@(l) isequal(l.Axes, a), lgs));
    if isempty(lg) || ismember(k, [3 7])         % B and E: no legends (note 9)
        na = copyobj(a, fig);
    else
        h = copyobj([lg(1) a], fig); na = h(2); new_lg{k} = h(1);
    end
    na.Units = 'inches'; na.PositionConstraint = 'innerposition'; na.Position = place{k, 2};
    new_ax(k) = na;
end
unthin(new_ax);
changes = restyle_titles(new_ax, [new_lg{:}]);

%% Type by role, then the rate axes, the spectra and the letters
canvas_type(fig, fig8_type_exceptions());
set([P.ax_A P.ax_B], 'LineWidth', st.axis_lw);   % bm_fig8_panels scales its axis lines with line_scale
drawnow;
for k = [3 7]                                    % B and E: rates as % of no stimulation
    a = new_ax(k);
    % Three labelled ticks: at 14 pt the nine labels 0 to 200 % ran into each other.
    keep = ismember(a.YTickLabel, {'0', '10%', '100%'});
    if any(keep); set(a, 'YTick', a.YTick(keep), 'YTickLabel', a.YTickLabel(keep)); end
    % The condition names short, on one line: at 1.6 in "No stimulation"
    % and "Stimulation" ran into each other even on two lines.
    a.XTickLabel = {'No stim.', 'Stim.'};
    a.XTickLabelRotation = 0;
end
% The y names on two lines no longer than the 1 in axes, so that the
% names of B and C (E and F) do not run into each other: "% of no
% stimulation" is said by the gray 100 % line and by the caption.
new_ax(3).YLabel.String = {'IED rate', '(%)'};
new_ax(7).YLabel.String = {'Event rate', '(%)'};
new_ax(4).YLabel.String = {'SOZ PSD', '(\muV^2/Hz)'};
new_ax(8).YLabel.String = {'PSD', '(a.u.^2/Hz)'};
for k = [4 8]                                    % C and F
    a = new_ax(k);
    psd_decades(a);
    % Decades only on the frequency axis, upright: five or six labels were
    % turned by MATLAB on a 1.6 in axis.
    xd = 10 .^ (ceil(log10(a.XLim(1)) - 1e-9):floor(log10(a.XLim(2)) + 1e-9));
    set(a, 'XTick', xd, 'XTickLabel', arrayfun(@(v) sprintf('%g', v), xd, 'UniformOutput', false), ...
        'XTickLabelRotation', 0, 'XMinorTick', 'off');
    % No key: on a 1.6 x 1 in axes it lay on the spectra wherever it was
    % put, and the colors are those of the traces' titles (purple no
    % stimulation, green stimulation).
    if ~isempty(new_lg{k}); delete(new_lg{k}); new_lg{k} = []; end
end
new_ax(4).XLabel.String = new_ax(8).XLabel.String;   % C: its own frequency name (Plan 3 review), as F
xl = xR - 1.45;                                  % the letters of the right column, left of the two-line y labels
letter(fig, 'A', 0.05, new_ax(1), 0.08); letter(fig, 'B', xl, new_ax(3)); letter(fig, 'C', xl, new_ax(4));
letter(fig, 'D', 0.05, new_ax(5), 0.08); letter(fig, 'E', xl, new_ax(7)); letter(fig, 'F', xl, new_ax(8));
letter(fig, 'G', 0.05, P.ax_A(1), 0.06);
letter(fig, 'H', 0.05, P.ax_C(1), 0.3);

%% Save
if ~isfolder(out_dir); mkdir(out_dir); end
fig.Position(3:4) = [W Ht];
drawnow; true_minus_ticks(fig);   % the true minus on every numeric tick label, last before export
save_figure_stable(out_dir, tag, fig, canvas_export_args());   % the whole canvas, untrimmed
prov = fullfile(out_dir, [tag '_provenance.md']);
write_provenance(prov, r, B, clinical_note, model_note, W, Ht, changes);
tbl = fullfile(out_dir, [tag '_table.md']);
bm_write_table(tbl, B, P.info, 'Figure 8 G and H');
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

function psd_decades(a)
% A labelled tick at every decade the spectrum's axis spans, and at least
% two: when it spans fewer, the limits are widened outward to whole
% decades. More than four decades: a tick at every other one, so that the
% labels do not crowd a short axis.
lo = floor(log10(a.YLim(1)) + 1e-9); hi = ceil(log10(a.YLim(2)) - 1e-9);
e = ceil(log10(a.YLim(1)) - 1e-9):floor(log10(a.YLim(2)) + 1e-9);
if numel(e) < 2
    a.YLim = 10 .^ [lo hi];
    e = lo:hi;
end
if numel(e) > 4; e = e(mod(e - e(end), 2) == 0); end
a.YTick = 10 .^ e;
a.YTickLabel = arrayfun(@(k) sprintf('10^{%d}', k), e, 'UniformOutput', false);
a.YMinorTick = 'off';
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
fprintf(fid, ['# Figure 8 provenance\n\nWritten by `scripts/paper/stim_fig8_compose.m` (through ' ...
    '`ied_stim_fig_combined`) on %s. Panels A to F are the axes of `fig_stim_combined` (mirror layout), copied, ' ...
    'placed on the common canvas and given its type by role; G and H are drawn from the burst_mode stage by ' ...
    '`bm_fig8_panels`.\n\nContains clinical recordings: do not add to git (`git add -f`) or publish outside ' ...
    'the manuscript.\n\n## Sources\n\n'], char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
src = r.sources;
fprintf(fid, '- Bursting network: `%s` (first network, seed %s)\n', B.index_file, mat2str(B.R{1}.seed));
fprintf(fid, '- Clinical traces: `%s`\n- Clinical rate: `%s`\n- Clinical PSD: `%s`\n- Model: `%s`\n', ...
    src.traces, src.rate, src.psd, src.model);
fprintf(fid, '- Clinical note: %s\n- Model note: %s\n', clinical_note, model_note);
fprintf(fid, '\n## Panel letters\n\n| Letter | Shows | Source |\n|---|---|---|\n');
m = r.letters;
for k = 1:size(m, 1)
    fprintf(fid, '| (%s) | %s | %s |\n', m{k, 1}, m{k, 2}, m{k, 3});
end
fprintf(fid, ['| (G) | bursting network: input; mean E and I rates; real part of the rightmost mode, gains frozen ' ...
    'and full Jacobian (1/s); after the step, dashed, the same network without the step (no-stimulation ' ...
    'continuation, same noise); at the right, one burst: the mode crosses the boundary before the rates take ' ...
    'off | burst_mode, bm_fig8_panels |\n' ...
    '| (H) | eigenvalues of M as rates: structure, quiet, before burst, burst, stimulation | same |\n']);
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
% case, no counts and no parentheses (polish list, 2026-11-22). A title
% that only states a statistic and a count ("Median, IQR over subjects
% (n = 12)") leaves the panel; the provenance file keeps every original.
% Returns the changes, n x 2 {was, is}.
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
for lg = lgs(:)'
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
