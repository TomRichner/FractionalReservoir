function files = stim_fig8_compose(r, B, out_dir, tag, clinical_note, model_note, visible)
%STIM_FIG8_COMPOSE Figure 8 in the new order: the bursting network's panels above the present figure's clinical and modular panels.
%
%   files = stim_fig8_compose(r, B, out_dir, tag, clinical_note, model_note, visible)
%
% r     fig_stim_combined's output, drawn with 'save' false (r.fig open,
%       the mirror layout: clinical A to C, modular model D to F)
% B     bm_load_runs' struct (the burst_mode stage; B.R{1} is drawn)
%
% The printed figure is 7 in wide (the present figure's width) and 4.75 in
% taller: the bursting panels (bm_fig8_panels, drawn on a canvas of their
% own) on top, every axes of r.fig (with its legend) below, both copied by
% normalized position. The present figure's letters A to F become D to I;
% A to C are new. Font sizes print as before, so the lower part looks
% exactly like the present figure.
%
% THE CANVAS IS SCALED TO THE SCREEN. 12.4 in is taller than the
% workstation's 1073-px screen, and a canvas taller than the screen is
% clipped by exportgraphics (fig_stim_combined warns of it). So the
% canvas is made s times the printed size (s = 0.84 there), every font,
% line and marker s times smaller, and the PNG written at 300 / s dpi:
% the same pixels as the printed size at 300 dpi (2100 px wide). The .fig
% and .svg are at the scaled size.
%
% Called by ied_stim_fig_combined; the output holds clinical recordings
% (never git add -f).
%
% See also: ied_stim_fig_combined, fig_stim_combined, bm_fig8_panels

old = r.fig;
guard = onCleanup(@() close_if(old));
old.Units = 'inches';
W = old.Position(3); Hm = old.Position(4);
Htop = 4.75;
Ht = Hm + Htop;
top = figure('Color', 'w', 'Units', 'inches', 'Position', [0.2 0.2 W Htop], 'Visible', 'off');
guard_top = onCleanup(@() close_if(top));
P = bm_fig8_panels(top, [0, 0.02, W, Htop - 0.3], B.R{1}, 'font_size', 9, 'annot_size', 8);
fs_letter = 12;
old_letters = findall(old, 'Tag', 'letter');
if ~isempty(old_letters); fs_letter = old_letters(1).FontSize; end
letter_in(top, P.ax_A(1), '(A)', 0.035, fs_letter, 0.03);
letter_in(top, P.ax_B(1), '(B)', P.ax_B(1).Position(1) - 0.5, fs_letter, 0.03);
letter_in(top, P.ax_C(1), '(C)', 0.035, fs_letter, 0.2);

mon = get(0, 'MonitorPositions'); ppi = get(0, 'ScreenPixelsPerInch');
s = min(1, (max(mon(:, 4)) / ppi - 0.8) / Ht);
vis = 'off'; if visible; vis = 'on'; end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [0.2 0.2 W * s Ht * s], 'Visible', vis);
move_into(top, fig, Hm / Ht, Htop / Ht, containers.Map({'(A)', '(B)', '(C)'}, {'(A)', '(B)', '(C)'}));
move_into(old, fig, 0, Hm / Ht, containers.Map({'(A)', '(B)', '(C)', '(D)', '(E)', '(F)'}, ...
    {'(D)', '(E)', '(F)', '(G)', '(H)', '(I)'}));
scale_marks(fig, s);
if ~isfolder(out_dir); mkdir(out_dir); end
for ext = {'.png', '.fig', '.svg'}
    ff = fullfile(out_dir, [tag ext{1}]); if isfile(ff); delete(ff); end
end
fig.PaperPositionMode = 'auto';
% 300 dpi at the printed size. Higher was tried and is wrong: at 450 and 538
% dpi on this canvas exportgraphics drew the y tick labels of the top panels
% displaced from their ticks (2026-11-21; at 359 dpi on the s-scaled canvas,
% that is 300 dpi printed, they are right).
exportgraphics(fig, fullfile(out_dir, [tag '.png']), 'Resolution', round(300 / s), 'Padding', 10);
savefig(fig, fullfile(out_dir, [tag '.fig']));
print(fig, fullfile(out_dir, [tag '.svg']), '-dsvg', '-vector');
prov = fullfile(out_dir, [tag '_provenance.md']);
write_provenance(prov, r, B, clinical_note, model_note, W, Ht);
tbl = fullfile(out_dir, [tag '_table.md']);
bm_write_table(tbl, B, P.info, 'Figure 8 A to C');
files = [fullfile(out_dir, strcat(tag, {'.png', '.fig', '.svg'})), {prov, tbl}];
if ~visible; close(fig); end
end

%% ------------------------------------------------------------------------
function letter_in(fig, ax, s, x_in, fs, raise)
u = ax.Units; ax.Units = 'inches'; p = ax.Position; ax.Units = u;
annotation(fig, 'textbox', 'Units', 'inches', 'Position', [x_in, p(2) + p(4) + raise, 0.4, 0.2], ...
    'String', s, 'FontSize', fs, 'FontName', 'Helvetica', 'LineStyle', 'none', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Margin', 0, 'Tag', 'letter');
end

function move_into(src, fig, y0, hf, letters)
% Copy every axes (with its legend) and every panel letter of src into
% fig, by normalized position, into the band [y0, y0 + hf] of fig's height;
% letters maps src's letter strings to fig's.
place = @(p) [p(1), y0 + p(2) * hf, p(3), p(4) * hf];
axs = findall(src, 'Type', 'axes');
lgs = findall(src, 'Type', 'legend');
for a = axs'
    a.Units = 'normalized';
    lg = lgs(arrayfun(@(l) isequal(l.Axes, a), lgs));
    if isempty(lg)
        na = copyobj(a, fig);
    else
        lg(1).Units = 'normalized'; lp = lg(1).Position; manual = strcmp(lg(1).Location, 'none');
        h = copyobj([lg(1) a], fig); na = h(2); nl = h(1);
        if manual
            nl.Units = 'normalized'; nl.Position = place(lp);
        end
    end
    na.Units = 'normalized';
    na.Position = place(a.Position);
end
for t = findall(src, 'Tag', 'letter')'
    t.Units = 'normalized';
    annotation(fig, 'textbox', place(t.Position), 'String', letters(char(t.String)), 'FontSize', t.FontSize, ...
        'FontName', t.FontName, 'LineStyle', 'none', 'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'left', 'Margin', 0, 'Tag', 'letter');
end
end

function scale_marks(fig, s)
% Every font size, line width and marker size times s, so that the canvas
% made s times smaller prints as the full-size one (the PNG is written at
% 450 / s dpi). Every size is read first and then set: an axes' labels
% would otherwise follow the axes and be scaled twice, depending on order.
if s >= 1; return; end
hf = findall(fig, '-property', 'FontSize');
sz = arrayfun(@(h) h.FontSize, hf);
for k = 1:numel(hf)
    hf(k).FontSize = sz(k) * s;
end
for h = findall(fig, 'Type', 'line')'
    h.LineWidth = h.LineWidth * s; h.MarkerSize = h.MarkerSize * s;
end
for h = findall(fig, 'Type', 'constantline')'
    h.LineWidth = h.LineWidth * s;
end
for h = findall(fig, 'Type', 'scatter')'
    h.SizeData = h.SizeData * s^2; h.LineWidth = h.LineWidth * s;
end
for h = findall(fig, 'Type', 'axes')'
    h.LineWidth = h.LineWidth * s;
end
drawnow;
end


function close_if(f)
if ~isempty(f) && isvalid(f); close(f); end
end

function write_provenance(file, r, B, clinical_note, model_note, W, Ht)
fid = fopen(file, 'w'); guard = onCleanup(@() fclose(fid));
fprintf(fid, ['# Figure 8 (new order) provenance\n\nWritten by `scripts/paper/stim_fig8_compose.m` (through ' ...
    '`ied_stim_fig_combined`) on %s. Panels D to I are the axes of the present figure (`fig_stim_combined`, ' ...
    'mirror layout), copied and moved down unchanged in size and font; A to C are drawn from the burst_mode stage ' ...
    'by `bm_fig8_panels`.\n\nContains clinical recordings: do not add to git (`git add -f`) or publish outside ' ...
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
fprintf(fid, '\n## Canvas\n\n%.2f x %.2f in (the printed size), PNG by print at 450 dpi.\n', W, Ht);
end
