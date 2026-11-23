function out = fig_supp9_grouped(cfg)
% FIG_SUPP9_GROUPED Supplementary Figure S9, memory details; saved data only.
%
%   out = FIG_SUPP9_GROUPED('run_dir', d, 'preset_name', p, 'out_dir', o)
%   out = FIG_SUPP9_GROUPED('run_dir', d, 'mc_run_dir', g, ...)   % the grid from another run
%
% What backs main Figure 7 and ties responsiveness to memory. Every panel
% is drawn by the function that owns it, invisibly and unsaved, and its
% axes are copied here (composite_tools); nothing is simulated. Memory only
% since Tom's review (Plan2_figure_feedback.md, S9; 2026-11-22), on the
% 9 in canvas:
%   (A) example reconstructions at four delays, two by two: the four trace
%       axes of fig_memory_capacity_example (its two summary axes are left
%       out; Figure 7A and 7B carry them for the ensemble).
%   (B) the responsive share of neurons over the gain by self-coupling
%       grid, one map per condition, in the layout of the capacity maps of
%       Figure 7E: fig_mc_gain_bias_maps with rows {'responsive'} (median
%       over seeds of 1 - frac_pinned - frac_silent, the quantity
%       'frac_on_slope' of mc_sweep_quantity). Self-coupling marked at its
%       ends and at 0, as on Figure 7.
%   (C) memory capacity against lambda_1 for every run of that grid:
%       fig_mc_vs_lambda1 with the bundle's registry arguments. Its
%       statistics boxes are not drawn; their text is in the table.
% The two panels that left S9: the shares of the leading direction against
% the longest SFA time constant are in fig_supp_leading_grouped, with the
% two other panels on what the leading direction is made of; the cumulative
% memory capacity is Figure 7B.
%
% WHICH RUN: as fig_main7_memory. run_dir holds mc_example and
% memory_capacity; the gain by self-coupling grid is read from mc_run_dir,
% by default the same run.
%
% Options: run_dir, mc_run_dir, preset_name (accepted for a uniform call),
% run_mode, figure_set (the bundle whose registry entries style the maps
% and the scatter), conditions (empty = the figure set's), out_dir, save,
% visible, fig_tag.
%
% out: figs, files, source, panels (letter -> source function).
%
% See also: composite_tools, fig_main7_memory, fig_mc_gain_bias_maps,
%           fig_mc_vs_lambda1, fig_memory_capacity_example, fig_supp_leading_grouped

arguments
    cfg.run_dir     (1,:) char    = ''
    cfg.mc_run_dir  (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''          % accepted for a uniform call
    cfg.run_mode    (1,:) char    = 'fastmedium'
    cfg.figure_set  (1,:) char    = '3cond'
    cfg.conditions                = {}
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = false
    cfg.verbose     (1,:) char    = 'minimal'   % accepted for a uniform call
    cfg.fig_tag     (1,:) char    = 'Fig_Supp9_Grouped'

end
setup_paths();
T = composite_tools();
st = manuscript_style();
if isempty(cfg.mc_run_dir); cfg.mc_run_dir = cfg.run_dir; end
bundle = mu7chain_gainbias_config(cfg.run_mode, cfg.figure_set);
if isempty(cfg.conditions); cfg.conditions = bundle.conditions; end
names = cellstr(cfg.conditions);
nc = numel(names);
quiet = {'save', false, 'visible', false};
map_entry = registry_entry(bundle.figures, 'fig_mc_maps_selfcoupling');
sct_entry = registry_entry(bundle.figures, 'fig_mc_vs_lambda1');

%% The source figures
S_ex  = fig_memory_capacity_example('run_dir', cfg.run_dir, 'conditions', names, quiet{:});
S_map = fig_mc_gain_bias_maps('run_dir', cfg.mc_run_dir, 'conditions', names, quiet{:}, ...
    map_entry.args{:}, 'rows', {'responsive'});
S_sct = fig_mc_vs_lambda1('run_dir', cfg.mc_run_dir, 'conditions', names, quiet{:}, sct_entry.args{:});
src_figs = [S_ex.figs(:); S_map.figs(:); S_sct.figs(:)];
cleanup = onCleanup(@() T.close(src_figs));

%% The canvas, laid out in pixels on canvas_w itself (fit has nothing to stretch)
W = round(st.canvas_w * get(groot, 'ScreenPixelsPerInch'));
% Rows from the bottom: the marker-size key, (C) the scatter panels with
% their two-line titles, (B) the maps with theirs, (A) the traces two by
% two with the scale bar under the last.
KEY_Y = 6; KEY_H = 70;
Y3 = 140; H3 = 170;                          % (C)
Y2 = Y3 + H3 + 48 + 72; MH = 180;            % (B): above the titles of (C) and the x labels of (B)
TR_H = 100; TR_GAP = 34;                     % (A): one trace (66 until the traces were no longer cut flat:
                                             % their range is now the readout's, up to about 4 times the input's), the title below
Y1 = Y2 + MH + 48 + 46;                      % bottom row of (A): above the titles of (B) and the scale bar
H = Y1 + 2 * TR_H + TR_GAP + 52;            % the titles of the top traces and the letter above them
fig = T.canvas(W, H, cfg.visible);

% (A): the trace axes of the example (titled 'delay = ...').
ex_axes = findall(S_ex.figs, 'Type', 'axes');
is_trace = arrayfun(@(a) startsWith(strjoin(cellstr(a.Title.String), ' '), 'delay'), ex_axes);
assert(nnz(is_trace) == 4, 'fig_supp9_grouped:NoTraces', ...
    'fig_memory_capacity_example drew %d axes titled ''delay = ...''; S9 lays out four.', nnz(is_trace));
tx = findall(ex_axes(is_trace), 'Type', 'text');            % the source's own panel letter
delete(tx(arrayfun(@(t) ~isempty(regexp(strjoin(cellstr(t.String), ' '), '^\([a-z]\)$', 'once')), tx)));
A = T.copy(ex_axes(is_trace), fig, false);
[~, order] = sort(arrayfun(@(a) -a.Position(2), A.axes)); ax_a = A.axes(order);   % top first: the delays in order
TR_X = [46 466]; TR_W = 350;                % the traces' hidden y tick labels still count in their inset
for k = 1:4
    row = ceil(k / 2); col = 2 - mod(k, 2);
    set(ax_a(k), 'Units', 'pixels', 'Position', [TR_X(col), Y1 + (2 - row) * (TR_H + TR_GAP), TR_W, TR_H]);
end

% (B): the responsive-share maps and their colour bar.
Bm = T.copy(S_map.figs, fig);
assert(numel(Bm.axes) == nc && isscalar(Bm.colorbars), 'fig_supp9_grouped:Maps', ...
    'Expected %d responsive-share maps and one colour bar.', nc);
[~, order] = sort(arrayfun(@(a) a.Position(1), Bm.axes)); maps = Bm.axes(order);
X0 = 82; MGAP = 14; CB_ROOM = 100;           % the maps' geometry of Figure 7
MW = (W - X0 - 2 * MGAP - 10 - 14 - CB_ROOM) / nc;
for c = 1:nc
    set(maps(c), 'Units', 'pixels', 'PositionConstraint', 'innerposition', ...
        'Position', [X0 + (c - 1) * (MW + MGAP), Y2, MW, MH]);
end
x_cb = X0 + nc * MW + (nc - 1) * MGAP + 10;
set(Bm.colorbars, 'Units', 'pixels', 'Position', [x_cb, Y2, 14, MH]);

% (C): the scatter panels, their colour bar and the marker-size key.
stats_text = scatter_statistics(S_sct.figs);
delete(findall(S_sct.figs, 'Tag', 'mc_sweep_scatter_stats'));
Cs = T.copy(S_sct.figs, fig);
is_key = arrayfun(@(a) strcmp(a.Tag, 'mc_sweep_scatter_size_legend'), Cs.axes);
sct = Cs.axes(~is_key);
[~, order] = sort(arrayfun(@(a) a.Position(1), sct)); sct = sct(order);
SX0 = 72; SGAP = 36; SCB_ROOM = 104;          % room for the colour bar's ticks and its 'Gain' (at 86 it touched the edge)
SW = (W - SX0 - (nc - 1) * SGAP - 12 - 14 - SCB_ROOM) / nc;
for c = 1:numel(sct)
    set(sct(c), 'Units', 'pixels', 'PositionConstraint', 'innerposition', ...
        'Position', [SX0 + (c - 1) * (SW + SGAP), Y3, SW, H3]);
end
if ~isempty(Cs.colorbars)
    set(Cs.colorbars, 'Units', 'pixels', 'Position', [SX0 + nc * SW + (nc - 1) * SGAP + 12, Y3, 14, H3]);
end

%% Style
T.style(fig);
% (A)
set(ax_a, 'TitleFontSizeMultiplier', 1);
% A key for the two traces (Plan 3 review: nothing said which was the input
% and which the reconstruction), under the bottom-left traces, off the data,
% beside the scale bar of the bottom-right ones.
hold(ax_a(3), 'on');
tr = findall(ax_a(3), 'Type', 'line');
in_col = [0 0 0];
rec = tr(arrayfun(@(l) ~isequal(l.Color, in_col) && numel(l.XData) > 2, tr));
hk = [plot(ax_a(3), NaN, NaN, '-', 'Color', in_col, 'LineWidth', 1.5), ...
    plot(ax_a(3), NaN, NaN, '-', 'Color', rec(1).Color, 'LineWidth', 1.5)];
lgA = legend(ax_a(3), hk, {'Input', 'Reconstruction'}, 'Orientation', 'horizontal', 'Box', 'off', ...
    'AutoUpdate', 'off', 'FontSize', st.tick_fs);
lgA.ItemTokenSize = [24 10];
% The scale bar's label as Figure 8's ("1 a.u."), whatever the source wrote.
au = findall(ax_a, 'Type', 'text');
set(au(arrayfun(@(t) any(strcmp(strtrim(cellstr(t.String)), '1 AU')), au)), 'String', '1 a.u.');
% (B): ticks and labels as Figure 7 styles the capacity maps.
GAIN_TICKS = [0.05 0.25 1 5];
for c = 1:nc
    ax = maps(c);
    title(ax, T.condition_title(names{c}, true), 'FontWeight', 'normal', 'FontSize', st.title_fs, ...
        'Color', st.condition_color(names{c}), 'Interpreter', 'tex');
    g_levels = ax.XTick;
    set(ax, 'XTick', log10(GAIN_TICKS), 'XTickLabel', compose('%g', GAIN_TICKS), ...
        'XTickLabelRotation', 0, 'XMinorTick', 'on');
    ax.XAxis.MinorTickValues = g_levels;
    s_levels = ax.YTick;
    % Self-coupling marked at its ends and at 0 (as Figure 7E and 7F; S9
    % note 4 with Figure 7 note 6), every level a minor tick.
    ticks = unique([s_levels([1 end]), 0]);
    set(ax, 'YTick', ticks(ticks >= ax.YLim(1) & ticks <= ax.YLim(2)), 'YMinorTick', 'on');
    ax.YAxis.MinorTickValues = s_levels;
    xlabel(ax, 'Gain');
    if c == 1
        ax.YTickLabel = compose('%g', ax.YTick);
        ylabel(ax, 'Self-Coupling, s');
    else
        ax.YTickLabel = {};
        ylabel(ax, '');
    end
end
Bm.colorbars.Label.String = 'Responsive Share';
% (C): two-line titles (on one line they ran into each other at 9 in).
for c = 1:numel(sct)
    title(sct(c), T.condition_title(names{c}, true), 'FontWeight', 'normal', 'Color', sct(c).Title.Color, ...
        'Interpreter', 'tex');
    xlabel(sct(c), 'Lyapunov Exp., \lambda_1 (s^{-1})', 'Interpreter', 'tex');
end
ylabel(sct(1), 'Memory Capacity');
set(sct, 'Box', 'off');                  % open axes, as every other panel (the source boxes them)
if ~isempty(Cs.colorbars); Cs.colorbars.Label.String = 'Gain'; end
clip_note = '';
if any(is_key)
    clip_note = respace_size_key(Cs.axes(is_key), W, KEY_Y, KEY_H);
end

%% The common canvas (the layout is on canvas_w already: fit moves nothing)
T.style(fig);                % again: the labels and titles set since take their roles' sizes
T.fit(fig);
% (A)'s key under the bottom-left traces, level with the scale bar's time label.
set(lgA, 'Units', 'pixels', 'Location', 'none');
lgA.Position(1:2) = [TR_X(1), Y1 - lgA.Position(4) - 6];

%% Panel letters
T.letter(ax_a(1), 'A', -TR_X(1), 22);
T.letter(maps(1), 'B', -X0 + 2, 50);
T.letter(sct(1), 'C', -SX0 + 2, 50);

%% Save
panels = struct( ...
    'A', 'fig_memory_capacity_example (its four trace axes): example reconstructions at several delays', ...
    'B', 'fig_mc_gain_bias_maps with rows {''responsive''}: responsive share of neurons over gain by self-coupling', ...
    'C', 'fig_mc_vs_lambda1 (fig_mc_sweep_scatter): memory capacity against lambda_1, every run of the grid');
sources = [{S_ex.source}, cellstr(string(S_map.source)), cellstr(string(S_sct.source))];
parts = cell(0, 2);
if cfg.save
    R = S_map.responsive_share;
    delays = arrayfun(@(a) ['- ' strjoin(cellstr(a.Title.String), ' ')], ax_a, 'UniformOutput', false);
    parts = {['Panel A, ' panels.A], strjoin([{'The delays drawn and the score of each reconstruction:', ''}, delays(:)'], newline)};
    for c = 1:nc
        parts(end + 1, :) = {sprintf('Panel B, map %d of %d (fig_mc_gain_bias_maps, rows {''responsive''})', c, nc), ...
            T.grid_table(sprintf('Responsive share (1 - pinned - silent), median over seeds, %s', names{c}), ...
            R.median(:, :, c), R.second_levels, R.first_levels, 'self-coupling s', 'gain')}; %#ok<AGROW>
    end
    parts(end + 1, :) = {['Panel C, ' panels.C], strjoin([{'The statistics the source figure prints on its panels (not drawn here):', ''}, stats_text], newline)};
    if ~isempty(clip_note)
        % The key's note on the clipped mark, left off the panel (rule 1, no
        % explanatory text on a panel; Plan 3 review) and kept for the caption.
        parts(end + 1, :) = {'Panel C, for the caption', sprintf(['Open squares: runs with %s ' ...
            '(the source key''s note, not drawn on the figure).'], clip_note)};
    end
end
caption = sprintf(['Supplementary Figure S9, memory (drawn by fig_supp9_grouped; run `%s`, grid from `%s`; conditions %s). ' ...
    '(A) Example reconstructions of the past input (black: target; colour: the readout) at four delays, with the ' ...
    'score of each; the example''s two summary axes are left out. (B) Responsive share of neurons, 1 - frac_pinned ' ...
    '- frac_silent (time-mean rate neither above 0.9 nor below 0.01), median over seeds, over gain by self-coupling ' ...
    's, one map per condition in the layout of Figure 7E; white 0, dark green 1; star: the preset; self-coupling ' ...
    'marked at -4, 0 and 6. (C) Memory capacity against lambda_1 with the task input on for every run of that grid, ' ...
    'coloured by gain, marker area by |s| (filled s >= 0, ring s < 0), open squares at the edge where |lambda_1| ' ...
    'exceeds the axis; the statistics boxes of the source figure are in this table, not on the panels. Canvas %d x ' ...
    '%d px (canvas_w of manuscript_style), type in points by role.'], ...
    cfg.run_dir, cfg.mc_run_dir, strjoin(names, ', '), getappdata(fig, 'composite_size'));
out_dir = '';
if cfg.save; out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath')); end
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Supplementary Figure S9: memory', caption, parts, sources);
out.panels = panels;
end

%% ------------------------------------------------------------------------
function e = registry_entry(figures, name)
for k = 1:numel(figures)
    if strcmp(figures{k}.name, name); e = figures{k}; return; end
end
error('fig_supp9_grouped:NoEntry', 'The bundle registers no figure %s.', name);
end

function note = respace_size_key(key, W, y, h)
% The marker-size key under the scatter panels, laid out for the 9 in
% canvas. The source spaces its reference markers for its own small labels
% on a wider sheet, with the heading left of the markers and the note on
% the clipped mark right of them; at tick_fs on 9 in the three do not fit
% in one line. Here the heading stays at the left and the markers and their
% values take the rest of the width; the clipped mark and its note are
% removed and the note returned (until the Plan 3 review they went under
% the heading, level with the values).
HEAD_PX = 250; LEFT = 16; RIGHT_PAD = 30;
tx = findall(key, 'Type', 'text');
is_val = arrayfun(@(t) strcmp(t.VerticalAlignment, 'top'), tx);
is_head = arrayfun(@(t) strcmp(t.HorizontalAlignment, 'right'), tx);
is_note = ~is_val & ~is_head;
n = nnz(is_val);
width = W - 2 * LEFT;
u = (width - HEAD_PX - RIGHT_PAD) / (n + 0.1);           % px per marker step: the first at 0.4 + 0.6 steps
key.Units = 'pixels';
key.Position = [LEFT, y, width, h];
x0 = 0.4 - HEAD_PX / u;
key.XLim = [x0, x0 + width / u];
T = composite_tools();
for k = 1:numel(tx); tx(k).String = T.sentence(tx(k).String); end          % 'Self-coupling s (R_W):'
% The clipped mark and its note ("|lambda_1| > 5, drawn at the edge") leave
% the figure (Plan 3 review, rule 1: no explanatory text on a panel); the
% note is returned for the table, where the caption takes it from.
note = '';
if any(is_note)
    note = cellfun(@(s) strtrim(strjoin(cellstr(s), ' ')), {tx(is_note).String}, 'UniformOutput', false);
    note = strjoin(note(~cellfun(@isempty, note)), '; ');
    delete(tx(is_note));
    delete(findall(key, 'Type', 'scatter', 'Marker', 'square'));
end
T.margin(ancestor(key, 'figure'), 2);                                      % the labels are the lowest thing drawn
end

function x = pixel_x(ax)
p = getpixelposition(ax, true);
x = p(1);
end

function lines = scatter_statistics(src)
% The text of the scatter panels' statistics boxes, one line per panel, left to right.
ax = findall(src, 'Type', 'axes');
ax = ax(arrayfun(@(a) ~isempty(findall(a, 'Tag', 'mc_sweep_scatter_stats')), ax));
[~, order] = sort(arrayfun(@(a) pixel_x(a), ax)); ax = ax(order);
lines = cell(1, numel(ax));
for k = 1:numel(ax)
    t = findall(ax(k), 'Tag', 'mc_sweep_scatter_stats');
    lines{k} = sprintf('- %s: %s', strjoin(cellstr(ax(k).Title.String), ' '), ...
        strjoin(cellstr(t(1).String), '; '));
end
end
