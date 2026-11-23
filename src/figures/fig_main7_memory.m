function out = fig_main7_memory(cfg)
% FIG_MAIN7_MEMORY Main Figure 7, memory: slow adaptive states hold a trace, over a wide range of connectivity; saved data only.
%
%   out = FIG_MAIN7_MEMORY('run_dir', d, 'preset_name', p, 'out_dir', o)
%   out = FIG_MAIN7_MEMORY('run_dir', d, 'mc_run_dir', g, ...)   % the grids from another run
%
% Every panel is drawn by the function that owns it, invisibly and unsaved,
% and its axes are copied here (composite_tools); nothing is simulated.
% The layout of Tom's review (Plan2_figure_feedback.md, Figure 7, notes 1
% to 6; 2026-11-22), on the 9 in canvas:
%   top row, four panels:
%   (A) reconstruction score against delay: the R^2 axes of
%       fig_memory_capacity.
%   (B) its cumulative sum against delay, which is how memory capacity is
%       calculated: the cumulative-capacity axes of fig_memory_capacity
%       (was S9 panel C).
%   (C) memory horizon, paired across conditions: the horizon axes of
%       fig_memory_capacity, the conditions named under it ("No
%       adaptation", "STS", "MTS"); narrower than the other three.
%   (D) memory capacity against input noise: the noise panel of
%       fig_main_mc_grouped (fig_mc_noise_sweeps), redrawn on an axis without
%       a break: linear from 0 (noise_scale 'linear', the default) or
%       logarithmic over the positive levels only ('log').
%   second row:
%   (E) memory capacity over gain by self-coupling, one map per condition:
%       the top row of fig_main_mc_grouped (fig_mc_gain_bias_maps).
%   third row:
%   (F) lambda_1 over the same grid: the bottom row of fig_main_mc_grouped.
% On both maps the self-coupling axis is marked at 0 as well as at its
% ends (note 6). One key for the conditions, in one line above the top row.
% lambda_1 against the longest SFA time constant (the former panel A) is in
% the supplementary figure on the leading direction, fig_supp_leading_grouped
% (note 1).
%
% WHICH RUN. run_dir holds the memory_capacity stage; the gain by
% self-coupling grid and the 1-D noise sweep are read from mc_run_dir, by
% default the same run (the A + C bundle, mu7chain_AC_config, keeps all of
% them in one run directory since FD-029). A run without one of them is an
% error of the function that reads it (resolve_data_file, resolve_run_dir).
%
% Options
%   run_dir, mc_run_dir   above
%   preset_name           accepted for a uniform call (the bundle passes it)
%   run_mode, figure_set  as fig_main_mc_grouped: the bundle whose registry
%                         entry styles the maps
%   conditions            names to draw, in order; empty = the figure set's
%   noise_scale           (D)'s x axis: 'linear' (default) or 'log' (see above)
%   out_dir, save, visible, fig_tag
%
% out: figs, files, source, panels (letter -> source function).
%
% See also: composite_tools, fig_main_mc_grouped, fig_memory_capacity,
%           fig_supp9_grouped, fig_supp_leading_grouped

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
    cfg.fig_tag     (1,:) char    = 'Fig_Main7_Grouped'
    cfg.noise_scale (1,:) char {mustBeMember(cfg.noise_scale, {'linear', 'log'})} = 'linear'

end
setup_paths();
T = composite_tools();
st = manuscript_style();
if isempty(cfg.mc_run_dir); cfg.mc_run_dir = cfg.run_dir; end
if isempty(cfg.conditions)
    bundle = mu7chain_gainbias_config(cfg.run_mode, cfg.figure_set);
    cfg.conditions = bundle.conditions;
end
names = cellstr(cfg.conditions);
nc = numel(names);
quiet = {'save', false, 'visible', false};

%% The source figures
S_mc  = fig_memory_capacity('run_dir', cfg.run_dir, 'conditions', names, quiet{:});
S_map = fig_main_mc_grouped('run_dir', cfg.mc_run_dir, 'run_mode', cfg.run_mode, ...
    'figure_set', cfg.figure_set, 'conditions', names, quiet{:});
src_figs = [S_mc.figs(:); S_map.figs(:)];
cleanup = onCleanup(@() T.close(src_figs));

%% The canvas, laid out in pixels on canvas_w itself (fit has nothing to stretch)
W = round(st.canvas_w * get(groot, 'ScreenPixelsPerInch'));
% Rows from the bottom: (F) the lambda_1 maps, (E) the capacity maps with
% the condition titles above them, the top row (A) to (D), and the
% conditions' key above the top row.
MAP_H = 180; Y_F = 72; Y_E = Y_F + MAP_H + 50;    % the gap holds the gain labels of (E), then the letter of (F)
TOP_H = 200; TOP_Y = Y_E + MAP_H + 48 + 70;    % the titles of (E), then the x labels of the top row
KEY_Y = TOP_Y + TOP_H + 34;
H = KEY_Y + 28;
fig = T.canvas(W, H, cfg.visible);

% (D), (E), (F) and the key: the blocks of fig_main_mc_grouped.
delete(findall(S_map.figs, 'Tag', 'main_mc_panel_letter'));
M = T.copy(S_map.figs, fig);
maps = gobjects(2, nc);
for r = 1:2
    for c = 1:nc
        maps(r, c) = findobj(M.axes, 'flat', 'Tag', sprintf('main_mc_r%d_c%d', r, c));
    end
end
ax_noise = M.axes(~ismember(M.axes, maps(:)));
assert(isscalar(ax_noise) && isscalar(M.legends) && numel(M.colorbars) == 2, 'fig_main7_memory:MapBlock', ...
    'Expected the noise panel, one legend and two colour bars beside the maps of fig_main_mc_grouped.');
is_lle_cb = arrayfun(@(cb) contains(strjoin(cellstr(cb.Label.String), ' '), 'lambda'), M.colorbars);
cb_mc = M.colorbars(~is_lle_cb); cb_l = M.colorbars(is_lle_cb);

% (A), (B), (C): the reconstruction, cumulative-capacity and horizon axes of
% fig_memory_capacity (its three axes left to right: cumulative capacity,
% R^2, horizon).
mc_axes = findall(S_mc.figs, 'Type', 'axes');
assert(numel(mc_axes) == 3, 'fig_main7_memory:MemoryAxes', ...
    'Expected three memory-capacity axes, found %d.', numel(mc_axes));
[~, order] = sort(arrayfun(@(a) pixel_x(a), mc_axes)); mc_axes = mc_axes(order);
A = T.copy(mc_axes(2), fig, false); ax_a = A.axes;
B = T.copy(mc_axes(1), fig, false); ax_b = B.axes;
C = T.copy(mc_axes(3), fig, false); ax_c = C.axes;

%% Geometry (pixels)
% The top row: four panels, each with room at its left for its y label and
% tick labels (MARGIN); (C) is narrower than the rest (Tom, 2026-11-23: "the
% memory horizon figure can be more narrow so that the other three can be
% wider"; 236 px until then), its conditions abbreviated under it.
MARGIN = [114 82 50 60]; GAP = 8; RIGHT = 36;   % (A): a y label of three lines; RIGHT: room for (D)'s
                                                % x label, whose sigma_u touched the edge at 20 (Plan 3 review)
WC = 170;
w_rest = (W - sum(MARGIN) - 3 * GAP - RIGHT - WC) / 3;
widths = [w_rest w_rest WC w_rest];
top = [ax_a ax_b ax_c ax_noise];
x = 0;
for k = 1:4
    x = x + MARGIN(k);
    set(top(k), 'Units', 'pixels', 'PositionConstraint', 'innerposition', 'Position', [x TOP_Y widths(k) TOP_H]);
    x = x + widths(k) + GAP;
end
% The maps: three columns, the colour bar right of each row.
X0 = 82; MGAP = 14; CB_ROOM = 104;
MW = (W - X0 - 2 * MGAP - 10 - 14 - CB_ROOM) / nc;
for r = 1:2
    for c = 1:nc
        set(maps(r, c), 'Units', 'pixels', 'PositionConstraint', 'innerposition', ...
            'Position', [X0 + (c - 1) * (MW + MGAP), (r == 1) * Y_E + (r == 2) * Y_F, MW, MAP_H]);
    end
end
x_cb = X0 + nc * MW + (nc - 1) * MGAP + 10;
set(cb_mc, 'Units', 'pixels', 'Position', [x_cb, Y_E, 14, MAP_H]);
set(cb_l, 'Units', 'pixels', 'Position', [x_cb, Y_F, 14, MAP_H]);

%% Style
T.style(fig);
% (A)
% Axis labels in sentence case (Tom, 2026-11-23, G1).
ylabel(ax_a, {'Input', 'reconstruction', 'score (R^2)'}, 'Interpreter', 'tex');   % three lines: on two it was taller than the panel
ax_a.YTick = [0 1];
ax_a.YLabel.PositionMode = 'auto';      % beside its tick labels (the copy kept the source's place, 0.5 wide)
% (B)
ylabel(ax_b, {'Cumulative', 'memory capacity'}, 'Interpreter', 'none');
ax_b.YTick = ax_b.YTick([1 end]);
% (C): ticks by condition NAME (the source labels them with the titles).
try
    horizon_names = grouped_condition_label(cellstr(ax_c.XTickLabel), 'name_of_title');
catch ME
    if ~strcmp(ME.identifier, 'grouped_condition_label:UnknownTitle'); rethrow(ME); end
    horizon_names = names;      % the source no longer labels by title: its order is the selection's
end
% The conditions in their colors, "No adaptation", "STS" and "MTS" (Tom,
% 2026-11-23, replacing the full names on two lines of 2026-11-22; the key
% above the top row spells STS and MTS out). The names say what the axis
% is, so it has no label of its own. The x range is narrowed so that the
% names clear each other on the narrower panel; the outer names reach a
% little past the plot box, under the y tick labels and into the gap.
ax_c.XLim = [0.8, numel(ax_c.XTick) + 0.2];
horizon_ticks(ax_c, horizon_names, st);
xlabel(ax_c, '');
ax_c.YTick = 0:6:ax_c.YLim(2);
ylabel(ax_c, 'Memory horizon (s)');
% (D): memory capacity labelled every 5, as (B) (G5, Tom 2026-11-23: every 2
% to 12 was seven labels on a narrow axes). The noise axis is as Tom set it.
ax_noise.YTick = 0:5:ax_noise.YLim(2);
% (E), (F): sentence-case condition titles over (E), in the colour the
% source gave them; the gain tick labels under (E) too (the source labels
% the bottom row only; polish list, 2026-11-22); self-coupling marked at
% its ends and at 0 (note 6).
for c = 1:nc
    title(maps(1, c), T.condition_title(names{c}, true), 'FontWeight', 'normal', 'FontSize', st.title_fs, ...
        'Color', maps(1, c).Title.Color, 'Interpreter', 'tex');
    maps(1, c).XTickLabel = maps(2, c).XTickLabel;
    for r = 1:2
        mark_zero_coupling(maps(r, c), c == 1);
    end
end
% The growth rate by its one name (G3, Tom 2026-11-23), on two lines as before.
cb_l.Label.String = {[extractBefore(st.label_growth, ',') ','], [st.label_lle ' (s^{-1})']};
cb_l.Label.Interpreter = 'tex';
% The key: the conditions, in one line above the top row, whose lines it
% names ((A), (B), (D); (C) names them under its axis, (E) over its maps).
M.legends.String = cellfun(@(n) T.condition_title(n), names, 'UniformOutput', false);
set(M.legends, 'Units', 'pixels', 'Orientation', 'horizontal', 'NumColumns', nc, 'ItemTokenSize', [22 18]);
% (D): the noise levels at their true values, no break (Tom, 2026-11-23; until
% then 0 stood left of a break and the positive levels at log spacing).
noise_axis(ax_noise, S_map.noise_levels, cfg.noise_scale);

%% The common canvas (the layout is on canvas_w already: fit moves nothing)
T.style(fig);                % again: the labels and titles set since take their roles' sizes
T.fit(fig);
for pass = 1:2               % the key settles its width only after a redraw
    drawnow;
    lp = M.legends.Position;
    M.legends.Position = [(W - lp(3)) / 2, KEY_Y, lp(3), lp(4)];
end

%% Panel letters, level with each other in the top row
T.letter(ax_a, 'A', -MARGIN(1), 8);
T.letter(ax_b, 'B', -MARGIN(2), 8);
T.letter(ax_c, 'C', -MARGIN(3), 8);
T.letter(ax_noise, 'D', -MARGIN(4), 8);
T.letter(maps(1, 1), 'E', -X0 + 2, 30);
T.letter(maps(2, 1), 'F', -X0 + 2, 8);
T.margin(fig, 4);            % the gain label of (F) is the lowest thing drawn: keep the crop off it

%% Save
panels = struct( ...
    'A', 'fig_memory_capacity (second axes): reconstruction score R^2 against delay', ...
    'B', 'fig_memory_capacity (first axes): cumulative memory capacity against delay', ...
    'C', 'fig_memory_capacity (third axes): memory horizon, paired across conditions', ...
    'D', 'fig_main_mc_grouped (noise panel; fig_mc_noise_sweeps, grid noise1d_wiener): memory capacity against input noise', ...
    'E', 'fig_main_mc_grouped (top row; fig_mc_gain_bias_maps): memory capacity over gain by self-coupling', ...
    'F', 'fig_main_mc_grouped (bottom row; fig_mc_gain_bias_maps): lambda_1 over the same grid');
sources = [{S_mc.source}, reshape(cellstr(string(S_map.source)), 1, [])];
parts = cell(0, 2);
if strcmp(cfg.noise_scale, 'linear')
    noise_note = 'The median over networks per condition; x is the noise level sigma_u on a linear axis from 0';
    noise_caption = 'the axis is linear from 0 to the largest level, every level at its value. ';
else
    noise_note = ['The median over networks per condition; x is the noise level sigma_u on a logarithmic axis, ' ...
        'the positive levels only (the zero-noise level is not drawn)'];
    noise_caption = 'the axis is logarithmic and shows the positive levels only (the zero-noise level is not drawn). ';
end
if cfg.save
    D = resolve_mc_gain_sweep('run_dir', cfg.mc_run_dir, 'grid', 'selfcoupling');
    grid_part = @(what, field, c) T.grid_table(sprintf('%s, median over seeds, %s', what, names{c}), ...
        mc_gain_grid(D, names{c}, field), D.second_levels, D.first_levels, ...
        D.second_axis.column, D.first_axis.column);
    parts = {['Panel A, ' panels.A], T.axes_table(ax_a, 'The mean over trials, per condition'); ...
             ['Panel B, ' panels.B], T.axes_table(ax_b, 'The cumulative memory capacity, mean over trials, per condition'); ...
             ['Panel C, ' panels.C], T.axes_table(ax_c, 'One gray line per trial across the conditions, the points per condition and its median bar'); ...
             ['Panel D, ' panels.D], T.axes_table(ax_noise, ...
                noise_note)};
    for c = 1:nc
        parts(end + 1, :) = {sprintf('Panel E, map %d of %d (fig_mc_gain_bias_maps)', c, nc), ...
            grid_part('Memory capacity', 'mc_total', c)}; %#ok<AGROW>
    end
    for c = 1:nc
        parts(end + 1, :) = {sprintf('Panel F, map %d of %d (fig_mc_gain_bias_maps)', c, nc), ...
            grid_part('lambda_1 with the task input on (1/s)', 'lle_input', c)}; %#ok<AGROW>
    end
end
caption = sprintf(['Main Figure 7, memory (drawn by fig_main7_memory; run `%s`, grids from `%s`; conditions %s). ' ...
    '(A) Reconstruction score R^2 of the past input against delay, mean and band over trials as ' ...
    'fig_memory_capacity draws them; dotted, the threshold of the memory horizon. (B) Its cumulative sum against ' ...
    'delay, mean and band over trials: the memory capacity is the value it reaches. (C) Memory horizon per trial, ' ...
    'the same trial joined across conditions, and the median per condition (bar); the conditions are named under ' ...
    'the axis, in their colors. (D) Memory capacity of the preset network against the input noise sigma_u: median ' ...
    'over networks and interquartile band; ' noise_caption '(E) Memory capacity (median over seeds) over gain by ' ...
    'self-coupling s, one map per condition; star: the preset. (F) lambda_1 with the task input on over the same ' ...
    'grid (median over seeds; symmetric-log colour scale with its darkest colour at zero). The self-coupling axis ' ...
    'of (E) and (F) is marked at -4, 0 and 6 (0 lies between grid levels). One key above the top row names the ' ...
    'conditions. Panels (D) to (F) are fig_main_mc_grouped''s, with its scales. Canvas %d x %d px (canvas_w of ' ...
    'manuscript_style), type in points by role.'], ...
    cfg.run_dir, cfg.mc_run_dir, strjoin(names, ', '), getappdata(fig, 'composite_size'));
out_dir = '';
if cfg.save; out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath')); end
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Main Figure 7: memory', caption, parts, sources);
out.panels = panels;
end

%% ------------------------------------------------------------------------
function x = pixel_x(ax)
p = getpixelposition(ax, true);
x = p(1);
end

function mark_zero_coupling(ax, labelled)
% The self-coupling axis of a map marked at its ends and at 0 (Tom,
% 2026-11-22, Figure 7 note 6). The axis is linear in s (mc_axis_spacing),
% so 0 sits at its value between the grid levels; every level keeps its
% minor tick.
yl = ax.YLim;
ends = ax.YTick([1 end]);
ticks = unique([ends, 0]);
ticks = ticks(ticks >= yl(1) & ticks <= yl(2));
ax.YTick = ticks;
if labelled
    ax.YTickLabel = compose('%g', ticks);
else
    ax.YTickLabel = {};
end
end

function horizon_ticks(ax, names, st)
% The horizon panel's conditions under their ticks, in their condition
% colors (Tom, 2026-11-23: "you can abbreviate it as No adaptation, STS,
% and MTS"): 'No adaptation', 'STS' and 'MTS', one line each where the
% labels clear each other by MIN_GAP_PX, otherwise 'No adaptation' on two
% lines. A condition whose tick name is neither single- nor
% multiple-timescale keeps its tick name.
MIN_GAP_PX = 6;
words = grouped_condition_label(names, 'tick');
assert(numel(words) == numel(ax.XTick), 'fig_main7_memory:HorizonTicks', ...
    'The horizon panel has %d ticks for %d conditions.', numel(ax.XTick), numel(words));
words = regexprep(words, '^Single-timescale$', 'STS');
words = regexprep(words, '^Multiple-timescale$', 'MTS');
xt = ax.XTick; yl = ax.YLim;
ax.XTickLabel = {};
h = gobjects(1, numel(words));
for k = 1:numel(words)
    h(k) = text(ax, xt(k), yl(1) - 0.04 * diff(yl), words{k}, 'Color', st.condition_color(names{k}), 'FontSize', st.tick_fs, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Clipping', 'off', 'Interpreter', 'none', ...
        'Tag', 'main7_horizon_tick');
end
% Do neighbours clear each other on one line? (extents in pixels)
drawnow;
set(h, 'Units', 'pixels');
e = cell2mat(arrayfun(@(t) t.Extent, h(:), 'UniformOutput', false));
clear_px = e(2:end, 1) - (e(1:end-1, 1) + e(1:end-1, 3));
set(h, 'Units', 'data');
if any(clear_px < MIN_GAP_PX)
    for k = 1:numel(h)
        parts = regexp(words{k}, '^(\S+?-|\S+) (\S+)$', 'tokens', 'once');
        if ~isempty(parts); h(k).String = parts; end
    end
end
ax.YLim = yl;
end

function noise_axis(ax, levels, scale)
% The noise panel's x axis at the true noise levels (Tom, 2026-11-23: "why
% is there a break in the x axis for the noise vs MC? I think you should not
% have a break ... we need a third tick"). The source places the levels at
% the positions of its broken axis (0, then log-spaced), one tick per level;
% every drawn object's x is mapped from those positions to its level.
%   'linear'  every level, 0 included, at its value, from 0 to the largest
%             level, labelled ticks at 0, 0.01 and 0.02 (those inside it);
%   'log'     the positive levels only (the zero-noise level is left out),
%             labelled ticks at the decades inside the range, written as
%             powers of ten, and the far end of the range when only one
%             decade falls inside (both ends when none does).
% No break marks either way.
xt = ax.XTick;
assert(numel(levels) == numel(xt), 'fig_main7_memory:NoiseLevels', ...
    'The noise panel has %d ticks and %d levels.', numel(xt), numel(levels));
levels = levels(:)';
kids = findall(ax, '-property', 'XData');
kids = kids(kids ~= ax);
for h = reshape(kids, 1, [])
    x = interp1(xt, levels, double(h.XData), 'linear', 'extrap');
    y = h.YData;
    if strcmp(scale, 'log')
        keep = x > 0;
        if isa(h, 'matlab.graphics.primitive.Patch')
            set(h, 'XData', x(keep), 'YData', y(keep));
            continue
        end
        y(~keep) = NaN;          % a line ends, a point is not drawn
        x(~keep) = NaN;
    end
    set(h, 'XData', x, 'YData', y);
end
setappdata(ax, 'noise_levels', levels);              % the level at each x, for whoever reads the figure back
switch scale
    case 'linear'
        ticks = [0 0.01 0.02];
        ticks = ticks(ticks <= max(levels));
        set(ax, 'XScale', 'linear', 'XLim', [0 max(levels)], 'XTick', ticks, ...
            'XTickLabel', compose('%g', ticks), 'XMinorTick', 'off', 'TickLabelInterpreter', 'tex');
    case 'log'
        pos = levels(levels > 0);
        lim = [min(pos) max(pos)];
        dec = 10 .^ (ceil(log10(lim(1)) - 1e-9):floor(log10(lim(2)) + 1e-9));
        lbl = arrayfun(@(d) sprintf('10^{%d}', round(log10(d))), dec, 'UniformOutput', false);
        if numel(dec) == 1          % the end of the range farther from the decade
            [~, far] = max(abs(log10(lim) - log10(dec)));
            dec = [dec lim(far)]; lbl = [lbl {sprintf('%g', lim(far))}];
        elseif isempty(dec)
            dec = lim; lbl = compose('%g', lim);
        end
        [dec, o] = sort(dec);
        set(ax, 'XScale', 'log', 'XLim', lim, 'XTick', dec, 'XTickLabel', lbl(o), ...
            'XMinorTick', 'off', 'TickLabelInterpreter', 'tex');
end
end
