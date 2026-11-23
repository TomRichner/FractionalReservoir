function out = fig_main7_memory(cfg)
% FIG_MAIN7_MEMORY Main Figure 7, memory: slow adaptive states hold a trace, over a wide range of connectivity; saved data only.
%
%   out = FIG_MAIN7_MEMORY('run_dir', d, 'preset_name', p, 'out_dir', o)
%   out = FIG_MAIN7_MEMORY('run_dir', d, 'mc_run_dir', g, ...)   % the grids from another run
%
% The figure plan of 2026-11-21 (figures_main_supp_appen.md, Figure 7)
% merges the former Figures 6 and 7. Every panel is drawn by the function
% that owns it, invisibly and unsaved, and its axes are copied here
% (composite_tools); nothing is simulated.
%   (A) lambda_1 against the longest SFA time constant: the first figure of
%       fig_sfa_EOC_allStd (the tau sweep; was 6A), with the display rules
%       of style_main6_grouped (lower limit -0.2, no mean markers).
%   (B) reconstruction score against delay: the R^2 axes of
%       fig_memory_capacity (was the right half of 6C; the cumulative
%       capacity half is in S9).
%   (C) memory horizon, paired across conditions: the horizon axes of
%       fig_memory_capacity (was 6D).
%   (D) memory capacity over gain by self-coupling, one map per condition:
%       the top row of fig_main_mc_grouped (fig_mc_gain_bias_maps; was 7A).
%   (E) memory capacity against input noise: the noise panel of
%       fig_main_mc_grouped (fig_mc_noise_sweeps; was 7C).
%   (F) lambda_1 over the same grid: the bottom row of fig_main_mc_grouped
%       (was 7B).
% One legend for the conditions, lower right (that of fig_main_mc_grouped).
%
% WHICH RUN. run_dir holds the tau sweep and the memory_capacity stage; the
% gain by self-coupling grid and the 1-D noise sweep are read from
% mc_run_dir, by default the same run (the A + C bundle, mu7chain_AC_config,
% keeps all of them in one run directory since FD-029). The paper's run of
% 2026-11-15 keeps them apart: data/mu7chain_AC_fastmedium and
% data/mu7chain_gainbias_fastmedium. A run without one of them is an error
% of the function that reads it (resolve_data_file, resolve_run_dir).
%
% Options
%   run_dir, mc_run_dir   above
%   preset_name           the bundle's preset (fig_sfa_EOC_allStd)
%   run_mode, figure_set  as fig_main_mc_grouped: the bundle whose registry
%                         entry styles the maps
%   conditions            names to draw, in order; empty = the figure set's
%   out_dir, save, visible, fig_tag
%
% out: figs, files, source, panels (letter -> source function).
%
% See also: composite_tools, fig_main_mc_grouped, fig_sfa_EOC_allStd,
%           fig_memory_capacity, fig_supp9_grouped

arguments
    cfg.run_dir     (1,:) char    = ''
    cfg.mc_run_dir  (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''
    cfg.run_mode    (1,:) char    = 'fastmedium'
    cfg.figure_set  (1,:) char    = '3cond'
    cfg.conditions                = {}
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = false
    cfg.verbose     (1,:) char    = 'minimal'   % accepted for a uniform call
    cfg.fig_tag     (1,:) char    = 'Fig_Main7_Grouped'

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
tau_args = [{'run_dir', cfg.run_dir, 'conditions', names}, quiet];
if ~isempty(cfg.preset_name); tau_args = [tau_args, {'preset_name', cfg.preset_name}]; end
S_tau = fig_sfa_EOC_allStd(tau_args{:});
S_mc  = fig_memory_capacity('run_dir', cfg.run_dir, 'conditions', names, quiet{:});
S_map = fig_main_mc_grouped('run_dir', cfg.mc_run_dir, 'run_mode', cfg.run_mode, ...
    'figure_set', cfg.figure_set, 'conditions', names, quiet{:});
src_figs = [S_tau.figs(:); S_mc.figs(:); S_map.figs(:)];
cleanup = onCleanup(@() T.close(src_figs));

%% The canvas: the map block of fig_main_mc_grouped below, the top row above
W = 1100;
H_MAPS = S_map.figs.Position(4);          % 685 px: its own canvas, kept pixel for pixel
assert(S_map.figs.Position(3) == W, 'fig_main7_memory:MapCanvas', ...
    'fig_main_mc_grouped drew on a canvas %g px wide; this figure is laid out for %d.', ...
    S_map.figs.Position(3), W);
TOP_Y = H_MAPS + 55; TOP_H = 170;
H = TOP_Y + TOP_H + 52;      % 962 px: under 0.9 of the workstation screen, so exported cropped as usual
fig = T.canvas(W, H, cfg.visible);

% (D), (E), (F) and the legend, where fig_main_mc_grouped put them.
delete(findall(S_map.figs, 'Tag', 'main_mc_panel_letter'));
M = T.copy(S_map.figs, fig);
maps = gobjects(2, nc);
for r = 1:2
    for c = 1:nc
        maps(r, c) = findobj(M.axes, 'flat', 'Tag', sprintf('main_mc_r%d_c%d', r, c));
    end
end
ax_noise = M.axes(~ismember(M.axes, maps(:)));
assert(isscalar(ax_noise) && isscalar(M.legends), 'fig_main7_memory:MapBlock', ...
    'Expected the noise panel and one legend beside the maps of fig_main_mc_grouped.');

% (A): the lambda_1 panel of the tau sweep (its first figure).
A = T.copy(S_tau.figs(1), fig);
ax_a = A.axes;
assert(isscalar(ax_a), 'fig_main7_memory:TauPanels', ...
    ['Expected one lambda_1 panel from the tau sweep, found %d: the sweep holds more than one of the ' ...
     'selected conditions. Select one.'], numel(ax_a));
T.place(ax_a, [95 TOP_Y 270 TOP_H]);

% (B), (C): the reconstruction and the horizon axes of fig_memory_capacity
% (its three axes left to right: cumulative capacity, R^2, horizon).
mc_axes = findall(S_mc.figs, 'Type', 'axes');
assert(numel(mc_axes) == 3, 'fig_main7_memory:MemoryAxes', ...
    'Expected three memory-capacity axes, found %d.', numel(mc_axes));
[~, order] = sort(arrayfun(@(a) pixel_x(a), mc_axes)); mc_axes = mc_axes(order);
B = T.copy(mc_axes(2), fig, false); ax_b = B.axes;
C = T.copy(mc_axes(3), fig, false); ax_c = C.axes;
T.place(ax_b, [495 TOP_Y 215 TOP_H]);
T.place(ax_c, [800 TOP_Y 275 TOP_H]);

%% Style
T.style(fig);
% (A), as style_main6_grouped draws it (Tom, 2026-11-17).
tau_name = tau_condition(S_tau.source, names);
set(ax_a, 'XScale', 'linear');
xlabel(ax_a, 'Longest SFA Time Constant (s)', 'Interpreter', 'none');
ylabel(ax_a, {'Lyapunov Exp.,', '\lambda_1 (s^{-1})'}, 'Interpreter', 'tex');   % two lines: at 15.4 pt one line ran into (A)
title(ax_a, T.condition_title(tau_name), 'FontWeight', 'normal', 'Interpreter', 'tex', ...
    'Color', st.condition_color(tau_name));
zero_ref = findall(ax_a, 'Type', 'constantline');
zero_ref = zero_ref(arrayfun(@(h) h.Value == 0, zero_ref));
set(zero_ref, 'Color', st.zeroline_color, 'LineStyle', '--', 'LineWidth', st.zeroline_lw);
delete(findall(ax_a, 'Type', 'line', 'Marker', 'o', 'LineStyle', 'none'));
ax_a.YLim(1) = -0.2;
% (B)
ylabel(ax_b, {'Input Reconstruction', 'Score (R^2)'}, 'Interpreter', 'tex');
ax_b.YTick = [0 1];
% (C): ticks by condition NAME (the source labels them with the titles).
try
    horizon_names = grouped_condition_label(cellstr(ax_c.XTickLabel), 'name_of_title');
catch ME
    if ~strcmp(ME.identifier, 'grouped_condition_label:UnknownTitle'); rethrow(ME); end
    horizon_names = names;      % the source no longer labels by title: its order is the selection's
end
% The conditions by name, on two lines and in their colors, where "1TS" and
% "MTS" stood (polish list, 2026-11-22: the abbreviations are defined
% nowhere on the figure). The names say what the axis is, so it has no
% label of its own; there is no room for one under two lines.
% Text objects, not tick labels: a tick label of two lines is set flush
% left of its first line, not centered under its tick.
horizon_ticks(ax_c, horizon_names, st);
xlabel(ax_c, '');
ax_c.YTick = 0:6:ax_c.YLim(2);
% (D), (F): sentence-case condition titles, in the colour the source gave them.
for c = 1:nc
    title(maps(1, c), T.condition_title(names{c}, true), 'FontWeight', 'normal', 'FontSize', st.title_fs, ...
        'Color', maps(1, c).Title.Color, 'Interpreter', 'tex');
end
M.legends.String = cellfun(@(n) T.condition_title(n), names, 'UniformOutput', false);
set(M.legends, 'Units', 'pixels');
for pass = 1:2          % the legend settles its width only after a redraw
    drawnow;
    lp = M.legends.Position;
    M.legends.Position = [W - 10 - lp(3), lp(2), lp(3), lp(4)];
end

% (D): the gain tick labels under the memory-capacity maps too. The source
% labels the bottom row only, which in this figure is two rows down, under
% (F) (polish list, 2026-11-22).
for c = 1:nc
    maps(1, c).XTickLabel = maps(2, c).XTickLabel;
end
% (E): the labelled noise levels as their values (the source rounds to one
% significant figure: 0.025 read "3 x 10^-2"), and a break drawn on the
% axis between 0 and the lowest positive level, which are not a decade apart.
relabel_noise_axis(ax_noise, S_map.noise_levels);

%% The common canvas: the layout above, made on W px, stretched across canvas_w
T.style(fig);                % again: the labels and titles set since take their roles' sizes
T.fit(fig);

%% Panel letters
T.letter(ax_a, 'A', -100, 24);       % (A) to (C) level with the title of (A)
T.letter(ax_b, 'B', -100, 24);
T.letter(ax_c, 'C', -85, 24);
T.letter(maps(1, 1), 'D', -95, 55);
T.letter(ax_noise, 'E', -70, 12);
T.letter(maps(2, 1), 'F', -95, 30);

%% Save
panels = struct( ...
    'A', 'fig_sfa_EOC_allStd (first figure): lambda_1 against the longest SFA time constant, tau sweep', ...
    'B', 'fig_memory_capacity (second axes): reconstruction score R^2 against delay', ...
    'C', 'fig_memory_capacity (third axes): memory horizon, paired across conditions', ...
    'D', 'fig_main_mc_grouped (top row; fig_mc_gain_bias_maps): memory capacity over gain by self-coupling', ...
    'E', 'fig_main_mc_grouped (noise panel; fig_mc_noise_sweeps, grid noise1d_wiener): memory capacity against input noise', ...
    'F', 'fig_main_mc_grouped (bottom row; fig_mc_gain_bias_maps): lambda_1 over the same grid');
sources = [{S_tau.source, S_mc.source}, reshape(cellstr(string(S_map.source)), 1, [])];
parts = cell(0, 2);
if cfg.save
    D = resolve_mc_gain_sweep('run_dir', cfg.mc_run_dir, 'grid', 'selfcoupling');
    grid_part = @(what, field, c) T.grid_table(sprintf('%s, median over seeds, %s', what, names{c}), ...
        mc_gain_grid(D, names{c}, field), D.second_levels, D.first_levels, ...
        D.second_axis.column, D.first_axis.column);
    parts = {['Panel A, ' panels.A], T.axes_table(ax_a, 'The median over repetitions at each level of the sweep'); ...
             ['Panel B, ' panels.B], T.axes_table(ax_b, 'The mean over trials, per condition'); ...
             ['Panel C, ' panels.C], T.axes_table(ax_c, 'One gray line per trial across the conditions, the points per condition and its mean bar')};
    for c = 1:nc
        parts(end + 1, :) = {sprintf('Panel D, map %d of %d (fig_mc_gain_bias_maps)', c, nc), ...
            grid_part('Memory capacity', 'mc_total', c)}; %#ok<AGROW>
    end
    parts(end + 1, :) = {['Panel E, ' panels.E], T.axes_table(ax_noise, ...
        'The median over networks per condition; x is the position of the noise level on the broken axis (0, then the log-spaced levels)')};
    for c = 1:nc
        parts(end + 1, :) = {sprintf('Panel F, map %d of %d (fig_mc_gain_bias_maps)', c, nc), ...
            grid_part('lambda_1 with the task input on (1/s)', 'lle_input', c)}; %#ok<AGROW>
    end
end
caption = sprintf(['Main Figure 7, memory (drawn by fig_main7_memory; run `%s`, grids from `%s`; conditions %s). ' ...
    '(A) The leading Lyapunov exponent lambda_1 against the longest SFA time constant in the tau sweep (%s): median ' ...
    'line, interquartile band and gray binned counts; the lower y limit is -0.2 and the mean markers of the source ' ...
    'figure are not drawn (Tom, 2026-11-17). (B) Reconstruction score R^2 of the past input against delay, mean ' ...
    'and band over trials as fig_memory_capacity draws them; dotted, the threshold of the memory horizon. (C) Memory ' ...
    'horizon per trial, the same trial joined across conditions, and the mean per condition; the conditions are named under the axis, in their colors. (D) Memory capacity ' ...
    '(median over seeds) over gain by self-coupling s, one map per condition; star: the preset. (E) Memory capacity ' ...
    'of the preset network against the input noise sigma_u: median over networks and interquartile band; the axis ' ...
    'is broken between 0 and the lowest positive level (break marks), and the labelled ticks read their levels. (F) lambda_1 with the task input on over the same grid ' ...
    '(median over seeds; symmetric-log colour scale with its darkest colour at zero). One legend names the ' ...
    'conditions. Panels (D) to (F) are fig_main_mc_grouped''s, with its scales. Canvas %d x %d px (canvas_w of manuscript_style), type in points by role.'], ...
    cfg.run_dir, cfg.mc_run_dir, strjoin(names, ', '), tau_name, getappdata(fig, 'composite_size'));
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

function name = tau_condition(tau_dir, requested)
% The condition the tau sweep drew: the selection intersected with what the
% sweep holds (fig_sfa_EOC_allStd's rule), which the caller asserted is one.
f = fullfile(tau_dir, 'param_space_summary.mat');
if isfile(f)
    S = load(f, 'summary_data'); conds = S.summary_data.conditions;
else
    psa = ParamSpaceAnalysis2.from_dir(tau_dir); conds = psa.conditions;
end
held = cellfun(@(c) c.name, conds, 'UniformOutput', false);
hit = requested(ismember(requested, held));
assert(isscalar(hit), 'fig_main7_memory:TauCondition', ...
    'The tau sweep holds %d of the selected conditions; panel (A) draws one.', numel(hit));
name = hit{1};
end

function horizon_ticks(ax, names, st)
% The horizon panel's conditions under their ticks: the condition's tick
% name ('No adaptation', 'Single-timescale', 'Multiple-timescale') on two
% lines, broken at its space or after its hyphen, centered, in its
% condition color.
words = grouped_condition_label(names, 'tick');
assert(numel(words) == numel(ax.XTick), 'fig_main7_memory:HorizonTicks', ...
    'The horizon panel has %d ticks for %d conditions.', numel(ax.XTick), numel(words));
xt = ax.XTick; yl = ax.YLim;
ax.XTickLabel = {};
for k = 1:numel(words)
    parts = regexp(words{k}, '^(\S+?-|\S+) ?(\S+)$', 'tokens', 'once');
    if isempty(parts); parts = words(k); end
    text(ax, xt(k), yl(1) - 0.04 * diff(yl), parts, 'Color', st.condition_color(names{k}), 'FontSize', st.tick_fs, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Clipping', 'off', 'Interpreter', 'none', ...
        'Tag', 'main7_horizon_tick');
end
ax.YLim = yl;
end

function relabel_noise_axis(ax, levels)
% The noise panel's x axis: every labelled tick reads its level (%g of the
% saved value), and the axis line is broken between the 0 and the first
% positive level, with a slanted stroke on each side of the gap.
assert(numel(levels) == numel(ax.XTick), 'fig_main7_memory:NoiseLevels', ...
    'The noise panel has %d ticks and %d levels.', numel(ax.XTick), numel(levels));
lbl = cellstr(ax.XTickLabel);
for k = 1:numel(lbl)
    if ~isempty(strtrim(lbl{k})); lbl{k} = sprintf('%g', levels(k)); end
end
set(ax, 'XTickLabel', lbl, 'TickLabelInterpreter', 'tex');
setappdata(ax, 'noise_levels', levels);              % the level at each tick, for whoever reads the figure back
if ~(levels(1) == 0 && numel(levels) > 1); return; end
drawnow;
xt = ax.XTick; yl = ax.YLim; xl = ax.XLim;
p = getpixelposition(ax);
ux = diff(xl) / p(3); uy = diff(yl) / p(4);          % axis units per pixel
xm = mean(xt(1:2)); half = 5 * ux;                   % a gap of 10 px, midway
args = {'Clipping', 'off', 'HandleVisibility', 'off', 'Tag', 'main7_noise_break'};
patch(ax, xm + half * [-1 1 1 -1], yl(1) + 3 * uy * [-1 -1 1 1], 'w', 'EdgeColor', 'none', args{:});
for x0 = xm + half * [-1 1]
    plot(ax, x0 + 3 * ux * [-1 1], yl(1) + 7 * uy * [-1 1], '-', 'Color', ax.XColor, 'LineWidth', ax.LineWidth, args{:});
end
set(ax, 'XLim', xl, 'YLim', yl);
end
