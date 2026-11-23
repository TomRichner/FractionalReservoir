function out = fig_supp9_grouped(cfg)
% FIG_SUPP9_GROUPED Supplementary Figure S9, memory details; saved data only.
%
%   out = FIG_SUPP9_GROUPED('run_dir', d, 'preset_name', p, 'out_dir', o)
%   out = FIG_SUPP9_GROUPED('run_dir', d, 'mc_run_dir', g, ...)   % the grid from another run
%
% The figure plan of 2026-11-21 (figures_main_supp_appen.md, S9): what backs
% main Figure 7 and ties responsiveness to memory. Every panel is drawn by
% the function that owns it, invisibly and unsaved, and its axes are copied
% here (composite_tools); nothing is simulated.
%   (A) example reconstructions at several delays: the four trace axes of
%       fig_memory_capacity_example (its two summary axes are left out;
%       Figure 7B and panel (B) here carry them for the ensemble).
%   (B) which state variables carry the leading perturbation, against the
%       longest SFA time constant: the second figure of fig_sfa_EOC_allStd
%       (was 6B), and beside it the cumulative memory capacity against
%       delay, the first axes of fig_memory_capacity (was the left half of
%       6C).
%   (C) the responsive share of neurons over the gain by self-coupling
%       grid, one map per condition, in the layout of the capacity maps of
%       Figure 7D: fig_mc_gain_bias_maps with rows {'responsive'} (median
%       over seeds of 1 - frac_pinned - frac_silent, the quantity
%       'frac_on_slope' of mc_sweep_quantity).
%   (D) memory capacity against lambda_1 for every run of that grid:
%       fig_mc_vs_lambda1 with the bundle's registry arguments. Its
%       statistics boxes are not drawn; their text is in the table.
%
% WHICH RUN: as fig_main7_memory. run_dir holds mc_example, the tau sweep
% and memory_capacity; the gain by self-coupling grid is read from
% mc_run_dir, by default the same run.
%
% Options: run_dir, mc_run_dir, preset_name, run_mode, figure_set (the
% bundle whose registry entries style the maps and the scatter),
% conditions (empty = the figure set's), out_dir, save, visible,
% fig_tag.
%
% out: figs, files, source, panels (letter -> source function).
%
% See also: composite_tools, fig_main7_memory, fig_mc_gain_bias_maps,
%           fig_mc_vs_lambda1, fig_memory_capacity_example

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
    cfg.fig_tag     (1,:) char    = 'Fig_Supp9_Grouped'

end
setup_paths();
T = composite_tools();
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
tau_args = [{'run_dir', cfg.run_dir, 'conditions', names}, quiet];
if ~isempty(cfg.preset_name); tau_args = [tau_args, {'preset_name', cfg.preset_name}]; end
S_tau = fig_sfa_EOC_allStd(tau_args{:});
S_mc  = fig_memory_capacity('run_dir', cfg.run_dir, 'conditions', names, quiet{:});
S_map = fig_mc_gain_bias_maps('run_dir', cfg.mc_run_dir, 'conditions', names, quiet{:}, ...
    map_entry.args{:}, 'rows', {'responsive'});
S_sct = fig_mc_vs_lambda1('run_dir', cfg.mc_run_dir, 'conditions', names, quiet{:}, sct_entry.args{:});
src_figs = [S_ex.figs(:); S_tau.figs(:); S_mc.figs(:); S_map.figs(:); S_sct.figs(:)];
cleanup = onCleanup(@() T.close(src_figs));
assert(numel(S_tau.figs) >= 2, 'fig_supp9_grouped:NoBlocks', ...
    'fig_sfa_EOC_allStd returned no leading-vector figure (the tau sweep saved no block fractions).');

%% The canvas
W = 1100; H = 958;           % under 0.9 of the workstation screen: exported cropped as usual
fig = T.canvas(W, H, cfg.visible);
Y1 = 700; H1 = 215;            % row 1: (A), (B)
Y2 = 385; MH = 190; MW = 165; GAP = 26; X0 = 105;   % row 2: (C), the maps' geometry of Figure 7
Y3 = 140; H3 = 150;            % row 3: (D); its size legend below

% (A): the trace axes of the example (titled 'delay = ...').
ex_axes = findall(S_ex.figs, 'Type', 'axes');
is_trace = arrayfun(@(a) startsWith(strjoin(cellstr(a.Title.String), ' '), 'delay'), ex_axes);
assert(any(is_trace), 'fig_supp9_grouped:NoTraces', ...
    'fig_memory_capacity_example drew no axes titled ''delay = ...''.');
tx = findall(ex_axes(is_trace), 'Type', 'text');            % the source's own panel letter
delete(tx(arrayfun(@(t) ~isempty(regexp(strjoin(cellstr(t.String), ' '), '^\([a-z]\)$', 'once')), tx)));
A = T.copy(ex_axes(is_trace), fig, false);
[~, order] = sort(arrayfun(@(a) -a.Position(2), A.axes)); ax_a = A.axes(order);   % top first
T.place(ax_a, [40 Y1 - 20 385 H1 + 20]);

% (B): leading-vector fractions (second figure of the tau sweep), then the
% cumulative capacity (first axes of fig_memory_capacity).
Bv = T.copy(S_tau.figs(2), fig);
assert(isscalar(Bv.axes), 'fig_supp9_grouped:TauPanels', ...
    'Expected one leading-vector panel from the tau sweep, found %d. Select one of its conditions.', ...
    numel(Bv.axes));
ax_b1 = Bv.axes;
T.place(ax_b1, [560 Y1 215 H1]);
mc_axes = findall(S_mc.figs, 'Type', 'axes');
assert(numel(mc_axes) == 3, 'fig_supp9_grouped:MemoryAxes', ...
    'Expected three memory-capacity axes, found %d.', numel(mc_axes));
[~, order] = sort(arrayfun(@(a) pixel_x(a), mc_axes)); mc_axes = mc_axes(order);

Bc = T.copy(mc_axes(1), fig, false);
ax_b2 = Bc.axes;
T.place(ax_b2, [885 Y1 190 H1]);

% (C): the responsive-share maps and their colour bar.
Cm = T.copy(S_map.figs, fig);
assert(numel(Cm.axes) == nc && isscalar(Cm.colorbars), 'fig_supp9_grouped:Maps', ...
    'Expected %d responsive-share maps and one colour bar.', nc);
[~, order] = sort(arrayfun(@(a) a.Position(1), Cm.axes)); maps = Cm.axes(order);
for c = 1:nc
    maps(c).Position = [X0 + (c - 1) * (MW + GAP), Y2, MW, MH];
end
x_cb = X0 + nc * MW + (nc - 1) * GAP + 14;
Cm.colorbars.Position = [x_cb, Y2, 14, MH];

% (D): the scatter panels, their colour bar and the marker-size legend.
stats_text = scatter_statistics(S_sct.figs);
delete(findall(S_sct.figs, 'Tag', 'mc_sweep_scatter_stats'));
Ds = T.copy(S_sct.figs, fig);
is_key = arrayfun(@(a) strcmp(a.Tag, 'mc_sweep_scatter_size_legend'), Ds.axes);
sct = Ds.axes(~is_key);
[~, order] = sort(arrayfun(@(a) a.Position(1), sct)); sct = sct(order);
T.place(sct, [95 Y3 880 H3], Ds.colorbars);
if any(is_key)
    Ds.axes(is_key).Position = [95 12 880 60];
end

%% Style
T.style(fig);
st = manuscript_style();
% (A)
set(ax_a, 'TitleFontSizeMultiplier', 1);
% (B)
tau_name = tau_condition(S_tau.source, names);
set(ax_b1, 'XScale', 'linear');          % as style_main6_grouped: the levels are evenly spaced
xlabel(ax_b1, {'Longest SFA Time', 'Constant (s)'}, 'Interpreter', 'none');
ylabel(ax_b1, {'Leading-Vector', 'Squared-Norm Fraction'}, 'Interpreter', 'none');
title(ax_b1, '');
ax_b1.YTick = [0 0.5 1];
% The state coordinates' palette and markers of style_main6_grouped (their
% own encoding, apart from the condition colours).
state_names = {'SFA', 'STD', 'x'};
state_colors = [213 94 0; 0 158 115; 142 68 173] / 255;
state_markers = {'o', 's', '^'};
for k = 1:numel(state_names)
    h = findall(ax_b1, 'Type', 'line', 'DisplayName', state_names{k});
    if k == 3 && isempty(h); h = findall(ax_b1, 'Type', 'line', 'DisplayName', 'Dendritic State'); end
    if isempty(h); continue; end
    set(h, 'Color', state_colors(k, :), 'MarkerFaceColor', state_colors(k, :), 'Marker', state_markers{k}, ...
        'MarkerSize', 5);
    if k == 3; set(h, 'DisplayName', 'Dendritic state'); end
end
set(Bv.legends, 'Location', 'none', 'Units', 'pixels');
ylabel(ax_b2, {'Cumulative', 'Memory Capacity'}, 'Interpreter', 'none');
ax_b2.YTick = ax_b2.YTick([1 end]);
% One legend for the conditions, from the lines of the cumulative-capacity
% axes (the source's own legend belongs to its R^2 axes, which is not
% drawn here), named by condition in the order drawn.
cond_lines = flipud(findall(ax_b2, 'Type', 'line'));
lg = gobjects(0);
if numel(cond_lines) == nc
    lg = legend(ax_b2, cond_lines, cellfun(@(n) T.condition_title(n), names, 'UniformOutput', false), ...
        'Box', 'off', 'FontSize', 14, 'NumColumns', 1, 'Location', 'none', 'Units', 'pixels');
end
% (C): ticks and labels as fig_main_mc_grouped styles the capacity maps.
GAIN_TICKS = [0.05 0.25 1 5];
for c = 1:nc
    ax = maps(c);
    title(ax, T.condition_title(names{c}, true), 'FontWeight', 'normal', 'FontSize', 14, ...
        'Color', st.condition_color(names{c}), 'Interpreter', 'tex');
    g_levels = ax.XTick;
    set(ax, 'XTick', log10(GAIN_TICKS), 'XTickLabel', compose('%g', GAIN_TICKS), ...
        'XTickLabelRotation', 0, 'XMinorTick', 'on');
    ax.XAxis.MinorTickValues = g_levels;
    s_levels = ax.YTick;
    set(ax, 'YTick', s_levels(mod(s_levels, 2) == 0), 'YMinorTick', 'on');
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
Cm.colorbars.Label.String = 'Responsive Share';
% (D)
for c = 1:numel(sct)
    title(sct(c), T.condition_title(names{c}), 'FontWeight', 'normal', 'Color', sct(c).Title.Color, ...
        'Interpreter', 'tex');
    xlabel(sct(c), 'Lyapunov Exp., \lambda_1 (s^{-1})', 'Interpreter', 'tex');
end
ylabel(sct(1), 'Memory Capacity');
if ~isempty(Ds.colorbars); Ds.colorbars.Label.String = 'Gain'; end

% The two keys of (B) go under it, in the free space right of the maps: the
% state coordinates' above, the conditions' below (inside the narrow
% leading-vector panel a key covers the curves). Legends settle their size
% only after a redraw; then they are placed.
x_keys = x_cb + 122;
set([Bv.legends(:); lg(:)], 'ItemTokenSize', [22 18]);
for pass = 1:2
    drawnow;
    y = Y2 + MH + 25;
    if ~isempty(Bv.legends)
        p = Bv.legends.Position;
        y = y - p(4);
        Bv.legends.Position = [x_keys, y, p(3), p(4)];
    end
    if ~isempty(lg)
        p = lg.Position;
        y = y - p(4) - 16;
        lg.Position = [x_keys, y, p(3), p(4)];
    end
end

%% Panel letters
T.letter(ax_a(1), 'A', -30, 22);
T.letter(ax_b1, 'B', -95, 12);
T.letter(maps(1), 'C', -95, 50);
T.letter(sct(1), 'D', -80, 25);

%% Save
panels = struct( ...
    'A', 'fig_memory_capacity_example (its four trace axes): example reconstructions at several delays', ...
    'B', ['fig_sfa_EOC_allStd (second figure): leading-vector squared-norm fractions against the longest SFA ' ...
          'time constant; fig_memory_capacity (first axes): cumulative memory capacity against delay'], ...
    'C', 'fig_mc_gain_bias_maps with rows {''responsive''}: responsive share of neurons over gain by self-coupling', ...
    'D', 'fig_mc_vs_lambda1 (fig_mc_sweep_scatter): memory capacity against lambda_1, every run of the grid');
sources = [{S_ex.source, S_tau.source, S_mc.source}, cellstr(string(S_map.source)), cellstr(string(S_sct.source))];
parts = cell(0, 2);
if cfg.save
    R = S_map.responsive_share;
    delays = arrayfun(@(a) ['- ' strjoin(cellstr(a.Title.String), ' ')], ax_a, 'UniformOutput', false);
    parts = {['Panel A, ' panels.A], strjoin([{'The delays drawn and the score of each reconstruction:', ''}, delays(:)'], newline); ...
        'Panel B, left (fig_sfa_EOC_allStd, second figure)', T.axes_table(ax_b1, 'The squared-norm fraction of the leading Lyapunov vector in each block of state variables'); ...
        'Panel B, right (fig_memory_capacity, first axes)', T.axes_table(ax_b2, 'The cumulative memory capacity, mean over trials, per condition')};
    for c = 1:nc
        parts(end + 1, :) = {sprintf('Panel C, map %d of %d (fig_mc_gain_bias_maps, rows {''responsive''})', c, nc), ...
            T.grid_table(sprintf('Responsive share (1 - pinned - silent), median over seeds, %s', names{c}), ...
            R.median(:, :, c), R.second_levels, R.first_levels, 'self-coupling s', 'gain')}; %#ok<AGROW>
    end
    parts(end + 1, :) = {['Panel D, ' panels.D], strjoin([{'The statistics the source figure prints on its panels (not drawn here):', ''}, stats_text], newline)};
end
caption = sprintf(['Supplementary Figure S9, memory (drawn by fig_supp9_grouped; run `%s`, grid from `%s`; conditions %s). ' ...
    '(A) Example reconstructions of the past input (black: target; colour: the readout) at four delays, with the ' ...
    'score of each; the example''s two summary axes are left out. (B) Left: the share of the squared norm of the ' ...
    'leading Lyapunov vector in the SFA, depression and dendritic coordinates against the longest SFA time ' ...
    'constant (%s; the state coordinates have their own colours and markers). Right: cumulative memory capacity ' ...
    'against delay, per condition. (C) Responsive share of neurons, 1 - frac_pinned - frac_silent (time-mean rate ' ...
    'neither above 0.9 nor below 0.01), median over seeds, over gain by self-coupling s, one map per condition in ' ...
    'the layout of Figure 7D; white 0, dark green 1; star: the preset. (D) Memory capacity against lambda_1 with the ' ...
    'task input on for every run of that grid, coloured by gain, marker area by |s| (filled s >= 0, ring s < 0), ' ...
    'open squares at the edge where |lambda_1| exceeds the axis; the statistics boxes of the source figure are in ' ...
    'this table, not on the panels. Canvas %d x %d px, text 14 pt.'], ...
    cfg.run_dir, cfg.mc_run_dir, strjoin(names, ', '), tau_name, W, H);
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

function name = tau_condition(tau_dir, requested)
f = fullfile(tau_dir, 'param_space_summary.mat');
if isfile(f)
    S = load(f, 'summary_data'); conds = S.summary_data.conditions;
else
    psa = ParamSpaceAnalysis2.from_dir(tau_dir); conds = psa.conditions;
end
held = cellfun(@(c) c.name, conds, 'UniformOutput', false);
hit = requested(ismember(requested, held));
assert(isscalar(hit), 'fig_supp9_grouped:TauCondition', ...
    'The tau sweep holds %d of the selected conditions; panel (B) draws one.', numel(hit));
name = hit{1};
end
