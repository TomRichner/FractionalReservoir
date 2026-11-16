function out = fig_mc_gain_bias_maps(cfg)
% FIG_MC_GAIN_BIAS_MAPS Memory capacity and lambda_1 on the grid of the sweep's two axes (gain by bias by default), one column per condition (block C).
%
%   out = FIG_MC_GAIN_BIAS_MAPS('run_dir', d)
%   out = FIG_MC_GAIN_BIAS_MAPS('data_file', f)
%   out = FIG_MC_GAIN_BIAS_MAPS('data', D, 'save', false)   % D from read_mc_gain_sweep, or a table
%
% Rows (Plan C, section 5):
%   1. memory capacity, median over seeds, one colour scale shared by every
%      condition, from 0 to the largest cell of any condition, white (none)
%      to dark red. The colour follows the SQUARE ROOT of memory capacity
%      (cfg.mc_scale = 'sqrt', the default; 'linear' is the alternative),
%      with the colour bar labelled in memory capacity itself: the floor of
%      a chaotic network (about 0.1) and of a silent one (about 0.5) then
%      differ visibly, while 5 against 14 keeps most of the colour range.
%      A log scale was not chosen: it gives the decade 0.1 to 1, which holds
%      only failed networks, as much colour as 1 to 10, which holds the
%      result.
%   2. lambda_1 with the memory-task input on, median over seeds, a
%      diverging scale centred on zero on a SYMMETRIC-LOG transform
%      (symlog_lle: linear within about +-cfg.lle_lin = 0.05 1/s,
%      logarithmic beyond, clipped at +-cfg.lle_clip = 10 1/s; the clipping
%      is stated in the colour bar's label). One scale serves both regimes:
%      without adaptation lambda_1 spans about -10 to +8, with chained MTS
%      -0.1 to +3, and a linear scale shared by both renders -0.1 as white.
%      cfg.lle_scale = 'linear' is the alternative: the same diverging map,
%      linear from -cfg.lle_clip to +cfg.lle_clip, values beyond clipped.
%   3. the interquartile range of memory capacity over seeds, shared scale,
%      only when cfg.show_iqr (default true) AND some cell has more than one
%      successful seed; with one seed the row would be blank, so it is
%      dropped and the caption says "one seed".
% x is the result's first axis (D.first_axis): the gain (level_of_chaos),
% or whatever property the sweep ran there, for example the self-coupling
% s; on a log axis for a log-spaced grid of positive levels and a linear
% one otherwise (cfg.first_scale, 'auto' by default; 'log' or 'linear' to
% force one). y is the result's second axis (D.second_axis): the bias S_c,
% or whatever property the sweep ran, for example the self-feedback k in
% units of R_W or the input noise sigma_u. Both axes' labels and their
% names in the caption come from mc_sweep_quantity. The ticks are the
% grid's own levels.
% AN AXIS OF A ZERO LEVEL THEN LOG-SPACED LEVELS (mc_gain_sweep_grid's
% 'zero_log', a noise amplitude; either axis, the first under 'auto') is
% drawn as its log ladder, every cell the same width, with the zero one log
% step below the lowest positive level and its tick labelled 0
% (mc_axis_spacing); the caption says the axis is broken there. Drawn at
% its values instead, such an axis would squeeze every level but the top
% one or two into a sliver beside the zero. Each cell is one grid point, drawn between the
% midpoints to its neighbours, so an irregular grid is drawn as it is.
%
% THE lambda_1 = 0 LINE is drawn on every row as the cell boundaries that
% separate a cell with lambda_1 > 0 from a neighbour with lambda_1 <= 0
% (the median map of the same condition), green dashed as everywhere in
% the manuscript. A staircase on the cell edges, not an interpolated
% contour: on a coarse grid an interpolated contour invents a position
% between levels and stops at the outermost cell centres.
%
% The preset's own point (D.preset_first, D.preset_second: gain 1, bias
% 0.2) is a white star, drawn where it is even when it is not a grid level. A
% cell with no successful seed is left grey (the axes colour).
%
% SATURATED AND SILENT CELLS: a cell whose median over seeds of the
% fraction of pinned neurons is >= cfg.flag_fraction (0.9) carries a small
% grey '+', one whose fraction of silent neurons is >= cfg.flag_fraction a
% small grey 'o', so a low memory capacity there is not read as chaos. 0.9
% rather than 0.5: a network with more than half its neurons silent can
% still hold a memory capacity of 5 to 9 (chained MTS at bias 0.55), so
% half silent is not a failed network. Flags are per condition.
%
% STYLE OPTIONS (FD-029, from Brian Lundstrom's feedback). Their defaults
% draw the figure as it was before them; a bundle's registry entry sets
% them for the paper's style.
%   ticks_top_right  (true) tick marks on the top and right edges too (the
%                    axes box); false draws ticks on the bottom and left only
%   show_zero_line   (true) the green dashed lambda_1 = 0 line
%   show_flags       (true) the '+' saturated and 'o' silent markers
%   show_iqr         (true) the spread row (above)
%   lle_colormap     'blue_white_red' (default): blue through WHITE at zero
%                    to red. Or a diverging map with a DARK center, so
%                    lambda_1 = 0 is the darkest band and white is never a
%                    value: 'guppy', 'redshift' or 'dark_magenta'
%                    (dark_center_colormaps says what each is and where it
%                    comes from)
%   mc_colormap      'white_orange_red' (default), or 'lle_center': white at
%                    zero, darkening to the center color of lle_colormap, so
%                    the two rows share one palette (needs a dark-centered
%                    lle_colormap)
%   missing_color    ([0.85 0.85 0.85], grey) the color of a cell with no
%                    successful job. With mc_colormap 'lle_center' on
%                    'redshift' the memory-capacity map is itself a gray
%                    scale, so give a color that is not gray there.
% A caption sentence about something switched off is dropped with it, and
% with a dark-centered lambda_1 map the caption says that the darkest color
% is lambda_1 = 0.
%
% WHICH CONDITIONS: cfg.conditions through select_conditions (empty = all the
% result holds, in its order); titles and colours by name (titles_by_name,
% mc_condition_colors). Any number of columns.
%
% The result is taken through resolve_mc_gain_sweep (data, data_file, or
% run_dir with cfg.grid naming which grid of the run directory).
%
% See also: fig_mc_gain_curves, fig_mc_robustness, mc_gain_bias_statistics, symlog_lle,
%           mc_gain_axis_scale

arguments
    cfg.verbose       (1,:) char    = 'minimal'
    cfg.run_dir       (1,:) char    = ''
    cfg.data_file     (1,:) char    = ''
    cfg.data                        = []
    cfg.grid          (1,:) char    = ''   % with run_dir: the grid's tag in the run directory
    cfg.out_dir       (1,:) char    = ''
    cfg.save          (1,1) logical = true
    cfg.visible       (1,1) logical = true
    cfg.preset_name   (1,:) char    = ''   % unused; accepted for a uniform call
    cfg.conditions                  = {}
    cfg.show_iqr      (1,1) logical = true
    cfg.flag_fraction (1,1) double  = 0.9
    cfg.mc_scale      (1,:) char    {mustBeMember(cfg.mc_scale, {'sqrt', 'linear'})} = 'sqrt'
    cfg.lle_scale     (1,:) char    {mustBeMember(cfg.lle_scale, {'symlog', 'linear'})} = 'symlog'
    cfg.lle_lin       (1,1) double  = 0.05 % symlog_lle's linear scale (1/s)
    cfg.lle_clip      (1,1) double  = 10   % |lambda_1| clipped here (1/s)
    cfg.fig_tag       (1,:) char    = 'Fig_MC_gain_bias_maps'
    cfg.first_scale   (1,:) char    {mustBeMember(cfg.first_scale, {'auto', 'log', 'linear'})} = 'auto'
    cfg.ticks_top_right (1,1) logical = true
    cfg.show_zero_line  (1,1) logical = true
    cfg.show_flags      (1,1) logical = true
    cfg.lle_colormap  (1,:) char    {mustBeMember(cfg.lle_colormap, ...
        {'blue_white_red', 'guppy', 'redshift', 'dark_magenta'})} = 'blue_white_red'
    cfg.mc_colormap   (1,:) char    {mustBeMember(cfg.mc_colormap, {'white_orange_red', 'lle_center'})} = 'white_orange_red'
    cfg.missing_color (1,3) double  {mustBeBetween(cfg.missing_color, 0, 1)} = [0.85 0.85 0.85]
end
DEFAULT_MISSING = [0.85 0.85 0.85];
if strcmp(cfg.mc_colormap, 'lle_center') && strcmp(cfg.lle_colormap, 'blue_white_red')
    error('fig_mc_gain_bias_maps:NoDarkCenter', ...
        ['mc_colormap ''lle_center'' darkens to the center color of a dark-centered lambda_1 map; ' ...
         'lle_colormap ''blue_white_red'' has a white center. Use guppy, redshift or dark_magenta.']);
end

setup_paths();
D = resolve_mc_gain_sweep('data', cfg.data, 'data_file', cfg.data_file, ...
    'run_dir', cfg.run_dir, 'grid', cfg.grid);
conds = select_conditions(D.conditions, cfg.conditions);
nc = numel(conds);
st = manuscript_style();
cond_titles = titles_by_name(conds);
cond_colors = mc_condition_colors(conds);

%% The grids
Q1 = mc_sweep_quantity(D.first_axis.column, D);    % the first axis (x): label and wording
Q2 = mc_sweep_quantity(D.second_axis.column, D);   % the second axis (y)
nb = numel(D.second_levels); ng = numel(D.first_levels);
MC = nan(nb, ng, nc); L = MC; IQ = MC; PIN = MC; SIL = MC; NOK = MC;
for c = 1:nc
    [MC(:, :, c), NOK(:, :, c)] = mc_gain_grid(D, conds{c}, 'mc_total');
    L(:, :, c)   = mc_gain_grid(D, conds{c}, 'lle_input');
    IQ(:, :, c)  = mc_gain_grid(D, conds{c}, 'mc_total', 'iqr');
    PIN(:, :, c) = mc_gain_grid(D, conds{c}, 'frac_pinned');
    SIL(:, :, c) = mc_gain_grid(D, conds{c}, 'frac_silent');
end
n_seeds = max([NOK(:); 0]);              % most successful seeds behind any cell
% A cell with fewer than two successful seeds has no spread to show: grey,
% not white, which would read as "the seeds agree".
IQ(NOK < 2) = NaN;
mc_max = max([MC(:); eps]);
switch cfg.mc_scale
    case 'sqrt'
        mc_f = @sqrt;
        mc_ticks = [0 0.1 0.5 1 2 5 10 15 20 30 50];
        mc_ticks = mc_ticks(mc_ticks <= mc_max * 1.0001);
        mc_cb = 'MC (median), square-root colour scale';
    otherwise
        mc_f = @(v) v;
        mc_ticks = [];
        mc_cb = 'MC (median)';
end
switch cfg.lle_scale
    case 'symlog'
        [~, l_tp, l_tl, l_lim] = symlog_lle(0, cfg.lle_lin, cfg.lle_clip);
        L_col = symlog_lle(L, cfg.lle_lin, cfg.lle_clip);
        l_cb = sprintf('\\lambda_1 (1/s), symmetric log, clipped at \\pm%g', cfg.lle_clip);
    otherwise
        l_lim = cfg.lle_clip * [-1 1];
        L_col = min(max(L, l_lim(1)), l_lim(2));   % NaN stays NaN (grey)
        l_tp = []; l_tl = {};
        l_cb = sprintf('\\lambda_1 (1/s), linear, clipped at \\pm%g', cfg.lle_clip);
end
iq_clim = [0, max([IQ(:); eps])];
[l_cmap, m_cmap] = row_colormaps(cfg.lle_colormap, cfg.mc_colormap);

rows = {struct('V', mc_f(MC), 'clim', mc_f([0 mc_max]), 'cmap', m_cmap, ...
               'label', 'Memory capacity', 'cb', mc_cb, ...
               'ticks', mc_f(mc_ticks), 'ticklabels', {num_labels(mc_ticks)})
        struct('V', L_col, 'clim', l_lim, 'cmap', l_cmap, ...
               'label', '\lambda_1, task input on', ...
               'cb', l_cb, ...
               'ticks', l_tp, 'ticklabels', {l_tl})};
if cfg.show_iqr && n_seeds > 1
    rows{end+1} = struct('V', IQ, 'clim', iq_clim, 'cmap', spread_cmap(), ...
                         'label', 'Spread over seeds', 'cb', 'MC IQR over seeds', ...
                         'ticks', [], 'ticklabels', {{}});
end
nr = numel(rows);

%% Draw
fig = figure('Color', 'w', 'Position', [60 60 270 * nc + 200, 260 * nr + 120]);
if ~cfg.visible; set(fig, 'Visible', 'off'); end
outer = tiledlayout(fig, nr, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
x_of = axis_transform(D.first_levels, cfg.first_scale);
% The second axis at its values, unless it is a zero then log-spaced axis
% (a noise amplitude), which is drawn as its log ladder (mc_axis_spacing).
[spacing2, y_of] = mc_axis_spacing(D.second_levels);
if ~strcmp(spacing2, 'zero_log'); y_of = @(v) v; end
broken = {};
if strcmp(spacing2, 'zero_log'); broken{end+1} = Q2.symbol; end
if strcmp(cfg.first_scale, 'auto') && strcmp(mc_axis_spacing(D.first_levels), 'zero_log')
    broken{end+1} = Q1.symbol;
end
x = x_of(D.first_levels); y = y_of(D.second_levels);
xe = edges(x); ye = edges(y);
xtl = arrayfun(@(v) sprintf('%.3g', v), D.first_levels, 'UniformOutput', false);
ytl = arrayfun(@(v) sprintf('%.3g', v), D.second_levels, 'UniformOutput', false);
ax_all = gobjects(nr, nc);
% A qualifier in parentheses goes to a second line, so that four long
% titles ("... Adaptation (product)") do not run together. When any title
% wraps, every title gets two lines (a blank second line where there is no
% qualifier), so the condition names stay on one baseline across the row.
titles = cell(1, nc);
for c = 1:nc
    titles{c} = strsplit(regexprep(cond_titles{c}, '\s+(\(.*\))$', '\n$1'), newline);
end
n_lines = max(cellfun(@numel, titles));
% The layout reserves too little height for a two-line title: exported, the
% wrapped second line then crosses the axes' top edge. One more blank line
% below every title lifts the text clear.
if n_lines > 1; n_lines = n_lines + 1; end
for c = 1:nc
    titles{c}(end+1:n_lines) = {' '};
end
for r = 1:nr
    R = rows{r};
    inner = tiledlayout(outer, 1, nc, 'TileSpacing', 'compact', 'Padding', 'tight');
    inner.Layout.Tile = r;
    for c = 1:nc
        ax = nexttile(inner); hold(ax, 'on');
        ax_all(r, c) = ax;
        draw_cells(ax, x, y, R.V(:, :, c));
        colormap(ax, R.cmap); clim(ax, R.clim);
        % Box 'on' mirrors the tick marks onto the top and right edges.
        set(ax, 'Color', cfg.missing_color, 'Layer', 'top', 'TickDir', 'out', ...
            'FontSize', st.tick_fs - 5, 'XTick', x, 'XTickLabel', xtl, ...
            'YTick', y, 'YTickLabel', ytl, 'Box', on_off(cfg.ticks_top_right), ...
            'XTickLabelRotation', 45);
        xlim(ax, xe([1 end])); ylim(ax, ye([1 end]));
        if cfg.show_flags
            mark_flags(ax, x, y, PIN(:, :, c) >= cfg.flag_fraction, SIL(:, :, c) >= cfg.flag_fraction);
        end
        if cfg.show_zero_line
            zero_boundary(ax, xe, ye, L(:, :, c), st);
        end
        plot(ax, x_of(D.preset_first), y_of(D.preset_second), 'p','MarkerSize', 11, ...
            'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
        if r == 1
            title(ax, titles{c}, 'FontSize', st.title_fs - 2, ...
                'FontWeight', 'normal', 'Color', cond_colors(c, :));
        end
        if c == 1
            ylabel(ax, {R.label, capitalized(Q2.label)}, 'FontSize', st.label_fs - 3);
        else
            set(ax, 'YTickLabel', []);
        end
        if r == nr
            xlabel(ax, capitalized(Q1.label), 'FontSize', st.label_fs - 3);
        else
            set(ax, 'XTickLabel', []);
        end
    end
    cb = colorbar(ax_all(r, nc));
    cb.Layout.Tile = 'east';
    cb.Label.String = R.cb; cb.FontSize = st.tick_fs - 5;
    if ~isempty(R.ticks)
        cb.Ticks = R.ticks; cb.TickLabels = R.ticklabels;
    end
end
if n_seeds <= 1
    seed_txt = 'One seed per cell';
else
    seed_txt = sprintf('Median over up to %d seeds', n_seeds);
end
if isequal(cfg.missing_color, DEFAULT_MISSING)
    grey_txt = 'Grey: no successful job';
else
    grey_txt = 'Background color: no successful job';
end
if nr == 3; grey_txt = [grey_txt ' (spread row: fewer than two)']; end
% Each sentence about a mark is there only when the mark is drawn.
line1 = [seed_txt '.'];
if cfg.show_zero_line
    line1 = [line1 ' Green dashed: edge between cells with \lambda_1 > 0 and \lambda_1 \leq 0.'];
end
if ~strcmp(cfg.lle_colormap, 'blue_white_red')
    line1 = [line1 ' The darkest color of the \lambda_1 row is \lambda_1 = 0.'];
end
line2 = sprintf('Star: preset (%s %.3g, %s %.3g).', Q1.symbol, D.preset_first, Q2.symbol, D.preset_second);
if cfg.show_flags
    line2 = [line2 sprintf(' + saturated, o silent (\\geq %.0f%% of neurons).', 100 * cfg.flag_fraction)];
end
caption = {line1, [line2 ' ' grey_txt '.']};
if ~isempty(broken)
    caption{end+1} = sprintf(['%s: 0, then log-spaced levels; cells evenly spaced, the axis broken ' ...
        'between 0 and the lowest positive level.'], capitalized(strjoin(broken, ' and ')));
end
title(outer, caption, 'FontSize', st.tick_fs - 5);

%% Save
out = struct('figs', fig, 'files', {{}}, 'source', D.source);
if cfg.save
    out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
    save_figure_stable(out_dir, cfg.fig_tag, fig);
    out.files = existing_outputs(out_dir, cfg.fig_tag);
end
end

%% ------------------------------------------------------------------------
function e = edges(v)
% Cell edges: midpoints, extended by half a step at each end.
v = v(:)';
if isscalar(v); e = v + [-0.5 0.5]; return; end
m = (v(1:end-1) + v(2:end)) / 2;
e = [v(1) - (m(1) - v(1)), m, v(end) + (v(end) - m(end))];
end

function draw_cells(ax, x, y, V)
xe = edges(x); ye = edges(y);
C = nan(numel(ye), numel(xe));
C(1:end-1, 1:end-1) = V;
[X, Y] = meshgrid(xe, ye);
surface(ax, X, Y, zeros(size(X)), C, 'FaceColor', 'flat', 'EdgeColor', 'none');
end

function mark_flags(ax, x, y, pinned, silent)
[X, Y] = meshgrid(x, y);
g = [0.35 0.35 0.35];
if any(pinned(:))
    plot(ax, X(pinned), Y(pinned), '+', 'Color', g, 'MarkerSize', 5, 'LineWidth', 0.8);
end
if any(silent(:))
    plot(ax, X(silent), Y(silent), 'o', 'Color', g, 'MarkerSize', 4, 'LineWidth', 0.8);
end
end

function zero_boundary(ax, xe, ye, Lc, st)
% Cell edges between a cell with lambda_1 > 0 and a neighbour with
% lambda_1 <= 0, both finite, drawn as one NaN-separated dashed line.
pos = Lc > 0; fin = isfinite(Lc);
X = []; Y = [];
[nb, ng] = size(Lc);
for b = 1:nb
    for g = 1:ng - 1       % neighbours along the first axis: a vertical edge
        if fin(b, g) && fin(b, g + 1) && pos(b, g) ~= pos(b, g + 1)
            X = [X, xe(g + 1), xe(g + 1), NaN]; Y = [Y, ye(b), ye(b + 1), NaN]; %#ok<AGROW>
        end
    end
end
for b = 1:nb - 1           % neighbours along the second axis: a horizontal edge
    for g = 1:ng
        if fin(b, g) && fin(b + 1, g) && pos(b, g) ~= pos(b + 1, g)
            X = [X, xe(g), xe(g + 1), NaN]; Y = [Y, ye(b + 1), ye(b + 1), NaN]; %#ok<AGROW>
        end
    end
end
if isempty(X); return; end
plot(ax, X, Y, '--', 'Color', st.zeroline_color, 'LineWidth', 2, 'Tag', 'lle_zero');
end

function f = axis_transform(levels, scale)
% Where a level of the first axis sits along x: log10 on a log axis, the
% value itself on a linear one. 'auto' follows the grid
% (mc_gain_axis_scale): linear when the levels are evenly spaced or not
% all positive, log otherwise.
% A zero then log-spaced axis (a noise amplitude) is drawn as its log
% ladder under 'auto' (mc_axis_spacing).
if strcmp(scale, 'auto')
    [spacing, f] = mc_axis_spacing(levels);
    if strcmp(spacing, 'zero_log'); return; end
    scale = mc_gain_axis_scale(levels);
end
if strcmp(scale, 'log') && any(levels <= 0)
    error('fig_mc_gain_bias_maps:LogOfNonPositive', ...
        'first_scale ''log'' needs positive levels of the first axis; they run from %g.', min(levels));
end
if strcmp(scale, 'linear')
    f = @(v) v;
else
    f = @log10;
end
end

function s = capitalized(s)
% The text with its first letter in upper case, for an axis label.
if ~isempty(s); s(1) = upper(s(1)); end
end

function c = num_labels(v)
c = arrayfun(@(x) sprintf('%g', x), v, 'UniformOutput', false);
end

function [l_cmap, m_cmap] = row_colormaps(lle_name, mc_name)
% The lambda_1 row's and the memory-capacity row's colormaps, by name.
if strcmp(lle_name, 'blue_white_red')
    l_cmap = lle_cmap();
    m_cmap = mc_cmap();      % mc_name is 'white_orange_red': checked by the caller
    return
end
[l_cmap, to_center] = dark_center_colormaps(lle_name);
if strcmp(mc_name, 'lle_center')
    m_cmap = to_center;
else
    m_cmap = mc_cmap();
end
end

function s = on_off(tf)
if tf; s = 'on'; else; s = 'off'; end
end

function c = mc_cmap()
% Memory capacity: white (none) through orange to dark red.
c = interpolated_cmap([1 1 1; 0.99 0.88 0.55; 0.93 0.55 0.25; 0.70 0.18 0.18; 0.30 0.02 0.20]);
end

function c = lle_cmap()
% lambda_1: diverging, blue (stable) through white (zero) to red (chaotic).
c = interpolated_cmap([0.02 0.19 0.38; 0.26 0.58 0.76; 1 1 1; 0.84 0.38 0.30; 0.40 0.00 0.05]);
end

function c = spread_cmap()
% Spread over seeds: white to dark purple. Not a grey scale: grey is the
% colour of a cell with no data (the axes colour).
c = interpolated_cmap([1 1 1; 0.80 0.76 0.89; 0.55 0.45 0.72; 0.25 0.10 0.40]);
end

function c = interpolated_cmap(anchors)
% 256 colours, linear between equally spaced anchor rows.
c = interp1(linspace(0, 1, size(anchors, 1)), anchors, linspace(0, 1, 256));
end
