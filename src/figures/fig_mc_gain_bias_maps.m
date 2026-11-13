function out = fig_mc_gain_bias_maps(cfg)
% FIG_MC_GAIN_BIAS_MAPS Memory capacity and lambda_1 on the gain-by-bias grid, one column per condition (block C).
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
%   3. the interquartile range of memory capacity over seeds, shared scale,
%      only when cfg.show_iqr (default true) AND some cell has more than one
%      successful seed; with one seed the row would be blank, so it is
%      dropped and the caption says "one seed".
% x is the gain (level_of_chaos) on a log axis, y the bias S_c; the ticks
% are the grid's own levels. Each cell is one grid point, drawn between the
% midpoints to its neighbours, so an irregular grid is drawn as it is.
%
% THE lambda_1 = 0 LINE is drawn on every row as the cell boundaries that
% separate a cell with lambda_1 > 0 from a neighbour with lambda_1 <= 0
% (the median map of the same condition), green dashed as everywhere in
% the manuscript. A staircase on the cell edges, not an interpolated
% contour: on a coarse grid an interpolated contour invents a position
% between levels and stops at the outermost cell centres.
%
% The preset's own point (D.preset_gain, D.preset_bias: gain 1, bias 0.2)
% is a white star, drawn where it is even when it is not a grid level. A
% cell with no successful seed is left grey (the axes colour).
%
% SATURATED AND SILENT CELLS: a cell whose median over seeds of the
% fraction of neurons pinned above 0.9 is >= cfg.flag_fraction (0.9)
% carries a small grey '+', one whose fraction silent below 0.01 is >=
% cfg.flag_fraction a small grey 'o', so a low memory capacity there is
% not read as chaos. 0.9 rather than 0.5: on the pilot the chained MTS at
% bias 0.55 has 57 % of its neurons silent and still a memory capacity of
% 5 to 9, so half silent is not a failed network; at bias 0.8 every
% neuron is silent. Flags are per condition.
%
% WHICH CONDITIONS: cfg.conditions through select_conditions (empty = all the
% result holds, in its order); titles and colours by name from
% manuscript_style. Works for three and four columns.
%
% The table is read by read_mc_gain_sweep (the one reader), located by
% resolve_mc_gain_sweep.
%
% See also: fig_mc_gain_curves, fig_mc_robustness, mc_gain_bias_statistics, symlog_lle

arguments
    cfg.verbose       (1,:) char    = 'minimal'
    cfg.run_dir       (1,:) char    = ''
    cfg.data_file     (1,:) char    = ''
    cfg.data                        = []
    cfg.out_dir       (1,:) char    = ''
    cfg.save          (1,1) logical = true
    cfg.visible       (1,1) logical = true
    cfg.preset_name   (1,:) char    = ''   % unused; accepted for a uniform call
    cfg.conditions                  = {}
    cfg.show_iqr      (1,1) logical = true
    cfg.flag_fraction (1,1) double  = 0.9
    cfg.mc_scale      (1,:) char    {mustBeMember(cfg.mc_scale, {'sqrt', 'linear'})} = 'sqrt'
    cfg.lle_lin       (1,1) double  = 0.05 % symlog_lle's linear scale (1/s)
    cfg.lle_clip      (1,1) double  = 10   % |lambda_1| clipped here (1/s)
    cfg.fig_tag       (1,:) char    = 'Fig_MC_gain_bias_maps'
end

setup_paths();
D = load_data(cfg);
[conds, ~] = select_conditions(D.conditions, cfg.conditions);
nc = numel(conds);
st = manuscript_style();

%% The grids
nb = numel(D.biases); ng = numel(D.gains);
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
% not white. White would read as "the seeds agree", and on a run in
% progress (seeds finishing unevenly) that was most of the row.
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
[~, l_tp, l_tl, l_lim] = symlog_lle(0, cfg.lle_lin, cfg.lle_clip);
L_col = symlog_lle(L, cfg.lle_lin, cfg.lle_clip);
iq_clim = [0, max([IQ(:); eps])];

rows = {struct('V', mc_f(MC), 'clim', mc_f([0 mc_max]), 'cmap', seq_cmap(), ...
               'label', 'Memory capacity', 'cb', mc_cb, ...
               'ticks', mc_f(mc_ticks), 'ticklabels', {num_labels(mc_ticks)})
        struct('V', L_col, 'clim', l_lim, 'cmap', div_cmap(), ...
               'label', '\lambda_1, task input on', ...
               'cb', sprintf('\\lambda_1 (1/s), symmetric log, clipped at \\pm%g', cfg.lle_clip), ...
               'ticks', l_tp, 'ticklabels', {l_tl})};
if cfg.show_iqr && n_seeds > 1
    rows{end+1} = struct('V', IQ, 'clim', iq_clim, 'cmap', gray_cmap(), ...
                         'label', 'Spread over seeds', 'cb', 'MC IQR over seeds', ...
                         'ticks', [], 'ticklabels', {{}});
end
nr = numel(rows);

%% Draw
fig = figure('Color', 'w', 'Position', [60 60 270 * nc + 200, 260 * nr + 120]);
if ~cfg.visible; set(fig, 'Visible', 'off'); end
outer = tiledlayout(fig, nr, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
x = log10(D.gains); y = D.biases;
xe = edges(x); ye = edges(y);
xtl = arrayfun(@(v) sprintf('%.3g', v), D.gains, 'UniformOutput', false);
ytl = arrayfun(@(v) sprintf('%.3g', v), D.biases, 'UniformOutput', false);
ax_all = gobjects(nr, nc);
% A qualifier in parentheses goes to a second line, so that four long
% titles ("... Adaptation (product)") do not run together. When any title
% wraps, every title gets two lines (a blank second line where there is no
% qualifier), so the condition names stay on one baseline across the row.
titles = cell(1, nc);
for c = 1:nc
    titles{c} = strsplit(regexprep(cond_title(st, conds{c}), '\s+(\(.*\))$', '\n$1'), newline);
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
        set(ax, 'Color', [0.85 0.85 0.85], 'Layer', 'top', 'TickDir', 'out', ...
            'FontSize', st.tick_fs - 5, 'XTick', x, 'XTickLabel', xtl, ...
            'YTick', y, 'YTickLabel', ytl, 'Box', 'on', 'XTickLabelRotation', 45);
        xlim(ax, xe([1 end])); ylim(ax, ye([1 end]));
        mark_flags(ax, x, y, PIN(:, :, c) >= cfg.flag_fraction, SIL(:, :, c) >= cfg.flag_fraction);
        zero_boundary(ax, xe, ye, L(:, :, c), st);
        plot(ax, log10(D.preset_gain), D.preset_bias, 'p', 'MarkerSize', 11, ...
            'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
        if r == 1
            title(ax, titles{c}, 'FontSize', st.title_fs - 2, ...
                'FontWeight', 'normal', 'Color', cond_color(st, conds{c}));
        end
        if c == 1
            ylabel(ax, {R.label, 'Bias S_c'}, 'FontSize', st.label_fs - 3);
        else
            set(ax, 'YTickLabel', []);
        end
        if r == nr
            xlabel(ax, 'Gain', 'FontSize', st.label_fs - 3);
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
grey_txt = 'Grey: no successful job';
if nr == 3; grey_txt = [grey_txt ' (spread row: fewer than two)']; end
title(outer, {sprintf('%s. Green dashed: edge between cells with \\lambda_1 > 0 and \\lambda_1 \\leq 0.', seed_txt), ...
    sprintf('Star: preset (gain %.3g, bias %.3g). + saturated, o silent (\\geq %.0f%% of neurons). %s.', ...
        D.preset_gain, D.preset_bias, 100 * cfg.flag_fraction, grey_txt)}, ...
    'FontSize', st.tick_fs - 5);

%% Save
out = struct('figs', fig, 'files', {{}}, 'source', D.source);
if cfg.save
    out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
    save_figure_stable(out_dir, cfg.fig_tag, fig);
    out.files = existing_outputs(out_dir, cfg.fig_tag);
end
end

%% ------------------------------------------------------------------------
function D = load_data(cfg)
if ~isempty(cfg.data)
    if isstruct(cfg.data) && isfield(cfg.data, 'gains'); D = cfg.data; else; D = read_mc_gain_sweep(cfg.data); end
else
    D = resolve_mc_gain_sweep(cfg.data_file, cfg.run_dir);
end
end

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
    for g = 1:ng - 1       % neighbours along gain: a vertical edge
        if fin(b, g) && fin(b, g + 1) && pos(b, g) ~= pos(b, g + 1)
            X = [X, xe(g + 1), xe(g + 1), NaN]; Y = [Y, ye(b), ye(b + 1), NaN]; %#ok<AGROW>
        end
    end
end
for b = 1:nb - 1           % neighbours along bias: a horizontal edge
    for g = 1:ng
        if fin(b, g) && fin(b + 1, g) && pos(b, g) ~= pos(b + 1, g)
            X = [X, xe(g), xe(g + 1), NaN]; Y = [Y, ye(b + 1), ye(b + 1), NaN]; %#ok<AGROW>
        end
    end
end
if isempty(X); return; end
plot(ax, X, Y, '--', 'Color', st.zeroline_color, 'LineWidth', 2, 'Tag', 'lle_zero');
end

function c = num_labels(v)
c = arrayfun(@(x) sprintf('%g', x), v, 'UniformOutput', false);
end

function c = seq_cmap()
k = [1 1 1; 0.99 0.88 0.55; 0.93 0.55 0.25; 0.70 0.18 0.18; 0.30 0.02 0.20];
c = interp1(linspace(0, 1, size(k, 1)), k, linspace(0, 1, 256));
end

function c = div_cmap()
k = [0.02 0.19 0.38; 0.26 0.58 0.76; 1 1 1; 0.84 0.38 0.30; 0.40 0.00 0.05];
c = interp1(linspace(0, 1, size(k, 1)), k, linspace(0, 1, 256));
end

function c = gray_cmap()
% White to dark purple, not grey: grey is the colour of a cell with no data
% (the axes colour), and a grey scale made the two indistinguishable.
k = [1 1 1; 0.80 0.76 0.89; 0.55 0.45 0.72; 0.25 0.10 0.40];
c = interp1(linspace(0, 1, size(k, 1)), k, linspace(0, 1, 256));
end

function t = cond_title(st, nm)
if isKey(st.condition_title, nm); t = st.condition_title(nm); else; t = nm; end
end

function c = cond_color(st, nm)
if isKey(st.condition_color, nm); c = st.condition_color(nm); else; c = [0.3 0.3 0.3]; end
end
