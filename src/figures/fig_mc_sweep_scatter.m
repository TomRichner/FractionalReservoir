function out = fig_mc_sweep_scatter(cfg)
% FIG_MC_SWEEP_SCATTER One quantity of the gain sweep of memory capacity against another, every run, one panel per condition (block C).
%
%   out = FIG_MC_SWEEP_SCATTER('run_dir', d, 'x', 'lle_input', 'y', 'mc_total')
%   out = FIG_MC_SWEEP_SCATTER('data_file', f, 'x', 'bias', 'y', 'rate_mean', 'color_by', 'gain')
%   out = FIG_MC_SWEEP_SCATTER('data', D, 'save', false)   % D from read_mc_gain_sweep, or a table
%
% Every run of run_mc_gain_sweep (every gain, level of the second axis and
% seed) is a point, x against y, colored by a third quantity, one panel per
% condition. The quantities are named, and defined once, by
% mc_sweep_quantity: the table's columns (gain, the result's second axis
% under its own name, bias or self_feedback, then mc_total, horizon,
% lle_input, lle_noinput, rate_E, rate_I, frac_pinned, frac_silent) and the
% derived rate_mean and
% frac_on_slope. Each brings its label, its axis scale (log for gain), its
% axis limits and an optional clip. fig_mc_vs_lambda1 is this figure with
% x = lle_input, y = mc_total, colored by gain.
%
% CLIPPED QUANTITIES (lambda_1 at +-5 1/s): a run beyond the clip on
% either axis is drawn at the axis edge as an open marker and counted in
% the panel; cfg.clip (a struct, one field per quantity name, e.g.
% struct('lle_input', 3)) replaces a quantity's clip, or clips one that
% has none. Axis limits are shared by every panel. Where x or y is
% lambda_1, its zero is a green dashed line, as everywhere in the
% manuscript.
%
% STATISTICS, in each panel and in out.stats (one row per condition), use
% the UNCLIPPED values of every run whose x and y are both finite (failed
% runs left out): Pearson r and Spearman rho over all runs; when x is
% lambda_1 (a quantity whose sign splits the regime), also Pearson r on the
% stable side (x < 0) and on the unstable side (x > 0) alone, each with its
% n; then the clipped and failed counts. Pearson and Spearman are the
% Statistics and Machine Learning Toolbox's corr (Spearman with tied ranks
% averaged, so a grid axis with few levels, like bias, gives many ties); a
% correlation over fewer than three points is NaN.
%
% OVERLAYS: cfg.overlay is a struct array of reference curves, each with
% fields x, y (vectors in the quantities' units), label, and conditions
% (cellstr of the panels it is drawn on; empty = every panel). A curve is
% drawn as a black line and named in that panel's legend.
%
% Conditions by cfg.conditions (select_conditions; empty = all the result
% holds, in its order), titles and colors by name (titles_by_name,
% mc_condition_colors). The result is taken through resolve_mc_gain_sweep
% (data, data_file or run_dir). The file tag defaults to
% Fig_MC_sweep_<y>_vs_<x>.
%
% out: figs, files, source, and stats (a table: condition, n, pearson_r,
% spearman_rho, [n_stable, r_stable, n_unstable, r_unstable when x is split
% at zero], n_clipped, n_failed).
%
% See also: mc_sweep_quantity, fig_mc_vs_lambda1, fig_mc_gain_bias_maps,
%           read_mc_gain_sweep, explore_bias_rate_lambda1

arguments
    cfg.verbose        (1,:) char    = 'minimal'
    cfg.run_dir        (1,:) char    = ''
    cfg.data_file      (1,:) char    = ''
    cfg.data                         = []
    cfg.out_dir        (1,:) char    = ''
    cfg.save           (1,1) logical = true
    cfg.visible        (1,1) logical = true
    cfg.preset_name    (1,:) char    = ''   % unused; accepted for a uniform call
    cfg.conditions                   = {}
    cfg.x              (1,:) char    = 'lle_input'
    cfg.y              (1,:) char    = 'mc_total'
    cfg.color_by       (1,:) char    = 'gain'
    cfg.clip           (1,1) struct  = struct()   % per quantity name, a clip replacing the registry's
    cfg.overlay        struct        = struct('x', {}, 'y', {}, 'label', {}, 'conditions', {})
    cfg.stats_position (1,2) double  = [0.56 0.70]   % normalized, top left corner of the statistics box
    cfg.fig_tag        (1,:) char    = ''
end

setup_paths();
D = resolve_mc_gain_sweep('data', cfg.data, 'data_file', cfg.data_file, 'run_dir', cfg.run_dir);
conds = select_conditions(D.conditions, cfg.conditions);
nc = numel(conds);
st = manuscript_style();
cond_titles = titles_by_name(conds);
cond_colors = mc_condition_colors(conds);
T = D.T;
quantity = @(name) mc_sweep_quantity(name, D, 'clip', clip_override(cfg.clip, name));
Qx = quantity(cfg.x);
Qy = quantity(cfg.y);
Qc = quantity(cfg.color_by);
fig_tag = cfg.fig_tag;
if isempty(fig_tag); fig_tag = sprintf('Fig_MC_sweep_%s_vs_%s', cfg.y, cfg.x); end

MARKER_SIZE      = 22;     % filled marker, an unclipped run
MARKER_SIZE_OPEN = 28;     % open marker, a clipped run at the axis edge
MARKER_ALPHA     = 0.85;
MAX_LEVEL_TICKS  = 8;      % a grid axis with more levels than this gets six colorbar ticks

ok_all = ~T.failed & isfinite(Qx.values) & isfinite(Qy.values);
in_set = ismember(T.condition, conds);
x_lim = Qx.limits(Qx.values(ok_all & in_set));
y_lim = Qy.limits(Qy.values(ok_all & in_set));

% The color variable and its scale: a grid axis takes its levels as ticks,
% any other quantity the range of the runs drawn.
to_scale = @(v) v;
if strcmp(Qc.scale, 'log'); to_scale = @log10; end
cvar = to_scale(Qc.values);
levels = Qc.levels;
if isempty(levels)
    finite_c = cvar(ok_all & in_set & isfinite(cvar));
    c_lim = [0 1];
    if ~isempty(finite_c); c_lim = [min(finite_c), max(finite_c)]; end
else
    c_lim = [min(to_scale(levels)), max(to_scale(levels))];
end
if diff(c_lim) <= 0; c_lim = c_lim(1) + [-0.5 0.5]; end
cmap = parula(256);
cb_label = Qc.label;
if strcmp(Qc.scale, 'log'); cb_label = [cb_label ' (log scale)']; end

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 10 * nc 11]);
if ~cfg.visible; set(fig, 'Visible', 'off'); end
tl = tiledlayout(fig, 1, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
stats = cell(nc, 1);
for c = 1:nc
    rows = strcmp(T.condition, conds{c});
    ok = rows & ok_all;
    x = Qx.values(ok);
    y = Qy.values(ok);
    g = cvar(ok);
    [xp, clip_x] = clip_to(x, Qx.clip);
    [yp, clip_y] = clip_to(y, Qy.clip);
    clipped = clip_x | clip_y;

    ax = nexttile(tl); hold(ax, 'on'); box(ax, 'on');
    if Qx.zero_line; xline(ax, 0, '--', 'Color', st.zeroline_color, 'LineWidth', st.zeroline_lw); end
    if Qy.zero_line; yline(ax, 0, '--', 'Color', st.zeroline_color, 'LineWidth', st.zeroline_lw); end
    scatter(ax, xp(~clipped), yp(~clipped), MARKER_SIZE, g(~clipped), 'filled', ...
        'MarkerFaceAlpha', MARKER_ALPHA, 'MarkerEdgeColor', 'none');
    scatter(ax, xp(clipped), yp(clipped), MARKER_SIZE_OPEN, g(clipped), 'LineWidth', 1);
    draw_overlays(ax, cfg.overlay, conds{c});
    colormap(ax, cmap); clim(ax, c_lim);
    set(ax, 'XScale', Qx.scale, 'YScale', Qy.scale);
    xlim(ax, x_lim); ylim(ax, y_lim);
    xlabel(ax, Qx.label);
    if c == 1; ylabel(ax, Qy.label); end
    title(ax, cond_titles{c}, 'Color', cond_colors(c, :), 'FontWeight', 'bold');
    set(ax, 'FontSize', st.tick_fs - 2, 'TickDir', 'out', 'TitleFontSizeMultiplier', 1);

    S = struct('condition', conds(c), 'n', numel(x), 'pearson_r', pearson(x, y), ...
        'spearman_rho', spearman(x, y));
    txt = {sprintf('Pearson r = %.2f, n = %d', S.pearson_r, S.n)
           sprintf('Spearman \\rho = %.2f', S.spearman_rho)};
    if Qx.zero_split
        s = x < 0; u = x > 0;
        S.n_stable = nnz(s);   S.r_stable = pearson(x(s), y(s));
        S.n_unstable = nnz(u); S.r_unstable = pearson(x(u), y(u));
        txt = [txt
               sprintf('%s<0: r = %.2f, n = %d', Qx.symbol, S.r_stable, S.n_stable)
               sprintf('%s>0: r = %.2f, n = %d', Qx.symbol, S.r_unstable, S.n_unstable)]; %#ok<AGROW>
    end
    S.n_clipped = nnz(clipped);
    S.n_failed = nnz(rows) - nnz(ok);
    stats{c} = S;
    txt{end + 1} = sprintf('clipped %d, failed %d', S.n_clipped, S.n_failed); %#ok<AGROW>
    text(ax, cfg.stats_position(1), cfg.stats_position(2), txt, 'Units', 'normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', st.tick_fs - 5, ...
        'BackgroundColor', [1 1 1 0.8], 'Margin', 2, 'Tag', 'mc_sweep_scatter_stats');
end
cb = colorbar(ax);
cb.Layout.Tile = 'east';
if ~isempty(levels)
    tick_levels = levels;   % a short grid (bias: 7 levels) shows every level
    if numel(levels) > MAX_LEVEL_TICKS
        tick_levels = levels(unique(round(linspace(1, numel(levels), 6))));
    end
    cb.Ticks = to_scale(tick_levels);
    cb.TickLabels = arrayfun(@(v) sprintf('%.3g', v), tick_levels, 'UniformOutput', false);
elseif strcmp(Qc.scale, 'log')
    cb.TickLabels = arrayfun(@(v) sprintf('%.3g', 10^v), cb.Ticks, 'UniformOutput', false);
end
cb.Label.String = cb_label;

Q2 = quantity(D.second_axis.column);   % the sweep's second axis: what its levels are called
caption = sprintf('%s against %s, every run (%d gains x %d %s x %d seeds)', Qy.symbol, ...
    Qx.symbol, numel(D.gains), numel(D.second_levels), Q2.plural, numel(D.seeds));
clipped_q = [Qx, Qy];
clipped_q = clipped_q(~cellfun(@isempty, {clipped_q.clip}));
if ~isempty(clipped_q)
    [~, first] = unique({clipped_q.symbol}, 'stable');
    parts = arrayfun(@(q) sprintf('|%s| > %g', q.symbol, q.clip), clipped_q(first), 'UniformOutput', false);
    caption = sprintf('%s; open markers: %s, drawn at the edge', caption, strjoin(parts, ' or '));
end
title(tl, caption, 'FontSize', st.title_fs - 2);

out = struct('figs', fig, 'files', {{}}, 'source', D.source, ...
    'stats', struct2table(vertcat(stats{:}), 'AsArray', true));
if cfg.save
    out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
end
end

%% ------------------------------------------------------------------------
function b = clip_override(clip, name)
% The caller's clip for this quantity, or [] for the registry's own.
b = [];
if isfield(clip, name); b = clip.(name); end
end

function [vp, clipped] = clip_to(v, bound)
% Values beyond +-bound drawn at the bound; an empty bound clips nothing.
if isempty(bound)
    vp = v; clipped = false(size(v));
    return
end
clipped = abs(v) > bound;
vp = min(max(v, -bound), bound);
end

function draw_overlays(ax, overlay, condition)
% The reference curves meant for this panel, black, named in a legend.
h = gobjects(0);
for k = 1:numel(overlay)
    o = overlay(k);
    if ~isempty(o.conditions) && ~any(strcmp(o.conditions, condition)); continue; end
    h(end + 1) = plot(ax, o.x, o.y, '-', 'Color', 'k', 'LineWidth', 1.2, ...
        'DisplayName', o.label, 'Tag', 'mc_sweep_scatter_overlay'); %#ok<AGROW>
end
if ~isempty(h)
    legend(ax, h, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
end
end

function r = pearson(x, y)
% Pearson correlation (corr); NaN with fewer than three points.
if numel(x) < 3; r = NaN; return; end
r = corr(x(:), y(:));
end

function rho = spearman(x, y)
% Spearman rank correlation (corr, tied ranks averaged); NaN with fewer
% than three points.
if numel(x) < 3; rho = NaN; return; end
rho = corr(x(:), y(:), 'Type', 'Spearman');
end
