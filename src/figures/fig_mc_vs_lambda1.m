function out = fig_mc_vs_lambda1(cfg)
% FIG_MC_VS_LAMBDA1 Memory capacity against lambda_1, every run of the gain-by-bias sweep, one panel per condition (block C).
%
%   out = FIG_MC_VS_LAMBDA1('run_dir', d)
%   out = FIG_MC_VS_LAMBDA1('data_file', f)
%   out = FIG_MC_VS_LAMBDA1('data', D, 'save', false)   % D from read_mc_gain_sweep, or a table
%
% Tom's question: how correlated is memory capacity with lambda_1, and is
% the relation simply linear? The maps (fig_mc_gain_bias_maps) show the two
% side by side on the gain-by-bias grid; this puts them on one pair of axes.
% One panel per condition, every run a point (every gain, bias and seed),
% x = lambda_1 (cfg.lle_field, default lle_input: with the memory task's
% input on), y = total memory capacity, colored by gain on a log scale
% (cfg.color_by = 'gain', the default) so the path along gain can be
% followed; cfg.color_by = 'bias' colors by the bias S_c instead.
%
% x IS LINEAR AND CLIPPED at +-cfg.lle_clip (default 5 1/s, as the linear
% lambda_1 map): a run beyond it is drawn at the axis edge as an open
% marker and counted in the panel. The y range, 0 to just above the
% largest memory capacity of the conditions drawn, is shared by every
% panel. The lambda_1 = 0 line is green dashed, as everywhere in the
% manuscript.
%
% STATISTICS, in each panel and in out.stats (one row per condition), use
% the UNCLIPPED values of every run whose lambda_1 and memory capacity are
% both finite (failed runs left out): Pearson r and Spearman rho over all
% runs, and Pearson r on the stable side (lambda_1 < 0) and on the unstable
% side (lambda_1 > 0) alone, each with its n; then the clipped and failed
% counts. If the relation is peaked near zero rather than monotone, the two
% one-sided r have opposite signs and the overall r understates it. Pearson
% and Spearman are the Statistics and Machine Learning Toolbox's corr
% (Spearman with tied ranks averaged); a correlation over fewer than three
% points is NaN.
%
% Conditions by cfg.conditions (select_conditions; empty = all the result
% holds, in its order), titles and colors by name (titles_by_name,
% mc_condition_colors). The result is taken through resolve_mc_gain_sweep
% (data, data_file or run_dir).
%
% out: figs, files, source, and stats (a table: condition, n, pearson_r,
% spearman_rho, n_stable, r_stable, n_unstable, r_unstable, n_clipped,
% n_failed).
%
% See also: fig_mc_gain_bias_maps, fig_mc_gain_curves, fig_mc_robustness,
%           read_mc_gain_sweep, explore_mc_vs_lambda1_scatter

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.run_dir     (1,:) char    = ''
    cfg.data_file   (1,:) char    = ''
    cfg.data                      = []
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.preset_name (1,:) char    = ''   % unused; accepted for a uniform call
    cfg.conditions                = {}
    cfg.lle_clip    (1,1) double  {mustBePositive} = 5   % |lambda_1| clipped here (1/s)
    cfg.lle_field   (1,:) char    {mustBeMember(cfg.lle_field, {'lle_input', 'lle_noinput'})} = 'lle_input'
    cfg.color_by    (1,:) char    {mustBeMember(cfg.color_by, {'gain', 'bias'})} = 'gain'
    cfg.fig_tag     (1,:) char    = 'Fig_MC_vs_lambda1'
end

setup_paths();
D = resolve_mc_gain_sweep('data', cfg.data, 'data_file', cfg.data_file, 'run_dir', cfg.run_dir);
conds = select_conditions(D.conditions, cfg.conditions);
nc = numel(conds);
st = manuscript_style();
cond_titles = titles_by_name(conds);
cond_colors = mc_condition_colors(conds);
T = D.T;
clip = cfg.lle_clip;

MARKER_SIZE      = 22;     % filled marker, an unclipped run
MARKER_SIZE_OPEN = 28;     % open marker, a clipped run at the axis edge
MARKER_ALPHA     = 0.85;
STATS_POSITION   = [0.56 0.70];   % normalized, top left corner of the statistics box

ok_all = ~T.failed & isfinite(T.(cfg.lle_field)) & isfinite(T.mc_total);
in_set = ismember(T.condition, conds);
mc_lim = [0, ceil(max([T.mc_total(ok_all & in_set); 0]) + 0.5)];

% The color variable and its scale
switch cfg.color_by
    case 'gain'
        cvar = log10(T.gain);
        levels = D.gains;
        to_scale = @log10;
        cb_label = 'gain (log scale)';
    case 'bias'
        cvar = T.bias;
        levels = D.biases;
        to_scale = @(v) v;
        cb_label = 'bias S_c';
end
c_lim = [min(to_scale(levels)), max(to_scale(levels))];
if diff(c_lim) <= 0; c_lim = c_lim(1) + [-0.5 0.5]; end
cmap = parula(256);

if strcmp(cfg.lle_field, 'lle_input')
    x_label = '\lambda_1 with input (1/s)';
else
    x_label = '\lambda_1 without input (1/s)';
end

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 10 * nc 11]);
if ~cfg.visible; set(fig, 'Visible', 'off'); end
tl = tiledlayout(fig, 1, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
stats = cell(nc, 1);
for c = 1:nc
    rows = strcmp(T.condition, conds{c});
    ok = rows & ok_all;
    x = T.(cfg.lle_field)(ok);
    y = T.mc_total(ok);
    g = cvar(ok);
    clipped = abs(x) > clip;
    xp = min(max(x, -clip), clip);

    ax = nexttile(tl); hold(ax, 'on'); box(ax, 'on');
    xline(ax, 0, '--', 'Color', st.zeroline_color, 'LineWidth', st.zeroline_lw);
    scatter(ax, xp(~clipped), y(~clipped), MARKER_SIZE, g(~clipped), 'filled', ...
        'MarkerFaceAlpha', MARKER_ALPHA, 'MarkerEdgeColor', 'none');
    scatter(ax, xp(clipped), y(clipped), MARKER_SIZE_OPEN, g(clipped), 'LineWidth', 1);
    colormap(ax, cmap); clim(ax, c_lim);
    xlim(ax, clip * [-1.05 1.05]); ylim(ax, mc_lim);
    xlabel(ax, x_label);
    if c == 1; ylabel(ax, 'Memory capacity'); end
    title(ax, cond_titles{c}, 'Color', cond_colors(c, :), 'FontWeight', 'bold');
    set(ax, 'FontSize', st.tick_fs - 2, 'TickDir', 'out', 'TitleFontSizeMultiplier', 1);

    s = x < 0; u = x > 0;
    S = struct('condition', conds(c), 'n', numel(x), 'pearson_r', pearson(x, y), ...
        'spearman_rho', spearman(x, y), 'n_stable', nnz(s), 'r_stable', pearson(x(s), y(s)), ...
        'n_unstable', nnz(u), 'r_unstable', pearson(x(u), y(u)), ...
        'n_clipped', nnz(clipped), 'n_failed', nnz(rows) - nnz(ok));
    stats{c} = S;
    txt = {sprintf('Pearson r = %.2f, n = %d', S.pearson_r, S.n)
           sprintf('Spearman \\rho = %.2f', S.spearman_rho)
           sprintf('\\lambda_1<0: r = %.2f, n = %d', S.r_stable, S.n_stable)
           sprintf('\\lambda_1>0: r = %.2f, n = %d', S.r_unstable, S.n_unstable)
           sprintf('clipped %d, failed %d', S.n_clipped, S.n_failed)};
    text(ax, STATS_POSITION(1), STATS_POSITION(2), txt, 'Units', 'normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', st.tick_fs - 5, ...
        'BackgroundColor', [1 1 1 0.8], 'Margin', 2, 'Tag', 'mc_lambda1_stats');
end
cb = colorbar(ax);
cb.Layout.Tile = 'east';
tick_levels = levels(unique(round(linspace(1, numel(levels), min(6, numel(levels))))));
cb.Ticks = to_scale(tick_levels);
cb.TickLabels = arrayfun(@(v) sprintf('%.3g', v), tick_levels, 'UniformOutput', false);
cb.Label.String = cb_label;
title(tl, sprintf(['Memory capacity against \\lambda_1, every run (%d gains x %d biases x %d seeds); ' ...
    'open markers: |\\lambda_1| > %g, drawn at the edge'], numel(D.gains), numel(D.biases), ...
    numel(D.seeds), clip), 'FontSize', st.title_fs - 2);

out = struct('figs', fig, 'files', {{}}, 'source', D.source, ...
    'stats', struct2table(vertcat(stats{:}), 'AsArray', true));
if cfg.save
    out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
    save_figure_stable(out_dir, cfg.fig_tag, fig);
    out.files = existing_outputs(out_dir, cfg.fig_tag);
end
end

%% ------------------------------------------------------------------------
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
