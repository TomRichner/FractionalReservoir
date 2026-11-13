% EXPLORE_MC_VS_LAMBDA1_SCATTER Memory capacity against lambda_1, every run of the main gain sweep.
%
% Tom asked how correlated memory capacity is with lambda_1, and whether
% the relation is simply linear. The block C maps show the two side by
% side on the gain-by-bias grid; this puts them on one pair of axes. One
% panel per condition, every run a point (11 gains x 7 biases x 2 seeds),
% x = lambda_1 with the task input on (lle_input), y = total memory
% capacity, colored by log10 gain so the path along gain can be followed.
%
% x is linear and clipped at +-5 1/s, as in Tom's linear map: a clipped
% run is drawn at the axis edge as an open marker and counted in the
% panel. The statistics in each panel use the UNCLIPPED values of every
% run with both values finite: Pearson r and Spearman rho over all runs,
% and Pearson r on the stable side (lambda_1 < 0) and on the unstable side
% (lambda_1 > 0) alone, each with its n. If the relation is peaked near
% zero rather than monotone, the two one-sided r have opposite signs and
% the overall r understates it.
%
% Output: figs/explore_C/mc_vs_lambda1/Fig_MC_vs_lambda1_scatter.{png,svg,fig}
%
% See also: fig_mc_gain_bias_maps, read_mc_gain_sweep, explore_mc_maps_linear_scales

root = fileparts(which('setup_paths'));
D = read_mc_gain_sweep(fullfile(root, 'data', 'mc_gain_sweep_main', 'mc_gain_sweep', ...
    'mc_gain_sweep_results.mat'));
out_dir = fullfile(root, 'figs', 'explore_C', 'mc_vs_lambda1');
fig_tag = 'Fig_MC_vs_lambda1_scatter';

conds = select_conditions(D.conditions, {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2', 'sfa3_std2'});
titles = titles_by_name(conds);
colors = mc_condition_colors(conds);
st = manuscript_style();
lle_clip = 5;                                   % 1/s, as Tom's linear map
T = D.T;

ok_all = ~T.failed & isfinite(T.lle_input) & isfinite(T.mc_total);
mc_lim = [0, ceil(max(T.mc_total(ok_all)) + 0.5)];
lg = log10(D.gains);
cmap = parula(256);

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 40 11]);
tl = tiledlayout(fig, 1, numel(conds), 'TileSpacing', 'compact', 'Padding', 'compact');
for c = 1:numel(conds)
    rows = strcmp(T.condition, conds{c});
    n_rows = nnz(rows);
    ok = rows & ok_all;
    x = T.lle_input(ok); y = T.mc_total(ok); g = log10(T.gain(ok));
    clipped = abs(x) > lle_clip;
    xp = min(max(x, -lle_clip), lle_clip);

    ax = nexttile(tl); hold(ax, 'on'); box(ax, 'on');
    xline(ax, 0, '--', 'Color', st.zeroline_color, 'LineWidth', st.zeroline_lw);
    scatter(ax, xp(~clipped), y(~clipped), 22, g(~clipped), 'filled', ...
        'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', 'none');
    scatter(ax, xp(clipped), y(clipped), 28, g(clipped), 'LineWidth', 1);
    colormap(ax, cmap); clim(ax, [min(lg), max(lg)]);
    xlim(ax, lle_clip * [-1.05 1.05]); ylim(ax, mc_lim);
    xlabel(ax, '\lambda_1 with input (1/s)');
    if c == 1; ylabel(ax, 'Memory capacity'); end
    title(ax, titles{c}, 'Color', colors(c, :), 'FontWeight', 'bold');
    set(ax, 'FontSize', st.tick_fs - 2, 'TickDir', 'out', 'TitleFontSizeMultiplier', 1);

    s = x < 0; u = x > 0;
    txt = {sprintf('Pearson r = %.2f, n = %d', pearson(x, y), numel(x))
           sprintf('Spearman \\rho = %.2f', spearman(x, y))
           sprintf('\\lambda_1<0: r = %.2f, n = %d', pearson(x(s), y(s)), nnz(s))
           sprintf('\\lambda_1>0: r = %.2f, n = %d', pearson(x(u), y(u)), nnz(u))
           sprintf('clipped %d, failed %d', nnz(clipped), n_rows - nnz(ok))};
    text(ax, 0.56, 0.70, txt, 'Units', 'normalized', 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', 'FontSize', st.tick_fs - 5, ...
        'BackgroundColor', [1 1 1 0.8], 'Margin', 2);
end
cb = colorbar(ax);
cb.Layout.Tile = 'east';
gt = D.gains(unique(round(linspace(1, numel(D.gains), 6))));
cb.Ticks = log10(gt);
cb.TickLabels = arrayfun(@(v) sprintf('%.3g', v), gt, 'UniformOutput', false);
cb.Label.String = 'gain (log scale)';
title(tl, sprintf(['Memory capacity against \\lambda_1, every run (%d gains x %d biases x %d seeds); ' ...
    'open markers: |\\lambda_1| > %g, drawn at the edge'], numel(D.gains), numel(D.biases), ...
    numel(D.seeds), lle_clip), 'FontSize', st.title_fs - 2);

if ~exist(out_dir, 'dir'); mkdir(out_dir); end
save_figure_stable(out_dir, fig_tag, fig);
fprintf('%s\n', fullfile(out_dir, [fig_tag '.png']));

%% ------------------------------------------------------------------------
function r = pearson(x, y)
% Pearson correlation; NaN with fewer than three points.
if numel(x) < 3; r = NaN; return; end
x = x - mean(x); y = y - mean(y);
r = sum(x .* y) / sqrt(sum(x.^2) * sum(y.^2));
end

function rho = spearman(x, y)
% Spearman rank correlation: Pearson on average ranks (ties share a rank).
rho = pearson(avg_rank(x), avg_rank(y));
end

function r = avg_rank(v)
[s, i] = sort(v(:));
r = zeros(numel(v), 1);
k = 1;
while k <= numel(s)
    j = k;
    while j < numel(s) && s(j + 1) == s(k); j = j + 1; end
    r(i(k:j)) = (k + j) / 2;
    k = j + 1;
end
end
