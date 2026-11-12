function fig = p3_summary_figure(T, visible)
% P3_SUMMARY_FIGURE Bars of the aberrant-activity rate and the mean E rate, per network, condition and ablation.
%
%   fig = P3_SUMMARY_FIGURE(T, visible)     T: the results table of run_p3_ablation
%
% One row of panels per network: left the aberrant events per minute
% (bursts for the bursting network, detector events for run13), right the
% mean E rate. Bars are coloured by regime (named beside the count); the
% dashed line is the intact chained MTS network's value. A saturated
% bursting network counts as ONE burst that never ends (0.6 per min).
%
% See also: run_p3_ablation

nets = unique(T.network, 'stable');
vis = 'off'; if visible, vis = 'on'; end
fig = figure('Color', 'w', 'Visible', vis, 'Units', 'pixels', 'Position', [40 40 1650 1150], ...
    'NumberTitle', 'off', 'Name', 'Pilot 3 summary');
tl = tiledlayout(fig, numel(nets), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
regimes = {'comparable', 'more', 'fewer', 'quiet', 'steady, no events', 'continuous high firing', 'saturation'};
colors = [0.20 0.60 0.30; 0.85 0.55 0.10; 0.45 0.65 0.85; 0.78 0.78 0.78; 0.40 0.40 0.45; 0.80 0.25 0.20; 0.45 0.05 0.05];
for i = 1:numel(nets)
    S = T(strcmp(T.network, nets{i}), :);
    labels = strcat(S.condition, {': '}, S.label);
    labels = strrep(labels, '_', ' ');
    ref = S(strcmp(S.condition, 'chain_mts') & strcmp(S.ablation, 'intact'), :);
    vals = {S.per_min, S.mean_rE};
    names = {'aberrant events per min', 'mean E rate'};
    for j = 1:2
        ax = nexttile(tl, (i - 1) * 2 + j); hold(ax, 'on');
        for k = 1:height(S)
            c = colors(strcmp(regimes, S.regime{k}), :);
            if isempty(c), c = [0.3 0.3 0.3]; end
            barh(ax, k, vals{j}(k), 'FaceColor', c, 'EdgeColor', 'none');
            note = '';
            if j == 1, note = [' (' S.regime{k} ')']; end
            text(ax, vals{j}(k), k, sprintf(' %.3g%s', vals{j}(k), note), 'FontSize', 7, 'VerticalAlignment', 'middle');
        end
        if j == 1, xline(ax, ref.per_min, '--k'); else, xline(ax, ref.mean_rE, '--k'); xlim(ax, [0 1.08]); end
        set(ax, 'YDir', 'reverse', 'YTick', 1:height(S), 'YTickLabel', labels, 'FontSize', 8, ...
            'Box', 'off', 'TickDir', 'out', 'TickLabelInterpreter', 'none');
        if j == 2, set(ax, 'YTickLabel', []); end
        ylim(ax, [0.4, height(S) + 0.6]);
        xlabel(ax, names{j});
        title(ax, sprintf('%s network: %s', nets{i}, names{j}), 'FontSize', 9, 'FontWeight', 'normal');
        if i == 1 && j == 2
            h = gobjects(1, numel(regimes));
            for q = 1:numel(regimes)
                h(q) = patch(ax, nan, nan, colors(q, :), 'EdgeColor', 'none');
            end
            legend(ax, h, regimes, 'Location', 'northeastoutside', 'Box', 'off', 'FontSize', 7);
        end
    end
end
title(tl, 'Pilot 3: unstimulated network with adaptation removed (dashed: intact chained MTS)', 'FontSize', 10);
end
