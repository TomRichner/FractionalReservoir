function fig = p1b_summary_figure(R, stims, conds, cond_labels, S, visible)
% P1B_SUMMARY_FIGURE Before against during, per stimulus case and condition.
%
%   fig = P1B_SUMMARY_FIGURE(R, stims, conds, cond_labels, S, visible)
%
% R.(cond){k} is the p1b_run_case output of condition cond under stims(k),
% or [] where that pair was not run. Six panels: bursts per s, mean E
% rate, largest population mean E rate, its peak-to-peak range, mean
% adaptation term (E), mean depression factor. Each condition
% is one colour, a small x offset apart; open marker = before (10-40 s),
% filled = during (50-100 s), joined by a line.

vis = 'off'; if visible, vis = 'on'; end
fig = figure('Color', 'w', 'Visible', vis, 'Units', 'pixels', ...
    'Position', [60 60 1300 760], 'NumberTitle', 'off', 'Name', 'pilot 1b summary');
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
cols = [0.80 0.40 0.00; 0.00 0.45 0.70; 0.00 0.60 0.50];
fields = {'burst_rate', 'mean_r_E', 'max_mean_r_E', 'p2p_mean_r_E', 'mean_adaptation', 'mean_depression'};
ylabs  = {'bursts per s', 'mean E rate', 'largest population mean E rate', ...
          'peak-to-peak of population mean E rate', 'mean (c/K)\Sigma_k a_k, E', 'mean depression factor'};
nK = numel(stims);
offs = linspace(-0.18, 0.18, numel(conds));
for f = 1:numel(fields)
    ax = nexttile(tl); hold(ax, 'on');
    h = gobjects(1, numel(conds));
    for c = 1:numel(conds)
        for k = 1:nK
            C = R.(conds{c}){k};
            if isempty(C), continue; end
            x = k + offs(c);
            yb = C.measures.before.(fields{f});
            yd = C.measures.during.(fields{f});
            plot(ax, [x x], [yb yd], '-', 'Color', cols(c, :), 'HandleVisibility', 'off');
            plot(ax, x, yb, 'o', 'Color', cols(c, :), 'MarkerSize', 7, 'HandleVisibility', 'off');
            h(c) = plot(ax, x, yd, 'o', 'Color', cols(c, :), 'MarkerFaceColor', cols(c, :), ...
                'MarkerSize', 7, 'DisplayName', cond_labels{c});
        end
    end
    hold(ax, 'off');
    set(ax, 'XTick', 1:nK, 'XTickLabel', {stims.short}, 'Box', 'off', 'TickDir', 'out');
    xlim(ax, [0.5 nK + 0.5]);
    ylabel(ax, ylabs{f});
    if f == 1
        ok = isgraphics(h);
        legend(ax, h(ok), 'Location', 'north', 'Box', 'off');
    end
end
title(tl, sprintf(['bursting_pairs, seed %s. Open: before stimulus (%g-%g s); filled: ' ...
    'during (%g-%g s); onset at %g s. Bursts exclude above-threshold stretches longer than %g s.'], ...
    mat2str(S.rng_seeds), S.win_before, S.win_during, S.t_stim, S.burst.sustained_max), ...
    'Interpreter', 'none', 'FontSize', 10);
end
