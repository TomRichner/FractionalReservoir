function fig = p2_frozen_figure(out, visible)
% P2_FROZEN_FIGURE Intact against frozen for one network of block B pilot 2.
%
% Columns: intact (the model's own run) and frozen from onset. Rows: the
% population mean rate; the mean adaptation term (c/K) sum a; the mean
% depression factor per route (SRNNCellTypePairs.std_route_factor); the
% mean input. Thin dashed lines are the same arm with NO stimulus. The
% vertical dotted line marks the freeze sample.
if nargin < 2, visible = 'off'; end
fig = figure('Visible', visible, 'Color', 'w', 'Units', 'pixels', 'Position', [60 60 1200 900]);
C = numel(out.cell_type_names);
colors = SRNNCellTypePairs.type_colors(C);
arms = {'intact', 'frozen'};
ctrl = {'nostim', 'nostim_frozen'};
titles = {'Adaptation intact', 'Adaptation frozen from onset'};
row_names = {'mean rate', 'mean adaptation term (c/K)\Sigma a', 'mean depression factor', 'mean input u'};
t = out.t;
tf = out.check.t_freeze;
ax = gobjects(4, 2);
tl = tiledlayout(fig, 4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for row = 1:4
    for col = 1:2
        A = out.arm.(arms{col});
        B = out.arm.(ctrl{col});
        ax(row, col) = nexttile(tl, (row - 1) * 2 + col);
        hold(ax(row, col), 'on');
        h = gobjects(0); lab = {};
        switch row
            case 1
                for q = C:-1:1
                    plot(t, B.rate(q, :), '--', 'Color', [colors(q, :) 0.55], 'LineWidth', 0.6);
                    h(q) = plot(t, A.rate(q, :), '-', 'Color', colors(q, :), 'LineWidth', 1.0);
                    lab{q} = sprintf('%s (mean after onset %.3f; no stimulus %.3f)', ...
                        out.cell_type_names{q}, A.rate_mean_after(q), B.rate_mean_after(q));
                end
            case 2
                for q = C:-1:1
                    if all(isnan(A.adapt(q, :))), continue; end
                    plot(t, B.adapt(q, :), '--', 'Color', [colors(q, :) 0.55], 'LineWidth', 0.6);
                    h(end + 1) = plot(t, A.adapt(q, :), '-', 'Color', colors(q, :), 'LineWidth', 1.0); %#ok<AGROW>
                    lab{end + 1} = out.cell_type_names{q}; %#ok<AGROW>
                end
            case 3
                rc = lines(size(A.dep, 1));
                for i = 1:size(A.dep, 1)
                    plot(t, B.dep(i, :), '--', 'Color', [rc(i, :) 0.55], 'LineWidth', 0.6);
                    h(i) = plot(t, A.dep(i, :), '-', 'Color', rc(i, :), 'LineWidth', 1.0);
                    lab{i} = A.dep_names{i};
                end
            case 4
                for q = C:-1:1
                    h(q) = plot(t, A.u(q, :), '-', 'Color', colors(q, :), 'LineWidth', 1.0);
                    lab{q} = sprintf('%s (mean over all %d; %d stimulated)', ...
                        out.cell_type_names{q}, out.n_per_type(q), out.stim_per_type(q));
                end
        end
        xline(ax(row, col), tf, 'k:', 'LineWidth', 1.0);
        box(ax(row, col), 'off'); set(ax(row, col), 'FontSize', 9);
        if row == 1, title(ax(row, col), titles{col}, 'FontWeight', 'normal', 'FontSize', 11); end
        if col == 1, ylabel(ax(row, col), row_names{row}, 'FontSize', 9); end
        if row == 4, xlabel(ax(row, col), 'time (s)'); end
        ok = isgraphics(h);
        if any(ok)
            legend(ax(row, col), h(ok), lab(ok), 'Location', 'best', 'Box', 'off', 'FontSize', 7);
        end
    end
    linkaxes(ax(row, :), 'xy');
end
linkaxes(ax(:), 'x'); xlim(ax(1, 1), [t(1) t(end)]);
title(tl, {out.label, 'Solid: with the stimulus. Dashed: the same arm with no stimulus. Dotted line: onset.'}, ...
    'FontSize', 10, 'Interpreter', 'none');
end
