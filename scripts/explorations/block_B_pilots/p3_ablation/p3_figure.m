function fig = p3_figure(runs, titles, fig_title, visible)
% P3_FIGURE One column per run: the rate of every neuron, and the population mean rates with the events marked.
%
%   fig = P3_FIGURE(runs, titles, fig_title, visible)
%
% runs     cell of p3_run_one results (one network), one per column
% titles   cell of column titles (first line; the counts are added)
% Row 1: rate of every neuron against time, in 0.05-s bins, black 0 to
% white 1. Bursting network: E above the blue line, I below. run13: by
% module (white lines), E before I within a module.
% Row 2: population mean rate of E and I; bursts shaded (bursting network)
% or events as triangles (run13, with the mean E rate of each of the five
% modules in grey behind); the burst threshold 0.2 dotted.
% Each row shares its y axis across the columns.
%
% See also: run_p3_ablation, p3_run_one

nc = numel(runs);
cE = [0.85 0.33 0.10]; cI = [0.00 0.45 0.74];
vis = 'off'; if visible, vis = 'on'; end
fig = figure('Color', 'w', 'Visible', vis, 'Units', 'pixels', ...
    'Position', [40 60 1750 620], 'NumberTitle', 'off', 'Name', fig_title);
tl = tiledlayout(fig, 2, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(2, nc);
for k = 1:nc
    R = runs{k}; M = R.meas; t = R.t;
    is_ied = strcmp(R.network, 'run13');
    Rall = zeros(numel(R.groups), numel(t));
    Rall(R.E_idx, :) = R.rE; Rall(R.I_idx, :) = R.rI;
    if is_ied
        is_I = false(numel(R.groups), 1); is_I(R.I_idx) = true;
        [~, order] = sortrows([R.groups, is_I, (1:numel(R.groups))']);
    else
        order = [R.E_idx(:); R.I_idx(:)];
    end
    % 0.05-s bins keep the saved .fig small
    step = max(1, round(0.05 / (t(2) - t(1))));
    nb = floor(numel(t) / step);
    img = squeeze(mean(reshape(Rall(order, 1:nb*step), numel(order), step, nb), 2));
    tb = mean(reshape(t(1:nb*step), step, nb), 1);

    ax(1, k) = nexttile(tl, k);
    imagesc(ax(1, k), tb, 1:numel(order), single(img), [0 1]);
    colormap(ax(1, k), gray); set(ax(1, k), 'YDir', 'reverse'); hold(ax(1, k), 'on');
    if is_ied
        edges = cumsum(accumarray(R.groups, 1)) + 0.5;
        yline(ax(1, k), edges(1:end-1), '-', 'Color', [1 1 1], 'LineWidth', 0.6);
        what = sprintf('%d events', M.count);
    else
        yline(ax(1, k), numel(R.E_idx) + 0.5, '-', 'Color', cI, 'LineWidth', 1);
        what = sprintf('%d bursts', M.count);
    end
    title(ax(1, k), {titles{k}, sprintf('%s, r_E %.2f, r_I %.2f', what, M.mean_rE, M.mean_rI), ...
        M.regime}, 'FontSize', 8, 'FontWeight', 'normal');
    if k == 1
        if is_ied, ylabel(ax(1, k), 'neuron, by module'); else, ylabel(ax(1, k), 'neuron (E above, I below)'); end
    else
        set(ax(1, k), 'YTickLabel', []);
    end
    set(ax(1, k), 'XTickLabel', []);

    ax(2, k) = nexttile(tl, nc + k); hold(ax(2, k), 'on');
    if ~is_ied
        for b = 1:M.n_bursts
            patch(ax(2, k), M.bursts.t_on(b) + [0 0 0 0] + [0 1 1 0] * (M.bursts.t_off(b) - M.bursts.t_on(b)), ...
                [0 0 1 1], [0.88 0.88 0.88], 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end
    if is_ied
        % the mean E rate of each module, behind: an event is local
        for g = 1:max(R.groups)
            plot(ax(2, k), t, mean(R.rE(R.groups(R.E_idx) == g, :), 1), 'Color', [0.72 0.72 0.72], ...
                'LineWidth', 0.3, 'HandleVisibility', 'off');
        end
    end
    hI = plot(ax(2, k), t, mean(R.rI, 1), 'Color', cI, 'LineWidth', 0.5);
    hE = plot(ax(2, k), t, mean(R.rE, 1), 'Color', cE, 'LineWidth', 0.5);
    yline(ax(2, k), 0.2, ':k', 'HandleVisibility', 'off');
    if is_ied && ~isempty(M.event_times)
        plot(ax(2, k), M.event_times, 0.95 * ones(size(M.event_times)), 'v', 'MarkerSize', 3, ...
            'MarkerFaceColor', [0.5 0.2 0.6], 'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');
    end
    ylim(ax(2, k), [0 1]); xlabel(ax(2, k), 'time (s)');
    if k == 1
        ylabel(ax(2, k), 'mean rate');
        legend(ax(2, k), [hE hI], {'E', 'I'}, 'Location', 'east', 'Box', 'off', 'FontSize', 7);
    else
        set(ax(2, k), 'YTickLabel', []);
    end
end
linkaxes(ax(:), 'x'); xlim(ax(1, 1), [runs{1}.t(1) runs{1}.t(end)]);
linkaxes(ax(1, :), 'y'); linkaxes(ax(2, :), 'y');
set(ax, 'FontSize', 8, 'Box', 'off', 'TickDir', 'out');
title(tl, fig_title, 'FontSize', 10, 'Interpreter', 'none');
end
