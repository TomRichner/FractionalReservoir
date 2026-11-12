function E = p6_examples_figure()
% P6_EXAMPLES_FIGURE Block B pilot 6, question 3: what the classes look like.
%
% Reads the traces written by p6_rebuild('show', ...) and p6_rebuild('verify',
% ...). One row per rebuilt network; for each of no adaptation and chained
% MTS, the population mean E and I rate and a rate image (60 E neurons above
% 60 I neurons). Also computes simple EVENT measures from the population
% rate after t = 5 s, which the sweep cannot (it stores no time series), and
% writes them, with the stored-against-rebuilt check, to p6_examples.csv.

P = p6_config();
tr_dir = fullfile(P.out_dir, 'traces');
nNet = numel(P.rebuild_idx);
nShow = 60;

f = figure('Visible', 'off', 'Color', 'w', 'Position', [30 30 1600 1250]);
tl = tiledlayout(f, nNet, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
rows = cell(0, 1);
for k = 1:nNet
    for c = 1:2
        cond = P.rebuild_conds{c};
        S = load(fullfile(tr_dir, sprintf('show_%d_%s.mat', P.rebuild_idx(k), cond)));
        V = load(fullfile(tr_dir, sprintf('verify_%d_%s.mat', P.rebuild_idx(k), cond)));
        t = S.t; late = t >= 5;
        mE = double(S.mean_rE); mI = double(S.mean_rI);
        nE = size(S.rE, 1); nI = size(S.rI, 1);
        pop = (nE * mE + nI * mI) / (nE + nI);

        % Event measures on the population rate, t >= 5 s
        ev = struct();
        ev.mean_E = mean(mE(late)); ev.mean_I = mean(mI(late));
        ev.cv_E = std(mE(late)) / mean(mE(late));
        ev.cv_I = std(mI(late)) / mean(mI(late));
        ev.frac_E_sat = mean(mean(S.rE(:, late) > 0.9, 1));      % share of E neuron-time saturated
        thr = median(pop(late)) + 0.1;                           % a population event: +0.1 above the median
        above = pop(late) > thr;
        ev.n_events = sum(diff([0 above]) == 1);
        ev.frac_time_event = mean(above);
        ev.peak_over_median = max(pop(late)) / median(pop(late));

        rows{end + 1, 1} = table(P.rebuild_idx(k), string(P.rebuild_class{k}), string(cond), ...
            V.LLE_stored, V.LLE_rebuilt, V.rate_stored, V.rate_rebuilt, ...
            ev.mean_E, ev.mean_I, ev.cv_E, ev.cv_I, ev.frac_E_sat, ev.n_events, ...
            ev.frac_time_event, ev.peak_over_median, 'VariableNames', ...
            {'config_idx', 'class_under_chained_MTS', 'condition', 'lambda1_stored', 'lambda1_rebuilt', ...
            'rate_stored', 'rate_rebuilt', 'mean_E_30s', 'mean_I_30s', 'cv_popE', 'cv_popI', ...
            'frac_E_saturated', 'n_pop_events', 'frac_time_in_event', 'peak_over_median'}); %#ok<AGROW>

        % traces
        ax = nexttile(tl, (k - 1) * 4 + (c - 1) * 2 + 1); hold(ax, 'on');
        plot(ax, t, mI, '-', 'Color', [0.1 0.4 0.8], 'LineWidth', 0.9);
        plot(ax, t, mE, '-', 'Color', [0.85 0.2 0.1], 'LineWidth', 0.9);
        set(ax, 'XLim', [0 30], 'YLim', [0 1.02], 'FontSize', 9, 'Box', 'off');
        title(ax, sprintf('net %d, %s: \\lambda_1 = %.2f /s, rate %.2f', P.rebuild_idx(k), ...
            strrep(P.cond_lbl{strcmp(P.conds, cond)}, '_', '\_'), V.LLE_stored, V.rate_stored), ...
            'FontSize', 9, 'FontWeight', 'normal');
        if c == 1, ylabel(ax, sprintf('%s\nmean rate', P.rebuild_class{k}), 'FontSize', 8); end
        if k == 1 && c == 1, legend(ax, {'I', 'E'}, 'Box', 'off', 'FontSize', 8, 'Location', 'east'); end
        if k == nNet, xlabel(ax, 'time (s)'); end

        % image: nShow E above nShow I
        iE = round(linspace(1, nE, min(nShow, nE))); iI = round(linspace(1, nI, min(nShow, nI)));
        ax = nexttile(tl, (k - 1) * 4 + (c - 1) * 2 + 2);
        imagesc(ax, t, 1:(numel(iE) + numel(iI)), [S.rE(iE, :); S.rI(iI, :)], [0 1]);
        colormap(ax, flipud(gray(256)));
        hold(ax, 'on'); plot(ax, [0 30], numel(iE) + [0.5 0.5], '-', 'Color', [0.85 0.2 0.1], 'LineWidth', 0.8);
        set(ax, 'FontSize', 9, 'YTick', [numel(iE) / 2, numel(iE) + numel(iI) / 2], 'YTickLabel', {'E', 'I'}, ...
            'TickLength', [0 0]);
        title(ax, sprintf('rates (white 0, black 1); events %d, CV_E %.2f, CV_I %.2f', ...
            ev.n_events, ev.cv_E, ev.cv_I), 'FontSize', 9, 'FontWeight', 'normal');
        if k == nNet, xlabel(ax, 'time (s)'); end
    end
end
title(tl, ['Six sampled networks rebuilt from their saved parameters, 30 s each: no adaptation (left pair) ' ...
    'and chained MTS (right pair). \lambda_1 and rate in the titles are the sweep''s stored values.'], 'FontSize', 10);

E = vertcat(rows{:});
writetable(E, fullfile(P.out_dir, 'p6_examples.csv'));
disp(E(:, [1 3 4 5 6 7 8 9 10 11 12 13]));

set(f, 'Visible', 'on');
savefig(f, fullfile(P.out_dir, 'p6_timeseries_examples.fig'));
set(f, 'Visible', 'off');
exportgraphics(f, fullfile(P.out_dir, 'p6_timeseries_examples.png'), 'Resolution', 105);
close(f);
end
