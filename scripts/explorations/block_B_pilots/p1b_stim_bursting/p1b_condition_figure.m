function fig = p1b_condition_figure(cases, S, title_str, visible)
% P1B_CONDITION_FIGURE One condition under every stimulus case: columns are cases.
%
%   fig = P1B_CONDITION_FIGURE(cases, S, title_str, visible)
%
% cases  1 x nK cell of p1b_run_case outputs, one condition
% Rows, as in pilot 0b plus the input: rate of every neuron (E above the
% blue line); population mean E and I rate with bursts shaded grey and
% sustained above-threshold stretches shaded orange, the burst threshold
% dotted; mean adaptation term over E; mean depression factor over the
% depressed routes; the input (mean over the stimulated neurons, and over
% all neurons). The stimulus onset is a red dashed line in every panel.
% Every row shares its y axis.

nK = numel(cases);
vis = 'off'; if visible, vis = 'on'; end
fig = figure('Color', 'w', 'Visible', vis, 'Units', 'pixels', ...
    'Position', [40 40 1500 1000], 'NumberTitle', 'off', 'Name', title_str);
tl = tiledlayout(fig, 5, nK, 'TileSpacing', 'compact', 'Padding', 'compact');
cE = [0.85 0.33 0.10]; cI = [0.00 0.45 0.74]; cOn = [0.80 0.10 0.10];
ax = gobjects(5, nK);
nE = size(cases{1}.r_E, 1);
for k = 1:nK
    C = cases{k}; B = C.bursts; M = C.measures;
    is_burst = B.duration <= S.burst.sustained_max;
    % 1 rates
    ax(1,k) = nexttile(tl, k);
    imagesc(ax(1,k), C.t, 1:(nE + size(C.r_I, 1)), double([C.r_E; C.r_I]), [0 1]);
    colormap(ax(1,k), gray); set(ax(1,k), 'YDir', 'reverse'); hold(ax(1,k), 'on');
    yline(ax(1,k), nE + 0.5, '-', 'Color', cI, 'LineWidth', 1);
    xline(ax(1,k), S.t_stim, '--', 'Color', cOn, 'LineWidth', 1.2);
    hold(ax(1,k), 'off');
    title(ax(1,k), sprintf('%s\nbursts/s %.2f -> %.2f, mean r_E %.3f -> %.3f', ...
        C.stim_label, M.before.burst_rate, M.during.burst_rate, ...
        M.before.mean_r_E, M.during.mean_r_E), 'Interpreter', 'tex', 'FontSize', 9);
    if k == 1, ylabel(ax(1,k), sprintf('neuron (E 1-%d, I below)', nE)); end
    % 2 mean rates
    ax(2,k) = nexttile(tl, nK + k); hold(ax(2,k), 'on');
    for b = 1:B.n
        col = [0.85 0.85 0.85];
        if ~is_burst(b), col = [1.00 0.88 0.72]; end
        patch(ax(2,k), [B.t_on(b) B.t_off(b) B.t_off(b) B.t_on(b)], [0 0 1 1], ...
            col, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    hI = plot(ax(2,k), C.t, C.mean_r_I, 'Color', cI, 'DisplayName', 'I');
    hE = plot(ax(2,k), C.t, C.mean_r_E, 'Color', cE, 'DisplayName', 'E');
    yline(ax(2,k), S.burst.thresh, ':k', 'HandleVisibility', 'off');
    xline(ax(2,k), S.t_stim, '--', 'Color', cOn, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    hold(ax(2,k), 'off');
    if k == 1
        ylabel(ax(2,k), 'mean rate');
        legend(ax(2,k), [hE hI], 'Location', 'northwest', 'Box', 'off');
    end
    % 3 adaptation
    ax(3,k) = nexttile(tl, 2*nK + k);
    plot(ax(3,k), C.t, C.mean_adaptation, 'Color', cE);
    xline(ax(3,k), S.t_stim, '--', 'Color', cOn, 'LineWidth', 1.2);
    if k == 1, ylabel(ax(3,k), 'mean (c/K)\Sigma_k a_k, E'); end
    % 4 depression
    ax(4,k) = nexttile(tl, 3*nK + k);
    plot(ax(4,k), C.t, C.mean_depression, 'Color', [0.3 0.3 0.3]);
    xline(ax(4,k), S.t_stim, '--', 'Color', cOn, 'LineWidth', 1.2);
    if k == 1
        ylabel(ax(4,k), sprintf('mean STD factor\n(%s)', strjoin(C.depression_routes, ', ')));
    end
    % 5 input
    ax(5,k) = nexttile(tl, 4*nK + k); hold(ax(5,k), 'on');
    hT = plot(ax(5,k), C.t, C.u_target, 'Color', [0.2 0.2 0.2], ...
        'DisplayName', 'mean over stimulated');
    hA = plot(ax(5,k), C.t, C.u_mean, '--', 'Color', [0.55 0.35 0.75], ...
        'DisplayName', 'mean over all 50');
    xline(ax(5,k), S.t_stim, '--', 'Color', cOn, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    hold(ax(5,k), 'off');
    if k == 1
        ylabel(ax(5,k), 'input u');
        legend(ax(5,k), [hT hA], 'Location', 'northwest', 'Box', 'off');
    end
    xlabel(ax(5,k), 'time (s)');
end
linkaxes(ax(:), 'x'); xlim(ax(1,1), [0 S.T_show]);
for r = 1:5, linkaxes(ax(r,:), 'y'); end
ylim(ax(2,1), [0 1]);
lo = min(cellfun(@(C) min(C.mean_depression), cases));
ylim(ax(4,1), [lo * 0.95, 1]);
hi = max(cellfun(@(C) max(C.mean_adaptation), cases));
ylim(ax(3,1), [0, hi * 1.05 + eps]);
ulo = min(cellfun(@(C) min([C.u_target C.u_mean]), cases));
uhi = max(cellfun(@(C) max([C.u_target C.u_mean]), cases));
ylim(ax(5,1), [min(0, ulo) - 0.02, uhi + 0.05]);
for a = ax(:)', set(a, 'Box', 'off', 'TickDir', 'out'); end
title(tl, title_str, 'Interpreter', 'none', 'FontSize', 10);
end
