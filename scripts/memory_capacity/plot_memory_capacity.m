function [fig1, fig2] = plot_memory_capacity(results_all, out_dir)
% PLOT_MEMORY_CAPACITY Regenerate the looped-MC figures from a results struct.
%
%   [fig1, fig2] = PLOT_MEMORY_CAPACITY(results_all)
%   [fig1, fig2] = PLOT_MEMORY_CAPACITY(results_all, out_dir)
%
% Rebuilds Fig1 (paired total-MC + horizon distributions) and Fig2 (per-delay
% R^2 and cumulative MC, mean +/- 95% CI) from the `results_all` struct saved by
% run_memory_capacity.m -- no simulation is re-run. This is the plotting half
% of the compute/plot split, so figure styling can be iterated freely:
%   S = load('..._results.mat'); plot_memory_capacity(S.results_all);
%
% Inputs:
%   results_all : struct saved by run_memory_capacity.m (settings, summary,
%                 conditions, MC_trials, H_trials, ...).
%   out_dir     : optional. If given (non-empty), save <run_tag>_Fig1/Fig2.{png,pdf}
%                 there. If omitted, the figures are only displayed.

    if nargin < 2; out_dir = ''; end

    % --- Unpack (no recompute; everything is already in results_all) ---
    condition_names = results_all.conditions;
    n_cond   = numel(condition_names);
    MC_trials = results_all.MC_trials;
    H_trials  = results_all.H_trials;
    n_trials  = size(MC_trials, 1);

    s = results_all.summary;
    MC_mean = s.MC_mean;  MC_ci = s.MC_ci95;
    H_mean  = s.H_mean;   H_ci  = s.H_ci95;
    R2_mean = s.R2_mean;  R2_ci = s.R2_ci95;

    % Cumulative curves (cheap to recompute from the per-delay means/CIs)
    R2_cum_mean  = cumsum(R2_mean, 2);
    R2_cum_ci_lo = cumsum(R2_ci.lo, 2);
    R2_cum_ci_hi = cumsum(R2_ci.hi, 2);

    cfg = results_all.settings;
    delay_s    = cfg.delay_s;
    input_type = cfg.input_type;
    n          = cfg.n;
    R2_thr     = cfg.R2_threshold_for_horizon;
    run_tag    = results_all.run_tag;

    % --- Style + palette (edit here to restyle) ---
    % Root defaults are set through with_graphics_defaults so they are RESTORED
    % when this function returns. Previously these were bare set(0,...) calls
    % that leaked into the rest of the session: after any memory-capacity plot,
    % DefaultTextInterpreter stayed 'none', which renders the sensitivity and
    % param-space sheets' '\lambda_1' / '\mu_{EE}' labels as literal backslash
    % text. Verified live before the fix. `style_cleanup` must stay in scope.
    style_cleanup = with_graphics_defaults( ...
        'DefaultAxesFontSize',      12, ...
        'DefaultTextInterpreter',   'none', ...
        'DefaultLegendInterpreter', 'none'); %#ok<NASGU>

    % Colors (edit here to restyle). Black / blue / green / red; the light CI
    % fill alpha (below) keeps overlapping bands from muddying.
    colors = [0.00 0.00 0.00;   % Baseline: black
              0.00 0.45 0.74;   % SFA:      blue
              0.20 0.60 0.20;   % STD:      green
              0.84 0.15 0.16];  % SFA+STD:  red
    if size(colors,1) < n_cond
        colors = lines(n_cond);   % fallback if more conditions than palette rows
    end

    xpos = 1:n_cond;

    % ================= Figure 1: total MC + horizon (paired) =================
    fig1 = figure('Color','w','Position',[100 100 1050 420]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    % 1A: Paired scatter with lines
    nexttile; hold on; grid off; box on;
    for k = 1:n_trials
        plot(xpos, MC_trials(k,:), '-', 'Color', [0 0 0 0.15], 'LineWidth', 1);
    end
    for i = 1:n_cond
        scatter(i*ones(n_trials,1), MC_trials(:,i), 20, ...
            'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor','none', ...
            'MarkerFaceAlpha',0.65);
    end
    for i = 1:n_cond
        plot(i, MC_mean(i), 'o', 'MarkerSize', 8, ...
            'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
        plot([i i], [MC_ci.lo(i) MC_ci.hi(i)], '-', 'Color', 'k', 'LineWidth', 2);
    end
    xlim([0.5 n_cond+0.5]);
    set(gca,'XTick',xpos,'XTickLabel',condition_names);
    ylabel('Total Memory Capacity (sum R^2)');
    title('Total MC (paired trials)');

    % 1B: Horizon distribution
    nexttile; hold on; grid off; box on;
    for k = 1:n_trials
        plot(xpos, H_trials(k,:), '-', 'Color', [0 0 0 0.15], 'LineWidth', 1);
    end
    for i = 1:n_cond
        scatter(i*ones(n_trials,1), H_trials(:,i), 20, ...
            'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor','none', ...
            'MarkerFaceAlpha',0.65);
    end
    for i = 1:n_cond
        plot(i, H_mean(i), 'o', 'MarkerSize', 8, ...
            'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
        plot([i i], [H_ci.lo(i) H_ci.hi(i)], '-', 'Color', 'k', 'LineWidth', 2);
    end
    xlim([0.5 n_cond+0.5]);
    set(gca,'XTick',xpos,'XTickLabel',condition_names);
    ylabel(sprintf('Memory horizon (s), R^2 > %.2f', R2_thr));
    title('Horizon (paired trials)');

    sgtitle(sprintf('Memory Capacity Comparison (%s input, n=%d, trials=%d)', ...
        input_type, n, n_trials));

    % ================= Figure 2: R^2(d) + cumulative (mean +/- CI) ===========
    fig2 = figure('Color','w','Position',[100 100 1050 420]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    % 2A: R^2(d)
    nexttile; hold on; grid off; box on;
    for i = 1:n_cond
        shaded_ci(delay_s, R2_ci.lo(i,:), R2_ci.hi(i,:), colors(i,:), 0.12);
        plot(delay_s, R2_mean(i,:), '-', 'Color', colors(i,:), 'LineWidth', 2);
    end
    xlabel('Delay (s)');
    ylabel('R^2(d)');
    title('Per-delay memory (mean ± 95% CI)');
    legend(condition_names,'Location','northeast');

    % 2B: Cumulative MC vs delay
    nexttile; hold on; grid off; box on;
    for i = 1:n_cond
        shaded_ci(delay_s, R2_cum_ci_lo(i,:), R2_cum_ci_hi(i,:), colors(i,:), 0.12);
        plot(delay_s, R2_cum_mean(i,:), '-', 'Color', colors(i,:), 'LineWidth', 2);
    end
    xlabel('Delay (s)');
    ylabel('Cumulative MC (sum_{j<=d} R^2(j))');
    title('Cumulative memory (mean ± 95% CI)');
    legend(condition_names,'Location','southeast');

    sgtitle(sprintf('Delay Profile of Memory (%s input, n=%d, trials=%d)', ...
        input_type, n, n_trials));

    % --- Optional save ---
    if ~isempty(out_dir)
        if ~isfolder(out_dir); mkdir(out_dir); end
        saveas(fig1, fullfile(out_dir, [run_tag '_Fig1_MC_Distributions.png']));
        print(fig1, fullfile(out_dir, [run_tag '_Fig1_MC_Distributions.pdf']), '-dpdf', '-painters');
        saveas(fig2, fullfile(out_dir, [run_tag '_Fig2_R2_Curves.png']));
        print(fig2, fullfile(out_dir, [run_tag '_Fig2_R2_Curves.pdf']), '-dpdf', '-painters');
        fprintf('Figures saved to:\n  %s\n', out_dir);
    end
end

%% ==================== Local helper ====================
function shaded_ci(x, lo, hi, rgb, alpha_fill)
% Draw shaded confidence interval [lo, hi] around x. HandleVisibility off so the
% CI patches don't appear in legends (only the mean lines should be labelled).
    x = x(:)'; lo = lo(:)'; hi = hi(:)';
    X = [x, fliplr(x)];
    Y = [hi, fliplr(lo)];
    fill(X, Y, rgb, 'FaceAlpha', alpha_fill, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end
