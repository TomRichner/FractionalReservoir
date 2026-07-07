function fig3 = plot_memory_capacity_combined(results_all, out_dir)
% PLOT_MEMORY_CAPACITY_COMBINED  Paper-ready 1x3 memory-capacity figure.
%
%   fig3 = PLOT_MEMORY_CAPACITY_COMBINED(results_all)
%   fig3 = PLOT_MEMORY_CAPACITY_COMBINED(results_all, out_dir)
%
% Assembles a single manuscript figure from the pieces of the two figures made
% by scripts/memory_capacity/plot_memory_capacity.m, laid out as a 1x3 strip:
%   (a) cumulative memory   (from Fig2 panel 2B)
%   (b) per-delay memory    (from Fig2 panel 2A)
%   (c) horizon paired trials (from Fig1 panel 1B)
%
% Everything is read from the saved `results_all` struct (see
% looped_memory_capacity.m) -- no simulation is re-run. Styling is kept local to
% this presentation copy so it can be tuned without touching the shared plotter.
%
% Inputs:
%   results_all : struct saved by looped_memory_capacity.m.
%   out_dir     : optional. If non-empty, save Fig_Memory_Capacity.{png,svg,fig}
%                 there (via save_some_figs_to_folder_2, which appends the figure
%                 number). If omitted, the figure is only displayed.

    if nargin < 2; out_dir = ''; end

    % Ensure AddLetters2Plots (src/plotting) is resolvable even if this function
    % is called standalone. Project root is 4 folders up from this file.
    if exist('AddLetters2Plots', 'file') ~= 2
        project_root = fileparts(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
        addpath(genpath(fullfile(project_root, 'src')));
    end

    % --- Unpack (no recompute; everything is already in results_all) ---
    condition_names = results_all.conditions;
    n_cond   = numel(condition_names);
    H_trials = results_all.H_trials;
    n_trials = size(H_trials, 1);

    s = results_all.summary;
    R2_mean = s.R2_mean;  R2_ci = s.R2_ci95;

    % Cumulative curves (recomputed from the per-delay means/CIs, exactly as in
    % plot_memory_capacity.m).
    R2_cum_mean  = cumsum(R2_mean, 2);
    R2_cum_ci_lo = cumsum(R2_ci.lo, 2);
    R2_cum_ci_hi = cumsum(R2_ci.hi, 2);

    cfg = results_all.settings;
    delay_s = cfg.delay_s;

    % --- Style + palette (edit here to restyle) ---
    set(0,'DefaultAxesFontSize',14);   % drives tick numbers + x/y labels
    set(0,'DefaultAxesLineWidth',1.0); % axis lines + tick marks
    set(0,'DefaultTextInterpreter','none');
    set(0,'DefaultLegendInterpreter','none');

    % Okabe-Ito colorblind-safe qualitative palette. Reddish-purple (instead of
    % bluish-green) keeps all four hues well separated -- sky-blue and bluish-
    % green were too similar -- and avoids green entirely.
    colors = [0.000 0.000 0.000;   % Baseline: black           #000000
              0.902 0.624 0.000;   % SFA:      orange          #E69F00
              0.337 0.706 0.914;   % STD:      sky blue        #56B4E9
              0.800 0.475 0.655];  % SFA+STD:  reddish purple  #CC79A7
    if size(colors,1) < n_cond
        colors = lines(n_cond);   % fallback if more conditions than palette rows
    end

    xpos = 1:n_cond;

    % ================= Figure 3: combined 1x3 strip =================
    fig3 = figure('Color','w','Position',[100 100 1350 380]);
    tl = tiledlayout(1,3,'Padding','compact','TileSpacing','loose');
    % Leave headroom above the tiles (so the panel letters sit above the axes
    % without clipping) and a margin below (so the x tick labels aren't cut off).
    tl.OuterPosition = [0 0.06 1 0.86];

    % (a) Cumulative MC vs delay
    nexttile; hold on; grid off; box off;
    for i = 1:n_cond
        shaded_ci(delay_s, R2_cum_ci_lo(i,:), R2_cum_ci_hi(i,:), colors(i,:), 0.12);
        plot(delay_s, R2_cum_mean(i,:), '-', 'Color', colors(i,:), 'LineWidth', 2);
    end
    xlim([0 10]);
    xlabel('Delay (s)');
    ylabel('Cumulative Memory Capacity');

    % (b) Per-delay R^2(d) -- carries the one legend for the figure
    nexttile; hold on; grid off; box off;
    for i = 1:n_cond
        shaded_ci(delay_s, R2_ci.lo(i,:), R2_ci.hi(i,:), colors(i,:), 0.12);
        plot(delay_s, R2_mean(i,:), '-', 'Color', colors(i,:), 'LineWidth', 2);
    end
    xlim([0 10]);
    xlabel('Delay (s)');
    ylabel('$R^2$', 'Interpreter', 'latex');
    legend(condition_names, 'Location', 'northeast');

    % (c) Horizon distribution (paired trials)
    nexttile; hold on; grid off; box off;
    % Faint paired-trial lines. Use a solid light-gray RGB rather than an
    % alpha'd color ([0 0 0 0.15]): the 4th (alpha) element is an undocumented
    % runtime-only transparency that is NOT saved in .fig files, so a saved/
    % reopened figure would otherwise show these as opaque black.
    for k = 1:n_trials
        plot(xpos, H_trials(k,:), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.75);
    end
    for i = 1:n_cond
        scatter(i*ones(n_trials,1), H_trials(:,i), 20, ...
            'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor','none', ...
            'MarkerFaceAlpha',0.65);
    end
    xlim([0.5 n_cond+0.5]);
    set(gca,'XTick',xpos,'XTickLabel',condition_names);
    ylabel('Memory Horizon (s)');

    % Panel letters (a)/(b)/(c) just above-left of each tile, outside the axes
    % (negative shifts push the label outside the plotting area). Uniform shift
    % for all three panels.
    AddLetters2Plots(fig3, {'(a)', '(b)', '(c)'}, ...
        'FontSize', 18, 'FontWeight', 'normal', ...
        'HShift', -0.04, 'VShift', -0.09);

    % --- Optional save ---
    % Save via save_some_figs_to_folder_2 (png + svg + fig). Pass this figure's
    % number explicitly so the open Fig1/Fig2 (from replot_memory_capacity) are
    % not also re-saved here.
    if ~isempty(out_dir)
        save_some_figs_to_folder_2(out_dir, 'Fig_Memory_Capacity', ...
            fig3.Number, {'png', 'svg', 'fig'});
        fprintf('Combined figure saved to:\n  %s\n', out_dir);
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
