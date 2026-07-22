%% Plot the Stability_Manuscript example memory-capacity figure.
% Loads mc_example_data.mat (written by compute_memory_capacity_example.m) and
% renders a single manuscript figure -- no simulation is re-run, so the look can
% be iterated quickly. Run compute_memory_capacity_example.m first.
%
% Layout (styled to match fig_memory_capacity/Fig_memory_capacity.m):
%   (a) Cumulative Memory Capacity vs delay (0-10 s), all 4 conditions.
%   (b) Per-delay R^2 vs delay (0-10 s), all 4 conditione t a few delays -- target
%   u(t-d) vs the trained readout -- each panel titled with the delay in seconds
%   and its R^2. Built from the saved mc_results (predictions/t_pred/R2_d).

clear; clc; close all;

% Assumes setup_paths.m has already been run (src/ + scripts/ on the MATLAB path).
% this_dir locates mc_example_data.mat and is where the final figure is saved.
this_dir = fileparts(mfilename('fullpath'));
% scripts/presentations/Stability_Manuscript/fig_memory_capacity_example -> root is 4 up
project_root = fileparts(fileparts(fileparts(fileparts(this_dir))));
out_dir = this_dir;   % write the final figure next to this script

%% Load precomputed results
data_file = fullfile(this_dir, 'mc_example_data.mat');
assert(isfile(data_file), ...
    'Missing %s -- run compute_memory_capacity_example.m first.', data_file);
S = load(data_file);
results         = S.results;
R2              = S.R2;
delay_s         = S.delay_s;
condition_names = S.condition_names;
n_cond          = numel(condition_names);

%% Style + palette (matched to Fig_memory_capacity.m / plot_memory_capacity_combined.m)
set(0,'DefaultAxesFontSize',14);    % tick numbers + x/y labels
set(0,'DefaultAxesLineWidth',1.0);  % axis lines + tick marks
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultLegendInterpreter','none');

% Okabe-Ito colorblind-safe palette (same as the combined MC figure):
%   Baseline black, SFA orange, STD sky blue, SFA+STD reddish purple.
colors = [0.000 0.000 0.000;
          0.902 0.624 0.000;
          0.337 0.706 0.914;
          0.800 0.475 0.655];
if size(colors,1) < n_cond
    colors = lines(n_cond);
end

% Short condition labels for the legend (match the combined figure).
if n_cond == 4
    display_names = {'Baseline', 'SFA', 'STD', 'SFA+STD'};
else
    display_names = condition_names;
end

xmax_s       = 7.5;             % (a)/(b) delay-axis limit (s)
recon_cond   = 4;               % reconstruct the SFA+STD condition
recon_delays = [1, 5, 10, 15];  % hold-delay indices to show (labeled in seconds): 0.3, 1.5, 3.0, 4.5 s at T_hold=0.3
recon_ylim   = [-0.6, 0.6];     % shared y-limits across all reconstruction panels
line_w       = 2;               % (a)/(b) curve width (matches the reference figure)

%% ==================== Combined figure ====================
% 24x2 tiled grid (fine rows so the gap below (a)/(b) is tunable): (a)/(b) span
% the top 8 rows; row 9 is left empty as a spacer; the reconstruction block
% (nested layout) fills rows 10-24. Explicit tile indices keep the spans clear.
fig = figure('Color', 'w', 'Position', [100 237 603 652]);   % compact size (hand-tuned)
tl = tiledlayout(24, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
% Leave headroom above the tiles so the (a)/(b) panel letters can sit clearly
% above the axes (and their y-labels) without clipping at the figure top. The
% reduced width leaves room at the right for the (c) amplitude scale-bar label.
tl.OuterPosition = [0 0 0.90 0.95];

% (a) Cumulative MC vs delay
ax_a = nexttile(1, [8 1]); hold on; grid off; box off;
for i = 1:n_cond
    plot(delay_s, cumsum(R2{i}), '-', 'Color', colors(i, :), 'LineWidth', line_w);
end
xlim([0 xmax_s]); ylim([0 11]);
xlabel('Delay (s)', 'FontSize', 15.4); ylabel({'Cumulative', 'Memory Capacity'}, 'FontSize', 15.4);
set(ax_a, 'XTick', [0 2.5 5], 'YTick', [0 5 10]);   % match plot_memory_capacity_combined.m
hold off;

% (b) Per-delay R^2 vs delay -- carries the one legend for the figure
ax_b = nexttile(2, [8 1]); hold on; grid off; box off;
for i = 1:n_cond
    plot(delay_s, R2{i}, '-', 'Color', colors(i, :), 'LineWidth', line_w);
end
xlim([0 xmax_s]); ylim([0 1.005]);
xlabel('Delay (s)', 'FontSize', 15.4); ylabel('$R^2$', 'Interpreter', 'latex', 'FontSize', 15.4);
set(ax_b, 'XTick', [0 2.5 5], 'YTick', [0 0.5 1]);   % match plot_memory_capacity_combined.m
lgd = legend(display_names, 'Location', 'northeast', 'Box', 'off');
% Nudge the legend further up and to the right (into the corner headroom).
lgd.Units = 'normalized';
lp = lgd.Position;
lgd.Position = [lp(1) + 0.02, lp(2) + 0.03, lp(3), lp(4)];
hold off;

% Reconstruction rows (SFA+STD): target u(t-d) vs trained readout, delay in s.
% The x-axis (time) is hidden on every row -- a 15 s scale bar on the bottom row
% conveys the timescale instead.
mcr        = results{recon_cond};
nR         = numel(recon_delays);
scalebar_s  = 15;             % length of the horizontal time scale bar (s)
scalebar_au = 1;              % length of the vertical amplitude scale bar (A.U.)
recon_win   = 60;             % time span shown in each reconstruction panel (s)
recon_ax   = gobjects(1, nR);

% Nested tiledlayout for the (c) rows so their vertical spacing can be tightened
% independently of (a)/(b). It occupies grid rows 10-24 of the outer 24x2 layout
% (row 9 above it is the empty spacer that opens the gap below (a)/(b)); 'tight'
% spacing packs the four panels closer than the outer 'compact' while still
% leaving room for each panel's title.
tl_c = tiledlayout(tl, nR, 1, 'TileSpacing', 'tight', 'Padding', 'none');
tl_c.Layout.Tile     = 19;       % top-left tile of the (c) block (row 10, col 1)
tl_c.Layout.TileSpan = [15 2];   % span rows 10-24 x 2 cols
for k = 1:nR
    d = recon_delays(k);
    p = mcr.predictions(d);
    t_delay = mcr.t_pred(p.t_indices);

    recon_ax(k) = nexttile(tl_c); hold on; grid off; box off;
    plot(t_delay, p.y_true, 'k-', 'LineWidth', 0.9);
    plot(t_delay, p.y_pred, '-', 'Color', colors(recon_cond, :), 'LineWidth', 1.4);
    ylim(recon_ylim);
    % Show only the first recon_win seconds of the test window.
    xlim([t_delay(1), t_delay(1) + recon_win]);
    title(sprintf('delay = %.1f seconds,  R^2 = %.2f', delay_s(d), mcr.R2_d(d)), ...
        'FontWeight', 'normal', 'Interpreter', 'tex');
    % Hide both axes entirely (no ticks, no black axis lines) -- the scale bars
    % convey time (x) and amplitude (y) instead.
    set(gca, 'XColor', 'none', 'YColor', 'none');
    if k == nR
        % L-shaped scale bar in the lower-right corner: a horizontal 15 s time bar
        % and a vertical 1 A.U. amplitude bar meeting at a right angle. Clipping is
        % off so the bars/labels can sit beneath and right of the axis limits.
        xl = get(gca, 'XLim');
        xc = xl(2) + 0.02 * diff(xl);                    % corner x (just right of the data)
        x1 = xc - scalebar_s;                            % left end of the time bar
        yb = recon_ylim(1) - 0.18 * diff(recon_ylim);    % corner y (below the trace)
        ya = yb + scalebar_au;                           % top of the amplitude bar

        % Horizontal time bar + vertical amplitude bar, sharing the corner (xc, yb).
        plot([x1 xc], [yb yb], 'k-', 'LineWidth', 2.5, 'Clipping', 'off');
        plot([xc xc], [yb ya], 'k-', 'LineWidth', 2.5, 'Clipping', 'off');

        % Time label centered below the horizontal bar.
        text((x1 + xc) / 2, yb - 0.06 * diff(recon_ylim), sprintf('%d s', scalebar_s), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontSize', 12, 'Clipping', 'off');
        % Amplitude label just right of the vertical bar, rotated -90 deg.
        text(xc + 0.03 * diff(xl), (yb + ya) / 2, sprintf('%g AU', scalebar_au), ...
            'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 12, 'Clipping', 'off');
    end
    hold off;
end

% Panel letters: (a)/(b) on the analysis panels via AddLetters2Plots. The
% remaining reconstruction rows are self-labeled by their delay/R^2 titles.
% Per-axis VShift: a modest upward nudge for (a)/(b) -- the two-line y-label on
% (a) frees up space, so the letters can sit lower than before.
AddLetters2Plots({ax_a, ax_b}, {'(a)', '(b)'}, ...
    'FontSize', 18, 'FontWeight', 'normal', 'HShift', -0.06, 'VShift', [-0.02 -0.02]);

% (c) is lettered manually: recon_ax(1) lives in the nested tiledlayout, whose
% reported position confuses AddLetters2Plots (it lands up by (a)). Place it in
% axis-normalized coordinates at the panel's top-left corner instead.
text(recon_ax(1), -0.025, 1.16, '(c)', 'Units', 'normalized', ...
    'FontSize', 18, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'Clipping', 'off');

%% ==================== Save the figure ====================
% Stable names (png/svg/fig), matching the other Stability_Manuscript scripts.
save_fig_stable(fig, out_dir, 'Fig_MC_Example');

% Log the git state alongside the figure.
capture_git_provenance(out_dir, project_root);

%% -------------------- Human-readable description --------------------
desc_path = fullfile(out_dir, 'README_fig_memory_capacity_example.txt');
fid = fopen(desc_path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Stability_Manuscript figure: Example memory capacity\n');
fprintf(fid, '====================================================\n\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now')));
fprintf(fid, 'By script: %s.m\n\n', mfilename);

fprintf(fid, 'HOW IT WAS MADE\n');
fprintf(fid, '  Two-step, no re-simulation at plot time: compute_memory_capacity_example.m\n');
fprintf(fid, '  runs the memory-capacity protocol for the 4 adaptation conditions and saves\n');
fprintf(fid, '  the per-condition mc_results structs to mc_example_data.mat (gitignored).\n');
fprintf(fid, '  This script loads that file and renders the figure, so the look can be\n');
fprintf(fid, '  iterated without re-running the sim. See git_provenance.txt for the commit.\n\n');

fprintf(fid, 'MODEL SETTINGS\n');
fprintf(fid, '  Match looped_memory_capacity.m (c_E = 0.5/3, sample_hold input, n = 300,\n');
fprintf(fid, '  f = 0.6, level_of_chaos = 2.5). See compute_memory_capacity_example.m.\n\n');

fprintf(fid, 'FIGURE PRODUCED (in this folder)\n');
fprintf(fid, '  Fig_MC_Example.png / .svg / .fig\n');
fprintf(fid, '    (a) Cumulative Memory Capacity vs delay (0-%d s), all 4 conditions.\n', xmax_s);
fprintf(fid, '    (b) Per-delay R^2 vs delay (0-%d s), all 4 conditions (legend).\n', xmax_s);
fprintf(fid, '    Below: SFA+STD input reconstruction (target vs readout) for hold-delays\n');
fprintf(fid, '    %s (delay indices), each titled with the delay in seconds and R^2;\n', mat2str(recon_delays));
fprintf(fid, '    all reconstruction panels share y-limits [%g, %g].\n', recon_ylim(1), recon_ylim(2));

clear cleanup;  % flush + close
fprintf('Description written: %s\n', desc_path);

fprintf('Rendered + saved the example memory-capacity figure.\n');

%% ==================== Local helpers ====================
function save_fig_stable(fig, out_dir, base)
% Save fig as <base>.{png,svg,fig} with a stable name: clear any prior <base>*
% outputs, save via save_some_figs_to_folder_2 (which suffixes the figure
% number), then rename to the fixed names.
    old = dir(fullfile(out_dir, [base '*']));
    for a = 1:numel(old)
        fp = fullfile(old(a).folder, old(a).name);
        if ~old(a).isdir && (endsWith(fp, '.png') || endsWith(fp, '.svg') || endsWith(fp, '.fig'))
            delete(fp);
        end
    end
    save_some_figs_to_folder_2(out_dir, base, fig.Number, {'png', 'svg', 'fig'});
    num = num2str(fig.Number);
    movefile(fullfile(out_dir, [base '_figure_' num '.png']), fullfile(out_dir, [base '.png']));
    movefile(fullfile(out_dir, [base '_figure_' num '.svg']), fullfile(out_dir, [base '.svg']));
    movefile(fullfile(out_dir, [base '_f_' num '.fig']),      fullfile(out_dir, [base '.fig']));
end
