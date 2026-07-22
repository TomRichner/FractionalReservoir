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

xmax_s       = 10;              % (a)/(b) delay-axis limit (s)
recon_cond   = 4;               % reconstruct the SFA+STD condition
recon_delays = [1, 5, 10, 15];  % hold-delay indices to show (labeled in seconds): 0.3, 1.5, 3.0, 4.5 s at T_hold=0.3
recon_ylim   = [-0.6, 0.6];     % shared y-limits across all reconstruction panels
line_w       = 2;               % (a)/(b) curve width (matches the reference figure)

%% ==================== Combined figure ====================
% 6x2 tiled grid: (a)/(b) span the top two rows; each reconstruction spans a
% full-width row below. Explicit tile indices keep the spans unambiguous.
fig = figure('Color', 'w', 'Position', [100 100 750 788]);   % 75% of the previous 1000x1050
tl = tiledlayout(6, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
% Leave headroom above the tiles so the (a)/(b) panel letters can sit clearly
% above the axes (and their y-labels) without clipping at the figure top.
tl.OuterPosition = [0 0 1 0.95];

% (a) Cumulative MC vs delay
ax_a = nexttile(1, [2 1]); hold on; grid off; box off;
for i = 1:n_cond
    plot(delay_s, cumsum(R2{i}), '-', 'Color', colors(i, :), 'LineWidth', line_w);
end
xlim([0 xmax_s]);
xlabel('Delay (s)'); ylabel('Cumulative Memory Capacity');
set(ax_a, 'XTick', [0 5 10], 'YTick', [0 4 8]);   % match plot_memory_capacity_combined.m
hold off;

% (b) Per-delay R^2 vs delay -- carries the one legend for the figure
ax_b = nexttile(2, [2 1]); hold on; grid off; box off;
for i = 1:n_cond
    plot(delay_s, R2{i}, '-', 'Color', colors(i, :), 'LineWidth', line_w);
end
xlim([0 xmax_s]);
xlabel('Delay (s)'); ylabel('$R^2$', 'Interpreter', 'latex');
set(ax_b, 'XTick', [0 5 10], 'YTick', [0 0.5 1]);   % match plot_memory_capacity_combined.m
legend(display_names, 'Location', 'northeast', 'Box', 'off');
hold off;

% Reconstruction rows (SFA+STD): target u(t-d) vs trained readout, delay in s.
% The x-axis (time) is hidden on every row -- a 15 s scale bar on the bottom row
% conveys the timescale instead.
mcr        = results{recon_cond};
tile_ids   = [5, 7, 9, 11];   % top-left tile of each full-width reconstruction row
nR         = numel(recon_delays);
scalebar_s = 15;              % length of the reconstruction scale bar (s)
recon_ax   = gobjects(1, nR);
for k = 1:nR
    d = recon_delays(k);
    p = mcr.predictions(d);
    t_delay = mcr.t_pred(p.t_indices);

    recon_ax(k) = nexttile(tile_ids(k), [1 2]); hold on; grid off; box off;
    plot(t_delay, p.y_true, 'k-', 'LineWidth', 0.9);
    plot(t_delay, p.y_pred, '-', 'Color', colors(recon_cond, :), 'LineWidth', 1.4);
    ylim(recon_ylim);
    title(sprintf('delay = %.1f seconds,  R^2 = %.2f', delay_s(d), mcr.R2_d(d)), ...
        'FontWeight', 'normal', 'Interpreter', 'tex');
    % Hide the x-axis entirely (no ticks, no black axis line).
    set(gca, 'XColor', 'none');
    if k == nR
        % 15 s scale bar in the lower-right corner.
        xl = get(gca, 'XLim');
        x2 = xl(2) - 0.03 * diff(xl);
        x1 = x2 - scalebar_s;
        yb = recon_ylim(1) + 0.05 * diff(recon_ylim);
        plot([x1 x2], [yb yb], 'k-', 'LineWidth', 2.5);
        text((x1 + x2) / 2, yb, sprintf('%d s', scalebar_s), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 12);
    end
    hold off;
end

% Panel letters: (a)/(b) on the analysis panels and (c) on the first
% reconstruction row, same style as the reference figure. The remaining
% reconstruction rows are self-labeled by their delay/R^2 titles.
% Per-axis VShift: push (a)/(b) further up into the headroom (clearing their
% y-labels); keep (c) where it is.
AddLetters2Plots({ax_a, ax_b, recon_ax(1)}, {'(a)', '(b)', '(c)'}, ...
    'FontSize', 18, 'FontWeight', 'normal', 'HShift', -0.06, 'VShift', [-0.075 -0.075 -0.04]);

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
