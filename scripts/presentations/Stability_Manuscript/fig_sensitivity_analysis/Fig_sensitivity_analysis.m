close all
clc

% Fig_sensitivity_analysis.m
% Stability_Manuscript presentation figure: combined LLE sensitivity.
% Regenerates the stacked 1D-sensitivity LLE figure from a saved run_all_<dt>
% run and writes the final figure into this presentation folder. No simulation
% is re-run -- edit data_root to point at a different run_all output.

this_dir     = fileparts(mfilename('fullpath'));
% .../Stability_Manuscript/fig_sensitivity_analysis -> project root is 4 up
project_root = fileparts(fileparts(fileparts(fileparts(this_dir))));

% Resolve the replot pipeline (scripts/run_all_analyses/replot) + src helpers
% (save_some_figs_to_folder_2, capture_git_provenance).
addpath(genpath(fullfile(project_root, 'scripts')));
addpath(genpath(fullfile(project_root, 'src')));

% Source run (a run_all_<dt> folder with 1D_sensitivity_* subdirs).
data_root = fullfile(project_root, 'data', 'param_space', 'run_all_jul_06_26_22_00');
out_dir   = this_dir;   % write the final figure next to this script

% Start from a clean slate: replot_sensitivity saves ALL open figures, so any
% stray figure lingering in the session (e.g. a previous combined figure) would
% pollute the per-param save and break assemble_sensitivity_figure.
close all force;

% LLE (growth-rate) y-axis range: re-bin the histograms to this range.
lle_range = [-1.75, 1.75];

% 1) Regenerate per-param LLE figs, 2) assemble the combined figure.
replot_dir = replot_sensitivity(data_root, lle_range);
assemble_sensitivity_figure(replot_dir, 'LLE');

% Re-open the saved combined figure and save it into the presentation folder
% (png + svg + fig), matching the memory-capacity figure's formats/naming.
combined_fig_path = fullfile(replot_dir, 'figures', 'sensitivity_LLE_combined.fig');
cf = openfig(combined_fig_path, 'visible');

% --- Clean up axis labels + fonts -----------------------------------------
% Enlarge tick numbers, relabel the ylabel (lambda_1 -> "Growth Rate", present
% only on the leftmost column), relabel each row's x-axis, and keep the
% condition titles only on the top row (removed from the lower two). The combined
% figure stacks rows alphabetically (f, level_of_chaos, n); rows are identified
% by their x-limits (f in [0.25,0.75], level_of_chaos in [0.5,3], n in [100,1000]):
%   f              -> "E:I ratio"     (E:I = f:(1-f)), ratio tick labels
%   level_of_chaos -> "Synaptic Gain" (g; edge of chaos ~1)
%   n              -> "Network Size"
ei_ticks  = [0.25, 1/3, 0.4, 0.5, 0.6, 2/3, 0.75];
ei_labels = {'1:3', '1:2', '2:3', '1:1', '3:2', '2:1', '3:1'};
tick_fs   = 14;   % larger tick numbers on both axes
ax_all = findobj(cf, 'Type', 'axes');
for a = 1:numel(ax_all)
    ax = ax_all(a);
    set(ax, 'FontSize', tick_fs);         % enlarge x & y tick numbers

    % ylabel (present only on the leftmost column): lambda_1 -> "Growth Rate"
    yl = get(ax, 'YLabel');
    if ~isempty(get(yl, 'String'))
        set(yl, 'String', 'Growth Rate', 'Interpreter', 'none', 'FontSize', 15);
    end

    xmax = max(xlim(ax));
    if xmax <= 1                          % f row -> E:I ratio
        xlabel(ax, 'E:I ratio', 'FontSize', 14);
        set(ax, 'XTick', ei_ticks, 'XTickLabel', ei_labels);
    elseif xmax <= 10                     % level_of_chaos row -> Synaptic Gain
        xlabel(ax, 'Synaptic Gain', 'FontSize', 14);
        title(ax, '');                    % condition titles only on the top row
    else                                  % n row -> Network Size
        xlabel(ax, 'Network Size', 'FontSize', 14);
        title(ax, '');                    % condition titles only on the top row
    end
end

% --- Vertical dividers between the 4 condition columns --------------------
% Draw 3 black lines spanning the plot height in the gaps between adjacent
% columns, to visually separate the conditions. Column boundaries come from the
% axes' normalized Positions; each line sits at the midpoint of the inter-column
% gap (so it falls in the gutter, not over any axes).
pos = cell2mat(get(ax_all, 'Position'));         % [left bottom width height]
[~, ~, col_of] = uniquetol(pos(:,1), 0.01);      % column index per axis (by left edge)
ncol      = max(col_of);
col_left  = accumarray(col_of, pos(:,1),          [ncol 1], @mean);
col_right = accumarray(col_of, pos(:,1)+pos(:,3), [ncol 1], @mean);
[col_left, ord] = sort(col_left);
col_right = col_right(ord);
y_bot = min(pos(:,2));                            % bottom of the bottom row
y_top = max(pos(:,2) + pos(:,4));                 % top of the top row
x_shift = 0.007;   % nudge dividers slightly left (normalized figure units)
for c = 1:ncol-1
    x_div = (col_right(c) + col_left(c+1)) / 2 - x_shift;
    annotation(cf, 'line', [x_div x_div], [y_bot y_top], ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
end

% --- Save ONLY the cleaned combined figure, with a STABLE name ------------
% save_some_figs_to_folder_2 suffixes filenames with the (run-dependent) figure
% number; save, then rename to fixed names so re-runs overwrite cleanly. First
% clear any prior Fig_Sensitivity_LLE outputs (stable or numbered) so nothing
% stale lingers.
fig_tag = 'Fig_Sensitivity_LLE';
old = dir(fullfile(out_dir, [fig_tag '*']));
for a = 1:numel(old)
    fp = fullfile(old(a).folder, old(a).name);
    if ~old(a).isdir && (endsWith(fp, '.png') || endsWith(fp, '.svg') || endsWith(fp, '.fig'))
        delete(fp);
    end
end
save_some_figs_to_folder_2(out_dir, fig_tag, cf.Number, {'png', 'svg', 'fig'});
num = num2str(cf.Number);
movefile(fullfile(out_dir, [fig_tag '_figure_' num '.png']), fullfile(out_dir, [fig_tag '.png']));
movefile(fullfile(out_dir, [fig_tag '_figure_' num '.svg']), fullfile(out_dir, [fig_tag '.svg']));
movefile(fullfile(out_dir, [fig_tag '_f_' num '.fig']),      fullfile(out_dir, [fig_tag '.fig']));

% Prep figures (per-param LLE + mean_rate, and the intermediate combined) exist
% only to build the final figure -- remove the whole replot folder so no extra
% figs are left behind in the data dir.
if isfolder(replot_dir)
    rmdir(replot_dir, 's');
end

% Log the git state alongside the figure so this presentation output can be
% traced back to an exact commit (+ working-tree diff if dirty).
capture_git_provenance(out_dir, project_root);

%% -------------------- Human-readable description --------------------
% Write a plain-text record of how this figure was produced: the source run,
% the swept-parameter subfolders used, the pipeline, and the output filenames.
fig_files = { ...
    'Fig_Sensitivity_LLE.png'; ...
    'Fig_Sensitivity_LLE.svg'; ...
    'Fig_Sensitivity_LLE.fig' };

sens_dirs = dir(fullfile(data_root, '1D_sensitivity_*'));
sens_dirs = sens_dirs([sens_dirs.isdir]);

desc_path = fullfile(out_dir, 'README_fig_sensitivity_analysis.txt');
fid = fopen(desc_path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Stability_Manuscript figure: LLE Sensitivity (combined)\n');
fprintf(fid, '=======================================================\n\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now')));
fprintf(fid, 'By script: %s.m\n\n', mfilename);

fprintf(fid, 'HOW IT WAS MADE\n');
fprintf(fid, '  Presentation replot -- no simulation is re-run. The script reloads the\n');
fprintf(fid, '  saved 1D-sensitivity PSA objects from a run_all_<dt> run and calls\n');
fprintf(fid, '  replot_sensitivity -> assemble_sensitivity_figure to rebuild the stacked\n');
fprintf(fid, '  LLE figure (rows = swept params, cols = adaptation conditions), then\n');
fprintf(fid, '  saves it here. See git_provenance.txt for the exact commit.\n\n');

fprintf(fid, 'SOURCE RUN\n');
fprintf(fid, '  %s\n', data_root);
fprintf(fid, '  1D_sensitivity subfolders used:\n');
for k = 1:numel(sens_dirs)
    fprintf(fid, '    %s\n', sens_dirs(k).name);
end
fprintf(fid, '\n');

fprintf(fid, 'FIGURES PRODUCED (in this folder)\n');
for k = 1:numel(fig_files)
    fprintf(fid, '  %s\n', fig_files{k});
end
fprintf(fid, '\n  Combined LLE sensitivity: one row per swept parameter\n');
fprintf(fid, '  (f, level_of_chaos, n), one column per adaptation condition.\n');
fprintf(fid, '  x-axes relabelled: f -> "E:I ratio" (E:I = f:(1-f), ticks\n');
fprintf(fid, '  1:3, 1:2, 2:3, 1:1, 3:2, 2:1, 3:1); level_of_chaos -> "Synaptic Gain";\n');
fprintf(fid, '  n -> "Network Size". ylabel lambda_1 -> "Growth Rate"; larger tick fonts.\n');
fprintf(fid, '  Condition titles kept only on the top row; vertical black dividers\n');
fprintf(fid, '  separate the 4 condition columns.\n');

clear cleanup;  % flush + close
fprintf('Description written: %s\n', desc_path);
