close all
clc

% Fig_param_space.m
% Stability_Manuscript presentation figure: parameter-space distributions.
% Recreates the param-space histogram figures (LLE + mean_rate distributions)
% from a saved run_all_<dt> run and combines them into a single 2x4 figure
% (row 1 = LLE, row 2 = mean firing rate; one column per adaptation condition).
% No simulation is re-run -- edit data_root to point at a different run_all run.

this_dir     = fileparts(mfilename('fullpath'));
% .../Stability_Manuscript/fig_param_space -> project root is 4 up
project_root = fileparts(fileparts(fileparts(fileparts(this_dir))));

% Resolve the replot pipeline (scripts/run_all_analyses/replot) + src helpers
% (save_some_figs_to_folder_2, capture_git_provenance).
setup_paths();

% Source run (a run_all_<dt> folder with a param_space_* subdir).
data_root = fullfile(project_root, 'data', 'param_space', 'run_all_jul_06_26_22_00');
out_dir   = this_dir;   % write the final figure next to this script

% Start from a clean slate: replot_param_space_analysis saves ALL open figures,
% so any stray figure lingering in the session would pollute the save.
close all force;

% 1) Regenerate the LLE + mean_rate distribution figures into a
%    replot_param_space_<dt>/figures/ folder under data_root.
replot_dir = replot_param_space_analysis(data_root);

% 2) Re-open both saved figures (invisible) and map by Name.
lle_fig = gobjects(0);
mr_fig  = gobjects(0);
fig_listing = dir(fullfile(replot_dir, 'figures', '*.fig'));
for k = 1:numel(fig_listing)
    f = openfig(fullfile(fig_listing(k).folder, fig_listing(k).name), 'invisible');
    switch get(f, 'Name')
        case 'LLE Distribution',       lle_fig = f;
        case 'mean_rate Distribution', mr_fig  = f;
        otherwise,                     close(f);
    end
end

% Row axes, sorted left-to-right (one per condition).
lle_ax = sort_axes_left_to_right(lle_fig);
mr_ax  = sort_axes_left_to_right(mr_fig);

% 3) Build the combined 2x4 figure: row 1 = LLE, row 2 = mean rate. Copy each
%    source axis into a subplot placeholder (same approach as
%    assemble_sensitivity_figure).
nCols = numel(lle_ax);
nRows = 2;
src   = {lle_ax, mr_ax};
combined = figure('Color', 'w', 'Position', [100 100 350*nCols 300*nRows]);
cax = gobjects(nRows, nCols);
for r = 1:nRows
    for c = 1:nCols
        ph = subplot(nRows, nCols, (r-1)*nCols + c, 'Parent', combined);
        target_pos = get(ph, 'Position');
        delete(ph);
        cax(r, c) = copyobj(src{r}(c), combined);
        set(cax(r, c), 'Position', target_pos);
    end
end
close(lle_fig);
close(mr_fig);

% 4) Clean up: fonts matched to the MC/sensitivity figures, condition titles
%    only on the top row (not bold), y-axes linked within each row.
tick_fs  = 14;    % tick numbers -- matches MC/sensitivity figures
label_fs = 15.4;  % axis labels  -- matches MC/sensitivity figures
title_fs = 14;
for r = 1:nRows
    for c = 1:nCols
        ax = cax(r, c);
        set(ax, 'FontSize', tick_fs);
        set(ax.YLabel, 'FontSize', label_fs);
        if r == 1
            xlabel(ax, 'Growth Rate', 'FontSize', label_fs);   % LLE (lambda_1) -> Growth Rate
            set(ax.Title, 'FontWeight', 'normal', 'FontSize', title_fs);  % titles, not bold
        else
            set(ax.XLabel, 'FontSize', label_fs);   % keep 'Mean Firing Rate'
            title(ax, '');   % condition titles only on the top row
        end
    end
    linkaxes(cax(r, :), 'y');   % shared probability axis within each row
end

% 5) Vertical gray dividers between the 4 condition columns (span both rows).
pos = cell2mat(get(cax(:), 'Position'));          % [left bottom width height]
[~, ~, col_of] = uniquetol(pos(:,1), 0.01);       % column index per axis (by left edge)
ncol      = max(col_of);
col_left  = accumarray(col_of, pos(:,1),          [ncol 1], @mean);
col_right = accumarray(col_of, pos(:,1)+pos(:,3), [ncol 1], @mean);
[col_left, ord] = sort(col_left);
col_right = col_right(ord);
y_bot = min(pos(:,2));
y_top = max(pos(:,2) + pos(:,4));
x_shift = 0.007;   % nudge dividers slightly left (normalized figure units)
for c = 1:ncol-1
    x_div = (col_right(c) + col_left(c+1)) / 2 - x_shift;
    annotation(combined, 'line', [x_div x_div], [y_bot y_top], ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
end

% 6) Save ONLY the combined figure, with a STABLE name (save_some_figs_to_folder_2
%    suffixes with the figure number; save then rename). Clearing Fig_ParamSpace*
%    first also removes the earlier separate LLE/mean_rate reproductions.
base = 'Fig_ParamSpace';
old = dir(fullfile(out_dir, [base '*']));
for a = 1:numel(old)
    fp = fullfile(old(a).folder, old(a).name);
    if ~old(a).isdir && (endsWith(fp, '.png') || endsWith(fp, '.svg') || endsWith(fp, '.fig'))
        delete(fp);
    end
end
save_some_figs_to_folder_2(out_dir, base, combined.Number, {'png', 'svg', 'fig'});
num = num2str(combined.Number);
movefile(fullfile(out_dir, [base '_figure_' num '.png']), fullfile(out_dir, [base '.png']));
movefile(fullfile(out_dir, [base '_figure_' num '.svg']), fullfile(out_dir, [base '.svg']));
movefile(fullfile(out_dir, [base '_f_' num '.fig']),      fullfile(out_dir, [base '.fig']));

% 7) Prep figures exist only to build the final -- remove the whole replot
%    folder so no extra figs are left behind in the data dir.
if isfolder(replot_dir)
    rmdir(replot_dir, 's');
end

% 8) Log the git state alongside the figure.
capture_git_provenance(out_dir, project_root);

%% -------------------- Human-readable description --------------------
ps_dirs = dir(fullfile(data_root, 'param_space_*'));
ps_dirs = ps_dirs([ps_dirs.isdir]);

desc_path = fullfile(out_dir, 'README_fig_param_space.txt');
fid = fopen(desc_path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Stability_Manuscript figure: Parameter-space distributions (combined)\n');
fprintf(fid, '====================================================================\n\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now')));
fprintf(fid, 'By script: %s.m\n\n', mfilename);

fprintf(fid, 'HOW IT WAS MADE\n');
fprintf(fid, '  Presentation replot -- no simulation is re-run. The script reloads the\n');
fprintf(fid, '  saved param-space PSA object from a run_all_<dt> run, calls\n');
fprintf(fid, '  replot_param_space_analysis (psa.plot for LLE + mean_rate), then copies\n');
fprintf(fid, '  the per-condition histogram axes into a single 2x4 figure:\n');
fprintf(fid, '    row 1 = LLE (growth-rate) distributions (green dashed zero line)\n');
fprintf(fid, '    row 2 = mean firing-rate distributions\n');
fprintf(fid, '    columns = No Adaptation, SFA, STD, SFA+STD\n');
fprintf(fid, '  Cleanups: LLE row xlabel "LLE (lambda_1)" -> "Growth Rate"; condition\n');
fprintf(fid, '  titles only on the top row (not bold); vertical gray column dividers;\n');
fprintf(fid, '  fonts matched to the MC/sensitivity figures; y-axes linked within each\n');
fprintf(fid, '  row. See git_provenance.txt for the exact commit.\n\n');

fprintf(fid, 'SOURCE RUN\n');
fprintf(fid, '  %s\n', data_root);
fprintf(fid, '  param_space subfolder(s) used:\n');
for k = 1:numel(ps_dirs)
    fprintf(fid, '    %s\n', ps_dirs(k).name);
end
fprintf(fid, '\n');

fprintf(fid, 'FIGURE PRODUCED (in this folder)\n');
fprintf(fid, '  Fig_ParamSpace.png / .svg / .fig\n');

clear cleanup;  % flush + close
fprintf('Description written: %s\n', desc_path);

%% ==================== Local helper ====================
function ax_sorted = sort_axes_left_to_right(fig)
% Return a figure's axes ordered left-to-right by their x-position.
    ax = findobj(fig, 'Type', 'axes');
    p = cell2mat(get(ax, 'Position'));
    [~, order] = sort(p(:, 1));
    ax_sorted = ax(order);
end
