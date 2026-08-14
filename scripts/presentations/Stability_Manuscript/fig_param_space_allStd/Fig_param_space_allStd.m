close all
clc

% Fig_param_space_allStd.m
% Stability_Manuscript presentation figure: parameter-space distributions, for the
% SRNNCellTypePairs "allStd" runs. Recreates the param-space histogram figures
% (LLE + mean_rate distributions) from a saved run_all_<dt> run and combines them
% into a single 2xN figure (row 1 = LLE, row 2 = mean firing rate; one column per
% adaptation condition). No simulation is re-run -- edit data_root below.
%
% Sibling of ../fig_param_space/Fig_param_space.m, which is left untouched. This
% is the most direct of the three ports: replot_param_space_analysis already
% works on SRNNCellTypePairs unchanged, because ParamSpaceAnalysis2.plot reads
% only 'metric' and 'pool_with' -- it takes no color_by, so the Pairs-vs-SRNNModel2
% colouring caveat that applies to plot_unit_histograms does not apply here.
%
% Note the LLE histogram range is fixed at [-1.5, 1.5] inside
% ParamSpaceAnalysis2.plot and is not settable from here. On this preset the
% param-space LLEs span roughly -10 to +4.8, so the outermost bins carry a large
% share of the distribution.
%
% See also: Fig_sensitivity_analysis_allStd, Fig_sfa_EOC_allStd

this_dir = fileparts(mfilename('fullpath'));
% Depth-independent project root (CLAUDE.md convention).
setup_paths();
project_root = fileparts(which('setup_paths'));

% Source run (a run_all_<dt> folder with a param_space_* subdir).
% Swap this one line to regenerate against the medium run.
data_root = fullfile(project_root, 'data', 'param_space', 'run_all_aug_14_26_12_14');
out_dir   = this_dir;   % write the final figure next to this script

% --- Presentation constants ------------------------------------------------
tick_fs  = 14;    % tick numbers -- matches MC/sensitivity figures
label_fs = 15.4;  % axis labels  -- matches MC/sensitivity figures
title_fs = 14;
% Probability y ticks, one entry per row. The two rows span different ranges
% (the LLE distributions peak near 0.8, the rate distributions near 0.4), so
% each gets its own sparse set rather than MATLAB's automatic ticks.
prob_yticks = {[0, 0.4, 0.8], [0, 0.2, 0.4]};
x_shift = 0.007;   % nudge column dividers slightly left (normalized figure units)

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
if isempty(lle_fig) || isempty(mr_fig)
    error('Fig_param_space_allStd:MissingFigure', ...
        ['Expected both an "LLE Distribution" and a "mean_rate Distribution" ' ...
         'figure in:\n  %s'], fullfile(replot_dir, 'figures'));
end

% Row axes, sorted left-to-right (one per condition).
lle_ax = sort_axes_left_to_right(lle_fig);
mr_ax  = sort_axes_left_to_right(mr_fig);

% 3) Build the combined 2xN figure: row 1 = LLE, row 2 = mean rate. Copy each
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
for r = 1:nRows
    for c = 1:nCols
        ax = cax(r, c);
        set(ax, 'FontSize', tick_fs);
        set(ax.YLabel, 'FontSize', label_fs);
        if r == 1
            % LLE is the x-axis here (row 1 is a distribution). 'tex' so the
            % symbol renders rather than printing the literal \lambda_1.
            xlabel(ax, '\lambda_1', 'Interpreter', 'tex', 'FontSize', label_fs);
            set(ax.Title, 'FontWeight', 'normal', 'FontSize', title_fs);  % titles, not bold
        else
            set(ax.XLabel, 'FontSize', label_fs);   % keep 'Mean Firing Rate'
            title(ax, '');   % condition titles only on the top row
        end
    end
    linkaxes(cax(r, :), 'y');   % shared probability axis within each row
    % After linkaxes, so the shared limits don't regenerate automatic ticks.
    set(cax(r, :), 'YTick', prob_yticks{r});
end

% 5) Vertical gray dividers between the condition columns (span both rows).
pos = cell2mat(get(cax(:), 'Position'));          % [left bottom width height]
[~, ~, col_of] = uniquetol(pos(:,1), 0.01);       % column index per axis (by left edge)
ncol      = max(col_of);
col_left  = accumarray(col_of, pos(:,1),          [ncol 1], @mean);
col_right = accumarray(col_of, pos(:,1)+pos(:,3), [ncol 1], @mean);
[col_left, ord] = sort(col_left);
col_right = col_right(ord);
y_bot = min(pos(:,2));
y_top = max(pos(:,2) + pos(:,4));
for c = 1:ncol-1
    x_div = (col_right(c) + col_left(c+1)) / 2 - x_shift;
    annotation(combined, 'line', [x_div x_div], [y_bot y_top], ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
end

% 6) Save ONLY the combined figure, with a STABLE name.
fig_tag = 'Fig_ParamSpace_allStd';
save_figure_stable(out_dir, fig_tag, combined);

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

desc_path = fullfile(out_dir, 'README_fig_param_space_allStd.txt');
fid = fopen(desc_path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Stability_Manuscript figure: Parameter-space distributions (allStd)\n');
fprintf(fid, '==================================================================\n\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now')));
fprintf(fid, 'By script: %s.m\n\n', mfilename);

fprintf(fid, 'HOW IT WAS MADE\n');
fprintf(fid, '  Presentation replot -- no simulation is re-run. The script reloads the\n');
fprintf(fid, '  saved param-space PSA object from a run_all_<dt> run, calls\n');
fprintf(fid, '  replot_param_space_analysis (psa.plot for LLE + mean_rate), then copies\n');
fprintf(fid, '  the per-condition histogram axes into a single 2x%d figure:\n', nCols);
fprintf(fid, '    row 1 = LLE distributions (green dashed zero line)\n');
fprintf(fid, '    row 2 = mean firing-rate distributions\n');
fprintf(fid, '    columns = No Adaptation, SFA, STD, SFA+STD\n');
fprintf(fid, '  Cleanups: LLE row xlabel "LLE (lambda_1)" -> \\lambda_1; condition\n');
fprintf(fid, '  titles only on the top row (not bold); vertical gray column dividers;\n');
fprintf(fid, '  fonts matched to the MC/sensitivity figures; y-axes linked within each\n');
fprintf(fid, '  row, with probability ticks at 0/0.4/0.8 (LLE row) and 0/0.2/0.4\n');
fprintf(fid, '  (rate row). See git_provenance.txt for the exact commit.\n\n');

fprintf(fid, 'SOURCE RUN\n');
fprintf(fid, '  %s\n', data_root);
fprintf(fid, '  param_space subfolder(s) used:\n');
for k = 1:numel(ps_dirs)
    fprintf(fid, '    %s\n', ps_dirs(k).name);
end
fprintf(fid, '\n');

fprintf(fid, 'FIGURE PRODUCED (in this folder)\n');
fprintf(fid, '  %s.png / .svg / .fig\n\n', fig_tag);

fprintf(fid, 'READING THIS FIGURE\n');
fprintf(fid, '  The LLE histogram range is fixed at [-1.5, 1.5] inside\n');
fprintf(fid, '  ParamSpaceAnalysis2.plot and cannot be set from this script. On this\n');
fprintf(fid, '  preset the param-space LLEs span roughly -10 to +4.8, so the outermost\n');
fprintf(fid, '  bins accumulate a large share of the distribution.\n');

clear cleanup;  % flush + close
fprintf('Description written: %s\n', desc_path);


%% ==================== Local helpers ====================
function ax_sorted = sort_axes_left_to_right(fig)
% Return a figure's axes ordered left-to-right by their x-position.
    ax = findobj(fig, 'Type', 'axes');
    p = cell2mat(get(ax, 'Position'));
    [~, order] = sort(p(:, 1));
    ax_sorted = ax(order);
end

function save_figure_stable(out_dir, fig_tag, fh)
% SAVE_FIGURE_STABLE Save one figure as <fig_tag>.{png,svg,fig} in out_dir.
%
% save_some_figs_to_folder_2 suffixes filenames with the (run-dependent) figure
% number; save, then rename to fixed names so re-runs overwrite cleanly. First
% clear any prior outputs for THIS tag (stable or numbered) so nothing stale
% lingers. Saving by fh.Number means other open figures are untouched.
%
% The renames are GUARDED because save_some_figs_to_folder_2 catches export
% failures, warns, and carries on (see its comment: roughly one run in four, a
% live figure reaches a state where every rasterizing path throws). An
% unguarded movefile would then abort this script before the replot folder is
% cleaned up.
    old = dir(fullfile(out_dir, [fig_tag '*']));
    for a = 1:numel(old)
        fp = fullfile(old(a).folder, old(a).name);
        if ~old(a).isdir && (endsWith(fp, '.png') || endsWith(fp, '.svg') || endsWith(fp, '.fig'))
            delete(fp);
        end
    end

    save_some_figs_to_folder_2(out_dir, fig_tag, fh.Number, {'png', 'svg', 'fig'});

    num = num2str(fh.Number);
    renames = { ...
        [fig_tag '_figure_' num '.png'], [fig_tag '.png']; ...
        [fig_tag '_figure_' num '.svg'], [fig_tag '.svg']; ...
        [fig_tag '_f_' num '.fig'],      [fig_tag '.fig']};
    for r = 1:size(renames, 1)
        src = fullfile(out_dir, renames{r, 1});
        if isfile(src)
            movefile(src, fullfile(out_dir, renames{r, 2}));
        else
            warning('Fig_param_space_allStd:ExportMissing', ...
                'Export did not produce %s; skipping that format for ''%s''.', ...
                renames{r, 1}, fig_tag);
        end
    end
end
