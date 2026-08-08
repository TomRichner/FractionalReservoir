close all
clc

% Fig_sensitivity_analysis.m
% Stability_Manuscript presentation figures: combined LLE + mean_rate sensitivity.
% Regenerates the stacked 1D-sensitivity figures from a saved run_all_<dt>
% run and writes the final figures into this presentation folder. No simulation
% is re-run -- edit data_root to point at a different run_all output.

this_dir     = fileparts(mfilename('fullpath'));
% .../Stability_Manuscript/fig_sensitivity_analysis -> project root is 4 up
project_root = fileparts(fileparts(fileparts(fileparts(this_dir))));

% Resolve the replot pipeline (scripts/run_all_analyses/replot) + src helpers
% (save_some_figs_to_folder_2, capture_git_provenance).
setup_paths();

% Source run (a run_all_<dt> folder with 1D_sensitivity_* subdirs).
data_root = fullfile(project_root, 'data', 'param_space', 'run_all_jul_06_26_22_00');
out_dir   = this_dir;   % write the final figures next to this script

% Start from a clean slate: replot_sensitivity saves ALL open figures, so any
% stray figure lingering in the session (e.g. a previous combined figure) would
% pollute the per-param save and break assemble_sensitivity_figure.
close all force;

% LLE y-axis range + histogram bin count. n_bins counts linspace EDGES, so both
% metrics draw n_bins+1 rows (n_bins-1 interior bins plus the two -inf/+inf
% overflow bins) -- 25 rows at n_bins = 24. mean_rate keeps its default [0,1]
% range but is given the same bin count so the two figures match vertically.
lle_range = [-1.75, 1.75];
n_bins    = 24;

% Regenerate the per-param figures. One call produces BOTH the LLE and the
% mean_rate figures for every swept parameter (replot_sensitivity plots both).
replot_dir = replot_sensitivity(data_root, lle_range, n_bins, n_bins);

% --- Shared cleanup constants ---------------------------------------------
% Enlarge tick numbers, enlarge the ylabel (present only on the leftmost
% column), relabel each row's x-axis, and keep the condition titles only on the
% top row (removed from the lower two). The combined figure stacks rows
% alphabetically (f, level_of_chaos, n); rows are identified by their x-limits
% (f in [0.25,0.75], level_of_chaos in [0.5,3], n in [100,1000]):
%   f              -> "E:I ratio"     (E:I = f:(1-f)), ratio tick labels
%   level_of_chaos -> "Synaptic Gain" (g; edge of chaos ~1)
%   n              -> "Network Size"
ei_ticks  = [0.25, 1/3, 0.4, 0.5, 0.6, 2/3, 0.75];
ei_labels = {'1:3', '1:2', '2:3', '1:1', '3:2', '2:1', '3:1'};
tick_fs   = 14;    % tick numbers -- matches MC figure DefaultAxesFontSize
label_fs  = 15.4;  % axis labels -- matches MC figure (14 * 1.1 label multiplier)
title_fs  = 20;    % condition titles (no adaptation, sfa only, ...) -- enlarged
letter_fs = 18;    % panel letters -- matches MC figure
row_shrink = 0.85; % shrink each row's height to open gaps between rows
top_headroom = 0.06; % shift the row stack down (normalized) to clear room above the top row for column headers
title_y   = 1.22;  % condition-title height above the top-row axes (normalized), reads as a column header
% Colormap ramps white (0 counts) -> 90% black (max), not pure black, so the
% blue median line stays visible over the darkest cells.
dark_cmap = repmat(linspace(1, 0.1, 256)', 1, 3);
median_alpha = 0.35;   % blue median line transparency (plot_sensitivity uses 0.55)
median_lw    = 3;      % blue median line width, 25% thinner (plot_sensitivity uses 4)
zeroline_lw  = 2;      % green dashed zero line width (plot_sensitivity uses 4)

% --- One combined figure per metric ----------------------------------------
% assemble_sensitivity_figure matches the per-param figs by their Name prefix
% ("<metric> Sensitivity - "), so the same cleanup serves both metrics. Only the
% ylabel, the zero line and the output name differ:
%   LLE       -> lambda_1, green dashed zero line kept (sign of lambda_1 matters)
%   mean_rate -> "Mean Firing Rate"; zero line dropped (range is [0,1], so
%                yline(0) just sits on the bottom axis and carries no meaning)
% clim_frac darkens imagesc: CLim is capped at total_reps*clim_frac (shared
% across panels of a metric so they stay comparable). Tuned per metric because
% the two use different bin counts/ranges, so counts concentrate differently.
metric_specs = struct( ...
    'name',      {'LLE',                 'mean_rate'}, ...
    'ylabel',    {'\lambda_1',           'Mean Firing Rate'}, ...
    'fig_tag',   {'Fig_Sensitivity_LLE', 'Fig_sensitivity_mean_rate'}, ...
    'zero_line', {true,                  false}, ...
    'clim_frac', {0.4,                   0.4});

for mi = 1:numel(metric_specs)
    spec = metric_specs(mi);

    % Assemble the per-param figs for this metric into one stacked figure, then
    % re-open it so it can be cleaned up and saved into the presentation folder
    % (png + svg + fig), matching the memory-capacity figure's formats/naming.
    assemble_sensitivity_figure(replot_dir, spec.name);
    combined_fig_path = fullfile(replot_dir, 'figures', ...
        sprintf('sensitivity_%s_combined.fig', spec.name));
    cf = openfig(combined_fig_path, 'visible');

    % --- Reshape the whole figure -----------------------------------------
    % Make it ~15% narrower and ~5% taller so the imagesc panels are less wide
    % (the added row spacing / raised headers had stretched their aspect ratio).
    % Axes use normalized positions, so they follow the new figure proportions.
    fig_pos = get(cf, 'Position');            % [left bottom width height]
    set(cf, 'Position', [fig_pos(1), fig_pos(2), fig_pos(3) * 0.85, fig_pos(4) * 1.05]);

    % --- Clean up axis labels + fonts -------------------------------------
    ax_all = findobj(cf, 'Type', 'axes');
    for a = 1:numel(ax_all)
        ax = ax_all(a);
        set(ax, 'FontSize', tick_fs);         % enlarge x & y tick numbers
        set(get(ax, 'Title'), 'FontWeight', 'normal', 'FontSize', title_fs);  % condition titles, not bold, enlarged
        box(ax, 'off');                       % drop the rectangle; keep x/y axes+ticks
        colormap(ax, dark_cmap);              % white -> 90% black (blue line stays visible)

        % Darken the histogram density: the panels are drawn with CLim = [0,
        % total_reps]; lower the ceiling so typical bin counts use more of the
        % dark range (kept shared across panels so they stay comparable).
        cl = get(ax, 'CLim');
        set(ax, 'CLim', [0, cl(2) * spec.clim_frac]);

        % Blue median line: more transparent + 25% thinner. (imagesc is Type
        % 'image', the zero line is 'constantline', so 'line' is the median.)
        ml = findobj(ax, 'Type', 'line');
        for m = 1:numel(ml)
            mc = get(ml(m), 'Color');
            if numel(mc) < 4; mc(4) = 1; end
            mc(4) = median_alpha;
            set(ml(m), 'Color', mc, 'LineWidth', median_lw);
        end
        % Green dashed zero line: thinner for LLE, dropped for mean_rate.
        zl = findobj(ax, 'Type', 'constantline');
        if spec.zero_line
            set(zl, 'LineWidth', zeroline_lw);
        else
            delete(zl);
        end

        % ylabel (present only on the leftmost column): use this metric's label.
        % Interpreter must be 'tex' so \lambda_1 renders as the symbol.
        yl = get(ax, 'YLabel');
        if ~isempty(get(yl, 'String'))
            set(yl, 'String', spec.ylabel, 'Interpreter', 'tex', 'FontSize', label_fs);
        end

        xmax = max(xlim(ax));
        if xmax <= 1                          % f row -> E:I ratio
            xlabel(ax, 'E:I ratio', 'FontSize', label_fs);
            set(ax, 'XTick', ei_ticks, 'XTickLabel', ei_labels);
        elseif xmax <= 10                     % level_of_chaos row -> Synaptic Gain
            xlabel(ax, 'Synaptic Gain', 'FontSize', label_fs);
            title(ax, '');                    % condition titles only on the top row
        else                                  % n row -> Network Size
            xlabel(ax, 'Network Size', 'FontSize', label_fs);
            set(ax, 'XTick', [100, 500, 1000]);   % fewer network-size ticks
            title(ax, '');                    % condition titles only on the top row
        end
    end

    % --- Open more vertical space between rows ----------------------------
    % Shrink each axis's height (widening the gaps below each row) and slide the
    % whole stack down by top_headroom, freeing space above the top row so the
    % condition titles can float well above it as column headers. Done before the
    % divider/letter code so those pick up the new positions.
    for a = 1:numel(ax_all)
        ax = ax_all(a);
        p  = get(ax, 'Position');             % [left bottom width height]
        new_h = p(4) * row_shrink;
        set(ax, 'Position', [p(1), p(2) + (p(4) - new_h) - top_headroom, p(3), new_h]);
    end

    % --- Lift the condition titles into column-header position ------------
    % Only the top-row axes still carry a title string (the loop above cleared
    % the lower rows). Raise each well above its axes and enlarge it so it
    % clearly labels the whole column, not just the first row.
    for a = 1:numel(ax_all)
        ax = ax_all(a);
        t  = get(ax, 'Title');
        if ~isempty(get(t, 'String'))
            set(t, 'Units', 'normalized', 'Position', [0.5, title_y, 0], ...
                'VerticalAlignment', 'bottom', 'FontSize', title_fs);
        end
    end

    % --- Vertical dividers between the 4 condition columns ----------------
    % Draw 3 black lines spanning the plot height in the gaps between adjacent
    % columns, to visually separate the conditions. Column boundaries come from
    % the axes' normalized Positions; each line sits at the midpoint of the
    % inter-column gap (so it falls in the gutter, not over any axes).
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

    % --- Panel letters (a), (b), ... up-and-left of each plot --------------
    % AddLetters2Plots sorts axes left-to-right, top-to-bottom, so (a) is the
    % top-left panel. Negative HShift/VShift push each label outside the axes, up
    % and to the left of its top-left corner.
    panel_letters = arrayfun(@(ch) sprintf('(%c)', ch), ...
        char('a' + (0:numel(ax_all)-1)), 'UniformOutput', false);
    AddLetters2Plots(cf, panel_letters, ...
        'FontSize', letter_fs, 'FontWeight', 'normal', 'HShift', -0.03, 'VShift', -0.04);

    % --- Save ONLY the cleaned combined figure, with a STABLE name ---------
    % save_some_figs_to_folder_2 suffixes filenames with the (run-dependent)
    % figure number; save, then rename to fixed names so re-runs overwrite
    % cleanly. First clear any prior outputs for THIS metric (stable or
    % numbered) so nothing stale lingers. Saving by cf.Number means the other
    % metric's figure, still open, is untouched.
    fig_tag = spec.fig_tag;
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
end

% Prep figures (per-param LLE + mean_rate, and the intermediate combined ones)
% exist only to build the final figures -- remove the whole replot folder so no
% extra figs are left behind in the data dir. Must come after BOTH metrics.
if isfolder(replot_dir)
    rmdir(replot_dir, 's');
end

% Log the git state alongside the figures so this presentation output can be
% traced back to an exact commit (+ working-tree diff if dirty).
capture_git_provenance(out_dir, project_root);

%% -------------------- Human-readable description --------------------
% Write a plain-text record of how these figures were produced: the source run,
% the swept-parameter subfolders used, the pipeline, and the output filenames.
fig_files = { ...
    'Fig_Sensitivity_LLE.png'; ...
    'Fig_Sensitivity_LLE.svg'; ...
    'Fig_Sensitivity_LLE.fig'; ...
    'Fig_sensitivity_mean_rate.png'; ...
    'Fig_sensitivity_mean_rate.svg'; ...
    'Fig_sensitivity_mean_rate.fig' };

sens_dirs = dir(fullfile(data_root, '1D_sensitivity_*'));
sens_dirs = sens_dirs([sens_dirs.isdir]);

desc_path = fullfile(out_dir, 'README_fig_sensitivity_analysis.txt');
fid = fopen(desc_path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Stability_Manuscript figures: LLE + mean firing rate Sensitivity (combined)\n');
fprintf(fid, '==========================================================================\n\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now')));
fprintf(fid, 'By script: %s.m\n\n', mfilename);

fprintf(fid, 'HOW THEY WERE MADE\n');
fprintf(fid, '  Presentation replot -- no simulation is re-run. The script reloads the\n');
fprintf(fid, '  saved 1D-sensitivity PSA objects from a run_all_<dt> run and calls\n');
fprintf(fid, '  replot_sensitivity (which plots BOTH metrics) -> assemble_sensitivity_figure\n');
fprintf(fid, '  once per metric to rebuild the stacked figures (rows = swept params,\n');
fprintf(fid, '  cols = adaptation conditions), then saves them here. See\n');
fprintf(fid, '  git_provenance.txt for the exact commit.\n\n');

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

fprintf(fid, '\nSHARED LAYOUT (both figures)\n');
fprintf(fid, '  One row per swept parameter (f, level_of_chaos, n), one column per\n');
fprintf(fid, '  adaptation condition. x-axes relabelled: f -> "E:I ratio"\n');
fprintf(fid, '  (E:I = f:(1-f), ticks 1:3, 1:2, 2:3, 1:1, 3:2, 2:1, 3:1);\n');
fprintf(fid, '  level_of_chaos -> "Synaptic Gain"; n -> "Network Size". Larger tick\n');
fprintf(fid, '  fonts. Condition titles kept only on the top row; vertical gray\n');
fprintf(fid, '  dividers separate the 4 condition columns. imagesc CLim capped at\n');
fprintf(fid, '  total_reps*0.4 (shared within a figure); colormap white -> 90%% black\n');
fprintf(fid, '  so the blue median line stays visible over the darkest cells. Panel\n');
fprintf(fid, '  letters (a)..(l) added up-and-left of each plot (AddLetters2Plots).\n');
fprintf(fid, '  Blue median line: alpha 0.35, 25%% thinner. Titles not bold; axis\n');
fprintf(fid, '  boxes removed (x/y axes + ticks kept).\n');

fprintf(fid, '\nPER-FIGURE DIFFERENCES\n');
fprintf(fid, '  Both metrics use n_bins = 24 (linspace edges), i.e. 25 plotted rows:\n');
fprintf(fid, '  23 interior bins plus the two -inf/+inf overflow bins. Matching the\n');
fprintf(fid, '  bin count gives the two figures the same vertical resolution.\n\n');
fprintf(fid, '  Fig_Sensitivity_LLE:       ylabel \\lambda_1; histogram range\n');
fprintf(fid, '                             [-1.75, 1.75]; green dashed zero line kept\n');
fprintf(fid, '                             (thinner) -- the sign of lambda_1 marks the\n');
fprintf(fid, '                             edge of chaos. The solid bands at the top and\n');
fprintf(fid, '                             bottom are the overflow bins (reps outside\n');
fprintf(fid, '                             the range).\n');
fprintf(fid, '  Fig_sensitivity_mean_rate: ylabel "Mean Firing Rate"; histogram range\n');
fprintf(fid, '                             [0, 1] (plot_sensitivity default; nothing can\n');
fprintf(fid, '                             fall outside it, so the overflow bins are\n');
fprintf(fid, '                             always empty); zero line removed -- at y=0 it\n');
fprintf(fid, '                             lands on the bottom axis and carries no\n');
fprintf(fid, '                             meaning for a rate.\n');

clear cleanup;  % flush + close
fprintf('Description written: %s\n', desc_path);
