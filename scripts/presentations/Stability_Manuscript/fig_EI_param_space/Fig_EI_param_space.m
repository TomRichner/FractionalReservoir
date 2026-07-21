close all
clc

% Fig_EI_param_space.m
% Stability_Manuscript presentation figure: parameter-space distributions with
% histograms color-coded by the excitatory fraction f (E:I ratio). A 2x5 layout:
% columns 1-4 are the four conditions (row 1 = LLE, row 2 = mean rate); each bar
% is a stack of per-network patches colored by f (blue = inhibition-dominated,
% red = excitation-dominated) via unit_histogram_patch. Column 5 embeds the
% f-gradient colorbar (upper-right cell) encoding f as an E:I ratio; the
% lower-right cell is left empty. No simulation is re-run.

this_dir     = fileparts(mfilename('fullpath'));
% .../Stability_Manuscript/EI_param_space -> project root is 4 up
project_root = fileparts(fileparts(fileparts(fileparts(this_dir))));

% Resolve src helpers (load_and_make_unit_histograms, unit_histogram_patch,
% blue_gray_red_colormap, save_some_figs_to_folder_2, capture_git_provenance).
addpath(genpath(fullfile(project_root, 'scripts')));
addpath(genpath(fullfile(project_root, 'src')));

% Source run + its param_space_* subdir (the colored builder takes the subdir).
data_root = fullfile(project_root, 'data', 'param_space', 'run_all_jul_06_26_22_00');
out_dir   = this_dir;   % write the final figures next to this script

ps_dirs = dir(fullfile(data_root, 'param_space_*'));
ps_dirs = ps_dirs([ps_dirs.isdir]);
assert(~isempty(ps_dirs), 'No param_space_* subdir found in %s', data_root);
ps_dir  = fullfile(ps_dirs(1).folder, ps_dirs(1).name);

% Start from a clean slate.
close all force;

% 1) Build the f-colored per-metric histogram figures (LLE + mean_rate), matching
%    the gray figure's LLE range [-1.5,1.5] and probability normalization. This
%    also creates a separate 'f Value Colorbar' figure.
[~, ~] = load_and_make_unit_histograms(ps_dir, ...
    'Metrics', {'lle', 'r'}, 'NormalizeMode', 'probability', 'LLERange', [-1.5, 1.5]);

% 2) Grab the metric figures (by Name, robust) + the colorbar figure.
lle_fig = findobj(0, 'Type', 'figure', 'Name', 'LLE Unit Histogram');
mr_fig  = findobj(0, 'Type', 'figure', 'Name', 'mean_rate Unit Histogram');
cb_fig  = findobj(0, 'Type', 'figure', 'Name', 'f Value Colorbar');

% Row axes, sorted left-to-right (one per condition).
lle_ax = sort_axes_left_to_right(lle_fig);
mr_ax  = sort_axes_left_to_right(mr_fig);

% 3) Build the combined 2x5 figure: columns 1-4 are the four conditions
%    (row 1 = LLE, row 2 = mean rate); column 5 holds the f-gradient colorbar in
%    the upper-right cell, with the lower-right cell left empty. Copy each source
%    axis into a subplot placeholder (same as Fig_param_space.m).
nCols = numel(lle_ax);        % 4 condition columns
nRows = 2;
nGrid = nCols + 1;            % 5-column layout; column 5 holds the colorbar
src   = {lle_ax, mr_ax};
combined = figure('Color', 'w', 'Position', [100 100 350*nGrid 300*nRows]);
cax = gobjects(nRows, nCols);
for r = 1:nRows
    for c = 1:nCols
        ph = subplot(nRows, nGrid, (r-1)*nGrid + c, 'Parent', combined);
        target_pos = get(ph, 'Position');
        delete(ph);
        cax(r, c) = copyobj(src{r}(c), combined);
        set(cax(r, c), 'Position', target_pos);
    end
end

% Colorbar into the upper-right cell (row 1, column 5 = subplot nGrid). The axes
% keeps its thin pbaspect [0.1 1 1], so it renders as a vertical bar centered in
% the cell; the lower-right cell (subplot 2*nGrid) is simply never filled.
cbax = gobjects(0);
if ~isempty(cb_fig)
    ph = subplot(nRows, nGrid, nGrid, 'Parent', combined);
    cb_target_pos = get(ph, 'Position');
    delete(ph);
    src_cb = findobj(cb_fig, 'Type', 'axes');
    cbax = copyobj(src_cb(1), combined);
    set(cbax, 'Position', cb_target_pos);
end
close(lle_fig);
close(mr_fig);
if ~isempty(cb_fig); close(cb_fig); end

% 4) Clean up: fonts matched to the MC/sensitivity figures, condition titles only
%    on the top row (not bold), y-axes linked within each row.
tick_fs  = 14;    % tick numbers -- matches MC/sensitivity figures
label_fs = 15.4;  % axis labels  -- matches MC/sensitivity figures
title_fs  = 20;   % condition titles -- matches sensitivity figure (column headers)
axes_lw   = 1.0;  % x/y axis line width (default 0.5)
letter_fs = 18;   % panel letters -- matches MC/sensitivity figures
row_shrink   = 0.85; % shrink each row's height to open the gap between rows
top_headroom = 0.06; % shift the row stack down (normalized) to clear room above the top row for column headers
title_y      = 1.22; % condition-title height above the top-row axes (normalized), reads as a column header
lle_yticks   = [0, 0.5];   % row 1 (Growth Rate) y ticks
rate_yticks  = [0, 0.3];  % row 2 (mean firing rate) y ticks
for r = 1:nRows
    for c = 1:nCols
        ax = cax(r, c);
        set(ax, 'FontSize', tick_fs, 'LineWidth', axes_lw);
        set(ax.YLabel, 'FontSize', label_fs);
        if r == 1
            xlabel(ax, 'Growth Rate', 'FontSize', label_fs);   % lambda_1 -> Growth Rate
            set(ax.Title, 'FontWeight', 'normal', 'FontSize', title_fs);  % titles, not bold
            set(findobj(ax, 'Type', 'constantline'), 'Color', [0 0.7 0]);  % zero line -> green
            set(ax, 'YTick', lle_yticks);    % just 0 and 0.5
        else
            set(ax.XLabel, 'FontSize', label_fs);   % keep 'Mean Firing Rate'
            title(ax, '');   % condition titles only on the top row
            set(ax, 'YTick', rate_yticks);   % just 0 and 0.25
        end
    end
    linkaxes(cax(r, :), 'y');   % shared probability axis within each row
end

% --- Open more space between the two rows + lift condition titles ---------
% Shrink each axis's height and slide the stack down: the top-fixed shrink
% widens the inter-row gap while the shift frees room above the top row. Then
% raise the condition titles into column-header position so they clearly label
% the whole column -- matching Fig_sensitivity_analysis.m.
for r = 1:nRows
    for c = 1:nCols
        ax = cax(r, c);
        p  = get(ax, 'Position');            % [left bottom width height]
        new_h = p(4) * row_shrink;
        set(ax, 'Position', [p(1), p(2) + (p(4) - new_h) - top_headroom, p(3), new_h]);
    end
end
for c = 1:nCols
    t = get(cax(1, c), 'Title');
    if ~isempty(get(t, 'String'))
        set(t, 'Units', 'normalized', 'Position', [0.5, title_y, 0], ...
            'VerticalAlignment', 'bottom', 'FontSize', title_fs);
    end
end
% Match the colorbar's vertical extent to the (now shrunk + shifted) top row, and
% nudge it left within its cell. The bar is centered in its Position box by
% pbaspect, so shifting the box left moves the bar left (increase cb_x_shift for
% more; normalized figure units).
cb_x_shift = 0.045;
if isgraphics(cbax)
    p = get(cbax, 'Position');
    new_h = p(4) * row_shrink;
    set(cbax, 'Position', [p(1) - cb_x_shift, p(2) + (p(4) - new_h) - top_headroom, p(3), new_h]);
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
x_shift = 0.012;   % nudge dividers slightly left (normalized figure units)
for c = 1:ncol-1
    x_div = (col_right(c) + col_left(c+1)) / 2 - x_shift;
    annotation(combined, 'line', [x_div x_div], [y_bot y_top], ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 2);
end

% Panel letters (a)..(h) on the 8 data panels only (reading order: row 1 left-
% to-right, then row 2). Pass the axes explicitly so the colorbar in column 5 is
% not lettered.
letter_axes = cell(1, numel(cax));
k = 0;
for r = 1:nRows
    for c = 1:nCols
        k = k + 1;
        letter_axes{k} = cax(r, c);
    end
end
panel_letters = arrayfun(@(ch) sprintf('(%c)', ch), ...
    char('a' + (0:numel(cax)-1)), 'UniformOutput', false);
AddLetters2Plots(letter_axes, panel_letters, ...
    'FontSize', letter_fs, 'FontWeight', 'normal', 'HShift', -0.03, 'VShift', -0.06);

% 6) Colorbar (now embedded in the upper-right cell): relabel the f gradient bar
%    with E:I ratios (E:I = f:(1-f)).
if isgraphics(cbax)
    ei_f   = [0.25, 1/3, 0.4, 0.5, 0.6, 2/3, 0.75];
    ei_lab = {'1:3', '1:2', '2:3', '1:1', '3:2', '2:1', '3:1'};
    ylim_cb = get(cbax, 'YLim');
    keep = ei_f >= ylim_cb(1) - 1e-6 & ei_f <= ylim_cb(2) + 1e-6;
    set(cbax, 'YTick', ei_f(keep), 'YTickLabel', ei_lab(keep), 'FontSize', tick_fs);
    ylabel(cbax, 'E:I ratio', 'FontSize', label_fs);
end

% 7) Save the single combined figure (colorbar now lives inside it).
save_fig_stable(combined, out_dir, 'Fig_EI_ParamSpace');

% 8) Log the git state alongside the figures.
capture_git_provenance(out_dir, project_root);

%% -------------------- Human-readable description --------------------
desc_path = fullfile(out_dir, 'README_fig_EI_param_space.txt');
fid = fopen(desc_path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Stability_Manuscript figure: Parameter-space distributions (E:I colored)\n');
fprintf(fid, '=======================================================================\n\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now')));
fprintf(fid, 'By script: %s.m\n\n', mfilename);

fprintf(fid, 'HOW IT WAS MADE\n');
fprintf(fid, '  Presentation replot -- no simulation is re-run. load_and_make_unit_histograms\n');
fprintf(fid, '  builds per-metric (LLE + mean_rate) 1x4 histograms where each bar is a stack\n');
fprintf(fid, '  of per-network patches colored by f (blue_gray_red_colormap; blue = low f /\n');
fprintf(fid, '  inhibition-dominated, red = high f / excitation-dominated), LLERange [-1.5,1.5],\n');
fprintf(fid, '  probability-normalized. Those axes are copied into a single 2x5 figure:\n');
fprintf(fid, '    row 1 = LLE ("Growth Rate", green dashed zero line)\n');
fprintf(fid, '    row 2 = mean firing rate\n');
fprintf(fid, '    columns 1-4 = No Adaptation, SFA, STD, SFA+STD\n');
fprintf(fid, '    column 5   = f-gradient colorbar (upper-right cell); lower-right cell empty\n');
fprintf(fid, '  Cleanups: condition titles raised into column-header position above the top\n');
fprintf(fid, '  row (not bold, enlarged to match the sensitivity figure); extra row spacing;\n');
fprintf(fid, '  y-ticks reduced (row 1: 0, 0.5; row 2: 0, 0.25); vertical gray column\n');
fprintf(fid, '  dividers; panel letters (a)..(h) on the 8 data panels only; fonts matched to\n');
fprintf(fid, '  the MC/sensitivity figures; y-axes linked within each row. The embedded\n');
fprintf(fid, '  colorbar encodes f as an E:I ratio\n');
fprintf(fid, '  (ticks 1:3, 1:2, 2:3, 1:1, 3:2, 2:1, 3:1). See git_provenance.txt.\n\n');

fprintf(fid, 'SOURCE RUN\n');
fprintf(fid, '  %s\n', data_root);
fprintf(fid, '  param_space subfolder used:\n');
fprintf(fid, '    %s\n\n', ps_dirs(1).name);

fprintf(fid, 'FIGURES PRODUCED (in this folder)\n');
fprintf(fid, '  Fig_EI_ParamSpace.png / .svg / .fig   (2x5: f-colored distributions + colorbar)\n');

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
