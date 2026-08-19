close all
clc

% Fig_sensitivity_analysis_allStd.m
% Stability_Manuscript presentation figures: combined LLE + mean_rate sensitivity,
% for the SRNNCellTypePairs "allStd" runs. Regenerates the stacked 1D-sensitivity
% figures from a saved run_all_<dt> run and writes the final figures into this
% presentation folder. No simulation is re-run -- edit data_root below to point at
% a different run_all output.
%
% Sibling of ../fig_sensitivity_analysis/Fig_sensitivity_analysis.m, which is left
% untouched. Three things differ, all forced by the move to SRNNCellTypePairs:
%
%   1. SEVEN swept parameters, not three. The Pairs sweeps add the four mu blocks
%      (mu_EE/EI/IE/II_relative) to f_E / level_of_chaos / n. Seven rows on one
%      sheet is ~2100 px tall, so each metric is split into two sheets: the three
%      "core" network parameters, and the four mu blocks (which are only really
%      comparable side by side anyway).
%   2. Rows are identified BY INDEX, not by guessing from the data. The original
%      inferred which parameter a row showed from max(xlim) -- <=1 is f, <=10 is
%      level_of_chaos, else n. That silently mislabels the Pairs sweeps:
%      mu_EE_relative (max 11) reads as "Network Size" and mu_EI_relative
%      (max -2.75) as "E:I ratio". Since assemble_sensitivity_figure takes an
%      ordered params list, row index -> parameter is already known exactly, so
%      the heuristic is replaced by a lookup. See row_style below.
%   3. f -> f_E, spanning [0.2, 0.8] rather than [0.25, 0.75], so the E:I ratio
%      ticks are different. Note effective_param(res,'f') on a Pairs run returns
%      the class default [0.5 0.5] rather than the swept value -- read f_E.
%
% See also: Fig_sfa_EOC_allStd, Fig_param_space_allStd, assemble_sensitivity_figure

this_dir     = fileparts(mfilename('fullpath'));
% Depth-independent project root (CLAUDE.md convention) -- setup_paths.m lives at
% the repo root, so this tolerates the script moving between subdirectories.
setup_paths();
project_root = fileparts(which('setup_paths'));

% Source run (a run_all_<dt> folder with 1D_sensitivity_* subdirs).
% Swap this one line to regenerate against the medium run.
data_root = fullfile(project_root, 'data', 'param_space', 'run_all_aug_14_26_17_25');
out_dir   = this_dir;   % write the final figures next to this script

% Start from a clean slate: replot_sensitivity saves ALL open figures, so any
% stray figure lingering in the session (e.g. a previous combined figure) would
% pollute the per-param save and break assemble_sensitivity_figure.
close all force;

% --- Which parameters go on which sheet ------------------------------------
% Order here IS the row order, and is what row_style is indexed by. Names must
% match the figure Name suffix that plot_sensitivity writes
% ("LLE Sensitivity - mu_EE_relative"), i.e. the raw property names.
sheets = struct( ...
    'tag',    {'core',                              'mu'}, ...
    'params', {{'f_E', 'level_of_chaos', 'n'}, ...
               {'mu_EE_relative', 'mu_EI_relative', 'mu_IE_relative', 'mu_II_relative'}});

% --- Per-parameter axis styling --------------------------------------------
% Replaces the original's max(xlim) heuristic. Keys are parameter names; each
% entry gives the x-axis label and, where the automatic ticks are unhelpful, an
% explicit tick vector (+ labels). Empty ticks means "leave MATLAB's alone".
% f_E: E:I = f_E:(1-f_E), so 0.2 -> 1:4. Ratios are clearer than the fraction.
% mu_*: RMT block means, indexed (post, pre) -- mu_EE is E<-E. Rendered in tex.
% The gain and mu axes are re-expressed as PERCENT OF THE PRESET DEFAULT below,
% so their raw ticks are left empty here.
row_style = containers.Map( ...
    {'f_E', 'level_of_chaos', 'n', ...
     'mu_EE_relative', 'mu_EI_relative', 'mu_IE_relative', 'mu_II_relative'}, ...
    {struct('xlabel', 'E:I ratio',     'xticks', [0.2 0.4 0.6 0.8], ...
                                       'xticklabels', {{'1:4','2:3','3:2','4:1'}}), ...
     struct('xlabel', 'Synaptic Gain', 'xticks', [], 'xticklabels', {{}}), ...
     struct('xlabel', 'Network Size',  'xticks', [100 500 1000], 'xticklabels', {{}}), ...
     struct('xlabel', '\mu_{EE}',      'xticks', [], 'xticklabels', {{}}), ...
     struct('xlabel', '\mu_{EI}',      'xticks', [], 'xticklabels', {{}}), ...
     struct('xlabel', '\mu_{IE}',      'xticks', [], 'xticklabels', {{}}), ...
     struct('xlabel', '\mu_{II}',      'xticks', [], 'xticklabels', {{}})});

% --- Axes shown as percent of the preset default ---------------------------
% Absolute mu_tilde values mean little on their own; what the sweep actually
% varies is the departure from the network the preset defines. Expressing these
% axes as (value/default - 1)*100 makes the four mu panels directly comparable
% and puts the preset's own network at 0%.
%
% The defaults come from the PRESET NAMED IN THE RUN'S OWN MANIFEST, not from a
% constant here, so this stays correct when data_root is re-pointed at a run
% built from a different preset. mu_tilde_relative is a C x C block indexed
% (post, pre); a 1 x C row is a presynaptic row broadcast down the columns.
%
% NOTE ON SIGN. mu_EI and mu_II have NEGATIVE defaults, so "+100%" means twice
% as inhibitory -- which in raw data coordinates is further LEFT. Those panels
% therefore get XDir reversed, so that on all four mu panels "further right"
% means "stronger synapse of this type" and the percent axis always ascends
% left-to-right. Without that the two inhibitory panels would read backwards.
% The defaults are resolved for EVERY swept parameter, not just the percent
% ones, because they also place the default marker (see default_mark_* below):
% on the f_E and Network Size rows there is no 0% tick to carry that
% information, which is exactly where the marker earns its keep.
default_value = preset_default_values(data_root, ...
    {'f_E', 'level_of_chaos', 'n', 'mu_EE_relative', 'mu_EI_relative', ...
     'mu_IE_relative', 'mu_II_relative'});
pct_params = {'level_of_chaos', 'mu_EE_relative', 'mu_EI_relative', ...
              'mu_IE_relative', 'mu_II_relative'};

% --- Histogram range + styling constants -----------------------------------
% Kept identical to the original figure by decision, so the two read the same.
% Worth knowing when reading the output: this preset's LLEs span roughly
% p1 = -10.0 to p99 = +3.7, well outside [-1.75, 1.75], so the solid bands at the
% top and bottom of each panel (the -inf/+inf overflow bins) carry a large share
% of the distribution. Widen lle_range here if that hides too much.
lle_range = [-1.75, 1.75];
n_bins    = 24;    % counts linspace EDGES, so 25 plotted rows (23 interior + 2 overflow)

tick_fs   = 14;    % tick numbers -- matches MC figure DefaultAxesFontSize
label_fs  = 15.4;  % axis labels -- matches MC figure (14 * 1.1 label multiplier)
title_fs  = 20;    % condition titles (no adaptation, sfa only, ...) -- enlarged
letter_fs = 18;    % panel letters -- matches MC figure
row_shrink = 0.85; % shrink each row's height to open gaps between rows
top_headroom = 0.06; % shift the row stack down (normalized) to clear room above the top row for column headers
title_y   = 1.22;  % condition-title height above the top-row axes (normalized), reads as a column header
% Colormap ramps white (0 counts) -> DARK GRAY (max), stopping well short of
% black. A ramp into black makes the densest cells read as a solid slab and
% swallows the overflow bands, which then dominate the panel; ending at 0.35 grey
% keeps the density gradient legible at the top end and leaves the median line
% clearly on top of it. Built inline rather than as a named colormap file -- it
% is one expression, and the interesting parameter is the end point, which is
% easier to see and tune here than behind a function call.
cmap_darkest = 0.35;   % grey level of the densest cell (0 = black, 1 = white)
dark_cmap = repmat(linspace(1, cmap_darkest, 256)', 1, 3);
% Median line: OPAQUE dark blue, not a transparent pure blue. Transparency was
% there to keep the line readable over near-black cells; with the colormap
% stopping at dark grey it is no longer needed, and an opaque line reads as one
% curve rather than as a wash whose apparent colour changes with the density
% underneath it.
median_color = [0 0 0.55];
median_lw    = 3;      % blue median line width, 25% thinner (plot_sensitivity uses 4)
zeroline_lw  = 2;      % green dashed zero line width (plot_sensitivity uses 4)

% Short tick rising off the x-axis at the preset's default value for that row --
% the "you are here" mark for the network the sweep departs from. Deliberately
% low-contrast: it is a reference, not data, and must not compete with the
% median line or the histogram.
default_mark_color = [0.80 0.60 0.60];  % reddish light gray
default_mark_lw    = 2;                 % points
default_mark_frac  = 0.05;              % height as a fraction of the y-range

% Regenerate the per-param figures. One call produces BOTH the LLE and the
% mean_rate figures for every swept parameter (replot_sensitivity plots both).
replot_dir = replot_sensitivity(data_root, lle_range, n_bins, n_bins);

% --- One combined figure per metric per sheet ------------------------------
% assemble_sensitivity_figure matches the per-param figs by their Name prefix
% ("<metric> Sensitivity - "), so the same cleanup serves both metrics. Only the
% ylabel, the zero line and the output name differ:
%   LLE       -> lambda_1, green dashed zero line kept (sign of lambda_1 matters)
%   mean_rate -> "Mean Firing Rate"; zero line dropped (range is [0,1], so
%                yline(0) just sits on the bottom axis and carries no meaning)
% clim_frac darkens imagesc: CLim is capped at total_reps*clim_frac (shared
% across panels of a metric so they stay comparable).
% yticks: [] leaves the automatic ticks alone. mean_rate is pinned to just the
% 0 and 1 endpoints -- the intermediate 0.5 adds clutter without information.
metric_specs = struct( ...
    'name',      {'LLE',                 'mean_rate'}, ...
    'ylabel',    {'\lambda_1',           'Mean Firing Rate'}, ...
    'fig_tag',   {'Fig_Sensitivity_LLE', 'Fig_sensitivity_mean_rate'}, ...
    'zero_line', {true,                  false}, ...
    'clim_frac', {0.8,                   0.8}, ...
    'yticks',    {[],                    [0, 1]});

made_tags = {};
for mi = 1:numel(metric_specs)
    spec = metric_specs(mi);

for si = 1:numel(sheets)
    sheet = sheets(si);

    % Assemble this sheet's rows into one stacked figure, then re-open it so it
    % can be cleaned up and saved into the presentation folder (png + svg + fig),
    % matching the memory-capacity figure's formats/naming. Passing params fixes
    % the row order, which is what makes the row_style lookup below exact.
    combined = assemble_sensitivity_figure(replot_dir, spec.name, sheet.params, sheet.tag);
    if isempty(combined)
        warning('Fig_sensitivity_allStd:NoSheet', ...
            'No %s figures for sheet ''%s''; skipping.', spec.name, sheet.tag);
        continue;
    end
    combined_fig_path = fullfile(replot_dir, 'figures', ...
        sprintf('sensitivity_%s_combined_%s.fig', spec.name, sheet.tag));
    cf = openfig(combined_fig_path, 'visible');

    % --- Reshape the whole figure -----------------------------------------
    % Make it ~15% narrower and ~5% taller so the imagesc panels are less wide
    % (the added row spacing / raised headers had stretched their aspect ratio).
    % Axes use normalized positions, so they follow the new figure proportions.
    fig_pos = get(cf, 'Position');            % [left bottom width height]
    set(cf, 'Position', [fig_pos(1), fig_pos(2), fig_pos(3) * 0.85, fig_pos(4) * 1.05]);

    % --- Map each axis to its ROW ------------------------------------------
    % assemble_sensitivity_figure lays the rows out with subplot(nRows,nCols,..)
    % in the order given by sheet.params, so grouping axes by their normalized
    % bottom edge and sorting top-down recovers row index exactly. This replaces
    % the original's max(xlim) guessing, which mislabels the mu sweeps.
    ax_all = findobj(cf, 'Type', 'axes');
    ax_pos = cell2mat(get(ax_all, 'Position'));       % [left bottom width height]
    [~, ~, row_of] = uniquetol(ax_pos(:,2), 0.01);    % group by bottom edge
    row_bottoms = accumarray(row_of, ax_pos(:,2), [max(row_of) 1], @mean);
    [~, top_down] = sort(row_bottoms, 'descend');     % row 1 = topmost
    rank_of_row = zeros(1, numel(top_down));          % group index -> row index
    rank_of_row(top_down) = 1:numel(top_down);
    row_idx = rank_of_row(row_of);

    % --- Clean up axis labels + fonts -------------------------------------
    for a = 1:numel(ax_all)
        ax = ax_all(a);
        set(ax, 'FontSize', tick_fs);         % enlarge x & y tick numbers
        set(get(ax, 'Title'), 'FontWeight', 'normal', 'FontSize', title_fs);  % condition titles, not bold, enlarged
        box(ax, 'off');                       % drop the rectangle; keep x/y axes+ticks
        colormap(ax, dark_cmap);              % white -> dark grey, not black

        % Darken the histogram density: the panels are drawn with CLim = [0,
        % total_reps]; lower the ceiling so typical bin counts use more of the
        % dark range (kept shared across panels so they stay comparable).
        cl = get(ax, 'CLim');
        set(ax, 'CLim', [0, cl(2) * spec.clim_frac]);

        % Median line: opaque dark blue + 25% thinner. Setting a 3-element Color
        % REPLACES plot_sensitivity's 4-element [0 0 1 0.55], which is what drops
        % the alpha. (imagesc is Type 'image', the zero line is 'constantline',
        % so 'line' is the median.)
        ml = findobj(ax, 'Type', 'line');
        set(ml, 'Color', median_color, 'LineWidth', median_lw);
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

        if ~isempty(spec.yticks)              % metric-specific y ticks
            set(ax, 'YTick', spec.yticks);
        end

        % x-axis: labelled from the KNOWN parameter for this row.
        this_param = sheet.params{row_idx(a)};
        rs = row_style(this_param);
        has_default = isKey(default_value, this_param);
        % A zero default has no percent scale (x/0), so such an axis stays in
        % raw units even if it is on the percent list.
        use_pct = has_default && ismember(this_param, pct_params) ...
                  && default_value(this_param) ~= 0;
        if use_pct
            apply_percent_axis(ax, default_value(this_param), rs.xlabel, label_fs);
        else
            xlabel(ax, rs.xlabel, 'Interpreter', 'tex', 'FontSize', label_fs);
            if ~isempty(rs.xticks)
                set(ax, 'XTick', rs.xticks);
                if ~isempty(rs.xticklabels)
                    set(ax, 'XTickLabel', rs.xticklabels);
                end
            end
        end

        % Short tick at the preset's default for this parameter. Drawn LAST so
        % it is not caught by the median-line restyle above, which selects on
        % Type 'line'.
        if has_default
            mark_default_value(ax, default_value(this_param), ...
                default_mark_color, default_mark_lw, default_mark_frac);
        end
        if row_idx(a) > 1
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
    % Draw lines spanning the plot height in the gaps between adjacent columns,
    % to visually separate the conditions. Column boundaries come from the axes'
    % normalized Positions; each line sits at the midpoint of the inter-column
    % gap (so it falls in the gutter, not over any axes).
    pos = cell2mat(get(ax_all, 'Position'));         % re-read: positions moved above
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
    fig_tag = sprintf('%s_%s', spec.fig_tag, sheet.tag);
    save_figure_stable(out_dir, fig_tag, cf);
    made_tags(end+1) = {fig_tag};
end
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
sens_dirs = dir(fullfile(data_root, '1D_sensitivity_*'));
sens_dirs = sens_dirs([sens_dirs.isdir]);

desc_path = fullfile(out_dir, 'README_fig_sensitivity_analysis_allStd.txt');
fid = fopen(desc_path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Stability_Manuscript figures: LLE + mean firing rate Sensitivity (allStd)\n');
fprintf(fid, '========================================================================\n\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now')));
fprintf(fid, 'By script: %s.m\n\n', mfilename);

fprintf(fid, 'HOW THEY WERE MADE\n');
fprintf(fid, '  Presentation replot -- no simulation is re-run. The script reloads the\n');
fprintf(fid, '  saved 1D-sensitivity PSA objects from a run_all_<dt> run and calls\n');
fprintf(fid, '  replot_sensitivity (which plots BOTH metrics) -> assemble_sensitivity_figure\n');
fprintf(fid, '  once per metric per sheet to rebuild the stacked figures (rows = swept\n');
fprintf(fid, '  params, cols = adaptation conditions), then saves them here. See\n');
fprintf(fid, '  git_provenance.txt for the exact commit.\n\n');

fprintf(fid, 'SOURCE RUN\n');
fprintf(fid, '  %s\n', data_root);
fprintf(fid, '  1D_sensitivity subfolders used:\n');
for k = 1:numel(sens_dirs)
    fprintf(fid, '    %s\n', sens_dirs(k).name);
end
fprintf(fid, '\n');

fprintf(fid, 'SHEETS\n');
fprintf(fid, '  This run is SRNNCellTypePairs, which sweeps SEVEN parameters rather than\n');
fprintf(fid, '  the three of the original SRNNModel2 figure. Seven rows on one sheet is\n');
fprintf(fid, '  unreadably tall, so each metric is split in two:\n');
for si = 1:numel(sheets)
    fprintf(fid, '    %-5s : %s\n', sheets(si).tag, strjoin(sheets(si).params, ', '));
end
fprintf(fid, '\n');

fprintf(fid, 'FIGURES PRODUCED (in this folder)\n');
for k = 1:numel(made_tags)
    fprintf(fid, '  %s.png / .svg / .fig\n', made_tags{k});
end

fprintf(fid, '\nSHARED LAYOUT (all sheets)\n');
fprintf(fid, '  One row per swept parameter, one column per adaptation condition.\n');
fprintf(fid, '  x-axes relabelled per parameter: f_E -> "E:I ratio"\n');
fprintf(fid, '  (E:I = f_E:(1-f_E), ticks 1:4, 2:3, 3:2, 4:1); level_of_chaos ->\n');
fprintf(fid, '  "Synaptic Gain"; n -> "Network Size"; mu_XY_relative -> \\mu_{XY}\n');
fprintf(fid, '  (RMT block means, indexed (post, pre)). Rows are identified by their\n');
fprintf(fid, '  index in the sheet''s parameter list, NOT by inspecting the data -- the\n');
fprintf(fid, '  original figure guessed from max(xlim), which mislabels the mu sweeps.\n\n');
fprintf(fid, '  PERCENT AXES. The Synaptic Gain and the four mu axes are shown as\n');
fprintf(fid, '  percent departure from the preset default, (value/default - 1)*100,\n');
fprintf(fid, '  rather than in raw mu_tilde_relative units -- absolute values mean\n');
fprintf(fid, '  little on their own, and this puts the preset''s own network at 0%%\n');
fprintf(fid, '  and makes the four mu panels directly comparable.\n\n');
fprintf(fid, '  SIGN: mu_EI and mu_II have NEGATIVE defaults, so "+100%%" means twice\n');
fprintf(fid, '  as inhibitory. In raw data coordinates that is further LEFT, so those\n');
fprintf(fid, '  two panels have their x-direction REVERSED; on all four mu panels\n');
fprintf(fid, '  rightward therefore means "stronger synapse of this type" and the\n');
fprintf(fid, '  percent axis ascends left-to-right. Only the ruler changes -- the\n');
fprintf(fid, '  underlying image is untouched.\n\n');
fprintf(fid, '  DEFAULT MARKER. Every row carries a short reddish-gray tick rising\n');
fprintf(fid, '  from the x-axis (%g of the y-range, %g pt) at the preset''s default\n', ...
    default_mark_frac, default_mark_lw);
fprintf(fid, '  for that parameter -- the network the sweep departs from. On the\n');
fprintf(fid, '  percent axes it sits at 0%%; on the f_E and Network Size rows there\n');
fprintf(fid, '  is no 0%% tick, which is where it earns its keep. A default lying\n');
fprintf(fid, '  outside the swept range is NOT drawn rather than clamped to the edge.\n\n');
fprintf(fid, '  WHERE THE DEFAULTS COME FROM. Not hardcoded, and not assumed from the\n');
fprintf(fid, '  preset struct''s field layout: the run''s run_manifest.mat names the\n');
fprintf(fid, '  preset and the model class, a model is constructed from that preset,\n');
fprintf(fid, '  and each value is read off the object so the CLASS''s own accessors\n');
fprintf(fid, '  resolve the aliases (f_E -> f(1), mu_EE_relative ->\n');
fprintf(fid, '  mu_tilde_relative(1,1) indexed (post,pre), and a 1 x C row broadcast\n');
fprintf(fid, '  down the columns). The run''s resolved_defaults cannot serve here:\n');
fprintf(fid, '  ParamSpaceAnalysis2 excludes grid axes from it, and each of these\n');
fprintf(fid, '  parameters is the axis of one of the sweeps. For this run:\n');
pk = keys(default_value);
for k = 1:numel(pk)
    fprintf(fid, '    %-16s default %g\n', pk{k}, default_value(pk{k}));
end
fprintf(fid, '\n');
fprintf(fid, '  Larger tick fonts. Condition titles kept only on the top row; vertical\n');
fprintf(fid, '  gray dividers separate the condition columns. imagesc CLim capped at\n');
fprintf(fid, '  total_reps*0.8 (shared within a figure); colormap ramps white -> grey\n');
fprintf(fid, '  %.2f, stopping well short of black -- a ramp into black makes the\n', cmap_darkest);
fprintf(fid, '  densest cells a solid slab and lets the overflow bands dominate the\n');
fprintf(fid, '  panel. Panel letters added up-and-left of each plot (AddLetters2Plots).\n');
fprintf(fid, '  Median line: OPAQUE dark blue [%.2f %.2f %.2f], %g pt (plot_sensitivity\n', ...
    median_color, median_lw);
fprintf(fid, '  draws it as transparent pure blue at 4 pt); the transparency existed\n');
fprintf(fid, '  only to survive near-black cells and is not needed now. Titles not\n');
fprintf(fid, '  bold; axis boxes removed.\n');

fprintf(fid, '\nPER-FIGURE DIFFERENCES\n');
fprintf(fid, '  Both metrics use n_bins = %d (linspace edges), i.e. %d plotted rows:\n', n_bins, n_bins + 1);
fprintf(fid, '  %d interior bins plus the two -inf/+inf overflow bins. Matching the\n', n_bins - 1);
fprintf(fid, '  bin count gives the two metrics the same vertical resolution.\n\n');
fprintf(fid, '  Fig_Sensitivity_LLE_*:       ylabel \\lambda_1; histogram range\n');
fprintf(fid, '                               [%.2f, %.2f]; green dashed zero line kept\n', lle_range(1), lle_range(2));
fprintf(fid, '                               (thinner) -- the sign of lambda_1 marks the\n');
fprintf(fid, '                               edge of chaos. The solid bands at the top and\n');
fprintf(fid, '                               bottom are the overflow bins (reps outside\n');
fprintf(fid, '                               the range). NOTE: this preset''s LLEs span\n');
fprintf(fid, '                               roughly p1 = -10.0 to p99 = +3.7, so those\n');
fprintf(fid, '                               overflow bands carry a large share of the\n');
fprintf(fid, '                               distribution. The range is kept identical to\n');
fprintf(fid, '                               the original figure by choice; widen\n');
fprintf(fid, '                               lle_range in the script to change that.\n');
fprintf(fid, '  Fig_sensitivity_mean_rate_*: ylabel "Mean Firing Rate"; histogram range\n');
fprintf(fid, '                               [0, 1] (plot_sensitivity default; nothing can\n');
fprintf(fid, '                               fall outside it, so the overflow bins are\n');
fprintf(fid, '                               always empty); y ticks at 0 and 1 only; zero\n');
fprintf(fid, '                               line removed -- at y=0 it lands on the bottom\n');
fprintf(fid, '                               axis and carries no meaning for a rate.\n');

clear cleanup;  % flush + close
fprintf('Description written: %s\n', desc_path);


%% ==================== Local helpers ====================
% There are none any more. preset_default_values, apply_percent_axis,
% mark_default_value and save_figure_stable were local subfunctions here until
% fig_sensitivity_medians needed the same four; they now live as standalone
% files in scripts/run_all_analyses/replot/ (on the path via setup_paths) and
% this script calls those. The bodies moved unchanged, EXCEPT that
% apply_percent_axis now signs the zero label ("+0%" rather than "0%"), which is
% why the mu and Synaptic Gain axes here read differently from the figures
% committed before 2026-08-18.
