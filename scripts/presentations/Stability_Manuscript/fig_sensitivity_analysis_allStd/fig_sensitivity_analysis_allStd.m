function out = fig_sensitivity_analysis_allStd(cfg)
% FIG_SENSITIVITY_ANALYSIS_ALLSTD Stacked 1-D sensitivity sheets, LLE and rate.
%
%   out = FIG_SENSITIVITY_ANALYSIS_ALLSTD()
%   out = FIG_SENSITIVITY_ANALYSIS_ALLSTD('run_dir', d)
%
% Rows are swept parameters, columns are adaptation conditions, each panel an
% imagesc of the full rep distribution with a median line. No simulation is
% re-run: the saved 1-D sensitivity PSA objects are reloaded, replot_sensitivity
% regenerates the per-parameter panels, and assemble_sensitivity_figure stacks
% them.
%
% TWO SHEETS PER METRIC. The Pairs sweeps cover SEVEN parameters, and seven rows
% on one sheet is about 2100 px tall. They are split by what the parameters MEAN
% rather than alphabetically: 'core' is the network-scale axes (f_E,
% level_of_chaos, n), 'mu' the four connectivity blocks, which are only really
% comparable side by side.
%
% ROWS ARE IDENTIFIED BY INDEX, not guessed from the data. The original inferred
% which parameter a row showed from max(xlim) -- <=1 is f, <=10 is
% level_of_chaos, else n -- which silently mislabels the Pairs sweeps
% (mu_EE_relative, max 11, reads as "Network Size"). assemble_sensitivity_figure
% takes an ordered parameter list, so row -> parameter is already known exactly.
%
% NO 'close all force'. Correct standalone -- replot_sensitivity saves ALL open
% figures -- but in a batch it destroys the previous figure's output before the
% caller can collect it. Explicitly named handles are saved instead.
%
% See also: resolve_run_dir, replot_sensitivity, assemble_sensitivity_figure,
%           preset_default_values, fig_sensitivity_medians

arguments
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
end

setup_paths();
out_dir      = default_out_dir(cfg.out_dir, mfilename('fullpath'));
project_root = fileparts(which('setup_paths'));
st           = manuscript_style();

run_dir = resolve_run_dir('run_dir', cfg.run_dir, 'preset_name', cfg.preset_name);
made_tags = {};

% stray figure lingering in the session (e.g. a previous combined figure) would
% pollute the per-param save and break assemble_sensitivity_figure.
% (no 'close all force' -- destroyed sibling figures in a batch; see header)

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
default_value = preset_default_values(run_dir, ...
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

tick_fs   = st.tick_fs;
label_fs  = st.label_fs;
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
% Median line: BRIGHT blue at 50% transparency. A 4-element Color is [r g b a],
% so the alpha travels with the colour rather than needing a separate handle
% property.
%
% This reverses an earlier choice, and the reasoning it reversed is still on the
% record: the line was made opaque dark blue because a transparent line reads as
% a wash whose apparent colour shifts with the histogram density underneath it,
% rather than as one curve. That effect is real and is the cost here. What is
% bought back is that the median no longer hides the cells it crosses, which
% matters most exactly where the distribution is densest and the line is
% therefore least separable from it. Chosen deliberately (TR, 2026-08-26); the
% dark-grey colormap ceiling above is what keeps a 50%-alpha line legible at all.
median_color = [0.10 0.50 1.00 0.50];   % bright blue, 50% alpha
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
replot_dir = replot_sensitivity(run_dir, lle_range, n_bins, n_bins);

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

        % Median line: bright blue at 50% alpha + 25% thinner. median_color is a
        % 4-element [r g b a], so this KEEPS an alpha channel where a 3-element
        % Color would drop it back to opaque. (imagesc is Type 'image', the zero
        % line is 'constantline', so 'line' is the median.)
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
    letters = panel_letters(numel(ax_all));   % (a).. and past (z) correctly
    AddLetters2Plots(cf, letters, ...
        'FontSize', letter_fs, 'FontWeight', 'normal', 'HShift', -0.03, 'VShift', -0.04);

    % --- Save ONLY the cleaned combined figure, with a STABLE name ---------
    fig_tag = sprintf('%s_%s', spec.fig_tag, sheet.tag);
    save_figure_stable(out_dir, fig_tag, cf);
    made_tags(end+1) = {fig_tag};
end
end

%% --- Clean up and record ----------------------------------------------------
% Prep figures (per-param LLE + mean_rate, and the intermediate combined ones)
% exist only to build the finals -- remove the whole replot folder so no extra
% figs are left in the data dir. Must come after BOTH metrics.
if isfolder(replot_dir)
    rmdir(replot_dir, 's');
end

out = struct('figs', gobjects(0), 'files', {{}}, 'source', run_dir);
if cfg.save
    for k = 1:numel(made_tags)
        out.files = [out.files, existing_outputs(out_dir, made_tags{k})];
    end
    capture_git_provenance(out_dir, project_root);

    axis_note = [ ...
        'x-axes are relabelled per parameter: f_E to E:I ratio; level_of_chaos ' ...
        'to Synaptic Gain; n to Network Size; mu_XY_relative to mu_{XY} (RMT ' ...
        'block means, indexed (post, pre)). The Synaptic Gain and the four mu ' ...
        'axes are shown as PERCENT DEPARTURE from the preset default, ' ...
        '(value/default - 1)*100, because absolute mu_tilde_relative values ' ...
        'mean little on their own; this puts the preset own network at 0 percent ' ...
        'and makes the four mu panels directly comparable. SIGN: mu_EI and ' ...
        'mu_II have NEGATIVE defaults, so +100 percent means twice as ' ...
        'inhibitory, which in raw coordinates is further LEFT -- those two ' ...
        'panels have their x-direction reversed so that on all four, rightward ' ...
        'means "stronger synapse of this type". Only the ruler changes; the ' ...
        'underlying image is untouched.'];

    default_note = [ ...
        'Every row carries a short reddish-grey tick rising from the x-axis at ' ...
        'the preset default for that parameter -- the network the sweep departs ' ...
        'from. On the percent axes it sits at 0 percent; on the f_E and Network ' ...
        'Size rows there is no 0 percent tick, which is where it earns its ' ...
        'keep. A default outside the swept range is NOT drawn rather than ' ...
        'clamped to the edge. The values come from preset_default_values, which ' ...
        'reads the run own run_manifest.mat, constructs a model from that ' ...
        'preset and reads each value off the object, so the CLASS accessors ' ...
        'resolve the aliases. resolved_defaults cannot serve here: ' ...
        'ParamSpaceAnalysis2 excludes grid axes from it, and each of these ' ...
        'parameters is the axis of one of the sweeps.'];

    metric_note = [ ...
        'Both metrics use the same bin count so they have the same vertical ' ...
        'resolution. The LLE sheets keep the green dashed zero line -- the sign ' ...
        'of lambda_1 marks the edge of chaos -- and the solid bands at top and ' ...
        'bottom are the overflow bins. On the aug_14 preset the LLEs spanned ' ...
        'roughly p1 = -10.0 to p99 = +3.7, so those bands carried a large share ' ...
        'of the distribution; re-check that against the run actually plotted. ' ...
        'The rate sheets use range [0, 1], where nothing can fall outside, so ' ...
        'their overflow bins are always empty and the zero line is dropped as ' ...
        'meaningless for a rate.'];

    write_figure_readme(out_dir, struct( ...
        'tag',    'fig_sensitivity_analysis_allStd', ...
        'title',  'Stability_Manuscript figures: LLE and mean firing rate sensitivity', ...
        'script', 'fig_sensitivity_analysis_allStd.m', ...
        'what',   ['One row per swept parameter, one column per adaptation ' ...
                   'condition. Each panel is an imagesc of the full ' ...
                   'across-reps distribution at each level, with the median ' ...
                   'overlaid as a bright blue line at 50% transparency, so it ' ...
                   'marks the central tendency without hiding the cells it ' ...
                   'crosses. Two sheets per metric: core (f_E, ' ...
                   'level_of_chaos, n) and mu (the four connectivity blocks).'], ...
        'how',    ['Presentation replot -- no simulation is re-run. ' ...
                   'replot_sensitivity reloads the saved PSA objects and ' ...
                   'regenerates both metrics; assemble_sensitivity_figure ' ...
                   'stacks them per sheet; the result is restyled and saved ' ...
                   'here, and the prep folder is deleted.'], ...
        'source', struct('run_dir', run_dir, 'preset', cfg.preset_name), ...
        ... % clim_frac is per-metric (it lives in the `spec` struct array), so
        ... % record the LLE sheet's value, which is the one worth quoting.
        'settings', struct('lle_range', lle_range, 'n_bins', n_bins, ...
                           'clim_frac_LLE', spec(1).clim_frac), ...
        'figures', {out.files}, ...
        'sections', struct( ...
            'heading', {'axis conventions', 'default marker', 'per-metric differences'}, ...
            'body',    {axis_note, default_note, metric_note})));
end
end
