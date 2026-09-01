function out = fig_EI_param_space(cfg)
% FIG_EI_PARAM_SPACE Param-space distributions, each bar coloured by E:I ratio.
%
%   out = FIG_EI_PARAM_SPACE()
%   out = FIG_EI_PARAM_SPACE('run_dir', d)
%
% A 2 x 5 sheet. Columns 1-4 are the four adaptation conditions (row 1 = LLE,
% row 2 = mean rate); each bar is a STACK of per-network patches coloured by the
% fraction excitatory (blue = inhibition-dominated, red = excitation-dominated).
% Column 5 embeds the colorbar, encoding the fraction as an E:I ratio.
%
% The sibling of fig_param_space_allStd, which shows the same distributions in
% plain grey. Both are paper figures: this one answers "does the E:I balance
% explain where in the distribution a network lands", which the grey sheet
% cannot.
%
% TWO THINGS CHANGED IN THE PORT.
%
% 1. IT WAS BUILT FROM A DIFFERENT RUN. The hardcoded data_root pointed at
%    run_all_jul_06_26_22_00 -- an old SRNNModel2 run -- while every sibling
%    figure pointed at an SRNNCellTypePairs run. The manuscript was therefore
%    showing two param-space figures computed from different models. It now
%    resolves the same run as everything else.
%
% 2. THE COLOUR AXIS IS NAMED, not defaulted. load_and_make_unit_histograms
%    colours by its ColorBy option, default 'f'. On SRNNCellTypePairs `f` is a
%    1 x C ROW, not a scalar, so the default breaks the per-network colour
%    assignment outright. The Pairs name for the same quantity is the scalar
%    alias f_E, which is exactly f(1) (SRNNCellTypePairs.m get.f_E), so the
%    colormap and the E:I tick labels carry over unchanged -- only the name
%    read off each result differs. It is resolved through
%    psa.effective_param, NOT res.config: effective_param('f') on a Pairs run
%    returns the class default [0.5 0.5] rather than the swept value, which
%    would have coloured every network identically.
%
% See also: fig_param_space_allStd, load_and_make_unit_histograms, resolve_run_dir

arguments
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    cfg.out_dir     (1,:) char    = ''
    cfg.color_by    (1,:) char    = ''      % '' -> per model class (f_E / f)
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
end

setup_paths();
out_dir      = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st           = manuscript_style();

run_dir = resolve_run_dir('run_dir', cfg.run_dir, 'preset_name', cfg.preset_name);

% The colour axis, per model class. get.f_E asserts exactly two cell types, so
% a three-type run must name its own axis rather than fall through to f_E.
color_by = cfg.color_by;
if isempty(color_by)
    [~, model_class] = srnn_param_preset(cfg.preset_name);
    if strcmp(model_class, 'SRNNCellTypePairs')
        color_by = 'f_E';
    else
        color_by = 'f';
    end
end

ps_dirs = dir(fullfile(run_dir, 'param_space_*'));
ps_dirs = ps_dirs([ps_dirs.isdir]);
assert(~isempty(ps_dirs), 'No param_space_* subdir found in %s', run_dir);
ps_dir  = fullfile(ps_dirs(1).folder, ps_dirs(1).name);

% Start from a clean slate.
% (no 'close all force' -- destroyed sibling figures in a batch; see header)

% 1) Build the f-colored per-metric histogram figures (LLE + mean_rate), matching
%    the gray figure's LLE range [-1.5,1.5] and probability normalization. This
%    also creates a separate 'f Value Colorbar' figure.
% CLim is FIXED at the swept range [0.2, 0.8] = 1:4 to 4:1, not derived from the
% data. Derived limits came from the min/max of the f_E values actually sampled,
% which was fine when the grid enumerated every level but is luck once the grid
% is randomly subsampled: a run that never draws f_E = 0.2 loses the 1:4 tick,
% the bar's span changes, and two runs' colours stop meaning the same thing. The
% axis is a property of the SWEEP, so it belongs in the script, not in the data.
EI_CLIM = [0.2, 0.8];
[~, ~] = load_and_make_unit_histograms(ps_dir, ...
    'Metrics', {'lle', 'r'}, 'NormalizeMode', 'probability', 'LLERange', [-1.5, 1.5], ...
    'ColorBy', color_by, 'CLim', EI_CLIM);

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
tick_fs  = st.tick_fs;
label_fs = st.label_fs;
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
letters = panel_letters(numel(cax));   % (a).. and past (z) correctly
AddLetters2Plots(letter_axes, letters, ...
    'FontSize', letter_fs, 'FontWeight', 'normal', 'HShift', -0.03, 'VShift', -0.06);

% 6) Colorbar (now embedded in the upper-right cell): relabel the f gradient bar
%    with E:I ratios (E:I = f:(1-f)).
if isgraphics(cbax)
    % Spans the FULL swept range, 0.2 to 0.8, not just 0.25-0.75. The keep
    % filter is retained as a guard, but with EI_CLIM fixed above it now keeps
    % every entry every time -- which is the point: the end ticks 1:4 and 4:1
    % are always drawn, rather than appearing only when the random subsample
    % happens to include a network at the extremes.
    ei_f   = [0.2, 0.25, 1/3, 0.4, 0.5, 0.6, 2/3, 0.75, 0.8];
    ei_lab = {'1:4', '1:3', '1:2', '2:3', '1:1', '3:2', '2:1', '3:1', '4:1'};
    ylim_cb = get(cbax, 'YLim');
    keep = ei_f >= ylim_cb(1) - 1e-6 & ei_f <= ylim_cb(2) + 1e-6;
    assert(all(keep), ['E:I colorbar lost %d tick(s): CLim is fixed at ' ...
        '[%g %g] so every label should fit.'], sum(~keep), EI_CLIM(1), EI_CLIM(2));
    set(cbax, 'YTick', ei_f(keep), 'YTickLabel', ei_lab(keep), 'FontSize', tick_fs);
    ylabel(cbax, 'E:I ratio', 'FontSize', label_fs);
end

if ~cfg.visible; set(combined, 'Visible', 'off'); end

%% --- Save -------------------------------------------------------------------
fig_tag = 'Fig_EI_ParamSpace';
out = struct('figs', combined, 'files', {{}}, 'source', ps_dir);
if cfg.save
    save_figure_stable(out_dir, fig_tag, combined);
    out.files = existing_outputs(out_dir, fig_tag);

end
end



