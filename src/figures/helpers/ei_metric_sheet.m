function combined = ei_metric_sheet(src_fig, cb_fig, spec, opts)
% EI_METRIC_SHEET One styled 1 x (N+1) sheet from a unit-histogram figure.
%
%   fig = EI_METRIC_SHEET(src_fig, cb_fig, spec, opts)
%
% src_fig   the '<field> Unit Histogram' figure load_and_make_unit_histograms
%           made (one axes per condition); cb_fig its colorbar figure (may be
%           empty); spec the sweep_metrics row; opts a struct from
%           manuscript_style-derived constants plus the colorbar labelling:
%             .tick_fs .label_fs .title_fs .axes_lw .letter_fs
%             .row_shrink .top_headroom .title_y .cb_x_shift
%             .xlabel   (char; '' keeps the registry label)
%             .yticks   (probability ticks, [] = automatic)
%             .zero_color
%             .cb_ticks .cb_labels .cb_ylabel .cb_clim
%
% Shared by fig_EI_param_space and fig_EI_weights_param_space, which used to
% carry two copies of this layout code for a fixed two-row (LLE, mean rate)
% sheet. With the measures now coming from the registry each figure makes ONE
% SHEET PER MEASURE, so the layout lives here once and takes the measure as
% data. The source figures are closed on return.
%
% See also: fig_EI_param_space, fig_EI_weights_param_space, sweep_metrics

src_ax = sort_axes_left_to_right(src_fig);
nCols  = numel(src_ax);
nGrid  = nCols + 1;                       % last column holds the colorbar
combined = figure('Color', 'w', 'Position', [100 100 350 * nGrid 300]);
cax = gobjects(1, nCols);
for c = 1:nCols
    ph = subplot(1, nGrid, c, 'Parent', combined);
    target_pos = get(ph, 'Position');
    delete(ph);
    cax(c) = copyobj(src_ax(c), combined);
    set(cax(c), 'Position', target_pos);
end
cbax = gobjects(0);
if ~isempty(cb_fig) && isgraphics(cb_fig)
    ph = subplot(1, nGrid, nGrid, 'Parent', combined);
    cb_target_pos = get(ph, 'Position');
    delete(ph);
    src_cb = findobj(cb_fig, 'Type', 'axes');
    cbax = copyobj(src_cb(1), combined);
    set(cbax, 'Position', cb_target_pos);
end
close(src_fig);
if ~isempty(cb_fig) && isgraphics(cb_fig); close(cb_fig); end

for c = 1:nCols
    ax = cax(c);
    set(ax, 'FontSize', opts.tick_fs, 'LineWidth', opts.axes_lw);
    set(ax.YLabel, 'FontSize', opts.label_fs);
    if isempty(opts.xlabel)
        xlabel(ax, spec.label, 'Interpreter', 'tex', 'FontSize', opts.label_fs);
    else
        xlabel(ax, opts.xlabel, 'FontSize', opts.label_fs);
    end
    set(ax.Title, 'FontWeight', 'normal', 'FontSize', opts.title_fs);
    if spec.zero_line
        set(findobj(ax, 'Type', 'constantline'), 'Color', opts.zero_color);
    end
    if ~isempty(opts.yticks); set(ax, 'YTick', opts.yticks); end
end
linkaxes(cax, 'y');

% Open headroom for the condition titles (they read as column headers).
for c = 1:nCols
    p = get(cax(c), 'Position');
    new_h = p(4) * opts.row_shrink;
    set(cax(c), 'Position', [p(1), p(2) + (p(4) - new_h) - opts.top_headroom, p(3), new_h]);
    t = get(cax(c), 'Title');
    if ~isempty(get(t, 'String'))
        set(t, 'Units', 'normalized', 'Position', [0.5, opts.title_y, 0], ...
            'VerticalAlignment', 'bottom', 'FontSize', opts.title_fs);
    end
end
if isgraphics(cbax)
    p = get(cbax, 'Position');
    new_h = p(4) * opts.row_shrink;
    set(cbax, 'Position', [p(1) - opts.cb_x_shift, p(2) + (p(4) - new_h) - opts.top_headroom, p(3), new_h]);
end

% Column dividers.
pos = cell2mat(get(cax(:), 'Position'));
[~, ~, col_of] = uniquetol(pos(:, 1), 0.01);
ncol      = max(col_of);
col_left  = accumarray(col_of, pos(:, 1),             [ncol 1], @mean);
col_right = accumarray(col_of, pos(:, 1) + pos(:, 3), [ncol 1], @mean);
[col_left, ord] = sort(col_left);
col_right = col_right(ord);
y_bot = min(pos(:, 2));
y_top = max(pos(:, 2) + pos(:, 4));
for c = 1:ncol - 1
    x_div = (col_right(c) + col_left(c + 1)) / 2 - 0.012;
    annotation(combined, 'line', [x_div x_div], [y_bot y_top], ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 2);
end

letter_axes = num2cell(cax);
AddLetters2Plots(letter_axes, panel_letters(nCols), ...
    'FontSize', opts.letter_fs, 'FontWeight', 'normal', 'HShift', -0.03, 'VShift', -0.06);

if isgraphics(cbax)
    ylim_cb = get(cbax, 'YLim');
    keep = opts.cb_ticks >= ylim_cb(1) - 1e-6 & opts.cb_ticks <= ylim_cb(2) + 1e-6;
    assert(all(keep), ['colorbar lost %d tick(s): CLim is fixed at [%g %g] so ' ...
        'every label should fit.'], sum(~keep), opts.cb_clim(1), opts.cb_clim(2));
    set(cbax, 'YTick', opts.cb_ticks(keep), 'YTickLabel', opts.cb_labels(keep), ...
        'FontSize', opts.tick_fs);
    ylabel(cbax, opts.cb_ylabel, 'FontSize', opts.label_fs);
end
end
