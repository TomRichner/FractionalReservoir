function out = fig_param_space_allStd(cfg)
% FIG_PARAM_SPACE_ALLSTD Param-space distributions, one 1 x N sheet per measure.
%
%   out = FIG_PARAM_SPACE_ALLSTD()
%   out = FIG_PARAM_SPACE_ALLSTD('run_dir', d)
%
% One sheet per registry measure with in_sheets set (sweep_metrics: lambda_1,
% mean rate, h_KS, D_KY, ...), one column per adaptation condition, tags
% Fig_ParamSpace_<stem>. No simulation is re-run: the saved param-space PSA
% object is reloaded and its per-condition histogram axes are copied into the
% sheets.
%
% NO 'close all force'. The original opened with it, and correctly so on its own
% terms -- replot_param_space_analysis saves ALL open figures, so a stray one
% would pollute the save. But run from a master script that renders figures in
% sequence, it DESTROYS the figures the previous entry just created, before the
% caller can collect or verify them. The requirement is met instead by saving an
% explicitly named handle, and by deleting the prep folder afterwards.
%
% WHICH CONDITIONS ARE DRAWN: cfg.conditions (FD-006). A cell of condition
% names, one column each, LEFT TO RIGHT IN THE ORDER GIVEN; empty (the
% default) draws every condition the param-space sample holds, in its order,
% exactly as before the option existed. Resolved by select_conditions inside
% replot_param_space_analysis against the run's own PSA object, so a name the
% run lacks is an error, not a missing column. Column titles come from
% srnn_condition_titles by NAME (psa.plot). The sheet is 1 x N, 350 px per
% column.
%
% THE lambda_1 AXIS: cfg.lle_scale. 'window' (the default) is the sheet as
% psa.plot draws it: 25 bins over the registry's dist_range (+-1.5 1/s) and
% EVERYTHING BEYOND GATHERED INTO THE TWO END BINS. That hides the
% distribution of a condition whose values lie outside the window: without
% adaptation three quarters of the sample sat in the left end bin (the
% values run from -10 to +8 1/s). 'symlog' redraws the lambda_1 sheet from
% the sample's own values on a symmetric-log axis, t = sign(x) log10(1 +
% |x| / 0.1): linear within about +-0.1 1/s, logarithmic beyond, wide enough
% for every value of every column, 40 bins of equal width in t (an edge at
% zero), probability per bin. Nothing is gathered. The other sheets are not
% touched. The values of each column are left in its axes' UserData
% (.values, .condition, .scale) for a caller's table. S5 (fig_supp5_grouped)
% asks for 'symlog' (2026-11-22).
%
% See also: resolve_run_dir, replot_param_space_analysis, fig_EI_param_space,
%           select_conditions

arguments
    cfg.verbose     (1,:) char    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (see verbose_level)
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.conditions                = {}          % names to draw, in order; empty = all the run holds (select_conditions)
    cfg.lle_scale   (1,:) char {mustBeMember(cfg.lle_scale, {'window', 'symlog'})} = 'window'
end

setup_paths();
out_dir      = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st           = manuscript_style();

run_dir = resolve_run_dir('run_dir', cfg.run_dir, 'preset_name', cfg.preset_name);

% --- Presentation constants ------------------------------------------------
tick_fs  = st.tick_fs;
label_fs = st.label_fs;
title_fs = st.title_fs;
x_shift = 0.007;   % nudge column dividers slightly left (normalized figure units)

% 1) Regenerate one distribution figure per registry measure with in_sheets
%    set (sweep_metrics) into a replot_param_space_<dt>/figures/ folder.
%    (no 'close all force' -- see the header note; it destroyed sibling figures)
replot_dir = replot_param_space_analysis(run_dir, 'conditions', cfg.conditions, 'verbose', cfg.verbose);
specs = sweep_metrics();
specs = specs([specs.in_sheets]);

% 2) Re-open the saved figures (invisible) and map them by Name.
by_name = containers.Map('KeyType', 'char', 'ValueType', 'any');
fig_listing = dir(fullfile(replot_dir, 'figures', '*.fig'));
for k = 1:numel(fig_listing)
    f = openfig(fullfile(fig_listing(k).folder, fig_listing(k).name), 'invisible');
    by_name(get(f, 'Name')) = f;
end

% 3) One 1 x N sheet per measure (one column per condition), styled alike.
figs = gobjects(1, numel(specs));
tags = cell(1, numel(specs));
for mi = 1:numel(specs)
    spec = specs(mi);
    key = sprintf('%s Distribution', spec.field);
    if ~isKey(by_name, key)
        error('Fig_param_space_allStd:MissingFigure', ...
            'Expected a "%s" figure in:\n  %s', key, fullfile(replot_dir, 'figures'));
    end
    src_fig = by_name(key);
    src_ax  = sort_axes_left_to_right(src_fig);
    nCols   = numel(src_ax);
    combined = figure('Color', 'w', 'Position', [100 100 350 * nCols 300]);
    cax = gobjects(1, nCols);
    for c = 1:nCols
        ph = subplot(1, nCols, c, 'Parent', combined);
        target_pos = get(ph, 'Position');
        delete(ph);
        cax(c) = copyobj(src_ax(c), combined);
        set(cax(c), 'Position', target_pos);
    end
    close(src_fig);

    for c = 1:nCols
        ax = cax(c);
        set(ax, 'FontSize', tick_fs);
        set(ax.YLabel, 'FontSize', label_fs);
        xlabel(ax, spec.label, 'Interpreter', 'tex', 'FontSize', label_fs);
        set(ax.Title, 'FontWeight', 'normal', 'FontSize', title_fs);
    end
    if strcmp(cfg.lle_scale, 'symlog') && strcmp(spec.field, 'LLE')
        symlog_histograms(cax, run_dir, cfg.conditions, spec);
    end
    % Each column's bars in its condition's colour (Plan 2, package A,
    % 2026-11-21: no condition in gray); the column's title names it.
    cond_rgb = condition_rgb_by_column(cax, reshape([cax.Position], 4, [])', [0.6 0.6 0.6], st);
    for c = 1:nCols
        set(findobj(cax(c), 'Type', 'histogram'), 'FaceColor', cond_rgb(c, :));
        set(findobj(cax(c), 'Type', 'bar'), 'FaceColor', cond_rgb(c, :));
    end
    linkaxes(cax, 'y');
    % The 300-px canvas leaves no room for the 20-pt x label, which printed cut
    % in half (on product runs too); extend the canvas downward by what it
    % needs (2026-11-12). Before the dividers, which are annotations and would
    % not follow.
    extend_canvas_below(combined, cax);

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
        x_div = (col_right(c) + col_left(c + 1)) / 2 - x_shift;
        annotation(combined, 'line', [x_div x_div], [y_bot y_top], ...
            'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
    end
    if ~cfg.visible; set(combined, 'Visible', 'off'); end
    figs(mi) = combined;
    tags{mi} = sprintf('Fig_ParamSpace_%s', spec.stem);
end
% Any leftover reopened figures (measures not in the sheet list) are closed.
for k = keys(by_name)
    if isgraphics(by_name(k{1})); close(by_name(k{1})); end
end

% 4) Save ONLY the combined figures, with STABLE names, then drop the prep folder.
out = struct('figs', figs, 'files', {{}}, 'source', run_dir);
if cfg.save
    for mi = 1:numel(specs)
        save_figure_stable(out_dir, tags{mi}, figs(mi));
        out.files = [out.files, existing_outputs(out_dir, tags{mi})];
    end
    if isfolder(replot_dir); rmdir(replot_dir, 's'); end
end
end

%% ------------------------------------------------------------------------
function symlog_histograms(cax, run_dir, conditions, spec)
% Replace the histograms of the lambda_1 sheet by ones of the sample's own
% values on a symmetric-log axis (see the header). cax: the sheet's axes,
% left to right, one per condition in the order replot_param_space_analysis
% drew them (select_conditions on the sample's own list).
C_LIN = 0.1;                                  % 1/s: linear within about this
N_BINS = 40;
sym = @(x) sign(x) .* log10(1 + abs(x) / C_LIN);
listing = dir(fullfile(run_dir, 'param_space_*'));
listing = listing([listing.isdir]);
values = {}; names = {};
for k = 1:numel(listing)
    d = fullfile(listing(k).folder, listing(k).name);
    if ~isfile(fullfile(d, 'psa_object.mat')); continue; end
    psa = ParamSpaceAnalysis2.from_dir(d);
    held = cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false);
    names = select_conditions(held, conditions);
    values = cellfun(@(n) grid_values(psa.results, n, spec.field), names, 'UniformOutput', false);
    break                                     % one joint sample per run, as the sheet's source
end
assert(numel(values) == numel(cax), 'fig_param_space_allStd:SymlogColumns', ...
    'The lambda_1 sheet has %d columns and the joint sample %d conditions to draw.', numel(cax), numel(values));
top = max([cellfun(@(v) max([abs(v(:)); 0]), values), 10]);
t_lim = sym(1.2 * top);
edges = linspace(-t_lim, t_lim, N_BINS + 1);
ticks = [-100 -10 -1 -0.1 0.1 1 10 100];
ticks = ticks(abs(sym(ticks)) <= t_lim);
for c = 1:numel(cax)
    ax = cax(c);
    delete(findobj(ax, 'Type', 'histogram'));
    delete(findobj(ax, 'Type', 'bar'));
    v = values{c};
    if ~isempty(v)
        hold(ax, 'on');
        h = histogram(ax, 'BinEdges', edges, 'BinCounts', histcounts(sym(v), edges) / numel(v), ...
            'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
        uistack(h, 'bottom');
        hold(ax, 'off');
    end
    set(ax, 'XLim', [-t_lim, t_lim], 'XTick', sym(ticks), 'XTickLabel', compose('%g', ticks), 'XTickLabelRotation', 0);
    ax.UserData = struct('values', v, 'condition', names{c}, 'scale', sprintf('symmetric log, linear within %g 1/s', C_LIN));
end
end

function vals = grid_values(results, cond_name, metric)
% Every successful value of the metric over the grid, one condition (the
% rule of ParamSpaceAnalysis2's own distribution plot, whose collector is
% private).
vals = [];
if ~isfield(results, cond_name); return; end
rc = results.(cond_name);
for k = 1:numel(rc)
    res = rc{k};
    if isstruct(res) && isfield(res, 'success') && res.success && isfield(res, metric) && ~isnan(res.(metric))
        vals(end + 1) = res.(metric); %#ok<AGROW>
    end
end
end
