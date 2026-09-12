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
% See also: resolve_run_dir, replot_param_space_analysis, fig_EI_param_space

arguments
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
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
replot_dir = replot_param_space_analysis(run_dir);
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
    linkaxes(cax, 'y');

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
