function out = fig_param_space_allStd(cfg)
% FIG_PARAM_SPACE_ALLSTD Param-space LLE and firing-rate distributions, 2 x N.
%
%   out = FIG_PARAM_SPACE_ALLSTD()
%   out = FIG_PARAM_SPACE_ALLSTD('run_dir', d)
%
% Row 1 = LLE distributions, row 2 = mean firing-rate distributions, one column
% per adaptation condition. No simulation is re-run: the saved param-space PSA
% object is reloaded and its per-condition histogram axes are copied into one
% combined figure.
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
project_root = fileparts(which('setup_paths'));
st           = manuscript_style();

run_dir = resolve_run_dir('run_dir', cfg.run_dir, 'preset_name', cfg.preset_name);

% --- Presentation constants ------------------------------------------------
tick_fs  = st.tick_fs;
label_fs = st.label_fs;
title_fs = st.title_fs;
% Probability y ticks, one entry per row. The two rows span different ranges
% (the LLE distributions peak near 0.8, the rate distributions near 0.4), so
% each gets its own sparse set rather than MATLAB's automatic ticks.
prob_yticks = {[0, 0.4, 0.8], [0, 0.2, 0.4]};
x_shift = 0.007;   % nudge column dividers slightly left (normalized figure units)

% Start from a clean slate: replot_param_space_analysis saves ALL open figures,
% so any stray figure lingering in the session would pollute the save.
% (no 'close all force' -- see the header note; it destroyed sibling figures)

% 1) Regenerate the LLE + mean_rate distribution figures into a
%    replot_param_space_<dt>/figures/ folder under data_root.
replot_dir = replot_param_space_analysis(run_dir);

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

if ~cfg.visible; set(combined, 'Visible', 'off'); end

%% --- Save -------------------------------------------------------------------
fig_tag = 'Fig_ParamSpace_allStd';
out = struct('figs', combined, 'files', {{}}, 'source', run_dir);
if cfg.save
    save_figure_stable(out_dir, fig_tag, combined);
    out.files = existing_outputs(out_dir, fig_tag);

    % The prep figures exist only to build the final one; remove the whole
    % replot folder so no extra figs are left behind in the data directory.
    if isfolder(replot_dir); rmdir(replot_dir, 's'); end

    capture_git_provenance(out_dir, project_root);

    write_figure_readme(out_dir, struct( ...
        'tag',    'fig_param_space_allStd', ...
        'title',  'Stability_Manuscript figure: parameter-space distributions', ...
        'script', 'fig_param_space_allStd.m', ...
        'what',   ['A 2 x N sheet. Row 1 is the distribution of the largest ' ...
                   'Lyapunov exponent over the whole parameter grid, with a ' ...
                   'green dashed line at zero; row 2 is the distribution of ' ...
                   'mean firing rate. One column per adaptation condition: ' ...
                   'No Adaptation, SFA, STD, SFA + STD.'], ...
        'how',    ['Presentation replot -- no simulation is re-run. ' ...
                   'replot_param_space_analysis regenerates psa.plot for LLE ' ...
                   'and mean_rate into a temporary replot folder; those axes ' ...
                   'are copied into a single combined figure, restyled, and ' ...
                   'the prep folder is deleted. Cleanups: LLE xlabel to ' ...
                   'lambda_1; condition titles only on the top row and not ' ...
                   'bold; vertical grey column dividers; y-axes linked within ' ...
                   'each row with sparse probability ticks.'], ...
        'source', struct('run_dir', run_dir, 'preset', cfg.preset_name), ...
        'figures', {out.files}, ...
        'notes',  ['THE LLE HISTOGRAM RANGE IS FIXED AT [-1.5, 1.5] inside ' ...
                   'ParamSpaceAnalysis2.plot and cannot be set from here. On ' ...
                   'the aug_14 preset the param-space LLEs spanned roughly ' ...
                   '-10 to +4.8, so the outermost bins accumulated a large ' ...
                   'share of the distribution. Check that against the run ' ...
                   'actually plotted before quoting the shape of the tails.']));
end
end
