function out = fig_supp5_grouped(cfg)
% FIG_SUPP5_GROUPED Supplementary Figure S5: the joint sample; saved data only.
%
%   out = FIG_SUPP5_GROUPED('run_dir', d, 'preset_name', p, 'out_dir', o)
%
% The figure plan of 2026-11-21 (figures_main_supp_appen.md, S5), for the
% medians and the share of contracting configurations that the Results
% quote. One column per condition; every panel is drawn by the function
% that owns it, invisibly and unsaved, and copied here (composite_tools).
%   (A) the distribution of lambda_1 over the jointly varied configurations
%       (the param_space sample): the lambda_1 sheet of
%       fig_param_space_allStd, asked for on its SYMMETRIC-LOG axis
%       ('lle_scale', 'symlog': linear within about +-0.1 1/s, logarithmic
%       beyond, every value at its place). Until 2026-11-22 it was that
%       sheet's window of +-1.5 1/s with everything beyond gathered into
%       the two end bins, which hid the distribution without adaptation
%       (three quarters of it in the left end bin). The table gives each
%       column's quartiles and the share below zero.
%   (B) the distribution of the mean firing rate over the same
%       configurations: its mean-rate sheet.
%   (C) local against accumulated growth rates: fig_local_vs_finite_lle
%       (the local rate pooled over time and networks, light, under the
%       finite-time lambda_1 per network, dark). NOTE that this function
%       draws the networks NEAR THE PRESET (every 1-D sweep at its level
%       nearest the preset's value), not the joint sample, which it reads
%       for its table only; the panel is the source's as it stands.
% The histograms of (A) and (B) are gray in the source; here each takes its
% condition's colour (manuscript_style, by name), as the plan's style rule
% asks, so the three rows agree.
%
% THE RUN DIRECTORY IS WRITTEN TO, briefly: fig_param_space_allStd makes
% its histograms in a prep folder <run_dir>/replot_param_space_<date> and,
% called without save, leaves it there. This function removes the prep
% folders that appeared during its call.
%
% Options: run_dir, preset_name, conditions (empty = all the run holds),
% out_dir, save, visible, fig_tag.
%
% out: figs, files, source, panels (letter -> source function).
%
% See also: composite_tools, composite_save, fig_param_space_allStd,
%           fig_local_vs_finite_lle, sweep_metrics

arguments
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''
    cfg.conditions                = {}
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = false
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.fig_tag     (1,:) char    = 'Fig_Supp5_Grouped'
end
setup_paths();
T = composite_tools();
st = manuscript_style();
quiet = {'save', false, 'visible', false, 'verbose', cfg.verbose};
ps_args = [{'run_dir', cfg.run_dir, 'conditions', cfg.conditions}, quiet];
if ~isempty(cfg.preset_name); ps_args = [ps_args, {'preset_name', cfg.preset_name}]; end

%% The source figures
prep_before = prep_folders(cfg.run_dir);
tidy = onCleanup(@() remove_new_prep(cfg.run_dir, prep_before));
S_ps  = fig_param_space_allStd(ps_args{:}, 'lle_scale', 'symlog');
cleanup_ps = onCleanup(@() T.close(S_ps.figs));
S_loc = fig_local_vs_finite_lle('run_dir', cfg.run_dir, 'conditions', cfg.conditions, quiet{:});
cleanup_loc = onCleanup(@() T.close(S_loc.figs));
specs = sweep_metrics();
specs = specs([specs.in_sheets]);
measures = {'LLE', 'mean_rate'};
[found, sheet] = ismember(measures, {specs.field});
assert(all(found) && numel(S_ps.figs) >= max(sheet), 'fig_supp5_grouped:Sheets', ...
    'fig_param_space_allStd returned no sheet for %s.', strjoin(measures(~found), ', '));

%% The canvas: three rows of one panel per condition
W = 1100; H = 890;             % under 0.9 of the workstation screen
fig = T.canvas(W, H, cfg.visible);
X0 = 100; GAP = 50; RH = 170;
ROW_Y = [650 380 80];          % 30 px more above (C) than between (A) and (B): its key sits there
ax = cell(1, 3);
src = {S_ps.figs(sheet(1)), S_ps.figs(sheet(2)), S_loc.figs};
nc = 0;
for r = 1:3
    Cp = T.copy(src{r}, fig);
    [~, order] = sort(arrayfun(@(a) a.Position(1), Cp.axes));
    ax{r} = Cp.axes(order);
    if r == 1
        nc = numel(ax{r});
        CW = (W - X0 - 15 - (nc - 1) * GAP) / nc;
        names = condition_names(ax{r}, cfg.conditions);
    end
    assert(numel(ax{r}) == nc, 'fig_supp5_grouped:Columns', ...
        'Row %d has %d panels, row 1 has %d: the source figures drew different conditions.', r, numel(ax{r}), nc);
    for c = 1:nc
        ax{r}(c).Position = [X0 + (c - 1) * (CW + GAP), ROW_Y(r), CW, RH];
    end
    if r == 3
        % One key for the two histograms of (C), made below: the source's
        % own (one per panel, inside it) go.
        delete(Cp.legends);
    end
end

%% Style
T.style(fig);
for r = 1:3
    for c = 1:nc
        title(ax{r}(c), '');
        if c > 1; ylabel(ax{r}(c), ''); end
        box(ax{r}(c), 'off');
    end
end
T.retitle(ax{1}, names, true);
% (A), (B): each histogram in its condition's colour.
for r = 1:2
    for c = 1:nc
        h = findall(ax{r}(c), 'Type', 'histogram');
        if isKey(st.condition_color, names{c})
            set(h, 'FaceColor', st.condition_color(names{c}), 'EdgeColor', 'none', 'FaceAlpha', 0.75);
        end
    end
    ylabel(ax{r}(1), 'Probability');
end
for c = 1:nc
    xlabel(ax{1}(c), 'Lyapunov Exp., \lambda_1 (s^{-1})', 'Interpreter', 'tex');
    xlabel(ax{2}(c), 'Mean Firing Rate');
    xlabel(ax{3}(c), 'Growth Rate (s^{-1})', 'Interpreter', 'tex');
end
ylabel(ax{3}(1), 'Density');
% The key of (C): light for the local rate, dark for the accumulated one,
% in neutral grays (it serves the three columns, each in its own colour),
% on one line ABOVE the first panel. Until 2026-11-22 it was the first
% panel's own legend in its upper right corner, where it covered that
% panel's median line and tallest bars.
hold(ax{3}(1), 'on');
key = [patch(ax{3}(1), NaN, NaN, [0.82 0.82 0.82], 'EdgeColor', 'none'), ...
       patch(ax{3}(1), NaN, NaN, [0.40 0.40 0.40], 'EdgeColor', 'none')];
hold(ax{3}(1), 'off');
legends = legend(ax{3}(1), key, {'Local', 'Accumulated'}, 'Box', 'off', 'FontSize', st.tick_fs, ...   % the caption says what is pooled
    'Orientation', 'horizontal', 'NumColumns', 2, 'Units', 'pixels', 'Location', 'none', 'AutoUpdate', 'off');
legends.ItemTokenSize = [22 18];
for pass = 1:2
    drawnow;
    p = legends.Position;
    q = ax{3}(1).Position;
    legends.Position = [q(1), q(2) + q(4) + 4, p(3), p(4)];
end

%% The common canvas: the layout above, made on W px, stretched across canvas_w
T.style(fig);                % again: the labels and titles set since take their roles' sizes
T.fit(fig);

%% Panel letters
T.letter(ax{1}(1), 'A', -85, 42);
T.letter(ax{2}(1), 'B', -85, 18);
T.letter(ax{3}(1), 'C', -85, 18);

%% Save
panels = struct( ...
    'A', 'fig_param_space_allStd (lambda_1 sheet): distribution of lambda_1 over the joint sample, per condition', ...
    'B', 'fig_param_space_allStd (mean-rate sheet): distribution of the mean firing rate over the joint sample', ...
    'C', 'fig_local_vs_finite_lle: local against accumulated growth rates, the networks near the preset');
sources = {S_ps.source, S_loc.source};
parts = cell(0, 2);
if cfg.save
    f = fieldnames(panels);
    what = {'Histogram of lambda_1 over the configurations (bin centres and probabilities)', ...
        'Histogram of the mean firing rate over the configurations (bin centres and probabilities)', ...
        'Densities of the local rate (pooled) and of the finite-time lambda_1 per network'};
    for r = 1:3
        for c = 1:nc
            body = T.axes_table(ax{r}(c), what{r});
            if r == 1; body = [values_summary(ax{r}(c)), body]; end
            parts(end + 1, :) = {sprintf('Panel %s, %s (%s)', f{r}, names{c}, strtok(panels.(f{r}), ' ')), body}; %#ok<AGROW>
        end
    end
end
caption = sprintf(['Supplementary Figure S5, the joint sample (drawn by fig_supp5_grouped; run `%s`; conditions %s), ' ...
    'one column per condition. (A) Distribution of the leading Lyapunov exponent lambda_1 over the jointly varied ' ...
    'configurations, on a symmetric-log axis (linear within about 0.1 1/s of zero, logarithmic beyond; 40 bins of ' ...
    'equal width on that axis, probability per bin; every value at its place, none gathered; dashed: ' ...
    'lambda_1 = 0). (B) Distribution of the mean firing rate over the same configurations. (C) Local against ' ...
    'accumulated growth rates for the networks near the preset (every 1-D sweep at its level nearest the ' ...
    'preset''s value): density of the local rate of the leading direction pooled over time and networks (light) ' ...
    'and of the finite-time lambda_1 per network (dark), with its median as a line; dotted: zero; the key above ' ...
    'the first panel is in gray and holds for the three columns. The histograms of (A) and (B) are gray in ' ...
    'fig_param_space_allStd and take the condition colours here. Canvas %d x %d px (canvas_w of manuscript_style), type in points by role.'], ...
    cfg.run_dir, strjoin(names, ', '), getappdata(fig, 'composite_size'));
out_dir = '';
if cfg.save; out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath')); end
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Supplementary Figure S5: the joint sample', caption, parts, sources);
out.panels = panels;
end

%% ------------------------------------------------------------------------
function txt = values_summary(ax)
% The values behind a symmetric-log lambda_1 panel (fig_param_space_allStd
% leaves them in the axes' UserData), as a paragraph for the table; '' for
% a panel drawn the old way. The x of the series below it is in the axis'
% own coordinate, sign(x) log10(1 + |x| / 0.1).
txt = '';
u = ax.UserData;
if ~(isstruct(u) && isfield(u, 'values')) || isempty(u.values); return; end
v = u.values(:);
q = prctile(v, [25 50 75]);
txt = sprintf(['lambda_1 over the %d configurations of `%s` (1/s): median %.3g, quartiles %.3g to %.3g, range %.3g to ' ...
    '%.3g; share below zero %.3g; share beyond the old window of +-1.5: %.3g below, %.3g above. Axis: %s; the ' ...
    'x of the series listed next is in that coordinate, sign(x) log10(1 + |x| / 0.1).\n\n'], ...
    numel(v), u.condition, q(2), q(1), q(3), min(v), max(v), mean(v < 0), mean(v < -1.5), mean(v > 1.5), u.scale);
end

function names = condition_names(axs, requested)
% The conditions of the columns: the selection when one was given (the
% sources draw it in the order given), else read from the titles painted.
if ~isempty(requested) && ~((ischar(requested) || isstring(requested)) && all(strlength(requested) == 0))
    names = cellstr(requested);
    assert(numel(names) == numel(axs), 'fig_supp5_grouped:Conditions', ...
        '%d conditions selected, %d columns drawn.', numel(names), numel(axs));
    return
end
names = grouped_condition_label(arrayfun(@(a) {a.Title.String}, axs), 'name_of_title');
end

function d = prep_folders(run_dir)
listing = dir(fullfile(run_dir, 'replot_param_space_*'));
d = {listing([listing.isdir]).name};
end

function remove_new_prep(run_dir, before)
now_there = prep_folders(run_dir);
added = setdiff(now_there, before);
for k = 1:numel(added)
    d = fullfile(run_dir, added{k});
    if isfolder(d); rmdir(d, 's'); end
end
end
