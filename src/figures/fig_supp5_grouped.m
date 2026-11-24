function out = fig_supp5_grouped(cfg)
% FIG_SUPP5_GROUPED Supplementary Figure S5: the joint sample; saved data only.
%
%   out = FIG_SUPP5_GROUPED('run_dir', d, 'preset_name', p, 'out_dir', o)
%
% The figure plan of 2026-11-21 (figures_main_supp_appen.md, S5), revised
% after Tom's review of 2026-11-22 (Plan2_figure_feedback.md, S5; Plan3.md,
% phase 0, question 12). One column per condition; every panel is drawn by
% the function that owns it, invisibly and unsaved, and copied here
% (composite_tools).
%   (A) the distribution of the mean firing rate over the jointly varied
%       configurations (the param_space sample), stacked by the E:I WEIGHT
%       ratio of each network: the mean-rate sheet of
%       fig_EI_weights_param_space, as it was shown before the flat
%       condition colors (note 1), with its E:I color bar.
%   (B) local against accumulated growth rates: fig_local_vs_finite_lle
%       (the local rate pooled over time and networks, light, under the
%       finite-time lambda_1 per network, dark). NOTE that this function
%       draws the networks NEAR THE PRESET (every 1-D sweep at its level
%       nearest the preset's value), not the joint sample, which it reads
%       for its table only; the panel is the source's as it stands.
% The histogram of lambda_1 over the joint sample, colored by E:I weight
% ratio, is main Figure 4B since Plan 3, so it is not repeated here (it was
% panel A, on a symmetric-log axis in flat condition colors, until then).
%
% Options: run_dir, preset_name, conditions (empty = all the run holds),
% out_dir, save, visible, fig_tag.
%
% out: figs, files, source, panels (letter -> source function).
%
% See also: composite_tools, composite_save, fig_EI_weights_param_space,
%           fig_local_vs_finite_lle, fig_main4_connectivity

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
ei_args = [{'run_dir', cfg.run_dir, 'conditions', cfg.conditions}, quiet];
if ~isempty(cfg.preset_name); ei_args = [ei_args, {'preset_name', cfg.preset_name}]; end

%% The source figures
S_ei = fig_EI_weights_param_space(ei_args{:});
cleanup_ei = onCleanup(@() T.close(S_ei.figs));
S_loc = fig_local_vs_finite_lle('run_dir', cfg.run_dir, 'conditions', cfg.conditions, quiet{:});
cleanup_loc = onCleanup(@() T.close(S_loc.figs));
specs = sweep_metrics();
specs = specs([specs.in_sheets]);
sheet = find(strcmp({specs.field}, 'mean_rate'));
assert(isscalar(sheet) && numel(S_ei.figs) >= sheet, 'fig_supp5_grouped:Sheets', ...
    'fig_EI_weights_param_space returned no mean-rate sheet.');

%% The canvas: two rows of one panel per condition, laid out on canvas_w itself
W = round(st.canvas_w * get(groot, 'ScreenPixelsPerInch'));
H = 590;
fig = T.canvas(W, H, cfg.visible);
X0 = 78; GAP = 46; RH = 170;
CB_W = 92;                     % the color bar of (A), its ticks and its label, on the right
ROW_Y = [345 72];
ax = cell(1, 2);
src = {S_ei.figs(sheet), S_loc.figs};
nc = 0; names = {}; bar_axes = gobjects(0);
for r = 1:2
    Cp = T.copy(src{r}, fig);
    a = Cp.axes;
    if r == 1
        % The sheet's color bar is an image axes of its own: kept aside,
        % replaced below by a color bar of the same ticks.
        is_bar = arrayfun(@(h) ~isempty(findall(h, 'Type', 'image')), a);
        bar_axes = a(is_bar);
        a = a(~is_bar);
    end
    [~, order] = sort(arrayfun(@(h) h.Position(1), a));
    ax{r} = a(order);
    if r == 1
        nc = numel(ax{r});
        CW = (W - X0 - CB_W - (nc - 1) * GAP) / nc;
        names = condition_names(ax{r}, cfg.conditions);
    end
    assert(numel(ax{r}) == nc, 'fig_supp5_grouped:Columns', ...
        'Row %d has %d panels, row 1 has %d: the source figures drew different conditions.', r, numel(ax{r}), nc);
    for c = 1:nc
        ax{r}(c).Position = [X0 + (c - 1) * (CW + GAP), ROW_Y(r), CW, RH];
    end
    if r == 2
        % One key for the two histograms of (B), made below: the source's
        % own (one per panel, inside it) go.
        delete(Cp.legends);
    end
end
assert(isscalar(bar_axes), 'fig_supp5_grouped:ColorBar', ...
    'Expected the one color bar of fig_EI_weights_param_space''s sheet, found %d.', numel(bar_axes));

%% Style
T.style(fig);
for r = 1:2
    for c = 1:nc
        title(ax{r}(c), '');
        ylabel(ax{r}(c), '');
        box(ax{r}(c), 'off');
    end
end
T.retitle(ax{1}, names, true);
ylabel(ax{1}(1), 'Probability');
ylabel(ax{2}(1), 'Density');
for c = 1:nc
    xlabel(ax{1}(c), st.label_rate);
    xlabel(ax{2}(c), 'Growth rate (s^{-1})', 'Interpreter', 'tex');
end
% The color bar of (A): the sheet's map, limits, ticks and label, as a
% color bar of the last panel at the right of the row.
cb = ei_colorbar(ax{1}(end), bar_axes, [W - CB_W + 10, ROW_Y(1), 13, RH]);
delete(bar_axes);
% The key of (B): light for the local rate, dark for the accumulated one,
% in neutral grays (it serves every column, each in its own colour), on one
% line ABOVE the first panel, where it covers no data.
hold(ax{2}(1), 'on');
key = [patch(ax{2}(1), NaN, NaN, [0.82 0.82 0.82], 'EdgeColor', 'none'), ...
       patch(ax{2}(1), NaN, NaN, [0.40 0.40 0.40], 'EdgeColor', 'none')];
hold(ax{2}(1), 'off');
legends = legend(ax{2}(1), key, {'Local', 'Accumulated'}, 'Box', 'off', 'FontSize', st.tick_fs, ...   % the caption says what is pooled
    'Orientation', 'horizontal', 'NumColumns', 2, 'Units', 'pixels', 'Location', 'none', 'AutoUpdate', 'off');
legends.ItemTokenSize = [22 18];
for pass = 1:2
    drawnow;
    p = legends.Position;
    q = ax{2}(1).Position;
    legends.Position = [q(1), q(2) + q(4) + 4, p(3), p(4)];
end
T.margin(fig, 6);
T.style(fig);                % again: the labels and titles set since take their roles' sizes
T.fit(fig);                  % laid out on canvas_w already: nothing moves

%% Panel letters, at the height of the titles of (A) and of the key of (B)
T.letter(ax{1}(1), 'A', 4 - X0, 44);
T.letter(ax{2}(1), 'B', 4 - X0, 8);

%% Save
panels = struct( ...
    'A', ['fig_EI_weights_param_space (mean-rate sheet): distribution of the mean firing rate over the joint sample, ' ...
          'stacked by E:I weight ratio, per condition'], ...
    'B', 'fig_local_vs_finite_lle: local against accumulated growth rates, the networks near the preset');
sources = {S_ei.source, S_loc.source};
parts = cell(0, 2);
if cfg.save
    psa = ParamSpaceAnalysis2.from_dir(S_ei.source);
    for c = 1:nc
        parts(end + 1, :) = {sprintf('Panel A, %s (fig_EI_weights_param_space)', names{c}), ...
            rate_summary(psa, names{c})}; %#ok<AGROW>
    end
    for c = 1:nc
        parts(end + 1, :) = {sprintf('Panel B, %s (fig_local_vs_finite_lle)', names{c}), ...
            T.axes_table(ax{2}(c), 'Densities of the local rate (pooled) and of the finite-time lambda_1 per network')}; %#ok<AGROW>
    end
end
caption = sprintf(['Supplementary Figure S5, the joint sample (drawn by fig_supp5_grouped; run `%s`; conditions %s), ' ...
    'one column per condition. (A) Distribution of the mean firing rate over the jointly varied configurations ' ...
    '(probability per bin), each bin stacked by the E:I weight ratio of its networks (the summed excitatory ' ...
    'against the summed inhibitory weight, color bar; ratios beyond 1:10 and 10:1 take the end colors). ' ...
    '(B) Local against accumulated growth rates for the networks near the preset (every 1-D sweep at its level ' ...
    'nearest the preset''s value): density of the local rate of the leading direction pooled over time and ' ...
    'networks (light) and of the finite-time lambda_1 per network (dark; its median is in the table below); dotted: ' ...
    'zero; the key above the first panel is in gray and holds for every column. The distribution of lambda_1 ' ...
    'over the joint sample is main Figure 4B. Canvas %d x %d px (canvas_w of manuscript_style), type in points by role.'], ...
    cfg.run_dir, strjoin(names, ', '), getappdata(fig, 'composite_size'));
out_dir = '';
if cfg.save; out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath')); end
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Supplementary Figure S5: the joint sample', caption, parts, sources);
out.panels = panels;
out.colorbar = cb;
end

%% ------------------------------------------------------------------------
function cb = ei_colorbar(ax, bar_axes, rect)
% A color bar of ax at rect (pixels) with the map, limits, ticks and label
% of the sheet's image axes bar_axes (fig_EI_weights_param_space draws its
% color bar as an image).
img = findall(bar_axes, 'Type', 'image');
colormap(ax, colormap(bar_axes));
clim(ax, [min(img(1).YData), max(img(1).YData)]);
cb = colorbar(ax, 'Units', 'pixels');
cb.Location = 'manual';
cb.Position = rect;
cb.Ticks = bar_axes.YTick;
cb.TickLabels = bar_axes.YTickLabel;
cb.Label.String = bar_axes.YLabel.String;
cb.Box = 'off';
end

function txt = rate_summary(psa, name)
% The mean firing rates of one condition over the joint sample, as a
% paragraph for the table: what the stacked histogram of (A) is made of.
v = [];
if isfield(psa.results, name)
    rc = psa.results.(name);
    for k = 1:numel(rc)
        r = rc{k};
        if isstruct(r) && isfield(r, 'success') && r.success && isfield(r, 'mean_rate') && isfinite(r.mean_rate)
            v(end + 1) = r.mean_rate; %#ok<AGROW>
        end
    end
end
if isempty(v); txt = sprintf('No successful configuration of %s.\n', name); return; end
q = prctile(v, [25 50 75]);
txt = sprintf(['Mean firing rate over the %d configurations of `%s`: median %.3g, quartiles %.3g to %.3g, ' ...
    'range %.3g to %.3g; share below 0.02 (quiet) %.3g, above 0.9 (saturated) %.3g. The color of each ' ...
    'network is its E:I weight ratio (fig_EI_weights_param_space).\n'], ...
    numel(v), name, q(2), q(1), q(3), min(v), max(v), mean(v < 0.02), mean(v > 0.9));
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
