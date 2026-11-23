function out = fig_main4_connectivity(cfg)
% FIG_MAIN4_CONNECTIVITY Main Figure 4: the sweeps (A), lambda_1 over the joint sample by E:I weight ratio (B), lambda_1 against rate (C).
%
%   out = FIG_MAIN4_CONNECTIVITY('run_dir', d, 'preset_name', p)
%
% The paper's Figure 4 as Tom laid it out on 2026-11-22
% (Plan2_figure_feedback.md, Figure 4, notes 2 to 7; Plan3.md, phase 0),
% three rows on the common canvas (canvas_w of manuscript_style):
%   (A) lambda_1 against each of six connectivity parameters, one varied
%       at a time (the 1-D sweeps), conditions overlaid: median over
%       repetitions and interquartile band (fig_sensitivity_medians, its
%       lambda_1 sheet). The y window is fixed and a curve that leaves it is
%       simply clipped: no off-scale triangles (note 2).
%   (B) the distribution of lambda_1 over the joint sample, one panel per
%       condition, each bin stacked by the E:I weight ratio of its networks
%       (fig_EI_weights_param_space, its lambda_1 sheet; note 6). As in that
%       figure, values beyond +-1.5 1/s are gathered into the two end bins.
%   (C) lambda_1 against mean population firing rate over the joint sample
%       ONLY, every dot one size (fig_lle_vs_rate with include_sweeps false;
%       note 5), colored by E:I weight ratio; the median of lambda_1 per bin
%       of rate in black, on top of the points over a white halo (note 3,
%       question 5), with a key saying so just above the multiple-timescale
%       panel (note 4).
% (B) and (C) share one E:I color bar at the right. The dendritic share of
% the leading perturbation (fig_leading_direction, the old panel C) has left
% this figure (note 7); that function is unchanged.
%
% HOW IT IS BUILT. With composite_tools, as the supplementary composites:
% each source function draws its panels invisibly and unsaved, the axes are
% copied onto one canvas and placed in pixels, laid out on canvas_w itself
% (not stretched), with the type by role (T.style). Nothing is simulated.
% Until Plan 3 this function joined fig_grouped_main(4) (panels A and B of
% that time, on a 13 in canvas) and fig_leading_direction's column with
% composite_figure; fig_grouped_main's group 4 is no longer part of it.
%
% OPTIONS
%   run_dir, preset_name   the run and its preset (the sweeps and the weight
%                    balances are read through the preset)
%   conditions       passed to every source figure (FD-006); empty = the
%                    preset's / the run's
%   fig_tag          the file name stem, default 'Fig_Main4_Connectivity'.
%                    The registry passes 'Fig_Main4_Grouped' so that the
%                    manuscript's figure_sources.json needs no repointing.
%
% Saves <tag>.png, .fig and <tag>_table.md (composite_save). out.figs,
% out.files, out.source, out.panels (letter -> source function).
%
% See also: fig_sensitivity_medians, fig_EI_weights_param_space,
%           fig_lle_vs_rate, composite_tools, composite_save

arguments
    cfg.verbose        (1,:) char    = 'minimal'
    cfg.data_file      (1,:) char    = ''     % unused; accepted for a uniform call
    cfg.out_dir        (1,:) char    = ''
    cfg.save           (1,1) logical = true
    cfg.visible        (1,1) logical = false
    cfg.run_dir        (1,:) char    = ''
    cfg.preset_name    (1,:) char    = ''
    cfg.conditions                   = {}
    cfg.fig_tag        (1,:) char    = 'Fig_Main4_Connectivity'
end

setup_paths();
T = composite_tools();
st = manuscript_style();
quiet = {'save', false, 'visible', false, 'verbose', cfg.verbose};
args = [{'run_dir', cfg.run_dir, 'conditions', cfg.conditions}, quiet];
if ~isempty(cfg.preset_name); args = [args, {'preset_name', cfg.preset_name}]; end

%% The source figures
specs = sweep_metrics();
specs = specs([specs.in_sheets]);
lle_sheet = find(strcmp({specs.field}, 'LLE'));
assert(isscalar(lle_sheet), 'fig_main4_connectivity:Sheets', 'sweep_metrics has no lambda_1 sheet.');
S_sens = fig_sensitivity_medians(args{:});
cleanup_sens = onCleanup(@() T.close(S_sens.figs));
S_ei = fig_EI_weights_param_space(args{:});
cleanup_ei = onCleanup(@() T.close(S_ei.figs));
S_rate = fig_lle_vs_rate('run_dir', cfg.run_dir, 'conditions', cfg.conditions, 'include_sweeps', false, quiet{:});
cleanup_rate = onCleanup(@() T.close(S_rate.figs));
assert(numel(S_sens.figs) >= lle_sheet && numel(S_ei.figs) >= lle_sheet, 'fig_main4_connectivity:Sheets', ...
    'A source figure returned no lambda_1 sheet.');

%% The canvas, laid out on canvas_w itself (pixels, origin at the bottom left)
W = round(st.canvas_w * get(groot, 'ScreenPixelsPerInch'));
H = 950;
fig = T.canvas(W, H, cfg.visible);
X0 = 72;                        % room for the y labels of the left column
% (A): 2 x 3; the first gap is the wider one, where the last E:I neuron
% count (400:100) meets the first network size (100).
A_GAP = [66 50]; A_H = 125;
A_W = (W - X0 - 8 - sum(A_GAP)) / 3;
A_Y = [768 568];                % bottom of the upper and the lower row
% (B) and (C): one panel per condition, the color bar at the right.
CB_W = 92; BC_GAP = 46;
B_Y = 340; B_H = 120;
C_Y = 62;  C_H = 150;

%% (A) the sweeps
src = S_sens.figs(lle_sheet);
G = T.grid(findall(src, 'Type', 'axes'));
assert(isequal(size(G), [2 3]), 'fig_main4_connectivity:Panels', ...
    'Expected the 2 x 3 panels of fig_sensitivity_medians, found %d x %d.', size(G, 1), size(G, 2));
y_label = G(1, 1).YLabel.String;
Cp = T.copy(reshape(G', [], 1), fig, true);      % row by row; the sheet's legend comes along
axA = Cp.axes;
for k = 1:6
    r = ceil(k / 3); c = k - 3 * (r - 1);
    axA(k).Position = [X0 + (c - 1) * A_W + sum(A_GAP(1:c - 1)), A_Y(r), A_W, A_H];
    ylabel(axA(k), '');
    axA(k).XTickLabelRotation = 0;
end
ylabel(axA(1), y_label, 'Interpreter', 'tex');
ylabel(axA(4), y_label, 'Interpreter', 'tex');
names = legend_conditions(Cp.legends, cfg.conditions);
delete(Cp.legends);

%% (B) lambda_1 over the joint sample by E:I weight ratio
Cb = T.copy(S_ei.figs(lle_sheet), fig);
is_bar = arrayfun(@(h) ~isempty(findall(h, 'Type', 'image')), Cb.axes);
bar_axes = Cb.axes(is_bar);
assert(isscalar(bar_axes), 'fig_main4_connectivity:ColorBar', ...
    'Expected the one color bar of fig_EI_weights_param_space''s sheet, found %d.', numel(bar_axes));
axB = Cb.axes(~is_bar);
[~, o] = sort(arrayfun(@(h) h.Position(1), axB)); axB = axB(o);
nc = numel(axB);
BC_W = (W - X0 - CB_W - (nc - 1) * BC_GAP) / nc;

%% (C) lambda_1 against rate, the joint sample
rate_axes = findall(S_rate.figs, 'Type', 'axes');
Cc = T.copy(rate_axes, fig, false);              % the axes only; the color bar is made below
axC = Cc.axes;
[~, o] = sort(arrayfun(@(h) h.Position(1), axC)); axC = axC(o);
assert(numel(axB) == numel(names) && numel(axC) == numel(names), 'fig_main4_connectivity:Columns', ...
    'The sources drew %d conditions in (A), %d in (B) and %d in (C).', numel(names), numel(axB), numel(axC));
for c = 1:nc
    axB(c).Position = [X0 + (c - 1) * (BC_W + BC_GAP), B_Y, BC_W, B_H];
    axC(c).Position = [X0 + (c - 1) * (BC_W + BC_GAP), C_Y, BC_W, C_H];
end

%% Style
T.style(fig);
for c = 1:nc
    for a = [axB(c), axC(c)]
        ylabel(a, ''); box(a, 'off');
    end
    xlabel(axB(c), [st.label_lle ' (1/s)'], 'Interpreter', 'tex');
    xlabel(axC(c), st.label_rate);
end
ylabel(axB(1), 'Probability');
ylabel(axC(1), [st.label_lle ' (1/s)'], 'Interpreter', 'tex');
T.retitle(axB, names, true);
T.retitle(axC, names, true);
% One E:I color bar for (B) and (C), the height of both rows.
cb = ei_colorbar(axC(end), bar_axes, [W - CB_W + 10, C_Y, 13, B_Y + B_H - C_Y]);
delete(bar_axes);
% The legend of (A): the conditions and the zero line, two columns above
% the block, right-aligned.
lines = flipud(findall(axA(1), 'Type', 'line'));
lines = lines(arrayfun(@(h) ~isempty(h.DisplayName), lines));
zero = findall(axA(1), 'Type', 'constantline');
zero = zero(arrayfun(@(h) h.Value == 0, zero));
assert(numel(lines) == nc && ~isempty(zero), 'fig_main4_connectivity:Legend', ...
    'Panel (A)''s first axes holds %d condition lines and %d zero lines.', numel(lines), numel(zero));
zero = zero(1);
zero.HandleVisibility = 'on';
lgA = legend(axA(1), [lines; zero], [cellfun(@(n) T.condition_title(n), names, 'UniformOutput', false), {'Zero growth'}], ...
    'NumColumns', 2, 'Box', 'off', 'Units', 'pixels', 'Location', 'none', 'AutoUpdate', 'off', 'Tag', 'main4_legend_A');
lgA.ItemTokenSize = [22 18];
% The key of the median in (C), just above the multiple-timescale panel.
mts = find(contains(cellfun(@(n) T.condition_title(n), names, 'UniformOutput', false), 'Multiple'), 1);
if isempty(mts); mts = nc; end
med = findall(axC(mts), 'Tag', 'rate_bin_median');
assert(isscalar(med), 'fig_main4_connectivity:Median', 'Panel (C) has no rate-bin median line.');
lgC = legend(axC(mts), med, {'Median per rate bin'}, 'Box', 'off', 'Units', 'pixels', 'Location', 'none', ...
    'AutoUpdate', 'off', 'Tag', 'main4_legend_C');
lgC.ItemTokenSize = [26 18];
T.margin(fig, 6);
T.style(fig);                % again: the labels and titles set since take their roles' sizes
T.fit(fig);                  % laid out on canvas_w already: nothing moves
for pass = 1:2
    drawnow;
    p = lgA.Position;
    lgA.Position = [W - 8 - p(3), H - 6 - p(4), p(3), p(4)];
    p = lgC.Position; q = axC(mts).Position;
    lgC.Position = [q(1) + q(3) - p(3), q(2) + q(4) + title_height(axC(mts)) + 2, p(3), p(4)];
end

%% Panel letters, at the height of the legend of (A) and of the titles of (B) and (C)
T.letter(axA(1), 'A', 4 - X0, H - 6 - lgA.Position(4) - (A_Y(1) + A_H) + 2);
T.letter(axB(1), 'B', 4 - X0, title_height(axB(1)) - 18);
T.letter(axC(1), 'C', 4 - X0, title_height(axC(1)) - 18);

%% Save
panels = struct( ...
    'A', 'fig_sensitivity_medians (lambda_1 sheet): the six 1-D sweeps, conditions overlaid', ...
    'B', 'fig_EI_weights_param_space (lambda_1 sheet): lambda_1 over the joint sample, stacked by E:I weight ratio', ...
    'C', 'fig_lle_vs_rate (include_sweeps false): lambda_1 against mean firing rate, the joint sample');
sources = {S_sens.source, S_ei.source, S_rate.source};
parts = cell(0, 2);
if cfg.save
    what = {'E:I neuron ratio', 'network size', 'E to E', 'I to E', 'E to I', 'I to I'};
    for k = 1:6
        parts(end + 1, :) = {sprintf('Panel A, %s (fig_sensitivity_medians)', what{k}), ...
            T.axes_table(axA(k), 'Median over repetitions of lambda_1 at each level, per condition')}; %#ok<AGROW>
    end
    psa = ParamSpaceAnalysis2.from_dir(S_ei.source);
    for c = 1:nc
        parts(end + 1, :) = {sprintf('Panel B, %s (fig_EI_weights_param_space)', names{c}), lle_summary(psa, names{c})}; %#ok<AGROW>
    end
    for c = 1:nc
        parts(end + 1, :) = {sprintf('Panel C, %s (fig_lle_vs_rate)', names{c}), ...
            T.axes_table(axC(c), 'lambda_1 against mean firing rate, one point per network of the joint sample, and the median per rate bin')}; %#ok<AGROW>
    end
end
caption = sprintf(['Main Figure 4, connectivity (drawn by fig_main4_connectivity; run `%s`; conditions %s). ' ...
    '(A) The leading Lyapunov exponent lambda_1 against each connectivity parameter, one varied at a time ' ...
    '(1-D sweeps): median over repetitions (line) and interquartile range (band), conditions overlaid; dashed: ' ...
    'zero growth; the short mark on the x axis: the preset''s value. The y window is fixed, so a curve that ' ...
    'leaves it is clipped. (B) Distribution of lambda_1 over the jointly varied configurations, one panel per ' ...
    'condition (probability per bin; values beyond +-1.5 1/s are counted in the end bins), each bin stacked by ' ...
    'the E:I weight ratio of its networks (color bar). (C) lambda_1 against mean population firing rate over the ' ...
    'same configurations, one point per network colored by its E:I weight ratio; black: the median of lambda_1 ' ...
    'in each bin of rate (bins of 0.1 holding at least three networks); dashed: zero growth. Canvas %d x %d px ' ...
    '(canvas_w of manuscript_style), type in points by role.'], ...
    cfg.run_dir, strjoin(names, ', '), getappdata(fig, 'composite_size'));
out_dir = '';
if cfg.save; out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath')); end
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Main Figure 4: connectivity', caption, parts, sources);
out.panels = panels;
out.colorbar = cb;
end

%% ------------------------------------------------------------------------
function h = title_height(ax)
% The height in pixels of the title above the plot box of ax (0 without one).
h = 0;
if isempty(ax.Title.String); return; end
u = ax.Title.Units; ax.Title.Units = 'pixels';
e = ax.Title.Extent;
ax.Title.Units = u;
p = getpixelposition(ax);
h = max(0, e(2) + e(4) - p(4));
end

function cb = ei_colorbar(ax, bar_axes, rect)
% A color bar of ax at rect (pixels) with the map, limits, ticks and label
% of the sheet's image axes bar_axes (fig_EI_weights_param_space draws its
% color bar as an image; fig_lle_vs_rate uses the same map and limits).
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

function txt = lle_summary(psa, name)
% lambda_1 of one condition over the joint sample, for the table: what the
% stacked histogram of (B) is made of, and how much of it the end bins hold.
v = [];
if isfield(psa.results, name)
    rc = psa.results.(name);
    for k = 1:numel(rc)
        r = rc{k};
        if isstruct(r) && isfield(r, 'success') && r.success && isfield(r, 'LLE') && isfinite(r.LLE)
            v(end + 1) = r.LLE; %#ok<AGROW>
        end
    end
end
if isempty(v); txt = sprintf('No successful configuration of %s.\n', name); return; end
q = prctile(v, [25 50 75]);
txt = sprintf(['lambda_1 over the %d configurations of `%s` (1/s): median %.3g, quartiles %.3g to %.3g, ' ...
    'range %.3g to %.3g; share below zero %.3g; share in the end bins (beyond +-1.5): %.3g below, %.3g above.\n'], ...
    numel(v), name, q(2), q(1), q(3), min(v), max(v), mean(v < 0), mean(v < -1.5), mean(v > 1.5));
end

function names = legend_conditions(legends, requested)
% The conditions the sweeps drew, in their order: from the sheet's legend
% (condition titles, by name), else the selection as given.
names = {};
if ~isempty(legends)
    names = grouped_condition_label(legends(1).String, 'name_of_title');
end
if isempty(names) && ~isempty(requested); names = cellstr(requested); end
assert(~isempty(names), 'fig_main4_connectivity:Conditions', ...
    'The conditions fig_sensitivity_medians drew could not be read from its legend. Pass ''conditions''.');
names = reshape(cellstr(names), 1, []);
end
