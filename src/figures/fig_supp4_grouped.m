function out = fig_supp4_grouped(cfg)
% FIG_SUPP4_GROUPED Supplementary Figure S4: sensitivity to each connectivity parameter; saved data only.
%
%   out = FIG_SUPP4_GROUPED('run_dir', d, 'preset_name', p, 'out_dir', o)
%
% The figure plan of 2026-11-21 (figures_main_supp_appen.md, S4): the
% summary of the sensitivity sheets, for the gain sweep and the stretches of
% positive growth that the Results on Figure 4 cite.
%   (A) lambda_1 against each of the SEVEN swept parameters, conditions
%       overlaid: median over repetitions and interquartile band, as main
%       Figure 4A.
%   (B) the mean firing rate for the same seven.
% Top row of each block: E:I neuron ratio, network size, all weights (the
% overall gain, level_of_chaos); bottom row: the four block means, E to E,
% I to E, E to I, I to I.
%
% WHERE THE PANELS COME FROM. Six of the seven are the panels of
% fig_sensitivity_medians (its lambda_1 sheet and its mean-rate sheet),
% drawn by it invisibly and unsaved and copied here (composite_tools). That
% function leaves the gain sweep out (six panels tile its 2 x 3 sheet), and
% this package may not edit it, so the seventh panel is drawn here from the
% same saved sweep by the same rule (ParamSpaceAnalysis2.collect_level_values,
% percentiles 25, 50, 75 of the successful repetitions at each level), on
% the percent axis of the block-mean panels (apply_percent_axis), with the
% preset mark (mark_default_value). Its colours, line width, band opacity,
% zero line, y window and ticks are READ FROM A COPIED PANEL, so a change of
% style in fig_sensitivity_medians reaches it. If that function gains an
% option for the panels it draws, the seventh panel should come from it too.
%
% Options: run_dir, preset_name, conditions (empty = the preset's), out_dir,
% save, visible, fig_tag.
%
% out: figs, files, source, panels (letter -> source function).
%
% See also: composite_tools, composite_save, fig_sensitivity_medians,
%           fig_sensitivity_analysis_allStd, sweep_metrics

arguments
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''
    cfg.conditions                = {}
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = false
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.fig_tag     (1,:) char    = 'Fig_Supp4_Grouped'
end
setup_paths();
T = composite_tools();
GAIN = 'level_of_chaos';
args = {'run_dir', cfg.run_dir, 'conditions', cfg.conditions, 'verbose', cfg.verbose, 'save', false, 'visible', false};
if ~isempty(cfg.preset_name); args = [args, {'preset_name', cfg.preset_name}]; end

%% The source figure: one sheet per measure, in the order of sweep_metrics
S = fig_sensitivity_medians(args{:});
cleanup = onCleanup(@() T.close(S.figs));
specs = sweep_metrics();
specs = specs([specs.in_sheets]);
measures = {'LLE', 'mean_rate'};
[found, sheet] = ismember(measures, {specs.field});
assert(all(found) && numel(S.figs) >= max(sheet), 'fig_supp4_grouped:Sheets', ...
    'fig_sensitivity_medians returned no sheet for %s.', strjoin(measures(~found), ', '));

%% The canvas: two blocks of 2 x 4 slots
W = 1100; H = 935;
fig = T.canvas(W, H, cfg.visible);
X0 = 95; CW = 196; GAP = 68; RH = 140;
ROW_Y = [745 530; 295 80];                 % ROW_Y(block, row)
slot = @(b, r, c) [X0 + (c - 1) * (CW + GAP), ROW_Y(b, r), CW, RH];

ax = cell(1, 2);        % ax{b}(1:7): ratio, size, gain, EE, IE, EI, II
legends = gobjects(0);
names = {};
for b = 1:2
    src = S.figs(sheet(b));
    G = T.grid(findall(src, 'Type', 'axes'));
    assert(isequal(size(G), [2 3]), 'fig_supp4_grouped:Panels', ...
        'Expected the 2 x 3 panels of fig_sensitivity_medians, found %d x %d.', size(G, 1), size(G, 2));
    y_label = G(1, 1).YLabel.String;
    order = [G(1, 1); G(1, 2); G(1, 3); G(2, 1); G(2, 2); G(2, 3)];   % ratio, size, EE, IE, EI, II
    Cp = T.copy(order, fig, b == 1);         % the legend comes with the first block only
    a = Cp.axes;
    if b == 1
        legends = Cp.legends;
        % The conditions drawn, from the legend the source made (titles, by name).
        names = legend_conditions(legends, cfg.conditions);
    end
    a(1).Position = slot(b, 1, 1);
    a(2).Position = slot(b, 1, 2);
    for k = 1:4; a(2 + k).Position = slot(b, 2, k); end
    gain = gain_panel(fig, cfg.run_dir, specs(sheet(b)), a(3), names, GAIN);
    gain.Position = slot(b, 1, 3);
    ax{b} = [a(1); a(2); gain; a(3:6)];
    set(ax{b}, 'Units', 'pixels');
    for k = 1:7; ylabel(ax{b}(k), ''); ax{b}(k).XTickLabelRotation = 0; end
    % A label of words on two lines: at label_fs 'Mean Firing Rate' was
    % longer than the 140-px panel it names.
    if ischar(y_label) && ~contains(y_label, '\') && numel(strsplit(y_label)) > 1
        w = strsplit(y_label);
        y_label = {strjoin(w(1:floor(end / 2)), ' '), strjoin(w(floor(end / 2) + 1:end), ' ')};
    end
    ylabel(ax{b}(1), y_label, 'Interpreter', 'tex');
    ylabel(ax{b}(4), y_label, 'Interpreter', 'tex');
end

%% Style
T.style(fig);
if ~isempty(legends)
    lg = legends(1);
    if ~isempty(names)
        lg.String = cellfun(@(n) T.condition_title(n), names, 'UniformOutput', false);
    end
    set(lg, 'Units', 'pixels', 'Location', 'none', 'Orientation', 'vertical', 'NumColumns', 1);
    lg.ItemTokenSize = [18 18];
    for pass = 1:2
        drawnow;
        p = lg.Position;
        lg.Position = [W - 10 - p(3), ROW_Y(1, 1) + RH - p(4), p(3), p(4)];
    end
end

% The x labels of (B)'s bottom row are the lowest thing drawn, and the
% cropped export ended on them (flush with the image edge, 2026-11-22):
% keep a margin under them.
T.margin(fig, 6);

%% The common canvas: the layout above, made on W px, stretched across canvas_w
T.style(fig);                % again: the labels and titles set since take their roles' sizes
T.fit(fig);

%% Panel letters
T.letter(ax{1}(1), 'A', -85, 18);
T.letter(ax{2}(1), 'B', -85, 18);

%% Save
panels = struct( ...
    'A', ['fig_sensitivity_medians (lambda_1 sheet, six panels) and, drawn here by the same rule, the gain sweep ' ...
          '(level_of_chaos, ''All Weights'')'], ...
    'B', 'fig_sensitivity_medians (mean-rate sheet, six panels) and the gain sweep drawn here');
sources = {S.source};
parts = cell(0, 2);
if cfg.save
    what = {'E:I neuron ratio', 'network size', 'all weights (overall gain)', 'E to E', 'I to E', 'E to I', 'I to I'};
    from = {'fig_sensitivity_medians', 'fig_sensitivity_medians', 'drawn by fig_supp4_grouped', ...
        'fig_sensitivity_medians', 'fig_sensitivity_medians', 'fig_sensitivity_medians', 'fig_sensitivity_medians'};
    letters = 'AB';
    for b = 1:2
        for k = 1:7
            parts(end + 1, :) = {sprintf('Panel %s, %s (%s)', letters(b), what{k}, from{k}), ...
                T.axes_table(ax{b}(k), sprintf('Median over repetitions of %s at each level, per condition (x in the sweep''s own units)', ...
                measures{b}))}; %#ok<AGROW>
        end
    end
end
caption = sprintf(['Supplementary Figure S4, sensitivity to each connectivity parameter (drawn by fig_supp4_grouped; ' ...
    'run `%s`; conditions %s). (A) The leading Lyapunov exponent lambda_1 against each swept parameter, one ' ...
    'parameter varied at a time: median over repetitions (line) and interquartile range (band), conditions ' ...
    'overlaid; dashed: lambda_1 = 0; the short mark on the x axis: the preset''s value. Top row: E:I neuron ' ...
    'ratio, network size, all weights (the overall gain); bottom row: the mean of each block of the ' ...
    'connectivity, read presynaptic to postsynaptic, as percent departure from the preset. The y window is ' ...
    'fixed, so a curve that leaves it is clipped. (B) The mean firing rate for the same seven sweeps. The gain ' ...
    'panel of each block is drawn by this function from the saved sweep, by the rule of fig_sensitivity_medians ' ...
    '(which omits it), in the style of the panels beside it. Canvas %d x %d px (canvas_w of manuscript_style), type in points by role.'], ...
    cfg.run_dir, strjoin(names, ', '), getappdata(fig, 'composite_size'));
out_dir = '';
if cfg.save; out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath')); end
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Supplementary Figure S4: sensitivity to each connectivity parameter', ...
    caption, parts, sources);
out.panels = panels;
end

%% ------------------------------------------------------------------------
function names = legend_conditions(legends, requested)
% The conditions the source drew, in its order: from its legend's titles,
% else the selection as given.
names = {};
if ~isempty(legends)
    try
        names = grouped_condition_label(legends(1).String, 'name_of_title');
    catch err
        if ~any(strcmp(err.identifier, {'grouped_condition_label:UnknownTitle', 'grouped_condition_label:AmbiguousTitle'}))
            rethrow(err);
        end
    end
end
if isempty(names) && ~isempty(requested); names = cellstr(requested); end
assert(~isempty(names), 'fig_supp4_grouped:Conditions', ...
    ['The conditions fig_sensitivity_medians drew could not be read from its legend. ' ...
     'Pass ''conditions''.']);
end

function ax = gain_panel(fig, run_dir, spec, donor, names, param)
% The sweep of the overall gain, drawn like the copied panel donor (a
% block-mean panel): the same colours, line width, band, zero line, y
% window and ticks, read from donor's own objects.
psa = sweep_of(run_dir, param);
x = psa.param_vectors{strcmp(psa.grid_params, param)};
x = x(:)';
% The donor's condition lines (in the order drawn) and one of its bands.
lines = flipud(findall(donor, 'Type', 'line'));
is_cond = arrayfun(@(h) ~isempty(h.DisplayName), lines);
cond_lines = lines(is_cond);
assert(numel(cond_lines) == numel(names), 'fig_supp4_grouped:DonorLines', ...
    'The copied sensitivity panel holds %d condition lines for %d conditions.', numel(cond_lines), numel(names));
band = findall(donor, 'Type', 'patch');
alpha = 0.18;
if ~isempty(band); alpha = band(1).FaceAlpha; end
mark = lines(~is_cond);
zero = findall(donor, 'Type', 'constantline');

ax = axes(fig, 'Units', 'pixels');
hold(ax, 'on');
if spec.zero_line && ~isempty(zero)
    yline(ax, 0, zero(1).LineStyle, 'Color', zero(1).Color, 'LineWidth', zero(1).LineWidth, ...
        'Alpha', zero(1).Alpha, 'HandleVisibility', 'off');
end
P = nan(numel(names), numel(x), 3);
for c = 1:numel(names)
    for k = 1:numel(x)
        v = ParamSpaceAnalysis2.collect_level_values(psa, param, k, names{c}, spec.field);
        if ~isempty(v); P(c, k, :) = prctile(v, [25 50 75]); end
    end
end
for c = 1:numel(names)
    ok = all(isfinite(squeeze(P(c, :, [1 3]))), 2)';
    fill(ax, [x(ok), fliplr(x(ok))], [P(c, ok, 1), fliplr(P(c, ok, 3))], cond_lines(c).Color, ...
        'FaceAlpha', alpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
for c = 1:numel(names)
    plot(ax, x, P(c, :, 2), '-', 'Color', cond_lines(c).Color, 'LineWidth', cond_lines(c).LineWidth, ...
        'DisplayName', cond_lines(c).DisplayName);
end
% The off-scale markers that fig_sensitivity_medians draws on the panels
% beside it (Plan2_canvas_standard.md: this panel, drawn here, had none):
% a triangle on the top or bottom edge where a median has left the y
% window, at the donor's marker size and fig_sensitivity_medians' default
% spacing.
donor_marks = findall(donor, 'Tag', 'offscale_marker');
marker_size = 5;
if ~isempty(donor_marks); marker_size = donor_marks(1).MarkerSize; end
cond_spec = struct('name', names(:)', 'color', arrayfun(@(h) h.Color, cond_lines(:)', 'UniformOutput', false));
draw_offscale_markers(ax, x, arrayfun(@(c) P(c, :, 2), 1:numel(names), 'UniformOutput', false), cond_spec, ...
    donor.YLim, marker_size, 0.035);
hold(ax, 'off');
xlim(ax, [min(x), max(x)]);
ylim(ax, donor.YLim);
set(ax, 'YTick', donor.YTick, 'FontSize', donor.FontSize, 'Box', 'off', 'LineWidth', donor.LineWidth, ...
    'TickDir', donor.TickDir);
default = preset_default_values(run_dir, {param});
if isKey(default, param) && default(param) ~= 0
    apply_percent_axis(ax, default(param), 'All Weights', donor.XLabel.FontSize, [-50 0 50]);
    if ~isempty(mark)
        frac = diff(mark(1).YData) / diff(donor.YLim);
        mark_default_value(ax, default(param), mark(1).Color, mark(1).LineWidth, frac);
    end
else
    xlabel(ax, 'All Weights');
end
end

function psa = sweep_of(run_dir, param)
% The 1-D sweep of run_dir over param; which parameter a folder swept is
% read off its PSA object, never guessed from the folder's name.
listing = dir(fullfile(run_dir, '1D_sensitivity_*'));
listing = listing([listing.isdir]);
for k = 1:numel(listing)
    d = fullfile(listing(k).folder, listing(k).name);
    if ~isfile(fullfile(d, 'psa_object.mat')); continue; end
    psa = ParamSpaceAnalysis2.from_dir(d);
    swept = setdiff(psa.grid_params, {'reps'}, 'stable');
    if ~isempty(swept) && strcmp(swept{1}, param); return; end
end
error('fig_supp4_grouped:NoGainSweep', ...
    'No 1D_sensitivity_* sweep over %s in:\n  %s\nThe sweeps stage of the bundle runs it.', param, run_dir);
end
