function out = fig_supp1_grouped(cfg)
% FIG_SUPP1_GROUPED Supplementary Figure S1: one neuron, and one against several timescales; saved data and the preset only.
%
%   out = FIG_SUPP1_GROUPED('run_dir', d, 'preset_name', p, 'out_dir', o)
%
% The figure plan of 2026-11-21 (figures_main_supp_appen.md, S1). Every
% panel is drawn by the function that owns it, invisibly and unsaved, and
% its axes are copied here (composite_tools).
%   (A) pulse responses of one unconnected neuron with no adaptation, with
%       SFA alone and with depression alone: fig_single_neuron_revised, the
%       rows named by cfg.pulse_rows (by default input, rate and synaptic
%       output; the source also draws the dendritic potential, the SFA
%       feedback and the depression factor), the columns cfg.pulse_conditions.
%   (B) SFA steady state and step response with one and with three
%       timescales: fig_SFA_steady_state.
%   (C) depression with one pool and with a reserve (the chain), steady
%       state and steps: fig_STD_steady_state. The numbers it writes beside
%       its ceiling lines are in the table, not on the panel.
%
% (A) reads the single_neuron stage of run_dir (an error of
% fig_single_neuron_revised when the run lacks it). (B) and (C) are computed
% from the preset by their own functions (analytic steady states and the
% step response of one synapse), as they always were; no network is
% simulated.
%
% Options
%   run_dir, preset_name
%   conditions         the network conditions (fig_STD_steady_state draws
%                      the depressing ones); empty = that function's default
%   pulse_conditions   the single-neuron conditions of (A), left to right
%   pulse_rows         the rows of (A), by the y label of the source's
%                      first column, top to bottom
%   out_dir, save, visible, fig_tag
%
% out: figs, files, source, panels (letter -> source function).
%
% See also: composite_tools, composite_save, fig_single_neuron_revised,
%           fig_SFA_steady_state, fig_STD_steady_state

arguments
    cfg.run_dir          (1,:) char    = ''
    cfg.preset_name      (1,:) char    = ''
    cfg.conditions                     = {}
    cfg.pulse_conditions (1,:) cell    = {'no_adaptation', 'sfa1_std0', 'sfa0_std1'}
    cfg.pulse_rows       (1,:) cell    = {'input u', 'raw rate r', 'synaptic output'}
    cfg.out_dir          (1,:) char    = ''
    cfg.save             (1,1) logical = true
    cfg.visible          (1,1) logical = false
    cfg.verbose          (1,:) char    = 'minimal'
    cfg.fig_tag          (1,:) char    = 'Fig_Supp1_Grouped'
end
setup_paths();
T = composite_tools();
quiet = {'save', false, 'visible', false};
preset = {};
if ~isempty(cfg.preset_name); preset = {'preset_name', cfg.preset_name}; end

%% The source figures
S_pulse = fig_single_neuron_revised('run_dir', cfg.run_dir, 'conditions', cfg.pulse_conditions, quiet{:});
S_sfa   = fig_SFA_steady_state(preset{:}, 'run_dir', cfg.run_dir, 'conditions', cfg.conditions, ...
    'verbose', cfg.verbose, quiet{:});
S_std   = fig_STD_steady_state(preset{:}, 'run_dir', cfg.run_dir, 'conditions', cfg.conditions, ...
    'verbose', cfg.verbose, quiet{:});
src_figs = [S_pulse.figs(:); S_sfa.figs(:); S_std.figs(:)];
cleanup = onCleanup(@() T.close(src_figs));

%% The canvas
W = 1100; H = 955;             % under 0.9 of the workstation screen: exported cropped as usual
fig = T.canvas(W, H, cfg.visible);

% (A): the named rows of the pulse grid.
G = T.grid(findall(S_pulse.figs, 'Type', 'axes'));
row_labels = arrayfun(@(a) strjoin(cellstr(a.YLabel.String), ' '), G(:, 1), 'UniformOutput', false);
[found, rows] = ismember(cfg.pulse_rows, row_labels);
assert(all(found), 'fig_supp1_grouped:NoSuchRow', ...
    'fig_single_neuron_revised draws no row labelled %s. Its rows: %s.', ...
    strjoin(cfg.pulse_rows(~found), ', '), strjoin(row_labels', ', '));
col_titles = arrayfun(@(a) a.Title.String, G(1, :), 'UniformOutput', false);
nr = numel(rows); ncol = size(G, 2);
pick = G(rows, :);
A = T.copy(pick(:), fig, false);
ax_a = reshape(A.axes, nr, ncol);          % copyobj keeps the order given
% An even grid: the rows left out leave no gap.
AX0 = 90; AW = 125; AGX = 65; AY0 = 660; AH = 68; AGY = 24;
for r = 1:nr
    for c = 1:ncol
        ax_a(r, c).Position = [AX0 + (c - 1) * (AW + AGX), AY0 + (nr - r) * (AH + AGY), AW, AH];
    end
end

% (B): the two SFA panels and their legend.
B = T.copy(S_sfa.figs, fig);
[~, order] = sort(arrayfun(@(a) a.Position(1), B.axes)); ax_b = B.axes(order);
assert(numel(ax_b) == 2, 'fig_supp1_grouped:SfaPanels', ...
    'Expected two panels from fig_SFA_steady_state, found %d.', numel(ax_b));
ax_b(1).Position = [740 730 135 165];
ax_b(2).Position = [935 730 135 165];

% (C): steady state (two panels, upper row) and the two step responses.
Cc = T.copy(S_std.figs, fig);
P = vertcat(Cc.axes.Position);
[~, order] = sortrows([-round(P(:, 2) / 10), P(:, 1)]);
ax_c = Cc.axes(order);
assert(numel(ax_c) == 4, 'fig_supp1_grouped:DepressionPanels', ...
    'Expected four panels from fig_STD_steady_state (two steady-state, two step responses), found %d.', numel(ax_c));
ax_c(1).Position = [90 365 290 180];
ax_c(2).Position = [480 365 290 180];
ax_c(3).Position = [90 70 425 180];
ax_c(4).Position = [640 70 425 180];

%% Style
T.style(fig);
% (A): the condition over each column, the quantity at the left of each row
% (the source's own words, on two lines where long), two y ticks per panel
% (each has its own range), time on the bottom row only.
for r = 1:nr
    for c = 1:ncol
        ax = ax_a(r, c);
        title(ax, '');
        if r == 1; title(ax, col_titles{c}, 'FontWeight', 'normal', 'FontSize', 14); end
        ax.YTick = ax.YTick([1 end]);
        ax.XTick = ax.XLim(1):10:ax.XLim(2);
        if r < nr
            ax.XTickLabel = {}; xlabel(ax, '');
        else
            ax.XTickLabelMode = 'auto'; ax.XTickLabelRotation = 0;
            xlabel(ax, 'time (s)');
        end
        if c == 1; ylabel(ax, two_lines(cfg.pulse_rows{r})); end
    end
end
% (B): the two panels show one quantity, named once; short titles tell the
% steady state from the step (the source's are too long for these panels).
title(ax_b(1), 'Steady state', 'FontWeight', 'normal');
title(ax_b(2), 'Step of the rate', 'FontWeight', 'normal');
ylabel(ax_b(2), '');
% (C): the values written beside the ceiling lines go to the table.
ceiling_text = {};
tx = findall(ax_c(2), 'Type', 'text');
for k = 1:numel(tx)
    s = strtrim(strjoin(cellstr(tx(k).String), ' '));
    if startsWith(s, 'ceiling')
        ceiling_text{end + 1} = sprintf('- %s (the dotted line in the colour [%.2f %.2f %.2f])', s, tx(k).Color); %#ok<AGROW>
        delete(tx(k));
    end
end
% Legends: the conditions' beside the steady-state panels with the keys of
% the reference lines under it; the SFA key under its two panels.
lgs = Cc.legends;
for k = 1:numel(lgs)
    lgs(k).String = sentence_conditions(lgs(k).String, T);
end
set([B.legends; lgs], 'Units', 'pixels', 'Location', 'none');
set(B.legends, 'NumColumns', 2, 'Orientation', 'horizontal');
for k = 1:numel(lgs); lgs(k).ItemTokenSize = [22 18]; end
if ~isempty(B.legends); B.legends.ItemTokenSize = [22 18]; end
for pass = 1:2
    drawnow;
    y = 545;
    for k = 1:numel(ax_c)
        own = lgs(arrayfun(@(l) isequal(l.Axes, ax_c(k)), lgs));
        for j = 1:numel(own)
            p = own(j).Position;
            y = y - p(4) - 6;
            own(j).Position = [790, y, p(3), p(4)];
        end
    end
    if ~isempty(B.legends)
        p = B.legends.Position;
        B.legends.Position = [W - 8 - p(3), 628, p(3), p(4)];
    end
end

%% Panel letters
T.letter(ax_a(1, 1), 'A', -75, 22);
T.letter(ax_b(1), 'B', -65, 28);
T.letter(ax_c(1), 'C', -75, 30);

%% Save
panels = struct( ...
    'A', sprintf('fig_single_neuron_revised (rows %s; conditions %s): pulse responses of one neuron', ...
        strjoin(cfg.pulse_rows, ', '), strjoin(cfg.pulse_conditions, ', ')), ...
    'B', 'fig_SFA_steady_state: SFA steady state and step response, one and three timescales', ...
    'C', 'fig_STD_steady_state: depression with one pool and with a reserve, steady state and step responses');
sources = {S_pulse.source, S_sfa.source, S_std.source};
parts = cell(0, 2);
if cfg.save
    for r = 1:nr
        for c = 1:ncol
            parts(end + 1, :) = {sprintf('Panel A, row %d, column %d (fig_single_neuron_revised)', r, c), ...
                T.axes_table(ax_a(r, c), sprintf('%s, %s', cfg.pulse_rows{r}, cfg.pulse_conditions{c}))}; %#ok<AGROW>
        end
    end
    parts(end + 1, :) = {'Panel B, left (fig_SFA_steady_state)', T.axes_table(ax_b(1), 'Steady-state SFA current against the rate')};
    parts(end + 1, :) = {'Panel B, right (fig_SFA_steady_state)', T.axes_table(ax_b(2), 'SFA current after a step of the rate to 1')};
    what = {'Steady-state depression factor against the rate', 'Steady-state synaptic output against the rate', ...
        'Depression factor through three steps of the rate', 'Synaptic output through the same steps'};
    for k = 1:4
        parts(end + 1, :) = {sprintf('Panel C, panel %d (fig_STD_steady_state)', k), T.axes_table(ax_c(k), what{k})}; %#ok<AGROW>
    end
    parts(end + 1, :) = {'Panel C, the ceilings of the synaptic output (written on the source panel, not drawn here)', ...
        strjoin([{'The values:', ''}, ceiling_text], newline)};
end
caption = sprintf(['Supplementary Figure S1, one neuron (drawn by fig_supp1_grouped; run `%s`, preset `%s`). (A) One ' ...
    'unconnected neuron driven by a pulse of input, with no adaptation, with SFA alone (one timescale) and with ' ...
    'depression alone (one pool): rows %s. SFA lowers the firing rate; depression leaves the rate and lowers the ' ...
    'synaptic output. Each panel has its own y range. (B) SFA with three timescales and with one: the steady-state ' ...
    'SFA current against the rate (the two lines coincide), and its approach after a step of the rate. (C) ' ...
    'Depression with one pool (single-timescale condition) and with a reserve pool behind it (the chain, ' ...
    'multiple-timescale condition): steady-state depression factor and synaptic output against the rate (dotted: ' ...
    'the ceiling of the output; gray: the undepressed output), and both through three steps of the rate (gray: ' ...
    'the rate; dotted: the steady state; dashed: the reserve). Canvas %d x %d px, text 14 pt.'], ...
    cfg.run_dir, cfg.preset_name, strjoin(cfg.pulse_rows, ', '), W, H);
out_dir = '';
if cfg.save; out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath')); end
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Supplementary Figure S1: one neuron, and one against several timescales', ...
    caption, parts, sources);
out.panels = panels;
end

%% ------------------------------------------------------------------------
function s = two_lines(s)
% A label of more than one word on two lines, broken at its first space.
k = find(s == ' ', 1);
if ~isempty(k) && numel(s) > 8
    s = {s(1:k - 1), s(k + 1:end)};
end
end

function strings = sentence_conditions(strings, T)
% Legend entries that are condition titles, in sentence case; the others as they are.
for k = 1:numel(strings)
    try
        name = grouped_condition_label(strings(k), 'name_of_title');
        strings{k} = T.condition_title(name{1});
    catch err
        if ~any(strcmp(err.identifier, {'grouped_condition_label:UnknownTitle', 'grouped_condition_label:AmbiguousTitle'}))
            rethrow(err);
        end
    end
end
end
