function out = fig_supp1_grouped(cfg)
% FIG_SUPP1_GROUPED Supplementary Figure S1: one neuron, and depression with one pool against a reserve; saved data and the preset only.
%
%   out = FIG_SUPP1_GROUPED('run_dir', d, 'preset_name', p, 'out_dir', o)
%
% The figure plan of 2026-11-21 (figures_main_supp_appen.md, S1), as revised
% after Tom's review of 2026-11-22 (Plan2_figure_feedback.md, S1): two
% panels. Each is drawn by the function that owns it, invisibly and unsaved,
% and its axes are copied here (composite_tools).
%   (A) pulse responses of one unconnected neuron with no adaptation, with
%       SFA alone and with depression alone: fig_single_neuron_revised, the
%       rows named by cfg.pulse_rows (by default input, rate and synaptic
%       output; the source also draws the dendritic potential, the SFA
%       feedback and the depression factor), the columns cfg.pulse_conditions.
%       Within a row the columns share ONE y axis that starts at 0, with
%       two ticks (0 and a round value under the row's largest value: 0.5
%       on the rate and synaptic-output rows), labelled on the first column.
%   (B) depression with one pool and with a reserve (the chain) through
%       three steps of the rate: the two step-response panels of
%       fig_STD_steady_state (depression factor, synaptic output), one above
%       the other, with one key that names only what they draw.
%
% WHAT LEFT THE FIGURE (the review): the old panel (B), the SFA current at
% steady state and after a step (fig_SFA_steady_state), and the two
% steady-state panels of the old (C). Both functions still draw them as
% figures of their own.
%
% LAYOUT. Made on the common canvas itself (canvas_w of manuscript_style, in
% pixels), not on a wider sheet squeezed onto it.
%
% (A) reads the single_neuron stage of run_dir (an error of
% fig_single_neuron_revised when the run lacks it). (B) is computed from the
% preset by its own function (the step response of one synapse), as it
% always was; no network is simulated.
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
%           fig_STD_steady_state, fig_SFA_steady_state

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
st = manuscript_style();
quiet = {'save', false, 'visible', false};
preset = {};
if ~isempty(cfg.preset_name); preset = {'preset_name', cfg.preset_name}; end

%% The source figures
S_pulse = fig_single_neuron_revised('run_dir', cfg.run_dir, 'conditions', cfg.pulse_conditions, quiet{:});
S_std   = fig_STD_steady_state(preset{:}, 'run_dir', cfg.run_dir, 'conditions', cfg.conditions, ...
    'verbose', cfg.verbose, quiet{:});
src_figs = [S_pulse.figs(:); S_std.figs(:)];
cleanup = onCleanup(@() T.close(src_figs));

%% (A): which rows of the pulse grid
G = T.grid(findall(S_pulse.figs, 'Type', 'axes'));
row_labels = arrayfun(@(a) strjoin(cellstr(a.YLabel.String), ' '), G(:, 1), 'UniformOutput', false);
[found, rows] = ismember(cfg.pulse_rows, row_labels);
assert(all(found), 'fig_supp1_grouped:NoSuchRow', ...
    'fig_single_neuron_revised draws no row labelled %s. Its rows: %s.', ...
    strjoin(cfg.pulse_rows(~found), ', '), strjoin(row_labels', ', '));
col_titles = arrayfun(@(a) a.Title.String, G(1, :), 'UniformOutput', false);
nr = numel(rows); ncol = size(G, 2);
pick = G(rows, :);

%% The canvas: the common width itself, so that T.fit has nothing to stretch
W = round(st.canvas_w * get(groot, 'ScreenPixelsPerInch'));
X0 = 98; RM = 16;              % left of the first plot box (a y label on two lines and the tick labels); right margin
AGX = 44; AH = 72; AGY = 20;   % (A): between columns (clear of the end labels of the time axes), a row's height, between rows
TOP = 46;                      % the column titles and the letter (A)
BELOW = 58;                    % tick labels and the axis name under a block
BETWEEN = 46;                  % from (A)'s axis name to (B)'s first plot box (its title, the letter)
BH = 165; BGY = 48;            % (B): a plot box's height; between the two (the lower one's title)
KEY_GAP = 14;                  % between (B)'s plot boxes and their key
AW = (W - X0 - RM - (ncol - 1) * AGX) / ncol;
H = TOP + nr * AH + (nr - 1) * AGY + BELOW + BETWEEN + 2 * BH + BGY + BELOW + 6;
fig = T.canvas(W, H, cfg.visible);

A = T.copy(pick(:), fig, false);
ax_a = reshape(A.axes, nr, ncol);          % copyobj keeps the order given
a_top = H - TOP;
for r = 1:nr
    for c = 1:ncol
        ax_a(r, c).Position = [X0 + (c - 1) * (AW + AGX), a_top - r * AH - (r - 1) * AGY, AW, AH];
    end
end

% (B): the two step responses of fig_STD_steady_state, its two lower rows
% (the full-width ones), copied without the source's legends.
src_ax = findall(S_std.figs, 'Type', 'axes');
P = zeros(numel(src_ax), 4);
for k = 1:numel(src_ax); P(k, :) = getpixelposition(src_ax(k), true); end
[~, order] = sortrows([-round(P(:, 2) / 10), P(:, 1)]);
src_ax = src_ax(order);
assert(numel(src_ax) == 4, 'fig_supp1_grouped:DepressionPanels', ...
    'Expected four panels from fig_STD_steady_state (two steady-state, two step responses), found %d.', numel(src_ax));
Bc = T.copy(src_ax(3:4), fig, false);
ax_b = Bc.axes;                             % depression factor, then synaptic output
b_top = a_top - nr * AH - (nr - 1) * AGY - BELOW - BETWEEN;
BW = W - X0 - RM;                           % narrowed below, once the key's width is known
ax_b(1).Position = [X0, b_top - BH, BW, BH];
ax_b(2).Position = [X0, b_top - 2 * BH - BGY, BW, BH];

%% Style
T.style(fig);
% (A): the condition over each column, the quantity at the left of each row
% (the source's own words, on two lines where long), time on the bottom row
% only. One y axis per row, from 0, labelled on the first column.
a_limits = zeros(nr, 2);
for r = 1:nr
    top = 0;
    for c = 1:ncol
        ln = findall(ax_a(r, c), 'Type', 'line');
        for k = 1:numel(ln); top = max(top, max(ln(k).YData(isfinite(ln(k).YData)))); end
    end
    tick = round_below(top);
    a_limits(r, :) = [0, max(1.12 * top, 1.05 * tick)];
    for c = 1:ncol
        ax = ax_a(r, c);
        title(ax, '');
        if r == 1; title(ax, pulse_title(cfg.pulse_conditions, c, col_titles{c}), 'FontWeight', 'normal', 'FontSize', st.title_fs); end
        set(ax, 'YLim', a_limits(r, :), 'YTick', [0, tick]);
        if c > 1; ax.YTickLabel = {}; end
        ax.XTick = ax.XLim(1):10:ax.XLim(2);
        if r < nr
            ax.XTickLabel = {}; xlabel(ax, '');
        else
            ax.XTickLabelMode = 'auto'; ax.XTickLabelRotation = 0;
            xlabel(ax, 'Time (s)');
        end
        if c == 1; ylabel(ax, two_lines(cap_first(cfg.pulse_rows{r}))); end
    end
    linkaxes(ax_a(r, :), 'y');
end
% (B): time on the lower plot only (the two share it); the long y label on two lines.
% Labels with a capital first letter (Tom, 2026-11-23, G1: the source's are
% lowercase); time ticks every 20 s, four on the 65 s axis (G5: they were
% every 10 s, seven).
xlabel(ax_b(1), ''); ax_b(1).XTickLabel = {};
ax_b(2).XTickLabelMode = 'auto';
xlabel(ax_b(2), 'Time (s)');
for k = 1:2; ax_b(k).XTick = 0:20:ax_b(k).XLim(2); end
ylabel(ax_b(1), cap_first(strjoin(cellstr(ax_b(1).YLabel.String), ' ')));
ylabel(ax_b(2), split_at(cap_first(strjoin(cellstr(ax_b(2).YLabel.String), ' ')), '\theta'));
% One key for the two plots, at their right, with what they draw and nothing
% else: the conditions, the reserve, the rate, the steady state. The lines
% are named by fig_STD_steady_state (Tag and DisplayName).
key = gobjects(0);
tags = {'std_step_condition', 'std_step_reserve', 'std_step_rate', 'std_step_steady'};
h_key = gobjects(0, 1); s_key = {};
for k = 1:numel(tags)
    h = findall(ax_b(1), 'Type', 'line', 'Tag', tags{k});
    if isempty(h); continue; end
    if k > 1; h = h(1); end                 % one entry per kind of reference line
    h_key = [h_key; h(:)]; %#ok<AGROW>
    s_key = [s_key, key_strings(get(h, {'DisplayName'})', T)]; %#ok<AGROW>
end
if ~isempty(h_key)
    key = legend(ax_b(1), h_key, s_key, 'Box', 'off', 'FontSize', st.tick_fs, 'Units', 'pixels', 'Location', 'none');
    key.ItemTokenSize = [22 18];
    for pass = 1:2
        drawnow;
        p = key.Position;
        BW = W - RM - p(3) - KEY_GAP - X0;
        ax_b(1).Position(3) = BW; ax_b(2).Position(3) = BW;
        key.Position = [W - RM - p(3), b_top - p(4), p(3), p(4)];
    end
end

%% The common canvas (the layout is already on it; fit leaves it as it is)
T.style(fig);                % again: the labels and titles set since take their roles' sizes
T.fit(fig);

%% Panel letters
T.letter(ax_a(1, 1), 'A', -(X0 - 6), 18);
T.letter(ax_b(1), 'B', -(X0 - 6), 18);

%% Save
panels = struct( ...
    'A', sprintf('fig_single_neuron_revised (rows %s; conditions %s): pulse responses of one neuron', ...
        strjoin(cfg.pulse_rows, ', '), strjoin(cfg.pulse_conditions, ', ')), ...
    'B', 'fig_STD_steady_state (its two step-response panels): depression with one pool and with a reserve, through steps of the rate');
sources = {S_pulse.source, S_std.source};
parts = cell(0, 2);
if cfg.save
    for r = 1:nr
        for c = 1:ncol
            parts(end + 1, :) = {sprintf('Panel A, row %d, column %d (fig_single_neuron_revised)', r, c), ...
                T.axes_table(ax_a(r, c), sprintf('%s, %s; y axis %g to %.3g, shared by the row', ...
                cfg.pulse_rows{r}, cfg.pulse_conditions{c}, a_limits(r, :)))}; %#ok<AGROW>
        end
    end
    what = {'Depression factor through three steps of the rate', 'Synaptic output through the same steps'};
    for k = 1:2
        parts(end + 1, :) = {sprintf('Panel B, plot %d (fig_STD_steady_state)', k), T.axes_table(ax_b(k), what{k})}; %#ok<AGROW>
    end
    if ~isempty(key)
        parts(end + 1, :) = {'Panel B, its key', sprintf('- %s', strjoin(strrep(cellstr(key.String), newline, ' '), [newline '- ']))};
    end
end
caption = sprintf(['Supplementary Figure S1, one neuron (drawn by fig_supp1_grouped; run `%s`, preset `%s`). (A) One ' ...
    'unconnected neuron driven by a pulse of input, with no adaptation, with SFA alone (one timescale) and with ' ...
    'depression alone (one pool): rows %s. SFA lowers the firing rate; depression leaves the rate and lowers the ' ...
    'synaptic output. Within a row the three columns share one y axis that starts at 0. The two single mechanisms ' ...
    'are not conditions of the paper and have colours of their own, named by the column titles: SFA alone purple ' ...
    '[0.56 0.27 0.68], depression alone bluish green [0 0.62 0.45]. (B) Depression with one pool (single-timescale ' ...
    'condition) and with a reserve pool behind it (the chain, multiple-timescale condition) through three steps of ' ...
    'the rate: the depression factor (above) and the synaptic output (below); gray: the rate; dotted: the steady ' ...
    'state of each step; dashed: the reserve. Canvas %d x %d px (canvas_w of manuscript_style), type in points by role.'], ...
    cfg.run_dir, cfg.preset_name, strjoin(cfg.pulse_rows, ', '), getappdata(fig, 'composite_size'));
out_dir = '';
if cfg.save; out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath')); end
out = composite_save(fig, cfg, out_dir, cfg.fig_tag, 'Supplementary Figure S1: one neuron, and depression with one pool against a reserve', ...
    caption, parts, sources);
out.panels = panels;
end

%% ------------------------------------------------------------------------
function t = pulse_title(names, c, fallback)
% The title over a column of (A), which is also the key to its colour: the
% mechanism in words where the column isolates one (the source writes
% 'SFA 1\tau', 'STD 1\tau'), else the source's own title.
t = fallback;
if c > numel(names); return; end
tok = regexp(names{c}, '^sfa(\d+)_std(c?)(\d+)$', 'tokens', 'once');
if strcmp(names{c}, 'no_adaptation')
    t = 'No adaptation';
elseif ~isempty(tok)
    n_sfa = str2double(tok{1}); n_std = str2double(tok{3});
    if n_sfa > 0 && n_std == 0
        t = 'SFA alone';
    elseif n_sfa == 0 && n_std > 0
        t = 'Depression alone';
    end
end
end

function s = cap_first(s)
% A label with its first letter capitalized (sentence case; the rest as it is).
if ~isempty(s); s(1) = upper(s(1)); end
end

function s = two_lines(s)
% A label of more than one word on two lines, broken at its first space.
k = find(s == ' ', 1);
if ~isempty(k) && numel(s) > 8
    s = {s(1:k - 1), s(k + 1:end)};
end
end

function s = split_at(s, mark)
% A label on two lines, broken before mark (as it is when mark is absent).
k = strfind(s, mark);
if ~isempty(k) && k(1) > 1
    s = {strtrim(s(1:k(1) - 1)), strtrim(s(k(1):end))};
end
end

function v = round_below(top)
% The largest of 1, 2, 5 times a power of ten that is not above top: the
% second tick of a row of (A) (0.5 for a row that reaches 0.55, 0.2 for one
% that reaches 0.25).
if ~(top > 0); v = 1; return; end
e = floor(log10(top));
m = [1 2 5] * 10^e;
v = m(find(m <= top * (1 + 1e-9), 1, 'last'));
end

function strings = key_strings(strings, T)
% Entries of (B)'s key: a condition title in sentence case on two lines,
% anything else with its first word capitalized.
for k = 1:numel(strings)
    try
        name = grouped_condition_label(strings(k), 'name_of_title');
        t = T.condition_title(name{1}, true);
        strings{k} = strjoin(cellstr(t), newline);
    catch err
        if ~any(strcmp(err.identifier, {'grouped_condition_label:UnknownTitle', 'grouped_condition_label:AmbiguousTitle'}))
            rethrow(err);
        end
        strings{k} = regexprep(T.sentence(strings{k}), '\s+', ' ');
    end
end
end
