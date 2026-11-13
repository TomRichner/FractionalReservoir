function out = fig_mc_gain_curves(cfg)
% FIG_MC_GAIN_CURVES Memory capacity and lambda_1 against gain at one level of the sweep's second axis (the bias by default), per condition (block C).
%
%   out = FIG_MC_GAIN_CURVES('run_dir', d)
%   out = FIG_MC_GAIN_CURVES('data_file', f)
%   out = FIG_MC_GAIN_CURVES('data_file', {f_std, f_3to2, f_2to3}, ...
%                            'row_labels', {'1 : 1', 'E : I = 3 : 2', 'E : I = 2 : 3'})
%   out = FIG_MC_GAIN_CURVES('data', {D1, D2}, 'save', false)
%   out = FIG_MC_GAIN_CURVES('run_dir', d, 'variants', {'EI3to2', 'EI2to3'})
%
% One ROW per result (the standard network, then the variant networks:
% E : I = 3 : 2 and 2 : 3 in block C), three columns:
%   (1) memory capacity against gain (log axis),
%   (2) lambda_1 with the memory-task input on against gain, zero line,
%   (3) memory capacity against lambda_1.
% Conditions overlaid in their manuscript colours: every seed a point, the
% median over seeds a line, the interquartile range a band (columns 1 and 2;
% in column 3 the line joins the per-gain medians of both quantities, in
% gain order).
%
% lambda_1 IS ON A SYMMETRIC-LOG AXIS (symlog_lle: linear within about
% +-cfg.lle_lin = 0.05 1/s, logarithmic beyond, clipped at +-cfg.lle_clip
% = 10 1/s), the y axis of column 2 and the x axis of column 3. Without
% adaptation lambda_1 runs from about -10 to +8 along gain; with chained
% MTS it sits at about -0.10 (the slowest adaptation mode) at every stable
% gain, and on a linear axis that cluster collapses onto the zero line. On
% the symmetric-log axis both are readable, and in column 3 a cluster (at
% least three per-gain medians with lambda_1 < 0 within 0.03 1/s of each
% other) is also marked by a bracket and its value, in the condition's
% colour, so the vertical stack of points reads as "lambda_1 does not
% move here" rather than as a measurement of distance to instability.
%
% THE LEVEL OF THE SECOND AXIS: cfg.second (default each result's own
% preset value on its second axis: bias 0.2, or self-feedback 0). Each row
% uses the nearest level its result holds, and the row title states both
% in the axis's own wording (mc_sweep_quantity), so a result swept only
% along gain (one level) is drawn at that level. cfg.bias is the same
% option under the name it had when the second axis was always the bias.
%
% WHICH RESULTS, one row each (resolve_mc_gain_sweep):
%   cfg.data       results in memory (one, or a cell of them)
%   cfg.data_file  results files (one, or a cell of them)
%   cfg.run_dir    a run directory: its main result, then the result of
%                  each variant network named in cfg.variants (the tags a
%                  bundle's mc_gain_sweep stage ran, in the order given). A
%                  variant the run directory does not hold is an error.
% ROW LABELS: cfg.row_labels, else each network's ratio of mean excitatory
% to mean inhibitory weight read from its preset ("E : I = 3 : 2"), else
% the file name.
%
% Rows share y limits per column, so the networks can be compared.
% Conditions by cfg.conditions (select_conditions, against each result).
%
% See also: fig_mc_gain_bias_maps, fig_mc_robustness, read_mc_gain_sweep

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.run_dir     (1,:) char    = ''
    cfg.data_file                 = ''      % char, or a cell of files (one row each)
    cfg.data                      = []      % read_mc_gain_sweep output / table, or a cell of them
    cfg.variants    (1,:) cell    = {}      % with run_dir: variant networks to add as rows
    cfg.row_labels                = {}      % one per row; empty = from each result's preset, or its file
    cfg.second                    = []      % level of the second axis; empty = each result's preset value
    cfg.bias                      = []      % cfg.second, by its earlier name
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.preset_name (1,:) char    = ''      % unused; accepted for a uniform call
    cfg.conditions                = {}
    cfg.lle_lin     (1,1) double  = 0.05    % symlog_lle's linear scale (1/s)
    cfg.lle_clip    (1,1) double  = 10      % |lambda_1| clipped here (1/s)
    cfg.fig_tag     (1,:) char    = 'Fig_MC_gain_curves'
end

setup_paths();
if ~isempty(cfg.second) && ~isempty(cfg.bias)
    error('fig_mc_gain_curves:TwoLevels', 'Give second or bias (its earlier name), not both.');
end
second_wanted = [cfg.second, cfg.bias];
Ds = load_all(cfg);
n_rows = numel(Ds);
st = manuscript_style();
labels = cfg.row_labels;
if isempty(labels); labels = cellfun(@row_label, Ds, 'UniformOutput', false); end
if ischar(labels); labels = {labels}; end
assert(numel(labels) == n_rows, 'fig_mc_gain_curves:RowLabels', ...
    '%d row labels for %d results.', numel(labels), n_rows);

sl = @(v) symlog_lle(v, cfg.lle_lin, cfg.lle_clip);
[~, l_tp, l_tl, l_lim] = symlog_lle(0, cfg.lle_lin, cfg.lle_clip);

fig = figure('Color', 'w', 'Position', [80 60 1150 300 * n_rows + 90]);
if ~cfg.visible; set(fig, 'Visible', 'off'); end
tl = tiledlayout(fig, n_rows, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(n_rows, 3);
legend_h = gobjects(0); legend_t = {};
for r = 1:n_rows
    D = Ds{r};
    conds = select_conditions(D.conditions, cfg.conditions);
    cond_colors = mc_condition_colors(conds);
    Q2 = mc_sweep_quantity(D.second_axis.column, D);
    b_want = second_wanted; if isempty(b_want); b_want = D.preset_second; end
    [~, ib] = min(abs(D.second_levels - b_want));
    b = D.second_levels(ib);
    g = D.gains;
    for k = 1:3; ax(r, k) = nexttile(tl); hold(ax(r, k), 'on'); end
    for c = 1:numel(conds)
        col = cond_colors(c, :);
        T = D.T(strcmp(D.T.condition, conds{c}) & Q2.values == b & ~D.T.failed, :);
        [mc_m, mc_lo, mc_hi] = per_gain(T, g, 'mc_total');
        [l_m, l_lo, l_hi]    = per_gain(T, g, 'lle_input');
        band(ax(r, 1), g, mc_lo, mc_hi, col);
        band(ax(r, 2), g, sl(l_lo), sl(l_hi), col);
        scatter(ax(r, 1), T.gain, T.mc_total, 14, col, 'filled', 'MarkerFaceAlpha', 0.45);
        scatter(ax(r, 2), T.gain, sl(T.lle_input), 14, col, 'filled', 'MarkerFaceAlpha', 0.45);
        scatter(ax(r, 3), sl(T.lle_input), T.mc_total, 14, col, 'filled', 'MarkerFaceAlpha', 0.45);
        h = plot(ax(r, 1), g, mc_m, '-', 'Color', col, 'LineWidth', st.line_lw);
        plot(ax(r, 2), g, sl(l_m), '-', 'Color', col, 'LineWidth', st.line_lw);
        plot(ax(r, 3), sl(l_m), mc_m, '-', 'Color', col, 'LineWidth', st.line_lw);
        mark_cluster(ax(r, 3), l_m, mc_m, col, c, sl);
        if ~ismember(conds{c}, legend_t)
            legend_h(end+1) = h; legend_t{end+1} = conds{c}; %#ok<AGROW>
        end
    end
    yline(ax(r, 2), 0, '--', 'Color', st.zeroline_color, 'LineWidth', st.zeroline_lw - 0.5);
    xline(ax(r, 3), 0, '--', 'Color', st.zeroline_color, 'LineWidth', st.zeroline_lw - 0.5);
    [xt, xtl] = gain_ticks(g);
    set(ax(r, 2), 'YTick', l_tp, 'YTickLabel', l_tl);
    set(ax(r, 3), 'XTick', l_tp, 'XTickLabel', l_tl);
    ylim(ax(r, 2), l_lim); xlim(ax(r, 3), l_lim);
    for k = 1:2
        set(ax(r, k), 'XScale', 'log', 'XTick', xt, 'XTickLabel', xtl);
        xlim(ax(r, k), [min(g) max(g)] .* [0.9 1.1]);
        xline(ax(r, k), D.preset_gain, ':', 'Color', [0.5 0.5 0.5]);
    end
    ylabel(ax(r, 1), 'Memory capacity', 'FontSize', st.label_fs - 3);
    ylabel(ax(r, 2), '\lambda_1 (1/s), symmetric log', 'FontSize', st.label_fs - 3);
    ylabel(ax(r, 3), 'Memory capacity', 'FontSize', st.label_fs - 3);
    level_txt = Q2.at(b);
    if abs(b - b_want) > 1e-9; level_txt = sprintf('%s (nearest to %.3g)', level_txt, b_want); end
    title(ax(r, 1), sprintf('%s, %s', labels{r}, level_txt), 'FontSize', st.title_fs - 3, ...
        'FontWeight', 'normal', 'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [0 1.02 0]);
    if r == n_rows
        xlabel(ax(r, 1), 'Gain', 'FontSize', st.label_fs - 3);
        xlabel(ax(r, 2), 'Gain', 'FontSize', st.label_fs - 3);
        xlabel(ax(r, 3), '\lambda_1 (1/s), task input on, symmetric log', 'FontSize', st.label_fs - 3);
    end
end
set(ax, 'FontSize', st.tick_fs - 3, 'TickDir', 'out', 'Box', 'off');
for k = 1:3; linkaxes(ax(:, k), 'y'); end
linkaxes(ax(:, 3), 'x');
lg = legend(ax(1, 3), legend_h, titles_by_name(legend_t), ...
    'Box', 'off', 'FontSize', st.tick_fs - 4);
lg.Layout.Tile = 'south';
lg.Orientation = 'horizontal';

out = struct('figs', fig, 'files', {{}}, 'source', {cellfun(@(d) d.source, Ds, 'UniformOutput', false)});
if cfg.save
    out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
    save_figure_stable(out_dir, cfg.fig_tag, fig);
    out.files = existing_outputs(out_dir, cfg.fig_tag);
end
end

%% ------------------------------------------------------------------------
function Ds = load_all(cfg)
% One result per row; the header's WHICH RESULTS says from where.
if ~isempty(cfg.data)
    Ds = cellfun(@(d) resolve_mc_gain_sweep('data', d), as_cell(cfg.data), 'UniformOutput', false);
elseif ~isempty(cfg.data_file)
    Ds = cellfun(@(f) resolve_mc_gain_sweep('data_file', char(f)), as_cell(cfg.data_file), ...
        'UniformOutput', false);
else
    Ds = [{resolve_mc_gain_sweep('run_dir', cfg.run_dir)}, ...
        cellfun(@(v) resolve_mc_gain_sweep('run_dir', cfg.run_dir, 'variant', v), cfg.variants, ...
        'UniformOutput', false)];
end
end

function c = as_cell(x)
if iscell(x); c = x; else; c = {x}; end
end

function s = row_label(D)
% The network's E : I ratio when its preset gives one, else the file name.
s = '';
if ~isempty(D.preset_name)
    s = ei_ratio_label(D.preset_name);
end
if isempty(s)
    [~, s] = fileparts(D.source);
end
end

function s = ei_ratio_label(preset_name)
% "E : I = p : q", the ratio of the mean excitatory to the mean inhibitory
% weight (|mu_E| : |mu_I|, as a ratio of small integers), read from the
% preset. Empty when the preset is not a two-type E, I network. A
% model_args override of the weights is not seen here; pass row_labels.
s = '';
d = srnn_param_preset(preset_name);
if ~(isfield(d, 'cell_type_names') && isequal(d.cell_type_names, {'E', 'I'})); return; end
mu = d.mu_tilde_relative;                 % (post, pre), or a presynaptic row
mu_by_pre = mean(mu, 1);
[p, q] = rat(abs(mu_by_pre(1) / mu_by_pre(2)));
s = sprintf('E : I = %d : %d', p, q);
end

function [m, lo, hi] = per_gain(T, g, field)
m = nan(size(g)); lo = m; hi = m;
for k = 1:numel(g)
    v = T.(field)(T.gain == g(k));
    v = v(isfinite(v));
    if isempty(v); continue; end
    m(k) = median(v); lo(k) = prctile(v, 25); hi(k) = prctile(v, 75);
end
end

function band(ax, x, lo, hi, col)
ok = isfinite(lo) & isfinite(hi);
if nnz(ok) < 2; return; end
x = x(ok); lo = lo(ok); hi = hi(ok);
fill(ax, [x, fliplr(x)], [hi, fliplr(lo)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
end

function mark_cluster(ax, l_m, mc_m, col, c, sl)
% A run of stable per-gain medians whose lambda_1 barely moves: bracket it
% and state the value, so the stack is not mistaken for a trend.
% The statement goes in a list in the axes' upper-right corner (positive
% lambda_1 with high memory capacity, where no curve runs: once lambda_1 is
% clearly positive the memory capacity has collapsed), one line per
% clustered condition in its colour, stacked downward. Not beside the
% bracket: two conditions that cluster at the same lambda_1 would overlap
% there. Not the lower-left corner, where it sat first: on the
% excitation-dominated network (E : I = 3 : 2) the no-adaptation curve runs
% there (lambda_1 near -10, memory capacity 0 to 3) and the list crossed it.
MIN_POINTS = 3;         % per-gain medians that make a cluster
HALF_WIDTH = 0.03;      % 1/s around their median lambda_1
k = find(isfinite(l_m) & isfinite(mc_m) & l_m < 0);
if numel(k) < MIN_POINTS; return; end
lv = l_m(k);
in = abs(lv - median(lv)) <= HALF_WIDTH;
if nnz(in) < MIN_POINTS; return; end
lc = median(lv(in)); y = mc_m(k(in));
x0 = sl(lc) - 0.15 - 0.08 * (c - 1);
plot(ax, x0 * [1 1], [min(y) max(y)], '-', 'Color', col, 'LineWidth', 1, 'HandleVisibility', 'off');
plot(ax, x0 + [0 0.05], min(y) * [1 1], '-', 'Color', col, 'LineWidth', 1, 'HandleVisibility', 'off');
plot(ax, x0 + [0 0.05], max(y) * [1 1], '-', 'Color', col, 'LineWidth', 1, 'HandleVisibility', 'off');
n_prev = numel(findobj(ax, 'Tag', 'lle_cluster'));
text(ax, 0.98, 0.98 - 0.07 * n_prev, sprintf('[ \\lambda_1 \\approx %.2f at %d gains', lc, nnz(in)), ...
    'Units', 'normalized', 'Color', col, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 8, 'Tag', 'lle_cluster');
end

function [xt, xtl] = gain_ticks(gains)
cand = [0.01 0.03 0.1 0.3 1 3 10 30 100];
xt = cand(cand >= min(gains) * 0.9 & cand <= max(gains) * 1.1);
if numel(xt) < 2; xt = gains([1 end]); end
xtl = arrayfun(@(v) sprintf('%g', v), xt, 'UniformOutput', false);
end

