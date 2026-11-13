function out = fig_mc_robustness(cfg)
% FIG_MC_ROBUSTNESS The distribution of memory capacity over the whole map of gain by the sweep's second axis, per condition (block C).
%
%   out = FIG_MC_ROBUSTNESS('run_dir', d)
%   out = FIG_MC_ROBUSTNESS('data_file', f)
%   out = FIG_MC_ROBUSTNESS('data', D, 'save', false)
%
% Tom's question: is chained MTS the most robust, that is, high memory
% capacity over more of the map than 1TS or no adaptation? Two panels:
%   (a) per condition, every job of the map (every grid point and seed) as a
%       jittered point, a box from the 25th to the 75th percentile, the
%       median as a thick bar, whiskers from the 5th to the 95th
%       percentile (the quantities mc_gain_bias_statistics reports);
%   (b) the share of the map with memory capacity at least x, against x,
%       conditions overlaid (one minus the empirical cumulative
%       distribution, over the same jobs).
% Failed jobs are left out; their count is in the subtitle.
%
% Conditions by cfg.conditions (select_conditions), titles and colours by
% name (titles_by_name, mc_condition_colors). The result is taken through
% resolve_mc_gain_sweep (data, data_file or run_dir).
%
% See also: mc_gain_bias_statistics, fig_mc_gain_bias_maps, fig_mc_gain_curves

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.run_dir     (1,:) char    = ''
    cfg.data_file   (1,:) char    = ''
    cfg.data                      = []
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.preset_name (1,:) char    = ''   % unused; accepted for a uniform call
    cfg.conditions                = {}
    cfg.fig_tag     (1,:) char    = 'Fig_MC_robustness'
end

setup_paths();
D = resolve_mc_gain_sweep('data', cfg.data, 'data_file', cfg.data_file, 'run_dir', cfg.run_dir);
conds = select_conditions(D.conditions, cfg.conditions);
nc = numel(conds);
st = manuscript_style();
cond_titles = titles_by_name(conds);
cond_colors = mc_condition_colors(conds);

BOX_HALF_WIDTH = 0.22;      % of the unit spacing between conditions
JITTER_WIDTH   = 0.28;      % the points spread a little wider than the box
QUANTILES      = [5 25 50 75 95];   % whisker, box, median, box, whisker

% Panel (a) is given 3/5 of the width and its tick labels are held
% horizontal: at four conditions MATLAB otherwise rotates the multi-line
% labels by 45 degrees.
fig = figure('Color', 'w', 'Position', [80 80 max(1000, 260 * nc + 560) 440]);
if ~cfg.visible; set(fig, 'Visible', 'off'); end
tl = tiledlayout(fig, 1, 5, 'TileSpacing', 'loose', 'Padding', 'compact');
ax_a = nexttile(tl, [1 3]); hold(ax_a, 'on');
ax_b = nexttile(tl, [1 2]); hold(ax_b, 'on');
s = RandStream('mt19937ar', 'Seed', 1);   % jitter only, fixed so the figure is repeatable
n_failed = 0;
h = gobjects(1, nc);
vmax = 0;
for c = 1:nc
    T = D.T(strcmp(D.T.condition, conds{c}), :);
    n_failed = n_failed + nnz(T.failed);
    v = T.mc_total(~T.failed);
    col = cond_colors(c, :);
    if isempty(v); continue; end
    vmax = max(vmax, max(v));
    jit = JITTER_WIDTH * (rand(s, numel(v), 1) - 0.5);
    scatter(ax_a, c + jit, v, 10, col, 'filled', 'MarkerFaceAlpha', 0.35);
    q = prctile(v, QUANTILES);
    plot(ax_a, [c c], q([1 2]), '-', 'Color', col, 'LineWidth', 1.2);
    plot(ax_a, [c c], q([4 5]), '-', 'Color', col, 'LineWidth', 1.2);
    rectangle(ax_a, 'Position', [c - BOX_HALF_WIDTH, q(2), 2 * BOX_HALF_WIDTH, max(q(4) - q(2), eps)], ...
        'EdgeColor', col, 'LineWidth', 1.5);
    plot(ax_a, c + BOX_HALF_WIDTH * [-1 1], q(3) * [1 1], '-', 'Color', col, 'LineWidth', st.median_lw);
    xs = sort(v);
    share = (numel(xs):-1:1)' / numel(xs);   % share with MC >= xs(k), so >= x for x in (xs(k-1), xs(k)]
    h(c) = stairs(ax_b, [0; xs], [share; 0], '-', 'Color', col, 'LineWidth', st.line_lw);
end
xlim(ax_a, [0.4, nc + 0.6]);
% Tick labels are the conditions' full titles, one word per line where it
% helps: the short names ('Baseline', 'SFA+STD') do not tell the chained
% MTS from the product one, nor match the other block C figures.
set(ax_a, 'XTick', 1:nc, 'XTickLabel', cellfun(@tick_label, cond_titles, 'UniformOutput', false), ...
    'TickLabelInterpreter', 'tex');
ylabel(ax_a, 'Memory capacity', 'FontSize', st.label_fs - 2);
Q2 = mc_sweep_quantity(D.second_axis.column, D);   % the second axis's wording
title(ax_a, sprintf('Every job: %d gains x %d %s x %d seed(s)', numel(D.gains), numel(D.second_levels), ...
    Q2.plural, numel(D.seeds)), 'FontWeight', 'normal', 'FontSize', st.title_fs - 2);
ylim(ax_a, [0, 1.05 * max(vmax, eps)]);
xlim(ax_b, [0, 1.05 * max(vmax, eps)]); ylim(ax_b, [0 1.02]);
xlabel(ax_b, 'Memory capacity x', 'FontSize', st.label_fs - 2);
ylabel(ax_b, 'Share of map with MC \geq x', 'FontSize', st.label_fs - 2);
ok = isgraphics(h);
legend(ax_b, h(ok), cond_titles(ok), ...
    'Location', 'northeast', 'Box', 'off', 'FontSize', st.tick_fs - 4);
set([ax_a ax_b], 'FontSize', st.tick_fs - 2, 'TickDir', 'out', 'Box', 'off');
set(ax_a, 'XTickLabelRotation', 0);
if n_failed > 0
    subtitle(ax_a, sprintf('%d failed job(s) left out', n_failed), 'FontSize', st.tick_fs - 5);
end

out = struct('figs', fig, 'files', {{}}, 'source', D.source);
if cfg.save
    out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
    save_figure_stable(out_dir, cfg.fig_tag, fig);
    out.files = existing_outputs(out_dir, cfg.fig_tag);
end
end

%% ------------------------------------------------------------------------
function t = tick_label(t)
% A condition's title broken over lines, so the tick labels of four
% conditions fit side by side.
t = regexprep(t, '-Timescale ', '-Timescale\\newline');
t = regexprep(t, '\s+(\(.*\))$', '\\newline$1');
t = strrep(t, 'No Adaptation', 'No\newlineAdaptation');
end
