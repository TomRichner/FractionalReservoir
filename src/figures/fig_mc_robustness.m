function out = fig_mc_robustness(cfg)
% FIG_MC_ROBUSTNESS The distribution of memory capacity over the whole gain-by-bias map, per condition (block C).
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
% Conditions by cfg.conditions (select_conditions), colours and titles by
% name from manuscript_style. The table is read by read_mc_gain_sweep.
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
if ~isempty(cfg.data)
    if isstruct(cfg.data) && isfield(cfg.data, 'gains'); D = cfg.data; else; D = read_mc_gain_sweep(cfg.data); end
else
    D = resolve_mc_gain_sweep(cfg.data_file, cfg.run_dir);
end
conds = select_conditions(D.conditions, cfg.conditions);
nc = numel(conds);
st = manuscript_style();

fig = figure('Color', 'w', 'Position', [80 80 max(900, 230 * nc + 520) 420]);
if ~cfg.visible; set(fig, 'Visible', 'off'); end
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'loose', 'Padding', 'compact');
ax_a = nexttile(tl); hold(ax_a, 'on');
ax_b = nexttile(tl); hold(ax_b, 'on');
s = RandStream('mt19937ar', 'Seed', 1);   % jitter only, fixed so the figure is repeatable
n_failed = 0;
h = gobjects(1, nc);
vmax = 0;
for c = 1:nc
    T = D.T(strcmp(D.T.condition, conds{c}), :);
    n_failed = n_failed + nnz(T.failed);
    v = T.mc_total(~T.failed);
    col = cond_color(st, conds{c});
    if isempty(v); continue; end
    vmax = max(vmax, max(v));
    jit = 0.28 * (rand(s, numel(v), 1) - 0.5);
    scatter(ax_a, c + jit, v, 10, col, 'filled', 'MarkerFaceAlpha', 0.35);
    q = prctile(v, [5 25 50 75 95]);
    plot(ax_a, [c c], q([1 2]), '-', 'Color', col, 'LineWidth', 1.2);
    plot(ax_a, [c c], q([4 5]), '-', 'Color', col, 'LineWidth', 1.2);
    rectangle(ax_a, 'Position', [c - 0.22, q(2), 0.44, max(q(4) - q(2), eps)], ...
        'EdgeColor', col, 'LineWidth', 1.5);
    plot(ax_a, c + [-0.22 0.22], q(3) * [1 1], '-', 'Color', col, 'LineWidth', st.median_lw);
    xs = sort(v);
    share = (numel(xs):-1:1)' / numel(xs);   % share with MC >= xs(k), so >= x for x in (xs(k-1), xs(k)]
    h(c) = stairs(ax_b, [0; xs], [share; 0], '-', 'Color', col, 'LineWidth', st.line_lw);
end
xlim(ax_a, [0.4, nc + 0.6]);
set(ax_a, 'XTick', 1:nc, 'XTickLabel', cellfun(@(n) cond_short(st, n), conds, 'UniformOutput', false));
ylabel(ax_a, 'Memory capacity', 'FontSize', st.label_fs - 2);
title(ax_a, sprintf('Every grid point and seed (%d x %d x %d)', numel(D.gains), numel(D.biases), ...
    numel(D.seeds)), 'FontWeight', 'normal', 'FontSize', st.title_fs - 2);
ylim(ax_a, [0, 1.05 * max(vmax, eps)]);
xlim(ax_b, [0, 1.05 * max(vmax, eps)]); ylim(ax_b, [0 1.02]);
xlabel(ax_b, 'Memory capacity x', 'FontSize', st.label_fs - 2);
ylabel(ax_b, 'Share of map with MC \geq x', 'FontSize', st.label_fs - 2);
ok = isgraphics(h);
legend(ax_b, h(ok), cellfun(@(n) cond_title(st, n), conds(ok), 'UniformOutput', false), ...
    'Location', 'northeast', 'Box', 'off', 'FontSize', st.tick_fs - 4);
set([ax_a ax_b], 'FontSize', st.tick_fs - 2, 'TickDir', 'out', 'Box', 'off');
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
function t = cond_title(st, nm)
if isKey(st.condition_title, nm); t = st.condition_title(nm); else; t = nm; end
end

function t = cond_short(st, nm)
if isKey(st.condition_short, nm); t = st.condition_short(nm); else; t = nm; end
end

function c = cond_color(st, nm)
if isKey(st.condition_color, nm); c = st.condition_color(nm); else; c = [0.3 0.3 0.3]; end
end
