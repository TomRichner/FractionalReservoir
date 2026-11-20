function files = fig_outliers(data_dir, fig_dir, opts)
%FIG_OUTLIERS The outlier figures, each with a sidecar .md that carries the numbers and a caption.
%
%   files = fig_outliers()
%   files = fig_outliers(data_dir, fig_dir, 'seed', [7 8])
%
% Reads the runs of explore_outliers (through outlier_table) and writes
% four PNGs, each making one point, with nothing written on the panels but
% axis labels and one short column title:
%
%   outliers_timeseries   rows: networks, columns: conditions. The rate of
%                         the rightmost eigenvalue of the frozen-gain
%                         connectivity against time (condition color), the
%                         bulk edge (gray), the K_eff prediction (thin
%                         black) and zero (dashed). The point: where the
%                         rightmost mode sits between the bulk and zero,
%                         moment by moment.
%   outliers_edge_vs_top  one panel per network: rightmost rate against
%                         bulk-edge rate, one dot per state, all
%                         conditions. On the diagonal there is no outlier;
%                         above it the rightmost mode is outside the bulk;
%                         above zero it is unstable.
%   outliers_dwell        one panel per network: survival curves of the
%                         dwell times beyond the bulk edge (solid) and
%                         beyond zero (dashed).
%   outliers_exposure     the proposed summary: mean exposure X by
%                         condition, one panel per network, one dot per
%                         seed.
%
% The first three use one seed (opts.seed; the first on disk by default);
% the last uses every seed on disk. Figures are made invisible and closed.
% Colors and names are manuscript_style's.
%
% data_dir   default data/explorations/outliers under the repository root
% fig_dir    default figs/explorations/outliers of the MAIN checkout
%            (../FractionalReservoir beside this worktree) when it exists,
%            so the figures are where Tom looks; else this repository's
%
% See also: explore_outliers, outlier_table, outlier_measures

arguments
    data_dir (1,:) char = ''
    fig_dir  (1,:) char = ''
    opts.seed double = []
    opts.resolution (1,1) double = 200
end

root = fileparts(which('setup_paths'));
if isempty(data_dir)
    data_dir = fullfile(root, 'data', 'explorations', 'outliers');
end
if isempty(fig_dir)
    main_checkout = fullfile(fileparts(root), 'FractionalReservoir');
    if isfolder(fullfile(main_checkout, 'figs'))
        fig_dir = fullfile(main_checkout, 'figs', 'explorations', 'outliers');
    else
        fig_dir = fullfile(root, 'figs', 'explorations', 'outliers');
    end
end
if ~isfolder(fig_dir); mkdir(fig_dir); end

[T, M] = outlier_table(data_dir, 'print', false);
st = manuscript_style();
cond_order = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
conds = cond_order(ismember(cond_order, T.condition));
nets = unique(T.network, 'stable');
if any(strcmp(nets, 'reference'))                       % the reference first
    nets = [{'reference'}; setdiff(nets(:), {'reference'}, 'stable')];
end
seed1 = opts.seed;
if isempty(seed1); seed1 = T.seed(1); else; seed1 = seed1(1); end
n_net = numel(nets); n_cond = numel(conds);
pick = @(net, cond, seed) find(strcmp(T.network, net) & strcmp(T.condition, cond) & T.seed == seed, 1);
net_label = containers.Map( ...
    {'reference', 'EEx1p5', 'EEx0p5', 'EI3to2'}, ...
    {'reference', 'E-to-E \times 1.5', 'E-to-E \times 0.5', 'E:I 3:2'});
gray = [0.55 0.55 0.55];
files = {};

%% 1. Time series
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [1 1 12 3.2 * n_net]);
tl = tiledlayout(fig, n_net, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
for w = 1:n_net
    for c = 1:n_cond
        ax = nexttile(tl); hold(ax, 'on');
        k = pick(nets{w}, conds{c}, seed1);
        if ~isempty(k)
            m = M{k};
            col = st.condition_color(conds{c});
            plot(ax, m.t, (m.K1 - 1) ./ m.tau_d, '-', 'Color', [0 0 0 0.35], 'LineWidth', 0.75);
            plot(ax, m.t, m.rate_edge, '-', 'Color', gray, 'LineWidth', 1.5);
            plot(ax, m.t, m.rate_top, '-', 'Color', col, 'LineWidth', 1.0);
            yline(ax, 0, '--', 'Color', 'k', 'LineWidth', 1.0);
            xlim(ax, [m.t(1) m.t(end)]);
            % The K_eff line can sit far from everything else; the view is
            % set by the two series the panel is about.
            v = [m.rate_edge, m.rate_top, 0];
            pad = 0.08 * max(range_(v), 0.5);
            ylim(ax, [min(v) - pad, max(v) + pad]);
        end
        set(ax, 'FontSize', st.tick_fs - 2, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
        if w == 1; title(ax, st.condition_title(conds{c}), 'FontSize', st.title_fs, 'FontWeight', 'normal'); end
        if c == 1; ylabel(ax, sprintf('%s\nrate (s^{-1})', net_label(nets{w})), 'FontSize', st.label_fs - 2); end
        if w == n_net; xlabel(ax, 'time (s)', 'FontSize', st.label_fs - 2); end
    end
end
files{end + 1} = save_fig(fig, fig_dir, 'outliers_timeseries', opts.resolution);
cols = {'lambda_1', 'kappa', 'rate_edge_med', 'rate_top_med', 'rate_top_p95', 'excess_med', ...
    'ratio_med', 'frac_detached', 'frac_unstable', 'K1_med', 'K_err_med', 'K_corr', 'dbar_E', 'dbar_I'};
write_sidecar(fig_dir, 'outliers_timeseries', ...
    ['Rows are the networks (reference, then the unbalanced one), columns the conditions. Colored line: the rate ' ...
    '(Re lambda - 1)/tau_d of the rightmost eigenvalue of the frozen-gain effective connectivity, every 10 ms over ' ...
    '20 s after 60 s of settling. Gray: the bulk edge, the radius predicted from the gains at each state times the ' ...
    'run''s null-calibrated factor kappa. Thin black: the rate the K_eff prediction gives the mean-induced outlier ' ...
    '(off the panel where it is far from both). Dashed: zero. The distance of the colored line above the gray one is ' ...
    'how far the rightmost mode is outside the bulk; above the dashed line it is locally unstable. Seed ' ...
    sprintf('%d', seed1) '. The y axes differ between panels.'], T(T.seed == seed1, :), cols);

%% 2. Rightmost rate against bulk-edge rate
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [1 1 5.2 * n_net 4.6]);
tl = tiledlayout(fig, 1, n_net, 'TileSpacing', 'compact', 'Padding', 'compact');
for w = 1:n_net
    ax = nexttile(tl); hold(ax, 'on');
    lo = Inf; hi = -Inf;
    for c = 1:n_cond
        k = pick(nets{w}, conds{c}, seed1);
        if isempty(k); continue; end
        m = M{k};
        scatter(ax, m.rate_edge, m.rate_top, 6, st.condition_color(conds{c}), 'filled', ...
            'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none');
        lo = min([lo, m.rate_edge, m.rate_top]); hi = max([hi, m.rate_edge, m.rate_top]);
    end
    if isfinite(lo)
        pad = 0.05 * max(hi - lo, 1);
        lims = [lo - pad, hi + pad];
        plot(ax, lims, lims, '-', 'Color', gray, 'LineWidth', 1.0);
        xline(ax, 0, '--', 'Color', 'k', 'LineWidth', 1.0);
        yline(ax, 0, '--', 'Color', 'k', 'LineWidth', 1.0);
        xlim(ax, lims); ylim(ax, lims);
    end
    axis(ax, 'square');
    set(ax, 'FontSize', st.tick_fs - 2, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
    title(ax, net_label(nets{w}), 'FontSize', st.title_fs, 'FontWeight', 'normal');
    xlabel(ax, 'bulk edge rate (s^{-1})', 'FontSize', st.label_fs - 2);
    if w == 1; ylabel(ax, 'rightmost rate (s^{-1})', 'FontSize', st.label_fs - 2); end
end
files{end + 1} = save_fig(fig, fig_dir, 'outliers_edge_vs_top', opts.resolution);
cols = {'rate_edge_med', 'rate_top_med', 'excess_med', 'excess_p95', 'ratio_med', 'frac_detached', ...
    'frac_unstable', 'frac_both', 'n_out_mean', 'n_out_max', 'n_unstable_mean', 'X'};
write_sidecar(fig_dir, 'outliers_edge_vs_top', ...
    ['One dot per sampled state; black, orange and blue are no adaptation, one timescale and multiple timescales ' ...
    '(manuscript_style). Horizontal: the rate of the bulk edge. Vertical: the rate of the rightmost eigenvalue of ' ...
    'the frozen-gain connectivity. On the gray diagonal the rightmost eigenvalue is the edge of the bulk: no ' ...
    'outlier. Above the diagonal it is outside the bulk by the vertical distance. The dashed lines are zero: a ' ...
    'cloud that is on the diagonal and at or left of the vertical dashed line has its outliers reined in and its ' ...
    'bulk at or inside the stability boundary. Seed ' sprintf('%d', seed1) '.'], T(T.seed == seed1, :), cols);

%% 3. Dwell times
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [1 1 5.2 * n_net 4.2]);
tl = tiledlayout(fig, 1, n_net, 'TileSpacing', 'compact', 'Padding', 'compact');
for w = 1:n_net
    ax = nexttile(tl); hold(ax, 'on');
    for c = 1:n_cond
        k = pick(nets{w}, conds{c}, seed1);
        if isempty(k); continue; end
        m = M{k};
        col = st.condition_color(conds{c});
        survival(ax, m.dwell.bulk, m.dt, '-', col);
        survival(ax, m.dwell.zero, m.dt, '--', col);
    end
    set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', st.tick_fs - 2, ...
        'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
    title(ax, net_label(nets{w}), 'FontSize', st.title_fs, 'FontWeight', 'normal');
    xlabel(ax, 'dwell time (s)', 'FontSize', st.label_fs - 2);
    if w == 1; ylabel(ax, 'fraction of excursions longer', 'FontSize', st.label_fs - 2); end
end
files{end + 1} = save_fig(fig, fig_dir, 'outliers_dwell', opts.resolution);
cols = {'n_exc_bulk', 'dwell_bulk_med', 'dwell_bulk_max', 'n_exc_zero', 'dwell_zero_med', 'dwell_zero_max', ...
    'n_exc_x', 'dwell_x_med', 'dwell_x_max', 'frac_same', 'base_overlap', 'identity_life_med', ...
    'turnover_per_s', 'N_eff', 'coh1_med'};
write_sidecar(fig_dir, 'outliers_dwell', ...
    ['Survival curves of excursion lengths: the fraction of excursions that last longer than the time on the ' ...
    'horizontal axis. Solid: runs of consecutive states with the rightmost eigenvalue beyond the bulk edge. Dashed: ' ...
    'runs with it beyond zero. Colors are the conditions. Sampling is every 10 ms, so the curves start there: a ' ...
    'shorter excursion can be missed, and two closer than that are counted as one. A run that touches either end of ' ...
    'the 20 s window is included at its observed length (the table of censored runs is in the data). A condition ' ...
    'that is beyond the edge or beyond zero throughout has a single run of 20 s. Seed ' sprintf('%d', seed1) '.'], ...
    T(T.seed == seed1, :), cols);

%% 4. The summary: exposure by condition
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [1 1 3.6 * n_net 3.8]);
tl = tiledlayout(fig, 1, n_net, 'TileSpacing', 'compact', 'Padding', 'compact');
for w = 1:n_net
    ax = nexttile(tl); hold(ax, 'on');
    for c = 1:n_cond
        sel = strcmp(T.network, nets{w}) & strcmp(T.condition, conds{c});
        X = T.X(sel);
        if isempty(X); continue; end
        col = st.condition_color(conds{c});
        plot(ax, c + linspace(-0.12, 0.12, numel(X))' .* (numel(X) > 1), X, 'o', ...
            'MarkerSize', 9, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', 'LineWidth', 1.0);
        plot(ax, c + [-0.25 0.25], mean(X) .* [1 1], '-', 'Color', col, 'LineWidth', 2.0);
    end
    set(ax, 'XTick', 1:n_cond, 'XTickLabel', cellfun(@(q) st.condition_short(q), conds, 'UniformOutput', false), ...
        'XLim', [0.4, n_cond + 0.6], 'FontSize', st.tick_fs - 2, 'LineWidth', st.axis_lw, ...
        'Box', 'off', 'TickDir', 'out');
    yl = ylim(ax); ylim(ax, [0, max(yl(2), 0.1)]);
    title(ax, net_label(nets{w}), 'FontSize', st.title_fs, 'FontWeight', 'normal');
    if w == 1; ylabel(ax, 'outlier exposure X (s^{-1})', 'FontSize', st.label_fs - 2); end
end
files{end + 1} = save_fig(fig, fig_dir, 'outliers_exposure', opts.resolution);
cols = {'lambda_1', 'X', 'frac_both', 'frac_detached', 'frac_unstable', 'ratio_med', ...
    'full_abscissa_med', 'full_abscissa_max', 'full_frac_unstable', 'full_n_net_mean', 'kin_edge', ...
    'full_minus_fast_med', 'full_states', 'shuffle_pct_med', 'uniform_minus_obs'};
write_sidecar(fig_dir, 'outliers_exposure', ...
    ['The proposed summary measure. Outlier exposure X is the time average of ' ...
    'max(0, Re lambda_1 - max(edge, 1))/tau_d: the part of the rightmost frozen-gain rate that is both outside the ' ...
    'bulk and beyond zero, averaged over the 20 s window (1/s). It is zero when the rightmost eigenvalue is inside ' ...
    'the bulk or the network is locally stable, and grows with how far out and how long the outlier is. One dot per ' ...
    'seed, the bar is their mean. It is blind to outliers that are far out but stable, to outliers away from the ' ...
    'real axis on the left, to what the slow variables do to a mode (compare the full-Jacobian columns below), and ' ...
    'to non-normal growth; and a network whose every gain has collapsed scores zero with nothing reined in, so it ' ...
    'is read beside lambda_1 and the ratio Re lambda_1 / edge.'], T, cols);
end

%% ------------------------------------------------------------------------
function survival(ax, dwell, dt, style, col)
% Empirical survival of the dwell times, as stairs from dt.
if isempty(dwell); return; end
d = sort(dwell(:));
frac = (numel(d):-1:1)' ./ numel(d);                   % fraction with dwell >= d(i)
stairs(ax, [dt; d], [1; frac], style, 'Color', col, 'LineWidth', 2.0);
end

function r = range_(v)
r = max(v) - min(v);
end

function file = save_fig(fig, fig_dir, name, resolution)
file = fullfile(fig_dir, [name '.png']);
exportgraphics(fig, file, 'Resolution', resolution);
close(fig);
fprintf('[fig_outliers] wrote %s\n', file);
end

function write_sidecar(fig_dir, name, caption, T, cols)
% The caption paragraph and the numbers of the figure's runs, as markdown.
fid = fopen(fullfile(fig_dir, [name '.md']), 'w');
closer = onCleanup(@() fclose(fid));
fprintf(fid, '# %s\n\n%s\n\n', name, caption);
fprintf(fid, ['Rates are in 1/s, dwell times in s, eigenvalues of the effective connectivity in units of W. ' ...
    'Definitions: `scripts/explorations/outliers/outlier_measures.m`.\n\n']);
head = [{'network', 'condition', 'seed'}, cols];
fprintf(fid, '| %s |\n', strjoin(head, ' | '));
fprintf(fid, '|%s\n', repmat('---|', 1, numel(head)));
for k = 1:height(T)
    cells = cell(1, numel(cols));
    for q = 1:numel(cols)
        cells{q} = sprintf('%.3g', T.(cols{q})(k));
    end
    fprintf(fid, '| %s | %s | %d | %s |\n', T.network{k}, T.condition{k}, T.seed(k), strjoin(cells, ' | '));
end
end
