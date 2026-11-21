function files = fig_outliers(data_dir, fig_dir, opts)
%FIG_OUTLIERS The outlier figures, each with a sidecar .md that carries the numbers and a caption.
%
%   files = fig_outliers()
%   files = fig_outliers(data_dir, fig_dir, 'seed', 7, 'dose_dir', '')
%
% Reads the runs of explore_outliers (through outlier_table) and writes six
% PNGs, each making one point, with nothing on the panels but axis labels
% and one short title per column:
%
%   outliers_complex_plane   rows: networks, columns: conditions. The whole
%                         frozen-gain spectrum at the calibration states in
%                         the plane of rates, the null-calibrated edge of
%                         the bulk (gray circle), the rightmost eigenvalue
%                         of each state (ringed; with its conjugate when it
%                         is complex), where the structure alone
%                         would put it (open diamond: uniform gains), zero
%                         (dashed). The point: is there an eigenvalue
%                         outside the bulk, and where is it relative to
%                         where the weights alone put it.
%   outliers_summary      top: detachment rho = Re lambda_1 / edge by
%                         condition (dot: median, bar: 5th to 95th
%                         percentile over time; one per seed), with the
%                         uniform-gain value (open diamond) and the edge
%                         (line at 1). Bottom: the rightmost rate with the
%                         slow variables frozen (open) and free (filled).
%   outliers_timeseries   the frozen-gain rightmost rate (condition color),
%                         the bulk edge (gray), the full Jacobian's
%                         rightmost real part (black) and zero (dashed)
%                         against time.
%   outliers_slow_variables   what the slow variables do to the rightmost
%                         frozen-gain mode: the full Jacobian's network
%                         eigenvalues at the dense states (dots), the
%                         frozen rightmost eigenvalue (open, with its
%                         conjugate when it is complex) and the full-
%                         Jacobian eigenvalue that carries its pattern
%                         (filled, with its conjugate).
%   outliers_dwell        survival curves of the time spent beyond zero:
%                         frozen-gain rightmost (dashed) and full Jacobian
%                         (solid).
%   outliers_dose         (when dose_dir has runs) detachment, responsive
%                         fraction of the E population and the E-to-I gain
%                         ratio against the structural detachment, across
%                         multiples of the E-to-E mean.
%
% WHAT CHANGED AFTER THE FIRST RUNS (2026-11-19). The exposure figure was
% dropped: exposure is largest exactly where the slow variables hold the
% outlier, and zero where a population is pinned (see the summary sidecar).
% The K_eff line left the time series (it is in the sidecar); the full
% Jacobian's rightmost real part took its place. The beyond-the-edge dwell
% curves were dropped: they are either chance-level blips of 10 to 40 ms
% (no outlier) or one run the length of the window (a standing outlier).
%
% The single-seed figures use opts.seed (the first on disk by default); the
% summary uses every seed. Figures are left OPEN (Tom looks at them in
% MATLAB), tagged 'fig_outliers'; a new call deletes only figures with that
% tag. Names are manuscript_style's, and so are the colors except that no
% adaptation is drawn in dark red (FR-004; see cond_color below).
%
% See also: explore_outliers, outlier_table, outlier_measures

arguments
    data_dir (1,:) char = ''
    fig_dir  (1,:) char = ''
    opts.seed double = []
    opts.dose_dir (1,:) char = 'default'
    opts.resolution (1,1) double = 200
end

root = fileparts(which('setup_paths'));
if isempty(data_dir)
    data_dir = fullfile(root, 'data', 'explorations', 'outliers');
end
if isempty(fig_dir)
    fig_dir = fullfile(root, 'figs', 'explorations', 'outliers');
end
if strcmp(opts.dose_dir, 'default')
    opts.dose_dir = fullfile(root, 'data', 'explorations', 'outliers_dose');
end
if ~isfolder(fig_dir); mkdir(fig_dir); end
delete(findall(groot, 'Type', 'figure', 'Tag', 'fig_outliers'));   % this function's own, nothing else

[T, M] = outlier_table(data_dir, 'print', false);
st = manuscript_style();
cond_order = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
conds = cond_order(ismember(cond_order, T.condition));
nets = unique(T.network, 'stable');
if any(strcmp(nets, 'reference'))                       % the reference first
    nets = [{'reference'}; setdiff(nets(:), {'reference'}, 'stable')];
end
seed1 = opts.seed;
if isempty(seed1); seed1 = max(T.seed); else; seed1 = seed1(1); end
n_net = numel(nets); n_cond = numel(conds);
pick = @(net, cond, seed) find(strcmp(T.network, net) & strcmp(T.condition, cond) & T.seed == seed, 1);
gray = [0.55 0.55 0.55];
fs_t = st.tick_fs - 2; fs_l = st.label_fs - 2;
files = {};
T1 = T(T.seed == seed1, :);

% The runs of the single-seed figures (spectra are not in the measures).
Rs = cell(n_net, n_cond);
for w = 1:n_net
    for c = 1:n_cond
        k = pick(nets{w}, conds{c}, seed1);
        if isempty(k); continue; end
        L = load(fullfile(data_dir, outlier_run_name(nets{w}, conds{c}, [seed1 seed1 + 1], '')), 'R');
        Rs{w, c} = L.R;
    end
end

%% 1. The complex plane
fig = new_fig([1 1 12 4.1 * n_net]);
tl = tiledlayout(fig, n_net, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
for w = 1:n_net
    axs = gobjects(1, n_cond); xr = [-22 5];
    for c = 1:n_cond
        ax = nexttile(tl); hold(ax, 'on'); axs(c) = ax;
        R = Rs{w, c};
        if ~isempty(R)
            m = M{pick(nets{w}, conds{c}, seed1)};
            col = cond_color(st, conds{c});
            tau = R.tau_d;
            ev = ([R.cal.ev_M] - 1) ./ tau;                    % n x n_cal, rates
            th = linspace(0, 2 * pi, 200);
            rad = m.kappa * median(R.R_eff) / tau;
            plot(ax, -1 / tau + rad * cos(th), rad * sin(th), '-', 'Color', gray, 'LineWidth', 1.5);
            scatter(ax, real(ev(:)), imag(ev(:)), 5, col, 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none');
            top = nan(1, size(ev, 2)) + 0i;
            for q = 1:size(ev, 2)
                [~, i1] = max(real(ev(:, q))); top(q) = ev(i1, q);
            end
            topc = with_conjugates(top);                       % a complex rightmost eigenvalue is one of a pair
            plot(ax, real(topc), imag(topc), 'o', 'MarkerSize', 6, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
            u = (median([R.cal.uniform_re]) - 1) / tau;
            plot(ax, u, 0, 'd', 'MarkerSize', 9, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.25);
            xline(ax, 0, '--', 'Color', 'k', 'LineWidth', 1.0);
            xr = [min(xr(1), -1 / tau - 1.15 * rad), max([xr(2), u + 3, max(real(top)) + 2])];
        end
        set(ax, 'FontSize', fs_t, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
        if w == 1; title(ax, st.condition_title(conds{c}), 'FontSize', st.title_fs, 'FontWeight', 'normal'); end
        if c == 1; ylabel(ax, sprintf('%s\nIm (s^{-1})', net_name(nets{w})), 'FontSize', fs_l); end
        if w == n_net; xlabel(ax, 'Re (s^{-1})', 'FontSize', fs_l); end
    end
    for c = 1:n_cond
        axis(axs(c), 'equal'); xlim(axs(c), xr); ylim(axs(c), [-1 1] * 0.44 * diff(xr));
    end
end
files{end + 1} = save_fig(fig, fig_dir, 'outliers_complex_plane', opts.resolution);
cols = {'lambda_1', 'kappa', 'agree_re', 'R_eff_med', 'rate_edge_med', 'rate_top_med', 'ratio_med', 'rho_uniform', ...
    'reined', 'frac_detached', 'K1_med', 'K_err_med', 'coh1_med', 'frac_on_E', 'frac_on_I', 'gain_ratio_EI', ...
    'n_unstable_mean', 'shuffle_pct_med'};
write_sidecar(fig_dir, 'outliers_complex_plane', ...
    ['Rows are the networks, columns the conditions. Dots: every eigenvalue of the frozen-gain effective connectivity ' ...
    '(each neuron''s outgoing weights times its current gain, slow variables held), drawn as a rate ' ...
    '(lambda - 1)/tau_d, at 12 states spread over 20 s after 60 s of settling, overlaid. Gray circle: the edge of the ' ...
    'bulk, centered on -1/tau_d = -10 1/s: the radius predicted from the gains (median over the window) times the ' ...
    'null-calibrated factor kappa (95th percentile of the rightmost real part of null matrices that have the same ' ...
    'entries and gains and no mean structure). Ringed dots: the rightmost eigenvalue of each state; when it is ' ...
    'complex it is one of a conjugate pair and both members are ringed. Open diamond: where ' ...
    'the rightmost eigenvalue would be if every neuron had the same gain (the rms of the observed gains, so the bulk ' ...
    'keeps its size): the outlier the weights alone make. Dashed: zero. The pile at -10 is the silent and saturated ' ...
    'neurons (zero gain). In the reference network the diamond is inside the circle (balanced means, no outlier) and ' ...
    'the ringed dots are the edge of the bulk. With E-to-E times 1.5 the diamond is far to the right (the structural ' ...
    'outlier, about five bulk radii out) and in every condition the actual rightmost eigenvalue is back at the bulk: ' ...
    'inside it without adaptation and with one timescale, where most excitatory neurons are pinned at the ceiling ' ...
    '(frac_on_E), and just outside it with multiple timescales, where they are not. Seed ' sprintf('%d', seed1) '.'], T1, cols);

%% 2. The summary: detachment, and the rightmost rate frozen and free
fig = new_fig([1 1 4.2 * n_net 7]);
tl = tiledlayout(fig, 2, n_net, 'TileSpacing', 'compact', 'Padding', 'compact');
for row = 1:2
    for w = 1:n_net
        ax = nexttile(tl); hold(ax, 'on');
        for c = 1:n_cond
            sel = find(strcmp(T.network, nets{w}) & strcmp(T.condition, conds{c}));
            col = cond_color(st, conds{c});
            off = linspace(-0.12, 0.12, numel(sel)) .* (numel(sel) > 1);
            for q = 1:numel(sel)
                k = sel(q); x = c + off(q);
                if row == 1
                    plot(ax, [x x], [T.ratio_p05(k) T.ratio_p95(k)], '-', 'Color', col, 'LineWidth', 2.0);
                    plot(ax, x, T.ratio_med(k), 'o', 'MarkerSize', 9, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', 'LineWidth', 1.0);
                    plot(ax, x, T.rho_uniform(k), 'd', 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.25);
                else
                    plot(ax, [x x] - 0.06, [T.rate_top_med(k) T.rate_top_p95(k)], '-', 'Color', col, 'LineWidth', 1.0);
                    plot(ax, x - 0.06, T.rate_top_med(k), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', col, 'LineWidth', 1.5);
                    plot(ax, [x x] + 0.06, [T.full_abscissa_p05(k) T.full_abscissa_p95(k)], '-', 'Color', col, 'LineWidth', 2.0);
                    plot(ax, x + 0.06, T.full_abscissa_med(k), 'o', 'MarkerSize', 9, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', 'LineWidth', 1.0);
                end
            end
        end
        if row == 1
            yline(ax, 1, '-', 'Color', gray, 'LineWidth', 1.5);
            set(ax, 'YScale', 'log', 'YLim', [0.4 8], 'YTick', [0.5 1 2 4 8]);
            title(ax, net_name(nets{w}), 'FontSize', st.title_fs, 'FontWeight', 'normal');
            if w == 1; ylabel(ax, 'detachment  Re \lambda_1 / bulk edge', 'FontSize', fs_l); end
        else
            yline(ax, 0, '--', 'Color', 'k', 'LineWidth', 1.0);
            if w == 1; ylabel(ax, 'rightmost rate (s^{-1})', 'FontSize', fs_l); end
        end
        set(ax, 'XTick', 1:n_cond, 'XTickLabel', cellfun(@(q) st.condition_short(q), conds, 'UniformOutput', false), ...
            'XLim', [0.4, n_cond + 0.6], 'FontSize', fs_t, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
    end
end
files{end + 1} = save_fig(fig, fig_dir, 'outliers_summary', opts.resolution);
cols = {'lambda_1', 'ratio_med', 'ratio_p05', 'ratio_p95', 'rho_uniform', 'reined', 'frac_on_E', 'frac_on_I', ...
    'rate_top_med', 'rate_top_p95', 'frac_unstable', 'full_abscissa_med', 'full_abscissa_p05', 'full_abscissa_p95', ...
    'full_frac_unstable', 'X_full', 'X', 'cont_re_med', 'cont_im_med', 'cont_ov_med', 'full_states'};
write_sidecar(fig_dir, 'outliers_summary', ...
    ['The summary. Columns are the networks; within a panel the conditions, two seeds each. TOP: detachment, the real ' ...
    'part of the rightmost eigenvalue of the frozen-gain connectivity divided by the null-calibrated edge of the bulk ' ...
    '(dot: median over 20 s; bar: 5th to 95th percentile; log axis). At 1 (gray line) the rightmost eigenvalue is the ' ...
    'edge of the bulk. Open diamonds: the same ratio with every gain made equal, which is what the weights alone give ' ...
    '(below 1 in the reference network: no structural outlier; near 5 with E-to-E times 1.5). The distance from the ' ...
    'diamond down to the dot is what the gains reined in; `reined` is that as a fraction (1: back at the edge). ' ...
    'BOTTOM: the rightmost rate with the slow variables frozen (open; bar from the median to the 95th percentile) and ' ...
    'with them free, the rightmost real part of the full Jacobian (filled; bar: 5th to 95th percentile). Dashed: zero. ' ...
    'Read the two rows with frac_on_E: detachment at or below 1 with frac_on_E near 0 is a pinned population, not a ' ...
    'controlled one. X is the exposure measure first proposed (mean of the frozen rightmost rate beyond both the ' ...
    'edge and zero); it is kept in the table to show why it was dropped: it is largest for multiple timescales at ' ...
    'E-to-E times 1.5, where the full Jacobian has no unstable eigenvalue at any sampled state. X_full is the mean ' ...
    'positive part of the full Jacobian''s rightmost real part. The continuation columns are the full-Jacobian ' ...
    'eigenvalue that carries the frozen rightmost pattern (median over 12 states).'], T, cols);

%% 3. Time series
fig = new_fig([1 1 12 3.2 * n_net]);
tl = tiledlayout(fig, n_net, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
for w = 1:n_net
    for c = 1:n_cond
        ax = nexttile(tl); hold(ax, 'on');
        k = pick(nets{w}, conds{c}, seed1);
        if ~isempty(k)
            m = M{k};
            col = cond_color(st, conds{c});
            plot(ax, m.t, m.rate_edge, '-', 'Color', gray, 'LineWidth', 1.5);
            plot(ax, m.t, m.rate_top, '-', 'Color', col, 'LineWidth', 1.0);
            plot(ax, m.full.t, m.full.abscissa, '-', 'Color', 'k', 'LineWidth', 1.25);
            yline(ax, 0, '--', 'Color', 'k', 'LineWidth', 1.0);
            xlim(ax, [m.t(1) m.t(end)]);
            v = [m.rate_edge, m.rate_top, m.full.abscissa, 0];
            pad = 0.08 * max(range_(v), 0.5);
            ylim(ax, [min(v) - pad, max(v) + pad]);
        end
        set(ax, 'FontSize', fs_t, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
        if w == 1; title(ax, st.condition_title(conds{c}), 'FontSize', st.title_fs, 'FontWeight', 'normal'); end
        if c == 1; ylabel(ax, sprintf('%s\nrate (s^{-1})', net_name(nets{w})), 'FontSize', fs_l); end
        if w == n_net; xlabel(ax, 'time (s)', 'FontSize', fs_l); end
    end
end
files{end + 1} = save_fig(fig, fig_dir, 'outliers_timeseries', opts.resolution);
cols = {'lambda_1', 'rate_edge_med', 'rate_top_med', 'rate_top_p95', 'rate_top_max', 'excess_med', 'ratio_med', ...
    'frac_detached', 'frac_unstable', 'full_abscissa_med', 'full_frac_unstable', 'full_minus_fast_med', ...
    'K1_med', 'K_err_med', 'K_corr', 'turnover_per_s', 'N_eff', 'full_states'};
write_sidecar(fig_dir, 'outliers_timeseries', ...
    ['Rows are the networks, columns the conditions. Colored: the rate (Re lambda - 1)/tau_d of the rightmost ' ...
    'eigenvalue of the frozen-gain effective connectivity, every 10 ms over 20 s after 60 s of settling. Gray: the ' ...
    'edge of the bulk (predicted radius at each state times the run''s kappa). Black: the rightmost real part of the ' ...
    'FULL Jacobian (slow variables free), every 100 ms (it stops early where the time budget ran out: full_states of ' ...
    '201). Dashed: zero. Colored above gray: an eigenvalue outside the bulk. Colored above zero with black at or below ' ...
    'zero: a mode that grows only if adaptation and depression are held still. The black line cannot go below the ' ...
    'slowest single-neuron adaptation rate (about -0.05 1/s with multiple timescales, -0.2 to -0.4 with one), which ' ...
    'is why it is flat where nothing is unstable. The K_eff prediction of the outlier (K1_med, in units of W, with its ' ...
    'median error where it predicts an outlier) is in the table. The y axes differ between panels. Seed ' ...
    sprintf('%d', seed1) '.'], T1, cols);

%% 4. What the slow variables do to the rightmost frozen-gain mode
fig = new_fig([1 1 12 3.6 * n_net]);
tl = tiledlayout(fig, n_net, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
for w = 1:n_net
    for c = 1:n_cond
        ax = nexttile(tl); hold(ax, 'on');
        R = Rs{w, c};
        if ~isempty(R) && ~isempty(R.dense)
            col = cond_color(st, conds{c});
            D = R.dense;
            lam = double([D.lam]); xf = double([D.xfrac]);
            net = lam(xf >= 0.1);
            scatter(ax, real(net), imag(net), 5, col, 'filled', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', 'none');
            xline(ax, 0, '--', 'Color', 'k', 'LineWidth', 1.0);
            fr = ([D.lam1_M] - 1) ./ R.tau_d;
            frc = with_conjugates(fr);                         % both members of a complex pair
            plot(ax, real(frc), imag(frc), 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.25);
            ct = [D.cont];
            plot(ax, [real(ct) real(ct)], [imag(ct) -imag(ct)], 'o', 'MarkerSize', 7, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
        end
        axis(ax, 'equal'); xlim(ax, [-14 10]); ylim(ax, [-9 9]);
        set(ax, 'FontSize', fs_t, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
        if w == 1; title(ax, st.condition_title(conds{c}), 'FontSize', st.title_fs, 'FontWeight', 'normal'); end
        if c == 1; ylabel(ax, sprintf('%s\nIm (s^{-1})', net_name(nets{w})), 'FontSize', fs_l); end
        if w == n_net; xlabel(ax, 'Re (s^{-1})', 'FontSize', fs_l); end
    end
end
files{end + 1} = save_fig(fig, fig_dir, 'outliers_slow_variables', opts.resolution);
cols = {'rate_top_med', 'cont_re_med', 'cont_im_med', 'cont_ov_med', 'dense_n_pos_med', 'dense_frac_unstable', ...
    'full_abscissa_med', 'full_frac_unstable', 'full_dwell_zero_med', 'full_dwell_zero_max', 'full_n_exc_zero', ...
    'X_full', 'kin_edge', 'lambda_1'};
write_sidecar(fig_dir, 'outliers_slow_variables', ...
    ['The right part of the spectrum of the FULL Jacobian (adaptation and depression variables free to respond) at ' ...
    'the same 12 states as the complex-plane figure, overlaid. Dots: eigenvalues whose eigenvector has at least 10% ' ...
    'of its squared norm in the dendritic variables (network modes; the thousands of single-neuron slow modes on the ' ...
    'negative real axis are left out). Open circles: the rightmost eigenvalue of the frozen-gain matrix at each state ' ...
    '(what the other figures follow), with its conjugate when it is complex. Filled ringed circles: the eigenvalue of the full Jacobian whose dendritic ' ...
    'pattern overlaps most with that frozen mode (cont_ov_med), with its conjugate. Dashed: zero. Without adaptation ' ...
    'there are no slow variables and the two coincide. With adaptation the frozen mode, real or nearly so and right ' ...
    'of zero, becomes a complex pair: the slow variables turn growth into a slow oscillation (Im of 4 to 6 1/s is ' ...
    '0.6 to 0.9 Hz) whose real part is near zero in the reference network and negative with E-to-E times 1.5. ' ...
    'The view is cut at -14 and +10 1/s. Seed ' sprintf('%d', seed1) '.'], T1, cols);

%% 5. Dwell times beyond zero
fig = new_fig([1 1 5.2 * n_net 4.2]);
tl = tiledlayout(fig, 1, n_net, 'TileSpacing', 'compact', 'Padding', 'compact');
for w = 1:n_net
    ax = nexttile(tl); hold(ax, 'on');
    for c = 1:n_cond
        for k = find(strcmp(T.network, nets{w}) & strcmp(T.condition, conds{c}))'
            m = M{k};
            col = cond_color(st, conds{c});
            survival(ax, m.dwell.zero, m.dt, '--', col, 1.25);
            survival(ax, m.dwell.full_zero, 0.1, '-', col, 2.25);
        end
    end
    set(ax, 'XScale', 'log', 'XLim', [0.01 25], 'YLim', [0 1.02], 'XTick', [0.01 0.1 1 10], 'FontSize', fs_t, ...
        'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
    title(ax, net_name(nets{w}), 'FontSize', st.title_fs, 'FontWeight', 'normal');
    xlabel(ax, 'time beyond zero (s)', 'FontSize', fs_l);
    if w == 1; ylabel(ax, 'fraction of excursions longer', 'FontSize', fs_l); end
end
files{end + 1} = save_fig(fig, fig_dir, 'outliers_dwell', opts.resolution);
cols = {'frac_unstable', 'n_exc_zero', 'dwell_zero_med', 'dwell_zero_max', 'full_frac_unstable', 'full_n_exc_zero', ...
    'full_dwell_zero_med', 'full_dwell_zero_max', 'frac_detached', 'n_exc_bulk', 'dwell_bulk_med', 'dwell_bulk_max', ...
    'frac_same', 'base_overlap', 'turnover_per_s', 'N_eff', 'coh1_med'};
write_sidecar(fig_dir, 'outliers_dwell', ...
    ['Survival curves of excursion lengths: the fraction of excursions beyond zero that last longer than the time on ' ...
    'the horizontal axis. Dashed: the rightmost eigenvalue of the frozen-gain matrix (10 ms resolution). Solid: the ' ...
    'rightmost eigenvalue of the full Jacobian (100 ms resolution). Colors are the conditions; both seeds are drawn. ' ...
    'A curve is absent where there is no excursion (with E-to-E times 1.5 the full Jacobian is never beyond zero in ' ...
    'any condition, and the frozen rightmost eigenvalue only with multiple timescales). A run that touches either end ' ...
    'of the 20 s window is drawn at its observed length, so a condition that is beyond zero throughout is one step at ' ...
    '20 s. Dwell beyond the EDGE of the bulk is in the table, not the figure: without a structural outlier it is ' ...
    'blips of 10 to 40 ms at the rate a 95th-percentile edge gives by chance, and with a standing outlier it is the ' ...
    'whole window.'], T, cols);

%% 6. Across multiples of the E-to-E mean
if isfolder(opts.dose_dir) && ~isempty(dir(fullfile(opts.dose_dir, 'run_*.mat')))
    Td = outlier_table(opts.dose_dir, 'print', false);
    Ta = [Td; T1];
    [~, order] = sort(Ta.structural_ratio); Ta = Ta(order, :);
    fig = new_fig([1 1 13 4.2]);
    tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    ycol = {'ratio_med', 'frac_on_E', 'gain_ratio_EI'};
    ylab = {'detachment  Re \lambda_1 / bulk edge', 'responsive fraction of E neurons', 'mean gain, E / I'};
    for p = 1:3
        ax = nexttile(tl); hold(ax, 'on');
        if p == 1
            u = Ta(strcmp(Ta.condition, conds{end}), :);
            plot(ax, u.structural_ratio, u.rho_uniform, '-d', 'Color', 'k', 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'LineWidth', 1.25);
            yline(ax, 1, '-', 'Color', gray, 'LineWidth', 1.5);
        end
        for c = 1:n_cond
            u = Ta(strcmp(Ta.condition, conds{c}), :);
            plot(ax, u.structural_ratio, u.(ycol{p}), '-o', 'Color', cond_color(st, conds{c}), 'LineWidth', 2.0, ...
                'MarkerSize', 8, 'MarkerFaceColor', cond_color(st, conds{c}), 'MarkerEdgeColor', 'w');
        end
        if p == 1; set(ax, 'YScale', 'log', 'YLim', [0.4 8], 'YTick', [0.5 1 2 4 8]); else; ylim(ax, [0 1.1]); end
        set(ax, 'FontSize', fs_t, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out');
        xlabel(ax, 'structural outlier / bulk radius', 'FontSize', fs_l);
        ylabel(ax, ylab{p}, 'FontSize', fs_l);
    end
    files{end + 1} = save_fig(fig, fig_dir, 'outliers_dose', opts.resolution);
    cols = {'structural_ratio', 'lambda_1', 'ratio_med', 'rho_uniform', 'reined', 'frac_detached', 'frac_on_E', ...
        'frac_on_I', 'gain_ratio_EI', 'rate_E', 'rate_I', 'K1_med', 'coh1_med', 'rate_edge_med', 'rate_top_med', ...
        'frac_unstable', 'N_eff'};
    write_sidecar(fig_dir, 'outliers_dose', ...
        ['The E-to-E mean at 1, 1.1, 1.25, 1.5 and 1.75 times the reference, placed on the horizontal axis by what it ' ...
        'does to the weights: the structural outlier (largest real part of the eigenvalues of the realized block-mean ' ...
        'matrix) over the structural bulk radius. Below 1 the outlier is inside the bulk. LEFT: detachment of the ' ...
        'rightmost frozen-gain eigenvalue from the calibrated edge (median over the window; log axis; gray line: the ' ...
        'edge). Open diamonds: the same with uniform gains, the weights alone. MIDDLE: the fraction of excitatory ' ...
        'neurons off the floor and ceiling of the nonlinearity (nonzero slope). RIGHT: mean gain of the excitatory ' ...
        'neurons over that of the inhibitory ones. Every condition holds the detachment at or below about 1.3 ' ...
        'whatever the weights do; they differ in how. Without adaptation the excitatory population is pinned as soon ' ...
        'as there is a structural outlier; with one timescale it is lost progressively; with multiple timescales it ' ...
        'stays responsive to a structural outlier of about five bulk radii while the E-to-I gain ratio falls in ' ...
        'proportion. Seed ' sprintf('%d', seed1) ' only; the three added levels (1.1, 1.25, 1.75) were sampled every ' ...
        '50 ms with 8 calibration states and no full-Jacobian pass, so their full-Jacobian columns are not given.'], Ta, cols);
end
end

%% ------------------------------------------------------------------------
function fig = new_fig(pos)
fig = figure('Color', 'w', 'Units', 'inches', 'Position', pos, 'Tag', 'fig_outliers');
end

function col = cond_color(st, name)
% The color every figure here draws a condition in. FR-004 (Tom,
% 2026-11-20): no adaptation is the dark red of the Figure 3 Ritz overlay
% (the value fig_eig_heatmap's draw_ritz substitutes for a colorless
% condition color), not manuscript_style's black, so that it cannot be
% confused with the marks that are black for another reason (the full
% Jacobian, the uniform-gain diamonds, zero lines). The shared map in
% manuscript_style is deliberately not changed here: the paper's figures
% change together on the next full figures run, and this helper goes then.
if strcmp(name, 'no_adaptation')
    col = [0.75 0.10 0.15];
else
    col = st.condition_color(name);
end
end

function name = net_name(net)
tok = regexp(net, '^EEx(\d+)p(\d+)$', 'tokens', 'once');
if ~isempty(tok)
    name = sprintf('E-to-E \\times %s.%s', tok{1}, tok{2});
elseif strcmp(net, 'EI3to2')
    name = 'E:I 3:2';
else
    name = net;
end
end

function z = with_conjugates(z)
% The eigenvalues z with the conjugate of each complex one appended, so that
% both members of a pair are marked. A real one (|Im| below 1e-9 of the
% larger of its modulus and 1) is kept once, on the real axis.
z = z(:).';
is_real = abs(imag(z)) <= 1e-9 * max(abs(z), 1);
z(is_real) = real(z(is_real));
z = [z, conj(z(~is_real))];
end

function survival(ax, dwell, dt, style, col, lw)
% Empirical survival of the dwell times, as stairs from dt.
if isempty(dwell); return; end
d = sort(dwell(:));
frac = (numel(d):-1:1)' ./ numel(d);                   % fraction with dwell >= d(i)
stairs(ax, [dt; d; d(end)], [1; frac; 0], style, 'Color', col, 'LineWidth', lw);
end

function r = range_(v)
r = max(v) - min(v);
end

function file = save_fig(fig, fig_dir, name, resolution)
file = fullfile(fig_dir, [name '.png']);
exportgraphics(fig, file, 'Resolution', resolution);
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
