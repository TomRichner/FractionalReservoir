function out = fig_example_timeseries(cfg)
% FIG_EXAMPLE_TIMESERIES Example SRNN time series at the paper's operating point.
%
%   out = FIG_EXAMPLE_TIMESERIES()
%   out = FIG_EXAMPLE_TIMESERIES('preset_name', p, 'out_dir', d)
%
% One realization of the network the rest of the paper analyses -- the same
% W, stimulus and setpoints (rng_seeds) under EVERY adaptation condition of
% the preset, one column per condition, rows sharing the time axis:
%
%   (a) dendritic state x of a subset of neurons (n_show_E E, n_show_I I;
%       25 + 25 = 50 per condition by default, TR 2026-09-14)
%   (b) firing rate r of the same neurons
%   (c) synaptic output r * prod(b) * prod(g) of the same neurons, route 1.
%       The preset must carry the same synaptic dynamics on every route
%       (SRNNCellTypePairs.routes_identical: route 1 is then the same array
%       as every other route); a per-route preset errors rather than
%       drawing a fallback nobody reads at 50 neurons
%   (d) the adaptation term (c/K) sum_k a_k of the same neurons, per neuron,
%       coloured by cell type (was the population mean of each a_k until
%       2026-09-14; the per-neuron term is what enters the rate)
%   (e) the depression factor of the same neurons (route 1, same rule as
%       (c)), from SRNNCellTypePairs.std_route_factor: prod(b) on a product
%       route, the releasable pool b_1 on a chained route (FD-006); a note
%       in each panel says which
%   (f) the leading LOCAL expansion rate (grey) and the ACCUMULATING
%       finite-time lambda_1 (condition colour), plot_local_and_finite_lle.
%
% WHAT CHANGED AND WHY (2026-09-14, Codex manuscript audit sec. 3 and 10).
% The figure used to be model.plot() under the most-adapted condition alone,
% with lya_method 'none': every neuron drawn, one panel per route for the
% synaptic output even when all four routes were identical, and no Lyapunov
% panel at all -- so the "here is what the model does" figure showed none of
% the quantities the Results are about. It now shows a readable subset, one
% synaptic-output trace per neuron where that is well defined, all three
% conditions side by side on one network, and the two Lyapunov quantities the
% paper keeps apart: local positive stretches are "transient expansion", the
% accumulated exponent's end value is what "stable" or "chaotic" refers to.
%
% The Lyapunov estimator is the sweeps' (top-K, K = 15, lya_dt 0.05,
% accumulation over [5, T]), so the printed lambda_1 is directly comparable
% with the sweep figures. Conditions without SFA or STD get an empty panel
% saying so rather than a missing row, so the columns align.
%
% cfg.condition, if given, restricts the figure to that one condition (the
% old single-column behaviour); empty draws all of them. cfg.conditions
% (FD-006) is the bundle-wide option: names drawn as columns in the order
% given, resolved by select_conditions against the preset's conditions;
% only the selected conditions are simulated. Give one or the other.
%
% cfg.model_args (2026-11-12) is a cell of name-value pairs handed to
% build_from_preset after this figure's own settings, so they override the
% preset for every condition drawn. It is for looking at a network
% interactively before a long run (another mu_tilde_relative, a stimulus
% in place of the preset's none) without adding a preset for each try:
%   fig_example_timeseries('preset_name', p, 'save', false, ...
%       'model_args', {'mu_tilde_relative', [6 -6; 6 -6]})
% Empty, the default, changes nothing. A figure made with it is not the
% preset's network; the registry never sets it.
%
% Two opt-in extras (2026-11-12), each a SEPARATE figure appended to
% out.figs after the main one and saved under its own tag. Both default
% false, so a default call and every registry entry make one figure, named
% as before. Both are drawn from the first condition's model: one seed
% gives every condition the same input and the same W.
%   cfg.show_stimulus   -> 'fig_example_stimulus': the external input u
%       (u_ex, which includes any intrinsic_drive). An image of u, neurons
%       by time with the types stacked in type order (E above I), zero
%       white; below it the trace of every stimulated neuron (u nonzero at
%       any sample) in its type's shades. The titles give the count and
%       fraction stimulated per type, the on intervals and max |u|.
%       out.stimulus holds the counts. With no input at all the figure
%       says so in words instead of drawing an empty axis.
%   cfg.show_W_spectrum -> 'fig_example_W_spectrum': eig of the W that is
%       simulated (get_params().W: level_of_chaos and any route scale
%       folded in), with the class's bulk radius R and outliers lambda_O,
%       and the spectral abscissa and radius in the title. R and lambda_O
%       describe the matrix before route scales, so they are left out (and
%       the title says why) when a scale or rescale_by_abscissa is on. If
%       the simulated W differs between the conditions drawn, the title
%       says so and the first is drawn. out.W_spectrum holds the numbers.
% For example, everything for one seed:
%   fig_example_timeseries('preset_name', p, 'save', false, ...
%       'show_stimulus', true, 'show_W_spectrum', true)
%
% Two options for the lambda row (2026-11-12), both off by default, so a
% default call draws and computes exactly what it did:
%   cfg.lambda_K    0 (default): the row is the local rate and the
%       accumulating lambda_1, from the sweeps' K = 15 estimator. K > 0:
%       the estimator is asked for the top K exponents (lya_K = K) and the
%       row draws ALL K accumulating finite-time exponents, lambda_1 thick
%       in the condition's colour and lambda_2..K thin beneath it, with the
%       final lambda_1 and how many of the K ended positive printed in the
%       corner. The local rate is not drawn in that form. Spectra nest
%       across K, so lambda_1 is the K = 15 value.
%   cfg.lambda_link   [] (default): the row's y axis is shared by every
%       column, like the other rows. A vector of column numbers: only those
%       columns share their y axis and every other column keeps its own
%       scale. [2 3 4] on a four-condition preset leaves the no-adaptation
%       column, whose exponents are an order of magnitude larger, on its
%       own scale so the three adapting columns can be read.
%
% See also: paper_config, srnn_param_preset, SRNNCellTypePairs,
%           plot_local_and_finite_lle, fig_lyapunov_spectrum

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.preset_name (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    cfg.out_dir     (1,:) char    = ''
    cfg.condition   (1,:) char    = ''      % '' -> every condition of the preset
    cfg.rng_seeds   (1,2) double  = [1 2]
    cfg.T_range     (1,2) double  = [0 20]
    cfg.n_show_E    (1,1) double  = 25
    cfg.n_show_I    (1,1) double  = 25
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''      % unused; accepted for a uniform call
    cfg.conditions                = {}      % names to draw, in order; empty = all the preset holds (select_conditions)
    cfg.model_args  (1,:) cell    = {}      % name-value overrides of the preset, for interactive looks (header)
    cfg.show_stimulus   (1,1) logical = false   % extra figure: the external input u (header)
    cfg.show_W_spectrum (1,1) logical = false   % extra figure: eigenvalues of the simulated W (header)
    cfg.lambda_K    (1,1) double {mustBeInteger, mustBeNonnegative} = 0   % > 0: the lambda row draws the top K exponents (header)
    cfg.lambda_link (1,:) double {mustBeInteger, mustBePositive} = []    % columns sharing the lambda row's y axis; [] = all (header)
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();

[~, ~, conditions] = srnn_param_preset(cfg.preset_name);
cond_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
if ~isempty(cfg.condition)
    if ~any(strcmp(cond_names, cfg.condition))
        error('fig_example_timeseries:NoSuchCondition', ...
            'Preset ''%s'' has no condition ''%s'' (has %s).', ...
            cfg.preset_name, cfg.condition, strjoin(cond_names, ', '));
    end
    cond_names = {cfg.condition};
    if ~isempty(cfg.conditions)
        error('fig_example_timeseries:BothSelections', ...
            'Give condition or conditions, not both.');
    end
end
cond_names = select_conditions(cond_names, cfg.conditions);
n_cond = numel(cond_names);
if any(cfg.lambda_link > n_cond)
    error('fig_example_timeseries:LambdaLink', ...
        'lambda_link names column %d; the figure has %d column(s).', max(cfg.lambda_link), n_cond);
end

%% Simulate: one network, every condition, the sweeps' Lyapunov estimator
T = cfg.T_range(2);
lya_T_interval = [min(5, T / 2), T];
lya_K = 15;                                   % the sweeps' K
if cfg.lambda_K > 0, lya_K = cfg.lambda_K; end
models = cell(1, n_cond);
for i = 1:n_cond
    models{i} = build_from_preset(cfg.preset_name, cond_names{i}, ...
        'rng_seeds', cfg.rng_seeds, 'T_range', cfg.T_range, 'fs', 400, ...
        'lya_method', 'topk', 'lya_K', lya_K, 'lya_K_auto', false, 'lya_dt', 0.05, ...
        'lya_T_interval', lya_T_interval, 'lya_warmup', min(5, T / 4), ...
        'verbose', cfg.verbose, cfg.model_args{:});
    models{i}.run();
    vprintf(cfg.verbose, 'verbose', '  [example_timeseries] %-14s lambda_1 = %+.4f\n', ...
        cond_names{i}, models{i}.lya_results.LLE);
end

%% Figure
row_names = {'x', 'r', 'synaptic output', 'SFA  (c/K)\Sigma a', 'depression factor', '\lambda'};
n_rows = numel(row_names);
fig = figure('Color', 'w', 'Position', [60 40 420 * n_cond + 80, 1150]);
tl = tiledlayout(fig, n_rows, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
tl.TileIndexing = 'columnmajor';
ax_all = gobjects(n_rows, n_cond);

for i = 1:n_cond
    m      = models{i};
    pd     = m.plot_data;
    params = m.get_params();
    t      = pd.t;
    names  = m.cell_type_names;
    C      = numel(names);
    tcol   = SRNNCellTypePairs.type_colors(C);
    % Subset: the first n_show of each type, and per-neuron shades of the type's hue.
    n_show = zeros(1, C);
    shades = cell(1, C);
    for q = 1:C
        want = cfg.n_show_I;
        if q == 1, want = cfg.n_show_E; end
        n_show(q) = min(want, size(pd.x.(names{q}), 1));
        shades{q} = SRNNCellTypePairs.neuron_colors(tcol(q, :), max(n_show(q), 1));
    end
    if ~SRNNCellTypePairs.routes_identical(params)
        error('fig_example_timeseries:RoutesDiffer', ...
            'Preset ''%s'' (%s) carries different synaptic dynamics on different routes; this figure draws route 1 per neuron and needs them identical.', ...
            cfg.preset_name, cond_names{i});
    end
    cond_col = st.condition_color(cond_names{i});
    ttl = cond_names{i};
    if isKey(st.condition_title, cond_names{i}), ttl = st.condition_title(cond_names{i}); end

    % (a) x, (b) r ------------------------------------------------------------
    for row = 1:2
        ax = nexttile(tl); ax_all(row, i) = ax; hold(ax, 'on');
        if row == 1, field = 'x'; else, field = 'r'; end
        for q = SRNNCellTypePairs.draw_order(C)
            D = pd.(field).(names{q});
            for k = 1:n_show(q)
                plot(ax, t, D(k, :), '-', 'Color', shades{q}(k, :), 'LineWidth', 0.5);
            end
        end
        if row == 1
            title(ax, ttl, 'FontWeight', 'normal', 'FontSize', st.title_fs);
            hl = gobjects(1, C);
            for q = 1:C
                hl(q) = plot(ax, NaN, NaN, '-', 'Color', tcol(q, :), 'LineWidth', 1.5);
            end
            legend(ax, hl, arrayfun(@(q) sprintf('%s (%d shown)', names{q}, n_show(q)), 1:C, ...
                'UniformOutput', false), 'Location', 'northeast', 'FontSize', 8, 'Box', 'off');
        end
        hold(ax, 'off');
    end

    % (c) synaptic output ------------------------------------------------------
    ax = nexttile(tl); ax_all(3, i) = ax; hold(ax, 'on');
    so = pd.synaptic_output;
    for q = SRNNCellTypePairs.draw_order(C)
        posts = fieldnames(so.(names{q}));
        D = so.(names{q}).(posts{1});             % route 1 == every route
        for k = 1:n_show(q)
            plot(ax, t, D(k, :), '-', 'Color', shades{q}(k, :), 'LineWidth', 0.5);
        end
    end
    hold(ax, 'off');

    % (d) adaptation term (c/K) sum_k a_k per neuron ---------------------------
    ax = nexttile(tl); ax_all(4, i) = ax; hold(ax, 'on');
    drew = false;
    for q = SRNNCellTypePairs.draw_order(C)
        A = pd.a.(names{q});                     % n_q x n_a x n_t, or []
        if isempty(A), continue; end
        Sa = params.c_eff(q) * reshape(sum(A, 2), size(A, 1), size(A, 3));   % n_q x n_t
        for k = 1:n_show(q)
            plot(ax, t, Sa(k, :), '-', 'Color', shades{q}(k, :), 'LineWidth', 0.5);
        end
        drew = true;
    end
    if ~drew
        text(ax, 0.5, 0.5, 'no SFA', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 10);
        set(ax, 'YTick', []);
    end
    hold(ax, 'off');

    % (e) the depression factor -----------------------------------------------
    % From the model (SRNNCellTypePairs.std_route_factor), never prod(b)
    % formed here: prod(b) on a product route, the releasable pool b1 on a
    % chained one. The panel says which.
    ax = nexttile(tl); ax_all(5, i) = ax; hold(ax, 'on');
    drew = false; is_chain = false; n_b_max = 0;
    for q = SRNNCellTypePairs.draw_order(C)
        posts = fieldnames(pd.b.(names{q}));
        B = pd.b.(names{q}).(posts{1});          % n_q x n_b x n_t, or []; route 1 == every route
        if isempty(B), continue; end
        post_q = find(strcmp(names, posts{1}), 1);
        Pb = SRNNCellTypePairs.std_route_factor(B, params, q, post_q);
        Pb = reshape(Pb, size(B, 1), size(B, 3));   % n_q x n_t
        is_chain = is_chain || params.std_chain(q, post_q);
        n_b_max = max(n_b_max, size(B, 2));
        for k = 1:n_show(q)
            plot(ax, t, Pb(k, :), '-', 'Color', shades{q}(k, :), 'LineWidth', 0.5);
        end
        drew = true;
    end
    if ~drew
        text(ax, 0.5, 0.5, 'no STD', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 10);
        set(ax, 'YTick', []);
    else
        ylim(ax, [0 1.02]);
        if is_chain
            note = 'b_1 (releasable pool)';
        elseif n_b_max > 1
            note = '\Pi_m b_m';
        else
            note = 'b (one pool)';
        end
        % just above the axes' top-right corner: 50 traces fill the inside
        text(ax, 1, 1.01, note, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', [0.3 0.3 0.3], 'Clipping', 'off');
    end
    hold(ax, 'off');

    % (f) local and finite-time lambda_1, or the top-K accumulating exponents ---
    ax = nexttile(tl); ax_all(6, i) = ax;
    if cfg.lambda_K > 0
        lambda_spectrum_row(ax, m.lya_results, cond_col, cfg.lambda_K);
    else
        plot_local_and_finite_lle(ax, m.lya_results, cond_col);
    end
    xlabel(ax, 'time (s)', 'FontSize', st.label_fs);
end

% Cosmetics shared across the grid
for row = 1:n_rows
    for i = 1:n_cond
        ax = ax_all(row, i);
        set(ax, 'FontSize', st.tick_fs - 3); box(ax, 'off');
        xlim(ax, cfg.T_range);
        if i == 1, ylabel(ax, row_names{row}, 'FontSize', st.label_fs - 2); end
        if row < n_rows, set(ax, 'XTickLabel', []); end
    end
    % Every row shares its y axis across the columns, except the lambda row
    % when lambda_link names the columns that do (the rest keep their own).
    link_cols = 1:n_cond;
    if row == n_rows && ~isempty(cfg.lambda_link)
        link_cols = unique(cfg.lambda_link);
    end
    if numel(link_cols) > 1
        linkaxes(ax_all(row, link_cols), 'y');
    end
end
title(tl, sprintf('%s   (seeds %s, T = %g s, %d E + %d I of %d neurons shown)', ...
    strrep(cfg.preset_name, '_', '\_'), mat2str(cfg.rng_seeds), T, cfg.n_show_E, cfg.n_show_I, models{1}.n), ...
    'FontWeight', 'normal', 'FontSize', 9);

if ~cfg.visible, set(fig, 'Visible', 'off'); end

%% Save
out = struct('figs', fig, 'files', {{}}, 'source', 'simulated inline', ...
    'lambda_1', cellfun(@(m) m.lya_results.LLE, models), 'cond_names', {cond_names});
if cfg.save
    save_figure_stable(out_dir, 'fig_example_timeseries', fig);
    out.files = existing_outputs(out_dir, 'fig_example_timeseries');
end

%% Opt-in extras: each its own figure, appended after the main one. Their
% tags are not prefixes of 'fig_example_timeseries' nor it of them, because
% save_figure_stable deletes <out_dir>/<tag>* before saving.
if cfg.show_stimulus
    [fig_s, out.stimulus] = stimulus_figure(models{1}, cfg, st);
    out.figs(end+1) = fig_s;
    if cfg.save
        save_figure_stable(out_dir, 'fig_example_stimulus', fig_s);
        out.files = [out.files, existing_outputs(out_dir, 'fig_example_stimulus')];
    end
end
if cfg.show_W_spectrum
    [fig_w, out.W_spectrum] = W_spectrum_figure(models, cond_names, cfg, st);
    out.figs(end+1) = fig_w;
    if cfg.save
        save_figure_stable(out_dir, 'fig_example_W_spectrum', fig_w);
        out.files = [out.files, existing_outputs(out_dir, 'fig_example_W_spectrum')];
    end
end
end

% =========================================================================
function lambda_spectrum_row(ax, lya_results, colour, K)
% The accumulating finite-time exponents of the top K directions against
% time (cfg.lambda_K): lambda_1 thick in the condition's colour, lambda_2..K
% thin beneath it, a zero line, and in the corner the final lambda_1 and how
% many of the K ended positive. Needs the top-K estimator's
% finite_LE_spectrum_t (n_t x K, NaN before accumulation starts).
if ~isfield(lya_results, 'finite_LE_spectrum_t')
    error('fig_example_timeseries:NoSpectrum', ...
        'lambda_K needs the top-K estimator''s finite_LE_spectrum_t; this result has none.');
end
F = lya_results.finite_LE_spectrum_t;
t = lya_results.t_lya(:);
n = min(numel(t), size(F, 1));
K = min(K, size(F, 2));
t = t(1:n); F = F(1:n, 1:K);
ok = isfinite(F(:, 1));
hold(ax, 'on');
% lambda_K first and lambda_1 last, so the leading exponents lie on top.
pale = colour + 0.6 * (1 - colour);
for k = K:-1:2
    plot(ax, t(ok), F(ok, k), '-', 'Color', pale, 'LineWidth', 0.5);
end
plot(ax, t(ok), F(ok, 1), '-', 'Color', colour, 'LineWidth', 2.0);
yline(ax, 0, ':', 'Color', [0.3 0.3 0.3]);
hold(ax, 'off');
ylabel(ax, '\lambda (s^{-1})');
% Limits from the settled part: the first fifth of the accumulation is the
% estimator converging, and it would set the scale for everything after.
i_ok = find(ok);
if ~isempty(i_ok)
    settled = i_ok(max(1, round(0.2 * numel(i_ok))):end);
    lo = min(F(settled, :), [], 'all'); hi = max(F(settled, :), [], 'all');
    lo = min(lo, 0); hi = max(hi, 0);
    pad = 0.08 * max(hi - lo, eps);
    ylim(ax, [lo - pad, hi + pad]);
    final = F(i_ok(end), :);
    % Just above the axes' top-right corner, as the depression panel's note:
    % inside, it would sit on the lambda_1 trace in one column or another.
    text(ax, 1, 1.01, sprintf('top %d:  \\lambda_1 = %+.3f,  %d > 0', K, final(1), nnz(final > 0)), ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'FontSize', 8, 'Color', colour, 'Clipping', 'off');
end
end

% =========================================================================
function [fig, info] = stimulus_figure(m, cfg, st)
% The external input u of the first model drawn (one seed: the same input
% under every condition). Top: image of u, neurons by time, types stacked
% in type order (E above I), zero white. Bottom: u of every stimulated
% neuron, in its type's shades. A neuron counts as stimulated when its u is
% nonzero at any sample, so a constant intrinsic_drive counts too; the
% title then says the input carries one.
t     = m.t_ex(:)';                      % 1 x nt
U     = m.u_ex;                          % n x nt
names = m.cell_type_names;
C     = numel(names);
tcol  = SRNNCellTypePairs.type_colors(C);
idx   = m.type_indices;
vis   = 'on'; if ~cfg.visible, vis = 'off'; end
fig   = figure('Color', 'w', 'Visible', vis, 'Position', [80 60 1000 720]);

stim = any(U ~= 0, 2);
n_stim = cellfun(@(ix) nnz(stim(ix)), idx);
n_type = cellfun(@numel, idx);
on = any(U ~= 0, 1);
info = struct('n_stimulated', nnz(stim), 'n_stimulated_per_type', n_stim, ...
    'n_per_type', n_type, 'type_names', {names}, 'on_intervals', on_intervals(t, on), ...
    'max_abs_u', max(abs(U(:))));

parts = arrayfun(@(q) sprintf('%s %d of %d (%.1f%%)', names{q}, n_stim(q), n_type(q), ...
    100 * n_stim(q) / max(n_type(q), 1)), 1:C, 'UniformOutput', false);
if ~any(stim)
    ax = axes(fig); axis(ax, 'off');
    text(ax, 0.5, 0.5, sprintf(['No external input: u = 0 for all %d neurons over [%g, %g] s.\n' ...
        '(the preset''s input; give model\\_args {''input\\_config'', ...} to add one)'], ...
        m.n, t(1), t(end)), 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', st.label_fs - 2);
    title(ax, 'External input u', 'FontWeight', 'normal', 'FontSize', st.title_fs);
    return;
end

const_drive = any(all(U ~= 0, 2) & max(U, [], 2) == min(U, [], 2));
iv = info.on_intervals;
iv_txt = strjoin(arrayfun(@(k) sprintf('%.2f-%.2f s', iv(k, 1), iv(k, 2)), ...
    1:size(iv, 1), 'UniformOutput', false), ', ');

tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
% (top) image of u, types stacked in type order
ax1 = nexttile(tl, [2 1]);
order = vertcat(idx{:});
umax = max(abs(U(:)));
imagesc(ax1, t, 1:m.n, U(order, :));
colormap(ax1, redwhiteblue_colormap(256));
clim(ax1, [-umax umax]);
cb = colorbar(ax1); cb.Label.String = 'u';
hold(ax1, 'on');
edges = cumsum(n_type);
for q = 1:C - 1
    yline(ax1, edges(q) + 0.5, 'k-', 'LineWidth', 1);
end
hold(ax1, 'off');
set(ax1, 'YTick', edges - n_type / 2, 'YTickLabel', names, 'FontSize', st.tick_fs - 3);
ylabel(ax1, 'neuron (by type)', 'FontSize', st.label_fs - 2);
set(ax1, 'XTickLabel', []);
title(ax1, sprintf('Stimulated: %s.  On: %s', strjoin(parts, ',  '), iv_txt), ...
    'FontWeight', 'normal', 'FontSize', st.title_fs - 3);

% (bottom) u of each stimulated neuron, in its type's shades
ax2 = nexttile(tl); hold(ax2, 'on');
hl = gobjects(1, C);
for q = SRNNCellTypePairs.draw_order(C)
    ks = idx{q}(stim(idx{q}));
    shades = SRNNCellTypePairs.neuron_colors(tcol(q, :), max(numel(ks), 1));
    for j = 1:numel(ks)
        plot(ax2, t, U(ks(j), :), '-', 'Color', shades(j, :), 'LineWidth', 0.75);
    end
end
for q = 1:C
    hl(q) = plot(ax2, NaN, NaN, '-', 'Color', tcol(q, :), 'LineWidth', 1.5);
end
hold(ax2, 'off');
legend(ax2, hl, arrayfun(@(q) sprintf('%s (%d)', names{q}, n_stim(q)), 1:C, ...
    'UniformOutput', false), 'Location', 'eastoutside', 'FontSize', 8, 'Box', 'off');
xlim(ax2, [t(1) t(end)]); box(ax2, 'off');
set(ax2, 'FontSize', st.tick_fs - 3);
xlabel(ax2, 'time (s)', 'FontSize', st.label_fs - 2);
ylabel(ax2, 'u of stimulated neurons', 'FontSize', st.label_fs - 2);
xlim(ax1, [t(1) t(end)]);

sup = sprintf('External input u  (seed %s, %d of %d neurons, |u| \\leq %.3g)', ...
    mat2str(cfg.rng_seeds), nnz(stim), m.n, umax);
if const_drive, sup = [sup '   includes a constant intrinsic drive']; end
title(tl, sup, 'FontWeight', 'normal', 'FontSize', 9);
end

% -------------------------------------------------------------------------
function iv = on_intervals(t, on)
% [start, end] times of each run of samples where on is true; a run ends one
% sample period after its last sample (the input is held over the sample).
iv = zeros(0, 2);
if ~any(on), return; end
dt = 0;
if numel(t) > 1, dt = t(2) - t(1); end
d = diff([false on false]);
t_start = t(d(1:end-1) == 1);            % first sample of each run
t_last  = t(d(2:end) == -1);             % last sample of each run
iv = [t_start(:), t_last(:) + dt];
iv(:, 2) = min(iv(:, 2), t(end));
end

% -------------------------------------------------------------------------
function [fig, info] = W_spectrum_figure(models, cond_names, cfg, st)
% Eigenvalues of the W that is simulated (params.W: level_of_chaos and any
% route scale folded in) for the first condition drawn, with the class's
% theoretical bulk radius R and outliers lambda_O. Those describe obj.W, the
% matrix before route scales, so they are marked only when the two agree.
m   = models{1};
W   = models{1}.get_params().W;
differs = cellfun(@(mk) ~isequal(mk.get_params().W, W), models(2:end));
lam = eig(full(W));
abscissa = max(real(lam));
radius   = max(abs(lam));
scaled   = ~isequal(W, m.W);
info = struct('eigenvalues', lam, 'abscissa', abscissa, 'radius', radius, ...
    'R', m.R, 'lambda_O', m.lambda_O, 'differs_across_conditions', any(differs), ...
    'route_scaled', scaled);

vis = 'on'; if ~cfg.visible, vis = 'off'; end
fig = figure('Color', 'w', 'Visible', vis, 'Position', [100 60 760 780]);
ax  = axes(fig); hold(ax, 'on');
h = plot(ax, real(lam), imag(lam), '.', 'Color', [0.15 0.15 0.15], 'MarkerSize', 7);
hs = h; ls = {sprintf('eig(W), n = %d', m.n)};
if ~scaled && ~m.rescale_by_abscissa
    th = linspace(0, 2 * pi, 400);
    hs(end+1) = plot(ax, m.R * cos(th), m.R * sin(th), '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1);
    ls{end+1} = sprintf('bulk radius R = %.3f', m.R);
    lo = m.lambda_O;
    lo(abs(lo) < 1e-9 * max(m.R, 1)) = 0;    % a balanced network's outliers are 0 to round-off
    hs(end+1) = plot(ax, real(lo), imag(lo), 'x', 'Color', [0.80 0.15 0.15], 'MarkerSize', 11, 'LineWidth', 1.5);
    ls{end+1} = ['outliers \lambda_O = ' strjoin(arrayfun(@(v) sprintf('%.3g', real(v)), lo(:)', ...
        'UniformOutput', false), ', ')];
end
hs(end+1) = xline(ax, abscissa, ':', 'Color', [0 0.45 0.70], 'LineWidth', 1.25);
ls{end+1} = sprintf('spectral abscissa %.3f', abscissa);
xline(ax, 0, '-', 'Color', [0.7 0.7 0.7]); yline(ax, 0, '-', 'Color', [0.7 0.7 0.7]);
hold(ax, 'off');
lim = 1.08 * max([radius, m.R, max(abs(m.lambda_O))]);
axis(ax, 'equal'); xlim(ax, [-lim lim]); ylim(ax, [-lim lim]);
grid(ax, 'on'); box(ax, 'off');
set(ax, 'FontSize', st.tick_fs - 3);
xlabel(ax, 'Re(\lambda)', 'FontSize', st.label_fs - 2);
ylabel(ax, 'Im(\lambda)', 'FontSize', st.label_fs - 2);
legend(ax, hs, ls, 'Location', 'southoutside', 'FontSize', 8, 'Box', 'off');

tl_txt = {sprintf('Eigenvalues of the simulated W   (%s, seeds %s)', ...
             strrep(cond_names{1}, '_', '\_'), mat2str(cfg.rng_seeds)), ...
         sprintf('spectral abscissa max Re \\lambda = %.4f,   spectral radius max |\\lambda| = %.4f', ...
             abscissa, radius), ...
         sprintf('level\\_of\\_chaos = %g', m.level_of_chaos)};
if scaled
    tl_txt{end+1} = 'route scale on: W differs from the drawn matrix; R and \lambda_O (of the drawn matrix) not shown';
end
if m.rescale_by_abscissa
    tl_txt{end+1} = 'rescale\_by\_abscissa on: R and \lambda_O do not apply; not shown';
end
if any(differs)
    tl_txt{end+1} = sprintf('simulated W DIFFERS across the conditions drawn (%s); this is %s', ...
        strrep(strjoin(cond_names([false differs]), ', '), '_', '\_'), strrep(cond_names{1}, '_', '\_'));
end
title(ax, tl_txt, 'FontWeight', 'normal', 'FontSize', 9);
end
