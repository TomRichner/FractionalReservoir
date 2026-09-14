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
%   (a) dendritic state x of a subset of neurons (n_show_E E, n_show_I I)
%   (b) firing rate r of the same neurons
%   (c) synaptic output r * prod(b) * prod(g) of the same neurons -- ONE
%       trace per neuron when every route carries the same synaptic dynamics
%       (SRNNCellTypePairs.routes_identical: route 1 is then the same array
%       as every other route), otherwise one trace per route with a legend
%   (d) population-mean SFA variable per timescale, per cell type
%   (e) depression prod(b) of the same neurons (route 1, same rule as (c))
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
% old single-column behaviour); empty draws all of them.
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
    cfg.n_show_E    (1,1) double  = 10
    cfg.n_show_I    (1,1) double  = 10
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''      % unused; accepted for a uniform call
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
end
n_cond = numel(cond_names);

%% Simulate: one network, every condition, the sweeps' Lyapunov estimator
T = cfg.T_range(2);
lya_T_interval = [min(5, T / 2), T];
models = cell(1, n_cond);
for i = 1:n_cond
    models{i} = build_from_preset(cfg.preset_name, cond_names{i}, ...
        'rng_seeds', cfg.rng_seeds, 'T_range', cfg.T_range, 'fs', 400, ...
        'lya_method', 'topk', 'lya_K', 15, 'lya_K_auto', false, 'lya_dt', 0.05, ...
        'lya_T_interval', lya_T_interval, 'lya_warmup', min(5, T / 4), ...
        'verbose', cfg.verbose);
    models{i}.run();
    vprintf(cfg.verbose, 'verbose', '  [example_timeseries] %-14s lambda_1 = %+.4f\n', ...
        cond_names{i}, models{i}.lya_results.LLE);
end

%% Figure
row_names = {'x', 'r', 'synaptic output', 'SFA  \langle a \rangle', 'depression  \Pi b', '\lambda'};
n_rows = numel(row_names);
fig = figure('Color', 'w', 'Position', [60 60 420 * n_cond + 80, 1000]);
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
    same_routes = SRNNCellTypePairs.routes_identical(params);
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
                plot(ax, t, D(k, :), '-', 'Color', shades{q}(k, :), 'LineWidth', 0.7);
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
    if same_routes
        for q = SRNNCellTypePairs.draw_order(C)
            posts = fieldnames(so.(names{q}));
            D = so.(names{q}).(posts{1});
            for k = 1:n_show(q)
                plot(ax, t, D(k, :), '-', 'Color', shades{q}(k, :), 'LineWidth', 0.7);
            end
        end
    else
        % per-route traces of the FIRST shown neuron of each type, with a legend
        hr = gobjects(0); lab = {};
        for q = 1:C
            posts = fieldnames(so.(names{q}));
            for p = 1:numel(posts)
                D = so.(names{q}).(posts{p});
                hr(end + 1) = plot(ax, t, D(1, :), '-', 'LineWidth', 0.8); %#ok<AGROW>
                lab{end + 1} = sprintf('%s\\rightarrow%s', names{q}, posts{p}); %#ok<AGROW>
            end
        end
        legend(ax, hr, lab, 'Location', 'northeast', 'FontSize', 8, 'Box', 'off');
        text(ax, 0.02, 0.95, 'routes differ: neuron 1 of each type, per route', ...
            'Units', 'normalized', 'FontSize', 8, 'VerticalAlignment', 'top');
    end
    hold(ax, 'off');

    % (d) SFA population means per timescale ---------------------------------
    ax = nexttile(tl); ax_all(4, i) = ax; hold(ax, 'on');
    ha = gobjects(0); lab = {};
    for q = 1:C
        A = pd.a.(names{q});                     % n_q x n_a x n_t
        if isempty(A), continue; end
        n_a = size(A, 2);
        tau = params.tau_a{q};
        for k = 1:n_a
            sh = 1 - 0.6 * (k - 1) / max(n_a - 1, 1);   % darker = slower
            colk = tcol(q, :) * sh;
            ha(end + 1) = plot(ax, t, reshape(mean(A(:, k, :), 1), 1, []), '-', 'Color', colk, 'LineWidth', 1.2); %#ok<AGROW>
            lab{end + 1} = sprintf('%s \\tau_a = %.3g s', names{q}, tau(k)); %#ok<AGROW>
        end
    end
    if isempty(ha)
        text(ax, 0.5, 0.5, 'no SFA', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 10);
        set(ax, 'YTick', []);
    else
        legend(ax, ha, lab, 'Location', 'northeast', 'FontSize', 7, 'Box', 'off', 'NumColumns', 2);
    end
    hold(ax, 'off');

    % (e) depression prod(b) ---------------------------------------------------
    ax = nexttile(tl); ax_all(5, i) = ax; hold(ax, 'on');
    drew = false;
    for q = SRNNCellTypePairs.draw_order(C)
        posts = fieldnames(pd.b.(names{q}));
        for p = 1:numel(posts)
            B = pd.b.(names{q}).(posts{p});      % n_q x n_b x n_t, or []
            if isempty(B), continue; end
            Pb = reshape(prod(B, 2), size(B, 1), size(B, 3));   % n_q x n_t
            for k = 1:n_show(q)
                plot(ax, t, Pb(k, :), '-', 'Color', shades{q}(k, :), 'LineWidth', 0.7);
            end
            drew = true;
            if same_routes, break; end           % route 1 is every route
        end
    end
    if ~drew
        text(ax, 0.5, 0.5, 'no STD', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 10);
        set(ax, 'YTick', []);
    else
        ylim(ax, [0 1.02]);
    end
    hold(ax, 'off');

    % (f) local and finite-time lambda_1 ---------------------------------------
    ax = nexttile(tl); ax_all(6, i) = ax;
    plot_local_and_finite_lle(ax, m.lya_results, cond_col);
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
    if n_cond > 1 && row ~= 4
        linkaxes(ax_all(row, :), 'y');
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
end
