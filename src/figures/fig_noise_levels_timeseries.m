function out = fig_noise_levels_timeseries(cfg)
% FIG_NOISE_LEVELS_TIMESERIES Time series of every condition at several levels of one parameter.
%
%   out = FIG_NOISE_LEVELS_TIMESERIES()
%   out = FIG_NOISE_LEVELS_TIMESERIES('level_name', p, 'levels', v, 'level_label', s, ...
%                                     'model_args', {...}, 'out_dir', d)
%
% One network, simulated under each adaptation condition (the COLUMNS) at
% each of a few values of one scalar model property (the BANDS of rows, one
% band per value, labeled with it). Within a band the rows are, for a
% sample of neurons:
%
%   x                  the dendritic state
%   r                  the firing rate
%   synaptic output    r * (depression factor) * (facilitation), route 1
%   depression factor  SRNNCellTypePairs.std_route_factor of route 1:
%                      prod(b) on a product route, the releasable pool b_1
%                      on a chained route (the panel says which); "no
%                      depression" in a condition without it
%
% Each variable keeps one y axis across the bands, so the levels can be
% compared by eye. r, the synaptic output and the depression factor share
% one axis across the whole figure. x shares one across the bands of each
% column (cfg.x_scale_per_column, default true; false shares it across the
% columns too): on the mu7chain network x spans about twenty times the range
% without adaptation that it spans with it, and a common axis flattened the
% adapting columns into a line in the first draw.
%
% It was written for release-noise plan step P2 (docs/Plan_release_noise.md,
% section 3, "Two time-series figures"): the four conditions with the Wiener
% process at three levels, and the same with the release (gamma) process at
% three levels. Nothing here names either noise. The parameter swept across
% the bands is cfg.level_name (default 'sigma_u_noise'), its values
% cfg.levels and its label cfg.level_label; anything to hold fixed in every
% panel goes in cfg.model_args. For example:
%
%   fig_noise_levels_timeseries('level_name', 'sigma_u_noise', ...
%       'levels', [0 0.025 0.05], 'level_label', '\sigma_u')
%   fig_noise_levels_timeseries('level_name', 'release_sigma', ...
%       'levels', [0 0.007 0.014], 'level_label', '\sigma_\gamma', ...
%       'model_args', {'sigma_u_noise', 0, 'release_rate', 0.01})
%
% THE SAME NETWORK AND THE SAME SEEDS IN EVERY PANEL: rng_seeds (W, the
% setpoints, the per-neuron ladders and depression draws, x0) and the noise
% seed, which is set explicitly (cfg.noise_seed, default the model's own
% default rng_seeds(1) + 224737) so it cannot depend on anything else. So
% between bands only the level differs; for the Wiener noise the path is the
% same path, scaled. The integrator is the stochastic one the sweeps use for
% a noisy preset (analysis_run_config(...).sde_solver at cfg.run_mode, with
% that cell's fs), at every level including zero, so the bands differ only in
% the level and not in the scheme. The series are stored at the integration
% rate (plot_freq = fs), not the class's 10 Hz plot decimation, which would
% hide the noise.
%
% Each run starts at -cfg.T_settle and the figure shows [0, cfg.T_show].
% The neurons drawn are the first n_show_E of type 1 and the first n_show_I
% of type 2 (as fig_example_timeseries chooses them), in per-neuron shades of
% the type's color. The preset must carry the same synaptic dynamics on
% every route (SRNNCellTypePairs.routes_identical), since route 1 is drawn
% for every neuron.
%
% No run directory: the figure simulates (lya_method 'none'; a few seconds a
% panel at n = 500).
%
% out.figs, out.files, out.source as every figure; out.levels and
% out.cond_names; out.series{l, i} holds what was drawn in band l, column i
% (t, and for the shown neurons x, r, synaptic_output and depression, rows
% E first then I; depression is [] without depression); out.x_std(l, i) is
% the temporal standard deviation of x over the shown window, averaged over
% ALL neurons of the network, and out.r_mean(l, i) the mean rate.
%
% See also: fig_example_timeseries, fig_noise_switch_rollouts, build_from_preset,
%           analysis_run_config, SRNNCellTypePairs.std_route_factor

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.preset_name (1,:) char    = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain'
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''      % unused; accepted for a uniform call
    cfg.conditions                = {}      % names to draw, in order; empty = all the preset holds
    cfg.level_name  (1,:) char    = 'sigma_u_noise'   % scalar model property (or alias) varied across the bands
    cfg.levels      (1,:) double  = [0 0.025 0.05]
    cfg.level_label (1,:) char    = ''      % TeX label of the level; '' -> the property name
    cfg.model_args  (1,:) cell    = {}      % name-value pairs applied to every panel, before the level
    cfg.rng_seeds   (1,2) double  = [1 2]
    cfg.noise_seed  double        = []      % [] -> rng_seeds(1) + 224737, the model's default
    cfg.run_mode    (1,:) char    = 'fast'  % analysis_run_config cell giving the integrator and fs
    cfg.T_settle    (1,1) double {mustBePositive} = 10
    cfg.T_show      (1,1) double {mustBePositive} = 5
    cfg.n_show_E    (1,1) double {mustBeInteger, mustBeNonnegative} = 10
    cfg.n_show_I    (1,1) double {mustBeInteger, mustBeNonnegative} = 10
    cfg.x_scale_per_column (1,1) logical = true   % x shares its axis across bands, not columns (header)
    cfg.tag         (1,:) char    = 'Fig_noise_levels_timeseries'
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();
if isempty(cfg.levels)
    error('fig_noise_levels_timeseries:NoLevels', 'levels is empty.');
end
level_label = cfg.level_label;
if isempty(level_label), level_label = strrep(cfg.level_name, '_', '\_'); end
noise_seed = cfg.noise_seed;
if isempty(noise_seed), noise_seed = cfg.rng_seeds(1) + 224737; end

[preset, ~, conditions] = srnn_param_preset(cfg.preset_name);
cond_names = select_conditions(cellfun(@(c) c.name, conditions, 'UniformOutput', false), cfg.conditions);
n_cond = numel(cond_names);
n_lev  = numel(cfg.levels);

% The sweeps' stochastic scheme and its sample rate, at every level.
rc = analysis_run_config('sensitivity', cfg.run_mode, preset);
solver = rc.sde_solver;
fs     = rc.model.fs;

%% Simulate: one network, every condition, every level
series = cell(n_lev, n_cond);
x_std  = nan(n_lev, n_cond);
r_mean = nan(n_lev, n_cond);
labels = cell(1, n_cond);                % depression-factor note per column
t0 = tic;
for l = 1:n_lev
    for i = 1:n_cond
        m = build_from_preset(cfg.preset_name, cond_names{i}, ...
            'rng_seeds', cfg.rng_seeds, 'noise_seed', noise_seed, ...
            'T_range', [-cfg.T_settle, cfg.T_show], 'fs', fs, 'ode_solver', solver, ...
            'plot_freq', fs, 'lya_method', 'none', 'verbose', cfg.verbose, ...
            cfg.model_args{:}, cfg.level_name, cfg.levels(l));
        m.run();
        [series{l, i}, x_std(l, i), r_mean(l, i), labels{i}] = ...
            sampled_series(m, cfg, cond_names{i});
        vprintf(cfg.verbose, 'verbose', '  [noise_levels] %s = %g  %-14s std(x) = %.4f  mean r = %.4f\n', ...
            cfg.level_name, cfg.levels(l), cond_names{i}, x_std(l, i), r_mean(l, i));
    end
end
vprintf(cfg.verbose, 'minimal', '  [noise_levels] %d conditions x %d levels of %s simulated in %.0f s\n', ...
    n_cond, n_lev, cfg.level_name, toc(t0));

%% Figure: an outer layout of bands, each band a rows-by-conditions layout
var_fields = {'x', 'r', 'synaptic_output', 'depression'};
var_labels = {'x', 'r', 'synaptic output', 'depression factor'};
n_var = numel(var_fields);
vis = 'on'; if ~cfg.visible, vis = 'off'; end
fig = figure('Color', 'w', 'Visible', vis, ...
    'Position', [40 20 330 * n_cond + 150, 112 * n_var * n_lev + 130]);
outer = tiledlayout(fig, n_lev, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax_all = gobjects(n_lev, n_var, n_cond);
for l = 1:n_lev
    band = tiledlayout(outer, n_var, n_cond, 'TileSpacing', 'tight', 'Padding', 'tight');
    band.Layout.Tile = l;
    title(band, sprintf('%s = %g', level_label, cfg.levels(l)), ...
        'FontWeight', 'bold', 'FontSize', st.title_fs - 2);
    for v = 1:n_var
        for i = 1:n_cond
            ax = nexttile(band, (v - 1) * n_cond + i); ax_all(l, v, i) = ax;
            draw_panel(ax, series{l, i}, var_fields{v});
            set(ax, 'FontSize', st.tick_fs - 5); box(ax, 'off');
            xlim(ax, [0 cfg.T_show]);
            if i == 1, ylabel(ax, var_labels{v}, 'FontSize', st.label_fs - 5); end
            if v < n_var || l < n_lev, set(ax, 'XTickLabel', []); end
            if v == n_var && l == n_lev, xlabel(ax, 'time (s)', 'FontSize', st.label_fs - 4); end
            if l == 1 && v == 1
                ttl = cond_names{i};
                if isKey(st.condition_title, ttl), ttl = st.condition_title(ttl); end
                title(ax, ttl, 'FontWeight', 'normal', 'FontSize', st.title_fs - 3, ...
                    'Color', st.condition_color(cond_names{i}));
            end
            if v == n_var && ~isempty(labels{i}) && l == 1
                % which factor this is, just above the panel's top-right corner
                text(ax, 1, 1.0, labels{i}, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'bottom', 'FontSize', 7, 'Color', [0.3 0.3 0.3], 'Clipping', 'off');
            end
        end
    end
end

% Shared y axes (header): r, the synaptic output and the depression factor
% across the whole figure; x across the bands of each column, or of all
% columns when cfg.x_scale_per_column is false.
for v = 1:n_var
    if strcmp(var_fields{v}, 'x') && cfg.x_scale_per_column
        groups = arrayfun(@(i) reshape(ax_all(:, v, i), [], 1), 1:n_cond, 'UniformOutput', false);
    else
        groups = {reshape(ax_all(:, v, :), [], 1)};
    end
    for g = 1:numel(groups)
        share_ylim(groups{g}, strcmp(var_fields{v}, 'depression'));
    end
end

% Legend for the cell types, under the last band.
ax = ax_all(n_lev, 1, 1); hold(ax, 'on');
tcol = SRNNCellTypePairs.type_colors(2);
h = gobjects(1, 2);
for q = 1:2
    h(q) = plot(ax, NaN, NaN, '-', 'Color', tcol(q, :), 'LineWidth', 1.5);
end
hold(ax, 'off');
lg = legend(ax, h, {sprintf('E (%d shown)', series{1}.n_shown(1)), ...
    sprintf('I (%d shown)', series{1}.n_shown(2))}, 'Orientation', 'horizontal', ...
    'FontSize', 8, 'Box', 'off');
lg.Layout.Tile = 'south';
title(outer, {sprintf('One network at %d levels of %s, every condition on the same seeds', n_lev, level_label), ...
    sprintf('%s   (n = %d, seeds %s, noise seed %d, %s at %g Hz, after %g s settling)', ...
    strrep(cfg.preset_name, '_', '\_'), series{1}.n, mat2str(cfg.rng_seeds), noise_seed, ...
    solver, fs, cfg.T_settle)}, 'FontWeight', 'normal', 'FontSize', 9);

%% Save
out = struct('figs', fig, 'files', {{}}, 'source', 'simulated inline', ...
    'levels', cfg.levels, 'level_name', cfg.level_name, 'cond_names', {cond_names}, ...
    'series', {series}, 'x_std', x_std, 'r_mean', r_mean, 'ode_solver', solver, 'fs', fs);
if cfg.save
    save_figure_stable(out_dir, cfg.tag, fig);
    out.files = existing_outputs(out_dir, cfg.tag);
end
end

% =========================================================================
function [S, x_std, r_mean, note] = sampled_series(m, cfg, cond_name)
% The shown window's series for the sampled neurons, from the model's own
% plot data and accessors (std_route_factor for the depression factor), plus
% two whole-network summaries.
pd     = m.plot_data;
params = m.get_params();
names  = m.cell_type_names;
C      = numel(names);
if ~SRNNCellTypePairs.routes_identical(params)
    error('fig_noise_levels_timeseries:RoutesDiffer', ...
        'Preset ''%s'' (%s) carries different synaptic dynamics on different routes; this figure draws route 1 per neuron and needs them identical.', ...
        cfg.preset_name, cond_name);
end
keep = pd.t >= -1e-9;
S = struct('t', pd.t(keep), 'n', m.n, 'n_shown', zeros(1, C), 'type', [], ...
    'x', [], 'r', [], 'synaptic_output', [], 'depression', []);
x_all = []; r_all = [];
note = '';
for q = 1:C
    want = cfg.n_show_I;
    if q == 1, want = cfg.n_show_E; end
    k = 1:min(want, size(pd.x.(names{q}), 1));
    S.n_shown(q) = numel(k);
    S.type = [S.type; q * ones(numel(k), 1)];
    X = pd.x.(names{q})(:, keep);
    R = pd.r.(names{q})(:, keep);
    x_all = [x_all; X]; r_all = [r_all; R]; %#ok<AGROW>
    S.x = [S.x; X(k, :)];
    S.r = [S.r; R(k, :)];
    posts = fieldnames(pd.synaptic_output.(names{q}));
    SO = pd.synaptic_output.(names{q}).(posts{1});       % route 1 == every route
    S.synaptic_output = [S.synaptic_output; SO(k, keep)];
    B = pd.b.(names{q}).(posts{1});                       % n_q x n_b x n_t, or []
    if ~isempty(B)
        post_q = find(strcmp(names, posts{1}), 1);
        D = SRNNCellTypePairs.std_route_factor(B, params, q, post_q);   % whole route, then rows
        D = reshape(D, size(B, 1), size(B, 3));
        S.depression = [S.depression; D(k, keep)];
        if params.std_chain(q, post_q)
            note = 'b_1 (releasable pool)';
        elseif size(B, 2) > 1
            note = '\Pi_m b_m';
        else
            note = 'b (one pool)';
        end
    end
end
x_std  = mean(std(x_all, 0, 2));
r_mean = mean(r_all, 'all');
end

% =========================================================================
function draw_panel(ax, S, field)
% One variable of one band and condition: every shown neuron, in its type's
% shades, type 1 drawn last so it lies on top.
D = S.(field);
if isempty(D)
    text(ax, 0.5, 0.5, 'no depression', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.4 0.4 0.4]);
    set(ax, 'YTick', []);
    return;
end
C = max(S.type);
tcol = SRNNCellTypePairs.type_colors(C);
hold(ax, 'on');
for q = SRNNCellTypePairs.draw_order(C)
    rows = find(S.type == q);
    shades = SRNNCellTypePairs.neuron_colors(tcol(q, :), max(numel(rows), 1));
    for k = 1:numel(rows)
        plot(ax, S.t, D(rows(k), :), '-', 'Color', shades(k, :), 'LineWidth', 0.5);
    end
end
hold(ax, 'off');
end

% =========================================================================
function share_ylim(axs, from_zero)
% One y range for the axes given, from the lines they hold; a panel with no
% line (a condition without the variable) is left alone. The depression
% factor's range starts at 0 and reaches at least 1.
has = arrayfun(@(a) ~isempty(findobj(a, 'Type', 'line')), axs);
axs = axs(has);
if isempty(axs), return; end
lo = inf; hi = -inf;
for a = axs'
    L = findobj(a, 'Type', 'line');
    yd = [L.YData];
    lo = min(lo, min(yd)); hi = max(hi, max(yd));
end
if from_zero, lo = 0; hi = max(hi, 1); end
pad = 0.04 * max(hi - lo, eps);
set(axs, 'YLim', [lo - pad, hi + pad]);
end
