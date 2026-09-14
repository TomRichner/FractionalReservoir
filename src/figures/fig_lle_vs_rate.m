function out = fig_lle_vs_rate(cfg)
% FIG_LLE_VS_RATE lambda_1 against mean firing rate, per condition, coloured by E:I weight balance.
%
%   out = FIG_LLE_VS_RATE('run_dir', d)
%
% THE CLAIM THIS TESTS (Codex audit sec. 4 of fig_to_do_in_future.md): mean
% firing rate and dynamical stability are related through the nonlinear
% operating point but are NOT directly correlated or interchangeable -- a quiet
% or saturated network can be stable, and networks with similar rates can have
% different lambda_1. "Merely increasing inhibition" is therefore an incomplete
% description of what changes the dynamics. The figure puts the two measures
% against each other, per adaptation condition, for every network of the joint
% param_space_* sample (and, with include_sweeps, of the 1-D sweeps), coloured
% by the realised E:I WEIGHT balance so the inhibition axis is visible too.
%
% Association is reported as Spearman's rank correlation with a bootstrap 95%
% interval, per condition. PEARSON IS DELIBERATELY NOT USED: the relation is
% visibly nonlinear (rate saturates at 1 and floors at 0 while lambda_1 does
% neither) and often multimodal, so a linear coefficient would misstate it. The
% quiet (mean rate < 0.02) and saturated (> 0.9) subgroups are outlined, and
% the table gives lambda_1's median inside each, which is where "a quiet or
% saturated network can be stable" is read off.
%
% Per job the fields read are LLE, mean_rate, config_idx and network_seed
% (the last two to rebuild the network for its weight balance, cached per grid
% point exactly as fig_EI_weights_param_space does). Failed jobs are skipped.
% Points with |lambda_1| > 5 are clamped to the axis edge and drawn as
% triangles so the range stays readable.
%
% See also: fig_EI_weights_param_space, fig_local_vs_finite_lle,
%           sweep_metrics, ParamSpaceAnalysis2.rebuild_model

arguments
    cfg.verbose        (1,:) char    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (see verbose_level)
    cfg.data_file      (1,:) char    = ''          % unused; accepted for a uniform call
    cfg.out_dir        (1,:) char    = ''
    cfg.save           (1,1) logical = true
    cfg.visible        (1,1) logical = true
    cfg.run_dir        (1,:) char    = ''
    cfg.preset_name    (1,:) char    = ''          % unused; the run records its preset
    cfg.include_sweeps (1,1) logical = true        % pool the 1D_sensitivity_* jobs too (smaller markers)
    cfg.n_boot         (1,1) double  = 1000
    cfg.quiet_rate     (1,1) double  = 0.02
    cfg.saturated_rate (1,1) double  = 0.9
    cfg.lle_clamp      (1,1) double  = 5
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();
if isempty(cfg.run_dir) || ~isfolder(cfg.run_dir)
    error('fig_lle_vs_rate:NoRunDir', ...
        'run_dir must name an existing run directory (got ''%s'').', cfg.run_dir);
end

% The colour axis: the same fixed E:I weight window fig_EI_weights_param_space
% pins (1:10 to 10:1), so the two sheets read on one scale.
W_CLIM  = [1/11, 10/11];
cb_ticks  = [1/11, 0.2, 1/3, 0.5, 2/3, 0.8, 10/11];
cb_labels = {'1:10', '1:4', '1:2', '1:1', '2:1', '4:1', '10:1'};

%% Load: joint sample (+ sweeps)
lle = struct(); rate = struct(); be = struct(); src = struct();
sweep_names = {};
ps = dir(fullfile(cfg.run_dir, 'param_space_*'));
ps = ps([ps.isdir]);
if isempty(ps)
    error('fig_lle_vs_rate:NoParamSpace', 'No param_space_* subdir found in %s.', cfg.run_dir);
end
dirs = {fullfile(ps(1).folder, ps(1).name)};
src_id = 1;
if cfg.include_sweeps
    sl = dir(fullfile(cfg.run_dir, '1D_sensitivity_*'));
    sl = sl([sl.isdir]);
    for k = 1:numel(sl); dirs{end + 1} = fullfile(sl(k).folder, sl(k).name); end %#ok<AGROW>
end
t_load = tic;
n_rebuilt = 0;
for d = 1:numel(dirs)
    if ~isfile(fullfile(dirs{d}, 'psa_object.mat')); continue; end
    psa = ParamSpaceAnalysis2.from_dir(dirs{d});
    if d > 1; sweep_names{end + 1} = strjoin(setdiff(psa.grid_params, {'reps'}, 'stable'), ','); end %#ok<AGROW>
    cache = containers.Map('KeyType', 'double', 'ValueType', 'double');   % per PSA: config_idx is per grid
    cn = cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false);
    for c = 1:numel(cn)
        name = cn{c};
        if ~isfield(psa.results, name); continue; end
        R = psa.results.(name);
        ok = cellfun(@(r) isstruct(r) && isfield(r, 'success') && r.success && isfield(r, 'LLE') && isfinite(r.LLE), R);
        R = R(ok);
        if isempty(R); continue; end
        if ~isfield(lle, name); lle.(name) = []; rate.(name) = []; be.(name) = []; src.(name) = []; end
        lle.(name)  = [lle.(name);  cellfun(@(r) r.LLE, R(:))];
        rate.(name) = [rate.(name); cellfun(@(r) r.mean_rate, R(:))];
        b = nan(numel(R), 1);
        for k = 1:numel(R)
            b(k) = ei_weight_fraction(psa, R{k}, cache, cfg.verbose);
        end
        be.(name)  = [be.(name); b];
        src.(name) = [src.(name); src_id * ones(numel(R), 1)];
    end
    n_rebuilt = n_rebuilt + cache.Count;
    src_id = src_id + 1;
end
cond_names = fieldnames(lle)';
n_cond = numel(cond_names);
if n_cond == 0
    error('fig_lle_vs_rate:NoJobs', 'No successful jobs found under %s.', cfg.run_dir);
end
vprintf(cfg.verbose, 'verbose', '[lle_vs_rate] %d networks rebuilt for the weight balance in %.1f s\n', n_rebuilt, toc(t_load));

%% Figure
fig = figure('Color', 'w', 'Position', [80 80 400 * n_cond + 80, 440]);
tl = tiledlayout(fig, 1, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
rows = cell(1, n_cond);
axs = gobjects(1, n_cond);
y_all = [];
for i = 1:n_cond; y_all = [y_all; lle.(cond_names{i})]; end %#ok<AGROW>
y_lim = [max(min(y_all), -cfg.lle_clamp), min(max(y_all), cfg.lle_clamp)];
y_lim = y_lim + 0.05 * diff(y_lim) * [-1 1];
y_lim(1) = min(y_lim(1), -0.25); y_lim(2) = max(y_lim(2), 0.25);
bin_edges = linspace(0, 1, 11);

for i = 1:n_cond
    name = cond_names{i};
    y = lle.(name); x = rate.(name); c = be.(name); s = src.(name);
    ax = nexttile(tl); hold(ax, 'on');
    clamped = abs(y) > cfg.lle_clamp;
    yc = max(min(y, cfg.lle_clamp), -cfg.lle_clamp);
    quiet = x < cfg.quiet_rate; sat = x > cfg.saturated_rate;
    sz = 30 * (s == 1) + 12 * (s > 1);
    % bulk
    plain = ~clamped & ~quiet & ~sat;
    scatter(ax, x(plain), yc(plain), sz(plain), c(plain), 'filled', 'MarkerFaceAlpha', 0.75);
    % subgroups outlined
    edge = quiet | sat;
    scatter(ax, x(edge & ~clamped), yc(edge & ~clamped), sz(edge & ~clamped) + 10, c(edge & ~clamped), ...
        'filled', 'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 1);
    % overflow
    up = clamped & y > 0; dn = clamped & y < 0;
    scatter(ax, x(up), yc(up), 40, c(up), '^', 'filled', 'MarkerEdgeColor', [0.1 0.1 0.1]);
    scatter(ax, x(dn), yc(dn), 40, c(dn), 'v', 'filled', 'MarkerEdgeColor', [0.1 0.1 0.1]);
    % running median in rate bins
    bm = nan(1, numel(bin_edges) - 1); bx = bm;
    for b = 1:numel(bm)
        in = x >= bin_edges(b) & x < bin_edges(b + 1);
        if b == numel(bm); in = in | x == 1; end
        if nnz(in) >= 3; bm(b) = median(yc(in)); bx(b) = median(x(in)); end
    end
    plot(ax, bx(~isnan(bm)), bm(~isnan(bm)), '-', 'Color', [0.15 0.15 0.15], 'LineWidth', 2.2);
    yline(ax, 0, ':', 'Color', [0 0.6 0], 'LineWidth', 1.2);
    xline(ax, cfg.quiet_rate, ':', 'Color', [0.5 0.5 0.5]);
    xline(ax, cfg.saturated_rate, ':', 'Color', [0.5 0.5 0.5]);
    hold(ax, 'off'); box(ax, 'off');
    colormap(ax, ei_colormap());
    set(ax, 'CLim', W_CLIM, 'FontSize', st.tick_fs, 'XLim', [-0.02 1.02], 'YLim', y_lim);
    [rho, ci] = spearman_boot(x, y, cfg.n_boot);
    title(ax, sprintf('%s\n\\rho_S = %+.2f [%+.2f, %+.2f], n = %d', st.condition_title(name), rho, ci(1), ci(2), numel(y)), ...
        'FontWeight', 'normal', 'FontSize', st.title_fs - 2);
    xlabel(ax, st.label_rate, 'FontSize', st.label_fs);
    if i == 1; ylabel(ax, [st.label_lle ' (1/s)'], 'FontSize', st.label_fs); end
    axs(i) = ax;
    rows{i} = sprintf('| %s | %+.2f [%+.2f, %+.2f] | %d | %d | %d | %s | %s | %s |', name, rho, ci(1), ci(2), ...
        numel(y), nnz(quiet), nnz(sat), med_txt(y(quiet)), med_txt(y(~quiet & ~sat)), med_txt(y(sat)));
end
cb = colorbar(axs(end));
cb.Ticks = cb_ticks; cb.TickLabels = cb_labels;
cb.Label.String = 'E:I weight ratio'; cb.Label.FontSize = st.label_fs;
cb.Layout.Tile = 'east';
title(tl, sprintf('\\lambda_1 vs mean rate%s. Line: median per rate bin. Outlined: quiet (< %.2g) and saturated (> %.2g).', ...
    tern(cfg.include_sweeps, ', joint sample + 1-D sweeps (small)', ', joint sample'), cfg.quiet_rate, cfg.saturated_rate), ...
    'FontWeight', 'normal', 'FontSize', 9);

hdr = '| Condition | Spearman rho [95% CI] | n | n quiet | n saturated | lambda_1 median quiet | mid | saturated |';
sep = '|---|---|---|---|---|---|---|---|';
vprintf(cfg.verbose, 'verbose', '%s\n%s\n', hdr, sep);
vprintf(cfg.verbose, 'verbose', '%s\n', rows{:});

if ~cfg.visible; set(fig, 'Visible', 'off'); end

fig_tag = 'Fig_LLE_vs_Rate';
out = struct('figs', fig, 'files', {{}}, 'source', cfg.run_dir);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    fid = fopen(fullfile(out_dir, [fig_tag '_table.md']), 'w');
    if fid > 0
        fprintf(fid, '# lambda_1 vs mean firing rate\n\nRun: `%s`. Sources: param_space sample%s. Spearman rank correlation with a %d-sample bootstrap 95%% CI; Pearson deliberately not used (nonlinear, saturating relation). Quiet: mean rate < %.2g; saturated: > %.2g.\n\n', ...
            cfg.run_dir, tern(cfg.include_sweeps, [' + 1-D sweeps (' strjoin(sweep_names, ', ') ')'], ''), ...
            cfg.n_boot, cfg.quiet_rate, cfg.saturated_rate);
        fprintf(fid, '%s\n%s\n', hdr, sep);
        fprintf(fid, '%s\n', rows{:});
        fclose(fid);
    end
end
end

%% ------------------------------------------------------------------------
function [rho, ci] = spearman_boot(x, y, n_boot)
n = numel(x);
if n < 3; rho = NaN; ci = [NaN NaN]; return; end
rho = corr(x(:), y(:), 'type', 'Spearman');
rng(0, 'twister');
bs = nan(n_boot, 1);
for b = 1:n_boot
    idx = randi(n, n, 1);
    if numel(unique(x(idx))) < 2 || numel(unique(y(idx))) < 2; continue; end
    bs(b) = corr(x(idx), y(idx), 'type', 'Spearman');
end
ci = prctile(bs(isfinite(bs)), [2.5 97.5]);
end

function s = med_txt(v)
if isempty(v); s = '-'; else; s = sprintf('%+.3f (n=%d)', median(v), numel(v)); end
end

function s = tern(c, a, b)
if c; s = a; else; s = b; end
end

function cm = ei_colormap()
% A diverging blue (inhibition-dominant) to red (excitation-dominant) map,
% white at 1:1; 256 rows.
n = 256;
top = [0.85 0.15 0.15]; mid = [0.96 0.96 0.96]; bot = [0.15 0.35 0.85];
h = n / 2;
cm = [interp1([0 1], [bot; mid], linspace(0, 1, h)'); interp1([0 1], [mid; top], linspace(0, 1, n - h)')];
end

function v = ei_weight_fraction(psa, res, cache, verbose)
% Copied from fig_EI_weights_param_space: the realised E:I weight balance of
% the network at this grid point, cached per config_idx (the network is shared
% by every condition run at a point). |sum W_E| / (|sum W_E| + |sum W_I|).
key = res.config_idx;
if isKey(cache, key)
    v = cache(key);
    return
end
m = psa.rebuild_model(res);
m.verbose = verbose;   % the model prints its build report only at 'verbose'
m.build();
ti = m.type_indices;
assert(numel(ti) >= 2, ...
    'ei_weight_fraction needs at least two cell types; got %d.', numel(ti));
S_E = full(sum(sum(m.W(:, ti{1}))));
S_I = full(sum(sum(m.W(:, ti{2}))));
denom = abs(S_E) + abs(S_I);
if denom == 0
    v = 0.5;
else
    v = abs(S_E) / denom;
end
cache(key) = v; %#ok<NASGU>  handle object: this mutates the caller's map
end
