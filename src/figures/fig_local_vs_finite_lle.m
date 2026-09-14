function out = fig_local_vs_finite_lle(cfg)
% FIG_LOCAL_VS_FINITE_LLE Local expansion vs long-interval stability, across trials.
%
%   out = FIG_LOCAL_VS_FINITE_LLE('run_dir', d)
%
% THE CLAIM THIS TESTS (Codex audit sec. 3 of fig_to_do_in_future.md): that
% multiple-timescale adaptation gives a regime that is slightly stable over an
% extended interval while still allowing transient epochs of local expansion.
% The evidence has to come from many trials, not one trace: the distribution
% of the LOCAL rate (the leading direction's instantaneous stretching rate,
% pooled over time and networks) against the distribution of the FINITE-TIME
% lambda_1 (one number per network, accumulated over the sweep's window). A
% regime that is "stable with transient expansion" has a finite-time
% distribution sitting at or below zero while a substantial share of the
% local-rate mass sits above it. "Transient expansion" is only ever said of
% the local positive intervals; "stable" / "chaotic" only of the finite-time
% window, whose length the sweep's lya_T_interval sets.
%
% No simulation is re-run. Two sets are read from the run directory:
%   near-default -- every 1D_sensitivity_* sweep at the level nearest the
%                   preset's default (preset_default_values), pooled over
%                   sweeps and reps: the reference network's own regime.
%   joint        -- the param_space_* sample, all grid points: the same
%                   measures across the whole explored parameter space.
% Per job the fields read are LLE, local_rate_lead (the stored leading
% local-rate series), frac_local_positive, mean_positive_excursion_s and
% p95_finite_0p2s (ParamSpaceAnalysis2.run_single_job, SRNNCellTypePairs.
% lya_summary). Failed jobs are skipped.
%
% Layout, ONE row, one column per condition (near-default set): density of
% the local rate pooled over time and jobs (light fill) with the finite-time
% lambda_1 across jobs on top (dark bars, median line); zero line; shared x.
% The text gives the share of local-rate samples > 0, the share of jobs with
% lambda_1 < 0 and the median lambda_1.
%
% A second row (boxcharts of frac_local_positive and
% mean_positive_excursion_s, near-default beside the joint sample) was
% removed on 2026-09-14 (TR): it was only about the single leading exponent,
% and a replacement built on the top-K local rates is future work. The joint
% sample is still read for the table.
%
% The table (<tag>_table.md, printed at 'verbose') gives per condition and set:
% n jobs, lambda_1 median [IQR], share lambda_1 < 0, pooled local-rate share
% > 0, and the medians of frac_local_positive, mean_positive_excursion_s and
% p95_finite_0p2s.
%
% See also: fig_lle_vs_rate, fig_sensitivity_medians, sweep_metrics,
%           SRNNCellTypePairs.lya_summary, preset_default_values

arguments
    cfg.verbose     (1,:) char    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (see verbose_level)
    cfg.data_file   (1,:) char    = ''          % unused; accepted for a uniform call
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = ''          % unused; the run records its preset
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();
if isempty(cfg.run_dir) || ~isfolder(cfg.run_dir)
    error('fig_local_vs_finite_lle:NoRunDir', ...
        'run_dir must name an existing run directory (got ''%s'').', cfg.run_dir);
end

%% Load
near  = load_near_default(cfg.run_dir, cfg.verbose);
joint = load_joint(cfg.run_dir, cfg.verbose);
cond_names = near.cond_names;
n_cond = numel(cond_names);
if n_cond == 0
    error('fig_local_vs_finite_lle:NoJobs', 'No successful jobs found under %s.', cfg.run_dir);
end

%% Figure
fig = figure('Color', 'w', 'Position', [80 80 380 * n_cond, 360]);
tl = tiledlayout(fig, 1, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');
ax_top = gobjects(1, n_cond);
rows = {};

% One x window for the top row: the local-rate mass is wide, the finite-time
% values narrow; the 1st-99th percentiles of the pooled local rates cover both.
all_local = [];
for i = 1:n_cond
    all_local = [all_local; near.local{i}(:)]; %#ok<AGROW>
end
x_lim = prctile(all_local(isfinite(all_local)), [1 99]);
x_lim = [min(x_lim(1), -0.5), max(x_lim(2), 0.5)];
edges = linspace(x_lim(1), x_lim(2), 61);

for i = 1:n_cond
    name = cond_names{i};
    col  = st.condition_color(name);
    light = 1 - 0.45 * (1 - col);
    loc  = near.local{i}; loc = loc(isfinite(loc));
    lle  = near.lle{i};

    ax = nexttile(tl); hold(ax, 'on');
    if ~isempty(loc)
        histogram(ax, loc, edges, 'Normalization', 'pdf', 'FaceColor', light, ...
            'EdgeColor', 'none', 'FaceAlpha', 0.9);
    end
    if ~isempty(lle)
        histogram(ax, lle, edges, 'Normalization', 'pdf', 'FaceColor', col, ...
            'EdgeColor', 'none', 'FaceAlpha', 0.75);
        xline(ax, median(lle), '-', 'Color', col, 'LineWidth', 2);
    end
    xline(ax, 0, ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.2);
    hold(ax, 'off'); box(ax, 'off');
    set(ax, 'FontSize', st.tick_fs, 'XLim', x_lim);
    title(ax, st.condition_title(name), 'FontWeight', 'normal', 'FontSize', st.title_fs);
    xlabel(ax, 'rate (1/s)', 'FontSize', st.label_fs);
    if i == 1; ylabel(ax, 'density', 'FontSize', st.label_fs); end
    txt = sprintf(['local > 0: %.0f%% of samples\n' ...
                   '\\lambda_1 < 0: %d of %d networks\n' ...
                   'median \\lambda_1 = %+.3f'], ...
        100 * mean(loc > 0), nnz(lle < 0), numel(lle), median(lle));
    text(ax, 0.97, 0.97, txt, 'Units', 'normalized', 'FontSize', 9, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    if i == n_cond   % the last column has room on the left; the first is full of mass
        legend(ax, {'local rate, pooled', 'finite-time \lambda_1 per network'}, ...
            'Location', 'northwest', 'FontSize', 8, 'Box', 'off');
    end
    ax_top(i) = ax;

    % Table rows for both sets (the joint sample is no longer drawn).
    sets = {near, joint};
    for s = 1:2
        j = find(strcmp(sets{s}.cond_names, name), 1);
        if isempty(j); continue; end
        rows{end + 1} = table_row(name, sets{s}.label, sets{s}, j); %#ok<AGROW>
    end
end
linkaxes(ax_top, 'y');
title(tl, sprintf('Local expansion vs finite-time \\lambda_1: %d near-default networks per condition, joint sample of %d', ...
    numel(near.lle{1}), numel(joint.lle{1})), 'FontWeight', 'normal', 'FontSize', 11);

hdr = '| Condition | set | n | lambda_1 median [IQR] | share lambda_1 < 0 | local rate share > 0 | median frac_local_positive | median mean excursion (s) | median p95 finite 0.2 s |';
sep = '|---|---|---|---|---|---|---|---|---|';
vprintf(cfg.verbose, 'verbose', '%s\n%s\n', hdr, sep);
vprintf(cfg.verbose, 'verbose', '%s\n', rows{:});

if ~cfg.visible; set(fig, 'Visible', 'off'); end

fig_tag = 'Fig_Local_vs_Finite_LLE';
out = struct('figs', fig, 'files', {{}}, 'source', cfg.run_dir);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    fid = fopen(fullfile(out_dir, [fig_tag '_table.md']), 'w');
    if fid > 0
        fprintf(fid, '# Local rate vs finite-time lambda_1\n\nRun: `%s`. Near-default = 1-D sweeps at the level nearest the preset default (%s); joint = param_space sample.\n\n', ...
            cfg.run_dir, strjoin(near.sweeps, ', '));
        fprintf(fid, '%s\n%s\n', hdr, sep);
        fprintf(fid, '%s\n', rows{:});
        fclose(fid);
    end
end
end

%% ------------------------------------------------------------------------
function s = table_row(name, label, S, j)
lle = S.lle{j}(:); loc = S.local{j}(:); loc = loc(isfinite(loc));
q = prctile(lle, [25 75]);
s = sprintf('| %s | %s | %d | %+.3f [%+.3f, %+.3f] | %.2f | %.2f | %.2f | %.2f | %+.2f |', ...
    name, label, numel(lle), median(lle), q(1), q(2), mean(lle < 0), mean(loc > 0), ...
    nanmedian_local(S.fpos{j}), nanmedian_local(S.exc{j}), nanmedian_local(S.p95{j}));
end

function m = nanmedian_local(v)
v = v(isfinite(v));
if isempty(v); m = NaN; else; m = median(v); end
end

function D = empty_set(label)
D = struct('label', label, 'cond_names', {{}}, 'lle', {{}}, 'local', {{}}, ...
    'fpos', {{}}, 'exc', {{}}, 'p95', {{}}, 'sweeps', {{}});
end

function D = add_jobs(D, cond, R)
% Append the successful jobs in R (a cell array of result structs) to set D.
ok = cellfun(@(r) isstruct(r) && isfield(r, 'success') && r.success && isfield(r, 'LLE'), R);
R = R(ok);
if isempty(R); return; end
j = find(strcmp(D.cond_names, cond), 1);
if isempty(j)
    D.cond_names{end + 1} = cond; j = numel(D.cond_names);
    D.lle{j} = []; D.local{j} = []; D.fpos{j} = []; D.exc{j} = []; D.p95{j} = [];
end
D.lle{j}  = [D.lle{j};  cellfun(@(r) r.LLE, R(:))];
D.fpos{j} = [D.fpos{j}; cellfun(@(r) field_or_nan(r, 'frac_local_positive'), R(:))];
D.exc{j}  = [D.exc{j};  cellfun(@(r) field_or_nan(r, 'mean_positive_excursion_s'), R(:))];
D.p95{j}  = [D.p95{j};  cellfun(@(r) field_or_nan(r, 'p95_finite_0p2s'), R(:))];
for k = 1:numel(R)
    if isfield(R{k}, 'local_rate_lead') && ~isempty(R{k}.local_rate_lead)
        D.local{j} = [D.local{j}; R{k}.local_rate_lead(:)];
    end
end
end

function v = field_or_nan(r, f)
if isfield(r, f) && ~isempty(r.(f)); v = r.(f); else; v = NaN; end
end

function D = load_near_default(run_dir, verbose)
% Every 1D_sensitivity_* sweep at the level nearest the preset default.
D = empty_set('near-default');
listing = dir(fullfile(run_dir, '1D_sensitivity_*'));
listing = listing([listing.isdir]);
if isempty(listing)
    error('fig_local_vs_finite_lle:NoSweeps', 'No 1D_sensitivity_* sweep in %s.', run_dir);
end
for k = 1:numel(listing)
    src = fullfile(listing(k).folder, listing(k).name);
    if ~isfile(fullfile(src, 'psa_object.mat')); continue; end
    psa = ParamSpaceAnalysis2.from_dir(src);
    swept = setdiff(psa.grid_params, {'reps'}, 'stable');
    if isempty(swept); continue; end
    g = swept{1};
    dv = preset_default_values(run_dir, {g});
    if ~isKey(dv, g)
        vprintf(verbose, 'verbose', '[local_vs_finite] %s: no default value, skipped\n', g);
        continue;
    end
    D.sweeps{end + 1} = g;
    cn = cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false);
    for c = 1:numel(cn)
        if ~isfield(psa.results, cn{c}); continue; end
        R = psa.results.(cn{c});
        ok = cellfun(@(r) isstruct(r) && isfield(r, 'success') && r.success, R);
        R = R(ok);
        if isempty(R); continue; end
        vals = cellfun(@(r) psa.effective_param(r, g), R);
        lv = unique(vals);
        [~, m] = min(abs(lv - dv(g)));
        D = add_jobs(D, cn{c}, R(vals == lv(m)));
    end
    vprintf(verbose, 'verbose', '[local_vs_finite] %-16s nearest level to default %.4g used\n', g, dv(g));
end
end

function D = load_joint(run_dir, verbose)
D = empty_set('joint');
ps = dir(fullfile(run_dir, 'param_space_*'));
ps = ps([ps.isdir]);
if isempty(ps)
    vprintf(verbose, 'minimal', '[local_vs_finite] no param_space_* in %s; joint set empty\n', run_dir);
    return;
end
psa = ParamSpaceAnalysis2.from_dir(fullfile(ps(1).folder, ps(1).name));
cn = cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false);
for c = 1:numel(cn)
    if isfield(psa.results, cn{c}); D = add_jobs(D, cn{c}, psa.results.(cn{c})); end
end
end
