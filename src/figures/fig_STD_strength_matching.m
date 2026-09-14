function out = fig_STD_strength_matching(cfg)
% FIG_STD_STRENGTH_MATCHING One- vs two-timescale STD, strength-matched at r_ref.
%
%   out = FIG_STD_STRENGTH_MATCHING('run_dir', d)
%
% The paper's single-timescale and multiple-timescale regimes differ in STD
% timescale COUNT, and until 2026-09-13 also in STD STRENGTH: both depression
% factors had rho = tau_rel/tau_rec = 0.125, so the two-timescale steady state
% 1/(1 + r/rho)^2 was the square of the one-timescale one. The matched presets
% (srnn_param_preset, dualStdScaled / dualStdUsage) make the steady-state
% synaptic output of the two-timescale routes equal to the one-timescale
% routes' at a reference rate r_ref = the occupied median rate of the
% multiple-timescale condition (0.25). This figure shows the match over the
% rates the network actually uses, so it is visible rather than asserted.
%
% Panels:
%   A. theta_ss(r) = s * r * prod_m b_m(r) against r for the single-timescale
%      route (reference), the unmatched dual route, the scaled dual route
%      (s = 1 + r_ref/rho, the primary variant) and the usage-matched dual
%      route (rho_u on both timescales, the control), with r_ref marked and
%      the occupied 5th-95th percentile band of the multiple-timescale rates.
%   B. The ratio theta_dual / theta_single against r, log y: exactly 1 at r_ref
%      for both matchings; the departure away from r_ref is what each matching
%      costs. The scaled variant's ratio at r -> 0 is its low-rate gain factor.
%   C. The occupied rates: per-network mean rate at the default point of the
%      run's 1-D sweeps (all three conditions, pooled across sweeps and reps)
%      and, when cfg.simulate is on, the per-neuron time-averaged rates of one
%      20-s run of each condition of cfg.preset_name.
%
% A markdown table beside the figure records theta_ss at r_ref and at the
% occupied percentiles, the low-rate gain ratio and the largest relative
% mismatch over the occupied band, for each variant.
%
% The three presets are options so a different reference/primary/control
% triple can be drawn; the curves come from the E->E route of each preset's
% conditions (all four routes are identical in the paper's presets).
%
% See also: fig_STD_steady_state, srnn_param_preset, test_route_scale,
%           docs/notes/STD_strength_matching_2026-09-13.md

arguments
    cfg.verbose          (1,:) char    = 'minimal'
    cfg.preset_name      (1,:) char    = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStdScaled_3cond_mu8p25'
    cfg.reference_preset (1,:) char    = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25'
    cfg.control_preset   (1,:) char    = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStdUsage_3cond_mu8p25'
    cfg.r_ref            (1,1) double  = 0.25
    cfg.run_dir          (1,:) char    = ''      % occupied rates from its 1D_sensitivity_* sweeps; '' skips
    cfg.simulate         (1,1) logical = true    % one 20-s run per condition for per-neuron rates
    cfg.rng_seeds        (1,2) double  = [1 2]
    cfg.out_dir          (1,:) char    = ''
    cfg.save             (1,1) logical = true
    cfg.visible          (1,1) logical = true
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st      = manuscript_style();

%% The four routes
single_route = std_route(cfg.reference_preset, 'sfa1_std1');
dual_unm     = std_route(cfg.reference_preset, 'sfa3_std2');
dual_scaled  = std_route(cfg.preset_name,      'sfa3_std2');
dual_usage   = std_route(cfg.control_preset,   'sfa3_std2');
V = {single_route, dual_unm, dual_scaled, dual_usage};
v_label = {'single timescale (reference)', 'dual, unmatched', ...
    sprintf('dual, route scale %.3g (primary)', dual_scaled.scale), ...
    sprintf('dual, usage \\rho_u = %.3g (control)', dual_usage.rho(1))};
v_style = {'-', ':', '-', '--'};
v_col   = {[0 0 0], [0.5 0.5 0.5], st.condition_color('sfa3_std2'), 0.6 * st.condition_color('sfa3_std2')};

r = linspace(0, 1, 2000);
theta = cellfun(@(v) theta_ss(v, r), V, 'UniformOutput', false);

%% Occupied rates
occ = struct('cond', {{}}, 'rates', {{}});
cond_names = {'no_adaptation', 'sfa1_std1', 'sfa3_std2'};
if ~isempty(cfg.run_dir)
    occ = occupied_rates(cfg.run_dir, cond_names, cfg.verbose);
end
band = [NaN NaN];
if ~isempty(occ.rates)
    v = occ.rates{strcmp(occ.cond, 'sfa3_std2')};
    if ~isempty(v), band = prctile(v, [5 95]); end
end
neuron_rates = struct();
if cfg.simulate
    for c = cond_names
        m = build_from_preset(cfg.preset_name, c{1}, 'rng_seeds', cfg.rng_seeds, ...
            'T_range', [0 20], 'fs', 400, 'lya_method', 'none', 'verbose', cfg.verbose);
        m.run();
        pd = m.plot_data;
        keep = pd.t >= 10;
        rr = [];
        for q = fieldnames(pd.r)'
            rr = [rr; mean(pd.r.(q{1})(:, keep), 2)]; %#ok<AGROW>
        end
        neuron_rates.(c{1}) = rr;
    end
end

%% Figure
fig = figure('Color', 'w', 'Position', [80 80 1250 400]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl); hold(ax1, 'on');
if all(isfinite(band))
    fill(ax1, [band(1) band(2) band(2) band(1)], [0 0 1 1], st.condition_color('sfa3_std2'), ...
        'FaceAlpha', 0.08, 'EdgeColor', 'none');
end
h = gobjects(1, 4);
for k = 1:4
    h(k) = plot(ax1, r, theta{k}, v_style{k}, 'Color', v_col{k}, 'LineWidth', 1.6);
end
xline(ax1, cfg.r_ref, '-', sprintf('r_{ref} = %.2g', cfg.r_ref), 'Color', [0.3 0.3 0.3], ...
    'LabelOrientation', 'horizontal', 'FontSize', 9);
hold(ax1, 'off'); box(ax1, 'off');
ylim(ax1, [0, 1.1 * max(cellfun(@max, theta))]);
xlabel(ax1, 'presynaptic rate r', 'FontSize', st.label_fs);
ylabel(ax1, 'steady-state synaptic output  s r \Pi b_m(r)', 'FontSize', st.label_fs);
title(ax1, 'A  steady state', 'FontWeight', 'normal', 'FontSize', st.title_fs);
legend(ax1, h, v_label, 'Location', 'northwest', 'FontSize', 8, 'Box', 'off');
set(ax1, 'FontSize', st.tick_fs);

ax2 = nexttile(tl); hold(ax2, 'on');
if all(isfinite(band))
    fill(ax2, [band(1) band(2) band(2) band(1)], [1e-2 1e-2 1e2 1e2], st.condition_color('sfa3_std2'), ...
        'FaceAlpha', 0.08, 'EdgeColor', 'none');
end
for k = 2:4
    plot(ax2, r(2:end), theta{k}(2:end) ./ theta{1}(2:end), v_style{k}, 'Color', v_col{k}, 'LineWidth', 1.6);
end
yline(ax2, 1, ':', 'Color', [0.3 0.3 0.3]);
xline(ax2, cfg.r_ref, '-', 'Color', [0.3 0.3 0.3]);
hold(ax2, 'off'); box(ax2, 'off');
set(ax2, 'YScale', 'log', 'FontSize', st.tick_fs);
ylim(ax2, [0.05 5]);
xlabel(ax2, 'presynaptic rate r', 'FontSize', st.label_fs);
ylabel(ax2, '\theta_{dual} / \theta_{single}', 'FontSize', st.label_fs);
title(ax2, 'B  ratio to the single-timescale route', 'FontWeight', 'normal', 'FontSize', st.title_fs);

ax3 = nexttile(tl); hold(ax3, 'on');
edges = linspace(0, 1, 41);
hh = gobjects(0); lab = {};
for i = 1:numel(cond_names)
    col = st.condition_color(cond_names{i});
    if ~isempty(occ.rates)
        v = occ.rates{strcmp(occ.cond, cond_names{i})};
        if ~isempty(v)
            hh(end + 1) = histogram(ax3, v, edges, 'Normalization', 'pdf', 'FaceColor', col, ...
                'FaceAlpha', 0.35, 'EdgeColor', 'none'); %#ok<AGROW>
            lab{end + 1} = sprintf('%s: network mean (n = %d)', st.condition_short(cond_names{i}), numel(v)); %#ok<AGROW>
        end
    end
    if isfield(neuron_rates, cond_names{i})
        [f, xi] = ksdensity(neuron_rates.(cond_names{i}), linspace(0, 1, 200), 'Support', [-1e-3 1 + 1e-3]);
        hh(end + 1) = plot(ax3, xi, f, '-', 'Color', col, 'LineWidth', 1.4); %#ok<AGROW>
        lab{end + 1} = sprintf('%s: neurons, one run', st.condition_short(cond_names{i})); %#ok<AGROW>
    end
end
xline(ax3, cfg.r_ref, '-', 'Color', [0.3 0.3 0.3]);
hold(ax3, 'off'); box(ax3, 'off');
xlabel(ax3, 'mean firing rate', 'FontSize', st.label_fs);
ylabel(ax3, 'density', 'FontSize', st.label_fs);
title(ax3, 'C  occupied rates', 'FontWeight', 'normal', 'FontSize', st.title_fs);
ylim(ax3, [0 25]);   % the per-neuron densities spike at 0 and 1; the interior is the point
if ~isempty(hh); legend(ax3, hh, lab, 'Location', 'northeast', 'FontSize', 8, 'Box', 'off'); end
set(ax3, 'FontSize', st.tick_fs);
title(tl, sprintf('STD strength matching at r_{ref} = %.2g (curves: E->E route of each preset)', cfg.r_ref), ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Table
lo = band(1); hi = band(2);
if ~all(isfinite(band)); lo = 0.05; hi = 0.35; end
hdr = sprintf('| Variant | tau_rec | tau_rel | scale | theta(r_ref) | theta(p5 = %.3f) | theta(p95 = %.3f) | low-rate gain ratio | max |log ratio| over [p5, p95] |', lo, hi);
rows = cell(1, 4);
rr = linspace(lo, hi, 400);
for k = 1:4
    v = V{k};
    ratio = theta_ss(v, rr) ./ theta_ss(single_route, rr);
    rows{k} = sprintf('| %s | %s | %s | %.4g | %.4f | %.4f | %.4f | %.3g | %.3f |', strrep(v_label{k}, '\rho_u', 'rho_u'), ...
        mat2str(v.tau_rec, 4), mat2str(v.tau_rel, 4), v.scale, theta_ss(v, cfg.r_ref), ...
        theta_ss(v, lo), theta_ss(v, hi), v.scale, max(abs(log(ratio))));
end
occ_rows = {};
if ~isempty(occ.rates)
    for i = 1:numel(occ.cond)
        v = occ.rates{i};
        occ_rows{end + 1} = sprintf('| %s | %d | %.3f | %.3f | %.3f | %.3f | %.3f |', occ.cond{i}, numel(v), ...
            median(v), prctile(v, 5), prctile(v, 25), prctile(v, 75), prctile(v, 95)); %#ok<AGROW>
    end
end
vprintf(cfg.verbose, 'verbose', '%s\n|---|---|---|---|---|---|---|---|---|\n%s\n', hdr, strjoin(rows, newline));

if ~cfg.visible; set(fig, 'Visible', 'off'); end
fig_tag = 'Fig_STD_Strength_Matching';
out = struct('figs', fig, 'files', {{}}, 'source', 'analytic + one run per condition', ...
    'occupied', occ, 'band', band);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    fid = fopen(fullfile(out_dir, [fig_tag '_table.md']), 'w');
    if fid > 0
        fprintf(fid, '# STD strength matching at r_ref = %.3g\n\nPrimary preset: `%s`\nReference: `%s`\nControl: `%s`\n\n', ...
            cfg.r_ref, cfg.preset_name, cfg.reference_preset, cfg.control_preset);
        fprintf(fid, '%s\n|---|---|---|---|---|---|---|---|---|\n%s\n', hdr, strjoin(rows, newline));
        if ~isempty(occ_rows)
            fprintf(fid, '\nOccupied network mean rates at the default point of the 1-D sweeps in `%s`:\n\n', cfg.run_dir);
            fprintf(fid, '| Condition | n | median | p5 | p25 | p75 | p95 |\n|---|---|---|---|---|---|---|\n%s\n', strjoin(occ_rows, newline));
        end
        fclose(fid);
    end
end
end

%% ------------------------------------------------------------------------
function v = std_route(preset_name, cond_name)
% tau_rec, tau_rel, rho and scale of the E->E route of one condition.
[~, ~, conditions] = srnn_param_preset(preset_name);
names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
sc = conditions{strcmp(names, cond_name)}.synapse_config;
route = sc.E.E;
v = struct('preset', preset_name, 'cond', cond_name, 'tau_rec', route.std.tau_rec(:)', ...
    'tau_rel', route.std.tau_rel(:)', 'scale', 1);
v.rho = v.tau_rel ./ v.tau_rec;
if isfield(route, 'scale') && ~isempty(route.scale); v.scale = route.scale; end
end

function th = theta_ss(v, r)
th = v.scale .* r .* prod(1 ./ (1 + (1 ./ v.rho(:)) * r(:)'), 1);
th = reshape(th, size(r));
end

function occ = occupied_rates(run_dir, cond_names, verbose)
% Per-network mean rate at the level nearest the preset default of every
% 1D_sensitivity_* sweep in run_dir, pooled over sweeps and reps, per condition.
occ = struct('cond', {cond_names}, 'rates', {cell(size(cond_names))});
listing = dir(fullfile(run_dir, '1D_sensitivity_*'));
if isempty(listing)
    warning('fig_STD_strength_matching:NoSweeps', 'No 1D_sensitivity_* sweep in %s; panel C has no network rates.', run_dir);
    return;
end
for i = 1:numel(cond_names); occ.rates{i} = []; end
for k = 1:numel(listing)
    src = fullfile(listing(k).folder, listing(k).name);
    psa = ParamSpaceAnalysis2.from_dir(src);
    g = psa.grid_params{1};
    dv = preset_default_values(run_dir, {g});
    if ~isKey(dv, g)
        vprintf(verbose, 'verbose', '  [matching] %s: no default value, skipped\n', g);
        continue;
    end
    for i = 1:numel(cond_names)
        if ~isfield(psa.results, cond_names{i}); continue; end
        R = psa.results.(cond_names{i});
        ok = cellfun(@(x) isstruct(x) && isfield(x, 'success') && x.success, R);
        R = R(ok);
        vals = cellfun(@(x) psa.effective_param(x, g), R);
        mr = cellfun(@(x) x.mean_rate, R);
        lv = unique(vals);
        [~, j] = min(abs(lv - dv(g)));
        occ.rates{i} = [occ.rates{i}; reshape(mr(vals == lv(j)), [], 1)];
    end
end
end
