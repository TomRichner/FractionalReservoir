function out = fig_sfa_EOC_allStd(cfg)
% FIG_SFA_EOC_ALLSTD SFA edge of chaos: lambda_1 against the slowest tau_a.
%
%   out = FIG_SFA_EOC_ALLSTD()
%   out = FIG_SFA_EOC_ALLSTD('run_dir', d)
%
% Replots the tau-sensitivity LLE panel -- how the largest Lyapunov exponent
% approaches 0 (the edge of chaos) as the slowest SFA adaptation timescale grows.
% No simulation is re-run: it reloads the saved tau sweep's PSA object (axis
% tau_a_EI, E and I together, since 2026-09-14; tau_a_E, E only, in older runs)
% and draws the single condition (the full-adaptation regime).
%
% THE WHOLE DISTRIBUTION IS SHOWN (2026-09-14). Per level of the sweep the
% reps' lambda_1 values are drawn as a vertical density strip (a 2-D
% histogram, level x lambda_1, white -> dark grey), with the median as a thick
% line inside a 25th-75th percentile band, the mean as a marker, and a black
% lambda_1 = 0 line. The y range is the 1st-99th percentile of ALL values,
% padded 10% and always including 0; values outside it are counted into two
% explicit OVERFLOW rows drawn inside the axes at the top and bottom and
% labelled with their counts, so nothing is dropped silently. Along the top,
% each level carries n reps and the share of reps with lambda_1 > 0.
%
% Until 2026-09-14 this figure restyled plot_sensitivity's panel with a fixed
% hist_range [-0.3 0.1] and y_view [-0.25 0.05]. On the medium run the tau
% LLEs spanned -0.26 .. +0.29 with a median of +0.008, so slightly over half
% the distribution was POSITIVE and 32% of it sat above the y_view ceiling in
% an overflow band that was itself off-screen: the panel read as almost
% entirely sub-zero when it was not. A table beside the figure
% (Fig_SFA_EOC_allStd_table.md) records, per level, n, median [IQR], mean,
% the share above zero and the overflow counts.
%
% THE RUN IS RESOLVED, NOT HARDCODED. This script used to carry
%   data_root = .../run_all_aug_14_26_17_25
% and three sibling figures carried the same line while a fourth pointed at a
% different, older run on a different model class -- so "regenerate the figures"
% silently built the manuscript from two runs. resolve_run_dir picks the newest
% run whose manifest names the requested preset, and errors if there is none.
%
% See also: resolve_run_dir, ParamSpaceAnalysis2/collect_level_values,
%           run_tau_sensitivity_analysis, paper_config

arguments
    cfg.verbose     (1,:) char    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (see verbose_level)
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
end

setup_paths();
out_dir      = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st           = manuscript_style();

run_dir = resolve_run_dir('run_dir', cfg.run_dir, 'preset_name', cfg.preset_name);

% --- Locate the tau sensitivity subfolder ----------------------------------
% By GLOB rather than by timestamped name, so re-pointing at another run needs
% no other edit.
tau_listing = dir(fullfile(run_dir, 'tau_sensitivity_*'));
tau_listing = tau_listing([tau_listing.isdir]);
if isempty(tau_listing)
    error('fig_sfa_EOC_allStd:NoTauDir', ...
        'No tau_sensitivity_* subfolder found in:\n  %s', run_dir);
end
if numel(tau_listing) > 1
    warning('fig_sfa_EOC_allStd:MultipleTauDirs', ...
        'Found %d tau_sensitivity_* subfolders; using the newest.', numel(tau_listing));
    [~, newest] = max([tau_listing.datenum]);
    tau_listing = tau_listing(newest);
end
tau_dir = fullfile(tau_listing.folder, tau_listing.name);

% --- Presentation constants ------------------------------------------------
fig_position = [457 500 420 340];
tick_fs   = st.tick_fs;
label_fs  = st.label_fs;
n_bins    = 40;        % lambda_1 bins between the 1st and 99th percentile
clim_frac = 0.8;       % darken: cap CLim at max count * clim_frac
% Colormap ramps white (0 counts) -> 90% black (max), not pure black, so the
% coloured median line stays visible over the darkest cells.
dark_cmap  = repmat(linspace(1, 0.1, 256)', 1, 3);
median_lw  = 2.5;
band_alpha = 0.25;

% --- Reload the tau PSA ------------------------------------------------------
psa = ParamSpaceAnalysis2.from_dir(tau_dir);

% The sweep's axis is tau_a_EI (E and I, 2026-09-14) or tau_a_E (older runs).
param = 'tau_a_E';
if isfield(psa.vector_param_lookup, 'tau_a_EI'); param = 'tau_a_EI'; end
if ~isfield(psa.vector_param_lookup, param)
    error('fig_sfa_EOC_allStd:NoTauAxis', ...
        'No tau_a_E / tau_a_EI vector parameter in %s (grid: %s).', tau_dir, strjoin(psa.grid_params, ', '));
end
lookup = psa.vector_param_lookup.(param);
x_tau  = cellfun(@(v) v(end), lookup);      % the slowest rung of every level
n_lev  = numel(x_tau);
if strcmp(param, 'tau_a_EI')
    x_label = 'slowest $\tau_a$ (E and I) (s)';
else
    x_label = 'slowest $\tau_a$ (E) (s)';
end
cond_names = cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false);
cond_use   = cond_names{1};                 % the tau sweep runs one condition
if numel(cond_names) > 1
    warning('fig_sfa_EOC_allStd:MultipleConditions', ...
        'The tau sweep carries %d conditions; drawing the first (%s).', numel(cond_names), cond_use);
end
cond_col = st.condition_color(cond_use);

% --- Collect every rep's lambda_1 per level ---------------------------------
vals = cell(1, n_lev);
for li = 1:n_lev
    vals{li} = ParamSpaceAnalysis2.collect_level_values(psa, param, li, cond_use, 'LLE');
end
all_vals = [vals{:}];
if isempty(all_vals)
    error('fig_sfa_EOC_allStd:NoData', 'No successful LLE values in %s.', tau_dir);
end
y_lo = prctile(all_vals, 1); y_hi = prctile(all_vals, 99);
pad  = 0.1 * max(y_hi - y_lo, eps);
y_lo = min(y_lo - pad, 0); y_hi = max(y_hi + pad, 0);
edges = linspace(y_lo, y_hi, n_bins + 1);
dy    = edges(2) - edges(1);
counts   = zeros(n_bins, n_lev);
n_above  = zeros(1, n_lev); n_below = zeros(1, n_lev);
n_rep    = zeros(1, n_lev); n_pos   = zeros(1, n_lev);
med = nan(1, n_lev); q1 = med; q3 = med; mu = med;
for li = 1:n_lev
    v = vals{li};
    n_rep(li)   = numel(v);
    n_pos(li)   = nnz(v > 0);
    n_above(li) = nnz(v > y_hi);
    n_below(li) = nnz(v < y_lo);
    inside = v(v >= y_lo & v <= y_hi);
    if ~isempty(inside); counts(:, li) = histcounts(inside, edges)'; end
    if ~isempty(v)
        med(li) = median(v); q1(li) = prctile(v, 25); q3(li) = prctile(v, 75); mu(li) = mean(v);
    end
end
% The two overflow rows are drawn as extra image rows of one bin height each.
img = [n_below; counts; n_above];
y_img = [y_lo - dy/2, y_hi + dy/2];        % centres of the first and last rows

% --- Draw --------------------------------------------------------------------
% x is the LEVEL INDEX (the sweep spaces the slowest tau linearly over its
% range, level_spacing 'linear'), ticked with the tau values in seconds.
cf = figure('Color', 'w', 'Position', fig_position);
ax = axes(cf); hold(ax, 'on');
imagesc(ax, 1:n_lev, y_img, img);
set(ax, 'YDir', 'normal');
colormap(ax, dark_cmap);
clim(ax, [0, max(1, max(img(:)) * clim_frac)]);
% Overflow rows: separated from the data by thin lines and labelled.
yline(ax, y_lo, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
yline(ax, y_hi, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
text(ax, 0.55, y_hi + dy/2, sprintf('> %.2f', y_hi), 'FontSize', 8, 'Color', [0.35 0.35 0.35], ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
text(ax, 0.55, y_lo - dy/2, sprintf('< %.2f', y_lo), 'FontSize', 8, 'Color', [0.35 0.35 0.35], ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
for li = 1:n_lev
    if n_above(li) > 0
        text(ax, li, y_hi + dy/2, sprintf('%d', n_above(li)), 'FontSize', 7, 'Color', [1 0.3 0.3], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
    end
    if n_below(li) > 0
        text(ax, li, y_lo - dy/2, sprintf('%d', n_below(li)), 'FontSize', 7, 'Color', [1 0.3 0.3], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
    end
end
% Zero line, band, median, mean.
yline(ax, 0, '-', 'Color', 'k', 'LineWidth', 1.5);
okb = isfinite(q1) & isfinite(q3);
xi = 1:n_lev;
fill(ax, [xi(okb), fliplr(xi(okb))], [q1(okb), fliplr(q3(okb))], cond_col, ...
    'FaceAlpha', band_alpha, 'EdgeColor', 'none');
plot(ax, xi, med, '-', 'Color', cond_col, 'LineWidth', median_lw);
plot(ax, xi, mu, 'o', 'Color', cond_col, 'MarkerFaceColor', 'w', 'MarkerSize', 4, 'LineWidth', 1);
% n reps and the share of reps above zero, along the top of the data area.
for li = 1:n_lev
    if n_rep(li) > 0
        text(ax, li, y_hi - 0.02 * (y_hi - y_lo), sprintf('%.0f%%', 100 * n_pos(li) / n_rep(li)), ...
            'FontSize', 7, 'Color', [0.2 0.2 0.2], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    end
end
text(ax, 0.55, y_lo + 0.03 * (y_hi - y_lo), sprintf('top row: share of reps with lambda_1 > 0 (n = %d per level)', max(n_rep)), 'FontSize', 7, 'Interpreter', 'none', ...
    'Color', [0.2 0.2 0.2], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
hold(ax, 'off');
xlim(ax, [0.5, n_lev + 0.5]);
ylim(ax, [y_lo - dy, y_hi + dy]);
set(ax, 'XTick', 1:n_lev, 'XTickLabel', arrayfun(@(v) sprintf('%.3g', v), x_tau, 'UniformOutput', false), ...
    'FontSize', tick_fs, 'Layer', 'top');
if n_lev > 7; set(ax, 'XTick', 1:2:n_lev, 'XTickLabel', arrayfun(@(v) sprintf('%.3g', v), x_tau(1:2:end), 'UniformOutput', false)); end
box(ax, 'off');
ylabel(ax, '$\lambda_1$ (1/s)', 'Interpreter', 'latex', 'FontSize', label_fs);
xlabel(ax, x_label, 'Interpreter', 'latex', 'FontSize', label_fs);
title(ax, sprintf('%s: all reps per level', st.condition_title(cond_use)), 'FontWeight', 'normal', 'FontSize', st.title_fs);

% --- Table ---------------------------------------------------------------------
hdr = '| level | slowest tau_a (s) | n | lambda_1 median [IQR] | mean | share > 0 | overflow above | overflow below |';
rows = cell(1, n_lev);
for li = 1:n_lev
    rows{li} = sprintf('| %d | %.3g | %d | %+.3f [%+.3f, %+.3f] | %+.3f | %.2f | %d | %d |', li, x_tau(li), n_rep(li), ...
        med(li), q1(li), q3(li), mu(li), n_pos(li) / max(1, n_rep(li)), n_above(li), n_below(li));
end
vprintf(cfg.verbose, 'verbose', '%s\n|---|---|---|---|---|---|---|---|\n%s\n', hdr, strjoin(rows, newline));

if ~cfg.visible; set(cf, 'Visible', 'off'); end

%% --- Second figure: WHERE the leading direction lives, vs max tau_a ----------
% The mechanistic form of the claim above. Per level of the sweep, the median
% over reps of the leading Lyapunov vector's norm^2 fractions in the x, SFA and
% STD blocks (stored per job by the top-K sweeps; see
% SRNNCellTypePairs.lya_summary). If the slowest SFA timescale sets lambda_1,
% the SFA fraction dominates until the STD mode takes over at the knee.
bf = gobjects(0);
has_blocks = isfield(psa.vector_param_lookup, param) && ...
    ~isempty(psa.results) && any(cellfun(@(c) has_field_in(psa.results, c, 'lead_frac_sfa'), ...
    cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false)));
if has_blocks
    bf = figure('Color', 'w', 'Position', [457 300 300 * numel(cond_names), 260]);
    blocks = {'lead_frac_sfa', 'lead_frac_std', 'lead_frac_x'};
    block_lab = {'SFA', 'STD', 'x'};
    block_col = [0.85 0.33 0.10; 0.00 0.45 0.74; 0.30 0.30 0.30];
    for ci = 1:numel(cond_names)
        ax_b = subplot(1, numel(cond_names), ci, 'Parent', bf); hold(ax_b, 'on');
        Y = nan(numel(x_tau), numel(blocks));
        for li = 1:numel(x_tau)
            for bi = 1:numel(blocks)
                v = ParamSpaceAnalysis2.collect_level_values(psa, param, li, cond_names{ci}, blocks{bi});
                if ~isempty(v); Y(li, bi) = median(v); end
            end
        end
        for bi = 1:numel(blocks)
            plot(ax_b, x_tau, Y(:, bi), '-o', 'Color', block_col(bi, :), 'LineWidth', 2, ...
                'MarkerFaceColor', block_col(bi, :), 'MarkerSize', 4, 'DisplayName', block_lab{bi});
        end
        hold(ax_b, 'off');
        set(ax_b, 'XScale', 'log', 'FontSize', tick_fs); box(ax_b, 'off');
        ylim(ax_b, [0 1]);
        xlabel(ax_b, x_label, 'Interpreter', 'latex', 'FontSize', label_fs);
        if ci == 1; ylabel(ax_b, 'leading-vector fraction', 'FontSize', label_fs); end
        if numel(cond_names) > 1
            title(ax_b, st.condition_title(cond_names{ci}), 'FontWeight', 'normal', 'FontSize', st.title_fs);
        end
        if ci == 1; legend(ax_b, 'Location', 'east', 'FontSize', 10); end
    end
    if ~cfg.visible; set(bf, 'Visible', 'off'); end
else
    warning('fig_sfa_EOC_allStd:NoBlockFractions', ...
        ['%s has no lead_frac_* fields (a Benettin-era run?); the leading-vector ' ...
         'block-fraction figure is skipped.'], tau_dir);
end

%% --- Save -------------------------------------------------------------------
fig_tag = 'Fig_SFA_EOC_allStd';
out = struct('figs', [cf, bf], 'files', {{}}, 'source', tau_dir);
if cfg.save
    save_figure_stable(out_dir, fig_tag, cf);
    out.files = existing_outputs(out_dir, fig_tag);
    fid = fopen(fullfile(out_dir, [fig_tag '_table.md']), 'w');
    if fid > 0
        fprintf(fid, '# lambda_1 vs the slowest SFA timescale\n\nAxis `%s`, condition `%s`, source `%s`.\n\n', param, cond_use, tau_dir);
        fprintf(fid, '%s\n|---|---|---|---|---|---|---|---|\n%s\n', hdr, strjoin(rows, newline));
        fclose(fid);
    end
    if isgraphics(bf)
        save_figure_stable(out_dir, 'Fig_SFA_EOC_blocks', bf);
        out.files = [out.files, existing_outputs(out_dir, 'Fig_SFA_EOC_blocks')];
    end
end
end

function tf = has_field_in(results, cond, field)
tf = false;
if ~isfield(results, cond); return; end
R = results.(cond); R = R(~cellfun(@isempty, R));
tf = ~isempty(R) && isfield(R{1}, field);
end



