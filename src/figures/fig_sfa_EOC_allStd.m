function out = fig_sfa_EOC_allStd(cfg)
% FIG_SFA_EOC_ALLSTD SFA edge of chaos: lambda_1 against the slowest tau_a.
%
%   out = FIG_SFA_EOC_ALLSTD()
%   out = FIG_SFA_EOC_ALLSTD('run_dir', d)
%   out = FIG_SFA_EOC_ALLSTD('run_dir', d, 'conditions', names)
%
% Replots the tau-sensitivity LLE panel -- how the largest Lyapunov exponent
% approaches 0 (the edge of chaos) as the slowest SFA adaptation timescale grows.
% No simulation is re-run: it reloads the saved tau sweep's PSA object (axis
% tau_a_EI, E and I together, since 2026-09-14; tau_a_E, E only, in older runs)
% and draws its condition(s): the fully-adapted regime(s) (see
% run_tau_sensitivity_analysis; one condition for older presets, the chained
% and the product MTS for the mu7chain presets).
%
% WHICH CONDITIONS ARE DRAWN: cfg.conditions (FD-006). Empty (the default):
% every condition the tau sweep holds, in its order. A cell of names: the
% selected conditions THAT THE TAU SWEEP HOLDS, in the selected order. This is
% the one figure where a selected condition the data lacks is NOT an error:
% the tau sweep holds only the fully-adapted conditions by design, so the
% paper's three-condition selection {no_adaptation, sfa1_std1, sfa3_stdc2}
% draws the chain alone, and the four-condition selection draws the chain and
% the product. The figure draws the INTERSECTION and errors
% (fig_sfa_EOC_allStd:NoSelectedCondition) only when it is empty. The table
% says which selected conditions the sweep does not hold.
%
% LAYOUT. One panel per drawn condition, side by side (420 px each). With ONE
% condition the figure is exactly what it was before the option existed: one
% panel, no title (the caption names the regime). With several, each panel is
% titled from srnn_condition_titles by condition NAME and coloured from
% manuscript_style by name; the density colour scale is shared so the panels
% read on one ruler. The leading-vector block figure has one subplot per
% drawn condition, as before.
%
% THE WHOLE DISTRIBUTION IS SHOWN (2026-09-14). Per level of the sweep the
% reps' lambda_1 values are drawn as a vertical density strip (a 2-D
% histogram, level x lambda_1, white -> dark grey), with the median as a thick
% line inside a 25th-75th percentile band, the mean as a marker, and a black
% lambda_1 = 0 line. The y range is FIXED at [-0.35 0.05] (ticks 0, -0.1,
% -0.2, -0.3) and the x axis is the slowest tau in seconds (ticks 0, 10, 20,
% 30); values outside the y range are counted in the table and NOT drawn, and
% the panel carries no annotation (TR, 2026-09-14; earlier that day it had two
% overflow rows with counts, a per-level share of reps with lambda_1 > 0
% along the top, and a caption, all removed).
%
% Until 2026-09-14 this figure restyled plot_sensitivity's panel with a fixed
% hist_range [-0.3 0.1] and y_view [-0.25 0.05]. On the medium run the tau
% LLEs spanned -0.26 .. +0.29 with a median of +0.008, so slightly over half
% the distribution was POSITIVE and 32% of it sat above the y_view ceiling in
% an overflow band that was itself off-screen: the panel read as almost
% entirely sub-zero when it was not. A table beside the figure
% (Fig_SFA_EOC_allStd_table.md) records, per condition and level, n, median
% [IQR], mean, the share above zero and the overflow counts.
%
% THE RUN IS RESOLVED, NOT HARDCODED. This script used to carry
%   data_root = .../run_all_aug_14_26_17_25
% and three sibling figures carried the same line while a fourth pointed at a
% different, older run on a different model class -- so "regenerate the figures"
% silently built the manuscript from two runs. resolve_run_dir picks the newest
% run whose manifest names the requested preset, and errors if there is none.
%
% See also: resolve_run_dir, ParamSpaceAnalysis2/collect_level_values,
%           run_tau_sensitivity_analysis, paper_config, select_conditions

arguments
    cfg.verbose     (1,:) char    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (see verbose_level)
    cfg.run_dir     (1,:) char    = ''
    cfg.preset_name (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.conditions                = {}          % names to draw, in order; empty = all the tau sweep holds (intersection, see header)
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
fig_position = [457 500 420 340];   % per panel; width grows with the panel count
tick_fs   = st.tick_fs;
label_fs  = st.label_fs;
n_bins    = 40;        % lambda_1 bins across y_range
clim_frac = 0.8;       % darken: cap CLim at max count * clim_frac
% Colormap ramps white (0 counts) -> 90% black (max), not pure black, so the
% coloured median line stays visible over the darkest cells.
dark_cmap  = repmat(linspace(1, 0.1, 256)', 1, 3);
median_lw  = 2.5;
band_alpha = 0.25;
y_range    = [-0.35 0.05];        % fixed lambda_1 window (TR, 2026-09-14); values outside are not drawn
y_ticks    = [-0.3 -0.2 -0.1 0];
x_ticks    = [0 10 20 30];         % slowest tau_a (s)

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
x_tau  = x_tau(:)';
n_lev  = numel(x_tau);
if strcmp(param, 'tau_a_EI')
    x_label = 'slowest $\tau_a$ (E and I) (s)';
else
    x_label = 'slowest $\tau_a$ (E) (s)';
end

% --- Which conditions: the selection INTERSECTED with the sweep (header) ----
held = cellfun(@(c) c.name, psa.conditions, 'UniformOutput', false);
requested = cfg.conditions;
if isempty(requested) || ((ischar(requested) || isstring(requested)) && all(strlength(requested) == 0))
    requested = {};                    % {}, [], '' and "" mean all, as in select_conditions
end
requested = cellstr(requested);
if isempty(requested)
    cond_draw = held;                  % the sweep's own order
    not_held  = {};
else
    % select_conditions over the requested list itself rejects duplicates and
    % malformed names; the intersection keeps the requested order.
    requested = select_conditions(requested, requested);
    in_sweep  = ismember(requested, held);
    cond_draw = requested(in_sweep);
    not_held  = requested(~in_sweep);
    if isempty(cond_draw)
        error('fig_sfa_EOC_allStd:NoSelectedCondition', ...
            ['None of the selected conditions {%s} is in the tau sweep, which holds {%s} ' ...
             '(%s).'], strjoin(requested, ', '), strjoin(held, ', '), tau_dir);
    end
end
n_draw = numel(cond_draw);

% --- Collect every rep's lambda_1 per condition and level --------------------
y_lo = y_range(1); y_hi = y_range(2);
edges = linspace(y_lo, y_hi, n_bins + 1);
y_ctr = (edges(1:end-1) + edges(2:end)) / 2;
S = struct('name', cond_draw, 'counts', [], 'n_rep', [], 'n_pos', [], 'n_above', [], ...
    'n_below', [], 'med', [], 'q1', [], 'q3', [], 'mu', []);
for ci = 1:n_draw
    vals = cell(1, n_lev);
    for li = 1:n_lev
        vals{li} = ParamSpaceAnalysis2.collect_level_values(psa, param, li, cond_draw{ci}, 'LLE');
    end
    if isempty([vals{:}])
        error('fig_sfa_EOC_allStd:NoData', 'No successful LLE values for %s in %s.', cond_draw{ci}, tau_dir);
    end
    % Fixed y range (TR, 2026-09-14): the paper's ruler for this panel. Values
    % outside it are counted in the table (n outside) and NOT drawn.
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
    S(ci).counts = counts; S(ci).n_rep = n_rep; S(ci).n_pos = n_pos;
    S(ci).n_above = n_above; S(ci).n_below = n_below;
    S(ci).med = med; S(ci).q1 = q1; S(ci).q3 = q3; S(ci).mu = mu;
end
% One colour scale for every panel, so they read on one ruler.
c_max = max(cellfun(@(c) max(c(:)), {S.counts}));

% --- Draw --------------------------------------------------------------------
% x is the slowest tau in seconds. The sweep spaces it linearly over its range
% (level_spacing 'linear'), so the image columns are uniform and imagesc can
% take the tau values directly; the ticks are fixed at 0, 10, 20, 30 s.
cf = figure('Color', 'w', 'Position', [fig_position(1:2), fig_position(3) * n_draw, fig_position(4)]);
for ci = 1:n_draw
    cond_col = st.condition_color(cond_draw{ci});
    if n_draw == 1
        ax = axes(cf); %#ok<LAXES>     % exactly the single-panel figure of before (runs once)
    else
        ax = subplot(1, n_draw, ci, 'Parent', cf);
    end
    hold(ax, 'on');
    imagesc(ax, x_tau, y_ctr, S(ci).counts);
    set(ax, 'YDir', 'normal');
    colormap(ax, dark_cmap);
    clim(ax, [0, max(1, c_max * clim_frac)]);
    % Zero line, band, median, mean.
    yline(ax, 0, '-', 'Color', 'k', 'LineWidth', 1.5);
    q1 = S(ci).q1; q3 = S(ci).q3;
    okb = isfinite(q1) & isfinite(q3);
    fill(ax, [x_tau(okb), fliplr(x_tau(okb))], [q1(okb), fliplr(q3(okb))], cond_col, ...
        'FaceAlpha', band_alpha, 'EdgeColor', 'none');
    plot(ax, x_tau, S(ci).med, '-', 'Color', cond_col, 'LineWidth', median_lw);
    plot(ax, x_tau, S(ci).mu, 'o', 'Color', cond_col, 'MarkerFaceColor', 'w', 'MarkerSize', 4, 'LineWidth', 1);
    hold(ax, 'off');
    dx = x_tau(2) - x_tau(1);
    xlim(ax, [min(0, x_tau(1) - dx/2), x_tau(end) + dx/2]);
    ylim(ax, y_range);
    set(ax, 'XTick', x_ticks, 'YTick', y_ticks, 'FontSize', tick_fs, 'Layer', 'top');
    box(ax, 'off');
    if ci == 1
        ylabel(ax, '$\lambda_1$ (1/s)', 'Interpreter', 'latex', 'FontSize', label_fs);
    end
    xlabel(ax, x_label, 'Interpreter', 'latex', 'FontSize', label_fs);
    if n_draw > 1
        title(ax, st.condition_title(cond_draw{ci}), 'FontWeight', 'normal', 'FontSize', st.title_fs);
    end
end

% --- Table ---------------------------------------------------------------------
hdr = '| level | slowest tau_a (s) | n | lambda_1 median [IQR] | mean | share > 0 | n above y range | n below y range |';
sep = '|---|---|---|---|---|---|---|---|';
tables = cell(1, n_draw);
for ci = 1:n_draw
    rows = cell(1, n_lev);
    for li = 1:n_lev
        rows{li} = sprintf('| %d | %.3g | %d | %+.3f [%+.3f, %+.3f] | %+.3f | %.2f | %d | %d |', li, x_tau(li), S(ci).n_rep(li), ...
            S(ci).med(li), S(ci).q1(li), S(ci).q3(li), S(ci).mu(li), ...
            S(ci).n_pos(li) / max(1, S(ci).n_rep(li)), S(ci).n_above(li), S(ci).n_below(li));
    end
    tables{ci} = sprintf('%s\n%s\n%s\n', hdr, sep, strjoin(rows, newline));
    if n_draw > 1
        vprintf(cfg.verbose, 'verbose', '%s:\n', cond_draw{ci});
    end
    vprintf(cfg.verbose, 'verbose', '%s', tables{ci});
end

if ~cfg.visible; set(cf, 'Visible', 'off'); end

%% --- Second figure: WHERE the leading direction lives, vs max tau_a ----------
% The mechanistic form of the claim above. Per level of the sweep, the median
% over reps of the leading Lyapunov vector's norm^2 fractions in the x, SFA and
% STD blocks (stored per job by the top-K sweeps; see
% SRNNCellTypePairs.lya_summary). If the slowest SFA timescale sets lambda_1,
% the SFA fraction dominates until the STD mode takes over at the knee.
% One subplot per DRAWN condition, in the drawn order.
bf = gobjects(0);
has_blocks = isfield(psa.vector_param_lookup, param) && ...
    ~isempty(psa.results) && any(cellfun(@(c) has_field_in(psa.results, c, 'lead_frac_sfa'), cond_draw));
if has_blocks
    bf = figure('Color', 'w', 'Position', [457 300 300 * n_draw, 260]);
    blocks = {'lead_frac_sfa', 'lead_frac_std', 'lead_frac_x'};
    block_lab = {'SFA', 'STD', 'x'};
    block_col = [0.85 0.33 0.10; 0.00 0.45 0.74; 0.30 0.30 0.30];
    for ci = 1:n_draw
        ax_b = subplot(1, n_draw, ci, 'Parent', bf); hold(ax_b, 'on');
        Y = nan(numel(x_tau), numel(blocks));
        for li = 1:numel(x_tau)
            for bi = 1:numel(blocks)
                v = ParamSpaceAnalysis2.collect_level_values(psa, param, li, cond_draw{ci}, blocks{bi});
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
        if n_draw > 1
            title(ax_b, st.condition_title(cond_draw{ci}), 'FontWeight', 'normal', 'FontSize', st.title_fs);
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
        if n_draw == 1
            % Byte-for-byte the single-condition table of before.
            fprintf(fid, '# lambda_1 vs the slowest SFA timescale\n\nAxis `%s`, condition `%s`, source `%s`.\n\n', param, cond_draw{1}, tau_dir);
            fprintf(fid, '%s', tables{1});
        else
            fprintf(fid, '# lambda_1 vs the slowest SFA timescale\n\nAxis `%s`, conditions %s (one panel each, left to right), source `%s`.\n\n', ...
                param, strjoin(cellfun(@(n) ['`' n '`'], cond_draw, 'UniformOutput', false), ', '), tau_dir);
        end
        if ~isempty(not_held)
            if n_draw == 1; fprintf(fid, '\n'); end   % after the table, not inside it
            fprintf(fid, ['Selected but not in the tau sweep (which holds only the fully-adapted ' ...
                'conditions), so not drawn: %s.\n\n'], strjoin(not_held, ', '));
        end
        if n_draw > 1
            for ci = 1:n_draw
                fprintf(fid, '## %s\n\n%s\n', cond_draw{ci}, tables{ci});
            end
        end
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
