function out = fig_eig_heatmap(cfg)
% FIG_EIG_HEATMAP Jacobian-eigenvalue occupancy, one panel per adaptation regime.
%
%   out = FIG_EIG_HEATMAP()
%   out = FIG_EIG_HEATMAP('data_file', f)
%   out = FIG_EIG_HEATMAP('density_scale', 'loglog')
%
% A Gaussian-smoothed 2-D DENSITY over the complex plane for each of the four
% adaptation regimes -- an "occupancy" heatmap showing how much time the
% instantaneous Jacobian's eigenvalues spend in each region, and in particular
% to the RIGHT of the imaginary axis (Re > 0, locally unstable). Panels share
% axis limits and one log-density colorbar.
%
% density_scale picks the colour transform of the smoothed density D:
%   'log'    -> log10(1 + D)               (default; the original figure)
%   'loglog' -> log10(1 + log10(1 + D))    (compresses the peaks further, so the
%                                           sparse outer eigenvalue cloud reads)
% The transform is applied HERE and the helper is handed plain values, so the
% shared plotter in src/plotting/ stays untouched.
%
% WHICH CONDITIONS ARE DRAWN: cfg.conditions (FD-006), a cell of condition
% names drawn as panels in the order given; empty (the default) draws every
% condition the .mat holds, in its own order, exactly as before. Resolved by
% select_conditions against the saved cond_names, so a name the run lacks is
% an error. The shared window limits and colour scale are computed from the
% DRAWN conditions only. Panel titles come from the condition NAME
% (titles_by_name), not from the titles run_eig_heatmap saved. Layout:
% floor(sqrt(N)) rows, so three conditions are 1 x 3 (the paper's set) and
% four are 2 x 2 (the review set; accepted, Tom 2026-11-12).
%
% THE RITZ OVERLAY (FD-033; 'ritz_overlay', off by default, and then the
% figure is the one it always was). The density is of the eigenvalues of the
% whole Jacobian at each sampled state: every state direction at once, most
% of which a perturbation never follows. The Lyapunov exponents belong to
% the few directions it does follow, and are real. The overlay joins them:
% for the SAME sampled states, the eigenvalues of the Jacobian projected on
% the leading Lyapunov subspace, H = Q' J Q (Ritz values; run_eig_heatmap
% with ritz_K > 0 saves H), drawn as points in the panel's condition color,
% all sampled states pooled as the density pools them. Options:
%   ritz_k         how many leading directions (0 = all the run tracked);
%                  the leading k x k block of the saved H is the top-k
%                  subspace, so any k up to the run's K can be drawn
%   ritz_scaling   'raw' (default) or 'type_sd': which of the run's state
%                  scalings; Ritz values depend on how x, the adaptation and
%                  the depression variables are scaled against one another
%   ritz_style     how the Ritz values are drawn. The eigenvalue density is
%                  a (log-scaled) count, and a plain scatter is not: a pile
%                  of a thousand points at one spot is one dot and a few
%                  dozen stray points look like many, so a scatter
%                  over-represents the outliers (Tom, 2026-11-19).
%                    'density'  (default) the Ritz values as a density of
%                               their own on the eigenvalue density's bins,
%                               the same log transform, in the condition's
%                               color, opacity following the count
%                    'marginal' very small transparent points, plus a
%                               histogram of Re along the bottom of the
%                               panel (linear counts, so the pile reads as
%                               a pile)
%                    'fine'     very small transparent points only
%                    'scatter'  the plain scatter, every point a visible dot
%   ritz_exponent_ticks  black ticks along the bottom edge at the top-k
%                  Lyapunov exponents of the same pass (the time averages,
%                  over the sampled window, of the local rates whose
%                  instantaneous values are the diagonal of H; default true)
%   ritz_inset     an inset on the real axis with the eigenvalue density
%                  recomputed at its scale, the Ritz values as a histogram
%                  of Re above it and the exponent ticks: the slow band of
%                  an adapting network, near -0.1 1/s, cannot be told from
%                  zero on a window tens of 1/s wide. Its range is chosen
%                  PER CONDITION from where that condition's Ritz values
%                  and exponents are (the 5th to 95th percentile of the real
%                  parts, the exponents and zero, padded), or is
%                  ritz_inset_xlim when that is given
%   NO TEXT ON THE PANELS (Tom, 2026-11-19; the repo's convention):
%   the counts and what the ticks are go in the table beside the figure,
%   which also carries a caption paragraph. ritz_annotate and
%   ritz_tick_label (both default false) print them on the panels for a
%   working copy. With the overlay the figure has ONE title line.
% The overlay is saved under its own tag, fig_eig_heatmap_ritz (with
% _typesd for the scaled version and _inset for the inset), with a table
% beside it, <tag>_table.md: per condition where the points sit, how many
% have positive real part, how many are complex and at what frequency, the
% share of the tracked subspace in dendritic coordinates, and how far the
% other scaling moves them. The tags share a prefix and save_figure_stable
% clears by prefix, so give each version its own out_dir. HOW TO READ IT:
% the tracked subspace is carried by the dynamics and is not an invariant
% subspace of the instantaneous Jacobian, so the points need not fall on the
% density (they lie in its numerical range); single points do not average
% to single exponents, only the sums over the top k agree with the local
% Lyapunov rates; and in a stable adapting network the least-contracting
% directions are a band of slow single-neuron adaptation modes, about one
% per neuron, so with k far below the network size the points mark that
% band and not the network's collective modes.
% A point right of zero is not
% an instability (H is not normal; every exponent can be negative), and the
% spread of the points depends strongly on the norm: on the 500-neuron
% multiple-timescale network the raw points spread over several 1/s where
% the type_sd points sit within a few tenths of -0.1 1/s (2026-11-19).
%
% SPECTRA_SOURCE (FD-034; 'eig_heatmap' by default, the figure as it was).
% 'local_lyapunov' draws the same figure from the states run_local_lyapunov
% sampled along its switching-input run ('spectra_n_samples';
% <run_dir>/local_lyapunov/local_lyapunov_spectra.mat): the density is of
% the Jacobian at states of that run, and the Ritz values are those of the
% top-K pass whose local rates the composite Figure 3 plots in its expansion
% row, so the diagonal of each H is the instantaneous value of those rates.
% Saved under fig_eig_heatmap_switch (with the overlay's suffixes), and the
% table gains a section on the shared pass: how the diagonal of H compares
% with the pass's interval rates, the first against the second half of the
% window, and the numbers by input step.
%
% PLOTTING HALF ONLY. run_eig_heatmap does the sampling and writes
% eig_heatmap_data.mat, so the look can be iterated without re-simulating.
%
% See also: run_eig_heatmap, manuscript_style

arguments
    cfg.verbose     (1,:) char    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (see verbose_level)
    cfg.data_file   (1,:) char    = ''    % '' -> search run_dir, then data/eig_heatmap
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''    % the run whose eig_heatmap data to plot
    cfg.preset_name (1,:) char    = ''    % unused; the preset is recorded in the .mat
    cfg.density_scale (1,:) char {mustBeMember(cfg.density_scale, {'log', 'loglog'})} = 'log'
    cfg.conditions                = {}    % names to draw, in order; empty = all the run holds (select_conditions)
    cfg.ritz_overlay (1,1) logical = false   % draw the Ritz values of the leading Lyapunov subspace (FD-033)
    cfg.ritz_k       (1,1) double {mustBeInteger, mustBeNonnegative} = 0   % leading directions drawn; 0 = all the run tracked
    cfg.ritz_scaling (1,:) char {mustBeMember(cfg.ritz_scaling, {'raw', 'type_sd'})} = 'raw'
    cfg.ritz_style   (1,:) char {mustBeMember(cfg.ritz_style, {'density', 'marginal', 'fine', 'scatter'})} = 'density'
    cfg.ritz_exponent_ticks (1,1) logical = true
    cfg.ritz_annotate (1,1) logical = false      % print the number of points and the share right of zero on each panel
    cfg.ritz_tick_label (1,1) logical = false    % print what the exponent ticks are on each panel
    cfg.ritz_inset   (1,1) logical = false
    cfg.ritz_inset_xlim double = []              % real-axis range of the inset (1/s); [] = chosen per condition
    cfg.spectra_source (1,:) char {mustBeMember(cfg.spectra_source, {'eig_heatmap', 'local_lyapunov'})} = 'eig_heatmap'   % which stage sampled the states (FD-034)
end

setup_paths();
out_dir      = default_out_dir(cfg.out_dir, mfilename('fullpath'));
st           = manuscript_style(); %#ok<NASGU>

% INSIDE THE RUN DIRECTORY ONLY. This used to load a hardcoded .mat sitting
% beside this file, with run_dir marked "unused" -- so on 2026-08-26 it plotted
% Aug 22 data while every other figure used the Aug 25 sweep, and nothing said
% so. The fix added the run directory but KEPT a data/eig_heatmap fallback, which
% is the same failure one step further out; that tier is now gone. Pass the
% standalone location AS run_dir, or pass data_file directly.
switch cfg.spectra_source
    case 'eig_heatmap'
        data_file = resolve_data_file(cfg.data_file, cfg.run_dir, ...
            {fullfile(cfg.run_dir, 'eig_heatmap')}, ...
            'eig_heatmap_data.mat', ...
            'Run run_eig_heatmap first');
    case 'local_lyapunov'
        data_file = resolve_data_file(cfg.data_file, cfg.run_dir, ...
            {fullfile(cfg.run_dir, 'local_lyapunov')}, ...
            'local_lyapunov_spectra.mat', ...
            'Run run_local_lyapunov with ''spectra_n_samples'' > 0 first');
end
D = load(data_file);
if isfield(D, 'cond_names')
    held = D.cond_names;
else
    held = D.condition_titles;   % a .mat from before cond_names was saved: titles stand in for names
end
[cond_names, idx] = select_conditions(held, cfg.conditions);
evals_by_cond    = D.evals_by_cond(idx);
condition_titles = titles_by_name(cond_names, D.condition_titles(idx));
lle_by_cond      = D.lle_by_cond(idx);
lle_window       = D.lle_window;
n_cond           = numel(cond_names);
if cfg.ritz_overlay
    % Read only when the overlay is asked for, so a run directory made
    % before the Ritz pass existed, or without it, draws as it always did.
    if ~isfield(D, 'ritz_by_cond')
        error('fig_eig_heatmap:NoRitz', ...
            ['%s holds no Ritz pass. Run run_eig_heatmap with ''ritz_K'' > 0 (a bundle: ' ...
             'cfg.eig_heatmap_ritz_K), or draw without ''ritz_overlay''.'], data_file);
    end
    ritz = cell(1, n_cond);
    for i = 1:n_cond
        ritz{i} = ritz_points(D.ritz_by_cond{idx(i)}, cfg.ritz_k, cfg.ritz_scaling, cond_names{i});
    end
end

%% ---- Heatmap / plotting parameters ----------------------------------------
grid_res   = 250;     % heatmap bins per axis
sigma_bins = 1.25;    % Gaussian smoothing width (bins)
keep_frac  = 0.999;   % fraction of eigenvalue density to keep inside the window

switch cfg.density_scale
    case 'log'
        scale_fn = @(D) log10(1 + D);
        cb_label = 'log_{10}(1 + eigenvalue density)';
    case 'loglog'
        scale_fn = @(D) log10(1 + log10(1 + D));
        cb_label = 'log_{10}(1 + log_{10}(1 + eigenvalue density))';
end


%% ---- Global, square, density-trimmed limits (shared for comparability) -----
% Trim to a SQUARE window (equal Re/Im span) that still contains keep_frac of
% the pooled eigenvalues, so the panels zoom in on where the density actually
% is instead of stretching to the extreme outliers.
all_evals = vertcat(evals_by_cond{:});
re = real(all_evals); im = imag(all_evals);

tail = (1 - keep_frac) / 2;                 % drop this fraction off each Re end
re_lo = quantile(re, tail);
re_hi = quantile(re, 1 - tail);
re_ctr = (re_lo + re_hi) / 2;
re_span = re_hi - re_lo;

im_half = quantile(abs(im), keep_frac);     % Im density is symmetric about 0
im_span = 2 * im_half;

% Common span = the larger of the two, so >= keep_frac is retained on each axis
% and the window is square. Re is centered on its own midpoint; Im on 0.
span = max(re_span, im_span);
re_lim = re_ctr + [-0.5, 0.5] * span;
im_lim = [-0.5, 0.5] * span;

re_edges = linspace(re_lim(1), re_lim(2), grid_res + 1);
im_edges = linspace(im_lim(1), im_lim(2), grid_res + 1);

D_by_cond = cell(n_cond, 1);
cmax = 0;
for i = 1:n_cond
    D_by_cond{i} = compute_eigenvalue_density( ...
        evals_by_cond{i}, re_edges, im_edges, sigma_bins);
    D_by_cond{i} = scale_fn(D_by_cond{i});               % colour transform, see density_scale
    cmax = max(cmax, max(D_by_cond{i}(:)));              % shared scaled-density max
end
clim = [0, max(cmax, 0.5)];   % the colorbar is ticked at 0, 0.25, 0.5

%% ---- Assemble the comparison figure ---------------------------------------
% Grid derived from the condition count, not hardcoded. This was `2, 2`, which
% fit the four adaptation regimes the preset had when it was written; the paper
% now runs SEVEN and `nexttile` throws "The layout does not have sufficient
% space" on the fifth. It went unnoticed because the figure was reading a
% four-condition .mat frozen beside its own .m -- fixing the data resolution is
% what surfaced this.
%
% floor(sqrt(.)) keeps 4 -> 2x2 exactly as before, and gives 7 -> 2x4 (one blank
% tile). Panel size is held roughly constant, so the window grows with the grid
% instead of squeezing panels.
n_rows = max(1, floor(sqrt(n_cond)));
n_cols = ceil(n_cond / n_rows);
fig = figure('Color', 'w', 'Position', [200, 150, 450*n_cols, 380*n_rows]);
tl  = tiledlayout(fig, n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

ax_panels = gobjects(n_cond, 1);
for i = 1:n_cond
    ax_panels(i) = nexttile(tl);
    plot_eigenvalue_heatmap_helper( ...
        ax_panels(i), D_by_cond{i}, re_edges, im_edges, clim, false);  % already scaled above
    title(ax_panels(i), condition_titles{i}, 'FontWeight', 'normal', 'FontSize', 14);

    % Finite-time lambda_1 over the last lle_window seconds, top-left.
    text(ax_panels(i), 0.03, 0.96, sprintf('\\lambda_1 = %+.3f', lle_by_cond(i)), ...
        'Units', 'normalized', 'Color', [0.25 0.25 0.25], 'FontSize', 12, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

    if cfg.ritz_overlay
        draw_ritz(ax_panels(i), ritz{i}, ritz_color(cond_names{i}), re_edges, im_edges, scale_fn, cfg);
    end
end

% Single shared colorbar for the whole layout (all panels share clim + colormap).
cb = colorbar(ax_panels(end));
cb.Layout.Tile = 'east';
cb.Label.String = cb_label;
cb.Ticks = [0 0.25 0.5];
cb.Box = 'off';

heading = {'Jacobian eigenvalue occupancy across adaptation regimes', ...
    sprintf('\\lambda_1 = finite-time largest Lyapunov exponent over the last %g s', lle_window)};
if cfg.ritz_overlay
    % One short line; what the points and lambda_1 are is in the table and the caption.
    heading = 'Jacobian eigenvalues (gray) and Ritz values of the leading Lyapunov directions (color)';
end
title(tl, heading, 'FontWeight', 'bold');

if cfg.ritz_overlay && cfg.ritz_inset
    % After the layout is complete (title, colorbar), so the plot boxes are final.
    for i = 1:n_cond
        draw_ritz_inset(fig, ax_panels(i), i, evals_by_cond{i}, ritz{i}, ritz_color(cond_names{i}), ...
            inset_range(ritz{i}, cfg.ritz_inset_xlim), scale_fn, cfg.ritz_exponent_ticks);
    end
end

if ~cfg.visible; set(fig, 'Visible', 'off'); end

%% --- Save -------------------------------------------------------------------
fig_tag = 'fig_eig_heatmap';
if strcmp(cfg.spectra_source, 'local_lyapunov'); fig_tag = [fig_tag '_switch']; end
if cfg.ritz_overlay
    fig_tag = [fig_tag '_ritz'];
    if ~strcmp(cfg.ritz_style, 'density'); fig_tag = [fig_tag '_' cfg.ritz_style]; end
    if strcmp(cfg.ritz_scaling, 'type_sd'); fig_tag = [fig_tag '_typesd']; end
    if cfg.ritz_inset; fig_tag = [fig_tag '_inset']; end
    table_text = ritz_table(D, idx, cond_names, condition_titles, ritz, cfg, re_lim, im_lim);
    vprintf(cfg.verbose, 'verbose', '%s\n', table_text);
end
out = struct('figs', fig, 'files', {{}}, 'source', data_file);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    if cfg.ritz_overlay
        fid = fopen(fullfile(out_dir, [fig_tag '_table.md']), 'w');
        if fid > 0
            fprintf(fid, '%s\n', table_text);
            fclose(fid);
        end
    end
end
end

%% ------------------------------------------------------------------------
function P = ritz_points(R, k, scaling, cond_name)
% The Ritz values of one condition to draw: the eigenvalues of the leading
% k x k block of the saved H at every valid sampled state, for the named
% scaling. P.z is k_used x states (complex); k is capped at the run's K.
sc = find(strcmp(R.scalings, scaling), 1);
if isempty(sc)
    error('fig_eig_heatmap:NoRitzScaling', ...
        'The run holds no ''%s'' Ritz pass for %s (it holds: %s).', scaling, cond_name, strjoin(R.scalings, ', '));
end
if k == 0; k = R.K; end
k_used = min(k, R.K);
states = find(R.valid(:))';
has_gram = isfield(R, 'JQ_gram');          % a run made before the residuals were saved has none
z = complex(zeros(k_used, numel(states)));
pair_res = nan(k_used, numel(states));
names = {'residual', 'share_diag_pos', 'share_ritz_pos', 'norm_diag', 'norm_lower', 'norm_upper_T', ...
    'coupling_TS', 'num_abscissa', 'max_re'};
S = cell2struct(repmat({nan(1, numel(states))}, numel(names), 1), names(:), 1);
for s = 1:numel(states)
    Hs = R.H(:, :, states(s), sc);
    if has_gram
        Dg = ritz_diagnostics(Hs, R.JQ_gram(:, :, states(s), sc), k_used);
        z(:, s) = Dg.values; pair_res(:, s) = Dg.pair_residual;
        for n = 1:numel(names); S.(names{n})(s) = Dg.(names{n}); end
    else
        z(:, s) = eig(Hs(1:k_used, 1:k_used));
    end
end
trace_k = arrayfun(@(s) trace(R.H(1:k_used, 1:k_used, s, sc)), states);
P = struct('z', z, 'pair_residual', pair_res, 'state', S, 'has_gram', has_gram, 'k', k_used, ...
    'k_requested', k, 'K', R.K, 'N', R.N, 'sc', sc, ...
    'states', states, 'exponents', R.LE_spectrum(1:k_used, sc), 'trace_k', trace_k, ...
    'x_share', mean(R.x_fraction(states, 1:k_used, sc), 2));
if isfield(R, 'full_num_abscissa')
    P.full = struct('max_re', R.full_max_re(states)', 'n_positive', R.full_n_positive(states)', ...
        'num_abscissa', R.full_num_abscissa(states, sc)');
else
    P.full = [];
end
end

function col = ritz_color(cond_name)
% The condition's color (manuscript_style), as every other figure uses it.
st = manuscript_style();
if st.condition_color.isKey(cond_name)
    col = st.condition_color(cond_name);
else
    col = [0.85 0.33 0.10];
end
end

function draw_ritz(ax, P, col, re_edges, im_edges, scale_fn, cfg)
% The Ritz values over the density, in the style asked for; the axes'
% limits are the density's and stay so (points outside are clipped and
% counted in the table). Every object is tagged ritz_*.
xl = xlim(ax); yl = ylim(ax);
z = P.z(:);
re_c = (re_edges(1:end - 1) + re_edges(2:end)) / 2;
im_c = (im_edges(1:end - 1) + im_edges(2:end)) / 2;
hold(ax, 'on');
switch cfg.ritz_style
    case 'scatter'
        % A thin white rim keeps a dark point readable on the dark part of the density.
        scatter(ax, real(z), imag(z), 10, col, 'filled', 'MarkerFaceAlpha', 0.45, ...
            'MarkerEdgeColor', 'w', 'MarkerEdgeAlpha', 0.5, 'LineWidth', 0.25, 'Tag', 'ritz_values');
    case {'fine', 'marginal'}
        scatter(ax, real(z), imag(z), 2.5, col, 'filled', 'MarkerFaceAlpha', 0.2, ...
            'MarkerEdgeColor', 'none', 'Tag', 'ritz_values');
        if strcmp(cfg.ritz_style, 'marginal')
            % Histogram of Re on the density's bins, LINEAR counts scaled to
            % the bottom fifth of the panel: where most of the points are.
            n = histcounts(real(z), re_edges);
            h = yl(1) + 0.2 * diff(yl) * n / max(n);
            xs = [re_edges(1:end - 1); re_edges(2:end)]; ys = [h; h];
            patch(ax, [xs(1); xs(:); xs(end)], [yl(1); ys(:); yl(1)], col, 'FaceAlpha', 0.55, ...
                'EdgeColor', col, 'LineWidth', 0.5, 'Tag', 'ritz_marginal');
        end
    case 'density'
        % The same binning, smoothing and log transform as the eigenvalue
        % density, normalized to its own maximum; a solid color whose
        % opacity is the transformed count. Black on a gray density would
        % not read, so a colorless condition color is drawn in dark red.
        if max(col) - min(col) < 0.1; col = [0.75 0.10 0.15]; end
        Dz = scale_fn(compute_eigenvalue_density(z, re_edges, im_edges, 1));
        A = Dz' / max(Dz(:));
        image(ax, 'XData', re_c, 'YData', im_c, 'CData', repmat(reshape(col, 1, 1, 3), numel(im_c), numel(re_c)), ...
            'AlphaData', 0.9 * A, 'Tag', 'ritz_density');
end
if cfg.ritz_exponent_ticks
    % One line object with NaN breaks: the Lyapunov exponents of the same
    % pass as ticks on the bottom edge. (One object of many points, so a
    % composite that looks for the two-point Re = 0 line does not find it.)
    lam = P.exponents(:)';
    len = 0.03 * diff(yl);
    xs = [lam; lam; nan(1, numel(lam))];
    ys = repmat([yl(1); yl(1) + len; NaN], 1, numel(lam));
    plot(ax, xs(:), ys(:), '-', 'Color', [0 0 0], 'LineWidth', 0.6, 'Tag', 'ritz_exponent_ticks');
    if cfg.ritz_tick_label
        text(ax, 0.02, 0.04, sprintf('| Lyapunov exponents (top %d)', P.k), 'Units', 'normalized', ...
            'FontSize', 8, 'Color', [0 0 0], 'VerticalAlignment', 'bottom', 'Tag', 'ritz_tick_label');
    end
end
if cfg.ritz_annotate
    n_pos_state = sum(real(P.z) > 0, 1);
    % Under the lambda_1 label of the stand-alone figure, top left.
    text(ax, 0.03, 0.87, sprintf('%d Ritz values\n(%d states x %d)\nRe > 0: %.1f%%,\nin %d of %d states', ...
        numel(z), size(P.z, 2), P.k, 100 * mean(real(z) > 0), nnz(n_pos_state > 0), size(P.z, 2)), ...
        'Units', 'normalized', 'FontSize', 8, 'Color', col, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', 'Tag', 'ritz_annotation');
end
hold(ax, 'off');
xlim(ax, xl); ylim(ax, yl);
end

function lim = inset_range(P, given)
% The inset's real-axis range: the one given, or chosen from where this
% condition's Ritz values and exponents are: the 5th to 95th percentile of
% the real parts, widened to hold every exponent and zero, padded by 15%
% and rounded outward to two significant digits.
if ~isempty(given)
    if numel(given) ~= 2 || given(2) <= given(1)
        error('fig_eig_heatmap:BadInsetRange', 'ritz_inset_xlim must be [lo hi] with lo < hi, or [].');
    end
    lim = given(:)'; return
end
re = real(P.z(:));
lo = min([quantile(re, 0.05); P.exponents(:); 0]);
hi = max([quantile(re, 0.95); P.exponents(:); 0]);
pad = 0.15 * max(hi - lo, 0.1);
lo = lo - pad; hi = hi + pad;
step = 10 ^ (floor(log10(hi - lo)) - 1);
lim = [floor(lo / step), ceil(hi / step)] * step;
end

function draw_ritz_inset(fig, ax, panel, evals, P, col, re_lim_inset, scale_fn, with_ticks)
% A small axes in the lower right corner of the panel's plot box, zoomed on
% the real axis near zero, where the Ritz values of an adapting network
% pile up. Lower half: the density of the full Jacobian's eigenvalues at
% that scale (its own gray scale, darkest = its own maximum) with the Ritz
% values as small points and the exponent ticks on the bottom edge. Upper
% half: a histogram of the real parts of the Ritz values in that range
% (linear counts), so the pile reads as a pile. The imaginary range is the
% 90th percentile of the points in range, so a few far points do not
% flatten it. Parented to the figure, tagged 'ritz_inset', with the panel
% index and its place within the plot box in UserData so a composite can
% put it back (fig_main3_stability).
rel = [0.50 0.13 0.47 0.30];                  % [x y w h] as fractions of the panel's plot box
box = plot_box_position(ax);
z = P.z(:);
in_x = real(z) >= re_lim_inset(1) & real(z) <= re_lim_inset(2);
im_half = 0.2;
if any(in_x); im_half = max(0.2, 1.5 * quantile(abs(imag(z(in_x))), 0.9)); end
re_edges = linspace(re_lim_inset(1), re_lim_inset(2), 91);
im_edges = linspace(-im_half, im_half, 31);
Dn = scale_fn(compute_eigenvalue_density(evals, re_edges, im_edges, 1));
ai = axes(fig, 'Units', 'normalized', ...
    'Position', [box(1) + rel(1) * box(3), box(2) + rel(2) * box(4), rel(3) * box(3), rel(4) * box(4)]);
% The image fills the lower half of the inset's y range; the histogram the upper half.
imagesc(ai, (re_edges(1:end - 1) + re_edges(2:end)) / 2, (im_edges(1:end - 1) + im_edges(2:end)) / 2, Dn');
axis(ai, 'xy'); colormap(ai, flipud(gray));
clim(ai, [0, max(max(Dn(:)), eps)]);
hold(ai, 'on');
st = manuscript_style();
top = 3.4 * im_half;                           % y of the top of the histogram band
plot(ai, [0 0], [-im_half top], '--', 'Color', st.zeroline_color, 'LineWidth', 1, 'Tag', 'ritz_inset_zero');
draw_col = col; if max(col) - min(col) < 0.1; draw_col = [0.75 0.10 0.15]; end
scatter(ai, real(z), imag(z), 3, draw_col, 'filled', 'MarkerFaceAlpha', 0.3, ...
    'MarkerEdgeColor', 'none', 'Tag', 'ritz_values_inset');
n = histcounts(real(z), re_edges);
if any(n)
    h = 1.2 * im_half + (top - 1.3 * im_half) * n / max(n);
    xs = [re_edges(1:end - 1); re_edges(2:end)]; ys = [h; h];
    patch(ai, [xs(1); xs(:); xs(end)], [1.2 * im_half; ys(:); 1.2 * im_half], draw_col, 'FaceAlpha', 0.6, ...
        'EdgeColor', draw_col, 'LineWidth', 0.5, 'Tag', 'ritz_inset_histogram');
end
plot(ai, re_lim_inset, [1 1] * 1.2 * im_half, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
if with_ticks
    lam = P.exponents(:)';
    xs = [lam; lam; nan(1, numel(lam))];
    ys = repmat([-im_half; -0.7 * im_half; NaN], 1, numel(lam));
    plot(ai, xs(:), ys(:), '-', 'Color', [0 0 0], 'LineWidth', 0.6, 'Tag', 'ritz_inset_exponent_ticks');
end
hold(ai, 'off');
set(ai, 'XLim', re_lim_inset, 'YLim', [-im_half top], 'Box', 'on', 'FontSize', 8, ...
    'TickDir', 'out', 'TickLength', [0.02 0.02], 'Layer', 'top', ...
    'XTick', unique([re_lim_inset(1), 0, re_lim_inset(2)]), ...
    'XTickLabel', arrayfun(@(v) sprintf('%g', v), unique([re_lim_inset(1), 0, re_lim_inset(2)]), 'UniformOutput', false), ...
    'YTick', [0 im_half], 'YTickLabel', {'0', sprintf('%.2g', im_half)}, ...
    'XTickLabelRotation', 0, ...
    'Tag', 'ritz_inset', 'UserData', struct('panel', panel, 'rel', rel));
end

function txt = ritz_table(D, idx, cond_names, titles, ritz, cfg, re_lim, im_lim)
% The numbers behind the overlay, as markdown. One row per condition.
other = setdiff({'raw', 'type_sd'}, {cfg.ritz_scaling});
L = {};
L{end + 1} = '# Ritz values of the Jacobian on the leading Lyapunov subspace';
L{end + 1} = '';
R1 = D.ritz_by_cond{idx(1)};
L{end + 1} = sprintf(['Drawn: the eigenvalues of the leading k x k block of H = Q''JQ at every sampled state, ' ...
    '`%s` norm, k = %d requested. Pass: top-%d, %g s of alignment before the first sampled state, window ' ...
    '[%g, %g] s, re-orthonormalized every %g s, basis seed %d. Same sampled states as the density.'], ...
    cfg.ritz_scaling, ritz{1}.k_requested, R1.K, R1.warmup_s, R1.interval(1), R1.interval(2), R1.lya_dt, R1.seed);
if min(R1.t) > R1.interval(1) + 1e-9
    % The window opens before the first sampled state (run_local_lyapunov's
    % spectra_start_s): the alignment ends where the window starts.
    L{end} = sprintf(['Drawn: the eigenvalues of the leading k x k block of H = Q''JQ at every sampled state, ' ...
        '`%s` norm, k = %d requested. Pass: top-%d, %g s of alignment over [%g, %g] s, window [%g, %g] s (the first ' ...
        'sampled state is at %g s), re-orthonormalized every %g s, basis seed %d. Same sampled states as the density.'], ...
        cfg.ritz_scaling, ritz{1}.k_requested, R1.K, R1.warmup_s, R1.interval(1) - R1.warmup_s, R1.interval(1), ...
        R1.interval(1), R1.interval(2), min(R1.t), R1.lya_dt, R1.seed);
end
L{end + 1} = '';
L{end + 1} = ['| Condition | state dim. | k drawn | states | points | Re > 0 | max Re (1/s) | median Re (1/s) | ' ...
    'Re > 0 per state: min, median, max | states with any | ' ...
    'complex | max freq. (Hz) | median freq. of complex (Hz) | outside the window | x share of the subspace | ' ...
    'lambda_1 | lambda_k | positive exponents | moved by the other scaling: median, max (1/s) |'];
L{end + 1} = '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|';
scale_lines = {};
for i = 1:numel(cond_names)
    P = ritz{i}; R = D.ritz_by_cond{idx(i)};
    z = P.z(:);
    is_cx = abs(imag(z)) > 1e-9 * max(1, abs(z));
    f = abs(imag(z)) / (2 * pi);
    outside = real(z) < re_lim(1) | real(z) > re_lim(2) | imag(z) < im_lim(1) | imag(z) > im_lim(2);
    if any(is_cx); f_med = sprintf('%.3g', median(f(is_cx))); else; f_med = 'none'; end
    moved = 'not in the run';
    if ismember(other{1}, R.scalings)
        Po = ritz_points(R, cfg.ritz_k, other{1}, cond_names{i});
        dist = arrayfun(@(s) matched_distance(P.z(:, s), Po.z(:, s)), 1:size(P.z, 2));
        moved = sprintf('%.3g, %.3g', median(dist), max(dist));
    end
    per_state = sum(real(P.z) > 0, 1);
    L{end + 1} = sprintf(['| %s | %d | %d | %d | %d | %d (%.1f%%) | %+.3f | %+.3f | %d, %g, %d | %d of %d | %d (%.1f%%) | %.3g | %s | ' ...
        '%d (%.1f%%) | %.3f | %+.4f | %+.4f | %d | %s |'], titles{i}, P.N, P.k, numel(P.states), numel(z), ...
        nnz(real(z) > 0), 100 * mean(real(z) > 0), max(real(z)), median(real(z)), ...
        min(per_state), median(per_state), max(per_state), nnz(per_state > 0), numel(per_state), nnz(is_cx), 100 * mean(is_cx), ...
        max(f), f_med, nnz(outside), 100 * mean(outside), median(P.x_share), P.exponents(1), P.exponents(end), ...
        nnz(P.exponents > 0), moved); %#ok<AGROW>
    scale_lines{end + 1} = sprintf('- %s: %s', titles{i}, strjoin(arrayfun(@(k) sprintf('%s %.3g', ...
        R.scales.kinds{k}, R.scales.sd(k)), 1:numel(R.scales.kinds), 'UniformOutput', false), ', ')); %#ok<AGROW>
end
L{end + 1} = '';
L{end + 1} = ['Columns. Re > 0, complex and outside the window: counts of points, with their share. Re > 0 per state: ' ...
    'each sampled state gives k points; the smallest, median and largest number of them right of zero in one state, ' ...
    'and in how many states at least one is. A point is ' ...
    'complex when |Im| exceeds 1e-9 of its modulus (or of 1); frequency is |Im| / 2 pi. x share: the share of the ' ...
    'drawn subspace''s squared norm in dendritic coordinates, median over states (in the drawn norm). lambda_1, ' ...
    'lambda_k: the largest and the k-th exponent of this pass over its window. Moved by the other scaling: per ' ...
    'state, the largest distance between a drawn point and its partner in the other scaling''s set (greedy ' ...
    'pairing); median and maximum over states.'];
L{end + 1} = '';
L{end + 1} = '`type_sd` scales, the pooled standard deviation in time of each kind of variable over the sampled window, per condition:';
L{end + 1} = '';
L = [L, scale_lines];
L{end + 1} = '';
L = [L, ritz_diagnostics_section(D, idx, cond_names, titles, cfg)];
L{end + 1} = '';
if isfield(D, 'state_step') && isfield(R1, 'local_rates_before')
    L = [L, shared_pass_section(D, idx, titles, ritz)];
    L{end + 1} = '';
end
L{end + 1} = '## Caption paragraph, to adapt';
L{end + 1} = '';
trajectory = 'a baseline simulation';
if isfield(D, 'settings') && isfield(D.settings, 'trajectory'); trajectory = D.settings.trajectory; end
L{end + 1} = sprintf(['Gray: density of the eigenvalues of the Jacobian at %d states sampled along ' trajectory ' ' ...
    '(darker = more). Color: at the same states, the eigenvalues of the Jacobian restricted to the %d directions the ' ...
    'Lyapunov computation tracks (Ritz values), drawn as a density on the same bins; these are instantaneous growth ' ...
    'and rotation rates of the tracked directions, not eigenvalues of the Jacobian. Black ticks on the bottom edge: ' ...
    'the Lyapunov exponents of those %d directions, averaged over [%g, %g] s. Inset: the real axis at the scale of ' ...
    'the Ritz values, with a histogram of their real parts (top), the eigenvalue density and the Ritz values ' ...
    '(bottom), the exponent ticks and zero (dashed). The share of Ritz values right of zero, and in how many states ' ...
    'they occur, is in the first table.'], numel(ritz{1}.states), ritz{1}.k, ritz{1}.k, R1.interval(1), R1.interval(2));
L{end + 1} = '';
L{end + 1} = '## Reading the overlay';
L{end + 1} = '';
L{end + 1} = ['- **What the points are.** The density pools the eigenvalues of the whole Jacobian, every state direction ' ...
    'at once. The points are the eigenvalues of the same Jacobians restricted to the k directions the Lyapunov ' ...
    'computation tracks: the directions a perturbation contracts least along the trajectory.'];
L{end + 1} = ['- **In a stable adapting network the points mark the slow band.** Its least-contracting directions are ' ...
    'a band of slow single-neuron adaptation modes, about one per neuron. With k far below the network size all k ' ...
    'directions lie in that band (a small x share in the table), so the points gather near the slowest adaptation ' ...
    'rate (about -0.1 1/s) and say little about the network''s collective modes. No affordable k passes the band ' ...
    'at 500 neurons.'];
L{end + 1} = ['- **The points are not eigenvalues of the Jacobian, and a point right of zero is not an instability.** ' ...
    'They would be eigenvalues only if the tracked subspace were invariant under the instantaneous Jacobian. It is ' ...
    'carried by the dynamics instead, and the residual column of the diagnostics says how far from invariant it is ' ...
    '(0 = invariant; on the 500-neuron networks it is 0.5 to 0.9, and single Ritz pairs are no closer). So each ' ...
    'point is an instantaneous rate, with rotation, of the tracked directions among themselves: a point of the ' ...
    'numerical range of the Jacobian, bounded by the numerical abscissa of H. Points right of zero say that, at that ' ...
    'state and in that norm, some tracked directions are growing; the Lyapunov exponents (the ticks) are the time ' ...
    'averages and can all be negative. Compare "Re > 0", "local rates > 0" and "positive exponents".'];
L{end + 1} = ['- **Do not read the points against the density one by one.** The density is of true eigenvalues, which ' ...
    'do not depend on the norm; the points do. In the raw norm the number of points right of zero follows, state by ' ...
    'state, the number of eigenvalues of the full Jacobian right of zero (fig_ritz_by_state); in the type_sd norm ' ...
    'almost none remain on the multiple-timescale network although those eigenvalues are still there.'];
L{end + 1} = ['- **Single points do not average to single exponents.** Only sums agree: at each state the real parts of ' ...
    'the k points sum to the trace of the block, the sum of the top-k instantaneous local Lyapunov rates.'];
L{end + 1} = ['- **The norm matters.** Ritz values depend on how x, the adaptation and the depression variables are ' ...
    'scaled against one another; the last column says by how much here. The edge of a top-k subspace inside a ' ...
    'continuous band also depends on the initial basis.'];
L{end + 1} = sprintf(['- **The black ticks on the bottom edge** are the Lyapunov exponents of the same k directions from the same ' ...
    'pass: the averages over the sampled window, [%g, %g] s, of the local rates whose instantaneous values are the ' ...
    'diagonal of H. They are k numbers per condition, not one set per state, and they are real. At the panel''s ' ...
    'scale they overlap where the exponents are close; the inset separates them.'], R1.interval(1), R1.interval(2));
L{end + 1} = ['- **How the points are drawn.** `ritz_style` = `' cfg.ritz_style '`. The eigenvalue density is a log-scaled ' ...
    'count; a plain scatter shows every point as a dot, so a pile of a thousand points looks like one and a few dozen ' ...
    'strays look like many. `density` draws the Ritz values as a log-scaled count on the same bins, `marginal` adds a ' ...
    'linear histogram of the real parts. The counts in the first table are the measure.'];
txt = strjoin(L, newline);
end

function L = ritz_diagnostics_section(D, idx, cond_names, titles, cfg)
% What the points are: per condition and per leading block k, medians over
% the sampled states of the residuals, of the structure of H and of the
% comparison with the full Jacobian at the same states (ritz_diagnostics).
L = {};
L{end + 1} = '## What the points are: diagnostics';
L{end + 1} = '';
R1 = D.ritz_by_cond{idx(1)};
if ~isfield(R1, 'JQ_gram')
    L{end + 1} = 'This run was made before the residuals were saved (no JQ_gram); rerun the stage for the diagnostics.';
    return
end
L{end + 1} = ['A Ritz value is an eigenvalue of the Jacobian only if the tracked subspace is invariant under it. ' ...
    'Medians over the sampled states (ranges in brackets), `' cfg.ritz_scaling '` norm, for the leading k directions:'];
L{end + 1} = '';
L{end + 1} = ['| Condition | k | residual | pair residual | local rates > 0 | Ritz Re > 0 | max Re Ritz (1/s) | ' ...
    'num. abscissa of H | full J: max Re | full J: eigenvalues Re > 0 | full J: num. abscissa | place | ' ...
    'diag. norm | lower (frame turning) | upper of T (feed-forward) | coupling (1/s) |'];
L{end + 1} = '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|';
med = @(v) median(v, 'omitnan');
for i = 1:numel(cond_names)
    R = D.ritz_by_cond{idx(i)};
    for k = unique(min([10, 25, R.K], R.K))
        P = ritz_points(R, k, cfg.ritz_scaling, cond_names{i});
        S = P.state; F = P.full;
        place = (S.max_re - F.max_re) ./ (F.num_abscissa - F.max_re);
        L{end + 1} = sprintf(['| %s | %d | %.3f [%.3f, %.3f] | %.3f | %.1f%% | %.1f%% | %+.3f [%+.3f, %+.3f] | %+.3f | ' ...
            '%+.3f | %g [%g, %g] | %+.2f | %.3f | %.3g | %.3g | %.3g | %.3g |'], titles{i}, P.k, ...
            med(S.residual), min(S.residual), max(S.residual), med(P.pair_residual(:)), ...
            100 * mean(S.share_diag_pos), 100 * mean(S.share_ritz_pos), med(S.max_re), min(S.max_re), max(S.max_re), ...
            med(S.num_abscissa), med(F.max_re), med(F.n_positive), min(F.n_positive), max(F.n_positive), ...
            med(F.num_abscissa), med(place), med(S.norm_diag), med(S.norm_lower), med(S.norm_upper_T), ...
            med(S.coupling_TS)); %#ok<AGROW>
    end
end
L{end + 1} = '';
L{end + 1} = ['Columns. **residual**: ||J Q - Q H||_F / ||J Q||_F, 0 for a subspace invariant under the Jacobian, near 1 ' ...
    'when the Jacobian sends the tracked directions almost entirely out of the subspace. **pair residual**: the same for ' ...
    'one Ritz pair, ||J Q y - theta Q y|| / ||J Q y||, median over pairs and states; small only where a point is nearly ' ...
    'an eigenvalue of the Jacobian. **local rates > 0**: share of the diagonal entries of H (the instantaneous local ' ...
    'Lyapunov rates of the k directions) that are positive, against **Ritz Re > 0**, the share of its eigenvalues; a ' ...
    'gap between the two is spread made by the off-diagonal part, not local growth. **num. abscissa of H**: the largest ' ...
    'eigenvalue of (H + H'')/2, the largest instantaneous growth rate inside the subspace; it bounds every Ritz real ' ...
    'part. **full J**: the Jacobian at the same states: its largest eigenvalue real part, how many of its eigenvalues ' ...
    'have positive real part, and its numerical abscissa in this norm, which bounds everything. **place**: where the ' ...
    'largest Ritz real part sits between the full spectral abscissa (0) and the full numerical abscissa (1); below 0 it ' ...
    'is left of the rightmost eigenvalue. **diag. norm, lower, upper of T**: H = T + S in Gram-Schmidt order, T upper ' ...
    'triangular with the local rates on its diagonal, S skew (the turning of the frame); the 2-norm of the diagonal and ' ...
    'the Frobenius norms of the strictly lower triangle of H (which is that of S) and of the strictly upper triangle of ' ...
    'T. **coupling**: the largest sqrt(|T_ij S_ji|), the scale by which feed-forward weight times frame turning moves ' ...
    'an eigenvalue of H off its diagonal.'];
end

function L = shared_pass_section(D, idx, titles, ritz)
% FD-034: the states were sampled along the local_lyapunov stage's own run
% and H comes from the pass whose local rates Figure 3 plots. Three tables:
% the diagonal of H against the pass's interval rates, the first against the
% second half of the window, and the numbers by input step.
L = {};
L{end + 1} = '## The shared pass: states of the switching run, H from the pass whose local rates are plotted';
L{end + 1} = '';
L{end + 1} = sprintf(['The %d states are evenly spaced over [%g, %g] s of %s, and H = Q''JQ comes from the top-%d pass ' ...
    'whose local rates the expansion row of Figure 3 sums (the same exponents bit for bit; run_local_lyapunov). ' ...
    'The diagonal of H is therefore the instantaneous value of those rates, and at each state the real parts of ' ...
    'the Ritz values sum to its trace.'], numel(D.state_times), D.state_times(1), D.state_times(end), ...
    D.settings.trajectory, D.ritz_by_cond{idx(1)}.K);
if D.state_times(1) > D.lya_T_interval(1) + 1e-9
    L{end} = [L{end} sprintf([' The pass accumulates over [%g, %g] s after %g s of alignment over [%g, %g] s; the ' ...
        'states start later, at %g s.'], D.lya_T_interval(1), D.lya_T_interval(2), D.settings.lya_warmup, ...
        D.lya_T_interval(1) - D.settings.lya_warmup, D.lya_T_interval(1), D.state_times(1))];
end
L{end + 1} = '';
L{end + 1} = ['| Condition | sum Re(Ritz) - trace(H): max abs (1/s) | diag(H) vs interval rate: correlation | median abs diff (1/s) | ' ...
    '95th pct. abs diff (1/s) | interval rates before vs after: median abs diff (1/s) | ' ...
    'expansion (bits/s), instantaneous vs interval: correlation | median instantaneous | median interval | median abs diff |'];
L{end + 1} = '|---|---|---|---|---|---|---|---|---|---|';
S = cell(1, numel(idx));
for i = 1:numel(idx)
    R = D.ritz_by_cond{idx(i)}; P = ritz{i};
    st = P.states;
    dg = cell2mat(arrayfun(@(s) diag(R.H(:, :, s, P.sc))', st(:), 'UniformOutput', false));   % states x K
    bef = R.local_rates_before(st, :); aft = R.local_rates_after(st, :);
    avg = mean(cat(3, bef, aft), 3, 'omitnan');             % the interval rate around the state
    ok = isfinite(avg);
    ad = abs(dg(ok) - avg(ok));
    both = isfinite(bef) & isfinite(aft);
    exp_inst = sum(max(dg, 0), 2) / log(2);
    exp_int  = mean([sum(max(bef, 0), 2), sum(max(aft, 0), 2)], 2, 'omitnan') / log(2);   % NaN rows sum to 0: mask below
    exp_int(~any(isfinite(aft), 2)) = sum(max(bef(~any(isfinite(aft), 2), :), 0), 2) / log(2);
    exp_int(~any(isfinite(bef), 2)) = sum(max(aft(~any(isfinite(bef), 2), :), 0), 2) / log(2);
    tr_err = max(abs(sum(real(P.z), 1) - P.trace_k));
    L{end + 1} = sprintf('| %s | %.2g | %.3f | %.3g | %.3g | %.3g | %.3f | %.3g | %.3g | %.3g |', titles{i}, tr_err, ...
        corr_of(dg(ok), avg(ok)), median(ad), quantile(ad, 0.95), median(abs(bef(both) - aft(both))), ...
        corr_of(exp_inst, exp_int), median(exp_inst), median(exp_int), median(abs(exp_inst - exp_int))); %#ok<AGROW>
    S{i} = struct('dg', dg, 'exp_inst', exp_inst, 'exp_int', exp_int, 't', R.t(st), 'step', D.state_step(st));
end
L{end + 1} = '';
L{end + 1} = ['Columns. **diag(H) vs interval rate**: every tracked direction at every state; the interval rate is the ' ...
    'pass''s local rate averaged over the two re-orthonormalization intervals that end and start at the state (one where ' ...
    'the run has only one). They are not the same quantity (an instant against an average over the interval), and ' ...
    '**before vs after** gives the scale: how much the interval rate itself changes from one interval to the next. ' ...
    '**expansion**: the sum of the positive rates over the tracked directions divided by log 2, what the expansion row ' ...
    'plots, from the diagonal of H (instantaneous) and as plotted for the intervals ending and starting at the state (their mean).'];
L{end + 1} = '';
L{end + 1} = 'First against second half of the window (the basis is aligned before the window; an unconverged basis would show as a drift):';
L{end + 1} = '';
L{end + 1} = '| Condition | half | states | median Re(Ritz) (1/s) | Ritz Re > 0 | local rates > 0 | residual, median | full J: eigenvalues Re > 0, median | x share, median |';
L{end + 1} = '|---|---|---|---|---|---|---|---|---|';
for i = 1:numel(idx)
    P = ritz{i}; t = S{i}.t; mid = (D.state_times(1) + D.state_times(end)) / 2;
    halves = {t <= mid, t > mid}; label = {'first', 'second'};
    for h = 1:2
        L{end + 1} = group_row(sprintf('| %s | %s', titles{i}, label{h}), P, S{i}, halves{h}); %#ok<AGROW>
    end
end
L{end + 1} = '';
L{end + 1} = sprintf('By input step (a new input every %g s; a state at a switch instant is counted with the step that ends there):', ...
    D.settings.step_duration_s);
L{end + 1} = '';
L{end + 1} = ['| Condition | step | states | median Re(Ritz) (1/s) | Ritz Re > 0 | local rates > 0 | residual, median | ' ...
    'full J: eigenvalues Re > 0, median | x share, median | expansion, instantaneous: median (bits/s) |'];
L{end + 1} = '|---|---|---|---|---|---|---|---|---|---|';
for i = 1:numel(idx)
    for q = unique(S{i}.step(:))'
        in = S{i}.step(:) == q;
        L{end + 1} = [group_row(sprintf('| %s | %d', titles{i}, q), ritz{i}, S{i}, in), ...
            sprintf(' %.3g |', median(S{i}.exp_inst(in)))]; %#ok<AGROW>
    end
end
end

function row = group_row(lead, P, S, in)
% One row of the shared-pass tables: the states selected by `in`.
in = in(:)';
z = P.z(:, in); dg = S.dg(in, :);
res = NaN; if P.has_gram; res = median(P.state.residual(in), 'omitnan'); end
n_pos = NaN; if ~isempty(P.full); n_pos = median(P.full.n_positive(in)); end
row = sprintf('%s | %d | %+.3f | %.1f%% | %.1f%% | %.3f | %g | %.3f |', lead, nnz(in), median(real(z(:))), ...
    100 * mean(real(z(:)) > 0), 100 * mean(dg(:) > 0), res, n_pos, median(P.x_share(in)));
end

function c = corr_of(a, b)
% Pearson correlation of two samples (no toolbox).
a = a(:) - mean(a(:)); b = b(:) - mean(b(:));
c = (a' * b) / sqrt((a' * a) * (b' * b));
end

function d = matched_distance(a, b)
% Largest distance after pairing two sets of complex numbers greedily.
a = a(:); b = b(:); d = 0;
for k = 1:numel(a)
    [dmin, j] = min(abs(b - a(k)));
    d = max(d, dmin); b(j) = [];
end
end









