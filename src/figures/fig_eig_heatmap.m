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
%   ritz_exponent_ticks  ticks along the bottom edge at the top-k Lyapunov
%                  exponents of the same pass (default true)
%   ritz_inset     an inset on the real axis near zero (ritz_inset_xlim,
%                  default [-1 0.5] 1/s) with the density recomputed at that
%                  scale: the slow band of an adapting network, near
%                  -0.1 1/s, cannot be told from zero on a window tens of
%                  1/s wide
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
    cfg.ritz_exponent_ticks (1,1) logical = true
    cfg.ritz_inset   (1,1) logical = false
    cfg.ritz_inset_xlim (1,2) double = [-1 0.5]   % real-axis range of the inset (1/s)
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
data_file = resolve_data_file(cfg.data_file, cfg.run_dir, ...
    {fullfile(cfg.run_dir, 'eig_heatmap')}, ...
    'eig_heatmap_data.mat', ...
    'Run run_eig_heatmap first');
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
        draw_ritz(ax_panels(i), ritz{i}, ritz_color(cond_names{i}), im_lim, cfg.ritz_exponent_ticks);
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
    heading{end + 1} = sprintf(['points: eigenvalues of the Jacobian projected on the top-%d Lyapunov ' ...
        'directions (Ritz values), %s norm'], ritz{1}.k_requested, strrep(cfg.ritz_scaling, '_', '\_'));
end
title(tl, heading, 'FontWeight', 'bold');

if cfg.ritz_overlay && cfg.ritz_inset
    % After the layout is complete (title, colorbar), so the plot boxes are final.
    for i = 1:n_cond
        draw_ritz_inset(fig, ax_panels(i), i, evals_by_cond{i}, ritz{i}, ritz_color(cond_names{i}), ...
            cfg.ritz_inset_xlim, scale_fn);
    end
end

if ~cfg.visible; set(fig, 'Visible', 'off'); end

%% --- Save -------------------------------------------------------------------
fig_tag = 'fig_eig_heatmap';
if cfg.ritz_overlay
    fig_tag = [fig_tag '_ritz'];
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
if k_used == R.K
    z = R.values(:, states, sc);
else
    z = complex(zeros(k_used, numel(states)));
    for s = 1:numel(states)
        z(:, s) = eig(R.H(1:k_used, 1:k_used, states(s), sc));
    end
end
trace_k = arrayfun(@(s) trace(R.H(1:k_used, 1:k_used, s, sc)), states);
P = struct('z', z, 'k', k_used, 'k_requested', k, 'K', R.K, 'N', R.N, 'sc', sc, ...
    'states', states, 'exponents', R.LE_spectrum(1:k_used, sc), 'trace_k', trace_k, ...
    'x_share', mean(R.x_fraction(states, 1:k_used, sc), 2));
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

function draw_ritz(ax, P, col, im_lim, with_ticks)
% The points, over the density and under nothing; the axes' limits are the
% density's and stay so (points outside are clipped and counted in the table).
xl = xlim(ax); yl = ylim(ax);
hold(ax, 'on');
% A thin white rim keeps a dark point readable on the dark part of the density.
scatter(ax, real(P.z(:)), imag(P.z(:)), 10, col, 'filled', 'MarkerFaceAlpha', 0.45, ...
    'MarkerEdgeColor', 'w', 'MarkerEdgeAlpha', 0.5, 'LineWidth', 0.25, 'Tag', 'ritz_values');
if with_ticks
    % One line object with NaN breaks: the Lyapunov exponents of the same
    % pass as ticks on the bottom edge. (One object of many points, so a
    % composite that looks for the two-point Re = 0 line does not find it.)
    lam = P.exponents(:)';
    len = 0.035 * diff(im_lim);
    xs = [lam; lam; nan(1, numel(lam))];
    ys = repmat([yl(1); yl(1) + len; NaN], 1, numel(lam));
    plot(ax, xs(:), ys(:), '-', 'Color', col, 'LineWidth', 0.75, 'Tag', 'ritz_exponent_ticks');
end
hold(ax, 'off');
xlim(ax, xl); ylim(ax, yl);
end

function draw_ritz_inset(fig, ax, panel, evals, P, col, re_lim_inset, scale_fn)
% A small axes in the lower right corner of the panel's plot box, zoomed on
% the real axis near zero: the density of the full Jacobian's eigenvalues at
% that scale (its own gray scale, darkest = its own maximum) and the Ritz
% values on it. Parented to the figure, tagged 'ritz_inset', with the panel
% index and its place within the plot box in UserData so a composite can
% put it back (fig_main3_stability).
rel = [0.50 0.10 0.47 0.25];                  % [x y w h] as fractions of the panel's plot box
box = plot_box_position(ax);
z = P.z(:);
in_x = real(z) >= re_lim_inset(1) & real(z) <= re_lim_inset(2);
im_half = max(0.5, 1.15 * max(abs(imag(z(in_x)))));
if isempty(im_half) || ~isfinite(im_half); im_half = 0.5; end
re_edges = linspace(re_lim_inset(1), re_lim_inset(2), 121);
im_edges = linspace(-im_half, im_half, 41);
Dn = scale_fn(compute_eigenvalue_density(evals, re_edges, im_edges, 1));
ai = axes(fig, 'Units', 'normalized', ...
    'Position', [box(1) + rel(1) * box(3), box(2) + rel(2) * box(4), rel(3) * box(3), rel(4) * box(4)]);
imagesc(ai, (re_edges(1:end - 1) + re_edges(2:end)) / 2, (im_edges(1:end - 1) + im_edges(2:end)) / 2, Dn');
axis(ai, 'xy'); colormap(ai, flipud(gray));
clim(ai, [0, max(max(Dn(:)), eps)]);
hold(ai, 'on');
st = manuscript_style();
plot(ai, [0 0], [-im_half im_half], '--', 'Color', st.zeroline_color, 'LineWidth', 1, 'Tag', 'ritz_inset_zero');
scatter(ai, real(z), imag(z), 8, col, 'filled', 'MarkerFaceAlpha', 0.45, ...
    'MarkerEdgeColor', 'w', 'MarkerEdgeAlpha', 0.5, 'LineWidth', 0.25, 'Tag', 'ritz_values_inset');
hold(ai, 'off');
set(ai, 'XLim', re_lim_inset, 'YLim', [-im_half im_half], 'Box', 'on', 'FontSize', 8, ...
    'TickDir', 'out', 'TickLength', [0.02 0.02], 'Layer', 'top', 'XTick', unique([re_lim_inset(1), 0, re_lim_inset(2)]), ...
    'YTick', [-im_half 0 im_half], 'YTickLabel', {sprintf('%.2g', -im_half), '0', sprintf('%.2g', im_half)}, ...
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
L{end + 1} = '';
L{end + 1} = ['| Condition | state dim. | k drawn | states | points | Re > 0 | max Re (1/s) | median Re (1/s) | ' ...
    'complex | max freq. (Hz) | median freq. of complex (Hz) | outside the window | x share of the subspace | ' ...
    'lambda_1 | lambda_k | positive exponents | moved by the other scaling: median, max (1/s) |'];
L{end + 1} = '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|';
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
    L{end + 1} = sprintf(['| %s | %d | %d | %d | %d | %d (%.1f%%) | %+.3f | %+.3f | %d (%.1f%%) | %.3g | %s | ' ...
        '%d (%.1f%%) | %.3f | %+.4f | %+.4f | %d | %s |'], titles{i}, P.N, P.k, numel(P.states), numel(z), ...
        nnz(real(z) > 0), 100 * mean(real(z) > 0), max(real(z)), median(real(z)), nnz(is_cx), 100 * mean(is_cx), ...
        max(f), f_med, nnz(outside), 100 * mean(outside), median(P.x_share), P.exponents(1), P.exponents(end), ...
        nnz(P.exponents > 0), moved); %#ok<AGROW>
    scale_lines{end + 1} = sprintf('- %s: %s', titles{i}, strjoin(arrayfun(@(k) sprintf('%s %.3g', ...
        R.scales.kinds{k}, R.scales.sd(k)), 1:numel(R.scales.kinds), 'UniformOutput', false), ', ')); %#ok<AGROW>
end
L{end + 1} = '';
L{end + 1} = ['Columns. Re > 0, complex and outside the window: counts of points, with their share. A point is ' ...
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
L{end + 1} = ['- **A point right of zero is not an instability.** The tracked directions are coupled through the ' ...
    'network, so H is not normal: its eigenvalues can leave the real axis and cross zero while every Lyapunov ' ...
    'exponent is negative. Compare the "Re > 0" and "positive exponents" columns.'];
L{end + 1} = ['- **The points need not fall on the density.** The tracked subspace is carried by the dynamics and is ' ...
    'not an invariant subspace of the instantaneous Jacobian. The points lie within its numerical range.'];
L{end + 1} = ['- **Single points do not average to single exponents.** Only sums agree: at each state the real parts of ' ...
    'the k points sum to the trace of the block, the sum of the top-k instantaneous local Lyapunov rates.'];
L{end + 1} = ['- **The norm matters.** Ritz values depend on how x, the adaptation and the depression variables are ' ...
    'scaled against one another; the last column says by how much here. The edge of a top-k subspace inside a ' ...
    'continuous band also depends on the initial basis.'];
L{end + 1} = ['- **Ticks on the bottom edge** are the top-k Lyapunov exponents of the same pass over its window ' ...
    '(at this scale they overlap where the exponents are close).'];
txt = strjoin(L, newline);
end

function d = matched_distance(a, b)
% Largest distance after pairing two sets of complex numbers greedily.
a = a(:); b = b(:); d = 0;
for k = 1:numel(a)
    [dmin, j] = min(abs(b - a(k)));
    d = max(d, dmin); b(j) = [];
end
end









