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
% PLOTTING HALF ONLY. run_eig_heatmap does the sampling and writes
% eig_heatmap_data.mat, so the look can be iterated without re-simulating.
%
% See also: run_eig_heatmap, manuscript_style

arguments
    cfg.data_file   (1,:) char    = ''    % '' -> search run_dir, then data/eig_heatmap
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''    % the run whose eig_heatmap data to plot
    cfg.preset_name (1,:) char    = ''    % unused; the preset is recorded in the .mat
    cfg.density_scale (1,:) char {mustBeMember(cfg.density_scale, {'log', 'loglog'})} = 'log'
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
evals_by_cond    = D.evals_by_cond;
condition_titles = D.condition_titles;
lle_by_cond      = D.lle_by_cond;
lle_window       = D.lle_window;
n_cond           = numel(condition_titles);

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
clim = [0, cmax];

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
fig = figure('Position', [200, 150, 450*n_cols, 380*n_rows]);
tl  = tiledlayout(fig, n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

ax_panels = gobjects(n_cond, 1);
for i = 1:n_cond
    ax_panels(i) = nexttile(tl);
    plot_eigenvalue_heatmap_helper( ...
        ax_panels(i), D_by_cond{i}, re_edges, im_edges, clim, false);  % already scaled above
    title(ax_panels(i), condition_titles{i}, 'FontWeight', 'normal', 'FontSize', 14);

    % Finite-time (Benettin) LLE over the last lle_window seconds, top-left.
    text(ax_panels(i), 0.03, 0.96, sprintf('LLE = %+.3f', lle_by_cond(i)), ...
        'Units', 'normalized', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
end

% Single shared colorbar for the whole layout (all panels share clim + colormap).
cb = colorbar(ax_panels(end));
cb.Layout.Tile = 'east';
cb.Label.String = cb_label;

title(tl, {'Jacobian eigenvalue occupancy across adaptation regimes', ...
    sprintf('LLE = finite-time Benettin exponent over the last %g s', lle_window)}, ...
    'FontWeight', 'bold');

if ~cfg.visible; set(fig, 'Visible', 'off'); end

%% --- Save -------------------------------------------------------------------
fig_tag = 'fig_eig_heatmap';
out = struct('figs', fig, 'files', {{}}, 'source', data_file);
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);


end
end









