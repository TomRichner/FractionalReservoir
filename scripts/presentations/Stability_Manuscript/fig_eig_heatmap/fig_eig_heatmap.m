% fig_eig_heatmap.m - Plot Jacobian-eigenvalue occupancy heatmaps (load + plot)
%
% Plotting half of the eig-heatmap figure. Loads the pooled eigenvalues computed
% by compute_eig_heatmap.m (eig_heatmap_data.mat) and renders a Gaussian-smoothed
% 2-D DENSITY over the complex plane for each of the four adaptation regimes --
% an "occupancy" heatmap showing how much time the eigenvalues spend in each
% region, and in particular to the RIGHT of the imaginary axis (Re > 0, locally
% unstable). Panels share axis limits and one colorbar (log density).
%
% Run compute_eig_heatmap.m first to (re)generate the data. This script only
% plots, so the figure look can be iterated without re-simulating.

close all; clear; clc;

setup_paths();
this_dir = fileparts(mfilename('fullpath'));

%% ---- Load computed eigenvalue data ----------------------------------------
data_file = fullfile(this_dir, 'eig_heatmap_data.mat');
if ~isfile(data_file)
    error('fig_eig_heatmap:NoData', ...
        'Missing %s. Run compute_eig_heatmap.m first to generate it.', data_file);
end
load(data_file, 'evals_by_cond', 'condition_titles', 'lle_by_cond', 'lle_window');
n_cond = numel(condition_titles);
if ~exist('lle_by_cond', 'var')
    error('fig_eig_heatmap:NoLLE', ...
        '%s has no LLE data. Re-run compute_eig_heatmap.m (it now computes Benettin LLE).', data_file);
end

%% ---- Heatmap / plotting parameters ----------------------------------------
grid_res   = 250;     % heatmap bins per axis
sigma_bins = 1.25;    % Gaussian smoothing width (bins)
keep_frac  = 0.999;   % fraction of eigenvalue density to keep inside the window

save_figs  = true;
save_types = {'fig', 'png', 'svg'};

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
    D_by_cond{i} = SRNNModel2.compute_eigenvalue_density( ...
        evals_by_cond{i}, re_edges, im_edges, sigma_bins);
    cmax = max(cmax, max(log10(1 + D_by_cond{i}(:))));  % shared log-density max
end
clim = [0, cmax];

%% ---- Assemble the 2x2 comparison figure -----------------------------------
fig = figure('Position', [200, 150, 900, 760]);
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax_panels = gobjects(n_cond, 1);
for i = 1:n_cond
    ax_panels(i) = nexttile(tl);
    SRNNModel2.plot_eigenvalue_heatmap_helper( ...
        ax_panels(i), D_by_cond{i}, re_edges, im_edges, clim, true);
    title(ax_panels(i), condition_titles{i}, 'FontWeight', 'normal', 'FontSize', 14);

    % Finite-time (Benettin) LLE over the last lle_window seconds, top-left.
    text(ax_panels(i), 0.03, 0.96, sprintf('LLE = %+.3f', lle_by_cond(i)), ...
        'Units', 'normalized', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
end

% Single shared colorbar for the whole layout (all panels share clim + colormap).
cb = colorbar(ax_panels(end));
cb.Layout.Tile = 'east';
cb.Label.String = 'log_{10}(1 + eigenvalue density)';

title(tl, {'Jacobian eigenvalue occupancy across adaptation regimes', ...
    sprintf('LLE = finite-time Benettin exponent over the last %g s', lle_window)}, ...
    'FontWeight', 'bold');

%% ---- Save -----------------------------------------------------------------
if save_figs
    save_some_figs_to_folder_2(this_dir, 'fig_eig_heatmap', fig.Number, save_types);
    fprintf('Saved figure to %s\n', this_dir);
end
