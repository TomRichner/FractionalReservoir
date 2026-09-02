function ax = plot_eigenvalue_heatmap_helper(ax, D, re_edges, im_edges, clim, use_log)
%PLOT_EIGENVALUE_HEATMAP_HELPER Render one eigenvalue-density panel.
%
%   ax = PLOT_EIGENVALUE_HEATMAP_HELPER(ax, D, re_edges, im_edges)
%   ax = PLOT_EIGENVALUE_HEATMAP_HELPER(ax, D, re_edges, im_edges, clim, use_log)
%
% Draws the smoothed density D as an image on the complex plane and overlays the
% Re = 0 stability line. Pass a shared clim across multiple panels for
% directly-comparable color scaling.
%
% Inputs:
%   ax        - target axes
%   D         - density matrix from compute_eigenvalue_density
%   re_edges  - real-axis bin edges
%   im_edges  - imag-axis bin edges
%   clim      - [lo hi] color limits (shared across panels)
%   use_log   - if true, display log10(1 + D) (default true)
%
% MOVED here from a static on SRNNModel2 (2026-09-02), with
% compute_eigenvalue_density. It draws into an axes from plain arrays and knows
% nothing about any model class.
%
% See also: compute_eigenvalue_density, fig_eig_heatmap

if nargin < 6 || isempty(use_log), use_log = true; end

re_centers = (re_edges(1:end-1) + re_edges(2:end)) / 2;
im_centers = (im_edges(1:end-1) + im_edges(2:end)) / 2;

if use_log
    plot_val = log10(1 + D);
else
    plot_val = D;
end

% D is (Re x Im); imagesc expects rows=y(Im), cols=x(Re) -> transpose
imagesc(ax, re_centers, im_centers, plot_val');
axis(ax, 'xy');
axis(ax, 'image');
colormap(ax, parula);
if nargin >= 5 && ~isempty(clim)
    caxis(ax, clim);
end

% Stability line at Re = 0
hold(ax, 'on');
yl = [im_centers(1), im_centers(end)];
plot(ax, [0, 0], yl, 'w--', 'LineWidth', 1.25);
hold(ax, 'off');

xlabel(ax, 'Re(\lambda)');
ylabel(ax, 'Im(\lambda)');
end
