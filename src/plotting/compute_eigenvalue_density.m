function D = compute_eigenvalue_density(eigenvalues, re_edges, im_edges, sigma_bins)
%COMPUTE_EIGENVALUE_DENSITY 2-D Gaussian-smoothed eigenvalue density.
%
%   D = COMPUTE_EIGENVALUE_DENSITY(eigenvalues, re_edges, im_edges)
%   D = COMPUTE_EIGENVALUE_DENSITY(eigenvalues, re_edges, im_edges, sigma_bins)
%
% Bins the (real, imag) eigenvalue coordinates onto the grid defined by
% re_edges x im_edges, then smooths with a small Gaussian kernel. conv2 is used
% (not imgaussfilt) to avoid an Image Processing Toolbox dependency.
%
% Inputs:
%   eigenvalues - vector of (complex) eigenvalues
%   re_edges    - bin edges along the real axis (1 x nRe+1)
%   im_edges    - bin edges along the imag axis (1 x nIm+1)
%   sigma_bins  - Gaussian smoothing width in bins (default 1.5)
%
% Output:
%   D - smoothed density matrix, size (numel(re_edges)-1) x (numel(im_edges)-1)
%       rows index the real axis, columns index the imaginary axis.
%
% MOVED here from a static on SRNNModel2 (2026-09-02). It takes arrays and
% returns an array -- no model state, no class coupling -- and its only caller,
% fig_eig_heatmap, applies it to SRNNCellTypePairs data. The class prefix was an
% accident of where it happened to be written, and it was one of the references
% keeping the paper pipeline dependent on SRNNModel2.
%
% See also: plot_eigenvalue_heatmap_helper, fig_eig_heatmap, run_eig_heatmap

if nargin < 4 || isempty(sigma_bins), sigma_bins = 1.5; end

H = histcounts2(real(eigenvalues(:)), imag(eigenvalues(:)), re_edges, im_edges);

if sigma_bins > 0
    % Build a normalized separable Gaussian kernel (+/- 3 sigma).
    half = max(1, ceil(3 * sigma_bins));
    gv = exp(-((-half:half).^2) / (2 * sigma_bins^2));
    gv = gv / sum(gv);
    K = gv(:) * gv(:)';          % 2-D isotropic Gaussian
    D = conv2(H, K, 'same');
else
    D = H;
end
end
