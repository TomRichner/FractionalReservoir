function edge = fresp_plane_edge(R_plane, sigma, w, tau_d, sigma_min)
%FRESP_PLANE_EDGE The rightmost point of the region R(z) >= |1 + tau_d z| of the complex plane: the predicted abscissa of the bulk.
%
%   edge = fresp_plane_edge(R_plane, sigma, w, tau_d)
%   edge = fresp_plane_edge(R_plane, sigma, w, tau_d, sigma_min)
%
% With the slow variables free the bulk of the full Jacobian's network
% eigenvalues is where the predicted radius R(z) (fresp_radius of the gains
% at z) is at least |1 + tau_d z|. This returns the largest real part of a
% grid point in that region, NaN when the region holds no grid point.
%
%   R_plane    n_sigma x n_w: R at z = sigma(a) + i w(b)
%   sigma, w   the grid (1/s, rad/s)
%   tau_d      the dendritic time constant (s)
%   sigma_min  grid points left of this are ignored (default -8 1/s): the
%              far slow band, where the gains have their own poles, is not
%              the bulk
%
% See also: fresp_run_one, fig_freqresp_plane

if nargin < 5 || isempty(sigma_min); sigma_min = -8; end
[SG, WW] = ndgrid(sigma, w);
inside = double(R_plane) >= abs(1 + tau_d * complex(SG, WW));
inside(sigma < sigma_min, :) = false;
rows = find(any(inside, 2));
edge = NaN;
if ~isempty(rows); edge = sigma(rows(end)); end
end
