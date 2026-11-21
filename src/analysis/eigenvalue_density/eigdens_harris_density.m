function H = eigdens_harris_density(alpha, f, mu_e, mu_i, sigma_e, sigma_i, r)
%EIGDENS_HARRIS_DENSITY The eigenvalue density of a sparse E-I random matrix, Harris et al. (2023) Eqs. 16, 18, 22 and 23.
%
%   H = eigdens_harris_density(alpha, f, mu_e, mu_i, sigma_e, sigma_i)
%   H = eigdens_harris_density(..., r)
%
% W = S .* (A D + u v'): connection probability alpha, a fraction f of
% excitatory columns with entries of mean mu_e and standard deviation
% sigma_e before the mask, the rest inhibitory (mu_i, sigma_i). The means
% and standard deviations are the paper's normalized ones times sqrt(N):
% pass them so that an i.i.d. matrix of unit radius has sigma = 1.
%
% Eq. 16   N sigma_sk^2 = alpha (1 - alpha) mu_k^2 + alpha sigma_k^2     k = e, i
% Eq. 18   R^2 = f N sigma_se^2 + (1 - f) N sigma_si^2
% Eq. 22   rho(z) = (1 / (2 pi N)) [ SP - DP H(DP |z|^2) ],  |z| <= R
% Eq. 23   H(x) = (x - Df N) / sqrt((x - Df N)^2 + N^2 (1 - Df^2))
% with the precisions P_sk = 1 / sigma_sk^2, SP their sum, DP = P_se - P_si,
% Df = 2 f - 1. In the variables v_k = N sigma_sk^2 (used here) N drops out.
% The density is uniform when v_e = v_i, i.e. (Eq. 20, g = 0) when
% (1 - alpha) mu_e^2 + sigma_e^2 = (1 - alpha) mu_i^2 + sigma_i^2.
%
% Output: struct with R, r, rho (the density, per unit area, integrating to
% 1 over the disk), rho_rel = rho pi R^2 (1 for a uniform disk), v_e, v_i,
% and g of Eq. 20.
%
% See also: eigdens_radial_theory (the same density from the general
% equation; eigdens_check_theory compares the two)

v_e = alpha * (1 - alpha) * mu_e^2 + alpha * sigma_e^2;
v_i = alpha * (1 - alpha) * mu_i^2 + alpha * sigma_i^2;
R = sqrt(f * v_e + (1 - f) * v_i);
if nargin < 7 || isempty(r)
    r = linspace(0, R, 400);
end
r = r(:)';
P_e = 1 / v_e; P_i = 1 / v_i;       % precisions divided by N
SP = P_e + P_i; DP = P_e - P_i; Df = 2 * f - 1;
x = DP .* r.^2;                     % (DP |z|^2) / N in the paper's variables
Hx = (x - Df) ./ sqrt((x - Df).^2 + (1 - Df^2));
rho = (SP - DP .* Hx) / (2 * pi);
rho(r > R) = 0;
H = struct('R', R, 'r', r, 'rho', rho, 'rho_rel', rho * pi * R^2, ...
    'v_e', v_e, 'v_i', v_i, 'g', 1 - v_i / v_e);
end
