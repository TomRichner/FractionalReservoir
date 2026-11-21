function [z, res, n_iter] = fresp_outlier_root(S, params, ws, z0)
%FRESP_OUTLIER_ROOT The population mode of the full Jacobian, from the mean part alone: the root of K(z) = 1 + tau_d z nearest z0.
%
%   [z, res, n_iter] = fresp_outlier_root(S, params, ws, z0)
%
% The mean part of M(z) is the C x C matrix K(z)(a, b) = mu(a, b) sum_{j in b} d(a, j; z)
% (fresp_radius). Where an eigenvalue kappa(z) of it lies outside the
% bulk, M(z) has an eigenvalue near it, so the full Jacobian has an
% eigenvalue near the root of
%
%   kappa(z) = 1 + tau_d z .
%
% That is the frequency-dependent version of "the outlier sits at the
% eigenvalue of K_eff": here the gains carry phase, and the root is in
% general complex even when the frozen outlier is real. Newton's method on
% the complex function h(z) = kappa(z) - 1 - tau_d z (analytic away from
% the neurons' own poles), with a finite-difference derivative; at each
% step kappa is the eigenvalue of K(z) nearest 1 + tau_d z.
%
%   z0   starting point (the continuation the outlier_spectra stage found,
%        or a guess)
%   z    the root; res = |h(z)|; n_iter the iterations used (50 at most)
%
% See also: fresp_radius, fresp_run_one

tau_d = params.tau_d;
z = z0;
for n_iter = 1:50
    hz = h(z);
    step = 1e-6 * (1 + abs(z));
    dh = (h(z + step) - hz) / step;
    dz = -hz / dh;
    z = z + dz;
    if abs(dz) < 1e-12 * (1 + abs(z)); break; end
end
res = abs(h(z));

    function v = h(zq)
        [~, K] = fresp_radius(fresp_gain(S, params, zq), ws);
        target = 1 + tau_d * zq;
        [~, i] = min(abs(K - target));
        v = K(i) - target;
    end
end
