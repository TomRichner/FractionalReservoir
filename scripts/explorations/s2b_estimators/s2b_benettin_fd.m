function [t_lya, loc] = s2b_benettin_fd(M, seg, v0, d0)
%S2B_BENETTIN_FD The stage's Benettin estimator, re-implemented line for line
% from SRNNCellTypePairs.benettin_algorithm_internal so that the initial
% perturbation can be chosen: a FINITE perturbation of size d0 (the class uses
% 1e-3), the perturbed state re-integrated with ode45 (the model's options)
% over each segment and compared with the stored trajectory.
if nargin < 4; d0 = 1e-3; end
dec = round(seg * M.fs);  dt = 1 / M.fs;  tau = dec * dt;
idx = 1:dec:(numel(M.t) - dec);
t_lya = M.t(idx);  loc = zeros(numel(idx), 1);
p = d0 * v0(:) / norm(v0);
for k = 1:numel(idx)
    [~, Sp] = ode45(M.rhs, t_lya(k) + (0:dt:tau), M.S(idx(k), :)' + p, M.ode_opts);
    delta = Sp(end, :)' - M.S(idx(k) + dec, :)';
    d = norm(delta);
    loc(k) = log(d / d0) / tau;
    p = d0 * delta / d;
end
end
