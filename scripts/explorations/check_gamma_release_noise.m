% CHECK_GAMMA_RELEASE_NOISE  Numerical checks behind section 12 of
% docs/EquationsParametersDocs/chained_MTS/Equations_chained_MTS.md, the
% proposed release noise written with a gamma process (not in the model code).
%
% Four parts, all serial, a few minutes in total:
%   1. Gamma increments: mean, variance, and two half steps against one full
%      step.
%   2. A linear first-order filter driven by matched Poisson and gamma
%      inputs: equal mean and variance, skewness and excess kurtosis.
%   3. The kick derivative kappa = dQ/dk of the unit-scale gamma quantile at
%      a fixed uniform: mean, variance, tail and bound against the shape per
%      step, and its small-shape limit.
%   4. A toy network of the simplified equations (S2), (S3), (S5), (S6) with
%      gamma release and no Wiener term (additive Wiener noise adds nothing
%      to any tangent), written out here and not taken from the model class:
%      lambda_1 from a Benettin run with shared uniforms against lambda_1
%      from the three tangent dynamics of section 12.5.
%
% Needs the Statistics Toolbox (gamrnd, poissrnd). Changes no model code and
% writes no file. Results are printed.

rng(20261114, 'twister');

%% 1. Gamma increments
fprintf('\n== 1. Gamma increments ==\n');
Lam = 8;            % shape per unit time (1/s)
q = 0.0125;         % scale (s)
dt = 0.0025;        % step (s)
n = 4e6;
full_step = gamrnd(Lam * dt, q, n, 1);
half_sum = gamrnd(Lam * dt / 2, q, n, 1) + gamrnd(Lam * dt / 2, q, n, 1);
fprintf('shape per step %.3g\n', Lam * dt);
fprintf('mean      : full %.4e  two halves %.4e  theory %.4e\n', ...
    mean(full_step), mean(half_sum), Lam * q * dt);
fprintf('variance  : full %.4e  two halves %.4e  theory %.4e\n', ...
    var(full_step), var(half_sum), Lam * q^2 * dt);
fprintf('3rd cumul.: full %.4e  two halves %.4e  theory %.4e\n', ...
    mean((full_step - mean(full_step)).^3), ...
    mean((half_sum - mean(half_sum)).^3), 2 * Lam * q^3 * dt);
[~, p_ks, ks_stat] = kstest2(full_step(1:5e5), half_sum(1:5e5));
fprintf('two-sample KS, 5e5 each: statistic %.2e, p = %.2f\n', ks_stat, p_ks);
fprintf('min increment %.3g (never negative: %d)\n', min(full_step), ...
    all(full_step >= 0));

%% 2. Linear filter driven by matched Poisson and gamma inputs
fprintf('\n== 2. tau dy = -y dt + dL, L Poisson (q dP) or gamma ==\n');
tau = 0.1;
dt = tau / 100;
n = 2e7;
decay = exp(-dt / tau);
fprintf(['events/tau   mean P     mean G     var P      var G      ', ...
    'skew P  skew G  ratio  (theory P, G)    exkurt P  exkurt G  (theory)\n']);
for events_per_tau = [1 10 100]
    Lam = events_per_tau / tau;
    q = 1 / events_per_tau;                 % mean drive Lam q = 1/tau
    yP = filter(1, [1 -decay], q * poissrnd(Lam * dt, n, 1) / tau);
    yG = filter(1, [1 -decay], gamrnd(Lam * dt, q, n, 1) / tau);
    yP = yP(2000:end);
    yG = yG(2000:end);
    sk_th = (2^1.5 / 3) / sqrt(events_per_tau);
    fprintf(['%8g   %.4f   %.4f   %.3e  %.3e  %.3f   %.3f   %.2f   ', ...
        '(%.3f, %.3f)   %.3f     %.3f     (%.3f, %.3f)\n'], events_per_tau, ...
        mean(yP), mean(yG), var(yP), var(yG), skewness(yP), skewness(yG), ...
        skewness(yG) / skewness(yP), sk_th, 2 * sk_th, ...
        kurtosis(yP) - 3, kurtosis(yG) - 3, ...
        1 / events_per_tau, 6 / events_per_tau);
end
clear yP yG full_step half_sum

%% 3. The kick derivative kappa = dQ/dk at a fixed uniform
% Q(U; k) is the quantile of Gamma(shape k, scale 1). At fixed U,
% dQ/dk = -(dF/dk) / f, with F the distribution function and f the density.
% Moments are integrals over x = Q with weight f(x), on a grid uniform in
% log x. Below x = exp(-700) both x and kappa are below 1e-300 and add
% nothing to any moment.
fprintf('\n== 3. kappa = dQ/dk of the unit-scale gamma quantile ==\n');
m_fun = @(x) x .* exp(x) .* expint(x);          % small-shape limit of k*kappa
c0 = integral(@(x) x .* exp(x) .* expint(x).^2, 1e-12, 500);
fprintf('small-shape limit: k*kappa -> m(x) = x e^x E1(x); ');
fprintf('integral E1 = %.6f; c0 = int x e^x E1^2 dx = %.4f\n', ...
    integral(@(x) expint(x), 1e-12, Inf), c0);
fprintf(['   shape k    E[kappa]   CV     c = k*Var   max k*kappa  ', ...
    'P(kappa>10)  P(kappa>100)  kurtosis   max|k*kappa - m|\n']);
for k = [1e-4 1e-3 1e-2 0.1 0.3 1 3 10 100]
    x_hi = k + 40 * sqrt(k) + 60;
    t = linspace(-700, log(x_hi), 6e5).';
    x = exp(t);
    h = 1e-4 * k;
    f = exp((k - 1) .* t - x - gammaln(k));
    dP = (gammainc(x, k + h) - gammainc(x, k - h)) / (2 * h);
    dS = (gammainc(x, k + h, 'upper') - gammainc(x, k - h, 'upper')) / (2 * h);
    lower_half = gammainc(x, k) < 0.5;
    kappa = dS ./ f;
    kappa(lower_half) = -dP(lower_half) ./ f(lower_half);
    w = f .* x;                                   % density in t
    ok = isfinite(kappa) & w > 0;
    kappa(~ok) = 0;
    w(~ok) = 0;
    m1 = trapz(t, w .* kappa);
    m2 = trapz(t, w .* kappa.^2);
    m4 = trapz(t, w .* (kappa - m1).^4);
    v = m2 - m1^2;
    in_bulk = x > 1e-6 & x < 30;
    fprintf(['%10.4g   %.5f   %6.2f   %.4f      %.4f      %.2e    %.2e   ', ...
        '%9.3g   %.2e\n'], k, m1, sqrt(v) / m1, k * v, ...
        max(k * kappa(x < k + 8 * sqrt(k) + 40)), ...
        trapz(t, w .* (kappa > 10)), trapz(t, w .* (kappa > 100)), ...
        m4 / v^2, max(abs(k * kappa(in_bulk) - m_fun(x(in_bulk)))));
end

% Monte Carlo cross-check of the mean and variance through the quantile
% function itself (gammaincinv), by a central difference in the shape.
for k = [0.1 1]
    U = rand(2e6, 1);
    h = 1e-4 * k;
    kap = (gammaincinv(U, k + h) - gammaincinv(U, k - h)) / (2 * h);
    fprintf('Monte Carlo, k = %g: E[kappa] = %.4f, k*Var = %.4f\n', ...
        k, mean(kap), k * var(kap));
end

% The range of shape per step the model would use. Matching the stationary
% spread of x under today's Wiener noise, sigma_u / sqrt(2 tau_d), with the
% paper network's 100 inputs of mean size 7F and spread 1.5F.
tau_d = 0.1; tau_rel = 0.25; sigma_u = 0.025; b_typ = 0.7;
F = 1 / sqrt(500 * 0.2 * 1.8);
row_w2 = 100 * (7^2 + 1.5^2) * F^2;
rs_over_nv = (sigma_u^2 / (2 * tau_d)) * 2 * tau_d / (row_w2 * tau_rel * b_typ);
fprintf(['matching sigma_u = %.3g: sum_j w_ij^2 = %.1f, so r_s / n_v = %.2e ', ...
    '(at b1 = %.1f)\n'], sigma_u, row_w2, rs_over_nv, b_typ);
for r_s = [0.005 0.02 0.1]
    n_v = r_s / rs_over_nv;
    fprintf('  r_s = %.3f -> n_v = %5.0f; shape per step at fs = 200, 400, 800:', ...
        r_s, n_v);
    fprintf(' %.4f', r_s * n_v * b_typ ./ ([200 400 800] * tau_rel));
    fprintf('\n');
end

%% 4. Toy network: lambda_1 by shared-uniform Benettin and by three tangents
fprintf('\n== 4. Toy network, N = 40, chain, no facilitation ==\n');
p = struct();
p.N = 40;
p.tau_d = 0.1; p.a0 = 0.2; p.c = 0.5;
p.tau_rec1 = 2; p.tau_rec2 = 5; p.tau_rel = 0.25; p.n2 = 2;
p.dt = 0.0025;
rng(7, 'twister');
W0 = randn(p.N) / sqrt(p.N);
% name, gain, K, r_s, n_v. The stable cases use one fast SFA timescale so
% that the leading exponent is set by the pools and not by a 10-s SFA mode.
% At high gain the network is chaotic under weak release only.
cases = { ...
    'high gain, weak  ', 30, 3, 0.02, 100; ...
    'high gain, medium', 30, 3, 0.05,  30; ...
    'high gain, strong', 30, 3, 0.10,  10; ...
    'low gain,  weak  ',  6, 1, 0.02, 100; ...
    'low gain,  strong',  6, 1, 0.10,  10};
t_acc = 1500;
fprintf('accumulation %g s after 20 s, step %g s; se from 20 blocks\n', ...
    t_acc, p.dt);
for ic = 1:size(cases, 1)
    p.W = cases{ic, 2} * W0;
    p.K = cases{ic, 3};
    if p.K == 3
        p.tau_a = [0.25 1.581 10];
    else
        p.tau_a = 0.25;
    end
    p.r_s = cases{ic, 4};
    p.n_v = cases{ic, 5};
    o = toy_run(p, 20, t_acc, 100 + ic);
    fprintf(['%s gain %g, K = %d, r_s = %.2f, n_v = %d: shape/step %.4f, ', ...
        'caps %d, mean rate %.3f, mean b1 %.3f, n_v*b1 %.1f\n'], ...
        cases{ic, 1}, cases{ic, 2}, p.K, p.r_s, p.n_v, o.shape_mean, ...
        o.n_caps, o.mean_r, o.mean_b1, p.n_v * o.mean_b1);
    fprintf(['      lambda_1: Benettin %+.5f | (c) pathwise %+.5f | ', ...
        '(b) mean %+.5f | (a) replay %+.5f\n'], ...
        o.lam_ben, o.lam(3), o.lam(2), o.lam(1));
    fprintf(['      (c) - Benettin %+.1e | (c) - (b) %+.5f +- %.5f | ', ...
        '(b) - (a) %+.5f +- %.5f\n'], o.lam(3) - o.lam_ben, ...
        o.lam(3) - o.lam(2), o.se_cb, o.lam(2) - o.lam(1), o.se_ba);
end

%% Local functions

function o = toy_run(p, t_warm, t_acc, seed)
% One run of the toy network with gamma release. The fiducial path, a
% Benettin copy driven by the same uniforms, and three tangent vectors
% (kick derivative 0, its mean, and the pathwise one) advance together.
    rng(seed, 'twister');
    N = p.N; K = p.K; dt = p.dt;
    q = p.tau_rel / p.n_v;
    ix = 1:N; ia = N + (1:N * K); ib1 = N + N * K + (1:N); ib2 = ib1 + N;
    z = [0.5 * randn(N, 1); zeros(N * K, 1); ones(2 * N, 1)];
    n_warm = round(t_warm / dt); n_acc = round(t_acc / dt);
    d0 = 1e-7;
    v0 = randn(numel(z), 1); v0 = v0 / norm(v0);
    V = repmat(v0, 1, 3);                 % (a), (b), (c)
    zb = z + d0 * v0;                     % Benettin copy
    logs = zeros(n_acc, 3); logs_ben = zeros(n_acc, 1);
    n_caps = 0; sum_r = 0; sum_b1 = 0; sum_shape = 0;
    for it = 1:(n_warm + n_acc)
        U = rand(N, 1);
        % fiducial step
        [zdet, k1, zp] = heun(z, p);
        shape = p.r_s * p.n_v * max(z(ib1), 0) * dt / p.tau_rel;
        gam = q * gammaincinv(U, shape);
        cap = p.tau_rel * max(zdet(ib1), 0);
        capped = gam > cap;
        dG = min(gam, cap);
        z_new = zdet;
        z_new(ix) = z_new(ix) + p.W * dG / p.tau_d;
        z_new(ib1) = z_new(ib1) - dG / p.tau_rel;
        % pathwise kick derivative d(dG)/d(b1) at fixed U
        hk = 1e-4 * shape;
        kappa = (gammaincinv(U, shape + hk) - gammaincinv(U, shape - hk)) ./ (2 * hk);
        kappa(~isfinite(kappa)) = 0;
        dG_db = p.r_s * dt * kappa;
        % tangents: derivative of the Heun map, then the kick's derivative
        J1 = jac_times(z, V, p);
        J2 = jac_times(zp, V + dt * J1, p);
        Vdet = V + (dt / 2) * (J1 + J2);
        kick = [zeros(N, 1), p.r_s * dt * V(ib1, 2), dG_db .* V(ib1, 3)];
        kick(capped, 3) = p.tau_rel * Vdet(ib1(capped), 3);
        V = Vdet;
        V(ix, :) = V(ix, :) + p.W * kick / p.tau_d;
        V(ib1, :) = V(ib1, :) - kick / p.tau_rel;
        % Benettin copy: same uniforms, its own pool
        [zbdet, ~, ~] = heun(zb, p);
        shape_b = p.r_s * p.n_v * max(zb(ib1), 0) * dt / p.tau_rel;
        dGb = min(q * gammaincinv(U, shape_b), p.tau_rel * max(zbdet(ib1), 0));
        zb = zbdet;
        zb(ix) = zb(ix) + p.W * dGb / p.tau_d;
        zb(ib1) = zb(ib1) - dGb / p.tau_rel;
        z = z_new;
        % renormalize
        nv = sqrt(sum(V.^2, 1));
        V = V ./ nv;
        dist = norm(zb - z);
        zb = z + (d0 / dist) * (zb - z);
        if it > n_warm
            j = it - n_warm;
            logs(j, :) = log(nv);
            logs_ben(j) = log(dist / d0);
            n_caps = n_caps + nnz(capped);
            hh = z(ix) - p.a0 - (p.c / K) * sum(reshape(z(ia), N, K), 2);
            sum_r = sum_r + mean(phi(hh));
            sum_b1 = sum_b1 + mean(z(ib1));
            sum_shape = sum_shape + mean(shape);
        end
        if it == n_warm
            zb = z + d0 * V(:, 3);       % align the copy with tangent (c)
        end
    end
    T = n_acc * dt;
    o.lam = sum(logs, 1) / T;
    o.lam_ben = sum(logs_ben) / T;
    nb = 20; len = floor(n_acc / nb);
    d = logs(1:nb * len, 3) - logs(1:nb * len, 2);
    o.se_cb = std(sum(reshape(d, len, nb), 1) / (len * dt)) / sqrt(nb);
    d = logs(1:nb * len, 2) - logs(1:nb * len, 1);
    o.se_ba = std(sum(reshape(d, len, nb), 1) / (len * dt)) / sqrt(nb);
    o.n_caps = n_caps;
    o.mean_r = sum_r / n_acc; o.mean_b1 = sum_b1 / n_acc;
    o.shape_mean = sum_shape / n_acc;
end

function [z_new, k1, zp] = heun(z, p)
    k1 = rhs(z, p);
    zp = z + p.dt * k1;
    z_new = z + (p.dt / 2) * (k1 + rhs(zp, p));
end

function r = phi(h)
    r = 1 ./ (1 + exp(-4 * h));
end

function dz = rhs(z, p)
% The deterministic part: (S2), (S3), (S5), (S6) and the drift of (S16).
    N = p.N; K = p.K;
    x = z(1:N); a = reshape(z(N + (1:N * K)), N, K);
    b1 = z(N + N * K + (1:N)); b2 = z(2 * N + N * K + (1:N));
    r = phi(x - p.a0 - (p.c / K) * sum(a, 2));
    refill = b2 .* (1 - b1) / p.tau_rec1;
    dx = (-x + p.W * (r .* b1)) / p.tau_d;
    da = (r - a) ./ p.tau_a;
    db1 = refill - b1 .* r / p.tau_rel;
    db2 = (1 - b2) / p.tau_rec2 - refill / p.n2;
    dz = [dx; da(:); db1; db2];
end

function JV = jac_times(z, V, p)
% The Jacobian of rhs times the columns of V (section 10, equation (20)).
    N = p.N; K = p.K; nc = size(V, 2);
    x = z(1:N); a = reshape(z(N + (1:N * K)), N, K);
    b1 = z(N + N * K + (1:N)); b2 = z(2 * N + N * K + (1:N));
    h = x - p.a0 - (p.c / K) * sum(a, 2);
    r = phi(h); dphi = 4 * r .* (1 - r);
    vx = V(1:N, :); va = reshape(V(N + (1:N * K), :), N, K, nc);
    vb1 = V(N + N * K + (1:N), :); vb2 = V(2 * N + N * K + (1:N), :);
    dr = dphi .* (vx - (p.c / K) * reshape(sum(va, 2), N, nc));
    drefill = (vb2 .* (1 - b1) - b2 .* vb1) / p.tau_rec1;
    dvx = (-vx + p.W * (b1 .* dr + r .* vb1)) / p.tau_d;
    dva = (reshape(dr, N, 1, nc) - va) ./ p.tau_a;
    dvb1 = drefill - (vb1 .* r + b1 .* dr) / p.tau_rel;
    dvb2 = -vb2 / p.tau_rec2 - drefill / p.n2;
    JV = [dvx; reshape(dva, N * K, nc); dvb1; dvb2];
end
