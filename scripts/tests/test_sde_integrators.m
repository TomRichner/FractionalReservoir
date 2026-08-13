% test_sde_integrators.m
% Verify src/sde_fixed_step.m: the three additive-noise schemes, their measured
% strong convergence orders, and the absolute-time noise indexing that Benettin
% depends on.
%
% The convergence measurement is the point of this file. A mistyped SRA1
% coefficient does not blow up -- it quietly degrades the scheme to strong
% order 1.0, which still converges and still looks plausible. The only thing
% that catches it is measuring the order, so that check is the acceptance
% criterion for SRA1 rather than a nice-to-have.
%
% Everything here runs on a small self-contained test system rather than on a
% model class, so a failure points at the integrator and nothing else.
%
% See also: sde_fixed_step, ode_rk4

fprintf('=== Testing sde_fixed_step ===\n\n');
all_passed = true;
rng(7);

%% A small, smooth, mildly nonlinear test system
n_sys = 4;
A = [-0.3 0.8 0.0 0.2; -0.6 -0.2 0.5 0.0; 0.1 -0.7 -0.4 0.6; 0.3 0.0 -0.5 -0.1];
tau = 0.1;
f = @(~, y) (-y + A * tanh(y) + 0.3) / tau;
y0 = [0.2; -0.1; 0.05; 0.3];

%% At sigma = 0 each scheme reduces to its deterministic parent, exactly
T_det = 0.2; fs_det = 800; t_det = (0:1/fs_det:T_det)';

[~, Y_sra1] = sde_fixed_step(f, t_det, y0, [], [], 'sra1');
Y_rk2 = ref_rk2(f, t_det, y0, 3/4);      % node 3/4, weights 1/3, 2/3
all_passed = check('sigma = 0: SRA1 equals RK2 with node 3/4 exactly', ...
    same(Y_sra1, Y_rk2)) && all_passed;

[~, Y_heun] = sde_fixed_step(f, t_det, y0, [], [], 'heun');
Y_trap = ref_trapezoid(f, t_det, y0);
all_passed = check('sigma = 0: Heun equals the trapezoidal corrector exactly', ...
    same(Y_heun, Y_trap)) && all_passed;

[~, Y_eul] = sde_fixed_step(f, t_det, y0, [], [], 'euler');
Y_fe = ref_forward_euler(f, t_det, y0);
all_passed = check('sigma = 0: Euler equals forward Euler exactly', ...
    same(Y_eul, Y_fe)) && all_passed;

% A noise struct whose sigma is zero must behave like no noise struct at all.
noise_zero = make_noise(zeros(n_sys, numel(t_det) - 1), ...
    zeros(n_sys, numel(t_det) - 1), 0, fs_det, 0, 1:n_sys);
[~, Y_sra1_z] = sde_fixed_step(f, t_det, y0, [], noise_zero, 'sra1');
all_passed = check('sigma = 0 with a noise struct matches no noise struct', ...
    same(Y_sra1_z, Y_sra1)) && all_passed;

%% Measured strong convergence order
% Reference: SRA1 on the finest grid. Coarse runs use the SAME Brownian path,
% rebuilt exactly by combining adjacent fine increments (see coarsen_noise).
T_conv = 0.25;
fs_fine = 51200;
sigma = 1.5;                 % large enough that noise error, not drift, dominates
n_paths = 16;
% Coarse step = refinement * h_fine. The smallest refinement sets how far the
% tested steps sit from the reference: too close and the reference's own error
% correlates with the tested run's (same scheme, same path) and cancels,
% flattering the apparent order. 8x separation keeps that below the noise.
refinements = [64 32 16 8];

fprintf('  measuring strong order over %d paths...\n', n_paths);
err = zeros(numel(refinements), 3, n_paths);
for p = 1:n_paths
    nsteps_fine = round(T_conv * fs_fine);
    xi1 = randn(n_sys, nsteps_fine);
    xi2 = randn(n_sys, nsteps_fine);

    nz_fine = make_noise(xi1, xi2, 0, fs_fine, sigma, 1:n_sys);
    t_fine = (0:1/fs_fine:T_conv)';
    [~, Y_ref] = sde_fixed_step(f, t_fine, y0, [], nz_fine, 'sra1');
    y_ref_end = Y_ref(end, :).';

    for ri = 1:numel(refinements)
        m = refinements(ri);
        [xi1_c, xi2_c] = coarsen_noise(xi1, xi2, 1 / fs_fine, m);
        fs_c = fs_fine / m;
        nz_c = make_noise(xi1_c, xi2_c, 0, fs_c, sigma, 1:n_sys);
        t_c = (0:1/fs_c:T_conv)';
        for si = 1:3
            names = {'euler', 'heun', 'sra1'};
            [~, Yc] = sde_fixed_step(f, t_c, y0, [], nz_c, names{si});
            err(ri, si, p) = norm(Yc(end, :).' - y_ref_end);
        end
    end
end

% Strong error is an L2 average over paths.
rms_err = sqrt(mean(err.^2, 3));
h_used = refinements(:) / fs_fine;
names = {'euler', 'heun', 'sra1'};
fprintf('    h          euler        heun         sra1\n');
for ri = 1:numel(refinements)
    fprintf('    %.2e   %.3e   %.3e   %.3e\n', h_used(ri), rms_err(ri, :));
end
slopes = zeros(1, 3);
for si = 1:3
    p = polyfit(log(h_used), log(rms_err(:, si)), 1);
    slopes(si) = p(1);
    fprintf('  %-6s measured strong order %.2f\n', names{si}, slopes(si));
end

% Euler and Heun are order 1.0 for additive noise and are checked two-sided:
% both bound the h^1.5 residual of the I(1,0) term the same way, and Heun's
% advantage is in the constant, not the rate.
all_passed = check('euler converges at strong order ~1.0', ...
    abs(slopes(1) - 1.0) < 0.2) && all_passed;
all_passed = check('heun converges at strong order ~1.0', ...
    abs(slopes(2) - 1.0) < 0.2) && all_passed;

% SRA1 is checked as a LOWER bound rather than against 1.5 exactly. Roessler
% guarantees strong order 1.5; that is a floor on the rate, not a prediction of
% it, and with a time-INDEPENDENT diffusion (our case) further terms cancel, so
% the observed rate lands between 1.5 and the scheme's deterministic order of
% 2 -- it measures ~1.7-1.9 here. Pinning it to 1.5 +/- a tolerance would be
% fitting a number the theory does not promise.
%
% What the test has to catch is a mistyped tableau coefficient, and that
% failure mode is unmistakable. Measured on this same system:
%     correct                          -> 1.83
%     dZ term dropped from I(1,0)      -> 1.03
%     B0(2,1) typed as 1 instead of 3/2 -> 0.90
% Both plausible typos collapse to strong order 1.0, so a 1.4 floor separates
% them decisively while leaving the true rate room to be what it is.
all_passed = check('sra1 converges at strong order >= 1.4 (guaranteed 1.5)', ...
    slopes(3) >= 1.4) && all_passed;
all_passed = check('sra1 converges strictly faster than the order-1.0 schemes', ...
    slopes(3) > max(slopes(1:2)) + 0.4) && all_passed;

% The half order really is a gain, not just a slope: at the coarsest step SRA1
% must beat Heun, which must beat Euler.
all_passed = check('SRA1 is more accurate than Heun at the same step', ...
    rms_err(1, 3) < rms_err(1, 2)) && all_passed;
all_passed = check('Heun is more accurate than Euler at the same step', ...
    rms_err(1, 2) < rms_err(1, 1)) && all_passed;

%% Absolute-time indexing: segments must reproduce the single-pass result
fs_seg = 400; T_seg = 1.0;
nsteps = round(T_seg * fs_seg);
nz = make_noise(randn(n_sys, nsteps), randn(n_sys, nsteps), 0, fs_seg, 0.8, 1:n_sys);
t_all = (0:1/fs_seg:T_seg)';
[~, Y_all] = sde_fixed_step(f, t_all, y0, [], nz, 'sra1');

split = round(numel(t_all) / 2);
t_a = t_all(1:split);
t_b = t_all(split:end);              % shares the boundary sample
[~, Y_a] = sde_fixed_step(f, t_a, y0, [], nz, 'sra1');
[~, Y_b] = sde_fixed_step(f, t_b, Y_a(end, :).', [], nz, 'sra1');
all_passed = check('two back-to-back segments reproduce the single pass', ...
    same(Y_b(end, :), Y_all(end, :))) && all_passed;

% ...and re-running one segment from a perturbed start still consumes the same
% increments, which is exactly what Benettin relies on.
pert = 1e-6 * randn(n_sys, 1);
[~, Y_b1] = sde_fixed_step(f, t_b, Y_a(end, :).' + pert, [], nz, 'sra1');
[~, Y_b2] = sde_fixed_step(f, t_b, Y_a(end, :).' + pert, [], nz, 'sra1');
all_passed = check('a perturbed re-integration is reproducible', ...
    same(Y_b1, Y_b2)) && all_passed;
% With a common noise path the additive noise cancels in the difference, so the
% separation stays of order the perturbation rather than of order sigma*sqrt(T).
sep = norm(Y_b1(end, :) - Y_b(end, :));
fprintf('  (separation after a 1e-6 perturbation over %.2f s: %.2e)\n', ...
    t_b(end) - t_b(1), sep);
all_passed = check('additive noise cancels between commonly-driven trajectories', ...
    sep < 1e-3) && all_passed;

%% Noise only reaches the states named by idx
idx_sub = [2 4];
nz_sub = make_noise(randn(2, nsteps), randn(2, nsteps), 0, fs_seg, 5.0, idx_sub);
[~, Y_sub] = sde_fixed_step(f, t_all, y0, [], nz_sub, 'euler');
[~, Y_none] = sde_fixed_step(f, t_all, y0, [], [], 'euler');
% States 1 and 3 are still coupled to the noisy ones through A, so they move --
% but only after the first step, whereas 2 and 4 move within it.
first_step_delta = abs(Y_sub(2, :) - Y_none(2, :));
all_passed = check('the first step moves only the states named by idx', ...
    all(first_step_delta(idx_sub) > 0) && ...
    all(first_step_delta(setdiff(1:n_sys, idx_sub)) == 0)) && all_passed;

%% Error paths
[threw, err_s] = capture(@() sde_fixed_step(f, [0 1], y0, [], [], 'euler'));
all_passed = check('a 2-point span errors', ...
    threw && contains(err_s.identifier, 'SpanTooShort')) && all_passed;

[threw, err_s] = capture(@() sde_fixed_step(f, t_all, y0, [], [], 'rk4'));
all_passed = check('an unknown scheme errors and lists the valid ones', ...
    threw && contains(err_s.identifier, 'UnknownScheme') && ...
    contains(err_s.message, 'sra1')) && all_passed;

% A span reaching past the pre-generated noise must error, not wrap or pad.
t_long = (0:1/fs_seg:(T_seg + 0.5))';
[threw, err_s] = capture(@() sde_fixed_step(f, t_long, y0, [], nz, 'euler'));
all_passed = check('a span past the end of the noise errors', ...
    threw && contains(err_s.identifier, 'NoiseOutOfRange')) && all_passed;

% A span starting before the noise begins likewise.
nz_late = nz; nz_late.t0 = 0.5;
[threw, err_s] = capture(@() sde_fixed_step(f, t_all, y0, [], nz_late, 'euler'));
all_passed = check('a span starting before the noise errors', ...
    threw && contains(err_s.identifier, 'NoiseOutOfRange')) && all_passed;

% A grid that does not match noise.fs would silently mis-pair increments.
t_wrong = (0:1/(fs_seg * 2):T_seg)';
[threw, err_s] = capture(@() sde_fixed_step(f, t_wrong, y0, [], nz, 'euler'));
all_passed = check('a grid that does not match noise.fs errors', ...
    threw && contains(err_s.identifier, 'GridMismatch')) && all_passed;

nz_bad = rmfield(nz, 'xi2');
[threw, err_s] = capture(@() sde_fixed_step(f, t_all, y0, [], nz_bad, 'euler'));
all_passed = check('a noise struct missing a field errors and names it', ...
    threw && contains(err_s.identifier, 'BadNoiseStruct') && ...
    contains(err_s.message, 'xi2')) && all_passed;

nz_bad = nz; nz_bad.idx = 1:2;
[threw, err_s] = capture(@() sde_fixed_step(f, t_all, y0, [], nz_bad, 'euler'));
all_passed = check('a row/idx count mismatch errors', ...
    threw && contains(err_s.identifier, 'BadNoiseStruct')) && all_passed;

%% Result
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');

%% Helpers
function nz = make_noise(xi1, xi2, t0, fs, sigma, idx)
nz = struct('xi1', xi1, 'xi2', xi2, 't0', t0, 'fs', fs, ...
    'sigma', sigma, 'idx', idx);
end

function [xi1_c, xi2_c] = coarsen_noise(xi1, xi2, h, m)
% Rebuild the SAME Brownian path on a step m times coarser.
%
% Increments simply add: dW_c = sum of the m fine increments. The second
% integral does not -- over two consecutive fine steps,
%   I_c = I_a + h*dW_a + I_b
% because the second sub-interval's area is measured from a base that has
% already risen by dW_a. Applying that pairwise, log2(m) times, is exact.
% The pair is then converted back to unit-variance normals by inverting the
% Kloeden-Platen identity I = (H/2)*(dW + dZ/sqrt(3)).
    assert(mod(log2(m), 1) == 0, 'coarsen_noise expects m to be a power of 2');
    dW = sqrt(h) * xi1;
    I10 = (h / 2) * (dW + sqrt(h) * xi2 / sqrt(3));
    hc = h;
    while size(dW, 2) > size(xi1, 2) / m
        a = 1:2:size(dW, 2);
        b = 2:2:size(dW, 2);
        I10 = I10(:, a) + hc * dW(:, a) + I10(:, b);
        dW = dW(:, a) + dW(:, b);
        hc = 2 * hc;
    end
    H = hc;
    xi1_c = dW / sqrt(H);
    xi2_c = sqrt(3) * (2 * I10 / H - dW) / sqrt(H);
end

function Y = ref_rk2(f, t, y0, node)
% Two-stage RK with c2 = node, b2 = 1/(2*node), b1 = 1 - b2.
b2 = 1 / (2 * node); b1 = 1 - b2;
y = y0(:); Y = zeros(numel(t), numel(y)); Y(1, :) = y.';
for k = 1:numel(t)-1
    h = t(k+1) - t(k);
    k1 = f(t(k), y);
    k2 = f(t(k) + node * h, y + node * h * k1);
    y = y + h * (b1 * k1 + b2 * k2);
    Y(k+1, :) = y.';
end
end

function Y = ref_trapezoid(f, t, y0)
y = y0(:); Y = zeros(numel(t), numel(y)); Y(1, :) = y.';
for k = 1:numel(t)-1
    h = t(k+1) - t(k);
    k1 = f(t(k), y);
    k2 = f(t(k) + h, y + h * k1);
    y = y + (h / 2) * (k1 + k2);
    Y(k+1, :) = y.';
end
end

function Y = ref_forward_euler(f, t, y0)
y = y0(:); Y = zeros(numel(t), numel(y)); Y(1, :) = y.';
for k = 1:numel(t)-1
    h = t(k+1) - t(k);
    y = y + h * f(t(k), y);
    Y(k+1, :) = y.';
end
end

function tf = same(a, b)
tf = isequal(size(a), size(b)) && max(abs(a(:) - b(:))) < 1e-14;
end

function [threw, err] = capture(fcn)
threw = false; err = [];
try
    fcn();
catch err
    threw = true;
end
end

function passed = check(name, condition)
if condition
    fprintf('  %s: PASS\n', name);
    passed = true;
else
    fprintf('  %s: FAIL\n', name);
    passed = false;
end
end
