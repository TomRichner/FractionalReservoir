% test_numerics_probe.m - SRNNNumericsProbe opens exactly the doors
% run_numerics_verification needs, and nothing else changes.
%
% The probe is a subclass of SRNNCellTypePairs that (1) lets a caller inject
% the Brownian path a run consumes, (2) keeps a handle to the path run()
% consumed after run() clears it, and (3) re-integrates a segment from a
% supplied state at an absolute time with the model's own integrator. The
% claims checked here:
%
%   1. With nothing injected the probe draws exactly the parent's path, and a
%      second probe with the same seeds draws the same one.
%   2. An injected path is what run() consumes: injecting the negated path
%      changes the trajectory, injecting the identical path reproduces it
%      bit-for-bit.
%   3. integrate_segment from a mid-run state reproduces the stored rows
%      bit-for-bit, deterministic AND stochastic. The stochastic case is the
%      one that matters: it proves the absolute-time noise indexing pairs a
%      re-integrated segment with the increments the original run used.
%   4. coarsen_noise (src/model/integrators) preserves the path: m = 1 is the
%      identity, and at m = 2 the coarse increments are the sums of the fine
%      ones. A coarse-step reshoot from a fine-reference state then lands
%      close to the fine reference -- the strong error, small but not zero.
%   5. set_noise rejects a path on the wrong grid.
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths has run.
%
% See also: SRNNNumericsProbe, run_numerics_verification, coarsen_noise,
%           test_sde_integrators

fprintf('=== Testing SRNNNumericsProbe ===\n\n');
all_passed = true;

common = {'n', 24, 'indegree', 12, 'n_cellTypes', 2, ...
    'cell_type_names', {'E', 'I'}, 'f', [0.5 0.5], ...
    'mu_tilde_relative', [3 -4], 'sigma_tilde_relative', [1 1], ...
    'tau_a', {log_ladder(0.25, 10, 3), []}, ...
    'T_range', [0 2], 'lya_method', 'none', 'store_full_state', true};

%% 1. Fall-through: no injection means the parent's draw
p1 = SRNNNumericsProbe(common{:}, 'fs', 200, 'ode_solver', 'sra1', 'sigma_u_noise', 0.05);
p1.build(); evalc('p1.run();');
nz1 = p1.consumed_noise;
all_passed = check('run() clears noise_increments as the parent does', ...
    isempty(p1.noise_increments)) && all_passed;
all_passed = check('consumed_noise keeps the path run() used', ...
    isstruct(nz1) && nz1.fs == 200 && size(nz1.xi1, 2) == numel(p1.t_out) - 1) && all_passed;

p1b = SRNNNumericsProbe(common{:}, 'fs', 200, 'ode_solver', 'sra1', 'sigma_u_noise', 0.05);
p1b.build(); evalc('p1b.run();');
all_passed = check('same seeds draw the same path (parent behaviour intact)', ...
    isequal(p1b.consumed_noise.xi1, nz1.xi1) && isequal(p1b.S_out, p1.S_out)) && all_passed;

%% 2. Injection is what run() consumes
p2 = SRNNNumericsProbe(common{:}, 'fs', 200, 'ode_solver', 'sra1', 'sigma_u_noise', 0.05);
nz_neg = nz1; nz_neg.xi1 = -nz1.xi1; nz_neg.xi2 = -nz1.xi2;
p2.set_noise(nz_neg);
p2.build(); evalc('p2.run();');
all_passed = check('an injected (negated) path changes the trajectory', ...
    isequal(p2.consumed_noise.xi1, -nz1.xi1) && ~isequal(p2.S_out, p1.S_out)) && all_passed;

p3 = SRNNNumericsProbe(common{:}, 'fs', 200, 'ode_solver', 'sra1', 'sigma_u_noise', 0.05);
p3.set_noise(nz1);
p3.build(); evalc('p3.run();');
all_passed = check('injecting the identical path reproduces the run bit-for-bit', ...
    isequal(p3.S_out, p1.S_out)) && all_passed;

p3.clear_noise(); evalc('p3.run();');
all_passed = check('clear_noise returns to the seeded draw', ...
    isequal(p3.S_out, p1.S_out) && isequal(p3.consumed_noise.xi1, nz1.xi1)) && all_passed;

%% 3. integrate_segment reproduces stored rows from a mid-run state
k = 101; L = 10;
p4 = SRNNNumericsProbe(common{:}, 'fs', 200, 'ode_solver', 'sra1', 'sigma_u_noise', 0);
p4.build(); evalc('p4.run();');
[~, seg] = p4.integrate_segment(p4.t_out(k:k+L), p4.S_out(k, :)');
all_passed = check('deterministic segment reproduces S_out rows bit-for-bit', ...
    isequal(seg, p4.S_out(k:k+L, :))) && all_passed;

p1.arm_noise();      % run() cleared the tensor; put the same path back
[~, seg] = p1.integrate_segment(p1.t_out(k:k+L), p1.S_out(k, :)');
all_passed = check('stochastic segment reproduces S_out rows bit-for-bit (absolute-time noise)', ...
    isequal(seg, p1.S_out(k:k+L, :))) && all_passed;
p1.disarm_noise();
all_passed = check('disarm_noise drops the tensor', isempty(p1.noise_increments)) && all_passed;

%% 4. coarsen_noise preserves the path
h = 1 / 200;
[c1, c2] = coarsen_noise(nz1.xi1, nz1.xi2, h, 1);
all_passed = check('coarsen_noise at m = 1 is the identity', ...
    max(abs(c1 - nz1.xi1), [], 'all') < 1e-12 && max(abs(c2 - nz1.xi2), [], 'all') < 1e-12) && all_passed;

[c1, ~] = coarsen_noise(nz1.xi1, nz1.xi2, h, 2);
dW_fine_sum = sqrt(h) * (nz1.xi1(:, 1:2:end) + nz1.xi1(:, 2:2:end));
all_passed = check('coarsen_noise at m = 2 sums the increments exactly', ...
    max(abs(sqrt(2 * h) * c1 - dW_fine_sum), [], 'all') < 1e-12) && all_passed;

% A coarse-step (fs 100) reshoot against the fine (fs 200) reference, on the
% same path. Two coarse steps (the shortest span sde_fixed_step accepts) from
% the reference state should land near the reference four fine steps later:
% not equal (that is the discretisation error) but far closer than the state moved.
[c1, c2] = coarsen_noise(nz1.xi1, nz1.xi2, h, 2);
p5 = SRNNNumericsProbe(common{:}, 'fs', 100, 'ode_solver', 'sra1', 'sigma_u_noise', 0.05);
p5.set_noise(struct('xi1', c1, 'xi2', c2, 't0', 0, 'fs', 100, 'sigma', 0, 'idx', []));
p5.build(); p5.arm_noise();
kk = 201;                                   % fine row 201 = t = 1.0 s, on the coarse grid
t_seg = p1.t_out(kk) + [0, 1, 2] / 100;
[~, seg] = p5.integrate_segment(t_seg, p1.S_out(kk, :)');
err   = norm(seg(end, :) - p1.S_out(kk + 4, :));
moved = norm(p1.S_out(kk + 4, :) - p1.S_out(kk, :));
fprintf('  (coarse reshoot error %.3g vs state movement %.3g over two coarse steps)\n', err, moved);
all_passed = check('a coarse-step reshoot on the coarsened path tracks the fine reference', ...
    err > 0 && err < 0.05 * moved) && all_passed;

%% 5. Grid checks
threw = false;
try
    p5.set_noise(struct('xi1', c1, 'xi2', c2, 't0', 0, 'fs', 200, 'sigma', 0, 'idx', []));
catch ME
    threw = strcmp(ME.identifier, 'SRNNNumericsProbe:NoiseGridMismatch');
end
all_passed = check('set_noise rejects a path at the wrong fs', threw) && all_passed;

threw = false;
try
    p5.set_noise(struct('xi1', c1, 'xi2', c2, 't0', -1, 'fs', 100, 'sigma', 0, 'idx', []));
catch ME
    threw = strcmp(ME.identifier, 'SRNNNumericsProbe:NoiseGridMismatch');
end
all_passed = check('set_noise rejects a path with the wrong t0', threw) && all_passed;

%% Summary
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');

function passed = check(name, condition)
if condition
    fprintf('  %s: PASS\n', name);
    passed = true;
else
    fprintf('  %s: FAIL\n', name);
    passed = false;
end
end
