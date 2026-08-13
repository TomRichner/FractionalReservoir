% test_SRNN2_noise.m
% Verify the additive Wiener noise on SRNNModel2's dendritic state:
%
%   dx_i = (-x_i + sum_j w_ij b_j r_j + u_i)/tau_d dt + sigma_u_noise/tau_d dW_i
%
% The claims checked here:
%   1. sigma_u_noise is INPUT-REFERRED -- the Dependent sigma_x_raw and
%      x_noise_std derive from it, and assigning either names the setter.
%   2. Turning noise on does not disturb W, the stimulus or x0: the draw takes
%      its own RNG substream and restores the stream, so sigma = 0 stays
%      bit-identical to the deterministic model.
%   3. sigma_u_noise means what the docs say -- with W = 0 and no stimulus the
%      dendritic state is an Ornstein-Uhlenbeck process whose stationary
%      variance matches sigma_x_raw^2 * tau_d / 2.
%   4. Noise reaches ONLY x. Nothing else in the state vector is driven
%      directly, which is what keeps the diffusion constant.
%   5. The realisation is reproducible and seed-controlled.
%   6. Asking for noise with a deterministic integrator is an error.
%   7. Benettin still works, and is reproducible, because the perturbed
%      re-integration sees the same noise path.
%
% See also: SRNNModel2, sde_fixed_step, test_sde_integrators

fprintf('=== Testing SRNNModel2 additive noise ===\n\n');
all_passed = true;

%% The input-referred parameterisation
m = SRNNModel2('sigma_u_noise', 0.05, 'tau_d', 0.1);
all_passed = check('sigma_x_raw = sigma_u_noise / tau_d', ...
    abs(m.sigma_x_raw - 0.5) < 1e-12) && all_passed;
all_passed = check('x_noise_std = sigma_u_noise / sqrt(2*tau_d)', ...
    abs(m.x_noise_std - 0.05 / sqrt(0.2)) < 1e-12) && all_passed;
all_passed = check('the default is noise-free', ...
    SRNNModel2().sigma_u_noise == 0) && all_passed;

[threw, err] = capture(@() SRNNModel2('sigma_x_raw', 1));
all_passed = check('assigning sigma_x_raw names sigma_u_noise', ...
    threw && strcmp(err.identifier, 'SRNNModel:RenamedProperty') && ...
    contains(err.message, 'sigma_u_noise')) && all_passed;
[threw, err] = capture(@() SRNNModel2('x_noise_std', 1));
all_passed = check('assigning x_noise_std names sigma_u_noise', ...
    threw && strcmp(err.identifier, 'SRNNModel:RenamedProperty') && ...
    contains(err.message, 'sigma_u_noise')) && all_passed;

%% Turning noise on leaves W, the stimulus and x0 bit-identical
plain = tiny('sigma_u_noise', 0, 'ode_solver', 'rk4');
noisy = tiny('sigma_u_noise', 0.4, 'ode_solver', 'sra1');
plain.build(); noisy.build();
all_passed = check('W is bit-identical with and without noise', ...
    max(abs(full(plain.W) - full(noisy.W)), [], 'all') == 0) && all_passed;
all_passed = check('S0 is bit-identical with and without noise', ...
    max(abs(plain.S0 - noisy.S0)) == 0) && all_passed;
all_passed = check('the external input is bit-identical too', ...
    max(abs(plain.u_ex - noisy.u_ex), [], 'all') == 0) && all_passed;

%% sigma = 0 allocates no noise at all
m0 = tiny('sigma_u_noise', 0, 'ode_solver', 'sra1');
m0.build(); evalc('m0.run();');
all_passed = check('sigma = 0 leaves the noise tensor empty', ...
    isempty(m0.noise_increments)) && all_passed;

%% The noise tensor is dropped after the run
mrun = tiny('sigma_u_noise', 0.4, 'ode_solver', 'sra1');
mrun.build(); evalc('mrun.run();');
all_passed = check('the noise tensor is cleared once the run is done', ...
    isempty(mrun.noise_increments)) && all_passed;

%% sigma_u_noise means what it says: an OU process with W = 0 and no stimulus
% level_of_chaos = 0 kills W, u_ex_scale = 0 kills the stimulus, so
% dx = -x/tau_d dt + sigma_x_raw dW exactly, with stationary variance
% sigma_x_raw^2 * tau_d / 2.
tau_d = 0.1; sigma_u = 0.2;
ou = SRNNModel2('n', 200, 'indegree', 20, 'level_of_chaos', 0, ...
    'u_ex_scale', 0, 'tau_d', tau_d, 'sigma_u_noise', sigma_u, ...
    'ode_solver', 'sra1', 'T_range', [0 40], 'fs', 400, ...
    'lya_method', 'none', 'store_full_state', true, 'x0_std', 0);
ou.build();
evalc('ou.run();');
x_block = (ou.N_sys_eqs - ou.n + 1):ou.N_sys_eqs;
settled = ou.t_out > 10 * tau_d;             % discard the OU transient
x_samples = ou.S_out(settled, x_block);
var_measured = var(x_samples(:));
var_expected = ou.sigma_x_raw^2 * tau_d / 2;
rel = abs(var_measured - var_expected) / var_expected;
fprintf('  (OU stationary var: measured %.5f, analytic %.5f, %.1f%% off)\n', ...
    var_measured, var_expected, 100 * rel);
all_passed = check('the stationary variance of x matches sigma_x_raw^2*tau_d/2', ...
    rel < 0.05) && all_passed;
all_passed = check('...and x_noise_std is its square root', ...
    abs(sqrt(var_expected) - ou.x_noise_std) < 1e-12) && all_passed;

%% Noise drives only x
% Over the very first step nothing but the dendritic block can have moved: the
% adaptation and depression states are driven by r, which depends on the state
% at the START of the step and is therefore identical in both runs.
%
% This needs a model that HAS other states, so SFA and STD are switched on
% explicitly -- the class defaults are n_a_E = n_b_E = 0, which would leave the
% state vector as nothing but x and make the check vacuous.
adapt = {'n_a_E', 3, 'n_a_I', 2, 'n_b_E', 1, 'n_b_I', 1};
a = tiny(adapt{:}, 'sigma_u_noise', 0, 'ode_solver', 'euler');
b = tiny(adapt{:}, 'sigma_u_noise', 0.5, 'ode_solver', 'euler');
a.build(); evalc('a.run();');
b.build(); evalc('b.run();');
xb = (a.N_sys_eqs - a.n + 1):a.N_sys_eqs;
other = setdiff(1:a.N_sys_eqs, xb);
step1 = abs(b.S_out(2, :) - a.S_out(2, :));
all_passed = check('the state vector really has non-x states to check', ...
    ~isempty(other)) && all_passed;
all_passed = check('the first step moves the dendritic block', ...
    all(step1(xb) > 0)) && all_passed;
all_passed = check('...and leaves every other state exactly alone', ...
    ~isempty(other) && max(step1(other)) == 0) && all_passed;

%% Reproducibility and seed control
r1 = run_traj(tiny('sigma_u_noise', 0.4, 'ode_solver', 'sra1'));
r2 = run_traj(tiny('sigma_u_noise', 0.4, 'ode_solver', 'sra1'));
all_passed = check('the same configuration reproduces the same trajectory', ...
    isequal(r1, r2)) && all_passed;

r3 = run_traj(tiny('sigma_u_noise', 0.4, 'ode_solver', 'sra1', 'noise_seed', 99));
all_passed = check('a different noise_seed gives a different realisation', ...
    ~isequal(r1, r3)) && all_passed;

% An explicit noise_seed pins the realisation even across different networks;
% the trajectories differ (different W) but the increments do not.
n1 = tiny('sigma_u_noise', 0.4, 'ode_solver', 'sra1', 'noise_seed', 7, 'rng_seeds', [11 12]);
n2 = tiny('sigma_u_noise', 0.4, 'ode_solver', 'sra1', 'noise_seed', 7, 'rng_seeds', [21 22]);
n1.build(); n2.build();
n1.build_noise_for_test(); n2.build_noise_for_test();
all_passed = check('an explicit noise_seed pins the draw across networks', ...
    isequal(n1.noise_increments.xi1, n2.noise_increments.xi1) && ...
    ~isequal(full(n1.W), full(n2.W))) && all_passed;

% Left empty, the seed follows rng_seeds, so a reps sweep varies the noise too.
d1 = tiny('sigma_u_noise', 0.4, 'ode_solver', 'sra1', 'rng_seeds', [11 12]);
d2 = tiny('sigma_u_noise', 0.4, 'ode_solver', 'sra1', 'rng_seeds', [21 22]);
d1.build(); d2.build();
d1.build_noise_for_test(); d2.build_noise_for_test();
all_passed = check('an empty noise_seed follows rng_seeds', ...
    ~isequal(d1.noise_increments.xi1, d2.noise_increments.xi1)) && all_passed;

%% Noise requires a stochastic integrator
for name = {'ode45', 'ode15s', 'rk4'}
    [threw, err] = capture(@() build_it('sigma_u_noise', 0.1, 'ode_solver', name{1}));
    all_passed = check(sprintf('sigma > 0 with ''%s'' errors', name{1}), ...
        threw && contains(err.message, 'sra1')) && all_passed;
end
[threw, err] = capture(@() build_it('sigma_u_noise', -1, 'ode_solver', 'sra1'));
all_passed = check('a negative sigma_u_noise errors', ...
    threw && contains(err.message, 'non-negative')) && all_passed;
% ...but a stochastic integrator at sigma = 0 is fine (the convergence tests
% and any deterministic-vs-stochastic comparison rely on it).
ok = build_it('sigma_u_noise', 0, 'ode_solver', 'heun');
all_passed = check('a stochastic integrator at sigma = 0 is allowed', ...
    ok.is_built) && all_passed;

%% Benettin survives the noise, and is reproducible
lya_args = {'T_range', [0 12], 'fs', 200, 'lya_method', 'benettin', ...
    'lya_T_interval', [6 12], 'lya_warmup', 4, 'ode_solver', 'sra1'};
b0 = tiny('sigma_u_noise', 0, lya_args{:});
b0.build(); evalc('b0.run();');
b1 = tiny('sigma_u_noise', 0.3, lya_args{:});
b1.build(); evalc('b1.run();');
b2 = tiny('sigma_u_noise', 0.3, lya_args{:});
b2.build(); evalc('b2.run();');
fprintf('  (LLE at sigma = 0: %+.4f; at sigma = 0.3: %+.4f)\n', ...
    b0.lya_results.LLE, b1.lya_results.LLE);
all_passed = check('Benettin returns a finite LLE under noise', ...
    isfinite(b1.lya_results.LLE)) && all_passed;
all_passed = check('...reproducibly, so the perturbed run saw the same path', ...
    b1.lya_results.LLE == b2.lya_results.LLE) && all_passed;
all_passed = check('...and the noise actually moved it', ...
    abs(b1.lya_results.LLE - b0.lya_results.LLE) > 1e-6) && all_passed;

%% QR survives the noise too
% Additive noise never enters the variational equation, so the only thing that
% changes is the trajectory the Jacobian is evaluated along.
q = SRNNModel2('n', 8, 'indegree', 4, 'n_a_E', 0, 'n_a_I', 0, ...
    'n_b_E', 0, 'n_b_I', 0, 'T_range', [0 4], 'fs', 200, ...
    'sigma_u_noise', 0.2, 'ode_solver', 'sra1', 'lya_method', 'qr', ...
    'lya_T_interval', [2 4], 'lya_warmup', 2, 'store_full_state', true);
q.build(); evalc('q.run();');
all_passed = check('QR returns a finite spectrum under noise', ...
    numel(q.lya_results.LE_spectrum) == q.N_sys_eqs && ...
    all(isfinite(q.lya_results.LE_spectrum))) && all_passed;

%% Mechanism: noise saturates the nonlinearity and drops the effective gain
% This is the route by which additive noise lowers the Lyapunov exponent -- it
% pushes x_eff away from the setpoint, where phi' is smaller. Reported as well
% as asserted, because it is the number that turns "the LLE fell" into a claim
% with a mechanism behind it.
% Run without adaptation so x_eff is exactly x and the measurement needs no
% unpacking; the mechanism does not depend on SFA being present.
gains = zeros(1, 2); sigmas = [0 0.6];
for k = 1:2
    g = tiny('sigma_u_noise', sigmas(k), 'ode_solver', 'sra1', ...
        'n_a_E', 0, 'n_a_I', 0, 'n_b_E', 0, 'n_b_I', 0, ...
        'T_range', [0 12], 'fs', 200, 'lya_method', 'none');
    g.build(); evalc('g.run();');
    x_only = g.S_out(g.t_out > 2, :);       % no adaptation => x_eff == x
    gains(k) = mean(g.activation_function_derivative(x_only(:)));
end
fprintf('  (mean phi'' at sigma = 0: %.4f; at sigma = 0.6: %.4f)\n', gains);
all_passed = check('noise lowers the mean effective gain', ...
    gains(2) < gains(1)) && all_passed;

%% Result
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');

%% Helpers
function m = tiny(varargin)
defaults = {'n', 24, 'indegree', 12, 'T_range', [0 3], 'fs', 200, ...
    'lya_method', 'none', 'store_full_state', true};
m = SRNN2TestAccess(defaults{:}, varargin{:});
end

function m = build_it(varargin)
m = tiny(varargin{:});
m.build();
end

function S = run_traj(m)
m.build();
evalc('m.run();');
S = m.S_out;
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
