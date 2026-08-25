% test_SRNNCellTypePairs_noise.m
% The SRNNCellTypePairs counterpart of test_SRNN2_noise. The two classes are
% duck-typed siblings that share no implementation, so the noise machinery is
% two parallel copies and has to be verified twice.
%
% This file checks the same core claims -- input-referred parameterisation, an
% undisturbed RNG stream, the analytic OU variance, noise reaching only x,
% seed control, the integrator requirement, and Benettin under noise -- on a
% THREE-type model, so nothing can silently depend on the E/I special case.
% The one structural difference from SRNNModel2 is where the dendritic block
% lives: this class names it directly as state_layout.x rather than deriving it
% as the trailing n entries.
%
% See also: SRNNCellTypePairs, test_SRNN2_noise, sde_fixed_step

fprintf('=== Testing SRNNCellTypePairs additive noise ===\n\n');
all_passed = true;

%% The input-referred parameterisation
m = pairs('sigma_u_noise', 0.05, 'tau_d', 0.1, 'ode_solver', 'sra1');
all_passed = check('sigma_x_raw = sigma_u_noise / tau_d', ...
    abs(m.sigma_x_raw - 0.5) < 1e-12) && all_passed;
all_passed = check('x_noise_std = sigma_u_noise / sqrt(2*tau_d)', ...
    abs(m.x_noise_std - 0.05 / sqrt(0.2)) < 1e-12) && all_passed;
all_passed = check('the default is noise-free', ...
    pairs().sigma_u_noise == 0) && all_passed;

[threw, err] = capture(@() pairs('sigma_x_raw', 1));
all_passed = check('assigning sigma_x_raw names sigma_u_noise', ...
    threw && strcmp(err.identifier, 'SRNNCellTypePairs:RenamedProperty') && ...
    contains(err.message, 'sigma_u_noise')) && all_passed;

%% Turning noise on leaves W, the stimulus and x0 bit-identical
plain = pairs('sigma_u_noise', 0, 'ode_solver', 'rk4');
noisy = pairs('sigma_u_noise', 0.4, 'ode_solver', 'sra1');
plain.build(); noisy.build();
all_passed = check('W is bit-identical with and without noise', ...
    max(abs(full(plain.W) - full(noisy.W)), [], 'all') == 0) && all_passed;
all_passed = check('S0 is bit-identical with and without noise', ...
    max(abs(plain.S0 - noisy.S0)) == 0) && all_passed;
all_passed = check('the external input is bit-identical too', ...
    max(abs(plain.u_ex - noisy.u_ex), [], 'all') == 0) && all_passed;

%% The noise tensor targets state_layout.x and is cleared after the run
nz = pairs('sigma_u_noise', 0.4, 'ode_solver', 'sra1');
nz.build();
nz.build_noise_for_test();
all_passed = check('the tensor targets state_layout.x', ...
    isequal(nz.noise_increments.idx, nz.cached_params.state_layout.x)) && all_passed;
all_passed = check('...with one row per neuron', ...
    size(nz.noise_increments.xi1, 1) == nz.n) && all_passed;
evalc('nz.run();');
all_passed = check('the noise tensor is cleared once the run is done', ...
    isempty(nz.noise_increments)) && all_passed;

%% sigma_u_noise means what it says: an OU process with W = 0 and no stimulus
tau_d = 0.1; sigma_u = 0.2;
ou = SRNNCellTypePairs('n', 180, 'indegree', 20, 'n_cellTypes', 3, ...
    'cell_type_names', {'A', 'B', 'C'}, 'f', [0.4 0.35 0.25], ...
    'mu_tilde_relative', [3 -4 -2], 'sigma_tilde_relative', [1 1 1], ...
    'tau_a', {[], [], []}, 'level_of_chaos', 0, 'u_ex_scale', 0, ...
    'tau_d', tau_d, 'sigma_u_noise', sigma_u, 'ode_solver', 'sra1', ...
    'T_range', [0 40], 'fs', 400, 'lya_method', 'none', ...
    'store_full_state', true, 'x0_std', 0);
ou.build();
evalc('ou.run();');
x_samples = ou.S_out(ou.t_out > 10 * tau_d, ou.cached_params.state_layout.x);
var_measured = var(x_samples(:));
var_expected = ou.sigma_x_raw^2 * tau_d / 2;
rel = abs(var_measured - var_expected) / var_expected;
fprintf('  (OU stationary var: measured %.5f, analytic %.5f, %.1f%% off)\n', ...
    var_measured, var_expected, 100 * rel);
all_passed = check('the stationary variance of x matches sigma_x_raw^2*tau_d/2', ...
    rel < 0.05) && all_passed;

%% Noise drives only x
% n_a is nonzero here so there are adaptation states to leave alone.
a = pairs('sigma_u_noise', 0, 'ode_solver', 'euler', 'tau_a', {srnn_sfa_timescales(2), srnn_sfa_timescales(1), srnn_sfa_timescales(1)});
b = pairs('sigma_u_noise', 0.5, 'ode_solver', 'euler', 'tau_a', {srnn_sfa_timescales(2), srnn_sfa_timescales(1), srnn_sfa_timescales(1)});
a.build(); evalc('a.run();');
b.build(); evalc('b.run();');
xb = a.cached_params.state_layout.x;
other = setdiff(1:a.N_sys_eqs, xb);
step1 = abs(b.S_out(2, :) - a.S_out(2, :));
all_passed = check('the state vector really has non-x states to check', ...
    ~isempty(other)) && all_passed;
all_passed = check('the first step moves the dendritic block', ...
    all(step1(xb) > 0)) && all_passed;
all_passed = check('...and leaves every other state exactly alone', ...
    ~isempty(other) && max(step1(other)) == 0) && all_passed;

%% Reproducibility and seed control
r1 = run_traj(pairs('sigma_u_noise', 0.4, 'ode_solver', 'sra1'));
r2 = run_traj(pairs('sigma_u_noise', 0.4, 'ode_solver', 'sra1'));
all_passed = check('the same configuration reproduces the same trajectory', ...
    isequal(r1, r2)) && all_passed;
r3 = run_traj(pairs('sigma_u_noise', 0.4, 'ode_solver', 'sra1', 'noise_seed', 99));
all_passed = check('a different noise_seed gives a different realisation', ...
    ~isequal(r1, r3)) && all_passed;

n1 = pairs('sigma_u_noise', 0.4, 'ode_solver', 'sra1', 'noise_seed', 7, 'rng_seeds', [11 12]);
n2 = pairs('sigma_u_noise', 0.4, 'ode_solver', 'sra1', 'noise_seed', 7, 'rng_seeds', [21 22]);
n1.build(); n2.build();
n1.build_noise_for_test(); n2.build_noise_for_test();
all_passed = check('an explicit noise_seed pins the draw across networks', ...
    isequal(n1.noise_increments.xi1, n2.noise_increments.xi1) && ...
    ~isequal(full(n1.W), full(n2.W))) && all_passed;

%% Noise requires a stochastic integrator
for name = {'ode45', 'ode15s', 'rk4'}
    [threw, err] = capture(@() build_it('sigma_u_noise', 0.1, 'ode_solver', name{1}));
    all_passed = check(sprintf('sigma > 0 with ''%s'' errors', name{1}), ...
        threw && contains(err.message, 'sra1')) && all_passed;
end

%% Benettin and QR both survive the noise
lya_args = {'T_range', [0 12], 'fs', 200, 'lya_method', 'benettin', ...
    'lya_T_interval', [6 12], 'lya_warmup', 4, 'ode_solver', 'sra1'};
b0 = pairs('sigma_u_noise', 0, lya_args{:}); b0.build(); evalc('b0.run();');
b1 = pairs('sigma_u_noise', 0.3, lya_args{:}); b1.build(); evalc('b1.run();');
b2 = pairs('sigma_u_noise', 0.3, lya_args{:}); b2.build(); evalc('b2.run();');
fprintf('  (LLE at sigma = 0: %+.4f; at sigma = 0.3: %+.4f)\n', ...
    b0.lya_results.LLE, b1.lya_results.LLE);
all_passed = check('Benettin returns a finite LLE under noise', ...
    isfinite(b1.lya_results.LLE)) && all_passed;
all_passed = check('...reproducibly, so the perturbed run saw the same path', ...
    b1.lya_results.LLE == b2.lya_results.LLE) && all_passed;

q = pairs('sigma_u_noise', 0.2, 'n', 9, 'indegree', 4, 'tau_a', {[], [], []}, ...
    'T_range', [0 4], 'fs', 200, 'ode_solver', 'sra1', 'lya_method', 'qr', ...
    'lya_T_interval', [2 4], 'lya_warmup', 2);
q.build(); evalc('q.run();');
all_passed = check('QR returns a finite spectrum under noise', ...
    numel(q.lya_results.LE_spectrum) == q.N_sys_eqs && ...
    all(isfinite(q.lya_results.LE_spectrum))) && all_passed;

%% Result
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');

%% Helpers
function p = pairs(varargin)
% Three cell types deliberately, so nothing can depend on the E/I special case.
defaults = {'n', 24, 'indegree', 12, 'n_cellTypes', 3, ...
    'cell_type_names', {'A', 'B', 'C'}, 'f', [0.4 0.35 0.25], ...
    'mu_tilde_relative', [3 -4 -2], 'sigma_tilde_relative', [1 1 1], ...
    'tau_a', {srnn_sfa_timescales(2), srnn_sfa_timescales(1), []}, 'T_range', [0 3], 'fs', 200, ...
    'lya_method', 'none', 'store_full_state', true};
p = SRNNPairsTestAccess(defaults{:}, varargin{:});
end

function p = build_it(varargin)
p = pairs(varargin{:});
p.build();
end

function S = run_traj(p)
p.build();
evalc('p.run();');
S = p.S_out;
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
