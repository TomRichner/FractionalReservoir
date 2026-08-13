% test_SRNN2_S_c_heterogeneity.m
% Verify per-neuron nonlinearity setpoints in SRNNModel2: mu_S_c_E/I and
% sigma_S_c_E/I are drawn at build() into the read-only S_c_vec, and the
% activation handles then carry a per-neuron centre.
%
% Three things have to hold and each is checked numerically rather than
% structurally:
%   1. With the new properties at their defaults NOTHING changes -- same W,
%      same S0, same phi, empty S_c_vec.
%   2. The draw has the requested per-population statistics, is reproducible,
%      and does not perturb the network or stimulus RNG.
%   3. phi with a vector centre equals the elementwise scalar-centre result,
%      and the analytic Jacobian still matches finite differences -- which is
%      what keeps the Lyapunov exponents meaningful.
%
% See also: SRNNModel2, test_SRNN2_activation, test_SRNNCellTypePairs_S_c_heterogeneity

fprintf('=== Testing SRNNModel2 per-neuron S_c ===\n\n');
all_passed = true;
xg = linspace(-1, 2, 101)';

%% Defaults are inert: no vector, and the scalar path is untouched
m_base = tiny_model();
m_base.build();
all_passed = check('S_c_vec is empty with the properties at their defaults', ...
    isempty(m_base.S_c_vec)) && all_passed;
all_passed = check('the handle is still the scalar-S_c function', ...
    same(m_base.activation_function(xg), ...
         SRNNModel2.logisticSigmoid(xg, m_base.S_c))) && all_passed;

%% Distinct means, no spread: exactly two values, in the right places
m = tiny_model('mu_S_c_E', 0.40, 'mu_S_c_I', 0.10);
m.build();
v = m.S_c_vec;
all_passed = check('S_c_vec is an n x 1 column', ...
    isequal(size(v), [m.n, 1])) && all_passed;
all_passed = check('mu_S_c_E applies to the E block', ...
    all(v(m.E_indices) == 0.40)) && all_passed;
all_passed = check('mu_S_c_I applies to the I block', ...
    all(v(m.I_indices) == 0.10)) && all_passed;
all_passed = check('sigma = 0 gives no spread at all', ...
    numel(unique(v)) == 2) && all_passed;

%% An empty mu falls back to the scalar S_c
m = tiny_model('S_c', 0.33, 'mu_S_c_I', 0.10);
m.build();
all_passed = check('an empty mu_S_c_E falls back to S_c', ...
    all(m.S_c_vec(m.E_indices) == 0.33) && ...
    all(m.S_c_vec(m.I_indices) == 0.10)) && all_passed;

%% Setting sigma alone adds spread around the shared S_c
m = tiny_model('S_c', 0.4, 'sigma_S_c_E', 0.05, 'sigma_S_c_I', 0.05, ...
    'n', 1000, 'indegree', 100);
m.build();
v = m.S_c_vec;
all_passed = check('sigma alone still centres on S_c', ...
    abs(mean(v) - 0.4) < 5e-3) && all_passed;

%% The requested statistics come out of the draw
m = tiny_model('n', 1000, 'indegree', 100, ...
    'mu_S_c_E', 0.40, 'sigma_S_c_E', 0.05, ...
    'mu_S_c_I', 0.10, 'sigma_S_c_I', 0.20);
m.build();
vE = m.S_c_vec(m.E_indices);
vI = m.S_c_vec(m.I_indices);
% 3 standard errors of the mean / of the sample SD, at n_E = n_I = 500.
all_passed = check('E sample mean matches mu_S_c_E', ...
    abs(mean(vE) - 0.40) < 3 * 0.05 / sqrt(numel(vE))) && all_passed;
all_passed = check('I sample mean matches mu_S_c_I', ...
    abs(mean(vI) - 0.10) < 3 * 0.20 / sqrt(numel(vI))) && all_passed;
all_passed = check('E sample SD matches sigma_S_c_E', ...
    abs(std(vE) - 0.05) < 3 * 0.05 / sqrt(2 * numel(vE))) && all_passed;
all_passed = check('I sample SD matches sigma_S_c_I', ...
    abs(std(vI) - 0.20) < 3 * 0.20 / sqrt(2 * numel(vI))) && all_passed;
all_passed = check('the two populations are drawn independently', ...
    abs(mean(vE) - mean(vI)) > 0.2) && all_passed;

%% Reproducibility and seed control
args = {'mu_S_c_E', 0.4, 'sigma_S_c_E', 0.05};
a = tiny_model(args{:}); a.build();
b = tiny_model(args{:}); b.build();
all_passed = check('the same configuration reproduces the same draw', ...
    isequal(a.S_c_vec, b.S_c_vec)) && all_passed;

c1 = tiny_model(args{:}, 'S_c_seed', 11); c1.build();
c2 = tiny_model(args{:}, 'S_c_seed', 12); c2.build();
all_passed = check('a different S_c_seed gives a different draw', ...
    ~isequal(c1.S_c_vec, c2.S_c_vec)) && all_passed;

d1 = tiny_model(args{:}, 'rng_seeds', [5 6]); d1.build();
d2 = tiny_model(args{:}, 'rng_seeds', [7 8]); d2.build();
all_passed = check('a different network seed also moves the setpoints', ...
    ~isequal(d1.S_c_vec, d2.S_c_vec)) && all_passed;

e1 = tiny_model(args{:}, 'rng_seeds', [5 6], 'S_c_seed', 99); e1.build();
e2 = tiny_model(args{:}, 'rng_seeds', [7 8], 'S_c_seed', 99); e2.build();
all_passed = check('an explicit S_c_seed pins the draw across networks', ...
    isequal(e1.S_c_vec, e2.S_c_vec) && ~isequal(full(e1.W), full(e2.W))) && all_passed;

%% The draw does not disturb W, the stimulus, or x0
plain = tiny_model('rng_seeds', [3 4]);              plain.build();
het   = tiny_model('rng_seeds', [3 4], args{:});     het.build();
all_passed = check('W is bit-identical with and without the draw', ...
    max(abs(full(plain.W) - full(het.W)), [], 'all') == 0) && all_passed;
all_passed = check('S0 is bit-identical with and without the draw', ...
    max(abs(plain.S0 - het.S0)) == 0) && all_passed;
all_passed = check('the external input is bit-identical too', ...
    max(abs(plain.u_ex - het.u_ex), [], 'all') == 0) && all_passed;

%% phi with a vector centre == elementwise phi with a scalar centre
for nm = {'logistic', 'piecewise'}
    m = tiny_model('activation', nm{1}, 'n', 30, 'indegree', 10, ...
        'mu_S_c_E', 0.4, 'sigma_S_c_E', 0.1, 'mu_S_c_I', 0.1, 'sigma_S_c_I', 0.1);
    m.build();
    cv = m.S_c_vec;
    x = linspace(-1.5, 2.5, m.n)';

    ref  = zeros(m.n, 1);
    dref = zeros(m.n, 1);
    for i = 1:m.n
        switch nm{1}
            case 'logistic'
                ref(i)  = SRNNModel2.logisticSigmoid(x(i), cv(i));
                dref(i) = SRNNModel2.logisticSigmoidDerivative(x(i), cv(i));
            case 'piecewise'
                ref(i)  = SRNNModel2.piecewiseSigmoid(x(i), m.S_a, cv(i));
                dref(i) = SRNNModel2.piecewiseSigmoidDerivative(x(i), m.S_a, cv(i));
        end
    end
    all_passed = check(sprintf('%s: vector centre == elementwise scalar centre', nm{1}), ...
        same(m.activation_function(x), ref) && ...
        same(m.activation_function_derivative(x), dref)) && all_passed;

    % The plotting path evaluates phi on an n x nt trajectory, so the same
    % must hold column by column under implicit expansion.
    X = x + [0, 0.25, -0.4];
    Y = m.activation_function(X);
    okmat = isequal(size(Y), size(X));
    for col = 1:size(X, 2)
        okmat = okmat && same(Y(:, col), m.activation_function(X(:, col)));
    end
    all_passed = check(sprintf('%s: broadcasts across an n x nt trajectory', nm{1}), ...
        okmat) && all_passed;
end

%% The analytic Jacobian still matches finite differences
m = SRNN2TestAccess('n', 40, 'indegree', 12, 'n_a_E', 2, 'n_b_E', 1, ...
    'n_a_I', 1, 'n_b_I', 1, 'f', 0.6, 'T_range', [0 0.5], 'fs', 200, ...
    'ode_solver', @ode_rk4, 'lya_method', 'none', ...
    'mu_S_c_E', 0.4, 'sigma_S_c_E', 0.1, 'mu_S_c_I', 0.1, 'sigma_S_c_I', 0.1);
m.build();
params = m.get_params();
params.u_interpolant = m.u_interpolant;
S = m.S0 + 0.05 * randn(size(m.S0));
J_analytic = full(SRNNModel2.compute_Jacobian_fast(S, params));
J_fd = fd_jacobian(@(s) SRNN2TestAccess.rhs(0.1, s, params), S, 1e-6);
jac_err = max(abs(J_analytic(:) - J_fd(:)));
fprintf('  (Jacobian error with heterogeneous S_c: %.3e)\n', jac_err);
all_passed = check('analytic Jacobian matches finite differences', ...
    jac_err < 2e-6) && all_passed;

%% Nonlinearities without a centre reject the request
[threw, err] = capture(@() build_it(tiny_model('activation', 'tanh', 'mu_S_c_E', 0.4)));
all_passed = check('tanh + mu_S_c errors, saying why', ...
    threw && contains(err.message, 'centre') && contains(err.message, 'tanh')) && all_passed;

[threw, err] = capture(@() build_it(tiny_model( ...
    'activation_custom', {@(x) 2*x, @(x) 2*ones(size(x))}, 'sigma_S_c_E', 0.1)));
all_passed = check('activation_custom + sigma_S_c errors', ...
    threw && contains(err.message, 'activation_custom')) && all_passed;

[threw, ~] = capture(@() build_it(tiny_model('sigma_S_c_E', -0.1)));
all_passed = check('a negative sigma_S_c errors', threw) && all_passed;

%% A heterogeneous model runs end to end and its LLE differs
m = tiny_model('n', 60, 'indegree', 20, 'n_a_E', 1, 'n_b_E', 1, ...
    'T_range', [0 2], 'fs', 200, 'ode_solver', @ode_rk4, ...
    'lya_method', 'benettin', 'lya_T_interval', [0.5 2]);
m.build(); m.run();
lle_homog = m.lya_results.LLE;

m = tiny_model('n', 60, 'indegree', 20, 'n_a_E', 1, 'n_b_E', 1, ...
    'T_range', [0 2], 'fs', 200, 'ode_solver', @ode_rk4, ...
    'lya_method', 'benettin', 'lya_T_interval', [0.5 2], ...
    'mu_S_c_E', 0.4, 'mu_S_c_I', 0.1, 'sigma_S_c_E', 0.1, 'sigma_S_c_I', 0.1);
m.build(); m.run();
all_passed = check('a heterogeneous model runs and its LLE moves', ...
    isfinite(m.lya_results.LLE) && ...
    abs(m.lya_results.LLE - lle_homog) > 1e-6) && all_passed;

%% Result
fprintf('\n========================================\n');
if all_passed
    fprintf('ALL TESTS PASSED!\n');
else
    fprintf('SOME TESTS FAILED!\n');
end
fprintf('========================================\n');

%% Helpers
function m = tiny_model(varargin)
% A small, quiet model: everything here is about the setpoints, not the network.
defaults = {'n', 40, 'indegree', 12, 'f', 0.5, 'T_range', [0 0.5], ...
    'fs', 200, 'ode_solver', @ode_rk4, 'lya_method', 'none'};
m = SRNNModel2(defaults{:}, varargin{:});
end

function build_it(m)
m.build();
end

function J = fd_jacobian(rhs, S, h)
% Central differences, column by column.
N = numel(S);
J = zeros(N, N);
for j = 1:N
    Sp = S; Sp(j) = Sp(j) + h;
    Sm = S; Sm(j) = Sm(j) - h;
    J(:, j) = (rhs(Sp) - rhs(Sm)) / (2 * h);
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
