% test_SRNNCellTypePairs_S_c_heterogeneity.m
% Verify per-neuron nonlinearity setpoints in SRNNCellTypePairs: mu_S_c and
% sigma_S_c are 1 x C rows, drawn at build() into the read-only S_c_vec, and
% the activation handles then carry a per-neuron centre.
%
% Uses THREE cell types throughout, so nothing here can accidentally depend on
% the two-type case the way an E/I-shaped API would.
%
% See also: SRNNCellTypePairs, test_SRNN2_S_c_heterogeneity, test_SRNNCellTypePairs

fprintf('=== Testing SRNNCellTypePairs per-neuron S_c ===\n\n');
all_passed = true;
xg = linspace(-1, 2, 101)';

%% Defaults are inert: no vector, and the scalar path is untouched
m_base = tiny_model();
m_base.build();
all_passed = check('S_c_vec is empty with the properties at their defaults', ...
    isempty(m_base.S_c_vec)) && all_passed;
all_passed = check('the handle is still the scalar-S_c function', ...
    same(m_base.activation_function(xg), ...
         SRNNCellTypePairs.logisticSigmoid(xg, m_base.S_c))) && all_passed;

%% Distinct means, no spread: one value per type, in the right blocks
m = tiny_model('mu_S_c', [0.40 0.10 -0.20]);
m.build();
v = m.S_c_vec;
all_passed = check('S_c_vec is an n x 1 column', ...
    isequal(size(v), [m.n, 1])) && all_passed;
per_type_ok = true;
for q = 1:m.n_cellTypes
    expected = [0.40 0.10 -0.20];
    per_type_ok = per_type_ok && all(v(m.type_indices{q}) == expected(q));
end
all_passed = check('each type gets its own mu_S_c', per_type_ok) && all_passed;
all_passed = check('sigma = 0 gives no spread at all', ...
    numel(unique(v)) == m.n_cellTypes) && all_passed;

%% A scalar broadcasts to every type; an empty mu falls back to S_c
m = tiny_model('mu_S_c', 0.25);
m.build();
all_passed = check('a scalar mu_S_c broadcasts to every type', ...
    isequal(m.mu_S_c, [0.25 0.25 0.25]) && all(m.S_c_vec == 0.25)) && all_passed;

m = tiny_model('S_c', 0.31, 'sigma_S_c', [0 0 0.1]);
m.build();
all_passed = check('an empty mu_S_c falls back to S_c', ...
    all(m.S_c_vec(m.type_indices{1}) == 0.31) && ...
    all(m.S_c_vec(m.type_indices{2}) == 0.31)) && all_passed;
all_passed = check('a scalar sigma_S_c is accepted and per-type sigma bites', ...
    std(m.S_c_vec(m.type_indices{3})) > 0 && ...
    std(m.S_c_vec(m.type_indices{1})) == 0) && all_passed;

%% The requested statistics come out of the draw
m = tiny_model('n', 900, 'indegree', 100, ...
    'mu_S_c', [0.40 0.10 -0.20], 'sigma_S_c', [0.05 0.20 0.10]);
m.build();
mu_want = [0.40 0.10 -0.20];
sd_want = [0.05 0.20 0.10];
mean_ok = true; sd_ok = true;
for q = 1:m.n_cellTypes
    vq = m.S_c_vec(m.type_indices{q});
    nq = numel(vq);
    mean_ok = mean_ok && abs(mean(vq) - mu_want(q)) < 3 * sd_want(q) / sqrt(nq);
    sd_ok   = sd_ok   && abs(std(vq)  - sd_want(q)) < 3 * sd_want(q) / sqrt(2*nq);
end
all_passed = check('every type''s sample mean matches its mu_S_c', mean_ok) && all_passed;
all_passed = check('every type''s sample SD matches its sigma_S_c', sd_ok) && all_passed;

%% Reproducibility and seed control
args = {'mu_S_c', [0.4 0.1 0.1], 'sigma_S_c', 0.05};
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

%% The draw does not disturb W, the stimulus, or x0
plain = tiny_model('rng_seeds', [3 4]);           plain.build();
het   = tiny_model('rng_seeds', [3 4], args{:});  het.build();
all_passed = check('W is bit-identical with and without the draw', ...
    max(abs(full(plain.W) - full(het.W)), [], 'all') == 0) && all_passed;
all_passed = check('S0 is bit-identical with and without the draw', ...
    max(abs(plain.S0 - het.S0)) == 0) && all_passed;
all_passed = check('the external input is bit-identical too', ...
    max(abs(plain.u_ex - het.u_ex), [], 'all') == 0) && all_passed;

%% phi with a vector centre == elementwise phi with a scalar centre
for nm = {'logistic', 'piecewise'}
    m = tiny_model('activation', nm{1}, ...
        'mu_S_c', [0.4 0.1 -0.2], 'sigma_S_c', [0.1 0.1 0.1]);
    m.build();
    cv = m.S_c_vec;
    x = linspace(-1.5, 2.5, m.n)';

    ref  = zeros(m.n, 1);
    dref = zeros(m.n, 1);
    for i = 1:m.n
        switch nm{1}
            case 'logistic'
                ref(i)  = SRNNCellTypePairs.logisticSigmoid(x(i), cv(i));
                dref(i) = SRNNCellTypePairs.logisticSigmoidDerivative(x(i), cv(i));
            case 'piecewise'
                ref(i)  = SRNNCellTypePairs.piecewiseSigmoid(x(i), m.S_a, cv(i));
                dref(i) = SRNNCellTypePairs.piecewiseSigmoidDerivative(x(i), m.S_a, cv(i));
        end
    end
    all_passed = check(sprintf('%s: vector centre == elementwise scalar centre', nm{1}), ...
        same(m.activation_function(x), ref) && ...
        same(m.activation_function_derivative(x), dref)) && all_passed;

    % The plotting path evaluates phi on an n x nt trajectory.
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
sc = struct();
sc.A.A.std = struct('tau_rec', 1, 'tau_rel', 0.25);
sc.B.A.stf = struct('tau_dec', 1, 'tau_fac', 0.5, 'G', 1.5);
sc.C.B.std = struct('tau_rec', 2, 'tau_rel', 0.5);
m = tiny_model('n', 30, 'indegree', 10, 'n_a', [2 1 0], ...
    'synapse_config', sc, ...
    'mu_S_c', [0.4 0.1 -0.2], 'sigma_S_c', [0.1 0.1 0.1]);
m.build();
params = m.get_params();
params.u_interpolant = m.u_interpolant;
S = m.S0 + 0.05 * randn(size(m.S0));
J_analytic = full(SRNNCellTypePairs.compute_Jacobian_fast(S, params));
J_fd = fd_jacobian(params, S, 1e-6);
jac_err = max(abs(J_analytic(:) - J_fd(:)));
fprintf('  (Jacobian error with heterogeneous S_c: %.3e)\n', jac_err);
all_passed = check('analytic Jacobian matches finite differences', ...
    jac_err < 2e-6) && all_passed;

%% Nonlinearities without a centre reject the request
[threw, err] = capture(@() build_it(tiny_model('activation', 'tanh', 'mu_S_c', 0.4)));
all_passed = check('tanh + mu_S_c errors, saying why', ...
    threw && contains(err.message, 'centre') && contains(err.message, 'tanh')) && all_passed;

[threw, err] = capture(@() build_it(tiny_model( ...
    'activation_custom', {@(x) 2*x, @(x) 2*ones(size(x))}, 'sigma_S_c', 0.1)));
all_passed = check('activation_custom + sigma_S_c errors', ...
    threw && contains(err.message, 'activation_custom')) && all_passed;

[threw, ~] = capture(@() build_it(tiny_model('sigma_S_c', [0.1 -0.1 0])));
all_passed = check('a negative sigma_S_c errors', threw) && all_passed;

[threw, ~] = capture(@() build_it(tiny_model('mu_S_c', [0.1 0.2])));
all_passed = check('a wrong-length mu_S_c errors', threw) && all_passed;

%% A heterogeneous model runs end to end and its LLE differs
run_args = {'T_range', [0 2], 'fs', 200, 'ode_solver', @ode_rk4, ...
    'lya_method', 'benettin', 'lya_T_interval', [0.5 2], 'n_a', [2 1 0]};
m = tiny_model(run_args{:});
m.build(); m.run();
lle_homog = m.lya_results.LLE;

m = tiny_model(run_args{:}, 'mu_S_c', [0.4 0.1 -0.2], 'sigma_S_c', [0.1 0.1 0.1]);
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
% Three cell types, small and quiet: this test is about the setpoints.
defaults = {'n', 30, 'indegree', 10, ...
    'n_cellTypes', 3, 'cell_type_names', {'A', 'B', 'C'}, ...
    'f', [1/3 1/3 1/3], ...
    'mu_tilde_relative',    [3 -4 -2; 2 -3 -1; 3 -2 -2], ...
    'sigma_tilde_relative', ones(3), ...
    'T_range', [0 0.5], 'fs', 200, 'ode_solver', @ode_rk4, 'lya_method', 'none'};
m = SRNNCellTypePairs(defaults{:}, varargin{:});
end

function build_it(m)
m.build();
end

function J = fd_jacobian(params, S, h)
% Central differences, column by column.
N = numel(S);
J = zeros(N, N);
for k = 1:N
    plus = S;  plus(k)  = S(k) + h;
    minus = S; minus(k) = S(k) - h;
    J(:, k) = (SRNNCellTypePairs.dynamics_fast(0.1, plus, params) - ...
        SRNNCellTypePairs.dynamics_fast(0.1, minus, params)) / (2 * h);
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
