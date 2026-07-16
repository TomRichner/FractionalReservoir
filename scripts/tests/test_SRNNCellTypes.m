% test_SRNNCellTypes.m - Verify generalized cell types, dynamics, and Jacobian.

clear; clc;
addpath(fileparts(fileparts(mfilename('fullpath'))));
setup_paths();
rng(0);

model = make_three_type_model();
model.build();
params = model.get_params();
params.u_interpolant = model.u_interpolant;

%% Layout and named interfaces
expected_states = model.n + sum(model.n_per_type .* (model.n_a + model.n_b));
assert(model.N_sys_eqs == expected_states);
assert(numel(model.S0) == expected_states);
assert(all(isfield(model.cell_indices, {'E', 'PV', 'SST'})));
assert(isequal(model.n_per_type, [6 4 2]));
assert(issparse(model.W));

%% Analytic Jacobian against central finite differences
S = random_state(model.S0, params);
for floor_on = [false true]
    params.std_zero_floor = floor_on;
    J_analytic = full(SRNNCellTypes.compute_Jacobian_fast(S, params));
    J_fd = finite_difference_jacobian(S, params, 1e-6);
    max_error = max(abs(J_analytic(:) - J_fd(:)));
    fprintf('Jacobian error (std_zero_floor=%d): %.3e\n', floor_on, max_error);
    assert(max_error < 1e-6);
end

%% End-to-end integration and named plot data
model.run();
assert(model.has_run);
assert(all(isfinite(model.S_out), 'all'));
assert(all(isfield(model.plot_data.r, {'E', 'PV', 'SST'})));
assert(size(model.plot_data.a.E, 2) == 2);
assert(isempty(model.plot_data.a.PV));
assert(size(model.plot_data.b.PV, 2) == 2);

%% Validation of required names and ragged timescales
assert_throws(@() SRNNCellTypes('n_cellTypes', 2, ...
    'cell_type_names', {'A', 'A'}, 'f', [0.5 0.5], ...
    'mu_tilde', [1 -1], 'sigma_tilde', [1 1]));
assert_throws(@() bad_tau_model());

fprintf('test_SRNNCellTypes: ALL TESTS PASSED\n');

function model = make_three_type_model()
model = SRNNCellTypes( ...
    'n', 12, 'indegree', 4, ...
    'n_cellTypes', 3, ...
    'cell_type_names', {'E', 'PV', 'SST'}, ...
    'f', [0.5 0.3 0.2], ...
    'mu_tilde', [0.1 -0.2 -0.08], ...
    'sigma_tilde', [0.01 0.02 0.015], ...
    'n_a', [2 0 1], ...
    'tau_a', {[0.25 1], [], 2}, ...
    'c', [0.05 0 0.03], ...
    'n_b', [1 2 0], ...
    'tau_b_rec', {1, [0.5 2], []}, ...
    'tau_b_rel', [0.25 0.4 0.3], ...
    'T_range', [0 0.2], 'fs', 100, ...
    'lya_method', 'none', 'store_full_state', true);
end

function S = random_state(S0, params)
S = S0;
for q = 1:params.n_cellTypes
    S(params.state_layout.a{q}) = 0.3 * rand(numel(params.state_layout.a{q}), 1);
    S(params.state_layout.b{q}) = 0.4 + ...
        0.6 * rand(numel(params.state_layout.b{q}), 1);
end
S(params.state_layout.x) = 0.5 * randn(params.n, 1);
end

function J = finite_difference_jacobian(S, params, h)
N = numel(S);
J = zeros(N, N);
for k = 1:N
    plus = S; minus = S;
    plus(k) = plus(k) + h;
    minus(k) = minus(k) - h;
    J(:, k) = (SRNNCellTypes.dynamics_fast(0, plus, params) - ...
        SRNNCellTypes.dynamics_fast(0, minus, params)) / (2 * h);
end
end

function bad_tau_model()
model = SRNNCellTypes('n_cellTypes', 2, ...
    'cell_type_names', {'A', 'B'}, 'f', [0.5 0.5], ...
    'mu_tilde', [1 -1], 'sigma_tilde', [1 1], ...
    'n_b', [2 0], 'tau_b_rec', {1, []});
model.build();
end

function assert_throws(f)
threw = false;
try
    f();
catch
    threw = true;
end
assert(threw, 'Expected operation to throw an error.');
end
