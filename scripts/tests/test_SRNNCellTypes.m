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
expected_states = model.n + sum(model.n_per_type .* ...
    (model.n_a + model.n_b + model.n_g));
assert(model.N_sys_eqs == expected_states);
assert(numel(model.S0) == expected_states);
assert(all(isfield(model.cell_indices, {'E', 'PV', 'SST'})));
assert(isequal(model.n_per_type, [6 4 2]));
assert(issparse(model.W));
for q = 1:model.n_cellTypes
    assert(all(model.S0(params.state_layout.g{q}) == 1));
end
legacy_layout = SRNNCellTypes.make_state_layout( ...
    model.n, model.n_per_type, model.n_a, model.n_b);
explicit_layout = SRNNCellTypes.make_state_layout( ...
    model.n, model.n_per_type, model.n_a, model.n_b, zeros(1, model.n_cellTypes));
assert(isequal(legacy_layout, explicit_layout));

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

%% Combined STF and STD synaptic output
params.std_zero_floor = false;
[~, ~, b_named, r_named, br_named, g_named] = ...
    SRNNCellTypes.unpack_and_compute_states(S', params);
for q = 1:model.n_cellTypes
    name = model.cell_type_names{q};
    nq = model.n_per_type(q);
    b_product = reshape(prod(b_named.(name), 2), nq, []);
    g_product = reshape(prod(g_named.(name), 2), nq, []);
    expected_br = g_product .* b_product .* r_named.(name);
    assert(max(abs(br_named.(name) - expected_br), [], 'all') < 1e-12);
end

%% End-to-end integration and named plot data
model.run();
assert(model.has_run);
assert(all(isfinite(model.S_out), 'all'));
assert(all(isfield(model.plot_data.r, {'E', 'PV', 'SST'})));
assert(size(model.plot_data.a.E, 2) == 2);
assert(isempty(model.plot_data.a.PV));
assert(size(model.plot_data.b.PV, 2) == 2);
assert(size(model.plot_data.g.E, 2) == 2);
assert(size(model.plot_data.g.SST, 2) == 1);
assert(all(isfinite(model.plot_data.br.E), 'all'));

%% STF defaults and disabled-state regression
default_stf = make_minimal_model('n_g', [1 0]);
assert(isequal(default_stf.tau_g_dec, {1, []}));
assert(isequal(default_stf.tau_g_fac, [0.25 0.25]));
assert(isequal(default_stf.G, [2 2]));
disabled_stf = make_minimal_model();
assert(disabled_stf.N_sys_eqs == disabled_stf.n);
disabled_stf.build();
assert(numel(disabled_stf.S0) == disabled_stf.n);
disabled_stf.run();
assert(all(isfinite(disabled_stf.S_out), 'all'));
assert(all(disabled_stf.plot_data.g.A == 1, 'all'));

%% Validation of required names and ragged timescales
assert_throws(@() SRNNCellTypes('n_cellTypes', 2, ...
    'cell_type_names', {'A', 'A'}, 'f', [0.5 0.5], ...
    'mu_tilde', [1 -1], 'sigma_tilde', [1 1]));
assert_throws(@() bad_tau_model());
assert_throws(@() make_minimal_model('n_g', [2 0], 'tau_g_dec', {1, []}));
assert_throws(@() make_minimal_model('n_g', [0.5 0]));
assert_throws(@() make_minimal_model('n_g', [1 0], 'tau_g_dec', {-1, []}));
assert_throws(@() make_minimal_model('tau_g_fac', [0 0.25]));
assert_throws(@() make_minimal_model('G', [0.9 2]));

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
    'n_g', [2 0 1], ...
    'tau_g_dec', {[0.3 1.5], [], 0.8}, ...
    'tau_g_fac', [0.2 0.25 0.4], ...
    'G', [2.5 2 1.8], ...
    'T_range', [0 0.2], 'fs', 100, ...
    'lya_method', 'none', 'store_full_state', true);
end

function S = random_state(S0, params)
S = S0;
for q = 1:params.n_cellTypes
    S(params.state_layout.a{q}) = 0.3 * rand(numel(params.state_layout.a{q}), 1);
    S(params.state_layout.b{q}) = 0.4 + ...
        0.6 * rand(numel(params.state_layout.b{q}), 1);
    if params.n_g(q) > 0
        S(params.state_layout.g{q}) = 1 + (params.G(q) - 1) .* ...
            rand(numel(params.state_layout.g{q}), 1);
    end
end
S(params.state_layout.x) = 0.5 * randn(params.n, 1);
end

function model = make_minimal_model(varargin)
model = SRNNCellTypes('n_cellTypes', 2, ...
    'n', 12, 'indegree', 4, ...
    'cell_type_names', {'A', 'B'}, 'f', [0.5 0.5], ...
    'mu_tilde', [0.1 -0.1], 'sigma_tilde', [0.01 0.01], ...
    'T_range', [0 0.05], 'fs', 100, 'lya_method', 'none', ...
    'store_full_state', true, varargin{:});
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
