% test_SRNNCellTypePairs.m - Pair-specific STD/STF dynamics and Jacobian.

clear; clc;
rng(0);

model = make_pair_model('none');
model.build();
params = model.get_params();
params.u_interpolant = model.u_interpolant;

%% Configuration compilation and state layout
assert(isequal(model.n_b_pairs, [2 0 1; 0 1 0; 0 0 0]));
assert(isequal(model.n_g_pairs, [0 1 2; 0 0 0; 0 1 0]));
assert(isequal(params.tau_b_rel{1, 1}, [0.2 0.2]));
assert(isequal(params.tau_g_fac{1, 3}, [0.2 0.3]));
expected = model.n + sum(model.n_per_type .* model.n_a) + ...
    sum(model.n_per_type .* sum(model.n_b_pairs + model.n_g_pairs, 2)');
assert(model.N_sys_eqs == expected);
assert(numel(model.S0) == expected);
for pre = 1:model.n_cellTypes
    for post = 1:model.n_cellTypes
        assert(all(model.S0(params.state_layout.b{pre, post}) == 1));
        assert(all(model.S0(params.state_layout.g{pre, post}) == 1));
    end
end

%% Analytic Jacobian against central finite differences
S = random_state(model.S0, params);
for floor_on = [false true]
    params.std_zero_floor = floor_on;
    J_analytic = full(SRNNCellTypePairs.compute_Jacobian_fast(S, params));
    J_fd = finite_difference_jacobian(S, params, 1e-6);
    max_error = max(abs(J_analytic(:) - J_fd(:)));
    fprintf('Pair Jacobian error (std_zero_floor=%d): %.3e\n', ...
        floor_on, max_error);
    assert(max_error < 2e-6);
end

%% Route readouts use independent products
params.std_zero_floor = false;
[~, ~, b, r, synaptic_output, g] = ...
    SRNNCellTypePairs.unpack_and_compute_states(S', params);
pre_name = 'E'; post_name = 'SST';
b_product = reshape(prod(b.(pre_name).(post_name), 2), ...
    model.n_per_type(1), []);
g_product = reshape(prod(g.(pre_name).(post_name), 2), ...
    model.n_per_type(1), []);
expected_output = b_product .* g_product .* r.(pre_name);
assert(max(abs(synaptic_output.(pre_name).(post_name) - expected_output), ...
    [], 'all') < 1e-12);
assert(isempty(b.E.PV));
assert(isempty(g.E.E));

%% Pair routing and dead-state triangular structure
routing = params;
routing.W = sparse(routing.n, routing.n);
E = routing.type_indices{1};
PV = routing.type_indices{2};
routing.W(PV, E) = 0.1;
J = full(SRNNCellTypePairs.compute_Jacobian_fast(S, routing));
row_x = routing.state_layout.x;
row_b_EE = routing.state_layout.b{1, 1};
row_g_EPV = routing.state_layout.g{1, 2};
assert(max(abs(J(row_x, row_b_EE)), [], 'all') < 1e-12);
assert(max(abs(J(row_x(setdiff(1:routing.n, PV)), row_g_EPV)), [], 'all') < 1e-12);
assert(any(abs(J(row_x(PV), row_g_EPV)) > 0, 'all'));
assert(all(diag(J(row_b_EE, row_b_EE)) < 0));

%% End-to-end integration and named outputs
model.run();
assert(model.has_run);
assert(all(isfinite(model.S_out), 'all'));
assert(isfield(model.plot_data.synaptic_output.E, 'PV'));
assert(size(model.plot_data.b.E.E, 2) == 2);
assert(size(model.plot_data.g.E.SST, 2) == 2);
old_visibility = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
visibility_cleanup = onCleanup(@() set(groot, ...
    'DefaultFigureVisible', old_visibility));
[figure_handle, ~] = model.plot(); close(figure_handle);
[figure_handle, celltype_axes] = model.plot_celltypes();
synaptic_PV_lines = findall(celltype_axes(4, 1), ...
    'Type', 'line', 'Tag', 'PostType_PV');
std_E_lines = findall(celltype_axes(6, 1), ...
    'Type', 'line', 'Tag', 'PostType_E');
std_SST_lines = findall(celltype_axes(6, 1), ...
    'Type', 'line', 'Tag', 'PostType_SST');
stf_PV_lines = findall(celltype_axes(7, 1), ...
    'Type', 'line', 'Tag', 'PostType_PV');
assert(numel(synaptic_PV_lines) == model.n_per_type(1));
assert(numel(std_E_lines) == model.n_per_type(1));
assert(numel(std_SST_lines) == model.n_per_type(1));
assert(numel(stf_PV_lines) == model.n_per_type(1));
std_E_colors = vertcat(std_E_lines.Color);
assert(size(unique(round(std_E_colors, 8), 'rows'), 1) > 1);
post_colors = lines(model.n_cellTypes);
expected_E_colors = 0.5 .* repmat(post_colors(1, :), ...
    model.n_per_type(1), 1) + 0.5 .* lines(model.n_per_type(1));
assert(max(abs(sortrows(std_E_colors) - sortrows(expected_E_colors)), ...
    [], 'all') < 1e-12);
close(figure_handle);

%% Uniform routes reproduce SRNNCellTypes recurrent dynamics
[old_model, pair_model] = make_parity_models();
old_model.build(); pair_model.build();
old_model.run(); pair_model.run();
old_x = old_model.S_out(:, old_model.get_params().state_layout.x);
pair_x = pair_model.S_out(:, pair_model.get_params().state_layout.x);
assert(max(abs(old_x - pair_x), [], 'all') < 2e-6);

%% Benettin and QR Lyapunov workflows
benettin_model = make_lya_model('benettin');
benettin_model.build(); benettin_model.run();
assert(isfinite(benettin_model.lya_results.LLE));
qr_model = make_lya_model('qr');
qr_model.build(); qr_model.run();
assert(numel(qr_model.lya_results.LE_spectrum) == qr_model.N_sys_eqs);
assert(all(isfinite(qr_model.lya_results.LE_spectrum)));

%% Validation failures
assert_throws(@() bad_unknown_type());
assert_throws(@() bad_missing_field());
assert_throws(@() bad_vector_length());
assert_throws(@() bad_G());

%% Empty configuration allocates no synaptic states
empty_model = SRNNCellTypePairs('n', 6, 'indegree', 3, ...
    'n_cellTypes', 2, 'cell_type_names', {'E', 'I'}, 'f', [0.5 0.5], ...
    'mu_tilde', [0.1 -0.1], 'sigma_tilde', [0.01 0.01], ...
    'T_range', [0 0.05], 'fs', 200, 'lya_method', 'none', ...
    'store_full_state', true);
assert(empty_model.N_sys_eqs == empty_model.n);
empty_model.build(); empty_model.run();
assert(all(isfinite(empty_model.S_out), 'all'));

fprintf('test_SRNNCellTypePairs: ALL TESTS PASSED\n');

function model = make_pair_model(lya_method)
synapse_config = struct();
synapse_config.E.E.std = struct('tau_rec', [0.3 1], 'tau_rel', 0.2);
synapse_config.E.PV.stf = struct('tau_dec', 0.5, 'tau_fac', 0.25, 'G', 2);
synapse_config.E.SST.std = struct('tau_rec', 0.8, 'tau_rel', 0.25);
synapse_config.E.SST.stf = struct( ...
    'tau_dec', [0.4 1.2], 'tau_fac', [0.2 0.3], 'G', [1.5 2]);
synapse_config.PV.PV.std = struct('tau_rec', 0.6, 'tau_rel', 0.3);
synapse_config.SST.PV.stf = struct('tau_dec', 0.9, 'tau_fac', 0.4, 'G', 1.8);
model = SRNNCellTypePairs( ...
    'n', 12, 'indegree', 4, ...
    'n_cellTypes', 3, 'cell_type_names', {'E', 'PV', 'SST'}, ...
    'f', [0.5 0.3 0.2], ...
    'mu_tilde', [0.1 -0.2 -0.08], ...
    'sigma_tilde', [0.01 0.02 0.015], ...
    'n_a', [2 0 1], 'tau_a', {[0.25 1], [], 2}, ...
    'c', [0.05 0 0.03], 'synapse_config', synapse_config, ...
    'T_range', [0 0.1], 'fs', 200, 'lya_method', lya_method, ...
    'store_full_state', true);
end

function S = random_state(S0, params)
S = S0;
for q = 1:params.n_cellTypes
    S(params.state_layout.a{q}) = 0.3 * rand(numel(params.state_layout.a{q}), 1);
end
for pre = 1:params.n_cellTypes
    npre = params.n_per_type(pre);
    for post = 1:params.n_cellTypes
        row_b = params.state_layout.b{pre, post};
        S(row_b) = 0.4 + 0.6 * rand(numel(row_b), 1);
        ng = params.n_g_pairs(pre, post);
        if ng > 0
            G = repmat(params.G{pre, post}, npre, 1);
            row_g = params.state_layout.g{pre, post};
            S(row_g) = 1 + (G(:) - 1) .* rand(numel(row_g), 1);
        end
    end
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
    J(:, k) = (SRNNCellTypePairs.dynamics_fast(0, plus, params) - ...
        SRNNCellTypePairs.dynamics_fast(0, minus, params)) / (2 * h);
end
end

function [old_model, pair_model] = make_parity_models()
common = {'n', 8, 'indegree', 4, 'n_cellTypes', 2, ...
    'cell_type_names', {'E', 'I'}, 'f', [0.5 0.5], ...
    'mu_tilde', [0.1 -0.1], 'sigma_tilde', [0.01 0.01], ...
    'n_a', [1 0], 'tau_a', {0.5, []}, 'c', [0.05 0], ...
    'T_range', [0 0.1], 'fs', 200, 'lya_method', 'none', ...
    'store_full_state', true, 'rng_seeds', [7 8]};
old_model = SRNNCellTypes(common{:}, ...
    'n_b', [1 0], 'tau_b_rec', {0.8, []}, 'tau_b_rel', [0.25 0.25], ...
    'n_g', [0 1], 'tau_g_dec', {[], 0.6}, ...
    'tau_g_fac', [0.25 0.3], 'G', [2 1.8]);
config = struct();
for post = {'E', 'I'}
    config.E.(post{1}).std = struct('tau_rec', 0.8, 'tau_rel', 0.25);
    config.I.(post{1}).stf = struct('tau_dec', 0.6, 'tau_fac', 0.3, 'G', 1.8);
end
pair_model = SRNNCellTypePairs(common{:}, 'synapse_config', config);
end

function model = make_lya_model(method)
config.E.I.std = struct('tau_rec', 0.8, 'tau_rel', 0.25);
config.I.E.stf = struct('tau_dec', 0.6, 'tau_fac', 0.3, 'G', 1.8);
model = SRNNCellTypePairs('n', 6, 'indegree', 3, ...
    'n_cellTypes', 2, 'cell_type_names', {'E', 'I'}, 'f', [0.5 0.5], ...
    'mu_tilde', [0.1 -0.1], 'sigma_tilde', [0.01 0.01], ...
    'synapse_config', config, 'T_range', [0 0.2], 'fs', 200, ...
    'lya_method', method, 'store_full_state', true, ...
    'store_decimated_state', false);
end

function bad_unknown_type()
config.X.E.std = struct('tau_rec', 1, 'tau_rel', 0.25);
make_minimal(config);
end

function bad_missing_field()
config.E.E.std = struct('tau_rec', 1);
make_minimal(config);
end

function bad_vector_length()
config.E.E.std = struct('tau_rec', [1 2], 'tau_rel', [0.2 0.3 0.4]);
make_minimal(config);
end

function bad_G()
config.E.E.stf = struct('tau_dec', 1, 'tau_fac', 0.2, 'G', 0.9);
make_minimal(config);
end

function model = make_minimal(config)
model = SRNNCellTypePairs('n_cellTypes', 2, ...
    'cell_type_names', {'E', 'I'}, 'f', [0.5 0.5], ...
    'mu_tilde', [0.1 -0.1], 'sigma_tilde', [0.01 0.01], ...
    'synapse_config', config);
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
