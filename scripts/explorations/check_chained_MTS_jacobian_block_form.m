% CHECK_CHAINED_MTS_JACOBIAN_BLOCK_FORM  The block form of the Jacobian in
% docs/EquationsParametersDocs/chained_MTS/Equations_chained_MTS.md,
% section 10, against SRNNCellTypePairs.compute_Jacobian_fast.
%
% The document writes the Jacobian of the simplified model (one set of
% constants for all neurons and routes, no spreads, depression on every
% route) as
%
%   J = D + U Phi' V^T                                   (equation (24))
%
% on the state z = [x; a_1; ...; a_K; b_1; b_2], N values each. This script
% assembles that J from the formulas of the document, coded here
% independently of the model class (the logistic nonlinearity and every
% block are written out below), and compares it with the code's sparse
% Jacobian:
%
%   1. One cell type, one route, n = 12: the code's state is a permutation
%      of z, so the comparison is elementwise, J_code = T J T^T. Four
%      conditions: no adaptation, single timescale, chained MTS, product
%      MTS. Two random states each.
%   2. Two cell types (E, I), four chain routes, n = 16: a neuron's two
%      routes hold two identical copies of its pools, so the code's state
%      is larger than z. With T the map that writes each pool into both
%      copies, the block form is the code's Jacobian on the subspace of
%      equal copies: J_code T = T J.
%   3. The reduction (27) to the dendritic states: each eigenvalue s of J
%      (chained MTS) makes (1 + s tau_d) I - W G(s) singular, G(s) the
%      diagonal matrix of the per-neuron gains (28).
%
% Assumes the session is set up (setup_paths). No pool, a few seconds.
% Prints PASS / FAIL per check and the largest difference of each.

tol = 1e-12;
all_passed = true;

tau_d = 0.1;
S_c = 0.2;
c = 0.5;
ladder = [0.25, sqrt(2.5), 10];
single = struct('tau_rec', 2, 'tau_rel', 0.25);
chain = struct('model', 'chain', 'tau_rec', [2 5], 'tau_rel', 0.25, 'n2', 2);
product = struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]);

conditions = { ...
    'no adaptation', [],        'none',    []; ...
    'single timescale', 0.25,   'single',  single; ...
    'chained MTS', ladder,      'chain',   chain; ...
    'product MTS', ladder,      'product', product};

%% 1. One cell type: elementwise
fprintf('1. One cell type, one route, n = 12 (elementwise, J_code = T J T^T)\n');
for k = 1:size(conditions, 1)
    [name, tau_a, kind, std_cfg] = conditions{k, :};
    args = {'n', 12, 'indegree', 6, 'n_cellTypes', 1, 'cell_type_names', {'all'}, 'f', 1, ...
        'mu_tilde_relative', 0.5, 'sigma_tilde_relative', 1, 'level_of_chaos', 1.3};
    [m, p] = build_network(args, {tau_a}, c, S_c, tau_d, std_route_config({'all'}, std_cfg));
    worst = 0;
    for trial = 1:2
        z = random_state(12, numel(tau_a), kind, 100 + 10 * k + trial);
        J = block_form_jacobian(z, full(p.W), tau_d, S_c, c, tau_a, kind, std_cfg);
        T = lift(p, numel(tau_a), kind);
        J_code = full(SRNNCellTypePairs.compute_Jacobian_fast(T * z, p));
        worst = max(worst, max(abs(T' * J_code * T - J), [], 'all') / max(abs(J), [], 'all'));
    end
    all_passed = report(sprintf('   %-17s N = %3d, largest relative difference %.1e', ...
        name, size(J, 1), worst), worst < tol) && all_passed;
    delete(m);
end

%% 2. Two cell types, four chain routes: the subspace of equal copies
fprintf('2. Two cell types, four chain routes, n = 16 (J_code T = T J)\n');
args = {'n', 16, 'indegree', 8, 'n_cellTypes', 2, 'cell_type_names', {'E', 'I'}, 'f', [0.5 0.5], ...
    'mu_tilde_relative', [3 -3], 'sigma_tilde_relative', [1 1], 'level_of_chaos', 1.3};
[m, p] = build_network(args, {ladder, ladder}, [c c], S_c, tau_d, std_route_config({'E', 'I'}, chain));
worst = 0;
for trial = 1:2
    z = random_state(16, 3, 'chain', 200 + trial);
    J = block_form_jacobian(z, full(p.W), tau_d, S_c, c, ladder, 'chain', chain);
    T = lift(p, 3, 'chain');
    J_code = full(SRNNCellTypePairs.compute_Jacobian_fast(T * z, p));
    worst = max(worst, max(abs(J_code * T - T * J), [], 'all') / max(abs(J), [], 'all'));
end
all_passed = report(sprintf('   %-17s N = %3d (code: %d), largest relative difference %.1e', ...
    'chained MTS', size(J, 1), size(T, 1), worst), worst < tol) && all_passed;

%% 3. The reduction to the dendritic states
% For each eigenvalue s of J, the smallest singular value of
% (1 + s tau_d) I - W G(s) over its largest, which is zero when s solves
% (27). Compared with the same ratio at points that are not eigenvalues.
fprintf('3. Reduction (27): (1 + s tau_d) I - W G(s) is singular at the eigenvalues of J\n');
W = full(p.W);
s_eig = eig(J);
ratio_at = @(s) singular_ratio(s, z, W, tau_d, S_c, c, ladder, chain);
at_eigs = arrayfun(ratio_at, s_eig);
off_eigs = arrayfun(ratio_at, s_eig + 0.05 * (1 + 1i));
all_passed = report(sprintf(['   %d eigenvalues: smallest over largest singular value at most %.1e ' ...
    '(median %.1e); at shifted points at least %.1e'], numel(s_eig), max(at_eigs), median(at_eigs), ...
    min(off_eigs)), max(at_eigs) < 1e-8 && min(off_eigs) > 1e3 * max(at_eigs)) && all_passed;
delete(m);

if all_passed
    fprintf('=== ALL block-form checks PASSED ===\n');
else
    fprintf('=== SOME block-form checks FAILED ===\n');
end

%% ------------------------------------------------------------------------
function J = block_form_jacobian(z, W, tau_d, S_c, c, tau_a, kind, s)
% Equation (24), J = D + U Phi' V^T, on z = [x; a_1 .. a_K; b_1; b_2].
% kind: 'none' | 'single' | 'chain' | 'product'; s holds the constants.
N = size(W, 1);
K = numel(tau_a);
n_b = pool_count(kind);
x = z(1:N);
a = reshape(z(N + (1:N * K)), N, K);
b = reshape(z(N + N * K + (1:N * n_b)), N, n_b);
I = eye(N);
O = zeros(N);

% Operating point (S2) and the slope of the logistic nonlinearity.
h = x - S_c;
if K > 0, h = h - (c / K) * sum(a, 2); end
r = 1 ./ (1 + exp(-4 * h));
Phi = diag(4 * r .* (1 - r));
R = diag(r);

% The synapse: output factor F (theta = r F), its derivatives Fb (one
% cell per pool), the pools' own block A_bb and their rate column U_b.
switch kind
    case 'none'
        F = I; Fb = {}; A_bb = zeros(0); U_b = zeros(0, N);
    case 'single'
        B1 = diag(b(:, 1));
        F = B1; Fb = {I};
        A_bb = -(I / s.tau_rec + R / s.tau_rel);
        U_b = -B1 / s.tau_rel;
    case 'chain'
        B1 = diag(b(:, 1)); B2 = diag(b(:, 2));
        F = B1; Fb = {I, O};
        A_bb = [-(B2 / s.tau_rec(1) + R / s.tau_rel), (I - B1) / s.tau_rec(1); ...
            B2 / (s.n2 * s.tau_rec(1)), -(I / s.tau_rec(2) + (I - B1) / (s.n2 * s.tau_rec(1)))];
        U_b = [-B1 / s.tau_rel; O];
    case 'product'
        B1 = diag(b(:, 1)); B2 = diag(b(:, 2));
        F = B1 * B2; Fb = {B2, B1};
        A_bb = blkdiag(-(I / s.tau_rec(1) + R / s.tau_rel(1)), -(I / s.tau_rec(2) + R / s.tau_rel(2)));
        U_b = [-B1 / s.tau_rel(1); -B2 / s.tau_rel(2)];
end

% D: what each variable does to itself and to the others directly.
T_a = diag(1 ./ tau_a(:));                       % K x K
D_xb = zeros(N, N * n_b);
for mm = 1:n_b
    D_xb(:, (mm - 1) * N + (1:N)) = W * R * Fb{mm} / tau_d;
end
D = [-I / tau_d, zeros(N, N * K), D_xb; ...
    zeros(N * K, N), -kron(T_a, I), zeros(N * K, N * n_b); ...
    zeros(N * n_b, N), zeros(N * n_b, N * K), A_bb];

% U Phi' V^T: everything that passes through the rate.
U = [W * F / tau_d; kron(1 ./ tau_a(:), I); U_b];
V = [I; -(c / max(K, 1)) * kron(ones(K, 1), I); zeros(N * n_b, N)];
J = D + U * Phi * V';
end

function ratio = singular_ratio(s, z, W, tau_d, S_c, c, tau_a, sc)
% Smallest over largest singular value of (1 + s tau_d) I - W G(s), with
% G(s) of (28) for the chain: g_i = b_1i (1 - r_i H_i(s) / tau_rel) *
% phi'_i / (1 + c phi'_i A(s)).
N = size(W, 1);
K = numel(tau_a);
x = z(1:N);
a = reshape(z(N + (1:N * K)), N, K);
b1 = z(N + N * K + (1:N));
b2 = z(2 * N + N * K + (1:N));
r = 1 ./ (1 + exp(-4 * (x - S_c - (c / K) * sum(a, 2))));
slope = 4 * r .* (1 - r);
A_s = mean(1 ./ (1 + s * tau_a));
d11 = b2 / sc.tau_rec(1) + r / sc.tau_rel;
d22 = 1 / sc.tau_rec(2) + (1 - b1) / (sc.n2 * sc.tau_rec(1));
a12 = (1 - b1) / sc.tau_rec(1);
a21 = b2 / (sc.n2 * sc.tau_rec(1));
H = (s + d22) ./ ((s + d11) .* (s + d22) - a12 .* a21);
g = b1 .* (1 - r .* H / sc.tau_rel) .* slope ./ (1 + c * slope * A_s);
sv = svd((1 + s * tau_d) * eye(N) - W * diag(g));
ratio = sv(end) / sv(1);
end

function n_b = pool_count(kind)
switch kind
    case 'none', n_b = 0;
    case 'single', n_b = 1;
    otherwise, n_b = 2;
end
end

function z = random_state(N, K, kind, seed)
% A state inside the variables' ranges: x about the setpoint, a and the
% pools in (0, 1).
stream = RandStream('mt19937ar', 'Seed', seed);
z = [0.2 + 0.5 * randn(stream, N, 1); rand(stream, N * K, 1); ...
    0.05 + 0.9 * rand(stream, N * pool_count(kind), 1)];
end

function sc = std_route_config(names, std_cfg)
% The same depression on every route between the named cell types.
sc = struct();
if isempty(std_cfg), return; end
for pre = 1:numel(names)
    for post = 1:numel(names)
        sc.(names{pre}).(names{post}).std = std_cfg;
    end
end
end

function [m, p] = build_network(args, tau_a, c, S_c, tau_d, synapse_config)
% A homogeneous network: no spreads, no input, no noise.
m = SRNNCellTypePairs(args{:}, 'tau_a', tau_a, 'c', c, 'S_c', S_c, 'tau_d', tau_d, ...
    'activation', 'logistic', 'synapse_config', synapse_config, 'rng_seeds', [3 4], ...
    'lya_method', 'none', 'verbose', 'near-none');
m.build();
p = m.get_params();
end

function T = lift(p, K, kind)
% The map from z = [x; a_1 .. a_K; b_1; b_2] to the code's state: row
% T(code index, z index) = 1. SFA rows are timescale-major within a cell
% type, a route's rows component-major, and every route leaving a cell
% type gets a copy of its neurons' pools.
N = p.n;
n_b = pool_count(kind);
layout = p.state_layout;
T = zeros(layout.N, N * (1 + K + n_b));
T(sub2ind(size(T), layout.x(:), (1:N)')) = 1;
for q = 1:p.n_cellTypes
    idx = p.type_indices{q}(:);
    for k = 1:K
        rows = layout.a{q}((k - 1) * numel(idx) + (1:numel(idx)));
        T(sub2ind(size(T), rows(:), N + (k - 1) * N + idx)) = 1;
    end
    for post = 1:p.n_cellTypes
        for mm = 1:n_b
            rows = layout.b{q, post}((mm - 1) * numel(idx) + (1:numel(idx)));
            T(sub2ind(size(T), rows(:), N + K * N + (mm - 1) * N + idx)) = 1;
        end
    end
end
end

function ok = report(msg, ok)
if ok, tag = 'PASS'; else, tag = 'FAIL'; end
fprintf('%s  [%s]\n', msg, tag);
end
