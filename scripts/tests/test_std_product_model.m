% test_std_product_model.m - the depression model SRNNCellTypePairs has
% today (independent components per route, combined as a PRODUCT), pinned
% down on routes with one component and with two.
%
% Written 2026-11-10, before a second depression model is added beside this
% one (docs/ISSUES.md, FD-004). The product stays the default and must stay
% bit-identical, so this states what "the product model" is in terms that
% do not depend on how the class computes it:
%
%   db_m/dt = (1 - b_m)/tau_rec_m - b_m r / tau_rel_m      (r: the raw rate)
%   theta   = r * prod_m b_m                               (per route)
%   dx_i/dt = (-x_i + sum_j W_ij theta_j^{type(j) -> type(i)} + u_i) / tau_d
%
% The network (make_net) has all four kinds of route at once: E->E with two
% components at DIFFERENT tau_rel/tau_rec ratios, E->I and I->I with one,
% I->E with none, and SFA on E.
%
% Checks:
%   1. Layout: n_b_pairs, the state count, every b at 1 in S0.
%   2. The right-hand side at random interior states equals a per-neuron,
%      loop-written evaluation of the equations above (round-off), with
%      std_zero_floor off and on. This is the check that depression is
%      driven by the raw rate and enters presynaptically, per route.
%   3. From a run: the stored synaptic output is rate times the product of
%      that route's components at every stored sample, exactly; a route
%      with no depression carries the rate itself.
%   4. Constant rate (unconnected neurons, W = 0, constant drive): each
%      component settles at b_m = 1/(1 + r tau_rec_m/tau_rel_m), on a
%      one-component route and on a two-component route whose components
%      have different ratios; theta settles at r times their product. The
%      approach to the steady state follows the closed-form exponential
%      with rate 1/tau_rec_m + r/tau_rel_m.
%   5. The analytic Jacobian against central finite differences, and
%      jacobian_times against the assembled Jacobian, at states from the
%      run and at random interior states, floor off and on.
%   6. The same two Jacobian checks on the manuscript's network preset at
%      n = 40 in each of its three conditions (no depression, one
%      component per route, two per route), at states from a short run.
%   7. Frozen trajectory checksums of the mixed network, noise off (rk4)
%      and on (sra1), captured from the code at commit e0e8e4d. If these
%      move, the default depression model's numbers moved.
%
% test_std_geomean covers the opt-in geometric-mean combine and, at equal
% ratios only, the product's steady state; test_jacobian_times and
% test_SRNNCellTypePairs cover the Jacobian on a three-type network with
% facilitation. This test adds the unequal-ratio steady state, the
% independent right-hand side, and the Jacobian on the paper's conditions.
%
% Prints PASS/FAIL per check and a final banner. About 5 s. Assumes
% setup_paths has run.
%
% See also: SRNNCellTypePairs.dynamics_fast, test_std_geomean,
%           test_jacobian_times, test_pairs_reproducibility

fprintf('=== Testing the product depression model (one and two components) ===\n\n');
all_passed = true;
rng_before = rng;
rng(20261110, 'twister');

%% 1. Layout
m = make_net();
m.build();
p = m.get_params();
p.u_interpolant = m.u_interpolant;
nE = m.n_per_type(1); nI = m.n_per_type(2);
all_passed = check('n_b_pairs is [2 1; 0 1] (pre, post)', isequal(m.n_b_pairs, [2 1; 0 1])) && all_passed;
N_expect = nE * 2 + (nE * 2 + nE * 1 + nI * 1) + m.n;      % a, b, x
all_passed = check(sprintf('N_sys_eqs = %d (a %d + b %d + x %d)', N_expect, nE * 2, nE * 3 + nI, m.n), ...
    m.N_sys_eqs == N_expect && numel(m.S0) == N_expect) && all_passed;
b_rows = [p.state_layout.b{1, 1}, p.state_layout.b{1, 2}, p.state_layout.b{2, 1}, p.state_layout.b{2, 2}];
all_passed = check('every depression variable starts at 1 and the I->E route has none', ...
    all(m.S0(b_rows) == 1) && isempty(p.state_layout.b{2, 1}) && numel(b_rows) == nE * 3 + nI) && all_passed;

%% 2. Right-hand side against the equations, written per neuron
for floor_on = [false true]
    p.std_zero_floor = floor_on;
    e = 0;
    for trial = 1:5
        S = random_state(m.S0, p);
        e = max(e, max(abs(SRNNCellTypePairs.dynamics_fast(0, S, p) - rhs_by_hand(S, p))));
    end
    all_passed = check(sprintf('dynamics_fast equals the per-neuron equations, std_zero_floor = %d (max abs %.1e)', ...
        floor_on, e), e < 1e-12) && all_passed;
end
p.std_zero_floor = false;

%% 3. Synaptic output from a run
m.run();
pd = m.plot_data;
rE = pd.r.E; rI = pd.r.I;
prod_of = @(B) reshape(prod(B, 2), size(B, 1), []);        % neurons x components x time
all_passed = check('stored b has two components on E->E and one on E->I and I->I', ...
    size(pd.b.E.E, 2) == 2 && size(pd.b.E.I, 2) == 1 && size(pd.b.I.I, 2) == 1) && all_passed;
all_passed = check('E->E (two components): synaptic output == r .* b_1 .* b_2 exactly', ...
    isequal(pd.synaptic_output.E.E, rE .* prod_of(pd.b.E.E))) && all_passed;
all_passed = check('E->I (one component): synaptic output == r .* b exactly', ...
    isequal(pd.synaptic_output.E.I, rE .* prod_of(pd.b.E.I))) && all_passed;
all_passed = check('I->I (one component): synaptic output == r .* b exactly', ...
    isequal(pd.synaptic_output.I.I, rI .* prod_of(pd.b.I.I))) && all_passed;
all_passed = check('I->E (no depression): synaptic output == r exactly, no b stored', ...
    isequal(pd.synaptic_output.I.E, rI) && isempty(pd.b.I.E)) && all_passed;
dep = min(pd.b.E.E(:));
all_passed = check(sprintf('depression is engaged in the run (min b on E->E = %.3f < 0.95)', dep), dep < 0.95) && all_passed;
all_passed = check('the two E->E components differ (different time constants)', ...
    max(abs(pd.b.E.E(:, 1, end) - pd.b.E.E(:, 2, end))) > 1e-3) && all_passed;

%% 4. Steady state and relaxation at constant rate
drive = [0.2; 0.4; 0.7];                 % logistic at S_c = 0.4: r about 0.31, 0.5, 0.77
one = struct('tau_rec', 2, 'tau_rel', 0.25);
two = struct('tau_rec', [1.5 6], 'tau_rel', [0.3 0.4]);        % ratios 5 and 15
for cfg = {one, two}
    s = cfg{1};
    out = constant_rate_run(s, drive);
    M = numel(s.tau_rec);
    r = out.r_end;
    b_ss = 1 ./ (1 + r .* (s.tau_rec ./ s.tau_rel));           % neurons x M
    e_b = max(abs(out.b_end - b_ss) ./ b_ss, [], 'all');
    e_th = max(abs(out.theta_end - r .* prod(b_ss, 2)) ./ (r .* prod(b_ss, 2)));
    all_passed = check(sprintf('%d component(s): each b_m -> 1/(1 + r tau_rec_m/tau_rel_m) (rel %.1e)', M, e_b), ...
        e_b < 1e-7) && all_passed;
    all_passed = check(sprintf('%d component(s): theta -> r * prod_m b_m at steady state (rel %.1e)', M, e_th), ...
        e_th < 1e-7) && all_passed;
    all_passed = check(sprintf('%d component(s): the rate is constant at phi(drive) (max dev %.1e)', M, out.r_dev), ...
        out.r_dev < 1e-9) && all_passed;
    % Closed form from b(t0) at constant r: b(t) = b_ss + (b(t0) - b_ss) exp(-(t - t0) k),
    % k = 1/tau_rec + r/tau_rel. t0 = 4 s, after x has settled (40 tau_d).
    k = 1 ./ s.tau_rec + r ./ s.tau_rel;
    e_t = 0;
    for j = 1:numel(out.t_check)
        dt = out.t_check(j) - out.t0;
        b_formula = b_ss + (out.b_t0 - b_ss) .* exp(-dt .* k);
        e_t = max(e_t, max(abs(out.b_check(:, :, j) - b_formula), [], 'all'));
    end
    all_passed = check(sprintf('%d component(s): b(t) follows the closed-form exponential, rate 1/tau_rec + r/tau_rel (max abs %.1e)', M, e_t), ...
        e_t < 1e-6) && all_passed;
end
% Two components at ONE ratio: the product is the square of the one-component factor.
sq = constant_rate_run(struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]), drive);
f1 = 1 ./ (1 + sq.r_end * 8);
e_sq = max(abs(sq.theta_end - sq.r_end .* f1 .^ 2) ./ (sq.r_end .* f1 .^ 2));
all_passed = check(sprintf('two components at one ratio: theta -> r / (1 + r tau_rec/tau_rel)^2 (rel %.1e)', e_sq), ...
    e_sq < 1e-7) && all_passed;

%% 5. Jacobian on the mixed network
for floor_on = [false true]
    p.std_zero_floor = floor_on;
    rows = unique(round([0.25 0.5 0.75 1] * size(m.S_out, 1)));
    states = [num2cell(m.S_out(rows, :)', 1), {random_state(m.S0, p), random_state(m.S0, p)}];
    [e_fd, e_jt] = jacobian_errors(states, p);
    all_passed = check(sprintf('analytic Jacobian vs central differences, std_zero_floor = %d, %d states (rel Frobenius %.1e)', ...
        floor_on, numel(states), e_fd), e_fd < 1e-7) && all_passed;
    all_passed = check(sprintf('jacobian_times == J * Y, std_zero_floor = %d (rel %.1e)', floor_on, e_jt), ...
        e_jt < 1e-12) && all_passed;
end
p.std_zero_floor = false;
% Structure: an x row sees b only through routes INTO its own type.
J = full(SRNNCellTypePairs.compute_Jacobian_fast(m.S_out(end, :)', p));
xE = p.state_layout.x(p.type_indices{1}); xI = p.state_layout.x(p.type_indices{2});
all_passed = check('dx_E/db is zero for E->I and I->I variables, dx_I/db is zero for E->E', ...
    ~any(J(xE, p.state_layout.b{1, 2}), 'all') && ~any(J(xE, p.state_layout.b{2, 2}), 'all') && ...
    ~any(J(xI, p.state_layout.b{1, 1}), 'all')) && all_passed;
all_passed = check('dx_E/db_EE, dx_I/db_EI and dx_I/db_II are not all zero', ...
    any(J(xE, p.state_layout.b{1, 1}), 'all') && any(J(xI, p.state_layout.b{1, 2}), 'all') && ...
    any(J(xI, p.state_layout.b{2, 2}), 'all')) && all_passed;
bEE = p.state_layout.b{1, 1};
all_passed = check('db/db is diagonal (components do not couple to each other) and negative', ...
    isdiag(J(bEE, bEE)) && all(diag(J(bEE, bEE)) < 0)) && all_passed;

%% 6. Jacobian on the manuscript's network, n = 40
P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0p025_dualStd_3cond_mu7revisedMedium';
[~, ~, conditions] = srnn_param_preset(P);
for i = 1:numel(conditions)
    name = conditions{i}.name;
    mp = [];
    evalc(['mp = build_from_preset(P, name, ''n'', 40, ''indegree'', 8, ''F_tracks_network'', true, ' ...
        '''sigma_u_noise'', 0, ''ode_solver'', ''rk4'', ''rng_seeds'', [1 2], ''T_range'', [0 3], ' ...
        '''lya_method'', ''none'', ''store_full_state'', true); mp.run();']);
    pp = mp.get_params();
    pp.u_interpolant = mp.u_interpolant;
    rows = unique(round([0.5 1] * size(mp.S_out, 1)));
    [e_fd, e_jt] = jacobian_errors(num2cell(mp.S_out(rows, :)', 1), pp);
    nb = unique(mp.n_b_pairs(:))';
    all_passed = check(sprintf('%-13s (components per route %s, N = %d): Jacobian vs differences (rel %.1e)', ...
        name, mat2str(nb), mp.N_sys_eqs, e_fd), e_fd < 1e-6) && all_passed;
    all_passed = check(sprintf('%-13s jacobian_times == J * Y (rel %.1e)', name, e_jt), e_jt < 1e-12) && all_passed;
end

%% 7. Frozen checksums
% Captured 2026-11-10 by running this file on commit e0e8e4d (R2026a,
% Windows), the code as it stood before FD-004. sum / norm of the final
% state and of the whole stored trajectory. Compared to a relative 1e-9:
% the run is 3 s on 30 neurons, far too short for round-off to grow that
% far, and a change of model moves these in the first few figures.
% DO NOT update these to make the test pass. If the default depression
% model was meant to change, say so in the commit that changes them.
frozen = {
    'rk4, noise off',   {'sigma_u_noise', 0,    'ode_solver', 'rk4'},  [40.2214538463209 7.03792608384453 25701.9724519229 190.012027429607]
    'sra1, noise 0.05', {'sigma_u_noise', 0.05, 'ode_solver', 'sra1'}, [39.8835459838761 7.06427860469557 25455.6415791316 191.364880541003]
    };
for i = 1:size(frozen, 1)
    mf = make_net(frozen{i, 2}{:});
    mf.build(); mf.run();
    got = [sum(mf.S_out(end, :)), norm(mf.S_out(end, :)), sum(mf.S_out(:)), norm(mf.S_out(:))];
    want = frozen{i, 3};
    ok = all(abs(got - want) <= 1e-9 * max(1, abs(want)));
    all_passed = check(sprintf('frozen checksums, %s', frozen{i, 1}), ok) && all_passed;
    if ~ok
        fprintf('      got  %s\n      want %s\n', mat2str(got, 15), mat2str(want, 15));
    end
end

rng(rng_before);
fprintf('\n');
if all_passed
    fprintf('=== ALL product depression model TESTS PASSED ===\n');
else
    fprintf('=== SOME product depression model TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function m = make_net(varargin)
% Two types, every kind of route: two components (unequal ratios), one, none.
sc = struct();
sc.E.E.std = struct('tau_rec', [1.5 6], 'tau_rel', [0.3 0.4]);
sc.E.I.std = struct('tau_rec', 2, 'tau_rel', 0.25);
sc.I.I.std = struct('tau_rec', 3, 'tau_rel', 1);
m = SRNNCellTypePairs( ...
    'n', 30, 'indegree', 10, 'n_cellTypes', 2, 'cell_type_names', {'E', 'I'}, ...
    'f', [0.6 0.4], 'mu_tilde_relative', [4 -4; 4 -4], 'sigma_tilde_relative', [1 1; 1 1], ...
    'level_of_chaos', 1.5, 'tau_a', {[0.25 2], zeros(1, 0)}, 'c', [0.3 0], ...
    'synapse_config', sc, 'rng_seeds', [11 12], 'T_range', [0 3], 'fs', 200, ...
    'ode_solver', 'rk4', 'lya_method', 'none', 'store_full_state', true, varargin{:});
m.input_config.intrinsic_drive = 0.3;
end

function dS = rhs_by_hand(S, p)
% The model equations, one neuron and one synapse at a time. Deliberately
% not vectorised and not sharing a line with dynamics_fast beyond the
% state layout and the nonlinearity.
L = p.state_layout;
n = p.n;
C = p.n_cellTypes;
type_of = zeros(n, 1);
for q = 1:C; type_of(p.type_indices{q}) = q; end
u = p.u_interpolant(0)';
x = S(L.x);
dS = zeros(size(S));

% Rates: r_i = phi(x_i - (c/K) sum_k a_ik)
x_eff = x;
for q = 1:C
    idx = p.type_indices{q};
    K = p.n_a(q);
    for jj = 1:numel(idx)
        for k = 1:K
            x_eff(idx(jj)) = x_eff(idx(jj)) - p.c_eff(q) * S(L.a{q}((k - 1) * numel(idx) + jj));
        end
    end
end
r = p.activation_function(x_eff);

% SFA: da_ik/dt = (r_i - a_ik) / tau_a_ik
for q = 1:C
    idx = p.type_indices{q};
    for jj = 1:numel(idx)
        for k = 1:p.n_a(q)
            row = L.a{q}((k - 1) * numel(idx) + jj);
            dS(row) = (r(idx(jj)) - S(row)) / p.tau_a_matrix{q}(jj, k);
        end
    end
end

% Depression: db/dt = (1 - b)/tau_rec - b r/tau_rel, driven by the raw rate
% of the presynaptic neuron; theta(j, post) = r_j * prod_m b_m.
theta = zeros(n, C);
for pre = 1:C
    idx = p.type_indices{pre};
    for post = 1:C
        M = p.n_b_pairs(pre, post);
        for jj = 1:numel(idx)
            j = idx(jj);
            factor = 1;
            for mm = 1:M
                row = L.b{pre, post}((mm - 1) * numel(idx) + jj);
                b = S(row);
                dS(row) = (1 - b) / p.tau_b_rec{pre, post}(mm) - b * r(j) / p.tau_b_rel{pre, post}(mm);
                factor = factor * b;
            end
            if M > 0 && p.std_zero_floor
                % The floor maps the factor at r = 1 to zero.
                p_min = prod(p.tau_b_rel{pre, post} ./ (p.tau_b_rec{pre, post} + p.tau_b_rel{pre, post}));
                factor = (factor - p_min) / (1 - p_min);
            end
            theta(j, post) = r(j) * factor;
        end
    end
end

% Dendrite: dx_i/dt = (-x_i + sum_j W_ij theta_j^{type(j) -> type(i)} + u_i) / tau_d
W = full(p.W);
for i = 1:n
    acc = 0;
    for j = 1:n
        if W(i, j) ~= 0
            acc = acc + W(i, j) * theta(j, type_of(i));
        end
    end
    dS(L.x(i)) = (-x(i) + acc + u(i)) / p.tau_d;
end
end

function S = random_state(S0, p)
S = S0;
for q = 1:p.n_cellTypes
    S(p.state_layout.a{q}) = 0.6 * rand(numel(p.state_layout.a{q}), 1);
end
for pre = 1:p.n_cellTypes
    for post = 1:p.n_cellTypes
        rows = p.state_layout.b{pre, post};
        S(rows) = 0.2 + 0.8 * rand(numel(rows), 1);
    end
end
S(p.state_layout.x) = 0.4 + 0.5 * randn(p.n, 1);
end

function [e_fd, e_jt] = jacobian_errors(states, p)
e_fd = 0; e_jt = 0;
for k = 1:numel(states)
    S = states{k};
    J = SRNNCellTypePairs.compute_Jacobian_fast(S, p);
    Jfd = SRNNCellTypePairs.finite_difference_jacobian(S, p, 1e-6);
    e_fd = max(e_fd, norm(full(J) - Jfd, 'fro') / norm(Jfd, 'fro'));
    Y = randn(numel(S), 6);
    e_jt = max(e_jt, norm(SRNNCellTypePairs.jacobian_times(S, Y, p) - J * Y, 'fro') / norm(J * Y, 'fro'));
end
end

function out = constant_rate_run(std_config, drive)
% Unconnected neurons (W = 0, no SFA, no noise) at a constant drive: x
% settles to the drive within a few tau_d, after which the rate is constant
% and each b obeys a linear ODE.
sc = struct(); sc.E.E.std = std_config;
n = numel(drive);
m = SRNNCellTypePairs('n', n, 'indegree', 1, 'n_cellTypes', 1, 'cell_type_names', {'E'}, ...
    'f', 1, 'mu_tilde_relative', 0, 'sigma_tilde_relative', 0, ...
    'tau_a', {zeros(1, 0)}, 'c', 0, 'x0_std', 0, 'synapse_config', sc, ...
    'T_range', [0 90], 'fs', 100, 'ode_solver', 'ode45', 'lya_method', 'none', ...
    'store_full_state', true);
m.input_config.amp = 0;
m.input_config.intrinsic_drive = drive;
m.build(); m.run();
p = m.get_params();
M = numel(std_config.tau_rec);
t = m.t_out;
b_at = @(row) reshape(m.S_out(row, p.state_layout.b{1, 1}), n, M);
r_all = p.activation_function(m.S_out(:, p.state_layout.x)')';     % time x n
[~, i0] = min(abs(t - 4));                 % 40 tau_d: x equals the drive to round-off
out.t0 = t(i0);
out.b_t0 = b_at(i0);
out.t_check = [4.5 5 6 8 15];
out.b_check = zeros(n, M, numel(out.t_check));
for j = 1:numel(out.t_check)
    [~, ij] = min(abs(t - out.t_check(j)));
    out.t_check(j) = t(ij);
    out.b_check(:, :, j) = b_at(ij);
end
out.r_end = r_all(end, :)';
out.r_dev = max(abs(r_all(i0:end, :) - r_all(end, :)), [], 'all');
out.b_end = b_at(size(m.S_out, 1));
out.theta_end = m.plot_data.synaptic_output.E.E(:, end);
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
