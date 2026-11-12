% test_std_chain_model.m - chained vesicle-pool depression in a network.
%
% SRNNCellTypePairs with routes that ask for std.model = 'chain' (FD-004,
% FD-005; docs/Plan_A_chained_vessicles.md, sections 1, 4 and 5). A chain
% route carries a releasable pool b1 and a reserve b2 per presynaptic
% neuron, with rate r:
%
%   db1/dt = b2 (1 - b1)/tau_rec1 - b1 r/tau_rel                      (C1)
%   db2/dt = (1 - b2)/tau_rec2 - (1/n2) b2 (1 - b1)/tau_rec1          (C2)
%   theta  = r b1                                                     (C3)
%
% and at a constant rate, with k1 = r tau_rec1/tau_rel,
% k2 = r tau_rec2/(n2 tau_rel), s = 1 + k1 + k2:
%
%   b1 = (s - sqrt(s^2 - 4 k2))/(2 k2),  b2 = 1 - k2 b1               (C4)
%
% The route model itself (route_pool) is checked alone in
% test_route_models. This file checks it IN THE CLASS: that the right
% numbers reach the right rows of the state vector and the right entries of
% the Jacobian, beside routes of the other model.
%
% THE NETWORK: 24 neurons, two cell types, SFA on both, and four routes of
% four kinds: E->E a chain, E->I one product component, I->E a chain with
% other parameters and a route scale, I->I undepressed.
%
% Checks (the numbers 1 to 9 are Plan A section 5's):
%   1. dynamics_fast equals a per-neuron loop over the equations, written
%      here from (C1), (C2) and the network equations, at random states.
%   2. Steady state: at the pools' (C4) values for the neurons' rates,
%      dynamics_fast's pool rows vanish; and unconnected neurons held at
%      six constant rates arrive at (C4) by integration.
%   3. The stored synaptic output is the rate times b1, exactly.
%   4. The analytic Jacobian equals finite_difference_jacobian at random
%      interior states and at states from a trajectory. The chain's 2 x 2
%      block is not symmetric, and this is the check that fails if the
%      class writes its entry (m, k) at (k, m): see "4b" below.
%   5. jacobian_times equals the assembled Jacobian times Y for K = 1,
%      several, and Y = eye(N). The chain routes take the full-A branch of
%      jacobian_times, and applying only the diagonal gives another answer.
%   6. n2 = 1e12 reproduces the single pool: b1 follows a one-component
%      product route's b to 1e-9 over a run, and b2 stays within 1e-10 of 1.
%   7. b1 and b2 stay in [0, 1] over a run whose rates swing from below
%      0.01 to above 0.99.
%   8. Each refusal raises its own error, at build: the zero floor,
%      combine, one or three tau_rec, a vector tau_rel, no n2, an invalid
%      n2, output-only stf on a chain route (the error names
%      stf.model = 'chain'; test_stf_chain_model has the release-coupled
%      form), an unknown std.model. Switching the floor on in params after
%      get_params raises the same error as at build.
%   9. Benettin, QR and top-K run on a chain network and on the same
%      network with two-component product routes; on each the three
%      leading exponents lie within 0.2 of each other, the bound for
%      finite-time exponents of a stable network. And Liouville: the sum
%      of all the chain network's exponents equals the mean trace of its
%      Jacobian over the window. And the TIGHT comparison, on a chain
%      network that is clearly chaotic (no small one is): the regression
%      network, n = 100, four chain routes, gain 2, where top-K and
%      Benettin agree to 1 % on lambda_1 = +1.27.
% And, from what the refactor behind the route interface found:
%  10. params: n_b_pairs, std_chain, std_n2, the nominal rows; both pools
%      start at 1.
%  11. routes_identical is true for identical chains, false between a chain
%      and a two-component product with the same time constants, and false
%      between chains that differ in n2 alone.
%  12. std_route_factor is b1 on a chain route (2-D and 3-D b), the product
%      on a product route, [] on an undepressed one; combine_std_pair,
%      apply_std_zero_floor_pair and std_zero_floor_pmin raise
%      ProductRouteOnly on a chain route and still serve a product route.
%  13. plot and plot_celltypes draw b1 for a chain route, not b1 b2, and
%      for a product route exactly prod(b), as before.
%  14. A chain model reproduces bit for bit: two fresh objects, run() twice
%      on one, and an object saved, loaded, built and run again.
%
% 4b, DONE ONCE BY HAND AND NOT KEPT: with entry (m, k) of A written at
% (k, m) in compute_Jacobian_fast (repmat and repelem exchanged), check 4
% fails with an error of the size of the Jacobian's entries and check 5
% fails with it. The commit message of this file records the numbers.
%
% Prints PASS/FAIL per check and a final banner. About 15 to 20 s, half of
% it the two QR runs of check 9. Assumes
% setup_paths has run.
%
% See also: route_pool, test_route_models, test_jacobian_times,
%           test_SRNNCellTypePairs, test_pairs_reproducibility

fprintf('=== Testing chained vesicle-pool depression in SRNNCellTypePairs ===\n\n');
all_passed = true;
rng_at_entry = rng;
restore_rng = onCleanup(@() rng(rng_at_entry));
old_visibility = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
restore_visibility = onCleanup(@() set(groot, 'DefaultFigureVisible', old_visibility));

sc = four_kinds();
m = small_network(sc);
m.build();
m.run();
p = m.get_params();
p.u_interpolant = m.u_interpolant;
lay = p.state_layout;
N = lay.N;
nE = p.n_per_type(1);
nI = p.n_per_type(2);
is_pool = false(N, 1);
is_pool([lay.b{1, 1}, lay.b{2, 1}]) = true;        % the two chain routes' rows

%% 10. params and the starting state
ok = isequal(p.n_b_pairs, [2 1; 2 0]) && isequal(p.std_chain, [true false; true false]) && ...
    isequaln(p.std_n2, [2 NaN; 1.5 NaN]) && isequal(p.tau_b_rec{1, 1}, [2 5]) && ...
    isequal(p.tau_b_rel{1, 1}, 0.25) && isequal(p.tau_b_rec{2, 1}, [1 3]) && ...
    isequal(p.route_scale, [1 1; 1.3 1]) && strcmp(p.route{1, 1}.P.model, 'pool') && ...
    strcmp(p.route{1, 2}.P.model, 'product') && ~p.route{1, 1}.A_is_diagonal && ...
    numel(lay.b{1, 1}) == 2 * nE && numel(lay.b{2, 1}) == 2 * nI && numel(lay.b{1, 2}) == nE && ...
    isempty(lay.b{2, 2}) && isempty([lay.g{:}]);
all_passed = check('10. params: n_b_pairs, std_chain, std_n2, nominal rows, the state layout', ok) && all_passed;
all_passed = check('10. both pools of every chain route start at 1', ...
    all(m.S0(is_pool) == 1) && all(m.S0(lay.b{1, 2}) == 1)) && all_passed;

%% 1. dynamics_fast against the per-neuron loop
rng(101, 'twister');
e_rhs = 0;
for trial = 1:5
    S = random_state(p, N);
    t = 0.3 + 0.7 * trial;
    e_rhs = max(e_rhs, rel_err(SRNNCellTypePairs.dynamics_fast(t, S, p), reference_rhs(t, S, m, p, sc)));
end
all_passed = check(sprintf('1. dynamics_fast equals the per-neuron loop of (C1), (C2) and the network (rel %.1e)', e_rhs), ...
    e_rhs <= 1e-13) && all_passed;

%% 2. steady state
% a. In the network: put every chain pool at (C4) for its neuron's rate.
e_ss = 0;
for trial = 1:4
    S = random_state(p, N);
    r = reference_rates(S, m, p);
    for pre = 1:2
        idx = p.type_indices{pre};
        cfg = sc.(m.cell_type_names{pre}).E.std;
        [b1, b2] = chain_steady_state(r(idx), cfg);
        S(lay.b{pre, 1}) = [b1; b2];
    end
    dS = SRNNCellTypePairs.dynamics_fast(1, S, p);
    e_ss = max(e_ss, max(abs(dS(is_pool))));
end
all_passed = check(sprintf('2. at (C4) for each neuron''s rate the pools'' derivatives vanish (max %.1e)', e_ss), ...
    e_ss <= 1e-14) && all_passed;

% b. By integration: six unconnected neurons, each held at a constant rate.
drive = [-0.2; 0.1; 0.3; 0.4; 0.6; 1.2];
u = SRNNCellTypePairs('n', 6, 'indegree', 1, 'n_cellTypes', 1, 'cell_type_names', {'E'}, 'f', 1, ...
    'mu_tilde_relative', 0, 'sigma_tilde_relative', 0, 'tau_a', {[]}, 'c', 0, ...
    'synapse_config', struct('E', struct('E', struct('std', sc.E.E.std))), ...
    'T_range', [0 110], 'fs', 20, 'ode_solver', 'rk4', 'lya_method', 'none', ...
    'store_full_state', true, 'verbose', 'near-none');
u.input_config.amp = 0;
u.input_config.intrinsic_drive = drive;
u.build();
u.run();
r_held = u.plot_data.r.E(:, end);
[b1, b2] = chain_steady_state(r_held, sc.E.E.std);
b_end = u.plot_data.b.E.E(:, :, end);
e_int = max(abs(b_end - [b1, b2]), [], 'all');
all_passed = check(sprintf(['2. unconnected neurons at constant rates %.3f to %.3f arrive at (C4) ' ...
    '(max %.1e)'], min(r_held), max(r_held), e_int), ...
    nnz(u.W) == 0 && e_int <= 1e-8 && max(r_held) - min(r_held) > 0.8) && all_passed;

%% 3. stored synaptic output
pd = m.plot_data;
ok = isequal(pd.synaptic_output.E.E, pd.r.E .* reshape(pd.b.E.E(:, 1, :), nE, [])) && ...
    isequal(pd.synaptic_output.I.E, pd.r.I .* reshape(pd.b.I.E(:, 1, :), nI, [])) && ...
    isequal(pd.synaptic_output.E.I, pd.r.E .* reshape(pd.b.E.I, nE, [])) && ...
    isequal(pd.synaptic_output.I.I, pd.r.I) && ...
    ~isequal(pd.synaptic_output.E.E, pd.r.E .* reshape(prod(pd.b.E.E, 2), nE, []));
all_passed = check('3. stored synaptic output is the rate times b1, exactly (and not times b1 b2)', ok) && all_passed;

%% 4. the analytic Jacobian against finite differences
rng(102, 'twister');
states = cell(1, 6);
for k = 1:3
    states{k} = random_state(p, N);
end
at = round(size(m.S_out, 1) * [0.3 0.6 1]);
for k = 1:3
    states{3 + k} = m.S_out(at(k), :)';
end
e_fd = zeros(1, 6);
e_fd_pool = zeros(1, 6);
asym = 0;
for k = 1:6
    J = full(SRNNCellTypePairs.compute_Jacobian_fast(states{k}, p));
    J_fd = SRNNCellTypePairs.finite_difference_jacobian(states{k}, p, 1e-6);
    e_fd(k) = max(abs(J - J_fd), [], 'all') / max(abs(J_fd), [], 'all');
    e_fd_pool(k) = max(abs(J(is_pool, is_pool) - J_fd(is_pool, is_pool)), [], 'all');
    B = J(is_pool, is_pool);
    asym = max(asym, max(abs(B - B'), [], 'all'));
end
all_passed = check(sprintf('4. analytic Jacobian equals finite differences at 3 random states (rel %.1e)', max(e_fd(1:3))), ...
    max(e_fd(1:3)) <= 1e-7) && all_passed;
all_passed = check(sprintf('4. and at 3 states from a trajectory (rel %.1e)', max(e_fd(4:6))), ...
    max(e_fd(4:6)) <= 1e-7) && all_passed;
all_passed = check(sprintf(['4. the pool block is not symmetric (|B - B''| up to %.2f), so a transposed 2 x 2 ' ...
    'cannot pass: block error %.1e'], asym, max(e_fd_pool)), asym > 0.05 && max(e_fd_pool) <= 1e-7) && all_passed;

%% 5. jacobian_times against the assembled Jacobian
rng(103, 'twister');
e_jt = zeros(1, 3);
e_diag_only = Inf;
p_diag = p;
for k = 1:numel(p_diag.route)
    p_diag.route{k}.A_is_diagonal = true;       % a model that understated its coupling
end
Ys = {randn(N, 1), randn(N, 7), eye(N)};
for k = 1:6
    S = states{k};
    J = SRNNCellTypePairs.compute_Jacobian_fast(S, p);
    for j = 1:3
        e_jt(j) = max(e_jt(j), rel_err(SRNNCellTypePairs.jacobian_times(S, Ys{j}, p), J * Ys{j}));
    end
    e_diag_only = min(e_diag_only, rel_err(SRNNCellTypePairs.jacobian_times(S, Ys{3}, p_diag), full(J)));
end
all_passed = check(sprintf('5. jacobian_times equals J * Y for K = 1 (rel %.1e), K = 7 (%.1e), Y = eye(N) (%.1e)', ...
    e_jt(1), e_jt(2), e_jt(3)), max(e_jt) <= 1e-13) && all_passed;
all_passed = check(sprintf(['5. the chain routes take the full-A branch: applying only the diagonal ' ...
    'is wrong by %.1e at least'], e_diag_only), ...
    ~p.route{1, 1}.A_is_diagonal && ~p.route{2, 1}.A_is_diagonal && p.route{1, 2}.A_is_diagonal && ...
    e_diag_only > 1e-3) && all_passed;

%% 6. a very large reserve is the single pool
big = sc;
big.E.E.std.n2 = 1e12;
big.I.E.std.n2 = 1e12;
one = sc;
one.E.E.std = struct('tau_rec', sc.E.E.std.tau_rec(1), 'tau_rel', sc.E.E.std.tau_rel);
one.I.E.std = struct('tau_rec', sc.I.E.std.tau_rec(1), 'tau_rel', sc.I.E.std.tau_rel);
m_big = small_network(big);
m_big.build();
m_big.run();
m_one = small_network(one);
m_one.build();
m_one.run();
bb = m_big.plot_data.b;
bo = m_one.plot_data.b;
e_b1 = max([max(abs(bb.E.E(:, 1, :) - bo.E.E), [], 'all'), max(abs(bb.I.E(:, 1, :) - bo.I.E), [], 'all')]);
e_b2 = max([max(abs(bb.E.E(:, 2, :) - 1), [], 'all'), max(abs(bb.I.E(:, 2, :) - 1), [], 'all')]);
e_x = max(abs(m_big.plot_data.x.E - m_one.plot_data.x.E), [], 'all');
depth = 1 - min(bo.E.E, [], 'all');
all_passed = check(sprintf(['6. n2 = 1e12: b1 follows the one-component route (max %.1e, x %.1e) while b ' ...
    'falls by %.2f; b2 stays at 1 (max %.1e)'], e_b1, e_x, depth, e_b2), ...
    isequal(m_big.W, m_one.W) && e_b1 <= 1e-9 && e_x <= 1e-9 && e_b2 <= 1e-10 && depth > 0.1) && all_passed;

%% 7. the pools stay in [0, 1] under large rate swings
m_swing = small_network(sc, 'T_range', [0 12], 'rng_seeds', [21 22]);
m_swing.input_config.amp = 4;
m_swing.input_config.n_steps = 8;
m_swing.input_config.no_stim_pattern = logical([1 0 0 1 0 0 1 0]);
m_swing.build();
m_swing.run();
ps = m_swing.plot_data;
pools = [ps.b.E.E(:); ps.b.I.E(:)];
rate_all = [ps.r.E(:); ps.r.I(:)];
all_passed = check(sprintf(['7. rates from %.3f to %.3f: b1 in [%.4f, %.4f], b2 in [%.4f, %.4f], ' ...
    'all inside [0, 1]'], min(rate_all), max(rate_all), ...
    min(ps.b.E.E(:, 1, :), [], 'all'), max(ps.b.E.E(:, 1, :), [], 'all'), ...
    min(ps.b.E.E(:, 2, :), [], 'all'), max(ps.b.E.E(:, 2, :), [], 'all')), ...
    all(pools >= 0 & pools <= 1) && min(rate_all) < 0.01 && max(rate_all) > 0.99 && ...
    min(ps.b.E.E(:, 1, :), [], 'all') < 0.2 && min(ps.b.E.E(:, 2, :), [], 'all') < 0.8) && all_passed;

%% 8. the refusals
chain = sc.E.E.std;
with = @(name, value) setfield(chain, name, value);
on_EE = @(std_cfg) setfield(sc, 'E', setfield(sc.E, 'E', struct('std', std_cfg)));
with_stf = sc;
with_stf.E.E.stf = struct('tau_dec', 1.5, 'tau_fac', 0.25, 'G', 3);
refusals = { ...
    'std_zero_floor on a network with a chain route', @() small_network(sc, 'std_zero_floor', true), 'ChainZeroFloor'
    'combine on a chain route',  @() small_network(on_EE(with('combine', 'geomean'))), 'ChainCombine'
    'one tau_rec',               @() small_network(on_EE(with('tau_rec', 2))), 'ChainTauRecCount'
    'three tau_rec',             @() small_network(on_EE(with('tau_rec', [2 5 9]))), 'ChainTauRecCount'
    'a vector tau_rel',          @() small_network(on_EE(with('tau_rel', [0.25 0.4]))), 'ChainTauRelNotScalar'
    'no n2',                     @() small_network(on_EE(rmfield(chain, 'n2'))), 'ChainMissingN2'
    'n2 not positive',           @() small_network(on_EE(with('n2', -1))), 'ChainInvalidN2'
    'output-only stf on a chain route', @() small_network(with_stf), 'ChainStfOutputOnly'
    'an unknown std.model',      @() small_network(on_EE(with('model', 'chian'))), 'UnknownStdModel'};
ok = true;
ids = cell(1, size(refusals, 1));
for k = 1:size(refusals, 1)
    make = refusals{k, 2};
    ids{k} = refusals{k, 3};
    ok_k = throws_id(@() build(make), ['SRNNCellTypePairs:' ids{k}]);
    if ~ok_k, fprintf('      %s did not raise SRNNCellTypePairs:%s at build\n', refusals{k, 1}, ids{k}); end
    ok = ok && ok_k;
end
all_passed = check(sprintf('8. %d refusals each raise their own identifier at build (%d distinct)', ...
    size(refusals, 1), numel(unique(ids))), ok && numel(unique(ids)) == size(refusals, 1) - 1) && all_passed;
p_floor = p;
p_floor.std_zero_floor = true;
S = states{4};
ok = throws_id(@() SRNNCellTypePairs.dynamics_fast(1, S, p_floor), 'SRNNCellTypePairs:ChainZeroFloor') && ...
    throws_id(@() SRNNCellTypePairs.compute_Jacobian_fast(S, p_floor), 'SRNNCellTypePairs:ChainZeroFloor') && ...
    throws_id(@() SRNNCellTypePairs.jacobian_times(S, Ys{1}, p_floor), 'SRNNCellTypePairs:ChainZeroFloor') && ...
    throws_id(@() SRNNCellTypePairs.std_route_factor(pd.b.E.E, p_floor, 1, 1), 'SRNNCellTypePairs:ChainZeroFloor');
all_passed = check('8. the floor switched on in params after get_params is refused too, not applied', ok) && all_passed;
explicit = on_EE(struct('model', 'product', 'tau_rec', [2 5], 'tau_rel', 0.25));
implicit = on_EE(struct('tau_rec', [2 5], 'tau_rel', 0.25));
pe = build(@() small_network(explicit)).get_params();
pim = build(@() small_network(implicit)).get_params();
all_passed = check('8. std.model = ''product'' is the same route as no model field', ...
    isequal(pe.route{1, 1}.P, pim.route{1, 1}.P) && ~pe.std_chain(1, 1) && isnan(pe.std_n2(1, 1))) && all_passed;

%% 11. routes_identical
same_chain = struct('std', sc.E.E.std);
all_chain = struct('E', struct('E', same_chain, 'I', same_chain), 'I', struct('E', same_chain, 'I', same_chain));
as_product = struct('std', struct('tau_rec', [2 5], 'tau_rel', 0.25));
mixed = all_chain;
mixed.I.I = as_product;
other_n2 = all_chain;
other_n2.I.I.std.n2 = 3;
all_product = struct('E', struct('E', as_product, 'I', as_product), 'I', struct('E', as_product, 'I', as_product));
identical = @(cfg) SRNNCellTypePairs.routes_identical(build(@() small_network(cfg)).get_params());
all_passed = check(['11. routes_identical: true for four identical chains and for four identical products; ' ...
    'false for a chain beside a product with the same time constants, and for chains differing in n2'], ...
    identical(all_chain) && identical(all_product) && ~identical(mixed) && ~identical(other_n2) && ...
    ~SRNNCellTypePairs.routes_identical(p)) && all_passed;
fields_chain = SRNNCellTypePairs.route_fields(p, 1, 1);
fields_product = SRNNCellTypePairs.route_fields(p, 1, 2);
all_passed = check('11. route_fields carries the route model''s signature', ...
    strcmp(fields_chain.signature.model, 'pool') && fields_chain.signature.n2 == 2 && ...
    strcmp(fields_product.signature.model, 'product') && fields_product.n_b == 1) && all_passed;

%% 12. std_route_factor, and the product-only delegations
b3 = pd.b.E.E;                                  % nE x 2 x nt
D3 = SRNNCellTypePairs.std_route_factor(b3, p, 1, 1);
D2 = SRNNCellTypePairs.std_route_factor(b3(:, :, end), p, 1, 1);
Dp = SRNNCellTypePairs.std_route_factor(pd.b.E.I, p, 1, 2);
ok = isequal(D3, b3(:, 1, :)) && isequal(size(D3), size(prod(b3, 2))) && isequal(D2, b3(:, 1, end)) && ...
    isequal(Dp, prod(pd.b.E.I, 2)) && isempty(SRNNCellTypePairs.std_route_factor(pd.b.I.I, p, 2, 2)) && ...
    throws_id(@() SRNNCellTypePairs.std_route_factor(b3(:, 1, :), p, 1, 1), 'SRNNCellTypePairs:RouteFactorShape');
all_passed = check('12. std_route_factor: b1 on a chain route (2-D and 3-D), the product on a product route, [] on none', ok) && all_passed;
id = 'SRNNCellTypePairs:ProductRouteOnly';
named = false;
try
    SRNNCellTypePairs.combine_std_pair(0.5, p, 1, 1);
catch err
    named = contains(err.message, 'std_route_factor');
end
ok = throws_id(@() SRNNCellTypePairs.combine_std_pair(0.5, p, 1, 1), id) && ...
    throws_id(@() SRNNCellTypePairs.apply_std_zero_floor_pair(0.5, p, 2, 1), id) && ...
    throws_id(@() SRNNCellTypePairs.std_zero_floor_pmin(p, 1, 1), id) && named && ...
    SRNNCellTypePairs.combine_std_pair(0.5, p, 1, 2) == 0.5 && ...
    SRNNCellTypePairs.apply_std_zero_floor_pair(0.5, p, 1, 2) == 0.5 && ...
    isequal(size(SRNNCellTypePairs.std_zero_floor_pmin(p, 1, 2)), [nE 1]);
all_passed = check('12. the three product-only delegations raise ProductRouteOnly on a chain route, naming std_route_factor', ok) && all_passed;

%% 13. the class's plots draw the factor
[fig, ax] = m.plot_celltypes();
std_row = 6;                                    % u, x, r, synaptic output, SFA, STD
drawn_EE = line_rows(ax(std_row, 1), 'PostType_E');
drawn_EI = line_rows(ax(std_row, 1), 'PostType_I');
drawn_IE = line_rows(ax(std_row, 2), 'PostType_E');
b1_EE = reshape(pd.b.E.E(:, 1, :), nE, []);
ok = isequal(drawn_EE, sortrows(b1_EE)) && isequal(drawn_IE, sortrows(reshape(pd.b.I.E(:, 1, :), nI, []))) && ...
    isequal(drawn_EI, sortrows(reshape(pd.b.E.I, nE, []))) && ...
    ~isequal(drawn_EE, sortrows(reshape(prod(pd.b.E.E, 2), nE, []))) && ...
    strcmp(ax(std_row, 1).YLabel.String, 'STD factor');
close(fig);
all_passed = check('13. plot_celltypes draws b1 for the chain routes (not b1 b2) and b for the product route', ok) && all_passed;
[fig, ax1] = m.plot();
mean_lines = findall(ax1(6), 'Type', 'line');
drawn = sortrows(vertcat(mean_lines.YData));
expected = sortrows([mean(b1_EE, 1); mean(reshape(pd.b.E.I, nE, []), 1); mean(reshape(pd.b.I.E(:, 1, :), nI, []), 1)]);
ok = isequal(drawn, expected);
close(fig);
all_passed = check('13. plot draws each route''s mean factor', ok) && all_passed;
m_prod = build(@() small_network(all_product));
m_prod.run();
[fig, axp] = m_prod.plot_celltypes();
drawn = line_rows(axp(std_row, 1), 'PostType_E');
ok = isequal(drawn, sortrows(reshape(prod(m_prod.plot_data.b.E.E, 2), nE, []))) && ...
    strcmp(axp(std_row, 1).YLabel.String, 'STD (product)');
close(fig);
all_passed = check('13. a product network draws prod(b) under the label it had', ok) && all_passed;

%% 14. reproducibility
A = small_network(sc, 'lya_method', 'topk', 'lya_K', 5, 'lya_K_auto', false, 'lya_dt', 0.1, ...
    'lya_T_interval', [2 4], 'lya_warmup', 1);
A.build();
A.run();
rng('shuffle');
B = small_network(sc, 'lya_method', 'topk', 'lya_K', 5, 'lya_K_auto', false, 'lya_dt', 0.1, ...
    'lya_T_interval', [2 4], 'lya_warmup', 1);
B.build();
B.run();
[same, what] = same_run(A, B);
all_passed = check(sprintf('14. two fresh chain models are identical%s (lambda_1 = %+.3f)', what, A.lya_results.LLE), ...
    same && isfinite(A.lya_results.LLE)) && all_passed;
S_first = A.S_out;
lya_first = untimed(A.lya_results);
A.run();
all_passed = check('14. run() twice on one object reproduces itself', ...
    isequal(A.S_out, S_first) && isequaln(untimed(A.lya_results), lya_first)) && all_passed;
file = [tempname '.mat'];
remove_file = onCleanup(@() delete_if_there(file));
save(file, 'A');
loaded = load(file);
L = loaded.A;
ok = isequal(L.S_out, S_first) && isequal(L.synapse_config, sc);
L.build();
L.run();
[same, what] = same_run(L, B);
all_passed = check(sprintf('14. saved, loaded, built and run again: identical%s', what), ok && same) && all_passed;

%% 9. the three Lyapunov estimators
% On 8 neurons and a 2.5 s window, because the ode45 QR costs about a
% second of wall time per simulated second even at this size.
lya = {'n', 8, 'indegree', 3, 'T_range', [0 5], 'lya_T_interval', [2.5 5], 'lya_warmup', 2, 'lya_dt', 0.25};
methods = {'benettin', 'qr', 'topk'};
lle = zeros(2, 3);
configs = {sc, product_twin(sc)};
for c = 1:2
    for k = 1:3
        extra = {};
        if strcmp(methods{k}, 'topk'), extra = {'lya_K', 5, 'lya_K_auto', false}; end
        mk = small_network(configs{c}, lya{:}, 'lya_method', methods{k}, extra{:});
        mk.build();
        mk.run();
        lle(c, k) = leading_exponent(mk);
    end
end
spread = max(lle, [], 2) - min(lle, [], 2);
fprintf('      chain    benettin %+.4f  qr %+.4f  topk %+.4f   spread %.4f\n', lle(1, :), spread(1));
fprintf('      product  benettin %+.4f  qr %+.4f  topk %+.4f   spread %.4f\n', lle(2, :), spread(2));
% Both networks are STABLE (no 12- to 24-neuron chain network tried was
% clearly chaotic), and on a stable network finite-time exponents from
% different initial bases differ by up to 0.2 over seconds
% (test_lyapunov_topk's header; docs/notes/Lyapunov_estimation_methods.md).
% So 0.2 is the bound for both, and the exact statement is the next check.
all_passed = check(sprintf(['9. Benettin, QR and top-K run on the chain network; their leading exponents ' ...
    'spread by %.3f, the product twin''s by %.3f (both < 0.2)'], spread(1), spread(2)), ...
    all(isfinite(lle(:))) && all(spread < 0.2)) && all_passed;
% Liouville: the sum of ALL exponents equals the time-mean of trace(J) over
% the window, an identity of the exact flow in any regime. It checks the
% chain's Jacobian inside the propagator end to end.
mk = small_network(sc, lya{:}, 'lya_dt', 0.1, 'lya_method', 'topk', 'lya_K', 0, 'lya_K_auto', false);
mk.build();
mk.run();
in_window = find(mk.t_out >= 2.5 & mk.t_out <= 5);
pk = mk.get_params();
trace_J = zeros(numel(in_window), 1);
for k = 1:numel(in_window)
    trace_J(k) = full(sum(diag(SRNNCellTypePairs.compute_Jacobian_fast(mk.S_out(in_window(k), :)', pk))));
end
e_liouville = abs(sum(mk.lya_results.LE_spectrum) - mean(trace_J)) / abs(mean(trace_J));
all_passed = check(sprintf('9. Liouville on the chain network: sum of all %d exponents %.3f, mean trace(J) %.3f (rel %.1e)', ...
    numel(mk.lya_results.LE_spectrum), sum(mk.lya_results.LE_spectrum), mean(trace_J), e_liouville), ...
    e_liouville < 1e-3) && all_passed;
% The tight comparison needs a CHAOTIC chain network, where the leading
% exponent is isolated and a finite-time estimate does not depend on the
% initial basis. The regression network (n = 100) with its four routes
% switched to the chain is one at gain 2, noise off: over [4, 12] s top-K
% gave +1.2748 and Benettin +1.2742 on 2026-11-12 (the settings of
% regression_rapid_chain_run, section 1d; docs/notes/
% Tolerance_tier_mutation_checks_2026-11-11.md, "A chaotic chain
% network"). About a second each. The bound is 1 % of the exponent; the
% two differ by 0.05 %.
chain4 = struct('std', sc.E.E.std);
sc4 = struct('E', struct('E', chain4, 'I', chain4), 'I', struct('E', chain4, 'I', chain4));
chaos = {'rng_seeds', [5 6], 'T_range', [0 12], 'fs', 200, 'lya_T_interval', [4 12], 'lya_warmup', 3, ...
    'verbose', 'near-none', 'sigma_u_noise', 0, 'ode_solver', 'rk4', 'level_of_chaos', 2, 'synapse_config', sc4};
preset = 'regression_pairs_noStim_noise0p025_dualStd_3cond';
mt = build_from_preset(preset, 'sfa3_std2', chaos{:}, 'lya_method', 'topk', 'lya_K', 10, ...
    'lya_K_auto', false, 'lya_dt', 0.05);
mt.run();
mb = build_from_preset(preset, 'sfa3_std2', chaos{:}, 'lya_method', 'benettin');
mb.run();
lle_topk = mt.lya_results.LLE;
lle_ben = mb.lya_results.LLE;
all_passed = check(sprintf(['9. a chaotic chain network (n = 100, four chain routes, gain 2): top-K %+.4f, ' ...
    'Benettin %+.4f, apart by %.2f %% (< 1 %%)'], lle_topk, lle_ben, 100 * abs(lle_topk - lle_ben) / lle_topk), ...
    all(p_chain_routes(mt)) && lle_topk > 1 && lle_ben > 1 && ...
    abs(lle_topk - lle_ben) <= 0.01 * lle_topk) && all_passed;

fprintf('\n');
if all_passed
    fprintf('=== ALL chain model TESTS PASSED ===\n');
else
    fprintf('=== SOME chain model TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function sc = four_kinds()
% Four routes of four kinds.
sc = struct();
sc.E.E.std = struct('model', 'chain', 'tau_rec', [2 5], 'tau_rel', 0.25, 'n2', 2);
sc.E.I.std = struct('tau_rec', 2, 'tau_rel', 0.25);
sc.I.E.std = struct('model', 'chain', 'tau_rec', [1 3], 'tau_rel', 0.3, 'n2', 1.5);
sc.I.E.scale = 1.3;
end

function tf = p_chain_routes(m)
% Which routes of a built model are chains (the flags of get_params).
p = m.get_params();
tf = p.std_chain(:)';
end

function sc = product_twin(sc)
% The same network with each chain route replaced by a two-component
% product route with the chain's time constants.
sc.E.E.std = struct('tau_rec', sc.E.E.std.tau_rec, 'tau_rel', sc.E.E.std.tau_rel);
sc.I.E.std = struct('tau_rec', sc.I.E.std.tau_rec, 'tau_rel', sc.I.E.std.tau_rel);
end

function m = small_network(sc, varargin)
m = SRNNCellTypePairs('n', 24, 'indegree', 8, 'n_cellTypes', 2, 'cell_type_names', {'E', 'I'}, ...
    'f', [0.6 0.4], 'mu_tilde_relative', [3 -4], 'sigma_tilde_relative', [1 1], ...
    'tau_a', {[0.25 2], 1}, 'c', [0.1 0.05], 'synapse_config', sc, ...
    'rng_seeds', [11 12], 'T_range', [0 4], 'fs', 200, 'ode_solver', 'rk4', ...
    'lya_method', 'none', 'store_full_state', true, 'verbose', 'near-none', varargin{:});
end

function m = build(make)
m = make();
m.build();
end

function S = random_state(p, N)
% An interior state: SFA in (0, 1), every route variable in (0.05, 0.95),
% x around the logistic's centre.
lay = p.state_layout;
S = zeros(N, 1);
ia = [lay.a{:}];
ib = [lay.b{:}];
S(ia) = rand(numel(ia), 1);
S(ib) = 0.05 + 0.9 * rand(numel(ib), 1);
S(lay.x) = 0.4 + 0.6 * randn(p.n, 1);
end

function [type_of, local] = neuron_types(p)
type_of = zeros(p.n, 1);
local = zeros(p.n, 1);
for q = 1:p.n_cellTypes
    idx = p.type_indices{q};
    type_of(idx) = q;
    local(idx) = 1:numel(idx);
end
end

function r = reference_rates(S, m, p)
% r_i = phi(x_i - (c/K) sum_k a_ik), phi the logistic 1/(1 + exp(-4 (x - S_c))).
lay = p.state_layout;
[type_of, local] = neuron_types(p);
r = zeros(p.n, 1);
for i = 1:p.n
    q = type_of(i);
    nq = p.n_per_type(q);
    K = numel(m.tau_a{q});
    sum_a = 0;
    for k = 1:K
        sum_a = sum_a + S(lay.a{q}((k - 1) * nq + local(i)));
    end
    x_eff = S(lay.x(i)) - (m.c(q) / max(1, K)) * sum_a;
    r(i) = 1 / (1 + exp(-4 * (x_eff - m.S_c)));
end
end

function dS = reference_rhs(t, S, m, p, sc)
% The whole right-hand side in scalar loops, from the equations: the
% network's (docs/EquationsParametersDocs/Equations_stability_paper.md) and
% (C1), (C2), (C3) for a chain route. It reads the configuration sc and the
% drawn weights m.W, and nothing the class compiled except the layout.
lay = p.state_layout;
names = m.cell_type_names;
C = numel(names);
n = p.n;
[type_of, local] = neuron_types(p);
r = reference_rates(S, m, p);
dS = zeros(lay.N, 1);

% SFA: da_ik/dt = (r_i - a_ik)/tau_k
for i = 1:n
    q = type_of(i);
    nq = p.n_per_type(q);
    for k = 1:numel(m.tau_a{q})
        row = lay.a{q}((k - 1) * nq + local(i));
        dS(row) = (r(i) - S(row)) / m.tau_a{q}(k);
    end
end

% The routes: each presynaptic neuron's variables on each route, and the
% synaptic output theta(j, post) it sends to type `post`.
theta = zeros(n, C);
scale = ones(C, C);
for j = 1:n
    pre = type_of(j);
    npre = p.n_per_type(pre);
    for post = 1:C
        F = 1;
        route = struct();
        if isfield(sc, names{pre}) && isfield(sc.(names{pre}), names{post})
            route = sc.(names{pre}).(names{post});
        end
        if isfield(route, 'scale'), scale(pre, post) = route.scale; end
        if isfield(route, 'std')
            cfg = route.std;
            rows = lay.b{pre, post};
            if isfield(cfg, 'model') && strcmp(cfg.model, 'chain')
                row1 = rows(local(j));
                row2 = rows(npre + local(j));
                b1 = S(row1);
                b2 = S(row2);
                dS(row1) = b2 * (1 - b1) / cfg.tau_rec(1) - b1 * r(j) / cfg.tau_rel;             % (C1)
                dS(row2) = (1 - b2) / cfg.tau_rec(2) - (1 / cfg.n2) * b2 * (1 - b1) / cfg.tau_rec(1);   % (C2)
                F = b1;                                                                           % (C3)
            else
                for k = 1:numel(cfg.tau_rec)
                    row = rows((k - 1) * npre + local(j));
                    tau_rel = cfg.tau_rel(min(k, numel(cfg.tau_rel)));
                    dS(row) = (1 - S(row)) / cfg.tau_rec(k) - S(row) * r(j) / tau_rel;
                    F = F * S(row);
                end
            end
        end
        theta(j, post) = r(j) * F;
    end
end

% dx_i/dt = (-x_i + sum_j w_ij theta_j + u_i)/tau_d
u = p.u_interpolant(t);
for i = 1:n
    post = type_of(i);
    total = 0;
    for j = 1:n
        w = full(m.W(i, j));
        if w ~= 0
            total = total + scale(type_of(j), post) * w * theta(j, post);
        end
    end
    dS(lay.x(i)) = (-S(lay.x(i)) + total + u(i)) / m.tau_d;
end
end

function [b1, b2] = chain_steady_state(r, cfg)
% (C4), in the quadratic formula's own form (r > 0).
k1 = r * cfg.tau_rec(1) / cfg.tau_rel;
k2 = r * cfg.tau_rec(2) / (cfg.n2 * cfg.tau_rel);
s = 1 + k1 + k2;
b1 = (s - sqrt(s .^ 2 - 4 * k2)) ./ (2 * k2);
b2 = 1 - k2 .* b1;
end

function rows = line_rows(ax, tag)
% The YData of the lines with this tag, one row per line, sorted.
found = findall(ax, 'Type', 'line', 'Tag', tag);
rows = sortrows(vertcat(found.YData));
end

function value = leading_exponent(m)
r = m.lya_results;
if isfield(r, 'LLE')
    value = r.LLE;
else
    value = max(r.LE_spectrum);
end
end

function [same, what] = same_run(A, B)
% Everything a run produces; `what` names the first field that differs.
fields = {'W', 'S0', 'u_ex', 't_ex', 't_out', 'S_out', 'plot_data', 'lya_results'};
same = true;
what = '';
for k = 1:numel(fields)
    a = A.(fields{k});
    b = B.(fields{k});
    if strcmp(fields{k}, 'lya_results')
        a = untimed(a);
        b = untimed(b);
    end
    if ~isequaln(a, b)
        same = false;
        what = sprintf(' (first difference: %s)', fields{k});
        return
    end
end
end

function r = untimed(r)
if isfield(r, 'seconds'), r = rmfield(r, 'seconds'); end
end

function delete_if_there(file)
if isfile(file), delete(file); end
end

function e = rel_err(A, B)
scale = max(abs(B(:)));
if isempty(scale) || scale == 0, scale = 1; end
e = max(abs(A(:) - B(:))) / scale;
end

function tf = throws_id(fn, id)
tf = false;
try
    fn();
catch err
    tf = strcmp(err.identifier, id);
end
end

function ok = check(label, ok)
if ok
    fprintf('  PASS  %s\n', label);
else
    fprintf('  FAIL  %s\n', label);
end
end
