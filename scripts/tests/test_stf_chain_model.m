% test_stf_chain_model.m - release-coupled facilitation in a network.
%
% SRNNCellTypePairs with routes that ask for stf.model = 'chain' (FD-018;
% docs/Plan_A_chained_vessicles.md, section 4.6). Such a route is served by
% the vesicle-pool model route_pool and carries, per presynaptic neuron
% with rate r, a release probability g and optionally a calcium trace h,
% beside no depression, a single pool b1 or the chain b1, b2:
%
%   db1/dt = b2 (1 - b1)/tau_rec1 - g b1 r/tau_rel
%   db2/dt = (1 - b2)/tau_rec2 - (1/n2) b2 (1 - b1)/tau_rec1
%   dg/dt  = (1 - g)/tau_dec + (G - g) d/tau_fac      d = h, or r without h
%   dh/dt  = (r - h)/tau_h
%   theta  = r g b1                                   r g without depression
%
% A route carries only the variables it has: none is clamped. Without the
% reserve the refill term is (1 - b1)/tau_rec1.
%
% The route model itself is checked alone in test_route_models (every
% variant: rhs against a per-neuron loop, linearize against central
% differences, the steady state, P.coupling). This file checks it IN THE
% CLASS: that the right numbers reach the right rows of the state vector
% and the right entries of the Jacobian, beside routes of the other model.
%
% THE NETWORK: 24 neurons, three cell types E, P, S, SFA on two of them,
% and nine routes of nine kinds, one of each pool variant the class can
% select, one product route and one with no synapse dynamics:
%   E->E  chain, g, h          P->E  chain, g (and a route scale)   S->E  chain
%   E->P  single pool, g       P->P  single pool, g, h              S->P  product: one b, two g
%   E->S  g, h                 P->S  g                              S->S  none
% E's routes use Tom's values (G = 5, tau_fac = tau_rel = 0.25,
% tau_dec = 1.5, tau_h = 0.5); P's use others, so nothing equals anything.
%
% Checks (the numbers 1 to 8 are Plan A section 4.6's):
%   1. dynamics_fast equals a per-neuron loop over the equations, written
%      here from the equations above and the network's, at random states.
%   2. Steady state: unconnected neurons held at six constant rates arrive,
%      by integration, at h = r, g = (1 + G k r)/(1 + k r) and the pools
%      of the facilitation-free route at the rate g r, for each of the
%      seven pool variants; the sustained output is that route's at g r.
%   3. The analytic Jacobian equals finite_difference_jacobian at random
%      interior states and at states from a trajectory; jacobian_times
%      equals the assembled Jacobian times Y for K = 1, several and
%      Y = eye(N); the pool routes take the full-A branch.
%   4. Limits. G = 1: g stays at 1 and the network is the one without
%      facilitation. tau_h small: the route with a calcium trace
%      approaches the one without. n2 = 1e12: the reserve stays at 1 and
%      b1 and g follow the single pool with coupled facilitation.
%   5. b1, b2 in [0, 1], g in [1, G], h in [0, 1] over a run whose rates
%      swing from below 0.01 to above 0.99.
%   6. Each refusal raises its own error, at build.
%   7. params (stf_chain, stf_tau_h, n_g_pairs, the nominal rows) and the
%      starting state; the stored synaptic output is r g b1 exactly;
%      plot_data's facilitation array has columns [g, h];
%      stf_route_factor is g on a coupled route (not g h) and the product
%      on an output-only one, and plot and plot_celltypes draw it; the
%      leading-vector fractions of lya_summary sum to 1 with h counted as
%      STF; route_fields and routes_identical tell every variant apart;
%      identical routes pass routes_identical and the reservoir's
%      redundancy check and give bit-identical synaptic output.
%   8. Benettin, QR and top-K run on such a network, and Liouville holds:
%      the sum of all exponents equals the mean trace of the Jacobian.
%      On a CHAOTIC network (n = 100, every route chain, g, h, gain 2)
%      top-K and Benettin agree to 2 % on lambda_1 = +1.85.
%   9. A model with coupled facilitation reproduces bit for bit.
%
% Prints PASS/FAIL per check and a final banner. About 35 s, most of it
% the QR run of check 8, the two n = 100 runs and the plots. Assumes
% setup_paths has run.
%
% See also: route_pool, test_route_models, test_std_chain_model,
%           test_jacobian_times

fprintf('=== Testing release-coupled facilitation in SRNNCellTypePairs ===\n\n');
all_passed = true;
rng_at_entry = rng;
restore_rng = onCleanup(@() rng(rng_at_entry));
old_visibility = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
restore_visibility = onCleanup(@() set(groot, 'DefaultFigureVisible', old_visibility));

sc = nine_kinds();
m = small_network(sc, 'lya_method', 'topk', 'lya_K', 5, 'lya_K_auto', false, 'lya_dt', 0.1, ...
    'lya_T_interval', [2 4], 'lya_warmup', 1);
m.build();
m.run();
p = m.get_params();
p.u_interpolant = m.u_interpolant;
lay = p.state_layout;
N = lay.N;
nE = p.n_per_type(1);
nP = p.n_per_type(2);
nS = p.n_per_type(3);
pd = m.plot_data;

%% 7a. params and the starting state
ok = isequal(p.n_b_pairs, [2 1 0; 2 1 0; 2 1 0]) && isequal(p.n_g_pairs, [2 1 2; 1 2 1; 0 2 0]) && ...
    isequal(p.stf_chain, [true true true; true true true; false false false]) && ...
    isequaln(p.stf_tau_h, [0.5 NaN 0.5; NaN 0.2 NaN; NaN NaN NaN]) && ...
    isequal(p.std_chain, [true false false; true false false; true false false]) && ...
    isequal(p.tau_g_dec{1, 1}, 1.5) && isequal(p.tau_g_fac{1, 1}, 0.25) && isequal(p.G{1, 1}, 5) && ...
    isequal(p.tau_g_dec{2, 3}, 0.9) && isequal(p.G{2, 2}, 2.86) && isequal(p.tau_g_dec{3, 2}, [1.5 4]) && ...
    isequal(p.tau_b_rec{1, 2}, 2) && isequal(p.tau_b_rec{1, 1}, [2 5]) && isempty(p.tau_b_rec{1, 3}) && ...
    strcmp(p.route{1, 3}.P.model, 'pool') && strcmp(p.route{3, 2}.P.model, 'product') && ...
    strcmp(p.route{3, 3}.P.model, 'product') && strcmp(p.route{3, 1}.P.model, 'pool') && ...
    numel(lay.g{1, 1}) == 2 * nE && numel(lay.g{1, 2}) == nE && numel(lay.g{2, 2}) == 2 * nP && ...
    numel(lay.b{1, 1}) == 2 * nE && isempty(lay.b{1, 3}) && isempty(lay.g{3, 1}) && ...
    ~p.route{1, 1}.A_is_diagonal && ~p.route{1, 2}.A_is_diagonal && ~p.route{1, 3}.A_is_diagonal && ...
    p.route{2, 3}.A_is_diagonal && p.route{3, 2}.A_is_diagonal;
all_passed = check('7. params: n_b_pairs, n_g_pairs, stf_chain, stf_tau_h, the nominal rows, the state layout', ok) && all_passed;
[is_g, is_h] = facilitation_rows(p);
ok = all(m.S0([lay.b{:}]) == 1) && all(m.S0(is_g) == 1) && all(m.S0(is_h) == 0) && ...
    nnz(is_h) == 2 * nE + nP && nnz(is_g) == 3 * nE + 3 * nP + 2 * nS;
all_passed = check('7. every pool and every g start at 1, every calcium trace at 0', ok) && all_passed;

%% 1. dynamics_fast against the per-neuron loop
rng(201, 'twister');
e_rhs = 0;
for trial = 1:5
    S = random_state(p, N);
    t = 0.3 + 0.7 * trial;
    e_rhs = max(e_rhs, rel_err(SRNNCellTypePairs.dynamics_fast(t, S, p), reference_rhs(t, S, m, p, sc)));
end
all_passed = check(sprintf('1. dynamics_fast equals the per-neuron loop of the equations and the network (rel %.1e)', e_rhs), ...
    e_rhs <= 1e-13) && all_passed;

%% 2. steady state, by integration
% Six unconnected neurons, each held at a constant rate, per pool variant.
drive = [-0.2; 0.1; 0.3; 0.4; 0.6; 1.2];
variants = {'chain, g, h', sc.E.E; 'single pool, g', sc.E.P; 'g, h', sc.E.S; 'chain, g', rmfield(sc.P.E, 'scale'); ...
    'single pool, g, h', sc.P.P; 'g', sc.P.S; 'chain', sc.S.E};
e_int = zeros(1, size(variants, 1));
e_out = zeros(1, size(variants, 1));
span = zeros(1, size(variants, 1));
pool = route_model('pool');
for v = 1:size(variants, 1)
    route = variants{v, 2};
    u = SRNNCellTypePairs('n', 6, 'indegree', 1, 'n_cellTypes', 1, 'cell_type_names', {'E'}, 'f', 1, ...
        'mu_tilde_relative', 0, 'sigma_tilde_relative', 0, 'tau_a', {[]}, 'c', 0, ...
        'synapse_config', struct('E', struct('E', route)), ...
        'T_range', [0 130], 'fs', 20, 'ode_solver', 'rk4', 'lya_method', 'none', ...
        'store_full_state', true, 'verbose', 'near-none');
    u.input_config.amp = 0;
    u.input_config.intrinsic_drive = drive;
    u.build();
    u.run();
    r_held = u.plot_data.r.E(:, end);
    span(v) = max(r_held) - min(r_held);
    [b_ss, g_ss, rho] = reference_steady_state(r_held, route);
    found = zeros(6, 0);
    if ~isempty(u.plot_data.b.E.E), found = u.plot_data.b.E.E(:, :, end); end
    if ~isempty(u.plot_data.g.E.E), found = [found, u.plot_data.g.E.E(:, :, end)]; end %#ok<AGROW>
    e_int(v) = max(abs(found - [b_ss, g_ss]), [], 'all');
    % The sustained output is the facilitation-free route's at the rate
    % rho = g r: rho times that route's steady factor at rho (rho itself
    % where the route has no depression).
    if isfield(route, 'std')
        plain = pool.steady_state(rho, pool.compile(route.std, [], struct('npre', 6)));
        expected = rho .* plain.F;
    else
        expected = rho;
    end
    e_out(v) = max(abs(u.plot_data.synaptic_output.E.E(:, end) - expected));
    ok_w = nnz(u.W) == 0;
end
all_passed = check(sprintf(['2. unconnected neurons at six constant rates arrive at h = r, g = (1 + G k r)/(1 + k r) and ' ...
    'the pools at the rate g r, for %d pool variants (max %.1e)'], size(variants, 1), max(e_int)), ...
    ok_w && max(e_int) <= 1e-8 && min(span) > 0.8) && all_passed;
all_passed = check(sprintf('2. the sustained output is the facilitation-free route''s at the rate g r (max %.1e)', max(e_out)), ...
    max(e_out) <= 1e-8) && all_passed;

%% 7b. stored synaptic output, and the facilitation array
g_of = @(pre, post, n) reshape(pd.g.(pre).(post)(:, 1, :), n, []);
b1_of = @(pre, post, n) reshape(pd.b.(pre).(post)(:, 1, :), n, []);
ok = isequal(pd.synaptic_output.E.E, pd.r.E .* (g_of('E', 'E', nE) .* b1_of('E', 'E', nE))) && ...
    isequal(pd.synaptic_output.E.P, pd.r.E .* (g_of('E', 'P', nE) .* b1_of('E', 'P', nE))) && ...
    isequal(pd.synaptic_output.E.S, pd.r.E .* g_of('E', 'S', nE)) && ...
    isequal(pd.synaptic_output.P.E, pd.r.P .* (g_of('P', 'E', nP) .* b1_of('P', 'E', nP))) && ...
    isequal(pd.synaptic_output.P.P, pd.r.P .* (g_of('P', 'P', nP) .* b1_of('P', 'P', nP))) && ...
    isequal(pd.synaptic_output.P.S, pd.r.P .* g_of('P', 'S', nP)) && ...
    isequal(pd.synaptic_output.S.E, pd.r.S .* b1_of('S', 'E', nS)) && ...
    isequal(pd.synaptic_output.S.S, pd.r.S) && ...
    ~isequal(pd.synaptic_output.E.S, pd.r.E .* reshape(prod(pd.g.E.S, 2), nE, []));
all_passed = check('7. stored synaptic output is r g b1 (r g without depression, r b1 without facilitation), exactly, and not r g h', ok) && all_passed;
h_EE = reshape(pd.g.E.E(:, 2, :), nE, []);
ok = isequal(size(pd.g.E.E, 2), 2) && isequal(size(pd.g.E.P, 2), 1) && isequal(size(pd.g.P.P, 2), 2) && ...
    all(h_EE(:, 1) == 0) && all(g_of('E', 'E', nE) >= 1, 'all') && max(h_EE, [], 'all') <= 1 && ...
    max(h_EE, [], 'all') > 0.05 && isempty(pd.g.S.E) && isempty(pd.g.S.S);
all_passed = check('7. plot_data.g has columns [g, h] on a route with a calcium trace, [g] without, [] without facilitation', ok) && all_passed;

%% 3. the analytic Jacobian against finite differences, and jacobian_times
rng(202, 'twister');
states = cell(1, 6);
for k = 1:3
    states{k} = random_state(p, N);
end
at = round(size(m.S_out, 1) * [0.3 0.6 1]);
for k = 1:3
    states{3 + k} = m.S_out(at(k), :)';
end
is_route = false(N, 1);
is_route([lay.b{:}, lay.g{:}]) = true;
e_fd = zeros(1, 6);
e_fd_route = zeros(1, 6);
bg_block = 0;
for k = 1:6
    J = full(SRNNCellTypePairs.compute_Jacobian_fast(states{k}, p));
    J_fd = SRNNCellTypePairs.finite_difference_jacobian(states{k}, p, 1e-6);
    e_fd(k) = max(abs(J - J_fd), [], 'all') / max(abs(J_fd), [], 'all');
    e_fd_route(k) = max(abs(J(is_route, is_route) - J_fd(is_route, is_route)), [], 'all');
    % The block that is new in kind: d(db1/dt)/dg, a facilitation column in a depression row.
    bg_block = max(bg_block, max(abs(J([lay.b{:}], [lay.g{:}])), [], 'all'));
    gb_block = max(abs(J([lay.g{:}], [lay.b{:}])), [], 'all');
end
all_passed = check(sprintf('3. analytic Jacobian equals finite differences at 3 random states (rel %.1e)', max(e_fd(1:3))), ...
    max(e_fd(1:3)) <= 1e-7) && all_passed;
all_passed = check(sprintf('3. and at 3 states from a trajectory (rel %.1e)', max(e_fd(4:6))), ...
    max(e_fd(4:6)) <= 1e-7) && all_passed;
all_passed = check(sprintf(['3. the (b, g) block is live (entries up to %.2f) and the (g, b) block is zero, so the routes'' ' ...
    'block is not symmetric: block error %.1e'], bg_block, max(e_fd_route)), ...
    bg_block > 0.05 && gb_block == 0 && max(e_fd_route) <= 1e-7) && all_passed;

rng(203, 'twister');
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
all_passed = check(sprintf('3. jacobian_times equals J * Y for K = 1 (rel %.1e), K = 7 (%.1e), Y = eye(N) (%.1e)', ...
    e_jt(1), e_jt(2), e_jt(3)), max(e_jt) <= 1e-13) && all_passed;
all_passed = check(sprintf('3. the coupled routes take the full-A branch: applying only the diagonal is wrong by %.1e at least', ...
    e_diag_only), e_diag_only > 1e-3) && all_passed;

%% 4. limits
% a. G = 1 on every coupled route: g cannot move, and the network is the
% one whose coupled routes have no facilitation at all.
flat = sc;
bare = sc;
for pair = {'E', 'E'; 'E', 'P'; 'E', 'S'; 'P', 'E'; 'P', 'P'; 'P', 'S'}'
    flat.(pair{1}).(pair{2}).stf.G = 1;
    bare.(pair{1}).(pair{2}) = rmfield(bare.(pair{1}).(pair{2}), 'stf');
end
m_flat = run_of(small_network(flat));
m_bare = run_of(small_network(bare));
pf = m_flat.plot_data;
pb = m_bare.plot_data;
[flat_g, ~] = facilitation_rows(m_flat.get_params());
g_flat = m_flat.S_out(:, flat_g & coupled_rows(m_flat.get_params()));
e_x = max(abs([pf.x.E; pf.x.P; pf.x.S] - [pb.x.E; pb.x.P; pb.x.S]), [], 'all');
e_b = max([max(abs(pf.b.E.E - pb.b.E.E), [], 'all'), max(abs(pf.b.E.P - pb.b.E.P), [], 'all'), ...
    max(abs(pf.b.P.E - pb.b.P.E), [], 'all'), max(abs(pf.b.P.P - pb.b.P.P), [], 'all')]);
depth = 1 - min(pb.b.E.E(:, 1, :), [], 'all');
all_passed = check(sprintf(['4. G = 1: g stays at 1 exactly, and x and the pools follow the network without ' ...
    'facilitation (x %.1e, pools %.1e) while b1 falls by %.2f'], e_x, e_b, depth), ...
    isequal(m_flat.W, m_bare.W) && all(g_flat(:) == 1) && e_x <= 1e-12 && e_b <= 1e-12 && depth > 0.1) && all_passed;

% b. tau_h -> 0: the calcium trace is the rate, and the route approaches
% the one without a trace. fs = 500 so that rk4 resolves tau_h = 2e-3.
no_h = sc;
for pair = {'E', 'E'; 'E', 'S'; 'P', 'P'}'
    no_h.(pair{1}).(pair{2}).stf = rmfield(no_h.(pair{1}).(pair{2}).stf, 'tau_h');
end
% The error is first order in tau_h (the trace lags the rate by tau_h):
% measured 2026-11-12, 1.6e-1 at tau_h = 2e-2 and 1.8e-2 at 2e-3 (and
% 8.7e-2, 9.3e-3 at 1e-2, 1e-3 with fs = 1000).
short = {'fs', 500, 'T_range', [0 2.5]};
m_no_h = run_of(small_network(no_h, short{:}));
e_h = zeros(1, 2);
tau_small = [2e-2 2e-3];
for k = 1:2
    quick = sc;
    for pair = {'E', 'E'; 'E', 'S'; 'P', 'P'}'
        quick.(pair{1}).(pair{2}).stf.tau_h = tau_small(k);
    end
    m_quick = run_of(small_network(quick, short{:}));
    pq = m_quick.plot_data;
    pn = m_no_h.plot_data;
    e_h(k) = max([max(abs(pq.g.E.E(:, 1, :) - pn.g.E.E), [], 'all'), max(abs(pq.g.E.S(:, 1, :) - pn.g.E.S), [], 'all'), ...
        max(abs(pq.g.P.P(:, 1, :) - pn.g.P.P), [], 'all'), max(abs(pq.b.E.E - pn.b.E.E), [], 'all'), ...
        max(abs(pq.x.E - pn.x.E), [], 'all')]);
end
moved = max(m_no_h.plot_data.g.E.E, [], 'all') - 1;
all_passed = check(sprintf(['4. tau_h = %g, %g: g, the pools and x approach the route without a calcium trace, ' ...
    'first order in tau_h (max %.1e, %.1e) while g rises by %.2f'], tau_small, e_h, moved), ...
    e_h(2) < e_h(1) / 5 && e_h(2) <= 4e-2 && moved > 0.5) && all_passed;

% c. n2 = 1e12: the reserve stays full and b1, g follow the single pool
% with coupled facilitation.
big = sc;
big.E.E.std.n2 = 1e12;
big.P.E.std.n2 = 1e12;
one = sc;
one.E.E.std = struct('tau_rec', sc.E.E.std.tau_rec(1), 'tau_rel', sc.E.E.std.tau_rel);
one.P.E.std = struct('tau_rec', sc.P.E.std.tau_rec(1), 'tau_rel', sc.P.E.std.tau_rel);
m_big = run_of(small_network(big));
m_one = run_of(small_network(one));
bb = m_big.plot_data;
bo = m_one.plot_data;
e_b1 = max([max(abs(bb.b.E.E(:, 1, :) - bo.b.E.E), [], 'all'), max(abs(bb.b.P.E(:, 1, :) - bo.b.P.E), [], 'all')]);
e_g = max([max(abs(bb.g.E.E - bo.g.E.E), [], 'all'), max(abs(bb.g.P.E - bo.g.P.E), [], 'all')]);
e_b2 = max([max(abs(bb.b.E.E(:, 2, :) - 1), [], 'all'), max(abs(bb.b.P.E(:, 2, :) - 1), [], 'all')]);
all_passed = check(sprintf(['4. n2 = 1e12: b1 and g follow the single pool with coupled facilitation (b1 %.1e, g %.1e); ' ...
    'b2 stays at 1 (max %.1e)'], e_b1, e_g, e_b2), ...
    isequal(m_big.W, m_one.W) && e_b1 <= 1e-9 && e_g <= 1e-9 && e_b2 <= 1e-10 && ...
    strcmp(m_one.get_params().route{1, 1}.P.model, 'pool') && m_one.get_params().n_b_pairs(1, 1) == 1) && all_passed;

%% 5. every variable stays in its range under large rate swings
m_swing = small_network(sc, 'T_range', [0 12], 'rng_seeds', [21 22]);
m_swing.input_config.amp = 4;
m_swing.input_config.n_steps = 8;
m_swing.input_config.no_stim_pattern = logical([1 0 0 1 0 0 1 0]);
m_swing.build();
m_swing.run();
ps = m_swing.get_params();
[sw_g, sw_h] = facilitation_rows(ps);
coupled = coupled_rows(ps);
S_sw = m_swing.S_out;
pools = S_sw(:, [ps.state_layout.b{:}]);
g_all = S_sw(:, sw_g & coupled);
h_all = S_sw(:, sw_h);
g_E = S_sw(:, ps.state_layout.g{1, 1}(1:ps.n_per_type(1)));             % E->E, G = 5
g_P = S_sw(:, ps.state_layout.g{2, 2}(1:ps.n_per_type(2)));             % P->P, G = 2.86
rate_all = [m_swing.plot_data.r.E(:); m_swing.plot_data.r.P(:); m_swing.plot_data.r.S(:)];
all_passed = check(sprintf(['5. rates from %.3f to %.3f: pools in [%.4f, %.4f], g in [%.4f, %.4f] (G = 5) and ' ...
    '[%.4f, %.4f] (G = 2.86), h in [%.4f, %.4f]'], min(rate_all), max(rate_all), min(pools(:)), max(pools(:)), ...
    min(g_E(:)), max(g_E(:)), min(g_P(:)), max(g_P(:)), min(h_all(:)), max(h_all(:))), ...
    all(pools(:) >= 0 & pools(:) <= 1) && all(g_all(:) >= 1) && all(g_E(:) <= 5) && all(g_P(:) <= 2.86) && ...
    all(h_all(:) >= 0 & h_all(:) <= 1) && min(rate_all) < 0.01 && max(rate_all) > 0.99 && ...
    max(g_E(:)) > 4 && max(h_all(:)) > 0.9 && min(pools(:)) < 0.1) && all_passed;

%% 6. the refusals
fac = sc.E.P.stf;                              % coupled, no calcium trace
fac_h = sc.E.E.stf;
on_EP = @(route) setfield(sc, 'E', setfield(sc.E, 'P', route));
std_one = sc.E.P.std;
std_two = struct('tau_rec', [2 5], 'tau_rel', 0.25);
output_only = struct('tau_dec', 1.5, 'tau_fac', 0.25, 'G', 5);
set_in = @(cfg, name, value) setfield(cfg, name, value);
refusals = { ...
    'a chain with output-only stf',           struct('std', sc.E.E.std, 'stf', output_only), 'ChainStfOutputOnly'
    'a chain with stf.model = ''product''',    struct('std', sc.E.E.std, 'stf', set_in(output_only, 'model', 'product')), 'ChainStfOutputOnly'
    'coupled stf beside a two-component product', struct('std', std_two, 'stf', fac), 'ChainStfProductStd'
    'coupled stf beside std.combine',         struct('std', set_in(std_one, 'combine', 'geomean'), 'stf', fac), 'ChainCombine'
    'a vector tau_dec',                       struct('std', std_one, 'stf', set_in(fac, 'tau_dec', [1.5 4])), 'ChainStfNotScalar'
    'G below one',                            struct('std', std_one, 'stf', set_in(fac, 'G', 0.9)), 'ChainStfInvalidG'
    'tau_h not positive',                     struct('stf', set_in(fac_h, 'tau_h', 0)), 'ChainStfInvalidTauH'
    'an unknown stf.model',                   struct('std', std_one, 'stf', set_in(fac, 'model', 'chian')), 'UnknownStfModel'
    'an unknown field of a coupled stf',      struct('std', std_one, 'stf', set_in(fac, 'tau_x', 1)), 'InvalidSynapseConfig'};
ok = true;
ids = cell(1, size(refusals, 1));
for k = 1:size(refusals, 1)
    route = refusals{k, 2};
    ids{k} = refusals{k, 3};
    ok_k = throws_id(@() build(@() small_network(on_EP(route))), ['SRNNCellTypePairs:' ids{k}]);
    if ~ok_k, fprintf('      %s did not raise SRNNCellTypePairs:%s at build\n', refusals{k, 1}, ids{k}); end
    ok = ok && ok_k;
end
named = false;
try
    build(@() small_network(on_EP(refusals{1, 2})));
catch err
    named = contains(err.message, 'stf.model = ''chain''');
end
only_E = struct('E', struct('E', struct('std', std_one, 'stf', fac)));
ok_floor = throws_id(@() build(@() small_network(sc, 'std_zero_floor', true)), 'SRNNCellTypePairs:ChainZeroFloor') && ...
    throws_id(@() build(@() small_network(only_E, 'std_zero_floor', true)), 'SRNNCellTypePairs:ChainZeroFloor');
all_passed = check(sprintf(['6. %d refusals each raise their own identifier at build (%d distinct); the zero floor beside a ' ...
    'coupled route with depression is refused; the chain''s output-only refusal names stf.model = ''chain'''], ...
    size(refusals, 1), numel(unique(ids))), ok && ok_floor && named && ...
    numel(unique(ids)) == size(refusals, 1) - 1) && all_passed;
explicit = sc;
explicit.S.P.stf.model = 'product';
explicit.E.P.std.model = 'product';
pe = build(@() small_network(explicit)).get_params();
all_passed = check('6. stf.model = ''product'' and std.model = ''product'' are the same routes as no model field', ...
    isequal(pe.route{3, 2}.P, p.route{3, 2}.P) && isequal(pe.route{1, 2}.P, p.route{1, 2}.P) && ~pe.stf_chain(3, 2)) && all_passed;
facil_only = struct('E', struct('S', struct('stf', fac_h)));
p_fo = build(@() small_network(facil_only, 'std_zero_floor', true)).get_params();
all_passed = check('6. the zero floor is not refused where the only pool route has no depression (it is not in force there)', ...
    strcmp(p_fo.route{1, 3}.P.model, 'pool') && p_fo.n_b_pairs(1, 3) == 0) && all_passed;

%% 7c. the display factor, and the class's plots
g3 = pd.g.E.E;                                  % nE x 2 x nt: [g, h]
F3 = SRNNCellTypePairs.stf_route_factor(g3, p, 1, 1);
F2 = SRNNCellTypePairs.stf_route_factor(g3(:, :, end), p, 1, 1);
Fp = SRNNCellTypePairs.stf_route_factor(pd.g.S.P, p, 3, 2);
ok = isequal(F3, g3(:, 1, :)) && isequal(size(F3), size(prod(g3, 2))) && isequal(F2, g3(:, 1, end)) && ...
    ~isequal(F3, prod(g3, 2)) && isequal(Fp, prod(pd.g.S.P, 2)) && size(pd.g.S.P, 2) == 2 && ...
    isequal(SRNNCellTypePairs.stf_route_factor(pd.g.E.P, p, 1, 2), pd.g.E.P) && ...
    isempty(SRNNCellTypePairs.stf_route_factor(pd.g.S.E, p, 3, 1)) && ...
    throws_id(@() SRNNCellTypePairs.stf_route_factor(g3(:, 1, :), p, 1, 1), 'SRNNCellTypePairs:RouteFactorShape') && ...
    isequal(SRNNCellTypePairs.std_route_factor(pd.b.E.P, p, 1, 2), pd.b.E.P) && ...
    isequal(SRNNCellTypePairs.std_route_factor(pd.b.E.E, p, 1, 1), pd.b.E.E(:, 1, :));
all_passed = check(['7. stf_route_factor: g on a coupled route (2-D and 3-D; not g h), the product on an output-only ' ...
    'route, [] on none; std_route_factor: b1 on every pool route'], ok) && all_passed;
[fig, ax] = m.plot_celltypes();
stf_row = 7;                                    % u, x, r, synaptic output, SFA, STD, STF, Lyapunov
drawn_EE = line_rows(ax(stf_row, 1), 'PostType_E');
drawn_ES = line_rows(ax(stf_row, 1), 'PostType_S');
drawn_SP = line_rows(ax(stf_row, 3), 'PostType_P');
ok = isequal(drawn_EE, sortrows(g_of('E', 'E', nE))) && isequal(drawn_ES, sortrows(g_of('E', 'S', nE))) && ...
    isequal(drawn_SP, sortrows(reshape(prod(pd.g.S.P, 2), nS, []))) && ...
    ~isequal(drawn_EE, sortrows(reshape(prod(pd.g.E.E, 2), nE, []))) && ...
    strcmp(ax(stf_row, 1).YLabel.String, 'STF factor') && strcmp(ax(stf_row - 1, 1).YLabel.String, 'STD factor');
close(fig);
all_passed = check('7. plot_celltypes draws g for the coupled routes (not g h) and prod(g) for the output-only route', ok) && all_passed;
[fig, ax1] = m.plot();
mean_lines = findall(ax1(7), 'Type', 'line');
drawn = sortrows(vertcat(mean_lines.YData));
expected = sortrows([mean(g_of('E', 'E', nE), 1); mean(g_of('E', 'P', nE), 1); mean(g_of('E', 'S', nE), 1); ...
    mean(g_of('P', 'E', nP), 1); mean(g_of('P', 'P', nP), 1); mean(g_of('P', 'S', nP), 1); ...
    mean(reshape(prod(pd.g.S.P, 2), nS, []), 1)]);
ok = isequal(drawn, expected) && strcmp(ax1(7).YLabel.String, 'STF factor');
close(fig);
all_passed = check('7. plot draws each route''s mean facilitation factor', ok) && all_passed;
old_stf = struct('tau_dec', [1.5 4], 'tau_fac', [0.25 1], 'G', [3 1.5]);
all_old = struct('E', struct('E', struct('stf', old_stf), 'P', struct('stf', old_stf)));
m_old = run_of(small_network(all_old));
[fig, axo] = m_old.plot_celltypes();
drawn = line_rows(axo(6, 1), 'PostType_E');      % no STD row on this network
ok = isequal(drawn, sortrows(reshape(prod(m_old.plot_data.g.E.E, 2), nE, []))) && ...
    strcmp(axo(6, 1).YLabel.String, 'STF (product)');
close(fig);
all_passed = check('7. an output-only network draws prod(g) under the label it had', ok) && all_passed;

%% 7d. the leading-vector fractions
S_lya = m.lya_summary();
total = S_lya.lead_frac_x + S_lya.lead_frac_sfa + S_lya.lead_frac_std + S_lya.lead_frac_stf;
e_unit = zeros(N, 1);
e_unit(find(is_h, 1)) = 1;
on_h = SRNNCellTypePairs.leading_vector_blocks(e_unit, lay);
all_passed = check(sprintf(['7. lya_summary''s leading-vector fractions sum to 1 (off by %.1e; STF %.3f); a vector on a ' ...
    'calcium trace alone is all STF'], abs(total - 1), S_lya.lead_frac_stf), ...
    abs(total - 1) <= 1e-12 && S_lya.lead_frac_stf > 0 && on_h.stf == 1 && on_h.std == 0 && on_h.x == 0) && all_passed;

%% 7e. route_fields, routes_identical, the reservoir's redundancy check
same = struct('std', sc.E.E.std, 'stf', sc.E.E.stf);
all_same = uniform_config(same);
other_h = all_same;
other_h.P.S.stf.tau_h = 0.6;
without_h = all_same;
without_h.P.S.stf = rmfield(without_h.P.S.stf, 'tau_h');
other_G = all_same;
other_G.S.E.stf.G = 4;
no_fac = all_same;
no_fac.E.S = rmfield(no_fac.E.S, 'stf');
one_pool_coupled = uniform_config(struct('std', std_one, 'stf', fac));
one_pool_output_only = uniform_config(struct('std', std_one, 'stf', output_only));
params_of = @(cfg) build(@() small_network(cfg)).get_params();
identical = @(cfg) SRNNCellTypePairs.routes_identical(params_of(cfg));
all_passed = check(['7. routes_identical: true for nine identical coupled routes; false when one differs in tau_h, ' ...
    'lacks the calcium trace, differs in G, or lacks facilitation'], ...
    identical(all_same) && identical(one_pool_coupled) && identical(one_pool_output_only) && ...
    ~identical(other_h) && ~identical(without_h) && ~identical(other_G) && ~identical(no_fac) && ...
    ~SRNNCellTypePairs.routes_identical(p)) && all_passed;
f_coupled = SRNNCellTypePairs.route_fields(params_of(one_pool_coupled), 1, 1);
f_output = SRNNCellTypePairs.route_fields(params_of(one_pool_output_only), 1, 1);
f_all = SRNNCellTypePairs.route_fields(p, 1, 1);
ok = ~isequal(f_coupled, f_output) && isequal(f_coupled.tau_dec, f_output.tau_dec) && isequal(f_coupled.G, f_output.G) && ...
    isequal(f_coupled.tau_rec, f_output.tau_rec) && f_coupled.n_g == f_output.n_g && ...
    strcmp(f_coupled.signature.model, 'pool') && strcmp(f_output.signature.model, 'product') && ...
    f_all.n_g == 2 && f_all.signature.tau_h == 0.5 && strcmp(f_all.signature.stf_model, 'chain');
all_passed = check(['7. route_fields tells a coupled route from an output-only one with the same numbers, ' ...
    'and carries tau_h'], ok) && all_passed;
redundant = true;
try
    SRNN_ESN_reservoir.assert_route_redundancy(params_of(all_same));
    SRNN_ESN_reservoir.assert_route_redundancy(params_of(one_pool_coupled));
catch
    redundant = false;
end
ok = redundant && ...
    throws_id(@() SRNN_ESN_reservoir.assert_route_redundancy(params_of(other_h)), 'SRNN_ESN_reservoir:AmbiguousSynapticOutput') && ...
    throws_id(@() SRNN_ESN_reservoir.assert_route_redundancy(params_of(without_h)), 'SRNN_ESN_reservoir:AmbiguousSynapticOutput');
all_passed = check(['7. the reservoir''s redundancy check passes identical coupled routes and refuses a type whose ' ...
    'routes differ in the calcium trace'], ok) && all_passed;
m_same = run_of(small_network(all_same));
so = m_same.plot_data.synaptic_output;
ok = isequal(so.E.E, so.E.P) && isequal(so.E.E, so.E.S) && isequal(so.P.E, so.P.S) && isequal(so.S.E, so.S.P) && ...
    max(m_same.plot_data.g.E.E(:, 1, :), [], 'all') > 1.5;
all_passed = check('7. identical coupled routes of one presynaptic type give bit-identical synaptic output', ok) && all_passed;

%% 9. reproducibility
rng('shuffle');
B = small_network(sc, 'lya_method', 'topk', 'lya_K', 5, 'lya_K_auto', false, 'lya_dt', 0.1, ...
    'lya_T_interval', [2 4], 'lya_warmup', 1);
B.build();
B.run();
[same_runs, what] = same_run(m, B);
all_passed = check(sprintf('9. two fresh models with coupled facilitation are identical%s (lambda_1 = %+.3f)', ...
    what, m.lya_results.LLE), same_runs && isfinite(m.lya_results.LLE)) && all_passed;
file = [tempname '.mat'];
remove_file = onCleanup(@() delete_if_there(file));
save(file, 'B');
loaded = load(file);
L = loaded.B;
ok = isequal(L.S_out, m.S_out) && isequal(L.synapse_config, sc);
L.build();
L.run();
[same_runs, what] = same_run(L, m);
all_passed = check(sprintf('9. saved, loaded, built and run again: identical%s', what), ok && same_runs) && all_passed;

%% 8. the three Lyapunov estimators
% On 6 neurons, two cell types and a 1.5 s window, because the ode45 QR
% integrates N^2 variational equations and is the slowest part of this
% file. E->E chain, g, h; I->E g; the other two routes none.
sl = struct();
sl.E.E = sc.E.E;
sl.I.E = sc.P.S;
lya = {'n', 6, 'indegree', 2, 'T_range', [0 3], 'lya_T_interval', [1.5 3], 'lya_warmup', 1, 'lya_dt', 0.25};
methods = {'benettin', 'qr', 'topk'};
lle = zeros(1, 3);
for k = 1:3
    extra = {};
    if strcmp(methods{k}, 'topk'), extra = {'lya_K', 5, 'lya_K_auto', false}; end
    mk = lyapunov_network(sl, lya{:}, 'lya_method', methods{k}, extra{:});
    mk.build();
    mk.run();
    lle(k) = leading_exponent(mk);
end
spread = max(lle) - min(lle);
fprintf('      benettin %+.4f  qr %+.4f  topk %+.4f   spread %.4f\n', lle, spread);
% The network is STABLE, with a nearly degenerate band of slow exponents,
% and on such a network finite-time exponents depend on the initial basis
% and converge slowly (test_lyapunov_topk's header;
% docs/notes/Lyapunov_estimation_methods.md). Measured 2026-11-12 on the
% 8-neuron version of this network (top-K -0.45, -0.47, -0.53): Benettin
% and top-K are 0.20 apart over a 2.5 s window, 0.07 over [10, 20] s and
% 0.04 over [20, 40] s. So this check asserts that the three run, agree
% in sign and stay within 0.5 (0.27 here); the exact statements are the
% two checks that follow.
all_passed = check(sprintf(['8. Benettin, QR and top-K run on a network with coupled facilitation: all negative, ' ...
    'spread %.3f over 1.5 s (< 0.5)'], spread), all(isfinite(lle)) && all(lle < 0) && spread < 0.5) && all_passed;
% Liouville: the sum of ALL exponents equals the time-mean of trace(J) over
% the window, an identity of the exact flow in any regime. It checks the
% coupled routes' Jacobian inside the propagator end to end.
mk = lyapunov_network(sl, lya{:}, 'lya_dt', 0.1, 'lya_method', 'topk', 'lya_K', 0, 'lya_K_auto', false);
mk.build();
mk.run();
in_window = find(mk.t_out >= 1.5 & mk.t_out <= 3);
pk = mk.get_params();
trace_J = zeros(numel(in_window), 1);
for k = 1:numel(in_window)
    trace_J(k) = full(sum(diag(SRNNCellTypePairs.compute_Jacobian_fast(mk.S_out(in_window(k), :)', pk))));
end
e_liouville = abs(sum(mk.lya_results.LE_spectrum) - mean(trace_J)) / abs(mean(trace_J));
all_passed = check(sprintf('8. Liouville: sum of all %d exponents %.3f, mean trace(J) %.3f (rel %.1e)', ...
    numel(mk.lya_results.LE_spectrum), sum(mk.lya_results.LE_spectrum), mean(trace_J), e_liouville), ...
    e_liouville < 1e-3) && all_passed;
% The tight comparison needs a CHAOTIC network, where the leading exponent
% is isolated: the regression network (n = 100) with all four routes a
% chain with coupled facilitation and a calcium trace, gain 2, noise off.
% Measured 2026-11-12 over [4, 12] s: top-K +1.8498, Benettin +1.8597
% (0.5 % apart; the bound is 2 %). About 2.5 s each.
route = sc.E.E;
sc4 = struct('E', struct('E', route, 'I', route), 'I', struct('E', route, 'I', route));
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
pt = mt.get_params();
all_passed = check(sprintf(['8. a chaotic network (n = 100, four routes with chain, g and h, gain 2): top-K %+.4f, ' ...
    'Benettin %+.4f, apart by %.2f %% (< 2 %%)'], lle_topk, lle_ben, 100 * abs(lle_topk - lle_ben) / lle_topk), ...
    all(pt.stf_chain(:)) && all(pt.n_g_pairs(:) == 2) && lle_topk > 1 && lle_ben > 1 && ...
    abs(lle_topk - lle_ben) <= 0.02 * lle_topk) && all_passed;

fprintf('\n');
if all_passed
    fprintf('=== ALL coupled facilitation TESTS PASSED ===\n');
else
    fprintf('=== SOME coupled facilitation TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function sc = nine_kinds()
% Nine routes of nine kinds (see the header).
chain_E = struct('model', 'chain', 'tau_rec', [2 5], 'tau_rel', 0.25, 'n2', 2);
chain_P = struct('model', 'chain', 'tau_rec', [1 3], 'tau_rel', 0.3, 'n2', 1.5);
fac_E = struct('model', 'chain', 'tau_dec', 1.5, 'tau_fac', 0.25, 'G', 5);
fac_P = struct('model', 'chain', 'tau_dec', 0.9, 'tau_fac', 0.4, 'G', 2.86);
with_h = @(cfg, tau_h) setfield(cfg, 'tau_h', tau_h);
sc = struct();
sc.E.E.std = chain_E;
sc.E.E.stf = with_h(fac_E, 0.5);
sc.E.P.std = struct('tau_rec', 2, 'tau_rel', 0.25);
sc.E.P.stf = fac_E;
sc.E.S.stf = with_h(fac_E, 0.5);
sc.P.E.std = chain_P;
sc.P.E.stf = fac_P;
sc.P.E.scale = 1.3;
sc.P.P.std = struct('tau_rec', 1.5, 'tau_rel', 0.3);
sc.P.P.stf = with_h(fac_P, 0.2);
sc.P.S.stf = fac_P;
sc.S.E.std = chain_E;
sc.S.P.std = struct('tau_rec', 2, 'tau_rel', 0.25);
sc.S.P.stf = struct('tau_dec', [1.5 4], 'tau_fac', [0.25 1], 'G', [3 1.5]);
end

function sc = uniform_config(route)
% The same route configuration on all nine routes.
sc = struct();
for pre = {'E', 'P', 'S'}
    for post = {'E', 'P', 'S'}
        sc.(pre{1}).(post{1}) = route;
    end
end
end

function m = small_network(sc, varargin)
m = SRNNCellTypePairs('n', 24, 'indegree', 8, 'n_cellTypes', 3, 'cell_type_names', {'E', 'P', 'S'}, ...
    'f', [0.5 0.25 0.25], 'mu_tilde_relative', [3 -4 -3], 'sigma_tilde_relative', [1 1 1], ...
    'tau_a', {[0.25 2], 1, []}, 'c', [0.1 0.05 0], 'synapse_config', sc, ...
    'rng_seeds', [11 12], 'T_range', [0 4], 'fs', 200, 'ode_solver', 'rk4', ...
    'lya_method', 'none', 'store_full_state', true, 'verbose', 'near-none', varargin{:});
end

function m = lyapunov_network(sc, varargin)
m = SRNNCellTypePairs('n_cellTypes', 2, 'cell_type_names', {'E', 'I'}, ...
    'f', [0.6 0.4], 'mu_tilde_relative', [3 -4], 'sigma_tilde_relative', [1 1], ...
    'tau_a', {[0.25 2], 1}, 'c', [0.1 0.05], 'synapse_config', sc, ...
    'rng_seeds', [11 12], 'fs', 200, 'ode_solver', 'rk4', ...
    'store_full_state', true, 'verbose', 'near-none', varargin{:});
end

function m = build(make)
m = make();
m.build();
end

function m = run_of(m)
m.build();
m.run();
end

function [is_g, is_h] = facilitation_rows(p)
% Which rows of the state vector are a facilitation variable that
% multiplies the output (g, and an output-only route's components) and
% which are a calcium trace, from params.stf_chain and params.stf_tau_h.
lay = p.state_layout;
is_g = false(lay.N, 1);
is_h = false(lay.N, 1);
for pre = 1:p.n_cellTypes
    npre = p.n_per_type(pre);
    for post = 1:p.n_cellTypes
        rows = lay.g{pre, post};
        if isempty(rows), continue; end
        if p.stf_chain(pre, post) && ~isnan(p.stf_tau_h(pre, post))
            is_g(rows(1:npre)) = true;
            is_h(rows(npre + 1:end)) = true;
        else
            is_g(rows) = true;
        end
    end
end
end

function tf = coupled_rows(p)
% The g-block rows of the release-coupled routes.
lay = p.state_layout;
tf = false(lay.N, 1);
for pre = 1:p.n_cellTypes
    for post = 1:p.n_cellTypes
        if p.stf_chain(pre, post), tf(lay.g{pre, post}) = true; end
    end
end
end

function S = random_state(p, N)
% An interior state: SFA in (0, 1), pools in (0.05, 0.95), g in (1, 2.8),
% calcium traces in (0, 1), x around the logistic's centre.
lay = p.state_layout;
S = zeros(N, 1);
ia = [lay.a{:}];
ib = [lay.b{:}];
[is_g, is_h] = facilitation_rows(p);
S(ia) = rand(numel(ia), 1);
S(ib) = 0.05 + 0.9 * rand(numel(ib), 1);
S(is_g) = 1 + 1.8 * rand(nnz(is_g), 1);
S(is_h) = rand(nnz(is_h), 1);
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
% the pool route's (the header). It reads the configuration sc and the
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
        b_rows = lay.b{pre, post};
        g_rows = lay.g{pre, post};
        coupled = isfield(route, 'stf') && isfield(route.stf, 'model') && strcmp(route.stf.model, 'chain');
        is_chain = isfield(route, 'std') && isfield(route.std, 'model') && strcmp(route.std.model, 'chain');
        if coupled || is_chain
            % A vesicle-pool route. Its variables, by name.
            use_factor = 1;                                % what multiplies b1 r/tau_rel
            if coupled
                cfg = route.stf;
                row_g = g_rows(local(j));
                g = S(row_g);
                if isfield(cfg, 'tau_h')
                    row_h = g_rows(npre + local(j));
                    h = S(row_h);
                    dS(row_g) = (1 - g) / cfg.tau_dec + (cfg.G - g) * h / cfg.tau_fac;
                    dS(row_h) = (r(j) - h) / cfg.tau_h;
                else
                    dS(row_g) = (1 - g) / cfg.tau_dec + (cfg.G - g) * r(j) / cfg.tau_fac;
                end
                use_factor = g;
                F = F * g;
            end
            if isfield(route, 'std')
                cfg = route.std;
                row1 = b_rows(local(j));
                b1 = S(row1);
                if is_chain
                    row2 = b_rows(npre + local(j));
                    b2 = S(row2);
                    dS(row1) = b2 * (1 - b1) / cfg.tau_rec(1) - use_factor * b1 * r(j) / cfg.tau_rel;
                    dS(row2) = (1 - b2) / cfg.tau_rec(2) - (1 / cfg.n2) * b2 * (1 - b1) / cfg.tau_rec(1);
                else
                    dS(row1) = (1 - b1) / cfg.tau_rec - use_factor * b1 * r(j) / cfg.tau_rel;
                end
                F = F * b1;
            end
        else
            % A product route: independent components, each a factor of the output.
            if isfield(route, 'std')
                cfg = route.std;
                for k = 1:numel(cfg.tau_rec)
                    row = b_rows((k - 1) * npre + local(j));
                    tau_rel = cfg.tau_rel(min(k, numel(cfg.tau_rel)));
                    dS(row) = (1 - S(row)) / cfg.tau_rec(k) - S(row) * r(j) / tau_rel;
                    F = F * S(row);
                end
            end
            if isfield(route, 'stf')
                cfg = route.stf;
                for k = 1:numel(cfg.tau_dec)
                    row = g_rows((k - 1) * npre + local(j));
                    dS(row) = (1 - S(row)) / cfg.tau_dec(k) + (cfg.G(k) - S(row)) * r(j) / cfg.tau_fac(k);
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

function [b_ss, g_ss, rho] = reference_steady_state(r, route)
% The steady state of a pool route at constant rates r (column), from the
% closed forms of Plan A section 4.6, the pools by the quadratic formula's
% own form. b_ss has the route's b columns, g_ss its g columns; rho is the
% effective rate g r the pools see.
b_ss = zeros(numel(r), 0);
g_ss = zeros(numel(r), 0);
rho = r;
if isfield(route, 'stf')
    cfg = route.stf;
    kappa = cfg.tau_dec / cfg.tau_fac;
    g = (1 + cfg.G * kappa * r) ./ (1 + kappa * r);
    rho = g .* r;
    g_ss = g;
    if isfield(cfg, 'tau_h'), g_ss = [g, r]; end
end
if isfield(route, 'std')
    cfg = route.std;
    k1 = rho * cfg.tau_rec(1) / cfg.tau_rel;
    if numel(cfg.tau_rec) == 2
        k2 = rho * cfg.tau_rec(2) / (cfg.n2 * cfg.tau_rel);
        s = 1 + k1 + k2;
        b1 = (s - sqrt(s .^ 2 - 4 * k2)) ./ (2 * k2);
        b_ss = [b1, 1 - k2 .* b1];
    else
        b_ss = 1 ./ (1 + k1);
    end
end
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
