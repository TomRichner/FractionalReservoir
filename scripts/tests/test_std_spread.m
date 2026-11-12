% test_std_spread.m - the per-neuron spread on the synapse parameters (FD-015).
%
% std_spread and stf_spread (1 x C, by PRESYNAPTIC cell type) give every
% neuron its own depression and facilitation parameters:
%   value = nominal * exp(spread * z),   G = 1 + (G_nominal - 1) * exp(spread * z)
% with z one row of standard normals per neuron and mechanism (std_draws,
% stf_draws), drawn at build() from a stream of its own and read by every
% route of that neuron, each quantity from a fixed column (route_model,
% "THE SPREAD"). The default is no spread. Plan A, sections 4.3, 4.6 and 5;
% the route models alone are checked in test_route_models (2g, 3d).
%
% Checks, for depression and for facilitation:
%   1. ZERO SPREAD IS NO SPREAD. A model that states std_spread = 0 and
%      stf_spread = 0 (and seeds) has the same params (every field, every
%      route's P) and the same run, bit for bit, as one that never mentions
%      them; nothing is drawn. A spread on a mechanism no route has draws
%      but changes no parameter and no number of the run.
%   2. THE DRAWS. They are randn(n, Q) from rng(seed, 'twister'), the seed
%      rng_seeds(1) + 32452843 (depression) or + 49979687 (facilitation)
%      unless std_seed / stf_seed is set; reproducible; moved by the seed
%      and by rng_seeds(1); each independent of whether the other
%      mechanism is spread; the leading columns do not depend on Q. W, the
%      setpoints, the SFA ladders, the stimulus, S0 and the global RNG
%      state after build() are untouched by turning either spread on,
%      ALSO with an input generator that does not seed the random stream:
%      the stock generator seeds it itself, so with it the stimulus and
%      x0 cannot show whether the draws put the stream back.
%   3. 1TS AND CHAINED MTS on one network seed have identical tau_rec_1
%      and tau_rel per neuron, equal to nominal * exp(spread(pre) * z) of
%      that neuron's own row; likewise a coupled single pool and a coupled
%      chain, and their facilitation parameters.
%   4. A NEURON'S ROUTES CARRY IDENTICAL VALUES: the per-neuron parameters
%      of a type's routes are equal, routes_identical is true (false once
%      the two types' spreads differ), the reservoir's redundancy check
%      passes, the synaptic outputs of a neuron's routes are bit-identical
%      along a run, and a memory-capacity run with the synaptic readout
%      works on a spread chain network (with facilitation too).
%   5. THE VALUES. Positive; G >= 1; at n = 1600 the log-mean and
%      log-spread of every jittered quantity match the nominal and the
%      type's spread within sampling error; the chain's ordering floor
%      holds on every neuron, tau_rec_1 is untouched by it, and
%      std_clipped counts the neurons moved (once, however many of their
%      routes are chains).
%   6. THE JACOBIANS with the spread on, spreads differing between types:
%      analytic against central finite differences, and jacobian_times
%      against the assembled Jacobian, on a network with a chain, a
%      one-component and a two-component product route side by side, on
%      the four-type network with every pool variant that has
%      facilitation, and on a product network with the zero floor.
%   7. PRODUCT MTS AND CHAINED MTS on one seed share tau_rec_1, tau_rel_1
%      and the standard normal behind tau_rec_2.
%   8. dynamics_fast on a spread network equals a per-neuron loop of the
%      equations with each neuron's own values, formed here from the
%      draws.
%   9. The class's surface: a scalar spread is broadcast; a negative,
%      non-finite or wrong-length spread and a bad seed error
%      InvalidParams; std_route_factor and stf_route_factor refuse an
%      array with fewer rows than the route has neurons once its
%      parameters differ per neuron (RouteFactorRows) and serve it
%      otherwise; a model changed after build() errors StaleSpreadDraws
%      until built again; the nominal params fields stay nominal; the
%      spreads are ordinary settable properties (so ParamSpaceAnalysis2
%      takes them as model defaults).
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths.
%
% See also: SRNNCellTypePairs, route_model, route_pool, route_product,
%           test_route_models, test_SRNNCellTypePairs_tau_a_spread

fprintf('=== Testing the per-neuron spread on the synapse parameters (std_spread, stf_spread) ===\n\n');
all_passed = true;
OFF_STD = 32452843;
OFF_STF = 49979687;

one   = struct('tau_rec', 2, 'tau_rel', 0.25);
two   = struct('tau_rec', [2 5], 'tau_rel', [0.25 0.5]);
chain = struct('model', 'chain', 'tau_rec', [2 5], 'tau_rel', 0.25, 'n2', 2);
fac   = struct('model', 'chain', 'tau_dec', 1.5, 'tau_fac', 0.25, 'G', 5);
fac_h = struct('model', 'chain', 'tau_dec', 1.5, 'tau_fac', 0.25, 'G', 5, 'tau_h', 0.5);
sc_one   = four_routes(struct('std', one));
sc_two   = four_routes(struct('std', two));
sc_chain = four_routes(struct('std', chain));
sc_one_f   = four_routes(struct('std', one, 'stf', fac_h));
sc_chain_f = four_routes(struct('std', chain, 'stf', fac_h));

%% 1. Zero spread is no spread
fprintf('1. zero spread\n');
for kind = {sc_chain_f, sc_two}
    plain = net(kind{1});
    zero = net(kind{1}, 'std_spread', 0, 'stf_spread', [0 0], 'std_seed', 7, 'stf_seed', 8);
    p_plain = plain.get_params();
    p_zero = zero.get_params();
    plain.run();
    zero.run();
    ok = isequaln(no_handles(p_plain), no_handles(p_zero)) && isempty(zero.std_draws) && isempty(zero.stf_draws) && ...
        isequal(zero.std_clipped, [0 0]) && isequal(zero.std_spread, [0 0]) && isequal(plain.std_spread, [0 0]) && ...
        isequal(plain.S_out, zero.S_out) && isequal(plain.S0, zero.S0) && isequaln(plain.plot_data, zero.plot_data);
    all_passed = check(sprintf('params (every route''s P) and a run are bit-identical with the spreads stated as zero (%s routes)', ...
        p_plain.route{1, 1}.P.model), ok) && all_passed;
end
% a spread on a mechanism no route has: drawn, and nothing else moves
plain = net(sc_two);
idle = net(sc_two, 'stf_spread', 0.3);
plain.run();
idle.run();
p_plain = plain.get_params();
p_idle = idle.get_params();
same_P = true;
for k = 1:4
    same_P = same_P && isequaln(p_plain.route{k}.P, p_idle.route{k}.P);
end
all_passed = check('stf_spread on a network without facilitation: draws exist, every route''s P and the run are unchanged', ...
    isequal(size(idle.stf_draws), [idle.n 4]) && same_P && isequal(plain.S_out, idle.S_out)) && all_passed;

%% 2. The draws
fprintf('2. the draws and their streams\n');
het = {'sigma_S_c', 0.05, 'tau_a_spread', 0.2};
rng(123, 'twister');
A = net(sc_chain_f, het{:});
state_A = rng;
rng(123, 'twister');
B = net(sc_chain_f, het{:}, 'std_spread', 0.25, 'stf_spread', [0.2 0.3]);
state_B = rng;
ok = isequal(A.W, B.W) && isequal(A.S_c_vec, B.S_c_vec) && isequal(A.tau_a_matrix, B.tau_a_matrix) && ...
    isequal(A.u_ex, B.u_ex) && isequal(A.S0, B.S0) && isequal(state_A, state_B) && ...
    ~isempty(A.S_c_vec) && ~isempty(A.tau_a_matrix);
all_passed = check('W, the setpoints, the SFA ladders, the stimulus, S0 and the RNG state after build() are untouched by the spreads', ok) && all_passed;
% The check above cannot see whether the draws put the global stream back:
% the stock input generator seeds the stream itself, so the stimulus and
% x0 never depend on what build_network left behind. With a generator that
% does NOT seed, they do (the tolerance tier's spread suite found this: a
% dropped restore passed every case until one had such an input).
unseeded = struct('generator', @(params, T, fs, ~, ~) deal(0.3 * randn(params.n, 1) .* ones(1, round(T * fs) + 1), ...
    (0:1 / fs:T)'), 'intrinsic_drive', 0);
A_u = net(sc_chain_f, het{:}, 'input_config', unseeded);
B_u = net(sc_chain_f, het{:}, 'input_config', unseeded, 'std_spread', 0.25, 'stf_spread', [0.2 0.3]);
C_u = net(sc_chain_f, het{:}, 'input_config', unseeded, 'rng_seeds', [3 2]);
all_passed = check('with an input generator that does not seed the stream, the stimulus and S0 are still untouched: the draws restore the global RNG state', ...
    isequal(A_u.u_ex, B_u.u_ex) && isequal(A_u.S0, B_u.S0) && ~isequal(A_u.u_ex, C_u.u_ex) && ...
    ~isequal(A_u.u_ex, A.u_ex)) && all_passed;

seed1 = B.rng_seeds(1);
keep = rng;
rng(seed1 + OFF_STD, 'twister'); want_std = randn(B.n, 5);
rng(seed1 + OFF_STF, 'twister'); want_stf = randn(B.n, 4);
rng(keep);
all_passed = check('std_draws and stf_draws are randn(n, Q) from rng(rng_seeds(1) + offset, ''twister''), Q = 5 and 4', ...
    isequal(B.std_draws, want_std) && isequal(B.stf_draws, want_stf) && ~isequal(want_std(:, 1:4), want_stf)) && all_passed;

B2 = net(sc_chain_f, het{:}, 'std_spread', 0.25, 'stf_spread', [0.2 0.3]);
B_std = net(sc_chain_f, 'std_spread', 0.25);                    % depression alone
B_stf = net(sc_chain_f, 'stf_spread', [0.2 0.3]);               % facilitation alone
B_seed = net(sc_chain_f, 'std_spread', 0.25, 'stf_spread', 0.2, 'std_seed', 99);
B_fseed = net(sc_chain_f, 'std_spread', 0.25, 'stf_spread', 0.2, 'stf_seed', 99);
B_net = net(sc_chain_f, 'std_spread', 0.25, 'stf_spread', 0.2, 'rng_seeds', [11 2]);
ok = isequal(B.std_draws, B2.std_draws) && isequal(B.stf_draws, B2.stf_draws) && ...
    isequaln(no_handles(B.get_params()), no_handles(B2.get_params())) && ...
    isequal(B_std.std_draws, B.std_draws) && isempty(B_std.stf_draws) && ...
    isequal(B_stf.stf_draws, B.stf_draws) && isempty(B_stf.std_draws) && ...
    ~isequal(B_seed.std_draws, B.std_draws) && isequal(B_seed.stf_draws, B.stf_draws) && ...
    isequal(B_fseed.std_draws, B.std_draws) && ~isequal(B_fseed.stf_draws, B.stf_draws) && ...
    isequal(B_seed.std_draws(:, 1:4), B_fseed.stf_draws) && ...
    ~isequal(B_net.std_draws, B.std_draws) && ~isequal(B_net.stf_draws, B.stf_draws);
all_passed = check(['same seed same draws; std_seed moves depression''s only, stf_seed facilitation''s only, ' ...
    'rng_seeds(1) both; each stream is the same whether or not the other is drawn'], ok) && all_passed;

% Q grows with the largest route; the leading columns do not move
sc_big = sc_chain;
sc_big.E.I.std = struct('tau_rec', [2 5 9], 'tau_rel', 0.3);
sc_big.I.I = struct('stf', struct('tau_dec', [0.4 1.2], 'tau_fac', [0.2 0.3], 'G', [1.5 2]));
big = net(sc_big, 'std_spread', 0.25, 'stf_spread', 0.2);
all_passed = check('a three-component product route makes Q = 7 and two facilitation components Q = 7; the leading columns are the same numbers', ...
    isequal(size(big.std_draws), [big.n 7]) && isequal(size(big.stf_draws), [big.n 7]) && ...
    isequal(big.std_draws(:, 1:5), B.std_draws) && isequal(big.stf_draws(:, 1:4), B.stf_draws)) && all_passed;

%% 3. 1TS and chained MTS share the releasable pool's parameters
fprintf('3. single-timescale and chained conditions on one seed\n');
spread = {'std_spread', [0.25 0.4], 'stf_spread', [0.3 0.2]};
m1 = net(sc_one, spread{:});
mc = net(sc_chain, spread{:});
p1 = m1.get_params();
pc = mc.get_params();
ok = true;
exact = true;
for pre = 1:2
    z = mc.std_draws(mc.type_indices{pre}, :);
    s = mc.std_spread(pre);
    for post = 1:2
        P1 = p1.route{pre, post}.P;
        Pc = pc.route{pre, post}.P;
        ok = ok && isequal(P1.tau_rest, Pc.tau_rec1) && isequal(P1.tau_drive, Pc.tau_rel);
        exact = exact && isequal(Pc.tau_rec1, 2 * exp(s * z(:, 1))) && isequal(Pc.tau_rel, 0.25 * exp(s * z(:, 2))) && ...
            isequal(Pc.n2, 2 * exp(s * z(:, 5)));
    end
end
all_passed = check('1TS and chained MTS: identical tau_rec_1 and tau_rel per neuron on every route', ...
    ok && isequal(m1.std_draws, mc.std_draws) && isequal(m1.W, mc.W)) && all_passed;
all_passed = check('each is nominal * exp(spread(pre) * z) of the neuron''s own row and its type''s spread; n2 from column 5', exact) && all_passed;

m1f = net(sc_one_f, spread{:});
mcf = net(sc_chain_f, spread{:});
p1f = m1f.get_params();
pcf = mcf.get_params();
ok = true;
for pre = 1:2
    z = mcf.stf_draws(mcf.type_indices{pre}, :);
    s = mcf.stf_spread(pre);
    for post = 1:2
        P1 = p1f.route{pre, post}.P;
        Pc = pcf.route{pre, post}.P;
        ok = ok && isequal(P1.tau_rec1, Pc.tau_rec1, pc.route{pre, post}.P.tau_rec1) && isequal(P1.tau_rel, Pc.tau_rel) && ...
            isequal(P1.tau_dec, Pc.tau_dec, 1.5 * exp(s * z(:, 1))) && isequal(P1.tau_fac, Pc.tau_fac, 0.25 * exp(s * z(:, 2))) && ...
            isequal(P1.tau_h, Pc.tau_h, 0.5 * exp(s * z(:, 3))) && isequal(P1.G, Pc.G, 1 + 4 * exp(s * z(:, 4))) && ...
            all(Pc.G >= 1);
    end
end
all_passed = check(['a coupled single pool and a coupled chain: identical pool and facilitation parameters per neuron; ' ...
    'tau_dec, tau_fac, tau_h and G - 1 from columns 1 to 4 of the facilitation draws'], ok) && all_passed;

%% 4. A neuron's routes carry identical values
fprintf('4. a neuron''s routes\n');
ok = true;
for pre = 1:2
    Pa = pcf.route{pre, 1}.P;
    Pb = pcf.route{pre, 2}.P;
    ok = ok && isequaln(Pa, Pb) && Pa.jittered;
end
same_spread = net(sc_chain_f, 'std_spread', 0.25, 'stf_spread', 0.25);
p_same = same_spread.get_params();
all_passed = check('the two routes of each presynaptic type have equal per-neuron parameters', ok) && all_passed;
all_passed = check('routes_identical: true with one spread for both types, false when the types'' spreads differ, true with none', ...
    SRNNCellTypePairs.routes_identical(p_same) && ~SRNNCellTypePairs.routes_identical(pcf) && ...
    SRNNCellTypePairs.routes_identical(net(sc_chain_f).get_params())) && all_passed;
f11 = SRNNCellTypePairs.route_fields(pcf, 1, 1);
f21 = SRNNCellTypePairs.route_fields(pcf, 2, 1);
f_plain = SRNNCellTypePairs.route_fields(p1, 1, 1);
all_passed = check('route_fields states the presynaptic type''s spread of each mechanism the route has', ...
    f11.std_spread == 0.25 && f11.stf_spread == 0.3 && f21.std_spread == 0.4 && f21.stf_spread == 0.2 && ...
    f_plain.std_spread == 0.25 && f_plain.stf_spread == 0) && all_passed;
all_passed = check('the reservoir''s redundancy check passes with the spread on (per presynaptic type)', ...
    no_error(@() SRNN_ESN_reservoir.assert_route_redundancy(pcf)) && ...
    no_error(@() SRNN_ESN_reservoir.assert_route_redundancy(pc))) && all_passed;
mcf.run();
[~, ~, ~, ~, syn] = SRNNCellTypePairs.unpack_and_compute_states(mcf.S_out, pcf);
plain_f = net(sc_chain_f);
plain_f.run();
all_passed = check('along a run the synaptic outputs of a neuron''s two routes are bit-identical, and the spread moved the run', ...
    isequal(syn.E.E, syn.E.I) && isequal(syn.I.E, syn.I.I) && ~isequal(mcf.S_out, plain_f.S_out) && ...
    min(syn.E.E(:)) < 0.9 * max(syn.E.E(:))) && all_passed;

mc_args = {'n', 60, 'indegree', 20, 'n_cellTypes', 2, 'cell_type_names', {'E', 'I'}, 'f', [0.6 0.4], ...
    'mu_tilde_relative', [3 -4; 3 -4], 'sigma_tilde_relative', ones(2, 2), 'level_of_chaos', 1.5, 'S_c', 0.35, ...
    'fs', 100, 'ode_solver', 'rk4', 'T_wash', 300, 'T_train', 1500, 'T_test', 900, 'd_max', 30, ...
    'input_type', 'sample_hold', 'T_hold', 0.1, 'lya_method', 'none', 'tau_a', {[0.25 1.58 10], []}};
for kind = {sc_chain, sc_chain_f}
    esn = SRNN_ESN_reservoir(mc_args{:}, 'synapse_config', kind{1}, 'std_spread', [0.25 0.25], 'stf_spread', 0.25);
    evalc('esn.build(); [MC, R2] = esn.run_memory_capacity(''readout_signal'', ''synaptic'');');
    esn0 = SRNN_ESN_reservoir(mc_args{:}, 'synapse_config', kind{1});
    evalc('esn0.build(); MC0 = esn0.run_memory_capacity(''readout_signal'', ''synaptic'');');
    all_passed = check(sprintf(['memory capacity, synaptic readout, spread chain network (%d route variables): ' ...
        'MC = %.4f, finite, and not the unspread network''s %.4f'], esn.cached_params.route{1}.P.ns, MC, MC0), ...
        isfinite(MC) && all(isfinite(R2)) && MC > 0 && MC ~= MC0 && isequal(size(esn.std_draws), [60 5]) && ...
        isequal(esn.W, esn0.W)) && all_passed;
end

%% 5. The values
fprintf('5. the drawn values at n = 1600\n');
sc_stat = struct();
sc_stat.E.E = struct('std', chain, 'stf', fac_h);
sc_stat.E.I = struct('std', struct('model', 'chain', 'tau_rec', [2 3], 'tau_rel', 0.25, 'n2', 2));
sc_stat.I.E = struct('std', two, 'stf', struct('tau_dec', [0.4 1.2], 'tau_fac', [0.2 0.3], 'G', [1.5 2]));
sc_stat.I.I = struct('std', one);
big = SRNNCellTypePairs('n', 1600, 'indegree', 100, 'n_cellTypes', 2, 'cell_type_names', {'E', 'I'}, ...
    'f', [0.5 0.5], 'mu_tilde_relative', [3 -4], 'sigma_tilde_relative', [1 1], 'synapse_config', sc_stat, ...
    'std_spread', [0.25 0.4], 'stf_spread', [0.3 0.15], 'T_range', [0 0.1], 'lya_method', 'none', 'verbose', 'near-none');
big.build();
pb = big.get_params();
EE = pb.route{1, 1}.P;
EI = pb.route{1, 2}.P;
IE = pb.route{2, 1}.P;
II = pb.route{2, 2}.P;
% {values, nominal, sigma}; of a ceiling, the excess over 1
Q = {EE.tau_rec1, 2, 0.25; EE.tau_rel, 0.25, 0.25; EE.n2, 2, 0.25; EI.n2, 2, 0.25; ...
    EE.tau_dec, 1.5, 0.3; EE.tau_fac, 0.25, 0.3; EE.tau_h, 0.5, 0.3; EE.G - 1, 4, 0.3; ...
    IE.tau_rest(:, 1), 2, 0.4; IE.tau_rest(:, 2), 5, 0.4; IE.tau_drive(:, 1), 0.25, 0.4; IE.tau_drive(:, 2), 0.5, 0.4; ...
    IE.tau_rest(:, 3), 0.4, 0.15; IE.tau_rest(:, 4), 1.2, 0.15; IE.tau_drive(:, 4), 0.3, 0.15; ...
    IE.target(:, 3) - 1, 0.5, 0.15; IE.target(:, 4) - 1, 1, 0.15; II.tau_rest, 2, 0.4; II.tau_drive, 0.25, 0.4};
nq = 800;
worst_mean = 0;
worst_std = 0;
positive = true;
for k = 1:size(Q, 1)
    v = Q{k, 1};
    sig = Q{k, 3};
    positive = positive && all(v > 0) && isequal(size(v), [nq 1]);
    worst_mean = max(worst_mean, abs(mean(log(v / Q{k, 2}))) / (sig / sqrt(nq)));
    worst_std = max(worst_std, abs(std(log(v / Q{k, 2})) / sig - 1) * sqrt(2 * nq));
end
all_passed = check(sprintf(['%d jittered quantities on four kinds of route: all positive; log-mean within 4 standard ' ...
    'errors of the nominal (worst %.2f) and log-spread within 4 of the type''s spread (worst %.2f)'], ...
    size(Q, 1), worst_mean, worst_std), positive && worst_mean <= 4 && worst_std <= 4) && all_passed;
all_passed = check('every ceiling is at least 1 (and above it: the excess is jittered, not G)', ...
    all(EE.G > 1) && all(IE.target(:, 3:4) > 1, 'all') && min(EE.G) < 4) && all_passed;
% the ordering floor: E->E has nominal log-span log(5/2), E->I log(3/2)
zE = big.std_draws(big.type_indices{1}, :);
gap_drawn = 0.25 * (zE(:, 3) - zE(:, 1));
want_EE = log(5 / 2) + gap_drawn < log(5 / 2) / 4;
want_EI = log(3 / 2) + gap_drawn < log(3 / 2) / 4;
gap_EE = log(EE.tau_rec2 ./ EE.tau_rec1);
gap_EI = log(EI.tau_rec2 ./ EI.tau_rec1);
ok = isequal(EE.clipped, want_EE) && isequal(EI.clipped, want_EI) && ...
    all(gap_EE >= log(5 / 2) / 4 - 1e-12) && all(gap_EI >= log(3 / 2) / 4 - 1e-12) && ...
    max(abs(gap_EE(want_EE) - log(5 / 2) / 4)) <= 1e-12 && ...
    isequal(EE.tau_rec1, 2 * exp(0.25 * zE(:, 1))) && isequal(EI.tau_rec1, EE.tau_rec1) && ...
    isequal(EE.tau_rec2(~want_EE), 5 * exp(0.25 * zE(~want_EE, 3))) && ...
    all(EE.tau_rec2(want_EE) > 5 * exp(0.25 * zE(want_EE, 3)));
all_passed = check(sprintf(['the ordering floor holds on every chain neuron (%d and %d of %d raised to it on the two ' ...
    'chain routes); only tau_rec_2 moved, tau_rec_1 is as drawn'], nnz(want_EE), nnz(want_EI), nq), ...
    ok && nnz(want_EE) > 5 && nnz(want_EI) > nnz(want_EE)) && all_passed;
all_passed = check(sprintf('std_clipped counts each moved neuron once: [%d %d]', big.std_clipped), ...
    isequal(big.std_clipped, [nnz(want_EE | want_EI), 0])) && all_passed;

%% 6. The Jacobians with the spread on
fprintf('6. Jacobians\n');
sc_kinds = struct();
sc_kinds.E.E.std = chain;
sc_kinds.E.I.std = one;
sc_kinds.I.E.std = two;
sc_kinds.I.I.std = struct('model', 'chain', 'tau_rec', [1 1.5], 'tau_rel', 0.3, 'n2', 1.5);
kinds = net(sc_kinds, 'std_spread', [0.25 0.4], 'stf_spread', 0.2);
mixed = pool_variants_network({'std_spread', [0.25 0.4 0.1 0.3], 'stf_spread', [0.3 0.2 0.25 0.15]});
floored = floor_network({'std_spread', [0.3 0.2 0.25], 'stf_spread', [0.2 0.3 0.1]});
labels = {'a chain, a one-component and a two-component product route', ...
    'every pool variant with facilitation, beside product routes (four types)', ...
    'product routes with the geometric mean, facilitation and the zero floor'};
models = {kinds, mixed, floored};
for q = 1:numel(models)
    m = models{q};
    m.run();
    p = m.get_params();
    p.u_interpolant = m.u_interpolant;
    n_jit = sum(cellfun(@(R) R.P.jittered, p.route(:)));
    N = size(m.S_out, 2);
    nt = size(m.S_out, 1);
    rng(40 + q, 'twister');
    e_fd = 0;
    e_jt = 0;
    for i = round(nt * [0.3 0.65 1])
        S = m.S_out(i, :)';
        J = SRNNCellTypePairs.compute_Jacobian_fast(S, p);
        J_fd = SRNNCellTypePairs.finite_difference_jacobian(S, p, 1e-6);
        Y = randn(N, 5);
        Z = SRNNCellTypePairs.jacobian_times(S, Y, p);
        e_fd = max(e_fd, max(abs(J(:) - J_fd(:))) / max(abs(J(:))));
        e_jt = max(e_jt, max(abs(Z - J * Y), [], 'all') / max(abs(Z(:))));
    end
    all_passed = check(sprintf('%s (%d spread routes, N = %d): analytic vs finite differences %.1e, jacobian_times vs J * Y %.1e', ...
        labels{q}, n_jit, N, e_fd, e_jt), n_jit >= 3 && e_fd <= 2e-6 && e_jt <= 1e-12) && all_passed;
end

%% 7. Product MTS and chained MTS
fprintf('7. product and chained multiple-timescale conditions on one seed\n');
mp = net(sc_two, spread{:});
pp = mp.get_params();
ok = isequal(mp.std_draws, mc.std_draws);
n_free = 0;
for pre = 1:2
    s = mc.std_spread(pre);
    z = mc.std_draws(mc.type_indices{pre}, :);
    for post = 1:2
        Pp = pp.route{pre, post}.P;
        Pc = pc.route{pre, post}.P;
        free = ~Pc.clipped;
        n_free = n_free + nnz(free);
        ok = ok && isequal(Pp.tau_rest(:, 1), Pc.tau_rec1) && isequal(Pp.tau_drive(:, 1), Pc.tau_rel) && ...
            isequal(Pp.tau_rest(:, 2), 5 * exp(s * z(:, 3))) && isequal(Pc.tau_rec2(free), Pp.tau_rest(free, 2)) && ...
            isequal(Pp.tau_drive(:, 2), 0.5 * exp(s * z(:, 4)));
    end
end
all_passed = check(sprintf(['product MTS and chained MTS: identical tau_rec_1 and tau_rel_1 per neuron, and tau_rec_2 ' ...
    'from the same standard normal (equal on the %d route neurons the floor did not move)'], n_free), ok) && all_passed;

%% 8. dynamics_fast against a per-neuron loop on each neuron's own values
fprintf('8. the right-hand side, neuron by neuron\n');
p = pcf;
p.u_interpolant = mcf.u_interpolant;
S = mcf.S_out(round(end / 2), :)';
dS = SRNNCellTypePairs.dynamics_fast(mcf.t_out(round(end / 2)), S, p);
[~, ~, ~, r] = SRNNCellTypePairs.unpack_and_compute_states(S', p);
names = mcf.cell_type_names;
e = 0;
for pre = 1:2
    rows = mcf.type_indices{pre};
    zs = mcf.std_draws(rows, :);
    zf = mcf.stf_draws(rows, :);
    ss = mcf.std_spread(pre);
    sf = mcf.stf_spread(pre);
    rate = r.(names{pre});
    for post = 1:2
        idx = p.route{pre, post}.idx;
        s = reshape(S(idx), numel(rows), 4);                 % [b1, b2, g, h]
        want = zeros(numel(rows), 4);
        for i = 1:numel(rows)
            tau_rec1 = 2 * exp(ss * zs(i, 1));
            tau_rel = 0.25 * exp(ss * zs(i, 2));
            tau_rec2 = 5 * exp(ss * zs(i, 3));
            if log(tau_rec2 / tau_rec1) < log(5 / 2) / 4
                tau_rec2 = tau_rec1 * exp(log(5 / 2) / 4);
            end
            n2 = 2 * exp(ss * zs(i, 5));
            tau_dec = 1.5 * exp(sf * zf(i, 1));
            tau_fac = 0.25 * exp(sf * zf(i, 2));
            tau_h = 0.5 * exp(sf * zf(i, 3));
            G = 1 + 4 * exp(sf * zf(i, 4));
            b1 = s(i, 1); b2 = s(i, 2); g = s(i, 3); h = s(i, 4);
            want(i, 1) = b2 * (1 - b1) / tau_rec1 - g * b1 * rate(i) / tau_rel;
            want(i, 2) = (1 - b2) / tau_rec2 - (1 / n2) * b2 * (1 - b1) / tau_rec1;
            want(i, 3) = (1 - g) / tau_dec + (G - g) * h / tau_fac;
            want(i, 4) = (rate(i) - h) / tau_h;
        end
        e = max(e, max(abs(reshape(dS(idx), numel(rows), 4) - want), [], 'all'));
    end
end
all_passed = check(sprintf('coupled chain network: the routes'' rows of dynamics_fast equal the per-neuron loop (max %.1e)', e), ...
    e <= 1e-13) && all_passed;

%% 9. The class's surface
fprintf('9. validation and display\n');
m = net(sc_chain, 'std_spread', 0.25, 'stf_spread', 0.1);
all_passed = check('a scalar spread is broadcast to every type', ...
    isequal(m.std_spread, [0.25 0.25]) && isequal(m.stf_spread, [0.1 0.1])) && all_passed;
bad = {{'std_spread', -0.1}, {'std_spread', [0.1 0.2 0.3]}, {'stf_spread', NaN}, {'stf_spread', [0.1 Inf]}, ...
    {'std_seed', [1 2]}, {'stf_seed', NaN}, {'std_seed', 'a'}};
ok = true;
for k = 1:numel(bad)
    ok = ok && throws_id(@() net(sc_chain, bad{k}{:}), 'SRNNCellTypePairs:InvalidParams');
end
all_passed = check('a negative, non-finite or wrong-length spread, and a seed that is not a finite scalar, error InvalidParams', ok) && all_passed;

pm = m.get_params();
all_passed = check('the nominal params fields stay nominal under a spread', ...
    isequal(pm.tau_b_rec{1, 1}, [2 5]) && isequal(pm.tau_b_rel{1, 1}, 0.25) && pm.std_n2(1, 1) == 2 && ...
    isequal(pm.std_spread, [0.25 0.25]) && isequal(pm.stf_spread, [0.1 0.1]) && ...
    isequal(pm.route{1, 1}.P.nominal.tau_rec, [2 5])) && all_passed;

m.run();
b = m.plot_data.b.E.E;
D = SRNNCellTypePairs.std_route_factor(b, pm, 1, 1);
plain = net(sc_chain);
plain.run();
p_plain = plain.get_params();
ok = isequal(D, b(:, 1, :)) && ...
    throws_id(@() SRNNCellTypePairs.std_route_factor(b(1:3, :, :), pm, 1, 1), 'SRNNCellTypePairs:RouteFactorRows') && ...
    isequal(SRNNCellTypePairs.std_route_factor(plain.plot_data.b.E.E(1:3, :, :), p_plain, 1, 1), plain.plot_data.b.E.E(1:3, 1, :));
gf = mcf.plot_data.g.E.E;
ok = ok && isequal(SRNNCellTypePairs.stf_route_factor(gf, pcf, 1, 1), gf(:, 1, :)) && ...
    throws_id(@() SRNNCellTypePairs.stf_route_factor(gf(1:3, :, :), pcf, 1, 1), 'SRNNCellTypePairs:RouteFactorRows');
all_passed = check(['std_route_factor and stf_route_factor: the whole route''s array is served; fewer rows error ' ...
    'RouteFactorRows on a spread route and are served on a route without a spread'], ok) && all_passed;
pf = floored.get_params();
bf = floored.plot_data.b.E.E;
k = pf.route{1, 1}.P.tau_rest ./ pf.route{1, 1}.P.tau_drive;
p_min = prod(1 ./ (1 + k), 2) .^ (1 / 2);                    % geometric mean of two, at r = 1
D = SRNNCellTypePairs.std_route_factor(bf, pf, 1, 1);
want = (prod(bf, 2) .^ (1 / 2) - p_min) ./ (1 - p_min);
all_passed = check(sprintf('the zero floor on a spread product route uses each neuron''s own p_min (max %.1e; p_min spans %.3f to %.3f)', ...
    max(abs(D - want), [], 'all'), min(p_min), max(p_min)), ...
    max(abs(D - want), [], 'all') <= 1e-13 && max(p_min) > 1.2 * min(p_min) && ...
    max(abs(SRNNCellTypePairs.std_zero_floor_pmin(pf, 1, 1) - p_min)) <= 1e-15) && all_passed;

m.n = 80;
stale = throws_id(@() m.get_params(), 'SRNNCellTypePairs:StaleSpreadDraws');
m.input_config.intrinsic_drive = [];                        % (sized for the old n by the first build)
m.build();
all_passed = check('a model whose n changed after build() errors StaleSpreadDraws from get_params until it is built again', ...
    stale && isequal(size(m.std_draws), [80 5]) && no_error(@() m.get_params())) && all_passed;

mc_class = meta.class.fromName('SRNNCellTypePairs');
props = mc_class.PropertyList;
ok = true;
for name = {'std_spread', 'stf_spread', 'std_seed', 'stf_seed'}
    pr = props(strcmp({props.Name}, name{1}));
    ok = ok && isscalar(pr) && strcmp(pr.SetAccess, 'public') && ~pr.Dependent;
end
for name = {'std_draws', 'stf_draws', 'std_clipped'}
    pr = props(strcmp({props.Name}, name{1}));
    ok = ok && isscalar(pr) && ~strcmp(pr.SetAccess, 'public');
end
all_passed = check('the spreads and seeds are ordinary public properties (model defaults of a sweep); the draws are read-only', ok) && all_passed;

fprintf('\n');
if all_passed
    fprintf('=== ALL std_spread TESTS PASSED ===\n');
else
    fprintf('=== SOME std_spread TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function sc = four_routes(route)
% One route configuration on all four routes of an E/I network.
sc = struct();
sc.E.E = route; sc.E.I = route; sc.I.E = route; sc.I.I = route;
end

function m = net(sc, varargin)
% A small E/I network with SFA on both types, built.
m = SRNNCellTypePairs('n', 60, 'indegree', 20, 'n_cellTypes', 2, 'cell_type_names', {'E', 'I'}, ...
    'f', [0.5 0.5], 'mu_tilde_relative', [3 -4], 'sigma_tilde_relative', [1 1], ...
    'tau_a', {[0.25 1 10], [0.25 1 10]}, 'c', [0.15 0.15], 'synapse_config', sc, ...
    'T_range', [0 3], 'fs', 200, 'ode_solver', 'rk4', 'lya_method', 'none', ...
    'store_full_state', true, 'verbose', 'near-none', varargin{:});
m.input_config.intrinsic_drive = [];                        % (sized for the old n by the first build)
m.build();
end

function m = pool_variants_network(extra)
% Four cell types with every variant of the pool model that has
% facilitation, beside product routes (the fmixed network of the tolerance
% tier's stf_chain suite), under a step staircase.
chain = @(tau_rec, tau_rel, n2) struct('model', 'chain', 'tau_rec', tau_rec, 'tau_rel', tau_rel, 'n2', n2);
fac = @(tau_dec, tau_fac, G) struct('model', 'chain', 'tau_dec', tau_dec, 'tau_fac', tau_fac, 'G', G);
fac_h = @(tau_dec, tau_fac, G, tau_h) struct('model', 'chain', 'tau_dec', tau_dec, 'tau_fac', tau_fac, ...
    'G', G, 'tau_h', tau_h);
sc = struct();
sc.E.E.std     = chain([2 5], 0.25, 2);
sc.E.E.stf     = fac_h(1.5, 0.25, 5, 0.5);
sc.E.PV.std    = struct('tau_rec', 2, 'tau_rel', 0.25);
sc.E.PV.stf    = fac(1.5, 0.25, 5);
sc.E.SST.stf   = fac_h(1.5, 0.25, 5, 0.5);
sc.E.VIP.std   = struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]);
sc.E.VIP.stf   = struct('tau_dec', [0.4 1.2], 'tau_fac', [0.2 0.3], 'G', [1.5 2]);
sc.PV.E.std    = chain([1 3], 0.3, 1.5);
sc.PV.E.stf    = fac(0.9, 0.4, 2.86);
sc.PV.E.scale  = 1.3;
sc.PV.PV.std   = struct('tau_rec', 0.6, 'tau_rel', 0.3);
sc.PV.PV.stf   = fac_h(0.9, 0.4, 2.86, 0.2);
sc.PV.SST.stf  = fac(0.9, 0.4, 2.86);
sc.SST.E.std   = chain([2 5], 0.25, 2);
sc.SST.PV.stf  = struct('tau_dec', 0.9, 'tau_fac', 0.4, 'G', 1.8);
sc.SST.SST.std = struct('tau_rec', 0.8, 'tau_rel', 0.25);
sc.SST.SST.stf = fac(2, 1, 1);
sc.VIP.E.std   = struct('tau_rec', [0.3 1], 'tau_rel', 0.2, 'combine', 'geomean');
sc.VIP.SST.std = chain([0.5 8], 0.2, 0.5);
sc.VIP.SST.stf = fac_h(3, 0.1, 8, 0.05);
sc.VIP.VIP.std = chain([1.5 4], 0.4, 1e6);
sc.VIP.VIP.stf = fac(0.5, 2, 1.5);
m = SRNNCellTypePairs('n', 40, 'indegree', 10, 'n_cellTypes', 4, 'cell_type_names', {'E', 'PV', 'SST', 'VIP'}, ...
    'f', [0.4 0.2 0.2 0.2], 'mu_tilde_relative', [3 -4 -2 -1.5], 'sigma_tilde_relative', [1 1 1 1], ...
    'tau_a', {[0.25 1], [], 2, 1}, 'c', [0.3 0 0.2 0.1], 'synapse_config', sc, ...
    'input_config', struct('n_steps', 6, 'step_density', struct(), 'amp', 0.5, ...
        'no_stim_pattern', false(1, 6), 'intrinsic_drive', 0, 'positive_only', false), ...
    'rng_seeds', [5 6], 'T_range', [0 3], 'fs', 200, 'ode_solver', 'rk4', 'lya_method', 'none', ...
    'store_full_state', true, 'verbose', 'near-none', extra{:});
m.input_config.intrinsic_drive = [];                        % (sized for the old n by the first build)
m.build();
end

function m = floor_network(extra)
% Three cell types with product routes only: the geometric mean, one and
% two depression components, facilitation, and the zero floor on.
sc = struct();
sc.E.E.std   = struct('tau_rec', [0.3 1], 'tau_rel', 0.2, 'combine', 'geomean');
sc.E.PV.stf  = struct('tau_dec', 0.5, 'tau_fac', 0.25, 'G', 2);
sc.E.SST.std = struct('tau_rec', 0.8, 'tau_rel', 0.25);
sc.E.SST.stf = struct('tau_dec', [0.4 1.2], 'tau_fac', [0.2 0.3], 'G', [1.5 2]);
sc.PV.PV.std = struct('tau_rec', [0.6 2], 'tau_rel', [0.3 0.4]);
sc.SST.PV.stf = struct('tau_dec', 0.9, 'tau_fac', 0.4, 'G', 1.8);
m = SRNNCellTypePairs('n', 30, 'indegree', 10, 'n_cellTypes', 3, 'cell_type_names', {'E', 'PV', 'SST'}, ...
    'f', [0.5 0.3 0.2], 'mu_tilde_relative', [3 -4 -2], 'sigma_tilde_relative', [1 1 1], ...
    'tau_a', {[0.25 1], [], 2}, 'c', [0.3 0 0.2], 'synapse_config', sc, 'std_zero_floor', true, ...
    'input_config', struct('n_steps', 6, 'step_density', struct(), 'amp', 0.5, ...
        'no_stim_pattern', false(1, 6), 'intrinsic_drive', 0, 'positive_only', false), ...
    'rng_seeds', [5 6], 'T_range', [0 3], 'fs', 200, 'ode_solver', 'rk4', 'lya_method', 'none', ...
    'store_full_state', true, 'verbose', 'near-none', extra{:});
m.input_config.intrinsic_drive = [];                        % (sized for the old n by the first build)
m.build();
end

function s = no_handles(s)
% A copy of a struct or cell with every function handle replaced by its
% text, so that two params structs compare with isequaln (two anonymous
% handles made separately are never isequal).
if isstruct(s)
    names = fieldnames(s);
    for k = 1:numel(s)
        for j = 1:numel(names)
            s(k).(names{j}) = no_handles(s(k).(names{j}));
        end
    end
elseif iscell(s)
    s = cellfun(@no_handles, s, 'UniformOutput', false);
elseif isa(s, 'function_handle')
    s = func2str(s);
end
end

function tf = no_error(fn)
try
    fn();
    tf = true;
catch
    tf = false;
end
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
