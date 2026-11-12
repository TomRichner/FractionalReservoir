% TEST_ESN_ROUTE_REDUNDANCY The 'synaptic' MC readout reads one route, safely.
%
% SRNNCellTypePairs gives a presynaptic neuron ONE SYNAPTIC OUTPUT PER TARGET
% TYPE: theta_j^{q->s} = r_j * prod(b^{q->s}) * prod(g^{q->s}). So "the synaptic
% output of neuron j" is only defined when all of j's outgoing routes share
% their synaptic dynamics, and SRNN_ESN_reservoir's 'synaptic' readout takes
% route 1 on exactly that assumption.
%
% THE CONFIG CHECK RESTS ON A NUMERICAL CLAIM, WHICH IS WHAT THIS TEST PROVES.
% assert_route_redundancy compares (tau_rec, tau_rel, tau_dec, tau_fac, G) and,
% since 2026-11-12, the route model and every parameter of it (sections 5 to 9:
% chain routes, whose output is r * b1), and
% concludes the trajectories coincide. That follows because
% db/dt = (1-b)/tau_rec - b*r_j/tau_rel depends only on the route's constants
% and the presynaptic neuron's own rate, from a common b = 1 -- identical ODE,
% identical driver, identical initial condition. Section 1 checks it by running
% the model rather than by trusting the argument.
%
% PER PRESYNAPTIC TYPE, NOT GLOBALLY. Neuron j of type q has only the routes
% q->1 ... q->C; what other types do cannot reach it. Section 3 pins that down:
% "E->* depressed, I->* not" is well defined and must be ACCEPTED, while STD on
% E->E alone must be refused.
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths has run.
%
% See also: SRNN_ESN_reservoir, run_memory_capacity, SRNNCellTypePairs

all_passed = true;
fprintf('\n=== test_esn_route_redundancy ===\n');

dual = struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]);
base = {'n', 60, 'indegree', 20, 'n_cellTypes', 2, ...
    'cell_type_names', {'E', 'I'}, 'f', [0.6 0.4], ...
    'mu_tilde_relative', [3 -4; 3 -4], 'sigma_tilde_relative', ones(2, 2), ...
    'level_of_chaos', 1.5, 'S_c', 0.35, 'fs', 100, 'ode_solver', 'rk4', ...
    'T_wash', 20, 'T_train', 60, 'T_test', 40, 'd_max', 8, ...
    'input_type', 'sample_hold', 'T_hold', 0.3};

%% ------------------------------------------------------------------------
%  1. Identical routes really do give bit-identical synaptic output.
%% ------------------------------------------------------------------------
fprintf('\n--- 1. identical routes are bit-identical ---\n');
sc_all = struct();
sc_all.E.E.std = dual;  sc_all.E.I.std = dual;
sc_all.I.E.std = dual;  sc_all.I.I.std = dual;

esn = SRNN_ESN_reservoir(base{:}, 'tau_a', {[0.25 1.58 10], []}, ...
    'synapse_config', sc_all);
evalc('esn.build(); esn.run_reservoir_esn()');
[~, ~, ~, ~, syn] = SRNNCellTypePairs.unpack_and_compute_states( ...
    esn.S_out, esn.cached_params);

d_E = max(abs(syn.E.E(:) - syn.E.I(:)));
d_I = max(abs(syn.I.E(:) - syn.I.I(:)));
all_passed = report(d_E == 0, sprintf('max|syn.E.E - syn.E.I| = %g (want exactly 0)', d_E)) && all_passed;
all_passed = report(d_I == 0, sprintf('max|syn.I.E - syn.I.I| = %g (want exactly 0)', d_I)) && all_passed;

% ...and the depression is real, so this is not passing on a trivial b == 1.
all_passed = report(min(syn.E.E(:)) < 0.9 * max(syn.E.E(:)), ...
    'depression is actually engaged (not a trivial b == 1 match)') && all_passed;

%% ------------------------------------------------------------------------
%  2. The check accepts identical routes and refuses mixed ones.
%% ------------------------------------------------------------------------
fprintf('\n--- 2. assert_route_redundancy ---\n');
all_passed = report(no_error(@() SRNN_ESN_reservoir.assert_route_redundancy( ...
    esn.cached_params)), 'accepts STD identical on all four routes') && all_passed;

sc_one = struct(); sc_one.E.E.std = dual;      % E->E only: E neurons ambiguous
esn_bad = SRNN_ESN_reservoir(base{:}, 'tau_a', {[], []}, 'synapse_config', sc_one);
evalc('esn_bad.build()');
all_passed = report(errors_with(@() SRNN_ESN_reservoir.assert_route_redundancy( ...
    esn_bad.cached_params), 'SRNN_ESN_reservoir:AmbiguousSynapticOutput'), ...
    'refuses STD on E->E only') && all_passed;

% Different TIME CONSTANTS on two routes is just as ambiguous as one missing.
sc_diff = struct();
sc_diff.E.E.std = dual;
sc_diff.E.I.std = struct('tau_rec', [1 3], 'tau_rel', [0.25 0.5]);
esn_diff = SRNN_ESN_reservoir(base{:}, 'tau_a', {[], []}, 'synapse_config', sc_diff);
evalc('esn_diff.build()');
all_passed = report(errors_with(@() SRNN_ESN_reservoir.assert_route_redundancy( ...
    esn_diff.cached_params), 'SRNN_ESN_reservoir:AmbiguousSynapticOutput'), ...
    'refuses two routes with different tau_rec') && all_passed;

% STF must be covered too: synaptic output is r * prod(b) * prod(g).
sc_stf = struct();
sc_stf.E.E.std = dual;  sc_stf.E.I.std = dual;
sc_stf.I.E.std = dual;  sc_stf.I.I.std = dual;
sc_stf.E.E.stf = struct('tau_dec', 1.0, 'tau_fac', 0.5, 'G', 2.0);
esn_stf = SRNN_ESN_reservoir(base{:}, 'tau_a', {[], []}, 'synapse_config', sc_stf);
evalc('esn_stf.build()');
all_passed = report(errors_with(@() SRNN_ESN_reservoir.assert_route_redundancy( ...
    esn_stf.cached_params), 'SRNN_ESN_reservoir:AmbiguousSynapticOutput'), ...
    'refuses STF on one route only (not just STD)') && all_passed;

%% ------------------------------------------------------------------------
%  3. The check is PER PRESYNAPTIC TYPE, so a per-type difference is fine.
%% ------------------------------------------------------------------------
fprintf('\n--- 3. per presynaptic type, not global ---\n');
% E->* depressed, I->* not. Every neuron still has ONE synaptic output, so this
% must be ACCEPTED -- a global "all routes identical" check would reject it.
sc_by_pre = struct();
sc_by_pre.E.E.std = dual;
sc_by_pre.E.I.std = dual;
esn_pre = SRNN_ESN_reservoir(base{:}, 'tau_a', {[], []}, 'synapse_config', sc_by_pre);
evalc('esn_pre.build()');
all_passed = report(no_error(@() SRNN_ESN_reservoir.assert_route_redundancy( ...
    esn_pre.cached_params)), 'accepts E->* depressed with I->* undepressed') && all_passed;

%% ------------------------------------------------------------------------
%  4. The paper's MC preset must be readable with 'synaptic'.
%% ------------------------------------------------------------------------
fprintf('\n--- 4. mc_pairs_dualStd is legal for the synaptic readout ---\n');
[preset, ~, conds] = srnn_param_preset('mc_pairs_dualStd');
ok_all = true;
for k = 1:numel(conds)
    c = rmfield(conds{k}, 'name');
    args = [namevalue(preset), namevalue(c), ...
        {'n', 40, 'indegree', 10, 'fs', 100, 'T_wash', 10, 'T_train', 30, ...
         'T_test', 20, 'd_max', 4, 'input_type', 'sample_hold', 'T_hold', 0.3}];
    m = SRNN_ESN_reservoir(args{:});
    evalc('m.build()');
    ok_all = ok_all && no_error(@() SRNN_ESN_reservoir.assert_route_redundancy(m.cached_params));
end
all_passed = report(ok_all, 'all 4 conditions pass the redundancy check') && all_passed;

%% ------------------------------------------------------------------------
%  5. Chain routes (std.model = 'chain', route_pool; FD-006).
%
%  A chain route carries a releasable pool b1 and a reserve b2, and its
%  synaptic output is r * b1 -- NOT r * b1 * b2. The reservoir class was
%  written when depression was always a product; these sections pin what it
%  does with a chain.
%% ------------------------------------------------------------------------
fprintf('\n--- 5. identical chain routes: the readout is rate times b1 ---\n');
chain = struct('model', 'chain', 'tau_rec', [2 4], 'tau_rel', 0.25, 'n2', 2);
sc_chain = struct();
sc_chain.E.E.std = chain;  sc_chain.E.I.std = chain;
sc_chain.I.E.std = chain;  sc_chain.I.I.std = chain;

esn_c = SRNN_ESN_reservoir(base{:}, 'tau_a', {[0.25 1.58 10], []}, ...
    'synapse_config', sc_chain);
evalc('esn_c.build(); esn_c.run_reservoir_esn()');
pc = esn_c.cached_params;
[~, ~, b_c, r_c, syn_c] = SRNNCellTypePairs.unpack_and_compute_states(esn_c.S_out, pc);

all_passed = report(no_error(@() SRNN_ESN_reservoir.assert_route_redundancy(pc)), ...
    'accepts identical chains on all four routes') && all_passed;
b1_E = reshape(b_c.E.E(:, 1, :), size(r_c.E));
b2_E = reshape(b_c.E.E(:, 2, :), size(r_c.E));
b1_I = reshape(b_c.I.E(:, 1, :), size(r_c.I));
all_passed = report(isequal(syn_c.E.E, r_c.E .* b1_E) && isequal(syn_c.I.E, r_c.I .* b1_I), ...
    'synaptic output equals rate .* b1 exactly, on E and on I') && all_passed;
all_passed = report(isequal(syn_c.E.E, syn_c.E.I) && isequal(syn_c.I.E, syn_c.I.I), ...
    'and is bit-identical across each neuron''s two routes') && all_passed;
% ...and the reserve is really in play, so r*b1 and r*b1*b2 are different
% arrays and the check above could have failed.
all_passed = report(min(b2_E(:)) < 0.999 && ~isequal(syn_c.E.E, r_c.E .* b1_E .* b2_E) ...
    && min(b1_E(:)) < 0.9, ...
    sprintf('the reserve is engaged (min b2 = %.3f), so r*b1 is not r*b1*b2', min(b2_E(:)))) && all_passed;

%% ------------------------------------------------------------------------
%  6. The check sees the route MODEL and every parameter of it, and its
%     message shows the difference.
%% ------------------------------------------------------------------------
fprintf('\n--- 6. chain routes that differ are refused, and the message says how ---\n');
amb = 'SRNN_ESN_reservoir:AmbiguousSynapticOutput';
plain = {'tau_a', {[], []}};

% a. A chain and a two-component product with the same time constants.
sc_mix = sc_chain;
sc_mix.E.I.std = struct('tau_rec', [2 4], 'tau_rel', 0.25);
msg = refusal(base, plain, sc_mix);
all_passed = report(strcmp(msg.id, amb), ...
    'refuses E->E chain beside E->I two-component product, same time constants') && all_passed;
all_passed = report(contains(msg.text, 'chained STD(') && contains(msg.text, 'product STD(') ...
    && contains(msg.text, 'model: ''pool'' vs ''product'''), ...
    'the message describes one route as a chain, the other as a product, and names the model') && all_passed;
all_passed = report(contains(msg.text, 'releasable-pool refill tau_rec1=2,') ...
    && contains(msg.text, 'reserve refill tau_rec2=4,') ...
    && contains(msg.text, 'release tau_rel=0.25,') && contains(msg.text, 'n2=2)'), ...
    'a chain is described as a chain: both refill times, the release time, n2') && all_passed;

% b. Two chains that differ ONLY in n2, by 1e-6 relative. Before 2026-11-12
%    the check did not look at n2 at all and this PASSED.
sc_n2 = sc_chain;
sc_n2.E.I.std.n2 = 2 * (1 + 1e-6);
msg = refusal(base, plain, sc_n2);
all_passed = report(strcmp(msg.id, amb), 'refuses two chains that differ only in n2 (1e-6 relative)') && all_passed;
all_passed = report(shows_difference(msg.text) && contains(msg.text, 'n2: 2 vs 2.000002') ...
    && contains(msg.text, 'relative difference 1e-06'), ...
    'the two route lines differ, and the message names n2 with both values') && all_passed;

% c. ...or only in the reserve's refill time.
sc_rec = sc_chain;
sc_rec.E.I.std.tau_rec = [2, 4 * (1 + 1e-6)];
msg = refusal(base, plain, sc_rec);
all_passed = report(strcmp(msg.id, amb), ...
    'refuses two chains that differ only in the reserve refill time (1e-6 relative)') && all_passed;
all_passed = report(shows_difference(msg.text) && contains(msg.text, 'tau_rec: [2 4] vs [2 4.000004]'), ...
    'the two route lines differ, and the message names tau_rec with both values') && all_passed;

% d. The same holds for PRODUCT routes: the old message rounded to three
%    digits and printed two identical lines under "non-identical".
sc_close = struct();
sc_close.E.E.std = dual;
sc_close.E.I.std = struct('tau_rec', [2 4], 'tau_rel', [0.25, 0.5 * (1 + 1e-6)]);
msg = refusal(base, plain, sc_close);
all_passed = report(strcmp(msg.id, amb) && shows_difference(msg.text) ...
    && contains(msg.text, 'tau_rel: [0.25 0.5] vs [0.25 0.5000005]'), ...
    'two product routes 1e-6 apart: refused, and the message shows it') && all_passed;

% e. Per presynaptic type still: E->* a chain, I->* a product.
fprintf('\n--- 7. E->* chain with I->* product is accepted ---\n');
sc_types = struct();
sc_types.E.E.std = chain;  sc_types.E.I.std = chain;
sc_types.I.E.std = dual;   sc_types.I.I.std = dual;
esn_t = SRNN_ESN_reservoir(base{:}, plain{:}, 'synapse_config', sc_types);
evalc('esn_t.build()');
all_passed = report(no_error(@() SRNN_ESN_reservoir.assert_route_redundancy( ...
    esn_t.cached_params)), 'accepts E->* chain with I->* two-component product') && all_passed;

% f. EVERY VARIANT of the pool model is described, not only the chain. A
%    pool route is also a single pool or has no depression at all (beside
%    release-coupled facilitation, stf.model = 'chain'). describe_route
%    read tau_rec(2) and n2 on every pool route, so two such routes that
%    differed made the check throw MATLAB:badsubscript where it must raise
%    AmbiguousSynapticOutput (found 2026-11-12, fixed with FD-015).
fprintf('\n--- 7b. non-identical pool routes with a single pool or no depression ---\n');
fac = struct('model', 'chain', 'tau_dec', 1.5, 'tau_fac', 0.25, 'G', 5);
fac_h = struct('model', 'chain', 'tau_dec', 1.5, 'tau_fac', 0.25, 'G', 5, 'tau_h', 0.5);
single = struct('tau_rec', 2, 'tau_rel', 0.25);
% single pool with g on both routes, G differing
sc_single = struct();
sc_single.E.E.std = single;  sc_single.E.E.stf = fac;
sc_single.E.I.std = single;  sc_single.E.I.stf = fac;
sc_single.E.I.stf.G = 4;
msg = refusal(base, plain, sc_single);
all_passed = report(strcmp(msg.id, amb) && shows_difference(msg.text) && ...
    contains(msg.text, 'single-pool STD(refill tau_rec=2, release tau_rel=0.25)') && ...
    contains(msg.text, 'release-coupled STF(tau_dec=1.5, tau_fac=0.25, G=5)') && contains(msg.text, 'G: 5 vs 4'), ...
    'two single-pool routes with coupled facilitation, G differing: refused, each described as a single pool') && all_passed;
% no depression: g alone against g with a calcium trace
sc_nodep = struct();
sc_nodep.E.E.stf = fac;
sc_nodep.E.I.stf = fac_h;
msg = refusal(base, plain, sc_nodep);
all_passed = report(strcmp(msg.id, amb) && shows_difference(msg.text) && ...
    contains(msg.text, 'calcium trace tau_h=0.5') && ~contains(msg.text, 'STD('), ...
    'two facilitating routes without depression, one with a calcium trace: refused, the trace is named') && all_passed;
% a single pool against a chain, both with g
sc_mix = struct();
sc_mix.E.E.std = single;  sc_mix.E.E.stf = fac;
sc_mix.E.I.std = chain;   sc_mix.E.I.stf = fac;
msg = refusal(base, plain, sc_mix);
all_passed = report(strcmp(msg.id, amb) && shows_difference(msg.text) && ...
    contains(msg.text, 'single-pool STD(') && contains(msg.text, 'chained STD('), ...
    'a single pool against a chain, both facilitated: refused, one described as each') && all_passed;
% a facilitating route without depression against one with no synapse dynamics
sc_bare = struct();
sc_bare.E.I.stf = fac_h;
msg = refusal(base, plain, sc_bare);
all_passed = report(strcmp(msg.id, amb) && shows_difference(msg.text) && ...
    contains(msg.text, 'no short-term plasticity'), ...
    'a facilitating pool route against a route with no synapse dynamics: refused') && all_passed;
% and identical pool routes of each variant are accepted
sc_ok = struct();
sc_ok.E.E.std = single;  sc_ok.E.E.stf = fac_h;
sc_ok.E.I.std = single;  sc_ok.E.I.stf = fac_h;
sc_ok.I.E.stf = fac;     sc_ok.I.I.stf = fac;
msg = refusal(base, plain, sc_ok);
all_passed = report(isempty(msg.id), ...
    'identical single pools with g and h from E, identical g alone from I: accepted') && all_passed;

%% ------------------------------------------------------------------------
%  8. Memory capacity runs on chain networks.
%% ------------------------------------------------------------------------
fprintf('\n--- 8. memory capacity on chain networks ---\n');
mc_base = [base, {'T_wash', 300, 'T_train', 1500, 'T_test', 900, 'd_max', 30, ...
    'T_hold', 0.1, 'lya_method', 'none'}];        % later name-value pairs win

esn_mc = SRNN_ESN_reservoir(mc_base{:}, 'tau_a', {[0.25 1.58 10], []}, ...
    'synapse_config', sc_chain);
evalc('esn_mc.build(); [MC_syn, R2_syn] = esn_mc.run_memory_capacity(''readout_signal'', ''synaptic'');');
all_passed = report(isfinite(MC_syn) && all(isfinite(R2_syn)) && MC_syn > 0, ...
    sprintf('all-chain network, synaptic readout: MC = %.4f, every R^2 finite', MC_syn)) && all_passed;
evalc('[MC_rate, R2_rate] = esn_mc.run_memory_capacity(''readout_signal'', ''rate'');');
all_passed = report(isfinite(MC_rate) && all(isfinite(R2_rate)) && MC_rate ~= MC_syn, ...
    sprintf('same network, rate readout: MC = %.4f (a different signal, a different MC)', MC_rate)) && all_passed;

% Genuinely per-route synapses: every route of a kind of its own. The
% synaptic readout is refused, and 'rate' is the way through.
sc_kinds = struct();
sc_kinds.E.E.std = chain;
sc_kinds.E.I.std = struct('tau_rec', 2, 'tau_rel', 0.25);
sc_kinds.I.E.std = struct('model', 'chain', 'tau_rec', [1 3], 'tau_rel', 0.3, 'n2', 1.5);
esn_k = SRNN_ESN_reservoir(mc_base{:}, 'tau_a', {[0.25 1.58 10], []}, ...
    'synapse_config', sc_kinds);
evalc('esn_k.build()');
all_passed = report(errors_with(@() esn_k.run_memory_capacity('readout_signal', 'synaptic'), amb), ...
    'per-route synapses: the synaptic readout is refused by run_memory_capacity') && all_passed;
evalc('[MC_k, R2_k] = esn_k.run_memory_capacity(''readout_signal'', ''rate'');');
all_passed = report(isfinite(MC_k) && all(isfinite(R2_k)) && MC_k > 0, ...
    sprintf('per-route synapses: the rate readout runs (MC = %.4f)', MC_k)) && all_passed;

%% ------------------------------------------------------------------------
%  9. plot_esn_timeseries draws a route's FACTOR (std_route_factor): b1 on
%     a chain route, prod(b) on a plain product route as it always did.
%
%  The panel multiplies the factors of a type's routes together, as it did
%  before chain routes existed. In esn_k type I has ONE depressed route, a
%  chain (I->E), so its trace is that route's b1 itself; type E has a chain
%  and a one-component product, so its trace is b1 times the product's b.
%% ------------------------------------------------------------------------
fprintf('\n--- 9. the time-series plot draws the factor ---\n');
old_visibility = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
restore_visibility = onCleanup(@() set(groot, 'DefaultFigureVisible', old_visibility));

% a. One chain route on I: the trace is b1 of that route, exactly.
[drawn, bk, ti] = depression_lines(esn_k);
nI = size(bk.I.E, 1);  nE = size(bk.E.E, 1);
b1_IE = reshape(bk.I.E(:, 1, ti), nI, []);
b12_IE = reshape(prod(bk.I.E(:, :, ti), 2), nI, []);
b1_EE = reshape(bk.E.E(:, 1, ti), nE, []);
b_EI = reshape(bk.E.I(:, 1, ti), nE, []);
want = sortrows([b1_IE; b1_EE .* b_EI]);
all_passed = report(isequal(drawn, want) && ~isequal(drawn, sortrows([b12_IE; b1_EE .* b_EI])), ...
    'chain route: the trace is b1 (not b1*b2); beside a product route, b1 times that route''s b') && all_passed;

% b. The all-chain network builds too; its two identical routes per type
%    multiply, as two identical product routes always have in this panel.
[drawn, bc, ti] = depression_lines(esn_mc);
b1_all = [reshape(bc.E.E(:, 1, ti), size(bc.E.E, 1), []); reshape(bc.I.E(:, 1, ti), size(bc.I.E, 1), [])];
all_passed = report(isequal(drawn, sortrows(b1_all .* b1_all)), ...
    'all-chain network: the plot builds, and its trace is made of b1 alone (b1 per route, two routes)') && all_passed;

% c. A plain product network draws what it drew before: prod(b) within a
%    route, multiplied across the type's routes.
esn_p = SRNN_ESN_reservoir(mc_base{:}, 'tau_a', {[0.25 1.58 10], []}, 'synapse_config', sc_all);
evalc('esn_p.build(); esn_p.run_memory_capacity(''readout_signal'', ''synaptic'');');
[drawn, bp, ti] = depression_lines(esn_p);
old_E = reshape(prod(bp.E.E, 2), size(bp.E.E, 1), []) .* reshape(prod(bp.E.I, 2), size(bp.E.I, 1), []);
old_I = reshape(prod(bp.I.E, 2), size(bp.I.E, 1), []) .* reshape(prod(bp.I.I, 2), size(bp.I.I, 1), []);
all_passed = report(isequal(drawn, sortrows([old_E(:, ti); old_I(:, ti)])), ...
    'product network: the trace is bit-identical to the prod(b) collapse it drew before') && all_passed;
clear restore_visibility

%% ------------------------------------------------------------------------
fprintf('\n');
if all_passed
    fprintf('test_esn_route_redundancy: ALL TESTS PASSED\n');
else
    fprintf(2, 'test_esn_route_redundancy: FAILURES ABOVE\n');
end

%% ------------------------------------------------------------------------
function ok = report(ok, msg)
if ok
    fprintf('  PASS  %s\n', msg);
else
    fprintf(2, '  FAIL  %s\n', msg);
end
end

function tf = no_error(fn)
try; fn(); tf = true; catch; tf = false; end
end

function tf = errors_with(fn, id)
try
    fn();
    tf = false;
catch ME
    tf = strcmp(ME.identifier, id);
end
end

function msg = refusal(base, extra, sc)
% Build a reservoir with synapse_config sc and return what
% assert_route_redundancy raises: msg.id and msg.text ('' if it passes).
m = SRNN_ESN_reservoir(base{:}, extra{:}, 'synapse_config', sc);
evalc('m.build()');
msg = struct('id', '', 'text', '');
try
    SRNN_ESN_reservoir.assert_route_redundancy(m.cached_params);
catch ME
    msg.id = ME.identifier;
    msg.text = ME.message;
end
end

function tf = shows_difference(text)
% The message's two route lines ('  E->E : ...' and '  E->I : ...') must not
% describe the two routes identically.
lines = strsplit(text, newline);
routes = lines(contains(lines, '->') & contains(lines, ' : '));
described = extractAfter(routes, ' : ');
tf = numel(described) == 2 && ~strcmp(described{1}, described{2});
end

function [drawn, b, ti] = depression_lines(esn)
% The rows drawn in plot_esn_timeseries's depression panel (sorted), the
% stored b struct, and the test-period sample indices they cover.
[fig, ax] = esn.plot_esn_timeseries(1);
closer = onCleanup(@() close(fig));
labels = arrayfun(@(a) get(get(a, 'YLabel'), 'String'), ax, 'UniformOutput', false);
panel = ax(strcmp(labels, 'depression'));
lines = findobj(panel, 'Type', 'line');
drawn = sortrows(vertcat(lines.YData));
b = esn.mc_results.b;
names = fieldnames(b);
posts = fieldnames(b.(names{1}));
ti = (esn.T_wash + esn.T_train + 1):size(b.(names{1}).(posts{1}), 3);
end

function nv = namevalue(s)
f = fieldnames(s);
nv = cell(1, 2*numel(f));
for i = 1:numel(f); nv{2*i-1} = f{i}; nv{2*i} = s.(f{i}); end
end
