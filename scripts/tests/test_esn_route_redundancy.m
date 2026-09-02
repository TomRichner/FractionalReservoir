% TEST_ESN_ROUTE_REDUNDANCY The 'synaptic' MC readout reads one route, safely.
%
% SRNNCellTypePairs gives a presynaptic neuron ONE SYNAPTIC OUTPUT PER TARGET
% TYPE: theta_j^{q->s} = r_j * prod(b^{q->s}) * prod(g^{q->s}). So "the synaptic
% output of neuron j" is only defined when all of j's outgoing routes share
% their synaptic dynamics, and SRNN_ESN_reservoir's 'synaptic' readout takes
% route 1 on exactly that assumption.
%
% THE CONFIG CHECK RESTS ON A NUMERICAL CLAIM, WHICH IS WHAT THIS TEST PROVES.
% assert_route_redundancy compares (tau_rec, tau_rel, tau_dec, tau_fac, G) and
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

function nv = namevalue(s)
f = fieldnames(s);
nv = cell(1, 2*numel(f));
for i = 1:numel(f); nv{2*i-1} = f{i}; nv{2*i} = s.(f{i}); end
end
