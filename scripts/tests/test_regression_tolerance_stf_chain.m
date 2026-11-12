% test_regression_tolerance_stf_chain.m
%
% The tolerance regression tier's STF_CHAIN suite (about 30 s, in the
% client; no pool is needed): networks with release-coupled facilitation
% (stf.model = 'chain', route_pool: a release probability g that
% multiplies the output and the rate at which the releasable pool is used,
% optionally driven by a calcium trace h), compared ELEMENT BY ELEMENT,
% within the tier's tolerances, with
% scripts/tests/fixtures/regression_tol_stf_chain.mat.
%
%   - pointwise: the right-hand side, the assembled Jacobian, the
%     matrix-free Jacobian product, the parameter view and every route's
%     depression and facilitation factor at stored states of the
%     regression network with four routes that carry a chain, g and h
%     (noise off and on, with and without the step staircase), of its
%     stable twin, and of a four-cell-type network with every variant of
%     the pool model that has facilitation (no depression, a single pool
%     or the chain; with and without the calcium trace) beside product
%     routes; among the states, pools depleted with g near its ceiling,
%     and a rested route with a nearly full reserve;
%   - stable networks end to end: trajectories, the three Lyapunov
%     estimators, transient gain, the finite-difference Jacobian, a sweep,
%     and memory capacity with the synaptic readout on four identical
%     coupled routes.
%
% Neither the anchor (test_regression_tolerance) nor the chain reference
% (test_regression_tolerance_chain) holds a coupled route, so this is the
% reference that says "coupled facilitation is unchanged". It was captured
% when coupled facilitation was first built (FD-018), before the
% per-neuron spread was built on top of it, and is the gate for that work
% and for any later edit to route_pool's facilitation.
%
% IF IT FAILS, the change moved a number by more than round-off: find the
% defect. Do not loosen a tolerance (regression_tol_settings). The
% reference is replaced only deliberately (make_regression_tol_reference's
% header says how), with the reason in the commit message. To bisect,
% capture a pre-edit reference with
% make_regression_tol_reference('<name>', 'suite', 'stf_chain').
% See docs/Plan_A_chained_vessicles.md, sections 4.6 and 7.3 (step 8), and
% docs/notes/Tolerance_tier_mutation_checks_2026-11-11.md for the
% deliberate errors this suite was shown to catch.

[ok, leaves] = regression_tol_check('stf_chain');
assert(ok, 'test_regression_tolerance_stf_chain:Differences', ...
    '%d leaf(s) of the tolerance tier''s stf_chain reference are outside tolerance or missing; see the list above.', ...
    nnz(ismember(leaves.status, {'FAIL', 'missing'})));
fprintf('PASS  tolerance regression tier, stf_chain suite, matches its reference within tolerance\n');
