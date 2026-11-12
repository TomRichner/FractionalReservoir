% test_regression_tolerance_chain.m
%
% The tolerance regression tier's CHAIN suite (about 30 s, in the client;
% no pool is needed): networks with chained vesicle-pool depression
% (std.model = 'chain', route_pool), compared ELEMENT BY ELEMENT, within
% the tier's tolerances, with scripts/tests/fixtures/regression_tol_chain.mat.
%
%   - pointwise: the right-hand side, the assembled Jacobian, the
%     matrix-free Jacobian product and the parameter view at stored states
%     of the regression network with four chain routes (noise off and on,
%     with and without the step staircase), of its stable twin, and of a
%     four-cell-type network in which one presynaptic type has a chain, a
%     two-component product, a one-component and an undepressed route;
%     among the states, both pools deeply depleted and a nearly full
%     reserve;
%   - stable networks end to end: trajectories, the three Lyapunov
%     estimators, transient gain, the finite-difference Jacobian, a sweep,
%     and memory capacity with the synaptic readout on four identical
%     chain routes.
%
% The anchor (test_regression_tolerance) holds no chain route, so this is
% the reference that says "the chain is unchanged". It was captured when
% the chain was first merged, before chained facilitation was built into
% route_pool, and is the gate for that work and for any later edit to the
% chain's arithmetic.
%
% IF IT FAILS, the change moved a chain number by more than round-off:
% find the defect. Do not loosen a tolerance (regression_tol_settings). The
% reference is replaced only deliberately (make_regression_tol_reference's
% header says how), with the reason in the commit message. To bisect,
% capture a pre-edit reference with
% make_regression_tol_reference('<name>', 'suite', 'chain').
% See docs/Plan_A_chained_vessicles.md, section 7.3, step 5, and
% docs/notes/Tolerance_tier_mutation_checks_2026-11-11.md for the
% deliberate errors this suite was shown to catch.

[ok, leaves] = regression_tol_check('chain');
assert(ok, 'test_regression_tolerance_chain:Differences', ...
    '%d leaf(s) of the tolerance tier''s chain reference are outside tolerance or missing; see the list above.', ...
    nnz(ismember(leaves.status, {'FAIL', 'missing'})));
fprintf('PASS  tolerance regression tier, chain suite, matches its reference within tolerance\n');
