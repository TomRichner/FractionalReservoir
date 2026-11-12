% test_regression_tolerance_spread.m
%
% The tolerance regression tier's SPREAD suite (about 30 s, in the client;
% no pool is needed): networks whose depression and facilitation
% parameters differ per neuron (std_spread, stf_spread; FD-015), compared
% ELEMENT BY ELEMENT, within the tier's tolerances, with
% scripts/tests/fixtures/regression_tol_spread.mat.
%
%   - pointwise: the standard normals drawn for each neuron, the value
%     each neuron has of every route parameter, which neurons the chain's
%     ordering rule moved, and with those values the right-hand side, the
%     assembled Jacobian, the matrix-free Jacobian product and every
%     route's depression and facilitation factor at stored states, on the
%     regression network with four spread chain routes (noise off and on,
%     with and without the step staircase), on its stable twin, on a
%     four-cell-type network with a spread on every kind of route (chain,
%     single pool, one and two product components, the geometric mean,
%     release-coupled and output-only facilitation) where the spreads
%     differ between the cell types and one type per mechanism has none,
%     and on a product network with the zero floor;
%   - stable networks end to end: trajectories, the three Lyapunov
%     estimators, transient gain, the finite-difference Jacobian, a sweep
%     with the spread as a model default, and memory capacity with the
%     synaptic readout on four identical spread chain routes.
%
% The other three references hold no spread: they say "with no spread
% nothing moved", this one says "the spread is unchanged". It was captured
% when the spread was first built.
%
% IF IT FAILS, the change moved a number by more than round-off: find the
% defect. Do not loosen a tolerance (regression_tol_settings). The
% reference is replaced only deliberately (make_regression_tol_reference's
% header says how), with the reason in the commit message. To bisect,
% capture a pre-edit reference with
% make_regression_tol_reference('<name>', 'suite', 'spread').
% See docs/Plan_A_chained_vessicles.md, sections 4.3 and 4.6, and
% docs/notes/Tolerance_tier_mutation_checks_2026-11-11.md for the
% deliberate errors this suite was shown to catch.

[ok, leaves] = regression_tol_check('spread');
assert(ok, 'test_regression_tolerance_spread:Differences', ...
    '%d leaf(s) of the tolerance tier''s spread reference are outside tolerance or missing; see the list above.', ...
    nnz(ismember(leaves.status, {'FAIL', 'missing'})));
fprintf('PASS  tolerance regression tier, spread suite, matches its reference within tolerance\n');
