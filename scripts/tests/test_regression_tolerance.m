% test_regression_tolerance.m
%
% The tolerance regression tier, layers 1 and 2 (about 20 s, in the client;
% no pool is needed): the right-hand side, both Jacobians and the packed
% physics evaluated at stored states, and small stable networks end to end
% (trajectories, three Lyapunov estimators, transient gain, a sweep, memory
% capacity), each compared ELEMENT BY ELEMENT, within a tolerance, with the
% anchor scripts/tests/fixtures/regression_tol_anchor.mat. The anchor was
% captured from the model class as it stood before Plan A touched it.
%
% This is the gate for a rewrite of the model's arithmetic, where the
% byte-exact tiers must differ (a reordered expression changes the last
% bit). It passes when nothing moved by more than round-off; the printout
% says how close the worst leaf of each layer was.
%
% IF IT FAILS, the change moved a number by more than round-off: find the
% defect. Do not loosen a tolerance (regression_tol_settings) and do not
% re-capture the anchor. To bisect a change, capture a pre-edit reference
% with make_regression_tol_reference('<name>') and compare with
% regression_tol_check('<name>'). Layer 3 (the chaotic tiers, structurally)
% is regression_tol_band_check, run on demand. See docs/ISSUES.md, FD-020,
% and docs/Plan_A_chained_vessicles.md, sections 2.2 and 6.

[ok, leaves] = regression_tol_check('anchor');
assert(ok, 'test_regression_tolerance:Differences', ...
    '%d leaf(s) of the tolerance tier''s anchor are outside tolerance or missing; see the list above.', ...
    nnz(ismember(leaves.status, {'FAIL', 'missing'})));
fprintf('PASS  tolerance regression tier matches the anchor within tolerance\n');
