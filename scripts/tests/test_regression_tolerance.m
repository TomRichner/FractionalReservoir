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

%
% SUPERSEDED LEAVES. The anchor is never re-captured. Leaves that a
% deliberate change of physics moved are listed, with the reason, in
% regression_tol_settings().anchor_superseded and are not compared against
% the anchor (regression_tol_check applies the list to the anchor only).
% Since F3 (commit cce78c9, 2026-11-12) that is 176 leaves of sfa1_std1.
% This test checks the count is what the settings say it is, so a pattern
% that grows to swallow other leaves, or a leaf that newly fails inside
% sfa1_std1, is seen; test_regression_tol_superseded checks the mechanism.

[ok, leaves, info] = regression_tol_check('anchor');
assert(ok, 'test_regression_tolerance:Differences', ...
    '%d leaf(s) of the tolerance tier''s anchor are outside tolerance or missing; see the list above.', ...
    nnz(ismember(leaves.status, {'FAIL', 'missing'})));
n_superseded = nnz(strcmp(leaves.status, 'superseded'));
assert(n_superseded == 176, 'test_regression_tolerance:Superseded', ...
    '%d leaves superseded, expected 176 (F3); a pattern in regression_tol_settings().anchor_superseded matches more or fewer leaves than it was written for.', ...
    n_superseded);
assert(all(contains(leaves.key(strcmp(leaves.status, 'superseded')), '.sfa1_std1')), ...
    'test_regression_tolerance:SupersededOutside', 'a superseded leaf lies outside sfa1_std1');
fprintf('PASS  tolerance regression tier matches the anchor within tolerance (%d leaves superseded, not compared: %s)\n', ...
    n_superseded, info.superseded(1).reason);
