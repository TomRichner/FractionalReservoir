% test_regression_rapid.m
%
% The ultrarapid regression (about 20 s): runs the model in the regression
% preset's three conditions with noise off and on, a small sweep through
% ParamSpaceAnalysis2 and a short memory-capacity run through
% SRNN_ESN_reservoir, and compares everything they produce, bit for bit, with
% scripts/tests/fixtures/regression_rapid_digest.mat.
%
% Run it after every edit to the model class, the reservoir class or the
% sweep driver. It is quick enough to be in run_all_tests('quick').
%
% If it fails because a change was MEANT to move results, the reference is
% replaced deliberately with make_regression_reference('rapid'), and the
% commit message says why. The full tier (every stage and figure, about 50
% minutes) is test_regression_bundle. See docs/ISSUES.md, FD-003.

[ok, diffs] = regression_rapid_check();
assert(ok, 'test_regression_rapid:Differences', ...
    '%d difference(s) from the ultrarapid regression reference; see the list above.', ...
    nnz(~diffs.ignored));
fprintf('PASS  ultrarapid regression matches the reference\n');
