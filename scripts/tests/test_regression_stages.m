% test_regression_stages.m
%
% The middle regression tier (a few minutes): calls three analysis stages
% directly on the 100-neuron regression network (Lyapunov spectrum,
% transient gain, the memory-capacity example), builds the figures that read
% them, and compares everything the stages save, and the list of figure
% files, with scripts/tests/fixtures/regression_stages_digest.mat.
%
% Run it before committing a change to an analysis stage, a figure, or
% anything the stages call in the model class. It is on run_all_tests' slow
% list.
%
% The three tiers: test_regression_rapid (seconds; the three classes
% directly), this one (minutes; stages and figures), test_regression_bundle
% (about 50 minutes; the whole bundle). If this fails because a change was
% MEANT to move results, the reference is replaced deliberately with
% make_regression_reference('stages'), and the commit message says why. See
% docs/ISSUES.md, FD-003.

[ok, diffs] = regression_stages_check();
assert(ok, 'test_regression_stages:Differences', ...
    '%d difference(s) from the stages regression reference; see the list above.', ...
    nnz(~diffs.ignored));
fprintf('PASS  stages regression matches the reference\n');
