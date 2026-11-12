% test_regression_rapid_chain.m
%
% The ultrarapid regression for chained vesicle-pool depression (about half
% a minute): runs chain networks (the regression network with its four
% routes switched to std.model = 'chain', noise off and on; a
% four-cell-type network with a chain beside every other kind of route; a
% clearly chaotic chain network), the three Lyapunov estimators, transient
% gain, a small sweep through ParamSpaceAnalysis2 and a short
% memory-capacity run through SRNN_ESN_reservoir, and compares everything
% they produce, bit for bit, with
% scripts/tests/fixtures/regression_rapid_chain_digest.mat.
%
% test_regression_rapid holds no chain route; this is its twin. Run both
% after every edit to the model class, the route models, the reservoir
% class or the sweep driver that should move no number.
%
% If it fails because a change was MEANT to move results, the reference is
% replaced deliberately with make_regression_reference('rapid_chain'), and
% the commit message says why. A rewrite of the chain's arithmetic (which
% changes the last bit) is gated by test_regression_tolerance_chain. See
% docs/Plan_A_chained_vessicles.md, section 7.3, step 5.

[ok, diffs] = regression_rapid_chain_check();
assert(ok, 'test_regression_rapid_chain:Differences', ...
    '%d difference(s) from the chain''s ultrarapid regression reference; see the list above.', ...
    nnz(~diffs.ignored));
fprintf('PASS  ultrarapid regression for chained depression matches the reference\n');
