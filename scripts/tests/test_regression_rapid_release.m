% test_regression_rapid_release.m
%
% The ultrarapid regression for release noise (about half a minute): runs
% the four conditions of the paper on the regression network with the
% gamma process in the synaptic output (FD-027), the chain with both noise
% sources and a seed schedule, the leak alone, a four-cell-type network
% with a different release in each type, top-K and Benettin, a sweep with
% release_sigma_all as a grid axis through ParamSpaceAnalysis2 and lean
% memory-capacity runs with the 'release_filtered' readout through
% SRNN_ESN_reservoir, and compares everything they produce, bit for bit,
% with scripts/tests/fixtures/regression_rapid_release_digest.mat.
%
% test_regression_rapid and test_regression_rapid_chain hold no run with
% spontaneous release; this is their twin. Run all three after every edit
% to the model class, the route models, the integrators, the reservoir
% class or the sweep driver that should move no number.
%
% If it fails because a change was MEANT to move results, the reference is
% replaced deliberately with make_regression_reference('rapid_release'),
% and the commit message says why. A rewrite of the release arithmetic
% (which changes the last bit) is gated by test_regression_tolerance_release.
% See docs/Plan_release_noise.md, section 4.

[ok, diffs] = regression_rapid_release_check();
assert(ok, 'test_regression_rapid_release:Differences', ...
    '%d difference(s) from the release noise''s ultrarapid regression reference; see the list above.', ...
    nnz(~diffs.ignored));
fprintf('PASS  ultrarapid regression for release noise matches the reference\n');
