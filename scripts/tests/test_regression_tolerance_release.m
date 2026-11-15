% test_regression_tolerance_release.m
%
% The tolerance regression tier's RELEASE suite (about 15 s, in the client;
% no pool is needed): runs with spontaneous release on (FD-027: a
% presynaptic neuron's release drive over a step is r dt + dR, dR a gamma
% increment; its mean r_s is a leak in the drift and the rest a kick after
% each step), compared ELEMENT BY ELEMENT, within the tier's tolerances,
% with scripts/tests/fixtures/regression_tol_release.mat.
%
%   - pointwise, for every route model with the leak on (no synapse
%     dynamics, one pool, the chained pool, the product of two pools, the
%     chain with release-coupled facilitation, and a four-cell-type
%     network with every other variant and a different release in each
%     presynaptic type): the right-hand side, the assembled Jacobian of
%     the mean dynamics, the matrix-free product, the parameter view and
%     the factors at stored states; the state after a kick by a fixed
%     pattern of increments (a zero draw, a small, a large and a huge
%     one), what each kick released on every route, the route models'
%     closed-form steady states, the increments a run draws and its kick
%     handle at one step, the generator's own output, and the tau_d filter
%     of the release readout on a fixed input;
%   - stable networks end to end: trajectories with the kicks, top-K, the
%     realized release of every route per step, the chain with both noise
%     sources, transient gain, the finite-difference Jacobian, Benettin, a
%     sweep with release_sigma_all as an axis, and lean memory capacity
%     with the 'release_filtered' readout.
%
% No other reference of the tier holds a run with release_rate or
% release_sigma above zero, so this is the one that says "the release
% noise is unchanged". The exact twin is test_regression_rapid_release.
%
% IF IT FAILS, the change moved a number by more than round-off: find the
% defect. Do not loosen a tolerance (regression_tol_settings). The
% reference is replaced only deliberately (make_regression_tol_reference's
% header says how), with the reason in the commit message. To bisect,
% capture a pre-edit reference with
% make_regression_tol_reference('<name>', 'suite', 'release').
% See docs/Plan_release_noise.md, section 4, and
% regression_tol_mutations_release for the deliberate errors this suite
% was shown to catch.

[ok, leaves] = regression_tol_check('release');
assert(ok, 'test_regression_tolerance_release:Differences', ...
    '%d leaf(s) of the tolerance tier''s release reference are outside tolerance or missing; see the list above.', ...
    nnz(ismember(leaves.status, {'FAIL', 'missing'})));
fprintf('PASS  tolerance regression tier, release suite, matches its reference within tolerance\n');
