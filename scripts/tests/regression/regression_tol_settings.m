function T = regression_tol_settings()
% REGRESSION_TOL_SETTINGS The tolerances of the tolerance regression tier, in one place.
%
%   T = regression_tol_settings()
%
% Every tolerance the tier applies is read from here, by
% regression_tol_compare (layers 1 and 2) and regression_tol_band_check
% (layer 3), and a reference stores the copy that was in force when it was
% captured. THE VALUES ARE PROPOSED, NOT FIXED: the mutation checks of
% Plan A step 0c (FD-020) decide them, by showing that every deliberate
% error fails the tier at these values while round-off does not.
%
% THE RULE FOR ONE FLOATING-POINT LEAF (an array; REF from the reference,
% NEW from this run). With scale = the largest |REF| in the leaf, element i
% passes when
%
%     |NEW_i - REF_i|  <=  rtol * max(|REF_i|, floor_frac * scale)
%
% and the leaf passes when every element does. The left side divided by
% the max(...) is the leaf's "relative difference"; its largest value over
% the leaf is what the comparison reports, pass or fail.
%
% WHY A FLOOR, AND WHY IT IS THE LEAF'S OWN SCALE. A purely element-wise
% relative test fails on round-off alone: an entry near zero is usually the
% difference of terms of the leaf's ordinary size (a right-hand side close
% to equilibrium, a Jacobian product that nearly cancels), so its round-off
% is eps times the size of those TERMS, not eps times its own size. The
% honest unit for such an entry is the leaf's scale. floor_frac = 1 makes
% the unit the leaf's scale for every element, so the rule reads "no
% element moved by more than rtol times the largest magnitude in the
% leaf"; a smaller floor_frac tightens the test on small entries, at the
% price of sitting closer to that cancellation round-off. A leaf that is
% identically zero in the reference has scale 0 and must stay exactly zero.
%
% WHAT IS COMPARED EXACTLY, whatever the tolerances: integer, logical and
% text leaves; every leaf's size and class; where the non-finite values
% are; and the DISCRETE leaves, floating-point arrays that hold counts or
% indices (the count of positive exponents, the K that lya_K_auto settled
% on, clip counts, sort orders). Those are named by the regular
% expressions in T.discrete, matched against the leaf's key. They can
% flip under round-off only if a case sits at a threshold, and the stable
% cases were chosen away from them (regression_tol_run).
%
% THE ROUND-OFF FLOOR these sit above was measured on 2026-11-11 by
% capturing at 14 computational threads and re-running at 1, 2, 4 and 8
% (the client is multithreaded, and the tier runs in the client). At this
% tier's sizes the thread count moves NOTHING in the pointwise layer and
% nothing in the stable layer's trajectories, top-K and Benettin results,
% transient gain or sweep: 1728 of 1749 float leaves are bit-identical.
% The exceptions are the three groups under T.override below, each with
% its measured floor. So the thread count is not what the two layer
% tolerances have to absorb; a reordered expression in the model is. That
% was measured the same day with a probe that was never committed: five
% harmless reorderings in SRNNCellTypePairs (the sum in dx/dt, the release
% term of db/dt, the division by tau_d in two Jacobian blocks and in
% jacobian_times), compared with a capture made just before. 452 of 1749
% float leaves moved, none of them an exact leaf, and the worst were:
%   pointwise   6.0e-16 (dS), 3.4e-16 (JY), 2.3e-16 (J)     rtol 1e-12
%   stable      2.4e-16 trajectories; 7.7e-15 top-K basis;
%               2.0e-12 Benettin local rates; 1.3e-12 lambda_1_drift
%                                                            rtol 1e-9
%   overrides   QR 6.1e-8; memory capacity 8.9e-10;
%               finite-difference Jacobian 1.7e-10
% so each tolerance sits 500 times or more above what a harmless
% reordering does (160 times for QR). The same probe is the reason for
% floor_frac = 1: purely element-wise, dS moved by 4e-12 and S by 9e-12,
% which a 1e-12 test on each element would have failed.
%
% Fields:
%   rtol.pointwise  layer 1, keys 'pointwise.*'
%   rtol.stable     layer 2, keys 'stable.*'
%   rtol.band       layer 3: value rows of the rapid / stages / full digests
%                   beyond this are LISTED by regression_tol_band_check;
%                   they do not fail it
%   floor_frac      see the rule above
%   discrete        cell of regular expressions: float leaves compared exactly
%   override        n x 2 cell {pattern, rtol}: a leaf whose key matches takes
%                   that rtol instead of its layer's (first match wins).
%                   Every row needs its reason written beside it.
%
% See also: regression_tol_compare, regression_tol_check,
%           regression_tol_band_check, make_regression_tol_reference
T = struct();
T.rtol = struct('pointwise', 1e-12, 'stable', 1e-9, 'band', 0.1);
T.floor_frac = 1;
T.discrete = {['(^|\.)(n_positive|n_positive_at_K|K_used|K|K_requested|retries|D_KY_resolved|' ...
    'sort_idx|tau_a_clipped|hold_len|config_idx|network_seed|i_peak|N|n_per_type|n_b|n_g|reps)$']};
T.override = { ...
    ... % The largest departure of the top-K basis from orthonormality: a
    ... % measurement OF round-off (about 1e-15), so it moves by tens of
    ... % percent of itself between thread counts. What must hold is that it
    ... % stays round-off: 1e3 keeps it below about 1e-12.
    '(^|\.)orth_defect_max$', 1e3; ...
    ... % A central difference with h = 1e-6 divides the right-hand side's
    ... % round-off by 2h: a reordered sum in dynamics_fast moved this leaf
    ... % by 1.7e-10 relative while moving dS itself by 6e-16.
    '(^|\.)finite_difference_J$', 1e-7; ...
    ... % The full-spectrum QR estimator integrates N^2 variational equations
    ... % with ode45, whose adaptive step control turns a last-bit difference
    ... % into one of up to 8e-8 relative (measured between 1, 2, 4, 8 and 14
    ... % threads on the reduced sfa3_std2 network).
    '^stable\.estimators\.\w+\.qr\.', 1e-5; ...
    ... % Memory capacity solves a ridge regression on reservoir states that
    ... % are nearly collinear in a stable network; the solve amplifies
    ... % round-off to as much as 6e-9 relative (same measurement).
    '^stable\.mc\.\w+\.(MC|R2)$', 1e-6};
end
