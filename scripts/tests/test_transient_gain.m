% test_transient_gain.m - SRNNCellTypePairs.transient_gain, leading_direction_at
% and excursion_samples against exact propagators and synthetic series.
%
% Checks:
%   1. SYNTHETIC frozen_x (a hand-made params struct, opts.Jxx given):
%      a symmetric (normal) J with every eigenvalue < 0 gives
%      G.worst(t) = exp(alpha t) -- nothing non-normal to amplify -- and
%      G.noise <= G.worst; the 2 x 2 non-normal [-1 100; 0 -2] peaks well
%      above 1 and matches norm(expm(J t)) at every report time. Both at
%      h = 0.25 ms so Heun's O(h^2) error sits below the 1e-6 tolerance.
%   2. THE ALL-BLOCKS NETWORK (n = 12, three types, SFA + dual STD + dual STF,
%      the make_pair_model of test_jacobian_times) at fs = 4000:
%      frozen_full through jacobian_times equals ||P_x expm(J t) P_x'|| to
%      1e-4 at every report time, frozen_x equals ||expm(J_xx t)||;
%      G(0) = 1 exactly for every variant; worst >= noise and worst >= every
%      direction reading; v_opt is a unit vector; frac_E in [0, 1];
%      a direction of the wrong length errors; an active horizon past the
%      trajectory errors.
%   3. PAPER PHYSICS at n = 40 (test_lyapunov_topk's networks A, chaotic
%      no_adaptation, and B, stable sfa3_std2; seed [1 2]): all three
%      variants over 1 s with the E/I difference, E/I sum and leading
%      Lyapunov directions -- finite, G(0) = 1, printed. The ordering of the
%      variants is PRINTED, NOT ASSERTED: on this 40-neuron stable net the
%      ACTIVE G_max (3.5) exceeds both frozen ones (1.5), i.e. the drift of
%      the state through regions of larger local expansion adds more than
%      adaptation's feedback removes over 1 s (lambda_1 = -0.13 but the
%      local rate fluctuates). Whether that holds at n = 500 is what the
%      stage measures; do not encode an expectation here.
%   4. leading_direction_at on A: unit norm, and its x part has |cos| > 0.9
%      with the x part of Q_final(:, 1) from a top-K run ending at the same
%      index (the trajectory's end).
%   5. excursion_samples on a synthetic +-1 series with known runs: the
%      onsets found are exactly the runs long enough AND preceded by enough
%      quiet; the quiet midpoints are exactly the long negative runs.
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths has run.
%
% See also: SRNNCellTypePairs.transient_gain, run_transient_gain,
%           test_jacobian_times, test_lyapunov_topk

fprintf('=== Testing transient_gain / leading_direction_at / excursion_samples ===\n\n');
all_passed = true;

%% 1. Synthetic frozen_x against expm
n = 6;
params = struct('n', n, 'type_indices', {{1:3, 4:6}}, ...
    'state_layout', struct('x', 1:n, 'N', n));
dt = 2.5e-4;
t_out = (0:dt:1)';
S_out = zeros(numel(t_out), n);
rng(3, 'twister');
A = randn(n); Js = -(A * A') / 5 - 0.5 * eye(n);              % symmetric, negative definite
alpha = max(eig(Js));
[G, info] = SRNNCellTypePairs.transient_gain(S_out, t_out, 1, params, 0.5, ...
    struct('variant', 'frozen_x', 'Jxx', Js, 'report_dt', 0.05));
err = max(abs(G.worst - exp(alpha * G.t)) ./ exp(alpha * G.t));
all_passed = check(sprintf('normal J: G.worst = exp(alpha t) (rel err %.1e)', err), err < 1e-6) && all_passed;
all_passed = check('normal J: G.noise <= G.worst, G(0) = 1, peak at t = 0', ...
    all(G.noise <= G.worst + 1e-12) && G.worst(1) == 1 && info.t_peak == 0) && all_passed;

J2 = [-1 100; 0 -2];
params2 = struct('n', 2, 'type_indices', {{1, 2}}, 'state_layout', struct('x', 1:2, 'N', 2));
[G2, info2] = SRNNCellTypePairs.transient_gain(zeros(numel(t_out), 2), t_out, 1, params2, 1, ...
    struct('variant', 'frozen_x', 'Jxx', J2, 'report_dt', 0.02));
ref = arrayfun(@(tt) norm(expm(J2 * tt)), G2.t);
err2 = max(abs(G2.worst - ref) ./ ref);
fprintf('  2x2 non-normal: G_max %.2f at t = %.2f s (expm peak %.2f)\n', info2.G_max, info2.t_peak, max(ref));
all_passed = check(sprintf('non-normal 2x2: G.worst = ||expm(J t)|| at every report time (rel err %.1e)', err2), ...
    err2 < 1e-6) && all_passed;
all_passed = check('non-normal 2x2: peak gain well above 1', info2.G_max > 10) && all_passed;

%% 2. The all-blocks network against expm of the full Jacobian
m = make_pair_model();
m.build();
evalc('m.run();');
p = m.get_params(); S = m.S_out; t = m.t_out; nn = m.n;
i0 = find(t >= 0.5, 1);
J = full(SRNNCellTypePairs.compute_Jacobian_fast(S(i0, :)', p));
ix = p.state_layout.x;
dirs = struct('ei_diff', [ones(p.n_per_type(1), 1); -ones(nn - p.n_per_type(1), 1)], ...
    'ei_sum', ones(nn, 1), 'rand', randn(nn, 1));
ok0 = true; ok_ord = true; ok_unit = true; err_full = 0; err_x = 0;
for var = {'frozen_x', 'frozen_full', 'active'}
    [Gv, iv] = SRNNCellTypePairs.transient_gain(S, t, i0, p, 0.25, ...
        struct('variant', var{1}, 'report_dt', 0.025, 'directions', dirs));
    ok0 = ok0 && Gv.worst(1) == 1 && Gv.noise(1) == 1;
    ok_ord = ok_ord && all(Gv.worst >= Gv.noise - 1e-12) && ...
        all(Gv.worst >= Gv.dir.ei_diff - 1e-12) && all(Gv.worst >= Gv.dir.ei_sum - 1e-12) && ...
        all(Gv.worst >= Gv.dir.rand - 1e-12);
    ok_unit = ok_unit && abs(norm(iv.v_opt) - 1) < 1e-12 && iv.frac_E >= 0 && iv.frac_E <= 1 && ...
        iv.participation >= 1 && iv.participation <= nn + 1e-9;
    switch var{1}
        case 'frozen_full'
            ref = arrayfun(@(tt) norm(subsref(expm(J * tt), substruct('()', {ix, ix}))), Gv.t);
            err_full = max(abs(Gv.worst - ref) ./ ref);
        case 'frozen_x'
            ref = arrayfun(@(tt) norm(expm(J(ix, ix) * tt)), Gv.t);
            err_x = max(abs(Gv.worst - ref) ./ ref);
    end
    fprintf('  all-blocks %-11s G_max %.4f at %.3f s, frac_E %.2f, participation %.1f\n', ...
        var{1}, iv.G_max, iv.t_peak, iv.frac_E, iv.participation);
end
all_passed = check(sprintf('frozen_full = ||P_x expm(J t) P_x''|| (rel err %.1e)', err_full), err_full < 1e-4) && all_passed;
all_passed = check(sprintf('frozen_x = ||expm(J_xx t)|| (rel err %.1e)', err_x), err_x < 1e-4) && all_passed;
all_passed = check('G(0) = 1 for every variant', ok0) && all_passed;
all_passed = check('worst >= noise and >= every direction reading', ok_ord) && all_passed;
all_passed = check('v_opt unit, frac_E in [0, 1], participation in [1, n]', ok_unit) && all_passed;
all_passed = check('a direction of the wrong length errors', throws_id(@() SRNNCellTypePairs.transient_gain( ...
    S, t, i0, p, 0.1, struct('directions', struct('bad', ones(3, 1)))), 'SRNNCellTypePairs:BadDirection')) && all_passed;
all_passed = check('an active horizon past the trajectory errors', throws_id(@() SRNNCellTypePairs.transient_gain( ...
    S, t, i0, p, 10, struct('variant', 'active')), 'SRNNCellTypePairs:HorizonExceedsTrajectory')) && all_passed;
all_passed = check('an unknown variant errors', throws_id(@() SRNNCellTypePairs.transient_gain( ...
    S, t, i0, p, 0.1, struct('variant', 'nope')), 'SRNNCellTypePairs:BadVariant')) && all_passed;

%% 3. Paper physics at n = 40
P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25';
base = {'n', 40, 'indegree', 8, 'F_tracks_network', true, 'sigma_u_noise', 0, ...
    'ode_solver', 'rk4', 'rng_seeds', [1 2], 'T_range', [0 12], ...
    'lya_T_interval', [4 12], 'lya_warmup', 4, 'lya_dt', 0.1, 'store_full_state', true, ...
    'lya_method', 'topk', 'lya_K', 5};
mA = run_preset(P, 'no_adaptation', base);
mB = run_preset(P, 'sfa3_std2', base);
fprintf('  A: lambda_1 %+.3f (chaotic);  B: lambda_1 %+.3f (stable)\n', mA.lya_results.LLE, mB.lya_results.LLE);
Gmax = struct();
for net = {'A', 'B'}
    if strcmp(net{1}, 'A'); mm = mA; else; mm = mB; end
    p = mm.get_params(); S = mm.S_out; t = mm.t_out; nn = mm.n;
    i0 = find(t >= 8, 1);
    nE = p.n_per_type(1);
    v_lyap = SRNNCellTypePairs.leading_direction_at(S, t, i0, p, 4);
    dirs = struct('ei_diff', [ones(nE, 1) / sqrt(nE); -ones(nn - nE, 1) / sqrt(nn - nE)], ...
        'ei_sum', ones(nn, 1), 'lyap', v_lyap);
    fin = true; g0 = true;
    for var = {'frozen_x', 'frozen_full', 'active'}
        [Gv, iv] = SRNNCellTypePairs.transient_gain(S, t, i0, p, 1, ...
            struct('variant', var{1}, 'report_dt', 0.02, 'directions', dirs));
        fin = fin && all(isfinite(Gv.worst)) && all(isfinite(Gv.noise)) && all(isfinite(Gv.dir.lyap));
        g0 = g0 && Gv.worst(1) == 1;
        Gmax.(net{1}).(var{1}) = iv.G_max;
        fprintf('  net %s %-11s G_max %8.3f at %.2f s | noise %.3f, ei_diff %.3f, ei_sum %.3f, lyap %.3f at peak | cos(v_opt, lyap) %.2f, frac_E %.2f, %.1f s\n', ...
            net{1}, var{1}, iv.G_max, iv.t_peak, Gv.noise(iv.i_peak), Gv.dir.ei_diff(iv.i_peak), ...
            Gv.dir.ei_sum(iv.i_peak), Gv.dir.lyap(iv.i_peak), iv.cos_opt.lyap, iv.frac_E, iv.seconds);
    end
    all_passed = check(sprintf('net %s: every variant finite with G(0) = 1', net{1}), fin && g0) && all_passed;
end
fprintf('  net B ordering (not asserted): frozen_x %.3f, frozen_full %.3f, active %.3f\n', ...
    Gmax.B.frozen_x, Gmax.B.frozen_full, Gmax.B.active);
all_passed = check('net A (no adaptation): frozen_full equals frozen_x (N = n)', abs(Gmax.A.frozen_full - Gmax.A.frozen_x) < 1e-9) && all_passed;

%% 4. leading_direction_at vs Q_final(:, 1) on A (both at the trajectory's end)
p = mA.get_params(); S = mA.S_out; t = mA.t_out;
[v_x, v_full] = SRNNCellTypePairs.leading_direction_at(S, t, numel(t), p, 6);
q1 = mA.lya_results.Q_final(:, 1);
qx = q1(p.state_layout.x); qx = qx / norm(qx);
c = abs(v_x' * qx);
all_passed = check('leading_direction_at returns unit vectors', abs(norm(v_x) - 1) < 1e-12 && abs(norm(v_full) - 1) < 1e-12) && all_passed;
all_passed = check(sprintf('x part aligned with Q_final(:,1) on the chaotic net (|cos| %.3f)', c), c > 0.9) && all_passed;
all_passed = check('a warm-up before the trajectory start warns and clamps', ...
    warns_id(@() SRNNCellTypePairs.leading_direction_at(S, t, 10, p, 1), 'SRNNCellTypePairs:LeadingDirectionWarmupClamped')) && all_passed;

%% 5. excursion_samples on a synthetic series
%          idx: 1-5 neg | 6-9 pos(4) | 10-11 neg(2) | 12-16 pos(5) | 17-36 neg(20) | 37-38 pos(2) | 39-45 neg(7) | 46-50 pos(5)
r = [-ones(5, 1); ones(4, 1); -ones(2, 1); ones(5, 1); -ones(20, 1); ones(2, 1); -ones(7, 1); ones(5, 1)];
tl = (0:numel(r) - 1)' * 0.05;
E = SRNNCellTypePairs.excursion_samples(r, tl, 4, 20);
% onsets: run 6-9 (len 4, preceded by 5 neg) yes; 12-16 (preceded by 2) no;
% 37-38 (len 2) no; 46-50 (len 5, preceded by 7) yes.  quiet: 17-36 -> mid 26.5 -> 27
all_passed = check(sprintf('onsets at segments 6 and 46 (%s)', mat2str(E.onset_idx')), isequal(E.onset_idx, [6; 46])) && all_passed;
all_passed = check(sprintf('quiet midpoint at segment 27 (%s)', mat2str(E.quiet_idx')), isequal(E.quiet_idx, 27)) && all_passed;
all_passed = check('times follow the index', isequal(E.onset_t, tl([6; 46])) && E.n_onset == 2 && E.n_quiet == 1) && all_passed;
E0 = SRNNCellTypePairs.excursion_samples(-ones(10, 1), (1:10)', 4, 20);
all_passed = check('no runs -> empty columns, zero counts', isempty(E0.onset_idx) && isempty(E0.quiet_idx) && E0.n_onset == 0 && E0.n_quiet == 0) && all_passed;

fprintf('\n');
if all_passed
    fprintf('=== ALL transient_gain TESTS PASSED ===\n');
else
    fprintf('=== SOME transient_gain TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function model = make_pair_model()
% The all-blocks network of test_SRNNCellTypePairs / test_jacobian_times,
% at fs = 4000 for 1 s so the Heun propagator can be compared with expm.
synapse_config = struct();
synapse_config.E.E.std = struct('tau_rec', [0.3 1], 'tau_rel', 0.2);
synapse_config.E.PV.stf = struct('tau_dec', 0.5, 'tau_fac', 0.25, 'G', 2);
synapse_config.E.SST.std = struct('tau_rec', 0.8, 'tau_rel', 0.25);
synapse_config.E.SST.stf = struct( ...
    'tau_dec', [0.4 1.2], 'tau_fac', [0.2 0.3], 'G', [1.5 2]);
synapse_config.PV.PV.std = struct('tau_rec', 0.6, 'tau_rel', 0.3);
synapse_config.SST.PV.stf = struct('tau_dec', 0.9, 'tau_fac', 0.4, 'G', 1.8);
model = SRNNCellTypePairs( ...
    'n', 12, 'indegree', 4, ...
    'n_cellTypes', 3, 'cell_type_names', {'E', 'PV', 'SST'}, ...
    'f', [0.5 0.3 0.2], ...
    'mu_tilde_relative', [0.1 -0.2 -0.08], ...
    'sigma_tilde_relative', [0.01 0.02 0.015], ...
    'tau_a', {[0.25 1], [], 2}, ...
    'c', [0.05 0 0.03], 'synapse_config', synapse_config, ...
    'T_range', [0 1], 'fs', 4000, 'ode_solver', 'rk4', 'lya_method', 'none', ...
    'store_full_state', true);
end

function m = run_preset(P, cond, base) %#ok<INUSD>  used inside evalc
m = [];
evalc('m = build_from_preset(P, cond, base{:});');
evalc('m.run();');
end

function ok = throws_id(fn, id)
ok = false;
try
    fn();
catch ME
    ok = strcmp(ME.identifier, id);
end
end

function ok = warns_id(fn, id)
lastwarn('');
evalc('fn();');
[~, wid] = lastwarn();
ok = strcmp(wid, id);
end

function passed = check(name, condition)
if condition
    fprintf('  %s: PASS\n', name);
    passed = true;
else
    fprintf('  %s: FAIL\n', name);
    passed = false;
end
end
