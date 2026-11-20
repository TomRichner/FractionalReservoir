% test_lyapunov_topk_ritz.m - The projected Jacobian H = Q' J Q that
% lyapunov_topk returns at requested states (opts.H_idx), and its
% eigenvalues, the Ritz values of the Jacobian on the tracked subspace
% (FD-033).
%
% Claims checked:
%   1. LINEAR SYSTEM, K = N. dY/dt = J Y with a constant J whose eigenvalues
%      are known and complex: Q is square, H is orthogonally similar to J,
%      so the Ritz values equal eig(J) to round-off at every requested
%      state, on the grid or inside a segment. On the assembled path and on
%      the matrix-free path.
%   2. LINEAR SYSTEM, K < N. J = U T U' with T in real Schur form, its
%      diagonal blocks in descending order of real part and a gap after the
%      fifth eigenvalue; the leading block is a complex pair at 15 Hz, ABOVE
%      the 10 Hz at which the finite-interval form log(Q_k' Q_{k+1} R_k) / tau
%      would alias at lya_dt = 0.05 s. After alignment the top-5 subspace is
%      the leading invariant subspace, so the Ritz values equal the five
%      leading eigenvalues, for a normal T and for a non-normal one. The
%      leading 2 x 2 block of H then carries the pair alone.
%   3. THE OPTION MOVES NOTHING ELSE. Every field lyapunov_topk returned
%      before the option existed is identical (isequaln, every field but the
%      wall time `seconds`) with H_idx given or not, and identical to the
%      copy of the function at the commit before FD-033 when that copy is on
%      the path (lyapunov_topk_head; a one-off gate, skipped otherwise).
%   4. ON THE MODEL (stable 40-neuron network, every Jacobian block live):
%      (a) H equals Q' (compute_Jacobian_fast Q) with the returned Q, Q is
%          orthonormal, and the Q at the end of the run equals Q_final;
%      (b) DIAGONAL: the trapezoid mean over a segment of H_jj at every
%          grid point of the segment equals that segment's local Lyapunov
%          rate log(R_jj) / tau of raw column j;
%      (c) TRACE: the same for the sum over the leading k columns.
%      (b) and (c) are identities of the exact flow; here the propagator is
%      Heun and the mean a trapezoid, both second order in the step. The
%      tolerance was MEASURED before it was fixed (see the constants below).
%   5. A requested state outside the propagated range gives H_valid false
%      and NaN, not an error; a bad H_idx errors :InvalidHIdx.
%   6. THE CLASS'S DOOR, SRNNCellTypePairs.projected_jacobian: it is the
%      same pass (H identical to a direct lyapunov_topk call on the class's
%      grid); K is capped at N; the block fractions of the subspace sum to
%      1; at K = N the Ritz values equal the eigenvalues of the assembled
%      Jacobian, WITH AND WITHOUT a state scaling (J_s = D^-1 J D is similar
%      to J), while at K < N the scaling moves them; a scaling of ones is
%      no scaling, to the bit; a bad scale errors :BadStateScale.
%   7. SRNNCellTypePairs.state_type_scales: one scale per kind of variable
%      (x, a1..a3, b1, b2 on this network), each the pooled standard
%      deviation of that kind's variables (root mean of their variances in
%      time over the rows), written into every variable of the kind; a kind that does not fluctuate gets 1 and is
%      named.
%   8. RESIDUALS AND A TIME-VARYING, NON-NORMAL SYSTEM WITH A KNOWN ANSWER
%      (ritz_diagnostics; H_gram_t). The Ritz values are eigenvalues of J
%      only when the tracked subspace is invariant under J, and the
%      residual ||J Q - Q H||_F / ||J Q||_F says whether it is.
%      (a) It is 0 (to the 1e-8 floor of a difference of squares) at K = N
%          and on the invariant subspaces of section 2.
%      (b) A ROTATING system, J(t) = R(t) A R(t)' with R = expm(Omega t),
%          Omega skew, A stable and non-normal. In the rotating frame the
%          tangent dynamics are dy/dt = (A - Omega) y, constant, so the
%          tracked subspace is R(t) V with V the leading invariant subspace
%          of B = A - Omega, NOT of A. The answer is then known in closed
%          form and is neither eig(A) nor eig(B): the Ritz values are
%          eig(V' A V), with residual ||A V - V V' A V||_F / ||A V||_F > 0,
%          and the Lyapunov exponents are the real parts of eig(B) (their
%          sum is checked: single exponents need longer to settle). This is
%          the case the model is in: a subspace carried by the dynamics,
%          not invariant under the instantaneous Jacobian.
%      (c) On the model: G = (J Q)' (J Q) against the assembled Jacobian;
%          the residual and the pair residuals from H and G alone equal the
%          ones computed directly from J, Q and the Ritz vectors; the
%          projected numerical abscissa bounds every Ritz real part and is
%          bounded by the full one.
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths has run.
%
% See also: lyapunov_topk, test_lyapunov_topk, run_eig_heatmap, fig_eig_heatmap

fprintf('=== Testing lyapunov_topk''s projected Jacobian (opts.H_idx) ===\n\n');
all_passed = true;

% Tolerances of the model identities 4(b) and 4(c), in 1/s. Measured
% 2026-11-19 on this test's network (fs 400, lya_dt 0.1, K 20, 20 segments):
% see the printed lines; the limits are the measured maxima with a margin of
% about ten.
TOL_DIAG  = 1e-4;
TOL_TRACE = 1e-4;

%% 1. Linear system, K = N
rng_state = rng; rng(7, 'twister');
N = 8;
[U, ~] = qr(randn(N));
T_full = schur_form([-0.2 + 2i * pi * 12, -0.6, -1.5 + 2i * pi * 2, -3, -5 + 1i], 0.8);
J = U * T_full * U';
rng(rng_state);
fs = 1000; t = (0:1 / fs:3)'; X = zeros(numel(t), N);
grid_fn = @simple_grid;
idx = [1, 137, 1001, 2444, numel(t)];          % on and off the 0.05 s grid, first and last row
o = struct('seed', 3, 'grid_fn', grid_fn, 'H_idx', idx);
r_asm = lyapunov_topk(X, t, fs, 0.05, [0 3], 0, N, @(S, p) J, [], o);
o_mf = o; o_mf.jac_times = @(S, Y, p) J * Y;
r_mf = lyapunov_topk(X, t, fs, 0.05, [0 3], 0, N, @(S, p) error('unused'), [], o_mf);
ev = eig(J);
err_asm = max(arrayfun(@(s) match_err(eig(r_asm.H_t(:, :, s)), ev), 1:numel(idx)));
err_mf  = max(arrayfun(@(s) match_err(eig(r_mf.H_t(:, :, s)), ev), 1:numel(idx)));
fprintf('  K = N linear: max |Ritz - eig(J)| assembled %.2e, matrix-free %.2e\n', err_asm, err_mf);
all_passed = check('K = N: Ritz values equal eig(J) at every requested state (both paths, 1e-10)', ...
    all(r_asm.H_valid) && all(r_mf.H_valid) && err_asm < 1e-10 && err_mf < 1e-10) && all_passed;
all_passed = check('H_idx and H_time are returned as requested', ...
    isequal(r_asm.H_idx, idx(:)) && isequal(r_asm.H_time, t(idx))) && all_passed;

%% 2. Linear system, K < N: the leading invariant subspace
lead = [-0.1 + 2i * pi * 15, -3, -4 + 2i * pi * 3];       % five eigenvalues, then a gap;
rest = [-12, -13, -14 + 2i, -15, -16, -17];               % gaps of 2.9 and 8 1/s, 7 s to align
t2 = (0:1 / fs:8)'; X2 = zeros(numel(t2), 12);
idx2 = [7003, 7250, 7777, numel(t2)];
names = {'normal', 'non-normal'}; couplings = [0, 0.8];
for v = 1:2
    rng_state = rng; rng(11, 'twister');
    T2 = schur_form([lead, rest], couplings(v));
    [U2, ~] = qr(randn(12));
    rng(rng_state);
    J2 = U2 * T2 * U2';
    o2 = struct('seed', 5, 'grid_fn', grid_fn, 'H_idx', idx2, 'jac_times', @(S, Y, p) J2 * Y);
    r2 = lyapunov_topk(X2, t2, fs, 0.05, [7 8], 7, 5, @(S, p) error('unused'), [], o2);
    want = eig(T2(1:5, 1:5));
    err5 = max(arrayfun(@(s) match_err(eig(r2.H_t(:, :, s)), want), 1:numel(idx2)));
    pair = [lead(1); conj(lead(1))];
    err2 = max(arrayfun(@(s) match_err(eig(r2.H_t(1:2, 1:2, s)), pair), 1:numel(idx2)));
    fprintf('  K = 5 of 12, %s: max |Ritz - leading eig| %.2e; leading 2 x 2 block vs the 15 Hz pair %.2e\n', ...
        names{v}, err5, err2);
    all_passed = check(sprintf('K < N, %s J: Ritz values equal the five leading eigenvalues (1e-6)', names{v}), ...
        err5 < 1e-6) && all_passed;
    all_passed = check(sprintf('K < N, %s J: the leading 2 x 2 block carries the 15 Hz pair (1e-6)', names{v}), ...
        err2 < 1e-6) && all_passed;
end

%% 3. The option moves nothing else
o_off = rmfield(o2, 'H_idx');
r_off = lyapunov_topk(X2, t2, fs, 0.05, [7 8], 7, 5, @(S, p) error('unused'), [], o_off);
all_passed = check('option off adds no field', ~any(isfield(r_off, {'H_t', 'H_idx', 'H_time', 'H_valid', 'Q_t'}))) && all_passed;
all_passed = check('every earlier field identical with the option on (linear, matrix-free)', ...
    same_but_seconds(r_off, r2)) && all_passed;
o_off_asm = rmfield(o, 'H_idx');
r_off_asm = lyapunov_topk(X, t, fs, 0.05, [0 3], 0, N, @(S, p) J, [], o_off_asm);
all_passed = check('every earlier field identical with the option on (linear, assembled)', ...
    same_but_seconds(r_off_asm, r_asm)) && all_passed;

%% 4. On the model
P = 'celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25';
m = [];
evalc(['m = build_from_preset(P, ''sfa3_std2'', ''n'', 40, ''indegree'', 8, ''F_tracks_network'', true, ' ...
    '''sigma_u_noise'', 0, ''ode_solver'', ''rk4'', ''rng_seeds'', [1 2], ''T_range'', [0 12], ' ...
    '''lya_method'', ''none'', ''store_full_state'', true);']);
evalc('m.run();');
params = m.get_params();
K = 20; lya_dt = 0.1; dec = round(lya_dt * m.fs);
jt = @(S, Y, p) SRNNCellTypePairs.jacobian_times(S, Y, p);
jf = @(S, p) SRNNCellTypePairs.compute_Jacobian_fast(S, p);
om = struct('seed', 9, 'grid_fn', grid_fn, 'jac_times', jt);
r0 = lyapunov_topk(m.S_out, m.t_out, m.fs, lya_dt, [4 12], 4, K, jf, params, om);
seg_start = round((r0.t_lya - m.t_out(1)) * m.fs) + 1;       % row of X at each segment start
k0 = find(r0.t_lya >= 6, 1); n_seg = 20;
rows = (seg_start(k0):seg_start(k0) + n_seg * dec)';
om_H = om; om_H.H_idx = [rows; numel(m.t_out)]; om_H.H_return_Q = true;
rH = lyapunov_topk(m.S_out, m.t_out, m.fs, lya_dt, [4 12], 4, K, jf, params, om_H);
all_passed = check('model: every earlier field identical with the option on', same_but_seconds(r0, rH)) && all_passed;
if exist('lyapunov_topk_head', 'file') == 2
    r_head = lyapunov_topk_head(m.S_out, m.t_out, m.fs, lya_dt, [4 12], 4, K, jf, params, om);
    om_asm = rmfield(om, 'jac_times');
    r_head_asm = lyapunov_topk_head(m.S_out, m.t_out, m.fs, lya_dt, [4 12], 4, 5, jf, params, om_asm);
    r_now_asm  = lyapunov_topk(m.S_out, m.t_out, m.fs, lya_dt, [4 12], 4, 5, jf, params, om_asm);
    all_passed = check('model: identical to the function before FD-033 (matrix-free and assembled)', ...
        same_but_seconds(r_head, r0) && same_but_seconds(r_head, rH) && ...
        same_but_seconds(r_head_asm, r_now_asm)) && all_passed;
else
    fprintf('  (lyapunov_topk_head not on the path: the comparison with the earlier function is skipped)\n');
end

% (a) H against the assembled Jacobian, Q orthonormal, last Q = Q_final
err_a = 0; err_orth = 0;
for s = [1, 7, round(numel(rows) / 2), numel(rows), numel(rows) + 1]
    Qs = rH.Q_t(:, :, s);
    Jd = full(jf(m.S_out(rH.H_idx(s), :)', params));
    Href = Qs' * (Jd * Qs);
    err_a = max(err_a, max(abs(rH.H_t(:, :, s) - Href), [], 'all') / max(abs(Href), [], 'all'));
    err_orth = max(err_orth, norm(Qs' * Qs - eye(K)));
end
Q_end = rH.Q_t(:, :, end);
fprintf('  model (a): max rel |H - Q''JQ| %.2e, max ||Q''Q - I|| %.2e\n', err_a, err_orth);
all_passed = check('model (a): H = Q'' J Q against the assembled Jacobian (1e-12), Q orthonormal (1e-12)', ...
    err_a < 1e-12 && err_orth < 1e-12) && all_passed;
all_passed = check('model (a): the basis at the last row is Q_final (unsorted)', ...
    isequal(Q_end(:, rH.sort_idx), rH.Q_final)) && all_passed;

% (b), (c): segment means of the diagonal and of the leading traces
raw = zeros(size(rH.local_LE_spectrum_t));
raw(:, rH.sort_idx) = rH.local_LE_spectrum_t;             % column j = basis vector j
d = zeros(numel(rows), K);
for s = 1:numel(rows); d(s, :) = diag(rH.H_t(:, :, s))'; end
w = [0.5; ones(dec - 1, 1); 0.5] / dec;                   % trapezoid over one segment
err_diag = 0; err_trace = 0; scale = 0;
for q = 1:n_seg
    seg_rows = (q - 1) * dec + (1:dec + 1);
    mean_d = w' * d(seg_rows, :);
    local  = raw(k0 + q - 1, :);
    err_diag  = max(err_diag, max(abs(mean_d - local)));
    err_trace = max(err_trace, max(abs(cumsum(mean_d) - cumsum(local))));
    scale = max(scale, max(abs(local)));
end
inst = max(abs(d(1:dec:end - 1, :) - raw(k0:k0 + n_seg - 1, :)), [], 'all');
fprintf(['  model (b): max |segment mean of H_jj - local rate| %.2e 1/s (rates up to %.3g 1/s);\n' ...
    '             the instantaneous H_jj at a segment start differs from that segment''s rate by up to %.2e\n'], ...
    err_diag, scale, inst);
fprintf('  model (c): max |segment mean of the leading-k trace - sum of the top-k local rates| %.2e 1/s\n', err_trace);
all_passed = check(sprintf('model (b): diagonal of H averages to the local rates (%.0e 1/s)', TOL_DIAG), ...
    err_diag < TOL_DIAG) && all_passed;
all_passed = check(sprintf('model (c): leading traces of H average to the sums of local rates (%.0e 1/s)', TOL_TRACE), ...
    err_trace < TOL_TRACE) && all_passed;

%% 5. Out of range, and argument errors
om_out = om; om_out.H_idx = [1, round(3 * m.fs) + 4];     % iteration starts at 2 s: row 1 precedes it
r_out = lyapunov_topk(m.S_out, m.t_out, m.fs, lya_dt, [4 12], 2, 5, jf, params, om_out);
all_passed = check('a state before the propagation starts: H_valid false and NaN, the other valid', ...
    isequal(r_out.H_valid, [false; true]) && all(isnan(r_out.H_t(:, :, 1)), 'all') && ...
    all(isfinite(r_out.H_t(:, :, 2)), 'all')) && all_passed;
bad = {0, numel(m.t_out) + 1, 2.5, [3 3]};
ok = true;
for b = 1:numel(bad)
    ob = om; ob.H_idx = bad{b};
    ok = ok && throws_id(@() lyapunov_topk(m.S_out, m.t_out, m.fs, lya_dt, [4 12], 2, 5, jf, params, ob), ...
        'lyapunov_topk:InvalidHIdx');
end
all_passed = check('H_idx out of range, non-integer or repeated errors :InvalidHIdx', ok) && all_passed;

%% 6. The class's door
rows6 = [round(6 * m.fs) + 1, round(7.013 * m.fs) + 1, numel(m.t_out)];
P1 = SRNNCellTypePairs.projected_jacobian(m.S_out, m.t_out, m.fs, rows6, K, params, struct('seed', 9));
r_direct = lyapunov_topk(m.S_out, m.t_out, m.fs, 0.05, [m.t_out(rows6(1)), m.t_out(end)], ...
    m.t_out(rows6(1)) - m.t_out(1), K, jf, params, ...
    struct('seed', 9, 'grid_fn', grid_fn, 'jac_times', jt, 'H_idx', rows6));
all_passed = check('projected_jacobian is the direct pass (H identical; all states valid)', ...
    all(P1.valid) && isequal(P1.H, r_direct.H_t) && isequal(P1.LE_spectrum, r_direct.LE_spectrum) && ...
    P1.warmup_s == m.t_out(rows6(1)) - m.t_out(1)) && all_passed;
all_passed = check('block fractions sum to 1 and x fractions lie in [0, 1]', ...
    max(abs(sum(P1.block_fractions, 2) - 1)) < 1e-12 && all(P1.x_fraction >= 0 & P1.x_fraction <= 1, 'all')) && all_passed;

Nm = params.state_layout.N;
Sc = SRNNCellTypePairs.state_type_scales(m.S_out, rows6(1):numel(m.t_out), params);
P_full   = SRNNCellTypePairs.projected_jacobian(m.S_out, m.t_out, m.fs, rows6, 10 * Nm, params, struct('seed', 9));
P_full_s = SRNNCellTypePairs.projected_jacobian(m.S_out, m.t_out, m.fs, rows6, Nm, params, ...
    struct('seed', 9, 'state_scale', Sc.d));
err_full = 0; err_full_s = 0;
for s = 1:numel(rows6)
    evJ = eig(full(jf(m.S_out(rows6(s), :)', params)));
    err_full   = max(err_full,   match_err(eig(P_full.H(:, :, s)), evJ) / max(abs(evJ)));
    err_full_s = max(err_full_s, match_err(eig(P_full_s.H(:, :, s)), evJ) / max(abs(evJ)));
end
fprintf('  K = N on the model: max rel |Ritz - eig(J)| raw %.2e, scaled %.2e (scales %s = %s)\n', ...
    err_full, err_full_s, strjoin(Sc.kinds, ', '), mat2str(Sc.sd, 3));
all_passed = check('K is capped at N', P_full.K == Nm && size(P_full.H, 1) == Nm) && all_passed;
all_passed = check('K = N: Ritz values equal eig(J), raw and scaled (1e-8 of the largest)', ...
    err_full < 1e-8 && err_full_s < 1e-8) && all_passed;
P1s = SRNNCellTypePairs.projected_jacobian(m.S_out, m.t_out, m.fs, rows6, K, params, ...
    struct('seed', 9, 'state_scale', Sc.d));
P1o = SRNNCellTypePairs.projected_jacobian(m.S_out, m.t_out, m.fs, rows6, K, params, ...
    struct('seed', 9, 'state_scale', ones(Nm, 1)));
moved = max(arrayfun(@(s) match_err(eig(P1s.H(:, :, s)), eig(P1.H(:, :, s))), 1:numel(rows6)));
fprintf('  K = %d: the scaling moves the Ritz values by up to %.3g 1/s\n', K, moved);
all_passed = check('K < N: the scaling changes the Ritz values, and is flagged', moved > 1e-6 && P1s.scaled && ~P1.scaled) && all_passed;
all_passed = check('a scaling of ones is no scaling, to the bit', isequal(P1o.H, P1.H)) && all_passed;
all_passed = check('a non-positive or mis-sized scale errors :BadStateScale', ...
    throws_id(@() SRNNCellTypePairs.projected_jacobian(m.S_out, m.t_out, m.fs, rows6, K, params, ...
        struct('state_scale', zeros(Nm, 1))), 'SRNNCellTypePairs:BadStateScale') && ...
    throws_id(@() SRNNCellTypePairs.projected_jacobian(m.S_out, m.t_out, m.fs, rows6, K, params, ...
        struct('state_scale', ones(3, 1))), 'SRNNCellTypePairs:BadStateScale')) && all_passed;

%% 7. Scales by kind of variable
L = params.state_layout; win = rows6(1):numel(m.t_out);
all_passed = check('kinds on this network: x, a1, a2, a3, b1, b2', ...
    isequal(Sc.kinds, {'x', 'a1', 'a2', 'a3', 'b1', 'b2'}) && isempty(Sc.constant)) && all_passed;
vx = m.S_out(win, L.x);
a2 = [];
for q = 1:params.n_cellTypes
    blk = reshape(L.a{q}, params.n_per_type(q), []);
    a2 = [a2, m.S_out(win, blk(:, 2))]; %#ok<AGROW>
end
all_passed = check('x and a2 scales are the pooled standard deviations, written to every variable of the kind', ...
    abs(Sc.sd(1) - sqrt(mean(var(vx)))) < 1e-15 && abs(Sc.sd(3) - sqrt(mean(var(a2)))) < 1e-15 && ...
    all(Sc.d(L.x) == Sc.sd(1)) && numel(Sc.d) == Nm && all(Sc.d > 0)) && all_passed;
S_rest = SRNNCellTypePairs.state_type_scales(repmat(m.S_out(end, :), 3, 1), 1:3, params);
all_passed = check('a state that does not fluctuate: every scale 1, every kind named', ...
    all(S_rest.d == 1) && isequal(S_rest.constant, Sc.kinds)) && all_passed;

%% 8. Residuals, and a rotating non-normal system with a known answer
res_full = max(arrayfun(@(s) getfield(ritz_diagnostics(r_asm.H_t(:, :, s), r_asm.H_gram_t(:, :, s)), 'residual'), 1:numel(idx))); %#ok<GFLD>
res_inv  = max(arrayfun(@(s) getfield(ritz_diagnostics(r2.H_t(:, :, s), r2.H_gram_t(:, :, s)), 'residual'), 1:numel(idx2))); %#ok<GFLD>
pair_inv = max(arrayfun(@(s) max(getfield(ritz_diagnostics(r2.H_t(:, :, s), r2.H_gram_t(:, :, s)), 'pair_residual')), 1:numel(idx2))); %#ok<GFLD>
fprintf('  residual: K = N %.1e; invariant subspace (non-normal J) %.1e, worst Ritz pair %.1e\n', res_full, res_inv, pair_inv);
all_passed = check('residual is zero at K = N and on an invariant subspace (1e-6)', ...
    res_full < 1e-6 && res_inv < 1e-6 && pair_inv < 1e-6) && all_passed;

rng_state = rng; rng(21, 'twister');
A_rot = schur_form([-0.2, -0.5, -1, -6, -7, -8], 0.8);          % stable, non-normal
W = randn(6); Om = 0.4 * (W - W') / norm(W - W');               % skew: the frame turns at up to 0.4 rad/s
rng(rng_state);
B_rot = A_rot - Om;
[Us, Ts] = schur(B_rot, 'real');
ev_B = ordeig(Ts);
[~, by_re] = sort(real(ev_B), 'descend');
sel = false(6, 1); sel(by_re(1:3)) = true;
gap_B = min(real(ev_B(sel))) - max(real(ev_B(~sel)));
Us = ordschur(Us, Ts, sel);
V = Us(:, 1:3);                                                 % leading invariant subspace of B
want_ritz = eig(V' * A_rot * V);
want_res  = norm(A_rot * V - V * (V' * A_rot * V), 'fro') / norm(A_rot * V, 'fro');
fs_r = 1000; t_r = (0:1 / fs_r:12)'; X_r = [t_r, zeros(numel(t_r), 5)];   % the "state" carries the time
jac_rot = @(S, p) expm(Om * S(1)) * A_rot * expm(-Om * S(1));
idx_r = [10003, 10511, 11250, numel(t_r)];
o_r = struct('seed', 2, 'grid_fn', grid_fn, 'H_idx', idx_r);
r_rot = lyapunov_topk(X_r, t_r, fs_r, 0.05, [10 12], 10, 3, jac_rot, [], o_r);
err_ritz = 0; err_res = 0; err_eigA = Inf; err_eigB = Inf;
top3 = @(M) maxk_real(eig(M), 3);
for s = 1:numel(idx_r)
    Dg = ritz_diagnostics(r_rot.H_t(:, :, s), r_rot.H_gram_t(:, :, s));
    err_ritz = max(err_ritz, match_err(Dg.values, want_ritz));
    err_res  = max(err_res, abs(Dg.residual - want_res));
    err_eigA = min(err_eigA, match_err(Dg.values, top3(A_rot)));
    err_eigB = min(err_eigB, match_err(Dg.values, top3(B_rot)));
end
err_le = abs(sum(r_rot.LE_spectrum) - sum(real(ev_B(sel))));   % the sum: single exponents need longer than 2 s
fprintf(['  rotating system (gap %.2f 1/s): |Ritz - eig(V''AV)| %.1e; residual %.4f against %.4f (error %.1e);\n' ...
    '    the Ritz values are %.2f from the leading eigenvalues of A and %.2f from those of B; sum of the exponents against that of Re eig(B) %.1e\n'], ...
    gap_B, err_ritz, Dg.residual, want_res, err_res, err_eigA, err_eigB, err_le);
all_passed = check('rotating non-normal system: Ritz values equal eig(V'' A V) (1e-5) and the residual its closed form (1e-5)', ...
    gap_B > 1 && err_ritz < 1e-5 && err_res < 1e-5) && all_passed;
all_passed = check('rotating system: a nonzero residual, and Ritz values that are eigenvalues of neither A nor B', ...
    want_res > 0.05 && err_eigA > 1e-2 && err_eigB > 1e-2) && all_passed;
all_passed = check('rotating system: the Lyapunov exponents sum to the real parts of the leading eig(A - Omega) (1e-4)', err_le < 1e-4) && all_passed;

% (c) on the model: G and the residuals against direct computation
err_G = 0; err_res_m = 0; err_pair_m = 0; bound_ok = true;
for s = [1, 7, round(numel(rows) / 2), numel(rows), numel(rows) + 1]
    Qs = rH.Q_t(:, :, s); Hs = rH.H_t(:, :, s);
    Jd = full(jf(m.S_out(rH.H_idx(s), :)', params));
    JQ = Jd * Qs;
    err_G = max(err_G, max(abs(rH.H_gram_t(:, :, s) - JQ' * JQ), [], 'all') / max(abs(JQ' * JQ), [], 'all'));
    Dg = ritz_diagnostics(Hs, rH.H_gram_t(:, :, s));
    err_res_m = max(err_res_m, abs(Dg.residual - norm(JQ - Qs * Hs, 'fro') / norm(JQ, 'fro')));
    [Yv, Tv] = eig(Hs);
    direct = vecnorm(JQ * Yv - Qs * Yv * Tv) ./ vecnorm(JQ * Yv);
    err_pair_m = max(err_pair_m, match_err(sort(Dg.pair_residual), sort(direct(:))));
    omega_J = max(eig((Jd + Jd') / 2));
    bound_ok = bound_ok && Dg.max_re <= Dg.num_abscissa + 1e-12 && Dg.num_abscissa <= omega_J + 1e-9 && Dg.trace_err < 1e-10;
end
fprintf('  model: G against (JQ)''(JQ) %.1e; residual from H and G against the direct one %.1e (last state: %.3f); pairs %.1e\n', ...
    err_G, err_res_m, Dg.residual, err_pair_m);
all_passed = check('model: G = (J Q)'' (J Q) (1e-12); residual and pair residuals from H and G equal the direct ones (1e-7)', ...
    err_G < 1e-12 && err_res_m < 1e-7 && err_pair_m < 1e-7) && all_passed;
all_passed = check('model: max Re Ritz <= projected numerical abscissa <= numerical abscissa of J; trace kept', bound_ok) && all_passed;

fprintf('\n');
if all_passed
    fprintf('=== ALL lyapunov_topk projected-Jacobian TESTS PASSED ===\n');
else
    fprintf('=== SOME lyapunov_topk projected-Jacobian TESTS FAILED ===\n');
end

%% ------------------------------------------------------------------------
function T = schur_form(ev, coupling)
% Real Schur form with the given eigenvalues in order (a complex entry
% stands for its conjugate pair, a 2 x 2 block) and Gaussian entries of
% standard deviation `coupling` above the blocks (0: a normal matrix).
n = sum(1 + (imag(ev) ~= 0));
T = triu(coupling * randn(n), 1);
c = 0;
for k = 1:numel(ev)
    if imag(ev(k)) == 0
        c = c + 1; T(c, c) = real(ev(k));
    else
        T(c + 1:c + 2, c + 1:c + 2) = [real(ev(k)), imag(ev(k)); -imag(ev(k)), real(ev(k))];
        c = c + 2;
    end
end
end

function v = maxk_real(ev, k)
% The k eigenvalues of largest real part.
[~, order] = sort(real(ev), 'descend');
v = ev(order(1:k));
end

function e = match_err(a, b)
% Largest distance after pairing two sets of complex numbers greedily.
a = a(:); b = b(:); e = 0;
if numel(a) ~= numel(b); e = Inf; return; end
for k = 1:numel(a)
    [dmin, j] = min(abs(b - a(k)));
    e = max(e, dmin); b(j) = [];
end
end

function tf = same_but_seconds(a, b)
% Every field of a is in b and identical there, wall time excepted.
tf = true;
names = setdiff(fieldnames(a), {'seconds'});
for k = 1:numel(names)
    tf = tf && isfield(b, names{k}) && isequaln(a.(names{k}), b.(names{k}));
end
end

function [idx_all, t_lya, acc_start] = simple_grid(t, dt, decimation, tau, interval, lya_warmup)
% The segment grid of the model classes' lyapunov_sample_grid, without the
% clamp warning: start lya_warmup before the window, whole segments only.
acc_start = max(0, interval(1));
i_start = max(1, round((acc_start - lya_warmup - t(1)) / dt) + 1);
idx_all = i_start:decimation:numel(t);
t_lya = t(idx_all); t_lya = t_lya(:);
keep = (idx_all + decimation <= numel(t)) & (t_lya(:).' + tau <= min(t(end), interval(2)) + dt * 1e-6);
idx_all = idx_all(keep); t_lya = t_lya(keep);
end

function ok = throws_id(fn, id)
ok = false;
try
    fn();
catch ME
    ok = strcmp(ME.identifier, id);
end
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
