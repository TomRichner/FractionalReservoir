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
