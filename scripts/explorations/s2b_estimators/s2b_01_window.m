% S2B_01_WINDOW  The stage's own window (T = [-8, 4] s, reduced network of 30
% neurons, seeds [1 2], no noise, ode45), both adapting conditions.
%  1. What panel B draws: the saved traces against a re-derivation (RK4 on
%     the variational equation): which raw QR column is the plotted "QR local
%     lambda_1", and is it the first vector?
%  2. Initial vector: six random tracers, the flow direction and e1, each
%     evolved alone (Benettin, infinitesimal). Do they give one trace?
%  3. Precision: the class's finite-perturbation Benettin (d0 = 1e-3, ode45)
%     against the linear tracer from the same initial vector; from e1 it must
%     reproduce the QR first column.
%  4. Segment length 5, 10, 20, 40, 100 ms; half the step (trajectory on an
%     800 Hz grid, RK4 step 1/400 instead of 1/200); an ode45 (1e-10)
%     integration of the variational equation for the first vector as the
%     reference.
[ddir, ~, root] = s2b_out_dirs();
D = load(fullfile(root, 'data', 'mu7chain_AC_fast', 'numerics_verification', 'numerics_verification_data.mat'));
conds = {'sfa1_std1', 'sfa3_stdc2'};
out = struct();
for c = 1:2
    cname = conds{c};
    Qs = D.results(strcmp(D.cond_names, cname)).qr(1);
    M = s2b_build(cname, [1 2]);
    rng(7, 'twister');
    V0 = [randn(M.N, 6), M.rhs(M.t(1), M.S(1, :)'), [1; zeros(M.N - 1, 1)]];
    o = struct('cname', cname, 'N', M.N, 'saved', Qs);
    % 1 and 2
    o.R100 = s2b_tangent(M, 0.1, 'V0', V0);          % full QR + tracers, the stage's QR segment
    o.R020 = s2b_tangent(M, 0.02, 'K', 3, 'V0', V0); % the stage's Benettin segment
    % 4, segment
    segs = [0.005 0.01 0.02 0.04 0.1];
    o.segs = segs;  o.Rseg = cell(1, numel(segs));
    for k = 1:numel(segs); o.Rseg{k} = s2b_tangent(M, segs(k), 'K', 1, 'V0', V0(:, 1)); end
    % 3, finite-perturbation Benettin as in the class
    [o.fd_t, o.fd_rand] = s2b_benettin_fd(M, 0.02, V0(:, 1));
    [~,      o.fd_e1]   = s2b_benettin_fd(M, 0.02, V0(:, 8));
    [~,      o.fd_e1_small] = s2b_benettin_fd(M, 0.02, V0(:, 8), 1e-5);
    % 4, half step
    M8 = s2b_build(cname, [1 2], 'fs', 800);
    o.R100_fs800 = s2b_tangent(M8, 0.1, 'K', 1, 'V0', V0(:, 1));
    % 4, reference: ode45 on the variational equation of ONE vector (e1),
    % states by pchip between the 800 Hz samples, tolerance 1e-10
    F = griddedInterpolant(M8.t, M8.S, 'pchip');
    v = V0(:, 8);  tt = M8.t(1):0.1:M8.t(end);  ref = zeros(numel(tt) - 1, 1);
    opt = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
    for k = 1:numel(tt) - 1
        [~, Y] = ode45(@(t, y) SRNNCellTypePairs.compute_Jacobian_fast(F(t)', M8.params) * y, tt(k:k + 1), v, opt);
        y = Y(end, :)';  ref(k) = log(norm(y)) / 0.1;  v = y / norm(y);
    end
    o.ref_t = tt(1:end - 1)';  o.ref_e1 = ref;
    out.(cname) = o;
    clear M M8 F
end
save(fullfile(ddir, 's2b_01_window.mat'), '-struct', 'out');

%% Print
w = @(t) t >= -4;
for c = 1:2
    o = out.(conds{c});  R = o.R100;  m = w(R.t);
    fprintf('\n=== %s, N = %d states ===\n', o.cname, o.N);
    [~, ord] = sort(mean(R.loc_qr(R.t >= 0, :)), 'descend');
    fprintf('1. raw QR column behind the plotted (sorted) lambda_1: %d; rms(saved - rederived) over t >= -4: %.2g s^-1\n', ...
        ord(1), rms(R.loc_qr(m, ord(1)) - o.saved.qr.local_LE_spectrum_t(m, 1)));
    sp = sort(mean(R.loc_qr(R.t >= 0, :)), 'descend');
    fprintf('   finite-time spectrum over [0, 4], first 8: %s\n', mat2str(sp(1:8), 3));
    fprintf('   saved Benettin: std over t >= -4 = %.3f; QR column 1: %.3f\n', ...
        std(o.saved.benettin.local_lya(o.saved.benettin.t_lya >= -4)), std(R.loc_qr(m, 1)));
    fprintf('2. std over t >= -4 of the local rate of each vector evolved alone (6 random | flow | e1): %s | %.3f | %.3f\n', ...
        mat2str(std(R.loc_v(m, 1:6)), 3), std(R.loc_v(m, 7)), std(R.loc_v(m, 8)));
    G = abs(R.V_end(:, 1:6)' * R.V_end(:, 1:6));  G = G(triu(true(6), 1));
    fprintf('   |cos| between the six random vectors at t = 4: min %.2f median %.2f max %.2f\n', min(G), median(G), max(G));
    fprintf('   |cos| with the flow direction at t = 4: random %s, QR q1 %.2f\n', mat2str(R.v_flow(end, 1:6), 2), R.q1_flow(end));
    fprintf('   norm of the projection of the random vectors on span(q1..qk) at t = 4, k = %s:\n', mat2str(R.ks));
    disp(round(squeeze(R.proj_v(end, 1:6, :)), 2));
    fprintf('   rotation of q1 per 0.1 s (rad), t >= -4: median %.3f max %.3f\n', median(R.rot_q1(m)), max(R.rot_q1(m)));
    m2 = w(o.fd_t);
    fprintf('3. finite (d0 = 1e-3) against linear, same initial vector, seg 0.02, t >= -4: rms difference random %.3g, e1 %.3g; e1 with d0 = 1e-5: %.3g s^-1\n', ...
        rms(o.fd_rand(m2) - o.R020.loc_v(m2, 1)), rms(o.fd_e1(m2) - o.R020.loc_v(m2, 8)), rms(o.fd_e1_small(m2) - o.R020.loc_v(m2, 8)));
    fprintf('   finite Benettin from e1: std %.3f (QR column 1 at 0.02 s: %.3f)\n', std(o.fd_e1(m2)), std(o.R020.loc_qr(m2, 1)));
    fprintf('4. segment (ms) -> std of the rate of random vector 1, t >= -4: ');
    for k = 1:numel(o.segs); fprintf('%g: %.3f  ', 1e3 * o.segs(k), std(o.Rseg{k}.loc_v(w(o.Rseg{k}.t), 1))); end
    fprintf('\n   segment (ms) -> std of QR column 1 (= e1): ');
    for k = 1:numel(o.segs); fprintf('%g: %.3f  ', 1e3 * o.segs(k), std(o.Rseg{k}.loc_qr(w(o.Rseg{k}.t), 1))); end
    fprintf('\n   half step: max |difference| of the 0.1 s rates, e1 %.2g, random %.2g s^-1\n', ...
        max(abs(o.R100_fs800.loc_qr(:, 1) - R.loc_qr(:, 1))), max(abs(o.R100_fs800.loc_v(:, 1) - R.loc_v(:, 1))));
    fprintf('   ode45 (1e-10) variational reference for e1: max |difference| from QR column 1 %.2g s^-1 (after t = -7: %.2g)\n', ...
        max(abs(o.ref_e1 - R.loc_qr(:, 1))), max(abs(o.ref_e1(o.ref_t >= -7) - R.loc_qr(R.t >= -7, 1))));
end
