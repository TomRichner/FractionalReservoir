% S2B_02_LONG  The same two networks (30 neurons, seeds [1 2]) run to 200 s.
% Where does the trajectory go after the stage's window, what are the leading
% exponents and their gaps, and how long do independent vectors take to
% agree on a direction?
[ddir] = s2b_out_dirs();
conds = {'sfa1_std1', 'sfa3_stdc2'};
out = struct();
for c = 1:2
    M = s2b_build(conds{c}, [1 2], 'T_range', [-8 200]);
    rng(7, 'twister');
    V0 = randn(M.N, 6);
    R = s2b_tangent(M, 0.1, 'V0', V0);
    o = struct('cname', conds{c}, 'N', M.N, 't', R.t, 'loc_qr', R.loc_qr, 'loc_v', R.loc_v, ...
        'proj_q1', R.proj_v(:, :, 1), 'q1_flow', R.q1_flow, 'v_flow', R.v_flow);
    idx = 1 + round((R.t - M.t(1)) * M.fs);
    o.speed = arrayfun(@(i) norm(M.rhs(M.t(i), M.S(i, :)')), idx);   % |dS/dt|
    o.x_ex = M.S(1:40:end, M.layout.x(1:3));  o.t_ex = M.t(1:40:end);
    % Jacobian at the last state: if the run ended on a fixed point, its
    % eigenvalues ARE the Lyapunov exponents (real parts)
    ev = eig(full(SRNNCellTypePairs.compute_Jacobian_fast(M.S(end, :)', M.params)));
    [~, k] = sort(real(ev), 'descend');  o.eig_end = ev(k);
    out.(conds{c}) = o;
    clear M R
end
save(fullfile(ddir, 's2b_02_long.mat'), '-struct', 'out');

for c = 1:2
    o = out.(conds{c});
    fprintf('\n=== %s ===\n', o.cname);
    fprintf('window      |dS/dt|    QR col1 mean  std    | random vectors: std of rate (min..max)   |cos| with q1 (min..max)   pairwise spread of the rate\n');
    for w0 = [-4 0 10 20 50 100 150 190]
        m = o.t >= w0 & o.t < w0 + 10;  if w0 == -4; m = o.t >= -4 & o.t < 4; end
        sv = std(o.loc_v(m, :));  pq = o.proj_q1(find(m, 1, 'last'), :);
        fprintf('[%4d,%4d)  %.2e   %+.4f      %.4f | %.4f .. %.4f                        %.2f .. %.2f              %.4f\n', w0, w0 + 10 - 2 * (w0 == -4), ...
            mean(o.speed(m)), mean(o.loc_qr(m, 1)), std(o.loc_qr(m, 1)), min(sv), max(sv), min(pq), max(pq), mean(std(o.loc_v(m, :), 0, 2)));
    end
    sp = sort(mean(o.loc_qr(o.t >= 100, :)), 'descend');
    fprintf('QR spectrum over [100, 200], first 8: %s\n', mat2str(sp(1:8), 3));
    fprintf('gaps lambda_1 - lambda_k, k = 2..6: %s\n', mat2str(sp(1) - sp(2:6), 2));
    fprintf('eigenvalues of J at t = 200, first 8 (real parts): %s\n', mat2str(real(o.eig_end(1:8))', 3));
    fprintf('   imaginary parts: %s\n', mat2str(imag(o.eig_end(1:8))', 3));
    fprintf('number of exponents within 0.05 s^-1 of lambda_1: %d of %d\n', nnz(sp >= sp(1) - 0.05), o.N);
end
