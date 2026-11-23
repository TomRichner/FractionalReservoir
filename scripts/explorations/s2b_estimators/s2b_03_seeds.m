% S2B_03_SEEDS  Other networks. For each seed pair and both adapting
% conditions (30 neurons; and 60 neurons, the size of the stage's production
% mode): the stage's window plus 56 s more. In the stage's window [-4, 4]:
% how much the local-rate trace depends on the vector (QR column 1 = e1
% against six random vectors evolved alone), whether the vectors have found a
% common direction, and what the trajectory is doing (decaying to a fixed
% point, or sustained).
[ddir] = s2b_out_dirs();
conds = {'sfa1_std1', 'sfa3_stdc2'};
cases = [repmat(30, 6, 1), (1:6)', (2:7)'; 60, 1, 2; 60, 2, 3];
rows = struct('cname', {}, 'n', {}, 'seeds', {}, 'N', {}, 'std_q1', {}, 'std_v', {}, 'rms_v_q1', {}, ...
    'cos_pair', {}, 'speed4', {}, 'speed60', {}, 'lam', {}, 'eig_end', {}, 't', {}, 'loc_q1', {}, 'loc_v', {}, 'speed', {});
for c = 1:2
    for k = 1:size(cases, 1)
        n = cases(k, 1);  seeds = cases(k, 2:3);
        M = s2b_build(conds{c}, seeds, 'n', n, 'T_range', [-8 60]);
        rng(7, 'twister');
        V0 = randn(M.N, 6);
        Rw = s2b_tangent(M, 0.1, 'K', 1, 'V0', V0, 't_stop', 4);     % the stage's window
        R  = s2b_tangent(M, 0.1, 'K', 6, 'V0', zeros(M.N, 0));         % to 60 s, six QR vectors
        m = Rw.t >= -4;
        r = struct('cname', conds{c}, 'n', n, 'seeds', seeds, 'N', M.N);
        r.std_q1 = std(Rw.loc_qr(m, 1));
        r.std_v  = std(Rw.loc_v(m, :));
        r.rms_v_q1 = rms(Rw.loc_v(m, :) - Rw.loc_qr(m, 1));
        G = abs(Rw.V_end' * Rw.V_end);  r.cos_pair = G(triu(true(6), 1))';
        sp = @(tq) norm(M.rhs(tq, M.S(1 + round((tq - M.t(1)) * M.fs), :)'));
        r.speed4 = sp(4);  r.speed60 = sp(60);
        r.lam = sort(mean(R.loc_qr(R.t >= 30, :)), 'descend');       % six leading, over [30, 60]
        ev = eig(full(SRNNCellTypePairs.compute_Jacobian_fast(M.S(end, :)', M.params)));
        r.eig_end = sort(real(ev), 'descend');  r.eig_end = r.eig_end(1:6)';
        r.t = Rw.t;  r.loc_q1 = Rw.loc_qr(:, 1);  r.loc_v = Rw.loc_v;
        idx = 1 + round((R.t - M.t(1)) * M.fs);
        r.speed = arrayfun(@(i) norm(M.rhs(M.t(i), M.S(i, :)')), idx);
        rows(end + 1) = r; %#ok<SAGROW>
        fprintf('%-10s n=%2d seeds [%d %d]: |dS/dt| at 4 s %.1e, at 60 s %.1e | std q1 %.3f | std random %.3f..%.3f | rms(random - q1) median %.3f | pair |cos| median %.2f | lam1-lam2 %.4f, lam1-lam6 %.4f (lam1 %+.3f)\n', ...
            conds{c}, n, seeds, r.speed4, r.speed60, r.std_q1, min(r.std_v), max(r.std_v), median(r.rms_v_q1), median(r.cos_pair), ...
            r.lam(1) - r.lam(2), r.lam(1) - r.lam(6), r.lam(1));
        clear M R Rw
    end
end
save(fullfile(ddir, 's2b_03_seeds.mat'), 'rows');
