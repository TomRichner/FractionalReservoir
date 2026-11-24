% S2B_04_WARMUP  Would a longer warm-up make the two estimators agree? Only
% where the activity is sustained (a network that decays to a fixed point has
% nothing left to compare after the warm-up). Multiple-timescale networks
% with sustained activity: 60 neurons seeds [1 2] and [2 3], 30 neurons seeds
% [4 5]; and single-timescale 60 neurons seeds [2 3] for contrast. Run to
% 120 s; in 8-s windows: the rms difference between the local rate of each
% of six random vectors and of the QR first vector (0.1 s segments), next to
% the swing (std) of the QR rate itself.
[ddir] = s2b_out_dirs();
cs = {'sfa3_stdc2', 60, [1 2]; 'sfa3_stdc2', 60, [2 3]; 'sfa3_stdc2', 30, [4 5]; 'sfa1_std1', 60, [2 3]};
out = struct('cname', {}, 'n', {}, 'seeds', {}, 't', {}, 'loc_q1', {}, 'loc_v', {}, 'proj_q1', {}, 'w0', {}, 'rms_w', {}, 'std_w', {}, 'lam', {});
w0 = -8:8:104;
for k = 1:size(cs, 1)
    M = s2b_build(cs{k, 1}, cs{k, 3}, 'n', cs{k, 2}, 'T_range', [-8 120]);
    rng(7, 'twister');
    R = s2b_tangent(M, 0.1, 'K', 6, 'V0', randn(M.N, 6));
    o = struct('cname', cs{k, 1}, 'n', cs{k, 2}, 'seeds', cs{k, 3}, 't', R.t, 'loc_q1', R.loc_qr(:, 1), ...
        'loc_v', R.loc_v, 'proj_q1', R.proj_v(:, :, 1), 'w0', w0, 'rms_w', [], 'std_w', [], 'lam', mean(R.loc_qr(R.t >= 60, :)));
    for j = 1:numel(w0)
        m = R.t >= w0(j) & R.t < w0(j) + 8;
        o.rms_w(j, :) = rms(R.loc_v(m, :) - R.loc_qr(m, 1));
        o.std_w(j) = std(R.loc_qr(m, 1));
    end
    out(k) = o;
    fprintf('\n%s n = %d seeds [%d %d]; six leading exponents over [60, 120]: %s\n', o.cname, o.n, o.seeds, mat2str(sort(o.lam, 'descend'), 3));
    fprintf('window start (s):          %s\n', sprintf('%6d', w0));
    fprintf('std of QR rate:            %s\n', sprintf('%6.2f', o.std_w));
    fprintf('rms(random - QR), median:  %s\n', sprintf('%6.2f', median(o.rms_w, 2)));
    fprintf('rms(random - QR), max:     %s\n', sprintf('%6.2f', max(o.rms_w, [], 2)));
    fprintf('|cos|(random, q1), min:    %s\n', sprintf('%6.2f', arrayfun(@(s) min(o.proj_q1(find(R.t < s + 8, 1, 'last'), :)), w0)));
    clear M R
end
save(fullfile(ddir, 's2b_04_warmup.mat'), 'out');
