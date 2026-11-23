function M = so_crossings(Q, S)
% SO_CROSSINGS Fractions of time beyond Re = 1, crossings, and how they line up with bursts.
%
%   M = so_crossings(Q, S)
%
% Q is a spectral series of so_run_network (fields t, lam1, lam_coh,
% lam_ref, R_eff, dbar, rate, bursts); S the settings (win_before,
% win_during, t_stim, tol_pre, tol_post).
%
% A CROSSING is a maximal run of consecutive samples with Re lam1 > 1 (the
% frozen-gain rate (Re lam1 - 1)/tau_d positive). A burst and a crossing
% COINCIDE when the crossing overlaps [burst on - tol_pre, burst off +
% tol_post]. Counts are over the baseline window win_before (a burst or
% crossing counts if it starts there) and the stimulus window win_during.
%
% Per window: frac_beyond (fraction of samples with Re lam1 > 1),
% frac_beyond_coh (the same for the most coherent mode), frac_beyond_ref
% (the followed structural mode), mean and max Re lam1, mean R_eff, mean
% gains dbar (2 x 2, post x pre), mean rates, n_bursts, n_cross,
% n_bursts_with_cross, n_cross_in_burst, lag (crossing onset - burst
% onset, s, for the first coinciding crossing of each burst), and the
% fraction beyond inside bursts and between them.

t = Q.t(:);
re1 = real(Q.lam1(:));
above = re1 > 1;
d = diff([false; above; false]);
c_on = t(find(d == 1)); c_off = t(find(d == -1) - 1);
b_on = Q.bursts.t_on(:); b_off = Q.bursts.t_off(:);

in_burst = false(size(t));
for b = 1:numel(b_on)
    in_burst = in_burst | (t >= b_on(b) & t <= b_off(b));
end

wins = struct('before', S.win_before, 'during', S.win_during, ...
    'transient', [S.t_stim, S.win_during(1)]);
names = fieldnames(wins);
M = struct('cross_on', c_on', 'cross_off', c_off');
for w = 1:numel(names)
    lim = wins.(names{w});
    k = t >= lim(1) & t < lim(2);
    r = struct();
    r.window = lim;
    r.frac_beyond = mean(above(k));
    r.frac_beyond_coh = mean(real(Q.lam_coh(k)) > 1);
    r.frac_beyond_ref = mean(real(Q.lam_ref(k)) > 1);
    r.mean_re1 = mean(re1(k));
    r.max_re1 = max(re1(k));
    r.min_re1 = min(re1(k));
    r.mean_coh1 = mean(Q.coh1(k));
    r.mean_R_eff = mean(Q.R_eff(k));
    r.max_R_eff = max(Q.R_eff(k));
    r.dbar = squeeze(mean(Q.dbar(k, :, :), 1));
    r.rate = mean(Q.rate(k, :), 1);
    r.frac_beyond_in_burst = mean(above(k & in_burst));
    r.frac_beyond_between = mean(above(k & ~in_burst));
    bsel = find(b_on >= lim(1) & b_on < lim(2));
    csel = find(c_on >= lim(1) & c_on < lim(2));
    r.n_bursts = numel(bsel);
    r.n_cross = numel(csel);
    hit_b = false(numel(bsel), 1); lag = nan(numel(bsel), 1);
    for i = 1:numel(bsel)
        b = bsel(i);
        ov = find(c_off >= b_on(b) - S.tol_pre & c_on <= b_off(b) + S.tol_post, 1);
        if ~isempty(ov)
            hit_b(i) = true;
            lag(i) = c_on(ov) - b_on(b);
        end
    end
    hit_c = false(numel(csel), 1);
    for i = 1:numel(csel)
        c = csel(i);
        hit_c(i) = any(c_off(c) >= b_on - S.tol_pre & c_on(c) <= b_off + S.tol_post);
    end
    r.n_bursts_with_cross = sum(hit_b);
    r.n_cross_in_burst = sum(hit_c);
    r.lag = lag';
    r.median_lag = median(lag, 'omitnan');
    r.cross_dur_median = median(c_off(csel) - c_on(csel));
    M.(names{w}) = r;
end
end
