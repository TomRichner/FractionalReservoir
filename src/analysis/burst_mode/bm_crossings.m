function M = bm_crossings(t, z, bursts, windows, tol_pre, tol_post)
%BM_CROSSINGS How a mode's real part (1/s) lines up with population bursts: crossings of the boundary, their lead, the share of time beyond it.
%
%   M = bm_crossings(t, z, bursts, windows, tol_pre, tol_post)
%
% t        sample times (s), uniform
% z        the mode's rate at each sample (complex or real, 1/s); the
%          boundary is Re z = 0
% bursts   bm_detect_bursts' struct (t_on, t_off)
% windows  struct of named [t0 t1] windows (s), e.g. before, during
% tol_pre, tol_post   s: a crossing COINCIDES with a burst when it overlaps
%          [burst on - tol_pre, burst off + tol_post]
%
% A CROSSING is a maximal run of consecutive samples with Re z > 0. Counts
% are over each window (a burst or a crossing counts where it starts).
% Per window:
%   frac_beyond          share of samples with Re z > 0
%   mean_re, max_re, min_re
%   n_bursts, n_cross
%   n_bursts_with_cross  bursts with a coinciding crossing
%   n_cross_in_burst     crossings that coincide with a burst
%   lead                 per burst, burst onset minus the onset of its
%                        first coinciding crossing (s; positive = the
%                        crossing came first); NaN without one
%   median_lead
% Ported from the pilot's so_crossings (2026-11-21), on rates instead of
% eigenvalues of M (Re lambda(M) > 1 there is Re z > 0 here), and with the
% lead signed as a lead (the pilot's lag was crossing minus burst onset).
%
% See also: bm_run_network, bm_detect_bursts

t = t(:); re = real(z(:));
above = re > 0;
d = diff([false; above; false]);
% d has one more element than t: an onset at k is d(k) = 1, an offset at k is d(k + 1) = -1
c_on = t(d(1:end - 1) == 1); c_off = t(d(2:end) == -1);
b_on = bursts.t_on(:); b_off = bursts.t_off(:);
names = fieldnames(windows);
M = struct('cross_on', c_on', 'cross_off', c_off');
for w = 1:numel(names)
    lim = windows.(names{w});
    k = t >= lim(1) & t < lim(2);
    r = struct('window', lim, 'frac_beyond', mean(above(k)), 'mean_re', mean(re(k)), ...
        'max_re', max(re(k)), 'min_re', min(re(k)));
    bsel = find(b_on >= lim(1) & b_on < lim(2));
    csel = find(c_on >= lim(1) & c_on < lim(2));
    r.n_bursts = numel(bsel);
    r.n_cross = numel(csel);
    lead = nan(1, numel(bsel));
    for i = 1:numel(bsel)
        b = bsel(i);
        ov = find(c_off >= b_on(b) - tol_pre & c_on <= b_off(b) + tol_post, 1);
        if ~isempty(ov); lead(i) = b_on(b) - c_on(ov); end
    end
    hit_c = false(1, numel(csel));
    for i = 1:numel(csel)
        c = csel(i);
        hit_c(i) = any(c_off(c) >= b_on - tol_pre & c_on(c) <= b_off + tol_post);
    end
    r.n_bursts_with_cross = nnz(isfinite(lead));
    r.n_cross_in_burst = nnz(hit_c);
    r.lead = lead;
    r.median_lead = median(lead, 'omitnan');
    M.(names{w}) = r;
end
end
