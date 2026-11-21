function s = bm_detect_bursts(t, mE, def)
%BM_DETECT_BURSTS Population bursts: the mean E rate above a threshold for long enough, close crossings merged.
%
%   s = bm_detect_bursts(t, mE, def)
%
% A BURST is the population mean E rate mE above def.thresh for at least
% def.min_dur s, threshold crossings separated by less than def.merge s
% being one burst. The rule of the block-B bursting pilots (pilot 0b, then
% p1b_detect_bursts in scripts/explorations/block_B_pilots/p1b_stim_bursting,
% which the stimulation pilot of 2026-11-21 used); ported here unchanged so
% that the burst_mode stage does not call into an exploration folder.
%
% t    1 x nt (or nt x 1), s, uniform
% mE   the same size, population mean E rate
% def  struct: thresh, min_dur, merge (s)
%
% s.n, s.t_on, s.t_off, s.duration (1 x nB, s), s.peak_t, s.peak (mean E
% rate at the peak).
%
% See also: burst_mode_stage_plan, bm_crossings

t = t(:)'; mE = mE(:)';
dt = t(2) - t(1);
above = mE > def.thresh;
d = diff([false, above, false]);
on = find(d == 1); off = find(d == -1) - 1;
k = 1;
while k < numel(on)
    if (on(k + 1) - off(k) - 1) * dt < def.merge
        off(k) = off(k + 1); on(k + 1) = []; off(k + 1) = [];
    else
        k = k + 1;
    end
end
dur = (off - on + 1) * dt;
ok = dur >= def.min_dur;
on = on(ok); off = off(ok); dur = dur(ok);
nB = numel(on);
peak_t = nan(1, nB); peak = nan(1, nB);
for b = 1:nB
    [peak(b), i] = max(mE(on(b):off(b)));
    peak_t(b) = t(on(b) + i - 1);
end
s = struct('n', nB, 't_on', t(on), 't_off', t(off), 'duration', dur, ...
    'peak_t', peak_t, 'peak', peak);
end
