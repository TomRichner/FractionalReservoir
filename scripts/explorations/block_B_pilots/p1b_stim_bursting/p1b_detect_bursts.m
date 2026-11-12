function s = p1b_detect_bursts(t, mE, def)
% P1B_DETECT_BURSTS Pilot 0b's burst definition, as a function.
%
%   s = P1B_DETECT_BURSTS(t, mE, def)
%
% A BURST is the population mean E rate mE above def.thresh for at least
% def.min_dur s, crossings separated by less than def.merge s being one.
% This is the definition of the local function detect_bursts in
% scripts/examples/bursting_chain_mts_example/bursting_chain_mts_figure.m
% (pilot 0b), which is not callable from outside that file; the arithmetic
% here is the same, line for line, without the burst-size part.
%
% t    1 x nt, s, uniform
% mE   1 x nt, population mean E rate
% def  struct: thresh, min_dur, merge (s)
%
% s.t_on, s.t_off, s.duration (1 x nB, s), s.peak_t, s.peak (mean E rate
% at the peak), s.n.

dt    = t(2) - t(1);
above = mE > def.thresh;
d   = diff([false, above, false]);
on  = find(d == 1); off = find(d == -1) - 1;
k = 1;
while k < numel(on)
    if (on(k+1) - off(k) - 1) * dt < def.merge
        off(k) = off(k+1); on(k+1) = []; off(k+1) = [];
    else
        k = k + 1;
    end
end
dur = (off - on + 1) * dt;
ok  = dur >= def.min_dur;
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
