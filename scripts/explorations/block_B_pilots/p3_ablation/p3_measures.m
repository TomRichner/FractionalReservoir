function M = p3_measures(R)
% P3_MEASURES The measures of pilot 3 for one run.
%
%   M = P3_MEASURES(R)      R from p3_run_one (t, rE, rI, x, groups, network)
%
% Aberrant activity, by each network's own definition, unchanged from the
% unstimulated pilots:
%   bursting  a BURST is the population mean E rate above 0.2 for at least
%             0.1 s, dips shorter than 0.2 s not splitting it; its size is
%             the fraction of all neurons with rate above 0.1 at its peak
%             (pilot 0b, bursting_chain_mts_figure).
%   run13     an EVENT of the IED pilot's detector v1 (ied_detect_events,
%             ied_event_details); its size is the fraction of neurons
%             recruited, its duration the waveform width.
% The burst definition is also applied to run13 (n_bursts there counts
% network-wide excursions), because detector v1 sets its thresholds from
% the segment itself and so cannot say "continuously high".
%
% M.count, M.per_min        aberrant events (bursts or IED events)
% M.size_median             fraction of neurons (see above)
% M.duration_median         s
% M.n_bursts, M.burst_*     the burst definition, both networks
% M.frac_time_burst         fraction of time inside a burst (mean E > 0.2)
% M.mean_rE, M.mean_rI      mean rates
% M.peak_mean_rE            largest population mean E rate
% M.frac_time_E_above_0p5   fraction of time the population mean E rate > 0.5
% M.frac_pinned             fraction of neurons with time-averaged rate > 0.9
% M.frac_silent             fraction of neurons with time-averaged rate < 0.01
% M.events / M.bursts       the details
%
% See also: p3_run_one

burst_thresh = 0.2; burst_min_dur = 0.1; burst_merge = 0.2; active_thresh = 0.1;
t = R.t;
T = t(end) - t(1);
mE = mean(R.rE, 1);
Rall = [R.rE; R.rI];

% Bursts (the definition of bursting_chain_mts_figure's detect_bursts).
dt = t(2) - t(1);
above = mE > burst_thresh;
d = diff([false, above, false]);
on = find(d == 1); off = find(d == -1) - 1;
k = 1;
while k < numel(on)
    if (on(k+1) - off(k) - 1) * dt < burst_merge
        off(k) = off(k+1); on(k+1) = []; off(k+1) = [];
    else
        k = k + 1;
    end
end
dur = (off - on + 1) * dt;
ok = dur >= burst_min_dur;
on = on(ok); off = off(ok); dur = dur(ok);
nB = numel(on);
frac_peak = nan(1, nB); peak_t = nan(1, nB);
for b = 1:nB
    [~, i] = max(mE(on(b):off(b)));
    ip = on(b) + i - 1;
    peak_t(b) = t(ip);
    frac_peak(b) = mean(Rall(:, ip) > active_thresh);
end
M = struct();
M.n_bursts = nB;
M.burst_per_min = nB * 60 / T;
M.burst_duration_median = median(dur);
M.burst_size_median = median(frac_peak);
M.frac_time_burst = mean(above);
M.bursts = struct('t_on', t(on), 't_off', t(off), 'peak_t', peak_t, ...
    'duration', dur, 'frac_active_peak', frac_peak);

switch R.network
    case 'bursting'
        M.count = nB;
        M.size_median = M.burst_size_median;
        M.duration_median = M.burst_duration_median;
        M.event_times = peak_t;
        M.events = [];
        M.long_group_excursions = NaN;
    case 'run13'
        det = ied_detect_events(t, R.x, Rall, R.groups, struct());
        events = ied_event_details(struct('t', t, 'x', R.x, 'groups', R.groups), det);
        M.count = det.metrics.event_count;
        M.size_median = median(events.Recruitment);
        M.duration_median = median(events.WaveformWidth, 'omitnan');
        M.event_times = [det.events.time];
        M.events = events;
        M.long_group_excursions = det.metrics.long_group_excursions;
end
M.per_min = M.count * 60 / T;

M.mean_rE = mean(R.rE, 'all');
M.mean_rI = mean(R.rI, 'all');
M.peak_mean_rE = max(mE);
M.frac_time_E_above_0p5 = mean(mE > 0.5);
rbar = mean(Rall, 2);
M.frac_pinned = mean(rbar > 0.9);
M.frac_silent = mean(rbar < 0.01);
end
