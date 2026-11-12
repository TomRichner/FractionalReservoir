function [A, companion, detail] = coactivation_excess(r, opts)
% COACTIVATION_EXCESS One measure of aberrant population activity, from per-neuron rates alone.
%
%   [A, companion, detail] = coactivation_excess(r)
%   [A, companion, detail] = coactivation_excess(r, Name, Value, ...)
%
% r   n x nt, the rate of every neuron (rows) over the analysis window
%     (columns), rates in [0, 1] as the model's sigmoids give them. No cell
%     types, modules or stimulus are needed. Pilot 4 used 100 Hz sampling
%     and windows of 35 to 100 s.
%
% THE MEASURE (block B pilot 4, 2026-11-12). Three steps:
%   1. Each neuron's own threshold: its median rate over the window plus
%      max(k * s_i, floor), where s_i = 1.4826 * (median absolute deviation
%      of its rate), the robust standard deviation. A neuron is "in
%      excursion" at time t when its rate is above that threshold.
%   2. Co-activation c(t): the fraction of all n neurons in excursion at t.
%   3. The chance level L = p + z * sqrt(p (1 - p) / n) + 1/n, where p is
%      the time average of c(t): what n INDEPENDENT neurons, each in
%      excursion a fraction p of the time, would reach only z standard
%      deviations out, and never fewer than two neurons together.
%   A = the time average of max(c(t) - L, 0): the fraction of the network
%   recruited into excursions BEYOND chance, averaged over time. It is 0
%   for a steady network, a silent one, a saturated one, and for neurons
%   that fluctuate independently; it grows with how often, how long and
%   how widely neurons leave their own baselines together.
%
% Name-value options (all dimensionless; the defaults are the pilot's)
%   k            3      robust standard deviations above the median
%   floor        0.1    least excursion, in rate units (a tenth of the
%                       rate's range); needed because a deterministic
%                       network between bursts has s_i near 0
%   z            5      standard deviations of the independent-neuron null
%   pinned_rate      0.9   a neuron whose time-mean rate is above this is pinned
%   saturated_frac   0.25  pinned fraction at or above which the network is 'saturated'
%   silent_rate      0.02  network mean rate below which it is 'silent'
%   aberrant_A       2.5e-4  A at or above which it is 'aberrant' (see the
%                       pilot's report for how unsure this is)
%
% companion   what must be read beside A, because A is 0 for a saturated
%             and for a silent network as it is for a healthy one:
%               mean_rate        mean over neurons and time
%               pinned_fraction  fraction of neurons with time-mean rate > pinned_rate
%               class            'silent' | 'saturated' | 'aberrant' | 'healthy',
%                                tested in that order
% detail      c (1 x nt), level L, p, threshold (n x 1), time_fraction
%             (share of the window with c > L), events_per_min (upward
%             crossings of L, crossings closer than 0.12 s merged; needs
%             opts.dt), q99 (99th percentile of c)
%
% WHAT IT CANNOT SEE: activity that is the network's majority state (a
% neuron above its own median more than half the time has that as its
% baseline), steady high firing (hence the companion), and events that
% recruit fewer neurons than the chance level L (about 5 % of a
% 350-neuron network here).
%
% Cost: two medians per neuron over the window; 0.1 s for 500 neurons x
% 10001 samples. It needs the window's rates in memory once (the model
% has them at the end of a run); nothing is stored.
%
% See also: run_p4_aberrant_measure, p4_candidates

arguments
    r double
    opts.k (1,1) double {mustBePositive} = 3
    opts.floor (1,1) double {mustBeNonnegative} = 0.1
    opts.z (1,1) double {mustBeNonnegative} = 5
    opts.pinned_rate (1,1) double = 0.9
    opts.saturated_frac (1,1) double = 0.25
    opts.silent_rate (1,1) double = 0.02
    opts.aberrant_A (1,1) double = 2.5e-4
    opts.dt double = []          % s per sample, for events_per_min only
end

n = size(r, 1);
med = median(r, 2);
s = 1.4826 * median(abs(r - med), 2);
threshold = med + max(opts.k * s, opts.floor);
c = mean(r > threshold, 1);
p = mean(c);
L = p + opts.z * sqrt(p * (1 - p) / n) + 1 / n;
A = mean(max(c - L, 0));

mean_per_neuron = mean(r, 2);
companion = struct('mean_rate', mean(mean_per_neuron), ...
    'pinned_fraction', mean(mean_per_neuron > opts.pinned_rate), 'class', '');
if companion.mean_rate < opts.silent_rate
    companion.class = 'silent';
elseif companion.pinned_fraction >= opts.saturated_frac
    companion.class = 'saturated';
elseif A >= opts.aberrant_A
    companion.class = 'aberrant';
else
    companion.class = 'healthy';
end

if nargout >= 3
    above = c > L;
    detail = struct('c', c, 'level', L, 'p', p, 'threshold', threshold, ...
        'time_fraction', mean(above), 'events_per_min', NaN, 'q99', prctile(c, 99), ...
        'event_on', [], 'event_off', []);
    if ~isempty(opts.dt)
        d = diff([false above false]);
        on = find(d == 1);  off = find(d == -1) - 1;
        j = 1;
        while j < numel(on)
            if (on(j + 1) - off(j) - 1) * opts.dt < 0.12
                off(j) = off(j + 1);  on(j + 1) = [];  off(j + 1) = [];
            else
                j = j + 1;
            end
        end
        detail.event_on = on;  detail.event_off = off;
        detail.events_per_min = numel(on) / (size(r, 2) * opts.dt) * 60;
    end
end
end
