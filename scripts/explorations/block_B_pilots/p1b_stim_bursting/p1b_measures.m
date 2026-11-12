function M = p1b_measures(C, S)
% P1B_MEASURES Pilot 1b's measures of one run, before and during the stimulus.
%
%   M = P1B_MEASURES(C, S)
%
% C  one case from p1b_run_case (t, rates, mean adaptation and depression,
%    bursts detected over the whole [0, T_show])
% S  the settings of p1b_stim_bursting: S.win_before, S.win_during (s),
%    S.t_stim (s), S.pinned_thresh
%
% M.before and M.during, each with
%   n_bursts         bursts (events no longer than S.burst.sustained_max)
%                    whose ONSET lies in the window
%   burst_rate       n_bursts / window length, per s
%   frac_above       fraction of the window's time with mean E rate above
%                    the burst threshold (separates many short bursts from
%                    one long stretch of high firing, which the burst count
%                    alone would call one burst)
%   longest_event    s, longest above-threshold event overlapping the window
%   n_sustained      events overlapping the window longer than
%                    S.burst.sustained_max (sustained high firing; these
%                    are NOT counted as bursts anywhere)
%   mean_r_E, mean_r_I, max_mean_r_E
%   sd_mean_r_E, p2p_mean_r_E   sd and max - min of the population mean
%                    E rate (the size of the fluctuation, threshold-free)
%   dominant_freq    Hz, largest spectral peak of the mean E rate, 0.05-10 Hz
%   mean_adaptation  mean (c/K) sum_k a_k over E
%   mean_depression  mean depression factor over the depressed routes
%   frac_pinned_t    time-mean fraction of all neurons with r > pinned_thresh
%   frac_pinned_n    fraction of neurons whose window-mean rate > pinned_thresh
%   mean_r_target    mean rate of the stimulated neurons (all of them when
%                    every neuron is stimulated)
% M.t_last_burst     s from stimulus onset to the ONSET of the last burst
%                    that starts at or after onset; NaN when none does
% M.n_bursts_after_onset  bursts starting in [t_stim, T_show]

wins = struct('before', S.win_before, 'during', S.win_during);
fn = fieldnames(wins);
rAll = [double(C.r_E); double(C.r_I)];
B = C.bursts;
% An above-threshold event longer than S.burst.sustained_max is SUSTAINED
% high firing, not a burst (a DC step that lifts the whole population
% above threshold makes one 60 s "event" under the bare definition).
is_burst = B.duration <= S.burst.sustained_max;
for w = 1:numel(fn)
    lim = wins.(fn{w});
    sel = C.t >= lim(1) & C.t <= lim(2);
    inb = is_burst & B.t_on >= lim(1) & B.t_on < lim(2);
    touch = B.t_on <= lim(2) & B.t_off >= lim(1);
    m = struct();
    m.window        = lim;
    m.n_bursts      = nnz(inb);
    m.burst_rate    = m.n_bursts / diff(lim);
    m.frac_above    = mean(C.mean_r_E(sel) > S.burst.thresh);
    m.longest_event = max([0, B.duration(touch)]);
    m.n_sustained   = nnz(touch & ~is_burst);
    m.mean_r_E      = mean(C.mean_r_E(sel));
    m.mean_r_I      = mean(C.mean_r_I(sel));
    m.max_mean_r_E  = max(C.mean_r_E(sel));
    % How much the population mean E rate moves, whatever the threshold
    % says: a steady state that sits AT the burst threshold crosses it on
    % every small wiggle, and the bare count then reports "bursts".
    y = C.mean_r_E(sel);
    m.sd_mean_r_E   = std(y);
    m.p2p_mean_r_E  = max(y) - min(y);
    dt = C.t(2) - C.t(1);
    Y  = abs(fft(y - mean(y)));
    fr = (0:numel(y) - 1) / (numel(y) * dt);
    band = fr >= 0.05 & fr <= 10;
    [~, i] = max(Y(band)); fb = fr(band);
    m.dominant_freq = fb(i);         % Hz, largest spectral peak in 0.05-10 Hz
    m.mean_adaptation = mean(C.mean_adaptation(sel));
    m.mean_depression = mean(C.mean_depression(sel));
    m.frac_pinned_t = mean(mean(rAll(:, sel) > S.pinned_thresh, 1));
    m.frac_pinned_n = mean(mean(rAll(:, sel), 2) > S.pinned_thresh);
    m.mean_r_target = mean(mean(rAll(C.target, sel), 2));
    M.(fn{w}) = m;
end
after = is_burst & B.t_on >= S.t_stim;
M.n_bursts_after_onset = nnz(after);
if any(after)
    M.t_last_burst = max(B.t_on(after)) - S.t_stim;
else
    M.t_last_burst = NaN;
end
end
