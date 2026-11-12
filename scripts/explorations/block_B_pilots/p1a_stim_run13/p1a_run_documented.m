function rows = p1a_run_documented(seed, A)
%P1A_RUN_DOCUMENTED The documented stimulation protocol of 2026-10-05, re-run for one seed (pilot 1a's regression check).
%
%   rows = p1a_run_documented(seed, A)
%
% ied_stim_experiment's protocol on the 'original' condition, rebuilt through
% run13_condition_model: one run of 10 s settling + 120 s, uniform input 0
% for 0..60 s and A for 60..120 s (dc_staircase_stimulus, no ramp, no input
% noise), detector v1 applied separately to 10..60 s (no stimulation) and
% 70..120 s (stimulation). Returns two rows (no_stim, stim) with the
% quantities of docs/IED_model/stim_2026_10_05_summary.csv that can be
% compared: event_count, events_per_min, events_per_min_fixed,
% median_recruitment, mean_rate, sfa_feedback (sum_q c_q sum a / n) and
% resources (mean E-output depression factor).

% Pin the generator type to the client's (see p1a_run_pair: the model seeds
% the noise and the setpoints with rng(seed), which keeps the session's type).
prior_rng = rng; rng(0, 'twister'); restore_rng = onCleanup(@() rng(prior_rng)); %#ok<NASGU>
ic = struct('dc_levels', [0 A], 'hold_dur', 60, 'ramp_dur', 0, ...
    'noise_intensity', 0, 'intrinsic_drive', [], 'generator', @dc_staircase_stimulus);
[m, groups, info] = run13_condition_model('original', 'duration', 120, 'seed', seed, 'input_config', ic);
m.run();
s = p1a_signals(m, info, groups);
[w0, det0] = p1a_window_measures(s, [10 60]);
w1 = p1a_window_measures(s, [70 120]);
% Fixed criterion: the no-stimulation excursion above each neuron's median,
% re-centred on the stimulated window's own median.
wf = p1a_window_measures(s, [70 120], (det0.threshold - det0.baseline) + p1a_baseline(s, [70 120]));
rows = table(["no_stim"; "stim"], [seed; seed], [A; A], ...
    [w0.event_count; w1.event_count], [w0.events_per_min; w1.events_per_min], ...
    [w0.events_per_min; wf.events_per_min], ...
    [w0.median_recruitment; w1.median_recruitment], [w0.rate; w1.rate], ...
    [w0.sfa_doc; w1.sfa_doc], [w0.dep; w1.dep], ...
    'VariableNames', {'segment', 'seed', 'amplitude', 'event_count', 'events_per_min', ...
    'events_per_min_fixed', 'median_recruitment', 'mean_rate', 'sfa_feedback', 'resources'});
end

function b = p1a_baseline(s, window)
% The detector's per-neuron baseline on a window: median of the 25-ms moving mean.
sel = s.t >= window(1) & s.t < window(2);
fs = 1 / median(diff(s.t(sel)));
b = median(movmean(double(s.x(:, sel)), max(1, round(.025 * fs)), 2), 2);
end
