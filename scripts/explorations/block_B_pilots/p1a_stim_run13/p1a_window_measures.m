function [w, det, events] = p1a_window_measures(s, window, threshold)
%P1A_WINDOW_MEASURES Detector v1 and the state means on one time window (pilot 1a).
%
%   [w, det, events] = p1a_window_measures(s, window)
%   [w, det, events] = p1a_window_measures(s, window, threshold)
%
% s       p1a_signals output
% window  [t0 t1] s; samples t0 <= t < t1, the selection ied_stim_experiment
%         uses
% threshold  optional n x 1 absolute threshold for ied_detect_events (the
%         fixed-criterion control); default: the window's own adaptive one
%
% w is a struct of scalars: event_count, events_per_min, median_recruitment
% (detector, fraction of all neurons), median_amplitude (WaveformAmplitude of
% ied_event_details, a.u.), median_fwhm (s), rate_E, rate_I, rate (all
% neurons, the detector's mean_rate), sfa_E, sfa_I, sfa_doc, dep,
% within_corr, between_corr. det is ied_detect_events' output (without the
% n x nt 'active' array), events ied_event_details' table.

if nargin < 3, threshold = []; end
sel = s.t >= window(1) & s.t < window(2);
t = s.t(sel);
det = ied_detect_events(t, s.x(:, sel), s.r(:, sel), s.groups, struct(), threshold);
[events, spatial] = ied_event_details(struct('t', t, 'x', s.x(:, sel), 'groups', s.groups), det);
w = struct( ...
    'event_count', det.metrics.event_count, ...
    'events_per_min', det.metrics.events_per_min, ...
    'median_recruitment', det.metrics.median_recruitment, ...
    'median_amplitude', median(events.WaveformAmplitude, 'omitnan'), ...
    'median_fwhm', median(events.WaveformWidth, 'omitnan'), ...
    'rate_E', mean(s.rate_E(sel)), 'rate_I', mean(s.rate_I(sel)), ...
    'rate', det.metrics.mean_rate, ...
    'sfa_E', mean(s.sfa_E(sel)), 'sfa_I', mean(s.sfa_I(sel)), 'sfa_doc', mean(s.sfa_doc(sel)), ...
    'dep', mean(s.dep(sel)), ...
    'within_corr', spatial.within_corr, 'between_corr', spatial.between_corr);
if isempty(events) || height(events) == 0
    w.median_amplitude = NaN; w.median_fwhm = NaN;
end
det = rmfield(det, 'active');
end
