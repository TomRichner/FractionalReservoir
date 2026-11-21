function E = ospec_excursions(t, x, threshold)
%OSPEC_EXCURSIONS The excursions of a series beyond a threshold: start, end, duration, and whether the window cut them off.
%
%   E = ospec_excursions(t, x)              beyond zero
%   E = ospec_excursions(t, x, threshold)
%
%   t   1 x n_s sample times (s), evenly spaced
%   x   1 x n_s series (for the stage: a rightmost real part in 1/s,
%       R.rate_top or real(R.full.net)). NaN counts as NOT beyond: in
%       R.full.net a NaN is a state whose rightmost network mode is to the
%       left of every eigenvalue computed (ospec_full_spectrum), which is
%       below zero whenever R.full.bound is.
%
% An excursion is a run of consecutive samples with x > threshold.
%   E.t_start, E.t_end   1 x m: times of its first and last sample
%   E.duration   1 x m: (number of samples) x dt, as the exploration
%                counted (outlier_measures). The resolution is dt: a
%                shorter excursion can be missed and two separated by less
%                are merged.
%   E.cut_start, E.cut_end   1 x m logical: the excursion contains the
%                first, the last sample of the series, so it began before
%                the window or was still going on when it ended, and its
%                duration is a LOWER BOUND
%   E.censored   cut_start | cut_end
%   E.dt, E.window_s (n_s x dt), E.fraction (of samples beyond), E.threshold
%
% See also: ospec_measures, ospec_load_run

if nargin < 3 || isempty(threshold); threshold = 0; end
t = t(:)'; x = x(:)';
n_s = numel(t);
dt = NaN;
if n_s > 1; dt = (t(end) - t(1)) / (n_s - 1); end
mask = x > threshold;                    % NaN > threshold is false
d = diff([false mask false]);
on = find(d == 1); off = find(d == -1) - 1;
E = struct('t_start', t(on), 't_end', t(off), 'duration', (off - on + 1) .* dt, ...
    'cut_start', on == 1, 'cut_end', off == n_s, 'censored', (on == 1) | (off == n_s), ...
    'dt', dt, 'window_s', n_s * dt, 'fraction', mean(mask), 'threshold', threshold);
end
