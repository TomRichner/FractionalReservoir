function v = p4_candidates(r, is_E, k)
% P4_CANDIDATES Every candidate measure of pilot 4 on one trace, as one struct of scalars.
%
%   v = p4_candidates(r, is_E, k)
%
% r      n x nt per-neuron rates over the analysis window
% is_E   n x 1 logical, the excitatory neurons (the population-mean
%        candidates use the mean E rate m(t), as the burst detector does;
%        the co-activation candidates use every neuron and ignore is_E)
% k      robust standard deviations above the median (default 3)
%
% rsd(.) below is 1.4826 * median absolute deviation.
%
%   pop_time       1.  fraction of time m(t) > median(m) + k rsd(m)
%   pop_mass       2.  time mean of max(m(t) - median(m) - k rsd(m), 0)
%   coact_time10   3a. fraction of time c(t) > 0.10, c(t) = fraction of
%                      neurons above median_i + max(k rsd_i, 0.1)
%   coact_q99      3b. 99th percentile of c(t)
%   coact_excess   3c. A of coactivation_excess (k, floor 0.1, z 5): the
%                      time mean of c(t) above the independent-neuron level
%   coact_time     3d. fraction of time c(t) is above that level
%   gated_mass     3e. time mean of [mean over neurons of max(r_i - thr_i, 0)]
%                      at the times c(t) is above that level
%   synchrony      4.  var_t(m) / mean over E neurons of var_t(r_i)
%   skewness       5a. skewness of m(t)
%   p99_over_median 5b. 99th percentile of m(t) / median(m)
% and the companions mean_rate, mean_rate_E, pinned_fraction, class.
%
% See also: coactivation_excess, run_p4_aberrant_measure

if nargin < 3, k = 3; end
r = double(r);
m = mean(r(is_E, :), 1);
mm = median(m);
ms = 1.4826 * median(abs(m - mm));
v = struct();
v.pop_time = mean(m > mm + k * ms);
v.pop_mass = mean(max(m - mm - k * ms, 0));

[A, comp, det] = coactivation_excess(r, 'k', k);
v.coact_time10 = mean(det.c > 0.10);
v.coact_q99 = det.q99;
v.coact_excess = A;
v.coact_time = det.time_fraction;
ex = mean(max(r - det.threshold, 0), 1);
v.gated_mass = mean(ex .* (det.c > det.level));

v.synchrony = var(m) / mean(var(r(is_E, :), 0, 2));
v.skewness = skewness(m);
v.p99_over_median = prctile(m, 99) / mm;

v.mean_rate = comp.mean_rate;
v.mean_rate_E = mean(m);
v.pinned_fraction = comp.pinned_fraction;
v.class = comp.class;
end
