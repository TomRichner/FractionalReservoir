function K = censored_survival(duration, censored)
% CENSORED_SURVIVAL Kaplan-Meier survival of durations, some of which are only lower bounds.
%
%   K = censored_survival(duration, censored)
%
%   duration   lengths (s), any shape; for the outlier figures the
%              excursions of a rightmost eigenvalue beyond zero
%              (ospec_excursions' E.duration)
%   censored   logical, same size: true where the length is a LOWER BOUND
%              (the window cut the excursion off at its start or its end,
%              E.censored). Empty or omitted: none.
%
% WHY. An excursion that was already under way when the window began, or
% still going when it ended, lasted AT LEAST its observed length. Counting
% it as complete makes long excursions look shorter than they are (a
% condition beyond zero for the whole window would read as "every
% excursion lasts exactly one window"); dropping it throws away exactly
% the long ones. The product-limit estimate uses it for what it is known
% to say: it was still going at its observed length. S drops only at the
% lengths of COMPLETE excursions, by the fraction of those still under
% observation that ended there:
%
%   S(t) = prod over complete lengths t_i <= t of (1 - d_i / n_i)
%
% d_i the number that ended at t_i, n_i the number, complete or cut, with
% length >= t_i (a cut one of exactly length t_i counts as still going).
% With nothing cut this is the ordinary empirical survival. An excursion
% cut at its START is treated like one cut at its end (its length is a
% lower bound either way). Beyond the longest observed length the estimate
% says nothing when that longest one was cut: K.open is then true, S stays
% above zero, and K.at returns NaN there.
%
% K, a struct:
%   t, S        the distinct complete lengths, ascending, and the survival
%               just AFTER each (column vectors; empty when none completed)
%   n, n_censored   number of excursions, and of those cut
%   t_max       the longest observed length (NaN when there are none)
%   open        the longest observed length was cut: S(end) > 0 and
%               nothing is known beyond t_max
%   at          @(tq) S at the times tq, the estimated fraction of
%               excursions LONGER than tq: 1 before the first complete
%               length, NaN beyond t_max when open (0 when not), NaN
%               everywhere when there are no excursions
%   quantile    @(p) the smallest complete length t with S(t) <= 1 - p
%               (p = 0.5: the median dwell time); NaN when the estimate
%               never gets that low (not reached within the window)
%
% See also: ospec_excursions, fig_outlier_dwell

if nargin < 2 || isempty(censored); censored = false(size(duration)); end
duration = duration(:); censored = logical(censored(:));
if numel(censored) ~= numel(duration)
    error('censored_survival:Size', 'duration and censored must have the same number of elements.');
end
if any(~isfinite(duration)) || any(duration < 0)
    error('censored_survival:Duration', 'Durations must be finite and not negative.');
end
K = struct('t', zeros(0, 1), 'S', zeros(0, 1), 'n', numel(duration), 'n_censored', nnz(censored), ...
    't_max', NaN, 'open', false);
if ~isempty(duration)
    K.t_max = max(duration);
    t = unique(duration(~censored));
    S = ones(numel(t), 1);
    s = 1;
    for i = 1:numel(t)
        n_risk = nnz(duration >= t(i));
        d = nnz(~censored & duration == t(i));
        s = s * (1 - d / n_risk);
        S(i) = s;
    end
    K.t = t; K.S = S;
    % A cut excursion as long as the longest complete one keeps S above zero.
    K.open = isempty(S) || S(end) > 0;
end
K.at = @(tq) survival_at(K, tq);
K.quantile = @(p) survival_quantile(K, p);
end

function S = survival_at(K, tq)
S = nan(size(tq));
if K.n == 0; return; end
S(:) = 1;
for i = 1:numel(K.t)
    S(tq >= K.t(i)) = K.S(i);
end
if K.open
    S(tq > K.t_max) = NaN;
end
end

function q = survival_quantile(K, p)
q = NaN;
i = find(K.S <= 1 - p + 1e-12, 1);
if ~isempty(i); q = K.t(i); end
end
