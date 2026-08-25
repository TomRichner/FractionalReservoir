function tau = srnn_sfa_timescales(K)
% SRNN_SFA_TIMESCALES The standard SFA timescale ladder, K values from 0.25 s to 10 s.
%
%   tau = SRNN_SFA_TIMESCALES(3)   % [0.25  1.5811  10]
%   tau = SRNN_SFA_TIMESCALES(1)   % 0.25   -- the FAST end
%   tau = SRNN_SFA_TIMESCALES(0)   % 1x0    -- no adaptation
%
% One place that knows the paper's adaptation timescales, so a condition asking
% for "the usual three" and a preset asking for "just the fast one" cannot drift
% apart.
%
% WHY THIS EXISTS RATHER THAN AN INLINE logspace. The model classes used to
% auto-fill tau_a from a COUNT with
%
%   logspace(log10(0.25), log10(10), n_a)
%
% and MATLAB's logspace(a,b,1) returns 10^b -- so asking for ONE timescale
% silently produced the SLOW 10 s end, not the fast 0.25 s one. That is wrong
% for the single-timescale conditions, and it is wrong in a way nothing reports:
% the model builds happily and runs a different experiment than intended.
%
% The deeper problem it exposed is that a COUNT CANNOT SAY WHICH TIMESCALES. At
% K = 1, 0.25, the geometric middle 1.58 and 10 are all defensible; no integer
% expresses the choice. So callers state the timescales, and this helper only
% supplies the standard ladder when that is what is wanted.
%
% K = 1 is special-cased to the fast end DELIBERATELY (TR's decision): the
% single-timescale conditions exist to ask whether fast adaptation alone
% reproduces what the spread does, so the fast end is the interesting one.
%
% See also: srnn_adaptation_conditions, srnn_param_preset

arguments
    K (1,1) double {mustBeInteger, mustBeNonnegative}
end

switch K
    case 0
        tau = zeros(1, 0);          % no adaptation; a 1x0 row, not []
    case 1
        tau = 0.25;                 % see above -- NOT logspace(...,1), which is 10
    otherwise
        tau = logspace(log10(0.25), log10(10), K);
end
end
