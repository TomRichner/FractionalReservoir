function values = log_ladder(lo, hi, n)
%LOG_LADDER n logarithmically spaced values from lo to hi, sane at n = 0 and 1.
%
%   values = LOG_LADDER(0.25, 10, 3)   % [0.25  1.5811  10]
%   values = LOG_LADDER(0.25, 10, 1)   % 0.25   -- the FIRST endpoint
%   values = LOG_LADDER(0.25, 10, 0)   % 1x0    -- empty, not []
%
% ENDPOINTS, NOT EXPONENTS. logspace takes log10 of its endpoints, so every
% caller in this repo wrapped it: logspace(log10(0.25), log10(10), 3). This takes
% the values you mean.
%
% WHY IT EXISTS. MATLAB's logspace(a, b, 1) returns 10^b -- the LAST endpoint,
% not the first. Asking for one value out of a range silently gives you the top
% of it, and nothing reports the substitution: the caller builds happily and runs
% a different experiment than intended. That is not hypothetical here. The SFA
% auto-fill did exactly this, so a single-timescale adaptation condition got the
% SLOW 10 s constant when it wanted the fast 0.25 s one, and the model ran
% without complaint.
%
% n = 1 RETURNS lo, THE FIRST ARGUMENT -- defined positionally, not as "the
% smaller one". So log_ladder(10, 0.25, 1) is 10. A descending range is
% something logspace already supports, and a rule about argument ORDER stays
% predictable where a rule about MAGNITUDE would silently flip.
%
% Picking the first endpoint rather than the last is the substantive choice, and
% it is the one a ladder wants: n = 1 is the degenerate case of "start here and
% spread upward", so it should give the start. For the SFA timescales that is
% also the scientifically interesting end -- the single-timescale conditions ask
% whether FAST adaptation alone reproduces what a spread does.
%
% n = 0 RETURNS zeros(1, 0), not []. Both are empty, but a 1x0 row concatenates
% and reshapes predictably where 0x0 does not, and SRNNCellTypePairs stores
% per-type timescales in a cell of rows.
%
% THIS REPLACED a helper named srnn_sfa_timescales(K), which hardcoded 0.25 and
% 10 in its own body and so could
% only ever speak for one ladder -- three different ladders already existed in
% the repo (0.25-10, 0.1-10, 0.5-15) and the other two got no protection from the
% n = 1 trap. It also wrote 0.25 twice, once as a literal for the n = 1 case and
% once inside log10() for the rest, so retuning the fast end would have missed
% one. Here the n = 1 case IS lo, derived rather than restated.
%
% See also: logspace, srnn_param_preset, srnn_adaptation_conditions

arguments
    lo (1,1) double {mustBePositive, mustBeFinite}
    hi (1,1) double {mustBePositive, mustBeFinite}
    n  (1,1) double {mustBeInteger, mustBeNonnegative}
end

switch n
    case 0
        values = zeros(1, 0);
    case 1
        values = lo;              % NOT logspace(...,1), which returns hi
    otherwise
        values = logspace(log10(lo), log10(hi), n);
end
end
