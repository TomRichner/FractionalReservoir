function names = full_adaptation_conditions(conditions)
% FULL_ADAPTATION_CONDITIONS Every condition at a preset's top adaptation count, chain first.
%
%   names = FULL_ADAPTATION_CONDITIONS(conditions)   % 1 x M cell of names, M >= 1
%
% The plural of full_adaptation_condition, which returns names{1}. It exists
% because a preset can now carry TWO fully-adapted regimes that differ only in
% the FORM of their depression, not in how much of it there is: the mu7chain
% bundle's sfa3_stdc2 (chained vesicle pools) and sfa3_std2 (the product of two
% independent pools). Both have three SFA timescales and two depression
% variables per route. A caller that wants one representative network takes
% the first; a caller that must cover both (the SFA-timescale sweep, Tom
% 2026-11-11) takes them all.
%
% RANKING. Each conforming name is parsed into [n_sfa n_std n_stf] and a chain
% flag. The top conditions are those with the greatest n_sfa + n_std + n_stf
% and, among those, the greatest n_sfa -- the rule full_adaptation_condition
% has always used. Within that top set the order is:
%   1. chained depression first (sfaX_stdcY), by an EXPLICIT rule rather than
%      by whichever the preset lists first: the chain is the paper's MTS, and a
%      four-condition preset lists the chain before the product only by
%      convention. Without this rule the strict ">" that picked the winner
%      returned the first-listed of a tie, so reordering a preset would have
%      silently swapped which network every figure drew.
%   2. otherwise the preset's own order.
%
% NAMES. Three exact patterns, each fully mandatory:
%   sfa<X>_std<Y>_stf<Z>   product or single-pool depression with facilitation
%   sfa<X>_std<Y>          product or single-pool depression
%   sfa<X>_stdc<Y>         chained depression; Y is the number of pools (2)
% Names matching none ('no_adaptation', or anything a future preset invents)
% are ignored. test_preset_conditions checks that a name matching a pattern
% tells the truth about its counts and, for stdc, that its routes really are
% chains.
%
% Errors, naming every condition, if nothing matches.
%
% See also: full_adaptation_condition, srnn_param_preset, build_from_preset

arguments
    conditions (1,:) cell
end

n = numel(conditions);
all_names = cell(1, n);
counts = nan(n, 3);
is_chain = false(1, n);
for i = 1:n
    all_names{i} = conditions{i}.name;
    [c, ch] = parse_condition_name(all_names{i});
    if isempty(c); continue; end
    counts(i, :) = c;
    is_chain(i) = ch;
end

ok = ~isnan(counts(:, 1))';
if ~any(ok)
    error('full_adaptation_condition:NoneMatch', ...
        ['No condition name matches sfa<X>_std<Y>[_stf<Z>] or sfa<X>_stdc<Y>, so ' ...
         'none can be identified as the fully-adapted regime.\n  Found: %s\n\n' ...
         'Pass the condition name explicitly if this preset uses a different ' ...
         'naming convention.'], strjoin(all_names, ', '));
end

total = sum(counts, 2)';
best_total = max(total(ok));
top = ok & total == best_total;
best_sfa = max(counts(top, 1));
top = top & counts(:, 1)' == best_sfa;

idx = find(top);
% Chain first, then preset order: a stable sort on "not a chain".
[~, order] = sort(~is_chain(idx));
names = all_names(idx(order));
end

%% ------------------------------------------------------------------------
function [counts, is_chain] = parse_condition_name(name)
% [n_sfa n_std n_stf] for a conforming name, [] otherwise; is_chain for stdc.
%
% EXACT PATTERNS, EACH FULLY MANDATORY, rather than one with an optional
% (?:_stf(\d+))? group or an optional 'c'. MATLAB's regexp does NOT return a
% token for a trailing optional group even when it matches: 'sfa1_std1_stf1'
% against that pattern yields two tokens, not three, and the facilitation count
% silently reads as zero. Named tokens behave the same way (the field is simply
% absent). Verified both. Exact patterns cannot express that failure.
counts = []; is_chain = false;
t3 = regexp(name, '^sfa(\d+)_std(\d+)_stf(\d+)$', 'tokens', 'once');
if ~isempty(t3)
    counts = cellfun(@str2double, t3);
    return
end
t2 = regexp(name, '^sfa(\d+)_std(\d+)$', 'tokens', 'once');
if ~isempty(t2)
    counts = [cellfun(@str2double, t2), 0];
    return
end
tc = regexp(name, '^sfa(\d+)_stdc(\d+)$', 'tokens', 'once');
if ~isempty(tc)
    counts = [cellfun(@str2double, tc), 0];
    is_chain = true;
end
end
