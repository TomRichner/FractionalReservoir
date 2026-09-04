function name = full_adaptation_condition(conditions)
% FULL_ADAPTATION_CONDITION The name of a preset's most-adapted regime.
%
%   name = FULL_ADAPTATION_CONDITION(conditions)
%   model = build_from_preset(p, full_adaptation_condition(cond_cell));
%
% Answers "which of this preset's conditions has everything switched on", for
% the figures and tables that want ONE representative fully-adapted network.
%
% WHY A LOOKUP RATHER THAN A LITERAL. Those callers used to hardcode
% 'sfa_and_std'. That worked only while every preset used the same name for the
% role, and it stopped working the moment one did not: the 3-condition preset
% called its full regime sfa3_std2, and four figures plus the tau sweep failed
% with NoSuchCondition. Since condition names now state their adaptation
% STRUCTURE, the role is spread over more names still, and which one a preset
% uses depends on its own timescales:
%
%   sfa3_std2        3cond, 4cond, 7cond, single_neuron_dualStd, mc_pairs_dualStd
%   sfa3_std1        default, overconnected, bursting_pairs
%   sfa1_std1_stf1   single_neuron_stf
%   sfa3_std0        sompolinsky_pairs (no routes at all)
%
% Four names for one role, so a literal cannot be right everywhere and this
% lookup is what keeps the callers correct for presets that do not exist yet.
%
% HOW. Parses each name as sfa<X>_std<Y> with an optional _stf<Z> and returns the
% one with the greatest X+Y+Z; ties break on the larger X. Names that do not
% match -- 'no_adaptation', or anything a future preset invents -- are ignored,
% which is what keeps the paper's naming convention from becoming a requirement.
% test_preset_conditions checks that a name matching the pattern tells the truth
% about its counts, so the arithmetic here is on values that have been verified.
%
% Errors, naming every condition, if nothing matches: a preset with no
% conforming name has no answer to give, and guessing (the last condition, say)
% would silently hand back whatever happened to be listed last.
%
% See also: srnn_param_preset, build_from_preset, srnn_condition_titles

arguments
    conditions (1,:) cell
end

best = -1; best_sfa = -1; name = '';
names = cell(1, numel(conditions));
for i = 1:numel(conditions)
    names{i} = conditions{i}.name;
    counts = parse_condition_name(names{i});
    if isempty(counts); continue; end
    total = sum(counts);
    if total > best || (total == best && counts(1) > best_sfa)
        best = total; best_sfa = counts(1); name = names{i};
    end
end

if isempty(name)
    error('full_adaptation_condition:NoneMatch', ...
        ['No condition name matches sfa<X>_std<Y>[_stf<Z>], so none can be ' ...
         'identified as the fully-adapted regime.\n  Found: %s\n\n' ...
         'Pass the condition name explicitly if this preset uses a different ' ...
         'naming convention.'], strjoin(names, ', '));
end
end

%% ------------------------------------------------------------------------
function counts = parse_condition_name(name)
% [n_sfa n_std n_stf] for a conforming name, [] otherwise.
%
% TWO PATTERNS, EACH FULLY MANDATORY, rather than one with an optional
% (?:_stf(\d+))? group. MATLAB's regexp does NOT return a token for a trailing
% optional group even when it matches: 'sfa1_std1_stf1' against that pattern
% yields two tokens, not three, and the facilitation count silently reads as
% zero. Named tokens behave the same way (the field is simply absent). Verified
% both. Two exact patterns cannot express that failure.
counts = [];
t3 = regexp(name, '^sfa(\d+)_std(\d+)_stf(\d+)$', 'tokens', 'once');
if ~isempty(t3)
    counts = cellfun(@str2double, t3);
    return
end
t2 = regexp(name, '^sfa(\d+)_std(\d+)$', 'tokens', 'once');
if ~isempty(t2)
    counts = [cellfun(@str2double, t2), 0];
end
end
