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
%   sfa3_stdc2       the mu7chain bundle (chained vesicle-pool depression),
%                    which ALSO carries sfa3_std2 and must return the chain
%
% Five names for one role, so a literal cannot be right everywhere and this
% lookup is what keeps the callers correct for presets that do not exist yet.
%
% HOW. Parses each name as sfa<X>_std<Y>, sfa<X>_std<Y>_stf<Z> or
% sfa<X>_stdc<Y> (three exact patterns; MATLAB's regexp drops the token of a
% trailing optional group, so optional groups cannot be used) and returns the
% one with the greatest X+Y+Z; ties break on the larger X, and then a CHAIN
% (stdc) wins over a product by an explicit rule, not by its place in the
% preset. Names that do not match -- 'no_adaptation', or anything a future
% preset invents -- are ignored, which is what keeps the paper's naming
% convention from becoming a requirement. test_preset_conditions checks that a
% name matching a pattern tells the truth about its counts, so the arithmetic
% here is on values that have been verified.
%
% WHAT WENT WRONG BEFORE THE THIRD PATTERN (2026-11-12, FD-007). With only the
% two product patterns, sfa3_stdc2 matched neither and was skipped: given the
% paper's three conditions (no_adaptation, sfa1_std1, sfa3_stdc2) this
% returned sfa1_std1, the SINGLE-timescale regime, with no error.
%
% The ranking lives in full_adaptation_conditions (plural), which returns
% every condition tied at the top, chain first; this returns the first of
% them.
%
% Errors, naming every condition, if nothing matches: a preset with no
% conforming name has no answer to give, and guessing (the last condition, say)
% would silently hand back whatever happened to be listed last.
%
% See also: full_adaptation_conditions, srnn_param_preset, build_from_preset,
%           srnn_condition_titles

arguments
    conditions (1,:) cell
end

names = full_adaptation_conditions(conditions);
name = names{1};
end
