function titles = srnn_condition_titles()
% SRNN_CONDITION_TITLES Display names for the adaptation regimes, in one place.
%
%   titles = SRNN_CONDITION_TITLES()   % containers.Map, name -> display string
%
% Every plot that labels a condition column reads this, so a regime is titled
% identically wherever it appears and adding one is a single edit.
%
% It replaced FIVE hardcoded copies of the same four-entry map -- four inside
% ParamSpaceAnalysis2's plotters and one in manuscript_style -- which meant
% adding a regime required finding all five, and any that were missed fell back
% to the raw snake_case name in one figure while reading properly in another.
%
% It lives in src/ rather than beside manuscript_style because
% ParamSpaceAnalysis2 uses it: src must not depend on scripts/.
%
% Callers should still guard with isKey() where a run may carry conditions this
% does not know about -- a saved run directory can name anything.
%
% NAMES STATE THE ADAPTATION STRUCTURE (2026-09-03). A condition is
% 'no_adaptation' or 'sfaX_stdY', X and Y being the number of SFA and depression
% timescales actually in use, plus '_stfZ' where facilitation is on. Before that
% the names were mechanism-based -- sfa_only, std_only, sfa_and_std -- and
% 'sfa_and_std' meant FOUR different things across the ten presets: 3 SFA + 2
% STD in the paper's networks, 3 + 1 in bursting_pairs and both SRNNModel2
% presets, 3 + 0 in sompolinsky_pairs, and 1 + 1 + STF in single_neuron_stf. Two
% run directories could carry the same condition folder name and different
% physics, and nothing caught it.
%
% Titles are MOSTLY structural -- 'SFA 3\tau + STD 1\tau' and so on -- with two
% deliberate exceptions. sfa1_std1 and sfa3_std2 are the 3-condition preset's
% adapting regimes, i.e. the paper's actual comparison, and are titled
% "Single-Timescale Adaptation" and "Multiple-Timescale Adaptation" because that
% is the claim being made; the counts are how it is implemented, not what it
% means.
%
% That those two words survive a rename is the point of keeping this map keyed by
% NAME rather than putting a title inside each condition: the wording changed
% here, in one file, without touching a preset, a condition, or a saved run. A
% title frozen into a run directory at compute time could not have moved at all.
%
% The cost, accepted: sfa3_std2 also appears in the 7-condition preset, next to
% sfa3_std1, so that sheet mixes these words with a structural label. The
% 7-condition set is exploratory and not in the paper.
%
% See also: srnn_param_preset, validate_preset_conditions, manuscript_style,
%           ParamSpaceAnalysis2/plot_sensitivity

current = { ...
    'no_adaptation',    'No Adaptation'; ...
    'sfa1_std0',        'SFA 1\tau'; ...
    'sfa3_std0',        'SFA 3\tau'; ...
    'sfa0_std1',        'STD 1\tau'; ...
    'sfa0_std2',        'STD 2\tau'; ...
    ... % THE PAPER'S COMPARISON, and titled in its own words rather than
    ... % structurally: these two are the 3-condition preset's adapting regimes,
    ... % and the manuscript's point is single vs multiple timescales, not the
    ... % counts themselves. TR's decision.
    ... %
    ... % They also appear in the 7-condition preset, where sfa3_std2 sits beside
    ... % sfa3_std1 -- so that sheet mixes these words with a structural label.
    ... % Accepted: the 7-condition set is exploratory and not in the paper.
    'sfa1_std1',        'Single-Timescale Adaptation'; ...
    'sfa3_std1',        'SFA 3\tau + STD 1\tau'; ...
    'sfa3_std2',        'Multiple-Timescale Adaptation'; ...
    ... % single_neuron_stf only -- the one preset with facilitation. Without the
    ... % suffix its regimes would be sfa0_std1 and sfa1_std1, and sfa1_std1 is
    ... % already the 3-condition preset's SFA+STD with NO facilitation.
    'sfa0_std1_stf1',   'STD 1\tau + STF'; ...
    'sfa1_std1_stf1',   'SFA 1\tau + STD 1\tau + STF'};

% LEGACY KEYS, for run directories saved before the rename. Their condition
% subfolders are still named this way, so a figure regenerated from an old run
% would otherwise fall back to raw snake_case in its panel titles. Kept
% deliberately and indefinitely: they cost one row each and old runs are the
% only record of experiments that are expensive to repeat.
%
% Note these names are AMBIGUOUS -- that is why they were replaced -- so the
% titles below are the ones they carried at the time, not a claim about how many
% timescales any particular old run actually used.
legacy = { ...
    'sfa_only',         'SFA Only'; ...
    'std_only',         'STD Only'; ...
    'sfa_and_std',      'SFA + STD'; ...
    'sfa_only_oneTS',   'SFA (1 \tau)'; ...
    'std_only_oneTS',   'STD (1 \tau)'};

all_rows = [current; legacy];
titles = containers.Map(all_rows(:, 1), all_rows(:, 2));
end
