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
% See also: srnn_param_preset, validate_preset_conditions, manuscript_style,
%           ParamSpaceAnalysis2/plot_sensitivity

% TWO TITLE VOCABULARIES, on purpose.
%
% The 4- and 7-regime sets name MECHANISMS -- "SFA Only", "STD Only",
% "SFA + STD" -- because those sets exist to separate SFA from STD.
%
% The 'single_multi' set names TIMESCALE STRUCTURE -- "Single-Timescale
% Adaptation" against "Multiple-Timescale Adaptation" -- because there the
% mechanisms are always present together and the only thing varying is how many
% timescales carry them. That is the comparison the paper is built on.
%
% This is why sfa3_std2 exists as a separate name for physics identical to
% sfa_and_std: one regime, two labels, so each set can say what its own
% comparison is about. sfa_and_std keeps "SFA + STD" for the sets that still
% use it.
%
% THE SINGLE_MULTI TITLES ARE ~3x LONGER than the mechanism ones (28 characters
% against 9). ParamSpaceAnalysis2's plotters set them as panel titles at
% FontSize 14 in figures sized 300 px per condition, where they fit but not by
% much. TR is aware and will handle any clipping in post; do NOT shorten these
% or adjust a figure to accommodate them.
titles = containers.Map( ...
    {'no_adaptation', 'sfa_only_oneTS', 'sfa_only', 'std_only_oneTS', ...
     'std_only', 'sfa3_std1', 'sfa_and_std', ...
     'sfa1_std1', 'sfa3_std2'}, ...
    {'No Adaptation', 'SFA (1 \tau)', 'SFA Only', 'STD (1 \tau)', ...
     'STD Only', 'SFA + STD (1 \tau)', 'SFA + STD', ...
     'Single-Timescale Adaptation', 'Multiple-Timescale Adaptation'});
end
