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
% See also: srnn_adaptation_conditions, manuscript_style,
%           ParamSpaceAnalysis2/plot_sensitivity

% sfa1_std1 belongs to the 'single_multi' set: SFA and STD both present, one
% timescale each. Read against its neighbours, "1 \tau each" is what
% distinguishes it from sfa3_std1, which has one STD timescale but the full SFA
% ladder.
%
% NOTE the 'single_multi' set's story is "none vs one timescale vs many", so a
% figure using it may want its legend to read "Single Timescale" / "Multi
% Timescale" instead. That would mean retitling sfa_and_std, which the 4- and
% 7-regime sets share -- so it is a deliberate decision, not a local tweak, and
% the accurate-and-consistent titles are what live here until it is made.
titles = containers.Map( ...
    {'no_adaptation', 'sfa_only_oneTS', 'sfa_only', 'std_only_oneTS', ...
     'std_only', 'sfa1_std1', 'sfa3_std1', 'sfa_and_std'}, ...
    {'No Adaptation', 'SFA (1 \tau)', 'SFA Only', 'STD (1 \tau)', ...
     'STD Only', 'SFA + STD (1 \tau each)', 'SFA + STD (1 \tau)', 'SFA + STD'});
end
