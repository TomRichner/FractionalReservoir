function names = run_mode_names()
% RUN_MODE_NAMES The run modes every stage of the pipeline must support.
%
%   names = RUN_MODE_NAMES()   % {'fast', 'medium', 'medium2', 'production'}
%
% ONE list, because there used to be four and they disagreed. analysis_run_config
% accepted five modes; run_memory_capacity, run_memory_capacity_example and
% run_eig_heatmap each had their own switch accepting three. A run_mode valid for
% the sweeps was therefore invalid for three later stages, and since
% run_all_paper_analyses wraps each stage, that surfaced as three stages silently
% producing nothing after the sweeps had already run.
%
% That is exactly what happened to single_multi_TS on 2026-09-03 with
% run_mode = 'fast2': the sweeps completed, memory capacity / mc_example /
% eig_heatmap all threw "Unknown run_mode", their two figures failed on missing
% data, and fig_memory_capacity fell back to a twelve-day-old .mat and reported
% SUCCESS. One stale figure looked exactly like a good one.
%
% 'fast2' HAS SINCE BEEN REMOVED (TR): it was 'fast' with more reps and a longer
% window on two of the three sweeps, which is a narrow enough difference that
% carrying a whole extra mode for it was not worth the coupling. Use 'fast'.
%
% Every consumer validates against THIS list, and test_run_modes asserts that
% each stage really accepts every name in it -- so a mode cannot be added to the
% sweeps again without the later stages being taught about it.
%
% See also: analysis_run_config, run_all_paper_analyses, test_run_modes

names = {'fast', 'medium', 'medium2', 'production'};
end
