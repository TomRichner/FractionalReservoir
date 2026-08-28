% REPRODUCE_PAPER_RUN Run the paper end to end: all compute, then all figures.
%
%   Open this file and press Run. Nothing else is needed -- not even a prior
%   setup_paths, which the first line does.
%
% Both entry points take their settings from paper_config() when called with no
% argument, so this script names nothing itself. **paper_config.m is the one file
% to edit** -- the preset (WHICH NETWORK) and the run mode (HOW MUCH COMPUTE)
% live there and nowhere else.
%
% WHAT IT PRODUCES
%   1. A dated run directory under data/param_space/run_all_<timestamp>/ holding
%      every sweep, the memory-capacity results, the eig-heatmap samples, a
%      run_manifest.mat, a git_provenance.txt and a parameters.md.
%   2. Every manuscript figure, written next to its own fig_*.m under
%      scripts/presentations/Stability_Manuscript/, plus the doc tables.
%
% COST at the default run mode ('medium'): about 3 h of compute and 8 min of
% figures on a 2026-era desktop (R2026a, Parallel Computing Toolbox, 16 workers).
% For the paper's final run, set run_mode to 'production' in paper_config --
% substantially longer, integrating with ode45 over a 50 s window and sampling
% 256 grid points instead of 64. 'fast' finishes in about a quarter hour and
% produces figures that are structurally correct but far too coarse to read.
%
% make_all_paper_figures resolves the newest run matching the preset, which is
% the one this script just produced. That is normally the same directory, but it
% is a lookup rather than a handoff: if another run of the same preset finishes
% in between, the figures follow that one instead. Pass the directory explicitly
% -- make_all_paper_figures(paper_config('run_dir', run_dir)) -- if that matters.
%
% ERROR ISOLATION. run_all_paper_analyses wraps each of its four stages, so one
% failure is reported and the queue continues rather than costing the whole run.
% make_all_paper_figures likewise reports per figure, and counts success as files
% actually on disk. Read the two summary tables at the end; do not assume that
% "it finished" means "everything worked".
%
% See also: paper_config, run_all_paper_analyses, make_all_paper_figures

setup_paths();

run_dir = run_all_paper_analyses();
results = make_all_paper_figures();

fprintf('\n========================================================\n');
fprintf('REPRODUCE PAPER RUN COMPLETE\n');
fprintf('  run_dir : %s\n', run_dir);
fprintf('  figures : %d of %d succeeded\n', sum([results.ok]), numel(results));
fprintf('========================================================\n');
