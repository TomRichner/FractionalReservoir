% REPRODUCE_PAPER_RUN Run the paper end to end: all compute, then all figures.
%
%   Open this file and press Run. Nothing else is needed -- not even a prior
%   setup_paths, which the first line does.
%
% This is the script form of the two-call sequence the project is built around:
%
%     cfg     = paper_config('preset_name', PRESET, 'run_mode', RUN_MODE);
%     run_dir = run_all_paper_analyses(cfg);            % hours
%     results = make_all_paper_figures(paper_config(... 'run_dir', run_dir));
%
% It exists so that repeating the paper's analysis is one action rather than a
% remembered pair of calls. The two entry points remain callable on their own;
% this adds no capability, only a documented starting point.
%
% WHAT IT PRODUCES
%   1. A dated run directory under data/param_space/run_all_<timestamp>/
%      holding every sweep, the memory-capacity results, the eig-heatmap
%      samples, a run_manifest.mat, a git_provenance.txt and a parameters.md.
%   2. Every manuscript figure, written next to its own fig_*.m under
%      scripts/presentations/Stability_Manuscript/, plus the doc tables.
%
% COST. At RUN_MODE = 'medium' the compute stage took 3.0 h and the figures
% 7.7 min on a 2026-era desktop (R2026a, Parallel Computing Toolbox, 16 workers).
% 'production' is substantially longer -- it integrates with ode45 over a 50 s
% window and samples 256 grid points instead of 64. Start with 'fast' if you
% only want to confirm the pipeline runs; it finishes in about a quarter hour and
% produces figures that are structurally correct but far too coarse to read.
%
% REPRODUCIBILITY. The run records the commit it ran at, so a run directory
% describes itself. One deliberate exception: the parameter-space stage draws its
% subset of grid points with rng('shuffle'), so rerunning picks a different
% sample of the joint grid. The 1-D sensitivity sweeps and the tau sweep are
% unaffected -- they cover every level in order, and each point's network seed is
% tied to its grid position rather than to execution order.
%
% ERROR ISOLATION. run_all_paper_analyses wraps each of its four stages, so one
% failure is reported and the queue continues rather than costing the whole run.
% make_all_paper_figures likewise reports per figure, and counts success as files
% actually on disk. Read the two summary tables at the end; do not assume that
% "it finished" means "everything worked".
%
% See also: paper_config, run_all_paper_analyses, make_all_paper_figures

setup_paths();

%% ------------------------------------------------------------------------
% The two knobs. They are orthogonal: the preset says WHICH NETWORK, the run
% mode says HOW MUCH COMPUTE.
%
% Both are named here explicitly rather than left to paper_config's defaults, so
% this script keeps meaning the same thing if those defaults are ever changed.
% To sweep a different network, prefer editing paper_config -- it is the single
% place the project names what the paper is built from, and both entry points
% follow it -- and change these only for a one-off.
PRESET   = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond';
RUN_MODE = 'medium';    % 'fast' | 'fast2' | 'medium' | 'medium2' | 'production'
%% ------------------------------------------------------------------------

fprintf('\n========================================================\n');
fprintf('REPRODUCE PAPER RUN\n');
fprintf('  preset   : %s\n', PRESET);
fprintf('  run_mode : %s\n', RUN_MODE);
fprintf('  start    : %s\n', datetime('now'));
fprintf('========================================================\n');
t_all = tic;

%% 1. All heavy compute
cfg     = paper_config('preset_name', PRESET, 'run_mode', RUN_MODE);
run_dir = run_all_paper_analyses(cfg);

%% 2. All figures, from the run just produced
% run_dir is passed EXPLICITLY rather than letting make_all_paper_figures resolve
% the newest matching run. They are normally the same directory, but only the
% explicit form guarantees the figures came from the sweep this script just ran
% -- an unrelated run finishing in between would otherwise be picked up silently.
fig_cfg = paper_config('preset_name', PRESET, 'run_mode', RUN_MODE, ...
                       'run_dir', run_dir);
results = make_all_paper_figures(fig_cfg);

%% Summary
n_ok = sum([results.ok]);
fprintf('\n========================================================\n');
fprintf('REPRODUCE PAPER RUN COMPLETE in %.2f h\n', toc(t_all)/3600);
fprintf('  run_dir : %s\n', run_dir);
fprintf('  figures : %d of %d succeeded\n', n_ok, numel(results));
fprintf('========================================================\n');
if n_ok < numel(results)
    fprintf(2, '  Some figures failed -- see the per-figure table above.\n');
end
