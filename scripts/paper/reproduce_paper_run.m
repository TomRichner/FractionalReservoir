% REPRODUCE_PAPER_RUN Rerun the manuscript's run end to end: all compute, then all figures.
%
%   Open this file and press Run. Nothing else is needed -- not even a prior
%   setup_paths, which the first line does.
%
% WHICH RUN. Manuscript5.md is built against figs/mu7revisedMedium (its
% figure_sources.json links figs/ there), and that figure set is the output of
% the mu7revisedMedium bundle. This script runs that bundle: every setting
% comes from mu7revisedMedium_config(), which states the preset, the run mode
% ('medium'), the side presets and the figure registry itself. When the
% manuscript moves to another bundle, change the one config name below.
%
% Until 2026-11-11 this script called both entry points with no argument, so
% it ran paper_config(), which names the older mu8p25 network: it reproduced a
% superseded run, not the manuscript's. paper_config is still what the two
% entry points use when called bare; it is not what the paper was built from.
%
% WHAT IT PRODUCES
%   1. data/mu7revisedMedium: every sweep, memory capacity, the eig-heatmap
%      samples, the numerics, spectrum, transient-gain, local-Lyapunov and
%      illustration stages, with run_manifest.mat, git_provenance.txt,
%      parameters.md and command_window.log.
%   2. figs/mu7revisedMedium: one folder per registry entry (main Figures 1 to
%      6 are the fig_mainN_grouped entries), with manifest.md, report.md and
%      git_provenance.txt at the root. OVERWRITTEN IN PLACE, and this is the
%      folder the manuscript's figs/ link reads.
%
% WHAT IT DOES NOT PRODUCE. Figure 7 is the clinical figure (Fig_Tica_Main.png,
% from the cscs_stim_stability pipeline through the Stim-Analysis-Figs repo).
% Figure 8 is the IED-like network: ied_stim_experiment then fig_ied_stim in
% scripts/explorations/IED_model, written to figs/IED_model/stim_20261005 (the
% record is docs/IED_model/).
%
% IT REFUSES A RUN DIRECTORY THAT ALREADY HOLDS A RUN. The manuscript's own
% run is in data/mu7revisedMedium (part of it tracked in git), so on a machine
% that has it run_all_paper_analyses stops with
% run_all_paper_analyses:TargetNotEmpty before computing or writing anything.
% That is deliberate: a reused directory accumulates sweep folders and the
% figures would pick between them. Move or delete data/mu7revisedMedium first
% if you mean to recompute. To redraw the figures from the run that is
% already there, without recomputing:
%
%   setup_paths(); results = make_all_paper_figures(mu7revisedMedium_config());
%
% COST. The analyses took 7.05 h and the figures 19 min on the 2026-09-16 run
% (R2026a, Parallel Computing Toolbox, n = 500 network at 'medium').
%
% WHAT "REPRODUCES" MEANS. One part of the run is not seeded: the joint
% parameter-space sample. ParamSpaceAnalysis2 draws it with rng('shuffle'),
% so its 128 grid points differ on every run, and whatever is drawn from the
% param_space stage matches the manuscript's statistically, not point for
% point.
%
% ERROR ISOLATION. run_all_paper_analyses wraps each stage, so one failure is
% reported and the rest continue. make_all_paper_figures likewise reports per
% figure, and counts success as files actually on disk. Read the two summary
% tables at the end; do not assume that "it finished" means "everything
% worked".
%
% See also: mu7revisedMedium_config, run_all_paper_analyses,
%           make_all_paper_figures, grouped_figure_registry, paper_config

setup_paths();

cfg     = mu7revisedMedium_config();      % the bundle the manuscript is built against
run_dir = run_all_paper_analyses(cfg);    % all compute -> data/mu7revisedMedium
results = make_all_paper_figures(cfg);    % all figures -> figs/mu7revisedMedium

fprintf('\n========================================================\n');
fprintf('REPRODUCE PAPER RUN COMPLETE\n');
fprintf('  preset  : %s (%s)\n', cfg.preset_name, cfg.run_mode);
fprintf('  run_dir : %s\n', run_dir);
fprintf('  figures : %d of %d succeeded\n', sum([results.ok]), numel(results));
fprintf('========================================================\n');
