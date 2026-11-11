% SINGLE_MULTI_TS_RUN Run the single- vs multi-timescale comparison end to end.
%
%   Open this file and press Run. Nothing else is needed -- not even a prior
%   setup_paths, which the first line does.
%
% Every setting comes from single_multi_TS_config(). **That is the file to
% edit**, and the one line most likely to need editing is run_mode: it is
% 'fast' for now and wants 'medium' or 'production' once the comparison is
% worth real compute.
%
%   preset    celltype_pairs_Sc0p2_noise0p025_dualStd_3cond
%   run_mode  fast
%   run_dir   data/single_multi_TS
%   fig_root  figs/single_multi_TS
%
% WHAT IT COMPARES: no adaptation vs adaptation on ONE timescale vs adaptation
% on MANY, SFA and STD always together.
%
%   no_adaptation   -                      -
%   sfa1_std1       1 tau_a (0.25 s)       1 STD pair (2 / 0.25 s)
%   sfa3_std2       3 tau_a (0.25-10 s)    2 STD pairs
%
% The contrast is TIMESCALE COUNT and only that, which is clean only because c
% is the total SFA budget divided by the number of timescales in use -- so
% one-timescale SFA adapts as hard in total as three-timescale SFA. Depression
% is deliberately not normalized that way, so its strength does grow with its
% timescale count; that asymmetry is intended, and
% docs/EquationsParametersDocs/Equations_stability_paper.md says why.
%
% NOT THE PAPER'S CONFIG YET. paper_config still names the 7-condition preset,
% and running this changes nothing about that. The network is IDENTICAL to the
% 4- and 7-condition presets -- only the regime set differs -- so results here
% are directly comparable with the exploratory seven-condition runs. The one
% thing to know when lining two runs up: sfa3_std2 here IS sfa_and_std there,
% byte-identical, named differently so each set can title its own comparison.
%
% IT CANNOT TOUCH THE PAPER'S OUTPUT. Both roots are named in the config, so
% nothing here writes into figs/paper or a dated run directory.
%
% RERUNNING REQUIRES DELETING data/single_multi_TS FIRST.
%
%   rmdir(fullfile(fileparts(which('setup_paths')), 'data', 'single_multi_TS'), 's')
%
% run_all_paper_analyses REFUSES a run directory that is not absent or empty,
% with TargetNotEmpty. A reused run directory accumulates: each sweep appends
% its own timestamped subfolder while the top-level manifest is replaced,
% leaving a manifest that describes one run in a directory holding several.
% figs/single_multi_TS is overwritten in place and needs no such care.
%
% A HANDOFF, NOT A LOOKUP. run_dir is named, so both halves address the same
% directory and there is nothing to race -- unlike reproduce_paper_run, where an
% empty run_dir means the figures resolve the newest matching run.
%
% ERROR ISOLATION. run_all_paper_analyses wraps each of its five stages, so one
% failure is reported and the queue continues. make_all_paper_figures reports
% per figure and counts success as files actually on disk. Read the two summary
% tables at the end; do not assume "it finished" means "everything worked".
%
% EXPECT SOME FIGURE LAYOUTS TO LOOK WRONG at three conditions -- anything that
% hardcodes seven columns will squash or blank. That is information, not a
% failure of this script.
%
% See also: single_multi_TS_config, paper_config, fast4_run,
%           run_all_paper_analyses, make_all_paper_figures

setup_paths();

cfg = single_multi_TS_config();

run_dir = run_all_paper_analyses(cfg);
results = make_all_paper_figures(cfg);

fprintf('\n========================================================\n');
fprintf('SINGLE vs MULTI TIMESCALE RUN COMPLETE\n');
fprintf('  preset  : %s (%s)\n', cfg.preset_name, cfg.run_mode);
fprintf('  run_dir : %s\n', run_dir);
fprintf('  figures : %d of %d succeeded\n', sum([results.ok]), numel(results));
fprintf('  to rerun: delete %s first\n', cfg.run_dir);
fprintf('========================================================\n');
