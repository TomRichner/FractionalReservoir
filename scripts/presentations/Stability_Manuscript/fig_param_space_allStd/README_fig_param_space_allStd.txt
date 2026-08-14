Stability_Manuscript figure: Parameter-space distributions (allStd)
==================================================================

Generated: 14-Aug-2026 13:58:40
By script: Fig_param_space_allStd.m

HOW IT WAS MADE
  Presentation replot -- no simulation is re-run. The script reloads the
  saved param-space PSA object from a run_all_<dt> run, calls
  replot_param_space_analysis (psa.plot for LLE + mean_rate), then copies
  the per-condition histogram axes into a single 2x4 figure:
    row 1 = LLE distributions (green dashed zero line)
    row 2 = mean firing-rate distributions
    columns = No Adaptation, SFA, STD, SFA+STD
  Cleanups: LLE row xlabel "LLE (lambda_1)" -> \lambda_1; condition
  titles only on the top row (not bold); vertical gray column dividers;
  fonts matched to the MC/sensitivity figures; y-axes linked within each
  row, with probability ticks at 0/0.4/0.8 (LLE row) and 0/0.2/0.4
  (rate row). See git_provenance.txt for the exact commit.

SOURCE RUN
  C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\param_space\run_all_aug_14_26_12_14
  param_space subfolder(s) used:
    param_space_test_refactor_nLevs_4_aug_14_26_13_35

FIGURE PRODUCED (in this folder)
  Fig_ParamSpace_allStd.png / .svg / .fig

READING THIS FIGURE
  The LLE histogram range is fixed at [-1.5, 1.5] inside
  ParamSpaceAnalysis2.plot and cannot be set from this script. On this
  preset the param-space LLEs span roughly -10 to +4.8, so the outermost
  bins accumulate a large share of the distribution.
