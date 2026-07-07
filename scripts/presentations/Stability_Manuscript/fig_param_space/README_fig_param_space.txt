Stability_Manuscript figure: Parameter-space distributions (combined)
====================================================================

Generated: 07-Jul-2026 15:52:48
By script: Fig_param_space.m

HOW IT WAS MADE
  Presentation replot -- no simulation is re-run. The script reloads the
  saved param-space PSA object from a run_all_<dt> run, calls
  replot_param_space_analysis (psa.plot for LLE + mean_rate), then copies
  the per-condition histogram axes into a single 2x4 figure:
    row 1 = LLE (growth-rate) distributions (green dashed zero line)
    row 2 = mean firing-rate distributions
    columns = No Adaptation, SFA, STD, SFA+STD
  Cleanups: condition titles only on the top row (not bold), vertical gray
  column dividers, fonts matched to the MC/sensitivity figures, y-axes
  linked within each row. See git_provenance.txt for the exact commit.

SOURCE RUN
  C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\param_space\run_all_jul_06_26_22_00
  param_space subfolder(s) used:
    param_space_test_refactor_nLevs_4_jul_07_26_10_44

FIGURE PRODUCED (in this folder)
  Fig_ParamSpace.png / .svg / .fig
