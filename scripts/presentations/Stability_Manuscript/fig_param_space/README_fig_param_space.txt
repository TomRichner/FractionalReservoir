Stability_Manuscript figure: Parameter-space distributions
==========================================================

Generated: 07-Jul-2026 15:39:36
By script: Fig_param_space.m

HOW IT WAS MADE
  Presentation replot -- no simulation is re-run. The script reloads the
  saved param-space PSA object from a run_all_<dt> run and calls
  replot_param_space_analysis (psa.plot for LLE + mean_rate) to rebuild the
  distribution histograms (1 row, one column per adaptation condition),
  then saves them here. See git_provenance.txt for the exact commit.

SOURCE RUN
  C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\param_space\run_all_jul_06_26_22_00
  param_space subfolder(s) used:
    param_space_test_refactor_nLevs_4_jul_07_26_10_44

FIGURES PRODUCED (in this folder)
  Fig_ParamSpace_LLE.png / .svg / .fig
  Fig_ParamSpace_MeanRate.png / .svg / .fig

  Fig_ParamSpace_LLE      = LLE (growth-rate) distribution histograms
  Fig_ParamSpace_MeanRate = mean firing-rate distribution histograms
  (one column per condition: No Adaptation, SFA, STD, SFA+STD).
