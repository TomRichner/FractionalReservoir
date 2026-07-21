Stability_Manuscript figure: LLE Sensitivity (combined)
=======================================================

Generated: 21-Jul-2026 15:10:02
By script: Fig_sensitivity_analysis.m

HOW IT WAS MADE
  Presentation replot -- no simulation is re-run. The script reloads the
  saved 1D-sensitivity PSA objects from a run_all_<dt> run and calls
  replot_sensitivity -> assemble_sensitivity_figure to rebuild the stacked
  LLE figure (rows = swept params, cols = adaptation conditions), then
  saves it here. See git_provenance.txt for the exact commit.

SOURCE RUN
  C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\param_space\run_all_jul_06_26_22_00
  1D_sensitivity subfolders used:
    1D_sensitivity_sensitivity_f_nLevs_15_jul_07_26_04_56
    1D_sensitivity_sensitivity_level_of_chaos_nLevs_15_jul_07_26_06_51
    1D_sensitivity_sensitivity_n_nLevs_15_jul_06_26_22_00

FIGURES PRODUCED (in this folder)
  Fig_Sensitivity_LLE.png
  Fig_Sensitivity_LLE.svg
  Fig_Sensitivity_LLE.fig

  Combined LLE sensitivity: one row per swept parameter
  (f, level_of_chaos, n), one column per adaptation condition.
  x-axes relabelled: f -> "E:I ratio" (E:I = f:(1-f), ticks
  1:3, 1:2, 2:3, 1:1, 3:2, 2:1, 3:1); level_of_chaos -> "Synaptic Gain";
  n -> "Network Size". ylabel lambda_1 -> "Growth Rate"; larger tick fonts.
  Condition titles kept only on the top row; vertical gray dividers
  separate the 4 condition columns. imagesc CLim capped at
  total_reps*0.5 (shared); colormap white -> 90% black so the blue
  median line stays visible over the darkest cells. Panel letters
  (a)..(l) added up-and-left of each plot (AddLetters2Plots).
  Blue median line: alpha 0.35, 25% thinner; green zero line thinner.
  Titles not bold; axis boxes removed (x/y axes + ticks kept).
  LLE histograms: range [-1.75, 1.75], 22 bins.
