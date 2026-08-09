Stability_Manuscript figures: LLE + mean firing rate Sensitivity (combined)
==========================================================================

Generated: 08-Aug-2026 19:40:13
By script: Fig_sensitivity_analysis.m

HOW THEY WERE MADE
  Presentation replot -- no simulation is re-run. The script reloads the
  saved 1D-sensitivity PSA objects from a run_all_<dt> run and calls
  replot_sensitivity (which plots BOTH metrics) -> assemble_sensitivity_figure
  once per metric to rebuild the stacked figures (rows = swept params,
  cols = adaptation conditions), then saves them here. See
  git_provenance.txt for the exact commit.

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
  Fig_sensitivity_mean_rate.png
  Fig_sensitivity_mean_rate.svg
  Fig_sensitivity_mean_rate.fig

SHARED LAYOUT (both figures)
  One row per swept parameter (f, level_of_chaos, n), one column per
  adaptation condition. x-axes relabelled: f -> "E:I ratio"
  (E:I = f:(1-f), ticks 1:3, 2:3, 1:1, 3:2, 3:1);
  level_of_chaos -> "Synaptic Gain"; n -> "Network Size". Larger tick
  fonts. Condition titles kept only on the top row; vertical gray
  dividers separate the 4 condition columns. imagesc CLim capped at
  total_reps*0.4 (shared within a figure); colormap white -> 90% black
  so the blue median line stays visible over the darkest cells. Panel
  letters (a)..(l) added up-and-left of each plot (AddLetters2Plots).
  Blue median line: alpha 0.35, 25% thinner. Titles not bold; axis
  boxes removed (x/y axes + ticks kept).

PER-FIGURE DIFFERENCES
  Both metrics use n_bins = 24 (linspace edges), i.e. 25 plotted rows:
  23 interior bins plus the two -inf/+inf overflow bins. Matching the
  bin count gives the two figures the same vertical resolution.

  Fig_Sensitivity_LLE:       ylabel \lambda_1; histogram range
                             [-1.75, 1.75]; green dashed zero line kept
                             (thinner) -- the sign of lambda_1 marks the
                             edge of chaos. The solid bands at the top and
                             bottom are the overflow bins (reps outside
                             the range).
  Fig_sensitivity_mean_rate: ylabel "Mean Firing Rate"; histogram range
                             [0, 1] (plot_sensitivity default; nothing can
                             fall outside it, so the overflow bins are
                             always empty); y ticks at 0 and 1 only; zero
                             line removed -- at y=0 it lands on the bottom
                             axis and carries no meaning for a rate.
