Stability_Manuscript figures: LLE + mean firing rate Sensitivity (allStd)
========================================================================

Generated: 14-Aug-2026 12:41:26
By script: Fig_sensitivity_analysis_allStd.m

HOW THEY WERE MADE
  Presentation replot -- no simulation is re-run. The script reloads the
  saved 1D-sensitivity PSA objects from a run_all_<dt> run and calls
  replot_sensitivity (which plots BOTH metrics) -> assemble_sensitivity_figure
  once per metric per sheet to rebuild the stacked figures (rows = swept
  params, cols = adaptation conditions), then saves them here. See
  git_provenance.txt for the exact commit.

SOURCE RUN
  C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\param_space\run_all_aug_14_26_12_04
  1D_sensitivity subfolders used:
    1D_sensitivity_sensitivity_f_E_nLevs_4_aug_14_26_12_05
    1D_sensitivity_sensitivity_level_of_chaos_nLevs_4_aug_14_26_12_06
    1D_sensitivity_sensitivity_mu_EE_relative_nLevs_4_aug_14_26_12_07
    1D_sensitivity_sensitivity_mu_EI_relative_nLevs_4_aug_14_26_12_07
    1D_sensitivity_sensitivity_mu_IE_relative_nLevs_4_aug_14_26_12_08
    1D_sensitivity_sensitivity_mu_II_relative_nLevs_4_aug_14_26_12_09
    1D_sensitivity_sensitivity_n_nLevs_4_aug_14_26_12_04

SHEETS
  This run is SRNNCellTypePairs, which sweeps SEVEN parameters rather than
  the three of the original SRNNModel2 figure. Seven rows on one sheet is
  unreadably tall, so each metric is split in two:
    core  : f_E, level_of_chaos, n
    mu    : mu_EE_relative, mu_EI_relative, mu_IE_relative, mu_II_relative

FIGURES PRODUCED (in this folder)
  Fig_Sensitivity_LLE_core.png / .svg / .fig
  Fig_Sensitivity_LLE_mu.png / .svg / .fig
  Fig_sensitivity_mean_rate_core.png / .svg / .fig
  Fig_sensitivity_mean_rate_mu.png / .svg / .fig

SHARED LAYOUT (all sheets)
  One row per swept parameter, one column per adaptation condition.
  x-axes relabelled per parameter: f_E -> "E:I ratio"
  (E:I = f_E:(1-f_E), ticks 1:4, 2:3, 3:2, 4:1); level_of_chaos ->
  "Synaptic Gain"; n -> "Network Size"; mu_XY_relative -> \mu_{XY}
  (RMT block means, indexed (post, pre)). Rows are identified by their
  index in the sheet's parameter list, NOT by inspecting the data -- the
  original figure guessed from max(xlim), which mislabels the mu sweeps.
  Larger tick fonts. Condition titles kept only on the top row; vertical
  gray dividers separate the condition columns. imagesc CLim capped at
  total_reps*0.4 (shared within a figure); colormap white -> 90% black
  so the blue median line stays visible over the darkest cells. Panel
  letters added up-and-left of each plot (AddLetters2Plots). Blue median
  line: alpha 0.35, 25% thinner. Titles not bold; axis boxes removed.

PER-FIGURE DIFFERENCES
  Both metrics use n_bins = 24 (linspace edges), i.e. 25 plotted rows:
  23 interior bins plus the two -inf/+inf overflow bins. Matching the
  bin count gives the two metrics the same vertical resolution.

  Fig_Sensitivity_LLE_*:       ylabel \lambda_1; histogram range
                               [-1.75, 1.75]; green dashed zero line kept
                               (thinner) -- the sign of lambda_1 marks the
                               edge of chaos. The solid bands at the top and
                               bottom are the overflow bins (reps outside
                               the range). NOTE: this preset's LLEs span
                               roughly p1 = -9.4 to p99 = +4.0, so those
                               overflow bands carry a large share of the
                               distribution. The range is kept identical to
                               the original figure by choice; widen
                               lle_range in the script to change that.
  Fig_sensitivity_mean_rate_*: ylabel "Mean Firing Rate"; histogram range
                               [0, 1] (plot_sensitivity default; nothing can
                               fall outside it, so the overflow bins are
                               always empty); y ticks at 0 and 1 only; zero
                               line removed -- at y=0 it lands on the bottom
                               axis and carries no meaning for a rate.
