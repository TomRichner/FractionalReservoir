Stability_Manuscript figures: LLE and mean firing rate sensitivity
==================================================================

Generated: 26-Aug-2026 17:45:18
By:        fig_sensitivity_analysis_allStd.m

WHAT IT SHOWS
  One row per swept parameter, one column per adaptation condition. Each panel
  is an imagesc of the full across-reps distribution at each level, with the
  median overlaid as a bright blue line at 50% transparency, so it marks the
  central tendency without hiding the cells it crosses. Two sheets per metric:
  core (f_E, level_of_chaos, n) and mu (the four connectivity blocks).

HOW IT WAS MADE
  Presentation replot -- no simulation is re-run. replot_sensitivity reloads
  the saved PSA objects and regenerates both metrics;
  assemble_sensitivity_figure stacks them per sheet; the result is restyled
  and saved here, and the prep folder is deleted.

SOURCE
  run_dir  C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\param_space\run_all_aug_25_26_22_08
  preset   celltype_pairs_Sc0p2_noise0p025_dualStd_7cond

PARAMETERS AS RUN
  lle_range      [-1.75 1.75]
  n_bins         24
  clim_frac_LLE  0.8

FIGURES PRODUCED (in this folder)
  Fig_Sensitivity_LLE_core.png
  Fig_Sensitivity_LLE_core.svg
  Fig_Sensitivity_LLE_core.fig
  Fig_Sensitivity_LLE_mu.png
  Fig_Sensitivity_LLE_mu.svg
  Fig_Sensitivity_LLE_mu.fig
  Fig_sensitivity_mean_rate_core.png
  Fig_sensitivity_mean_rate_core.svg
  Fig_sensitivity_mean_rate_core.fig
  Fig_sensitivity_mean_rate_mu.png
  Fig_sensitivity_mean_rate_mu.svg
  Fig_sensitivity_mean_rate_mu.fig

AXIS CONVENTIONS
  x-axes are relabelled per parameter: f_E to E:I ratio; level_of_chaos to
  Synaptic Gain; n to Network Size; mu_XY_relative to mu_{XY} (RMT block
  means, indexed (post, pre)). The Synaptic Gain and the four mu axes are
  shown as PERCENT DEPARTURE from the preset default, (value/default - 1)*100,
  because absolute mu_tilde_relative values mean little on their own; this
  puts the preset own network at 0 percent and makes the four mu panels
  directly comparable. SIGN: mu_EI and mu_II have NEGATIVE defaults, so +100
  percent means twice as inhibitory, which in raw coordinates is further LEFT
  -- those two panels have their x-direction reversed so that on all four,
  rightward means "stronger synapse of this type". Only the ruler changes; the
  underlying image is untouched.

DEFAULT MARKER
  Every row carries a short reddish-grey tick rising from the x-axis at the
  preset default for that parameter -- the network the sweep departs from. On
  the percent axes it sits at 0 percent; on the f_E and Network Size rows
  there is no 0 percent tick, which is where it earns its keep. A default
  outside the swept range is NOT drawn rather than clamped to the edge. The
  values come from preset_default_values, which reads the run own
  run_manifest.mat, constructs a model from that preset and reads each value
  off the object, so the CLASS accessors resolve the aliases.
  resolved_defaults cannot serve here: ParamSpaceAnalysis2 excludes grid axes
  from it, and each of these parameters is the axis of one of the sweeps.

PER-METRIC DIFFERENCES
  Both metrics use the same bin count so they have the same vertical
  resolution. The LLE sheets keep the green dashed zero line -- the sign of
  lambda_1 marks the edge of chaos -- and the solid bands at top and bottom
  are the overflow bins. On the aug_14 preset the LLEs spanned roughly p1 =
  -10.0 to p99 = +3.7, so those bands carried a large share of the
  distribution; re-check that against the run actually plotted. The rate
  sheets use range [0, 1], where nothing can fall outside, so their overflow
  bins are always empty and the zero line is dropped as meaningless for a
  rate.

