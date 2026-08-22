Stability_Manuscript figure: sensitivity medians, conditions overlaid
=====================================================================

Generated: 22-Aug-2026 10:59:59
By:        fig_sensitivity_medians.m

WHAT IT SHOWS
  One subplot per swept parameter, four median lines per subplot (one per
  adaptation condition), 2 x 3 on a single sheet per metric. The compact
  counterpart of the allStd sheets, which tile the conditions as columns
  instead.

HOW IT WAS MADE
  Medians are computed straight off the saved PSA objects via
  ParamSpaceAnalysis2.collect_level_values -- no replot detour, because
  harvesting median lines out of saved .fig files would be fragile.

SOURCE
  run_dir  C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\param_space\run_all_aug_21_26_17_36
  preset   celltype_pairs_Sc0p2_noise0p025_dualStd

FIGURES PRODUCED (in this folder)
  Fig_Sensitivity_LLE_medians.png
  Fig_Sensitivity_LLE_medians.svg
  Fig_Sensitivity_LLE_medians.fig
  Fig_sensitivity_mean_rate_medians.png
  Fig_sensitivity_mean_rate_medians.svg
  Fig_sensitivity_mean_rate_medians.fig

PANELS
  Panels, row-major over 2 x 3: f_E (E:I neuron ratio), n (Network Size), then
  the four connectivity blocks mu_EE, mu_EI, mu_IE, mu_II. level_of_chaos is
  deliberately dropped -- seven parameters do not tile 2 x 3, and the gain
  sweep is the least surprising of the seven.

E:I AXIS
  E:I NEURON RATIO. The f_E sweep varies the fraction excitatory with mu_EI
  and mu_IE held fixed, so what the axis really reports is the E:I neuron
  COUNT. The ticks are therefore spelled as counts rather than as reduced
  ratios, which would hide the network size the counts are drawn from.

PERCENT AXES
  PERCENT AXES. The four mu axes are shown as percent departure from the
  preset default, (value/default - 1)*100. mu_EI and mu_II have NEGATIVE
  defaults, so +100 percent means twice as inhibitory; those two panels have
  their x-direction reversed so that on all four, rightward means stronger
  synapse of this type. Only the ruler changes -- the plotted data is
  untouched. Every panel also carries a short reddish-grey tick at the preset
  default, resolved by preset_default_values from the run own run_manifest.mat
  rather than hardcoded.

STATISTIC
  STATISTIC. Median across reps at each swept level, per condition
  (ParamSpaceAnalysis2.collect_level_values then prctile) -- exactly the blue
  median line of the allStd sheets. Failed or NaN reps are excluded first.
  Condition colours come from manuscript_style and are matched BY NAME, so a
  run declaring its conditions in a different order cannot silently recolour
  the figure.

CLIPPING
  CLIPPING. The LLE panel keeps the same y window as the allStd sheets and the
  green dashed zero line (the sign of lambda_1 marks the edge of chaos). On
  the aug_14 preset several medians -- No Adaptation above all -- left that
  window and were CLIPPED rather than rescaling every panel around them.
  Re-check against the run actually plotted. The rate panel uses [0, 1], where
  nothing can fall outside, so its zero line is dropped as meaningless for a
  rate.

