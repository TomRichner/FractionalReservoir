Stability_Manuscript figure: parameter-space distributions, E:I coloured
========================================================================

Generated: 22-Aug-2026 11:53:14
By:        fig_EI_param_space.m

WHAT IT SHOWS
  A 2 x 5 sheet. Columns 1-4 are the four adaptation conditions, row 1 the LLE
  distribution (green dashed zero line) and row 2 the mean firing rate. Each
  bar is a stack of per-network patches coloured by the fraction excitatory
  (blue = inhibition-dominated, red = excitation-dominated). Column 5 holds
  the colorbar, labelled as an E:I ratio; the lower-right cell is empty.

HOW IT WAS MADE
  Presentation replot -- no simulation is re-run.
  load_and_make_unit_histograms builds per-metric histograms whose bars are
  stacks of per-network patches, coloured through blue_gray_red_colormap and
  normalized as probability with the LLE range fixed at [-1.5, 1.5]. Those
  axes are copied into a single combined figure.

SOURCE
  run_dir                C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\param_space\run_all_aug_21_26_17_36
  param_space_subfolder  param_space_nLevs_3_aug_21_26_17_40
  preset                 celltype_pairs_Sc0p2_noise0p025_dualStd
  color_by               f_E

FIGURES PRODUCED (in this folder)
  Fig_EI_ParamSpace.png
  Fig_EI_ParamSpace.svg
  Fig_EI_ParamSpace.fig

READING THIS FIGURE
  COLOUR AXIS. The bars are coloured by 'f_E'. On SRNNCellTypePairs the scalar
  fraction excitatory is the alias f_E (exactly f(1)); the property f is a 1 x
  C row there and would break the colouring. The value is read through
  psa.effective_param, not res.config, because effective_param('f') on a Pairs
  run returns the class default rather than the swept value.

  This figure used to be built from run_all_jul_06_26_22_00, an older
  SRNNModel2 run, while its sibling grey sheet used a SRNNCellTypePairs run --
  two param-space figures in one manuscript computed from different models.
  Both now resolve the same run.

