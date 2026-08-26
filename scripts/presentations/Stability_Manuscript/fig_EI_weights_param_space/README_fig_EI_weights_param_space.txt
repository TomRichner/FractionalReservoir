Stability_Manuscript figure: parameter-space distributions, E:I WEIGHT coloured
==============================================================================

Generated: 26-Aug-2026 17:39:21
By:        fig_EI_weights_param_space.m

WHAT IT SHOWS
  The same sheet as fig_EI_param_space -- one column per adaptation condition,
  row 1 the LLE distribution (green dashed zero line), row 2 the mean firing
  rate, each bar a stack of per-network patches -- but coloured by the balance
  of the synaptic WEIGHTS rather than by the fraction of neurons that are
  excitatory. The colorbar in the final column is labelled as an E:I weight
  ratio.

HOW IT WAS MADE
  Presentation replot; no simulation is re-run, but each grid point's NETWORK
  is rebuilt. A stored result records only scalars, never W, so
  psa.rebuild_model(res) reconstructs the constructor call the sweep made --
  including rng_seeds = [network_seed, network_seed + 1] -- and build()
  redraws the identical W. The colour value is sum(W(:,E)) / (sum(W(:,E)) +
  |sum(W(:,I))|), summed over presynaptic columns. Rebuilds are cached per
  config_idx, since the network belongs to the grid point and is shared by
  every condition run there.

SOURCE
  run_dir                C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\param_space\run_all_aug_25_26_22_08
  param_space_subfolder  param_space_nLevs_4_aug_26_26_00_33
  preset                 celltype_pairs_Sc0p2_noise0p025_dualStd_7cond
  clim                   [0.0909091 0.909091]

FIGURES PRODUCED (in this folder)
  Fig_EI_Weights_ParamSpace.png
  Fig_EI_Weights_ParamSpace.svg
  Fig_EI_Weights_ParamSpace.fig

READING THIS FIGURE
  WHY NOT f_E. Its sibling fig_EI_param_space colours by the fraction of
  neurons that are excitatory, which told the whole story while that was the
  only thing varying. The joint grid now also sweeps mu_EE / mu_EI / mu_IE /
  mu_II over -75% to +200%, so two networks with the same f_E can sit at
  opposite ends of the E:I balance: an 80% excitatory network with inhibitory
  synapses three times as strong is inhibition-dominated.

  MEASURED on the medium run of 2026-08-26, 64 grid points: the weight
  fraction spans [0.074, 0.980] where f_E spans [0.2, 0.8], and the two
  correlate at only 0.675. Roughly 1:12 to 48:1 in weights against a
  construction-capped 1:4 to 4:1 in neuron counts.

  COLOUR LIMITS are fixed at 1:10 to 10:1 ([1/11, 10/11]) rather than taken
  from the data, so the bar means the same thing in every run -- with a
  randomly subsampled grid, derived limits depend on which networks happened
  to be drawn. They were briefly 1:19 to 19:1, which covered nearly the whole
  measured spread but spent most of the colormap on a few extreme networks and
  flattened the contrast across the bulk near 1:1. Networks beyond 1:10 or
  10:1 clamp to the end colours and so read as "at least 10:1" rather than
  showing how far past it they go; load_and_make_unit_histograms warns with
  the count.

  The ratio is invariant to level_of_chaos and to rescale_by_abscissa: both
  are scalars multiplying all of W and cancel between numerator and
  denominator. What moves this axis is f and the four mu blocks.

