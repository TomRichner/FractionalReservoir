Stability_Manuscript figure: sensitivity MEDIANS, collapsed across conditions
============================================================================

Generated: 18-Aug-2026 21:59:46
By script: Fig_sensitivity_medians.m

WHAT THIS IS
  The compact counterpart of ../fig_sensitivity_analysis_allStd. That figure
  tiles the four adaptation conditions as COLUMNS and shows the full rep
  distribution per panel as an imagesc histogram (28 panels over two sheets
  per metric). Here the conditions are OVERLAID in one axes per swept
  parameter, as median lines only, so all six parameters fit one 2 x 3 sheet
  and the conditions can be compared directly rather than across columns.
  No simulation is re-run.

SOURCE RUN
  C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\param_space\run_all_aug_14_26_17_25
  1D_sensitivity subfolders used:
    f_E             
    n               
    mu_EE_relative  
    mu_EI_relative  
    mu_IE_relative  
    mu_II_relative  

FIGURES PRODUCED (in this folder)
  Fig_Sensitivity_LLE_medians.png / .svg / .fig
  Fig_sensitivity_mean_rate_medians.png / .svg / .fig

PANELS (row-major, 2 x 3)
  (a) f_E               x-axis: E:I neuron ratio
  (b) n                 x-axis: Network Size
  (c) mu_EE_relative    x-axis: \mu_{EE}
  (d) mu_EI_relative    x-axis: \mu_{EI}
  (e) mu_IE_relative    x-axis: \mu_{IE}
  (f) mu_II_relative    x-axis: \mu_{II}

  level_of_chaos ("Synaptic Gain") is deliberately DROPPED: the allStd run
  sweeps seven parameters, which do not tile 2 x 3, and the gain sweep is the
  least surprising of the seven (it simply scales W).

E:I NEURON RATIO. The f_E sweep varies the fraction excitatory with mu_EI and
  mu_IE held fixed, so what the axis really reports is the E:I neuron COUNT.
  This run has n = 500, so the ticks read 100:400 / 250:250 / 400:100 rather
  than the reduced ratios 1:4 / 1:1 / 4:1, which would hide the network size
  the counts come from.

PERCENT AXES. The four mu axes are shown as percent departure from the preset
  default, (value/default - 1)*100. mu_EI and mu_II have NEGATIVE defaults, so
  "+100%" means twice as inhibitory; those two panels have their x-direction
  REVERSED so that on all four rightward means "stronger synapse of this type".
  Only the ruler changes -- the plotted data is untouched.

DEFAULT MARKER. Every panel carries a short reddish-gray tick rising from the
  x-axis (0.05 of the y-range, 2 pt) at the preset default for that parameter.
  Resolved by preset_default_values() from the run's own run_manifest.mat, not
  hardcoded. For this run:
    f_E              default 0.5
    n                default 500
    mu_EE_relative   default 5.5
    mu_EI_relative   default -5.5
    mu_IE_relative   default 5.5
    mu_II_relative   default -5.5

CONDITION COLOURS (Okabe-Ito, colorblind-safe; same as fig_memory_capacity_example)
  no_adaptation  No Adaptation  [0.000 0.000 0.000]
  sfa_only       SFA Only       [0.902 0.624 0.000]
  std_only       STD Only       [0.337 0.706 0.914]
  sfa_and_std    SFA + STD      [0.800 0.475 0.655]
  Conditions are matched to colours BY NAME, so a run declaring them in a
  different order cannot silently recolour the figure.

STATISTIC. Median across reps at each swept-parameter level, per condition
  (ParamSpaceAnalysis2.collect_level_values -> prctile), i.e. exactly the blue
  median line of the allStd sheets. Failed / NaN reps are excluded first. The
  script computes prctile(vals, pcts) with pcts = 50 and plots only the median;
  adding e.g. an IQR band means extending pcts and setting band_pcts, whose
  shading branch is already written.

PER-METRIC DIFFERENCES
  Fig_Sensitivity_LLE_medians: ylabel \lambda_1; y window [-1.75, 1.75], kept identical to the
    allStd sheets; green dashed zero line kept (the sign of lambda_1 marks the
    edge of chaos). NOTE: this preset's LLEs span roughly p1 = -10.0 to
    p99 = +3.7, so several medians -- No Adaptation above all -- leave the
    window and are CLIPPED rather than rescaling every panel around them.
  Fig_sensitivity_mean_rate_medians: ylabel "Mean Firing Rate"; y window [0, 1] (nothing can
    fall outside it); y ticks at 0 and 1 only; zero line removed.
