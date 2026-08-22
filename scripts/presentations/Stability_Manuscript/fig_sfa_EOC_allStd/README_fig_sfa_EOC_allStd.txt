Stability_Manuscript figure: SFA edge of chaos (tau_a sensitivity)
==================================================================

Generated: 22-Aug-2026 11:53:20
By:        fig_sfa_EOC_allStd.m

WHAT IT SHOWS
  Single panel: the largest Lyapunov exponent lambda_1 against the maximum SFA
  adaptation timescale tau_a (the last and largest of the log-spaced tau_a_E
  vector), under SFA + STD.

HOW IT WAS MADE
  Presentation replot -- no simulation is re-run. The saved tau_a_E_max PSA
  object is reloaded via ParamSpaceAnalysis2.from_dir and plot_sensitivity
  rebuilds the panel, which is then restyled: imagesc CLim capped at
  total_reps*clim_frac; colormap white to 90 percent black so the blue median
  stays visible; axis box removed; the condition title dropped; x relabelled
  to max tau_a (s).

SOURCE
  run_dir        C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\param_space\run_all_aug_21_26_17_36
  tau_subfolder  tau_sensitivity_tau_timescales_tau_a_E_max_nLevs_7_aug_21_26_17_40
  preset         celltype_pairs_Sc0p2_noise0p025_dualStd

PARAMETERS AS RUN
  lle_range     [-0.3 0.1]
  y_view        [-0.25 0.05]
  y_ticks       [-0.2 -0.1 0 0.1]
  clim_frac     0.8
  median_alpha  0.35

FIGURES PRODUCED (in this folder)
  Fig_SFA_EOC_allStd.png
  Fig_SFA_EOC_allStd.svg
  Fig_SFA_EOC_allStd.fig

READING THIS FIGURE
  THE HISTOGRAM RANGE AND Y-VIEW CROP REAL DATA, and are kept identical to the
  original figure by choice. On the aug_14 preset the tau LLEs spanned about
  -0.26 to +0.29 with a MEDIAN of +0.008, so slightly over half the
  distribution was positive and about a third of it sat above the y_view
  ceiling -- in the +inf overflow band, which is itself off-screen. The panel
  therefore reads as almost entirely sub-zero when it is not. Widen lle_range,
  y_view and y_ticks together to show the full spread. Re-check those numbers
  against the run actually plotted; they are preset-dependent.

