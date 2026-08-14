Stability_Manuscript figure: SFA Edge of Chaos (tau_a sensitivity, allStd)
=========================================================================

Generated: 14-Aug-2026 14:03:51
By script: Fig_sfa_EOC_allStd.m

HOW IT WAS MADE
  Presentation replot -- no simulation is re-run. The script reloads the
  saved tau_a_E_max PSA object via ParamSpaceAnalysis2.from_dir and calls
  plot_sensitivity('metric','LLE') to rebuild the single-condition
  (SFA+STD) LLE-vs-tau_a panel, then restyles and saves it here. See
  git_provenance.txt for the exact commit.

SOURCE RUN
  C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\param_space\run_all_aug_14_26_12_14
  tau sensitivity subfolder used:
    tau_sensitivity_tau_timescales_tau_a_E_max_nLevs_11_aug_14_26_13_24

FIGURES PRODUCED (in this folder)
  Fig_SFA_EOC_allStd.png
  Fig_SFA_EOC_allStd.svg
  Fig_SFA_EOC_allStd.fig

  Single panel: largest Lyapunov exponent (lambda_1) vs the maximum
  SFA adaptation timescale tau_a (the last, largest of 3 log-spaced
  tau_a_E elements). x-axis relabelled tau_a_E(last) -> "max \tau_a (s)";
  ylabel kept as \lambda_1 (latex, matching the xlabel); condition title
  ("SFA + STD") removed.
  imagesc CLim capped at total_reps*0.8; colormap white -> 90% black so
  the blue median line stays visible. Blue median: alpha 0.35, width 3;
  green zero line width 2; axis box removed. LLE histogram range
  [-0.3, 0.1], view cropped to [-0.25, 0.05].

READING THIS PANEL -- IMPORTANT
  The histogram range and y-view are kept identical to the original
  figure by choice. On this preset that crops real data: the tau LLEs
  span about -0.26 to +0.29 with a MEDIAN of +0.008, so slightly over
  half the distribution is positive and 32% of it sits above the view
  ceiling of 0.05 -- in the +inf overflow band, which is itself
  off-screen. The panel therefore reads as almost entirely sub-zero
  when it is not. To show the full spread, widen lle_range, y_view and
  y_ticks together at the top of the script.
