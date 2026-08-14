Stability_Manuscript figure: SFA Edge of Chaos (tau_a sensitivity, allStd)
=========================================================================

Generated: 14-Aug-2026 12:41:33
By script: Fig_sfa_EOC_allStd.m

HOW IT WAS MADE
  Presentation replot -- no simulation is re-run. The script reloads the
  saved tau_a_E_max PSA object via ParamSpaceAnalysis2.from_dir and calls
  plot_sensitivity('metric','LLE') to rebuild the single-condition
  (SFA+STD) LLE-vs-tau_a panel, then restyles and saves it here. See
  git_provenance.txt for the exact commit.

SOURCE RUN
  C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\param_space\run_all_aug_14_26_12_04
  tau sensitivity subfolder used:
    tau_sensitivity_tau_timescales_tau_a_E_max_nLevs_7_aug_14_26_12_10

FIGURES PRODUCED (in this folder)
  Fig_SFA_EOC_allStd.png
  Fig_SFA_EOC_allStd.svg
  Fig_SFA_EOC_allStd.fig

  Single panel: largest Lyapunov exponent (lambda_1) vs the maximum
  SFA adaptation timescale tau_a (the last, largest of 3 log-spaced
  tau_a_E elements). x-axis relabelled tau_a_E(last) -> "max \tau_a (s)";
  ylabel kept as \lambda_1 (latex, matching the xlabel); condition title
  ("SFA + STD") removed.
  imagesc CLim capped at total_reps*0.4; colormap white -> 90% black so
  the blue median line stays visible. Blue median: alpha 0.35, width 3;
  green zero line width 2; axis box removed. LLE histogram range
  [-0.3, 0.1], view cropped to [-0.25, 0.05].

READING THIS PANEL -- IMPORTANT
  The histogram range and y-view are kept identical to the original
  figure by choice. On this preset that crops real data: the tau LLEs
  span about -0.23 to +0.26 with a MEDIAN of +0.002, so roughly half the
  distribution is positive, while the view stops at 0.05 and the +inf
  overflow band sits off-screen. The panel therefore looks entirely
  sub-zero when it is not. To show the full spread, widen lle_range and
  y_view (and y_ticks) at the top of the script.
