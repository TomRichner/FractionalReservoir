Stability_Manuscript figure: SFA Edge of Chaos (tau_a sensitivity, allStd2)
=========================================================================

Generated: 14-Aug-2026 12:55:36
By script: Fig_sfa_EOC_allStd2.m

HOW IT WAS MADE
  Presentation replot -- no simulation is re-run. The script reloads the
  saved tau_a_E_max PSA object via ParamSpaceAnalysis2.from_dir and calls
  plot_sensitivity('metric','LLE') to rebuild the single-condition
  (SFA+STD) LLE-vs-tau_a panel, then restyles and saves it here. See
  git_provenance.txt for the exact commit.

SOURCE RUN
  C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\param_space\run_all_aug_14_26_00_48
  tau sensitivity subfolder used:
    tau_sensitivity_tau_timescales_tau_a_E_max_nLevs_11_aug_14_26_02_44

FIGURES PRODUCED (in this folder)
  Fig_SFA_EOC_allStd2.png
  Fig_SFA_EOC_allStd2.svg
  Fig_SFA_EOC_allStd2.fig

  Single panel: largest Lyapunov exponent (lambda_1) vs the maximum
  SFA adaptation timescale tau_a (the last, largest of 3 log-spaced
  tau_a_E elements). x-axis relabelled tau_a_E(last) -> "max \tau_a (s)";
  ylabel kept as \lambda_1 (latex, matching the xlabel); condition title
  ("SFA + STD") removed.
  imagesc CLim capped at total_reps*0.4; colormap white -> 90% black so
  the blue median line stays visible. Blue median: alpha 0.35, width 3;
  green zero line width 2; axis box removed. LLE histogram range
  [-0.3, 0.1], view cropped to [-0.25, 0.05].

READING THIS PANEL
  The histogram range and y-view are kept identical to the original
  figure by choice, and that crops part of the data: the tau LLEs span
  about -0.17 to +0.21 with a MEDIAN of -0.017. The median sits inside
  the view, but the upper tail (roughly the top 5%) lands in the +inf
  overflow band, which is itself above the y_view ceiling of 0.05 and so
  not drawn. To show the full spread, widen lle_range, y_view and
  y_ticks together at the top of the script.

COMPARING WITH Fig_SFA_EOC_allStd
  This run predates commit e429331 and sweeps tau_a_E(last) over the
  ORIGINAL [5, 60] s; the newer allStd run uses the retargeted [1, 30] s.
  The x-axes therefore cover different ranges and the two panels are not
  directly superimposable.
