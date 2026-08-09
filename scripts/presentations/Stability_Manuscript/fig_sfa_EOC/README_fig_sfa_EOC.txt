Stability_Manuscript figure: SFA Edge of Chaos (tau_a sensitivity)
=================================================================

Generated: 08-Aug-2026 21:39:02
By script: Fig_sfa_EOC.m

HOW IT WAS MADE
  Presentation replot -- no simulation is re-run. The script reloads the
  saved tau_a_E_max PSA object (psa_object.mat) and calls
  ParamSpaceAnalysis2.plot_sensitivity('metric','LLE') to rebuild the
  single-condition (SFA+STD) LLE-vs-tau_a panel, then restyles and saves
  it here. See git_provenance.txt for the exact commit.

SOURCE RUN
  C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\param_space\run_all_jul_06_26_22_00
  tau sensitivity subfolder used:
    tau_sensitivity_tau_timescales_tau_a_E_max_nLevs_15_jul_07_26_08_49

FIGURES PRODUCED (in this folder)
  Fig_SFA_EOC.png
  Fig_SFA_EOC.svg
  Fig_SFA_EOC.fig

  Single panel: largest Lyapunov exponent (lambda_1) vs the maximum
  SFA adaptation timescale tau_a (the last, largest of 3 log-spaced
  tau_a_E elements, swept 5..60 s). As tau_a grows lambda_1 rises
  toward 0 (the green dashed edge-of-chaos line) but stays negative.
  x-axis relabelled tau_a_E(last) -> "max \tau_a (s)"; ylabel kept as
  \lambda_1 (latex, matching the xlabel); condition title ("SFA + STD")
  removed.
  imagesc CLim capped at total_reps*0.4; colormap white -> 90% black so
  the blue median line stays visible. Blue median: alpha 0.35, width 3;
  green zero line width 2; axis box removed. LLE range [-0.3, 0.1].
