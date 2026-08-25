Stability_Manuscript figure: Jacobian eigenvalue occupancy
==========================================================

Generated: 25-Aug-2026 16:42:02
By:        fig_eig_heatmap.m

WHAT IT SHOWS
  A Gaussian-smoothed 2-D density over the complex plane for each adaptation
  regime -- how much time the instantaneous Jacobian's eigenvalues spend in
  each region, and in particular to the RIGHT of the imaginary axis, where the
  local dynamics are unstable. Panels share axis limits and one log-density
  colorbar.

HOW IT WAS MADE
  Plotting half only. run_eig_heatmap samples the Jacobian at fixed intervals
  through each run, pools the eigenvalues, and saves them; this renders that
  file, so the look can be iterated without re-simulating. All four regimes
  share rng_seeds, so W is identical and they are comparable.

SOURCE
  data_file  C:\Users\m218089\Desktop\github_repos\FractionalReservoir\scripts\presentations\Stability_Manuscript\fig_eig_heatmap\eig_heatmap_data.mat
  preset     celltype_pairs_Sc0p2_noise0p025_dualStd
  run_mode   fast

PARAMETERS AS RUN
  preset_name     celltype_pairs_Sc0p2_noise0p025_dualStd
  run_mode        fast
  T_range         [0 40]
  fs              200
  n_samples       40
  lle_window      20
  lya_T_interval  [20 40]
  model_class     SRNNCellTypePairs
  level_of_chaos  1
  n               500
  ode_solver      sra1

FIGURES PRODUCED (in this folder)
  fig_eig_heatmap.png
  fig_eig_heatmap.svg
  fig_eig_heatmap.fig

MEASURED EXPONENTS
  Finite-time Benettin exponents over the last 20 s of each run: No adaptation
  +2.6461, SFA only +2.4944, STD only -0.3045, SFA + STD -0.1149.

THE GAIN
  SYNAPTIC GAIN IS THE PRESET'S (level_of_chaos = 1), not the 3.0 the original
  script set. That 3.0 existed because the original was DETERMINISTIC and
  needed high gain to make the eigenvalues wander at all -- its own comment
  said so. The paper's preset is stochastic, and the noise is what moves the
  state, and therefore the Jacobian, around. Using the preset's own gain keeps
  this panel showing the same network as every other figure. If the cloud
  turns out to be static at this gain, that is a real result about the network
  rather than a reason to raise the gain silently.

