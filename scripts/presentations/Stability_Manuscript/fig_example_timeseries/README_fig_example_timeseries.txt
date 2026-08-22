Stability_Manuscript figure: example time series
================================================

Generated: 22-Aug-2026 10:45:37
By:        fig_example_timeseries.m

WHAT IT SHOWS
  One realization of the network the rest of the paper analyses: stimulus,
  dendritic state x, firing rate r, synaptic output, and the
  adaptation/depression states, under the SFA + STD condition. Simulated
  inline -- there is no saved run behind it.

HOW IT WAS MADE
  Model built from the named preset under one adaptation condition, then
  model.plot(). The integrator is chosen from the preset: a preset carrying
  sigma_u_noise > 0 requires a stochastic scheme (sra1), because this figure
  does not go through analysis_run_config, which is what selects one for the
  sweeps.

SOURCE
  preset       celltype_pairs_Sc0p2_noise0p025_dualStd
  model_class  SRNNCellTypePairs
  condition    sfa_and_std
  rng_seeds    [1 2]
  ode_solver   sra1

PARAMETERS AS RUN
  n                     500
  indegree              100
  f                     [0.5 0.5]
  level_of_chaos        1
  tau_d                 0.1
  activation            piecewise
  S_a                   0.8
  S_c                   0.2
  c                     [0.166667 0]
  sigma_u_noise         0.025
  ode_solver            sra1
  fs                    400
  T_range               [0 20]
  mu_tilde_relative     [5.5 -5.5;5.5 -5.5]
  sigma_tilde_relative  [1.5 1.5;1.5 1.5]
  F_tracks_network      false
  cell_type_names       {E, I}
  n_a                   [3 0]
  rng_seeds             [1 2]
  x0_std                0.1
  R_theoretical         3.83333
  N_sys_eqs             3250
  spectral_radius       4.14766

FIGURES PRODUCED (in this folder)
  fig_example_timeseries.png
  fig_example_timeseries.svg
  fig_example_timeseries.fig

READING THIS FIGURE
  PORTED from SRNNModel2. This figure used to show a different model from
  every other figure in the paper -- logistic nonlinearity, n = 300, no noise
  -- while the sweeps were SRNNCellTypePairs, piecewise, n = 500, with noise.
  It now takes the same preset as the sweeps.

