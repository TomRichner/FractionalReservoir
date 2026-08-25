Stability_Manuscript figure: stimulus engages adaptation (bursting)
===================================================================

Generated: 25-Aug-2026 16:37:46
By:        fig_stim_engages_adaptation.m

WHAT IT SHOWS
  A small, sparse network driven by a DC staircase. Below the step it is
  quiet; above it, adaptation and depression are engaged and the population
  BURSTS. Figure 1 is the model's time-series summary; figure 2 the power
  spectrum of the mean dendritic potential, one trace per DC level; figure 3
  (optional) population synchrony per level.

HOW IT WAS MADE
  Built from the bursting_pairs preset with a DC-staircase stimulus plus a
  white-noise probe. The staircase sets the operating point; the noise probes
  how the network filters its input, which the per-level PSD reads out. Each
  level skips its first settling seconds before the analysis window.

SOURCE
  preset     bursting_pairs
  rng_seeds  [42 42]
  dc_levels  [0 0.4]
  hold_dur   50

PARAMETERS AS RUN
  n                     50
  indegree              10
  f                     [0.7 0.3]
  level_of_chaos        1
  tau_d                 0.1
  activation            piecewise
  S_a                   0.9
  S_c                   0.5
  c                     [0.5 0]
  sigma_u_noise         0
  ode_solver            rk4
  fs                    400
  T_range               [0 100]
  mu_tilde_relative     [3 -2;3 -2]
  sigma_tilde_relative  [1 1;1 1]
  F_tracks_network      true
  cell_type_names       {E, I}
  n_a                   [3 0]
  rng_seeds             [42 42]
  x0_std                0.1
  R_theoretical         1.97203
  N_sys_eqs             225
  spectral_radius       3.27169

FIGURES PRODUCED (in this folder)
  bursting_timeseries.png
  bursting_timeseries.svg
  bursting_timeseries.fig
  bursting_psd.png
  bursting_psd.svg
  bursting_psd.fig
  bursting_synchrony.png
  bursting_synchrony.svg
  bursting_synchrony.fig

TRACES COME FROM PLOT_DATA
  TRACES COME FROM model.plot_data, not from the raw state vector. The
  original indexed S_out by hand -- x_cols = (nE*n_a_E + nI*n_a_I + nE*n_b_E +
  nI*n_b_I) + (1:n) -- which encodes SRNNModel2's state packing and does not
  hold for SRNNCellTypePairs, whose b-states are per ROUTE rather than per
  population. It also broke silently whenever a condition changed an
  adaptation count. plot_data is keyed by cell type and route, so the same
  code serves either class and cannot go out of step with the packing.
  plot_deci = 1 keeps it at the full sample rate, which the PSD requires.

SYNCHRONY
  Population synchrony on the adapting population (n = 35), measured over the
  settled part of each hold. chi2 = var_t(mean_i M) / mean_i var_t(M_i): 1 is
  perfect synchrony, 1/n perfect asynchrony. The reported rho_bar = (chi2 -
  1/n)/(1 - 1/n) removes that floor. Computed on x rather than r because r
  clips at the nonlinearity ceiling in the burst regime and saturation
  distorts variance. Values: DC = [0 0.4], rho_bar(x) = [0.386 0.198],
  rho_bar(r) = [0.583 0.324].

