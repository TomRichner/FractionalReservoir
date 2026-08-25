Stability_Manuscript figure 1 panel A: chaos onset in a random network
======================================================================

Generated: 25-Aug-2026 16:35:27
By:        fig_introductory_concepts.m

WHAT IT SHOWS
  Two figures. statetraces: membrane-potential time series at three levels of
  synaptic gain. eigenspectra: the Jacobian eigenvalue disk at the same three
  gains, with the imaginary axis marking the stability boundary. Three
  networks share ONE underlying random weight matrix and differ only in gain,
  so the comparison isolates regulation from architecture.

HOW IT WAS MADE
  Reproduces Sompolinsky, Crisanti & Sommers (1988) on this project's model
  class. Purely random Dale-free Gaussian connectivity, tanh, fully connected,
  no adaptation and no external input, so the spectral radius is R = gamma
  exactly and chaos onset sits at gamma = 1. The Jacobian at x = 0 is J = (-I
  + W)/tau_d, since tanh'(0) = 1 and gamma is already folded into W.

SOURCE
  preset    sompolinsky_pairs
  net_seed  0
  gammas    [0.9 1.6 2.5]

PARAMETERS AS RUN
  n                     200
  indegree              200
  f                     1
  level_of_chaos        2.5
  tau_d                 1
  activation            tanh
  S_a                   0.9
  S_c                   0.4
  c                     0
  sigma_u_noise         0
  ode_solver            rk4
  fs                    200
  T_range               [0 60]
  mu_tilde_relative     0
  sigma_tilde_relative  1
  F_tracks_network      true
  cell_type_names       {all}
  n_a                   0
  rng_seeds             [0 1]
  x0_std                1
  R_theoretical         2.5
  N_sys_eqs             200
  spectral_radius       2.60265

FIGURES PRODUCED (in this folder)
  statetraces/panelA_bottom_traces.png
  statetraces/panelA_bottom_traces.svg
  statetraces/panelA_bottom_traces.fig
  eigenspectra/panelA_eigenspectrum.png
  eigenspectra/panelA_eigenspectrum.svg
  eigenspectra/panelA_eigenspectrum.fig

MEASURED EXPONENTS
  gamma = [0.9 1.6 2.5] -> LLE = [-0.0997 -0.00109 0.0989]

ONE CELL TYPE
  ONE CELL TYPE. This is a single undifferentiated population: zero-mean
  Gaussian weights taking both signs, so it is Dale-free and no E/I split
  exists to name. The type is called 'all' for that reason. Until 2026-08-23
  the preset used TWO statistically identical types named A and B, purely
  because SRNNCellTypePairs could not build a one-type model (build_network
  configured its RMTBlocks generator piecemeal where set_types is required to
  change the number of types). That is fixed, and the weight matrix is
  BIT-IDENTICAL across the change: the two types had identical zero-mean
  statistics, so the per-block scaling was uniform and the underlying draw
  never depended on how it was partitioned. This figure is therefore
  unchanged, measured exponents included. The spectral radius is R = gamma
  exactly either way.

