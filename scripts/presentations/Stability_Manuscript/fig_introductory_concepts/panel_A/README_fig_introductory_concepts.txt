Stability_Manuscript figure 1 panel A: chaos onset in a random network
======================================================================

Generated: 22-Aug-2026 11:15:37
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
  f                     [0.5 0.5]
  level_of_chaos        2.5
  tau_d                 1
  activation            tanh
  S_a                   0.9
  S_c                   0.4
  c                     [0 0]
  sigma_u_noise         0
  ode_solver            rk4
  fs                    200
  T_range               [0 60]
  mu_tilde_relative     [0 0;0 0]
  sigma_tilde_relative  [1 1;1 1]
  F_tracks_network      true
  cell_type_names       {A, B}
  n_a                   [0 0]
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

WHY TWO CELL TYPES
  WHY TWO CELL TYPES, NAMED A AND B. SRNNCellTypePairs cannot build a
  one-cell-type model: build_W assigns the RMTBlocks generator piecemeal where
  set_types is the only supported way to change the number of types, so a
  scalar f is expanded back to two and the 1x1 mu_tilde then fails validation.
  Two types with IDENTICAL zero-mean blocks is the same network and builds
  today. They are named A and B rather than E and I because they are
  statistically indistinguishable and the weights take both signs -- calling
  them excitatory and inhibitory would be a lie. The traces concatenate both
  types, which together are the whole network.

