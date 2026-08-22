Stability_Manuscript figure: single-neuron adaptation mechanisms (sfa_std)
==========================================================================

Generated: 22-Aug-2026 11:12:19
By:        fig_adaptation_methods.m

WHAT IT SHOWS
  One unconnected neuron driven by a step in external input, under each
  mechanism combination as a column. Rows: the stimulus, the dendritic state
  x, the firing rate r, the synaptic output, and the active mechanism states.

HOW IT WAS MADE
  Built from the named preset with n_a and synapse_config overridden per
  column, so every column shares the preset's timescales and differs only in
  which mechanisms are on. Traces are read from model.plot_data directly
  rather than harvested from model.plot() axes.

SOURCE
  preset   celltype_pairs_Sc0p2_noise0p025_dualStd
  variant  sfa_std
  columns  {none, sfa, std, sfa+std}

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
  T_range               [0 25]
  mu_tilde_relative     [5.5 -5.5;5.5 -5.5]
  sigma_tilde_relative  [1.5 1.5;1.5 1.5]
  F_tracks_network      false
  cell_type_names       {E, I}
  n_a                   [3 0]
  rng_seeds             [1 2]
  x0_std                0.1
  R_theoretical         3.83333
  N_sys_eqs             1750
  spectral_radius       4.14766

FIGURES PRODUCED (in this folder)
  Fig_single_neuron_SFA_STD.png
  Fig_single_neuron_SFA_STD.svg
  Fig_single_neuron_SFA_STD.fig

WHY N = 2
  WHY n = 2 AND NOT 1. SRNNCellTypePairs enforces n >= n_cellTypes and rejects
  indegree = 0 (it requires 0 < indegree <= n), and it cannot build a
  one-cell-type model at all -- build_W assigns the RMTBlocks generator
  piecemeal where set_types is the only supported way to change the number of
  types, so a scalar f is expanded back to two types and the 1x1 mu_tilde then
  fails validation. Two neurons with identically zero weights is the smallest
  expressible unconnected network. Only the E neuron is plotted, and with W ==
  0 the second neuron cannot influence it.

STF MAPPING
  This variant shows SFA and STD only. The facilitation columns are in the
  sfa_std_stf variant, which uses the single_neuron_stf preset.

