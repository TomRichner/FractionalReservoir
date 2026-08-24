Stability_Manuscript figure: single-neuron adaptation mechanisms (sfa_std)
==========================================================================

Generated: 23-Aug-2026 22:40:26
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
  preset   single_neuron_dualStd
  variant  sfa_std
  columns  {none, sfa, std, sfa+std}

PARAMETERS AS RUN
  n                     1
  indegree              1
  f                     1
  level_of_chaos        1
  tau_d                 0.1
  activation            piecewise
  S_a                   0.8
  S_c                   0.2
  c                     0.166667
  sigma_u_noise         0.025
  ode_solver            sra1
  fs                    400
  T_range               [0 25]
  mu_tilde_relative     0
  sigma_tilde_relative  0
  F_tracks_network      true
  cell_type_names       {E}
  n_a                   3
  rng_seeds             [1 2]
  x0_std                0
  R_theoretical         0
  N_sys_eqs             6
  spectral_radius       0

FIGURES PRODUCED (in this folder)
  Fig_single_neuron_SFA_STD.png
  Fig_single_neuron_SFA_STD.svg
  Fig_single_neuron_SFA_STD.fig

THE NETWORK
  THE NETWORK IS ONE NEURON. Both variants run n = 1 with a single cell type
  and W = 0 -- no recurrence of any kind, so every trace is the mechanism
  responding to the external step alone. Check PARAMETERS AS RUN above: it is
  read off the built object, so it reports what was actually simulated. Until
  2026-08-23 variant A named the paper's n = 500 network preset and overrode
  nothing, so it plotted neuron 1 of a fully recurrent chaotic network while
  claiming to show one unconnected neuron; the mechanisms were invisible. The
  single-neuron networks now live in their own presets (single_neuron_dualStd,
  single_neuron_stf), which is what stops that recurring. Variant A inherits
  the paper's c, three SFA timescales and dual depression, so the cartoon
  explains the network figures rather than some other model.

STF MAPPING
  This variant shows SFA and STD only. The facilitation columns are in the
  sfa_std_stf variant, which uses the single_neuron_stf preset.

