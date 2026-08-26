Stability_Manuscript figure: single-neuron adaptation mechanisms (sfa_std_stf)
==============================================================================

Generated: 26-Aug-2026 14:50:39
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
  preset   single_neuron_stf
  variant  sfa_std_stf
  columns  {none, sfa, std, stf, sfa+std, std+stf, sfa+std+stf}

PARAMETERS AS RUN
  n                     1
  indegree              1
  f                     1
  level_of_chaos        1
  tau_d                 0.1
  activation            piecewise
  S_a                   1
  S_c                   0.5
  c                     1
  sigma_u_noise         0
  ode_solver            rk4
  fs                    400
  T_range               [0 25]
  mu_tilde_relative     0
  sigma_tilde_relative  0
  F_tracks_network      true
  cell_type_names       {E}
  n_a                   1
  rng_seeds             [1 2]
  x0_std                0
  R_theoretical         0
  N_sys_eqs             4
  spectral_radius       0

FIGURES PRODUCED (in this folder)
  Fig_single_neuron_SFA_STD_STF.png
  Fig_single_neuron_SFA_STD_STF.svg
  Fig_single_neuron_SFA_STD_STF.fig

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
  STF PARAMETERS map EXACTLY from the deleted test_single_neuron_stf.m. That
  model carried a release probability p resting at p0 with gain p/p0 and dp/dt
  = (p0-p)/tau_f + kappa*(1-p)*r. Substituting the gain variable u = p/p0
  gives du/dt = (1-u)/tau_f + kappa*(1/p0 - u)*r, which IS this class's dg/dt
  = (1-g)/tau_dec + (G-g)*r/tau_fac with tau_dec = tau_f, tau_fac = 1/kappa
  and G = 1/p0. Here tau_dec = 6, tau_fac = 2.5, G = 2.857, i.e. the old tau_f
  = 6, kappa = 0.4, p0 = 0.35. WHAT DOES NOT MAP: the old STD depleted in
  proportion to p (db/dt = (1-b)/tau_rec - (p*b*r)/tau_rel, the
  Tsodyks-Markram coupling), where this class's depletion is independent of g.
  No tau_rel reproduces that -- it is a constant factor at rest but
  ACCELERATES as p rises. tau_rel is carried literally, so this figure will
  NOT match the archived PNG (about 2.9x stronger depression at rest). That is
  expected. Whether the decoupling is correct at all is recorded in
  UserNotes.md.

