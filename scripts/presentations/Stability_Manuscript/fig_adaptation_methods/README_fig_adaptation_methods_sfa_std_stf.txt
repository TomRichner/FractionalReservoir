Stability_Manuscript figure: single-neuron adaptation mechanisms (sfa_std_stf)
==============================================================================

Generated: 22-Aug-2026 11:28:13
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
  n                     2
  indegree              1
  f                     [0.5 0.5]
  level_of_chaos        1
  tau_d                 0.1
  activation            piecewise
  S_a                   1
  S_c                   0.5
  c                     [1 0]
  sigma_u_noise         0
  ode_solver            rk4
  fs                    400
  T_range               [0 25]
  mu_tilde_relative     [0 0;0 0]
  sigma_tilde_relative  [0 0;0 0]
  F_tracks_network      true
  cell_type_names       {E, I}
  n_a                   [1 0]
  rng_seeds             [1 2]
  x0_std                0
  R_theoretical         0
  N_sys_eqs             5
  spectral_radius       0

FIGURES PRODUCED (in this folder)
  Fig_single_neuron_SFA_STD_STF.png
  Fig_single_neuron_SFA_STD_STF.svg
  Fig_single_neuron_SFA_STD_STF.fig

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

