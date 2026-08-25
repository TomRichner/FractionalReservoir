Stability_Manuscript figure: SFA steady state, one timescale vs several
=======================================================================

Generated: 25-Aug-2026 16:36:47
By:        fig_SFA_steady_state.m

WHAT IT SHOWS
  Conceptual, ANALYTIC figure (no simulation). Setting da_k/dt = 0 in da_k/dt
  = (-a_k + r)/tau_k gives a_k = r for every timescale, so c*sum(a) = c*n_a*r.
  Left: c*sum(a) against r, for the preset's timescale count and for one.
  Right: the transient response to a step r: 0 -> 1.

SOURCE
  preset  celltype_pairs_Sc0p2_noise0p025_dualStd_7cond

PARAMETERS AS RUN
  tau_a            [0.25 1.58114 10]
  n_a              3
  c_budget         0.5
  c_per_timescale  0.166667
  tau_single       0.25
  c_single         0.5

FIGURES PRODUCED (in this folder)
  Fig_SFA_steady_state.png
  Fig_SFA_steady_state.svg
  Fig_SFA_steady_state.fig

THE RESULT
  THE TWO LINES IN THE LEFT PANEL ARE IDENTICAL: both are 0.5*r exactly.
  Because a_k = r independently of tau_k, the steady-state SFA current depends
  only on the PRODUCT c*n_a, and the model dividing c by n_a holds that
  product fixed. The dashed line is drawn over the solid one to show both were
  computed; they coincide to machine precision. This is why the model
  normalises c by its timescale count -- adding timescales does not silently
  move the operating point.

THE TRANSIENT
  WHERE THE COUNT SHOWS UP: THE TRANSIENT. tau_a (n_a = 3) = [0.25 1.581 10]
  against a single tau_a = 0.25, which is the FAST component of that set.
  srnn_sfa_timescales(1) returns that fast end deliberately: the
  logspace(log10(0.25), log10(10), n_a) auto-fill it replaced returned 10^b at
  n_a = 1, the slow end, so the comparison would have been against the slowest
  timescale alone. Several timescales give a multi-exponential approach --
  fast partial adaptation then a slow tail -- against a single exponential.
  Both settle to 0.5000.

CONTRAST WITH STD
  CONTRAST WITH STD (fig_STD_steady_state). SFA enters the dynamics as a SUM,
  so normalising c by the count makes that count invisible at steady state.
  STD enters as a PRODUCT, so two timescales square the depression however the
  taus are chosen -- no normalisation of tau would make dual STD match single
  STD, and none is applied.

PROVENANCE
  PARAMETERS COME FROM THE PRESET. c IS the total adaptation budget and is
  read directly; the model divides it by the number of timescales in use. It
  used to be stored already divided, so this figure multiplied it back by n_a
  -- which after that change tripled it. tau_a is condition-owned rather than
  a preset field, so it is read off a BUILT model.

