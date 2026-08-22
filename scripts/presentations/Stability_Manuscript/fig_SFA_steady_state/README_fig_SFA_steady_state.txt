Stability_Manuscript figure: SFA steady state, one timescale vs several
=======================================================================

Generated: 22-Aug-2026 11:49:23
By:        fig_SFA_steady_state.m

WHAT IT SHOWS
  Conceptual, ANALYTIC figure (no simulation). Setting da_k/dt = 0 in da_k/dt
  = (-a_k + r)/tau_k gives a_k = r for every timescale, so c*sum(a) = c*n_a*r.
  Left: c*sum(a) against r, for the preset's timescale count and for one.
  Right: the transient response to a step r: 0 -> 1.

SOURCE
  preset  celltype_pairs_Sc0p2_noise0p025_dualStd

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
  only on the PRODUCT c*n_a, and splitting c as a budget holds that product
  fixed. The dashed line is drawn over the solid one to show both were
  computed; they coincide to machine precision. This is why the project splits
  c over its timescales -- adding timescales does not silently move the
  operating point.

THE TRANSIENT
  WHERE THE COUNT SHOWS UP: THE TRANSIENT. tau_a (n_a = 3) = [0.25 1.581 10]
  against a single tau_a = 0.25, which is the FAST component of that set. The
  model rule logspace(log10(0.25), log10(10), n_a) is deliberately not used at
  n_a = 1: logspace(a, b, 1) returns 10^b, the slow end, so the comparison
  would be against the slowest timescale alone. Several timescales give a
  multi-exponential approach -- fast partial adaptation then a slow tail --
  against a single exponential. Both settle to 0.5000.

CONTRAST WITH STD
  CONTRAST WITH STD (fig_STD_steady_state). SFA enters the dynamics as a SUM,
  so a budget-split c makes the count invisible at steady state. STD enters as
  a PRODUCT, so two timescales square the depression however the taus are
  chosen -- there is no budget-split of tau that would make dual STD match
  single STD.

PROVENANCE
  PARAMETERS COME FROM THE PRESET. The c BUDGET is recovered as c(1)*n_a(1),
  because the preset stores the already-divided value, not the budget. tau_a
  is not in the preset at all -- it is the class default, auto-filled per n_a
  -- so it is read off a BUILT model.

