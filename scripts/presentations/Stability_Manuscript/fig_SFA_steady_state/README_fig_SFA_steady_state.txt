Stability_Manuscript figure: SFA steady state, 1 vs 3 timescales
===============================================================

Generated: 19-Aug-2026 18:06:13
By script: Fig_SFA_steady_state.m

WHAT IT SHOWS
  Conceptual, ANALYTIC figure (no simulation). Setting da/dt = 0 in
    da_k/dt = (-a_k + r)/tau_k
  gives a_k = r for EVERY timescale, so c*sum(a) = c*n_a*r.

  Left  : c*sum(a) vs r, for n_a = 3 (c = 0.1667) and n_a = 1 (c = 0.5).
  Right : the transient response to a step r: 0 -> 1.

THE RESULT: THE TWO LINES IN THE LEFT PANEL ARE IDENTICAL
  Both are 0.5*r exactly. Because a_k = r independent of tau_k, the
  steady-state SFA current depends only on the PRODUCT c*n_a -- and
  splitting c as a budget (0.5/3 over three timescales vs 0.5 over one)
  holds that product fixed. The dashed line is drawn over the solid one
  to show both were computed; they coincide to machine precision.

  This is why the project uses c_E = 0.15/3: adding timescales does not
  silently move the operating point.

CONTRAST WITH STD (see fig_STD_steady_state)
  SFA enters the dynamics as a SUM, so a budget-split c makes the count
  invisible at steady state. STD enters as a PRODUCT, so n_b = 2 squares
  the depression no matter how the taus are chosen -- there is no
  budget-split of tau that would make dual STD match single STD.

WHERE THE COUNT DOES SHOW UP: THE TRANSIENT
  tau_a (n_a=3) = [0.25 1.5811 10]
  tau_a (n_a=1) = 0.25   (= tau_a(1), the FAST component of the n_a=3 set)
  The model rule logspace(log10(0.25), log10(10), n_a) is deliberately
  NOT used at n_a = 1: logspace(a, b, 1) returns 10^b, the slow end, so
  the comparison would be against the SLOWEST timescale alone and would
  be a statement about the slow tail rather than about the count.
  Three timescales give a multi-exponential approach -- fast partial
  adaptation followed by a slow tail -- against a single fast exponential.
  Same destination, different route. Both settle to 0.5000.

FIGURES PRODUCED (in this folder)
  Fig_SFA_steady_state.png
  Fig_SFA_steady_state.svg
  Fig_SFA_steady_state.fig
