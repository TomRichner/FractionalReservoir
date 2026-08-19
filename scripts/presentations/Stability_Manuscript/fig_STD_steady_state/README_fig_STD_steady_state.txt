Stability_Manuscript figure: dual-timescale STD steady state
===========================================================

Generated: 19-Aug-2026 15:04:04
By script: Fig_STD_steady_state.m

WHAT IT SHOWS
  Conceptual, ANALYTIC figure (no simulation). Setting db/dt = 0 in
    db_k/dt = (1 - b_k)/tau_rec_k - b_k*r/tau_rel_k
  gives b_k(r) = 1/(1 + r*tau_rec_k/tau_rel_k); the synapse multiplies
  the timescales, so it sees prod_k b_k(r).

  Left  : prod(b) vs r, with one component b_k dashed for reference.
  Right : prod(b)*r vs r -- what the recurrent sum actually receives,
          again with one component b_k*r dashed for reference. Drawn
          with equal pixel scaling on both axes (daspect [1 1 1]) so
          the faint green identity line y = r is a literal 45 degrees.
          That line is the UNDEPRESSED synapse (b == 1); the vertical
          gap down to a curve is what depression costs at that rate.

PRESET (taus are read from it, not hardcoded)
  celltype_pairs_Sc0p2_noise0p025_dualStd
  tau_rec = [2 4]
  tau_rel = [0.25 0.5]
  tau_rec/tau_rel = [8 8]  (the only combination that sets the steady state)

NUMBERS WORTH QUOTING
  prod(b) at r = 0.3 : 0.0865
  prod(b) at r = 1   : 0.0123   (a 81x reduction in synaptic gain)
  output prod(b)*r peaks at r = 0.125, value 0.0312
  Beyond that peak the output DECREASES: depression outruns the rate
  driving it, so firing harder delivers less to the network.
  A SINGLE timescale does not do this -- b_k*r rises monotonically to
  an asymptote of 0.1250 (= tau_rel/tau_rec), reaching 0.1111 at r = 1,
  i.e. 3.6x the product's peak. The turnover is specific to n_b > 1.

FIGURES PRODUCED (in this folder)
  Fig_STD_steady_state.png
  Fig_STD_steady_state.svg
  Fig_STD_steady_state.fig
