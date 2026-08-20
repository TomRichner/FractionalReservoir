Stability_Manuscript figure: dual-timescale STD steady state
===========================================================

Generated: 20-Aug-2026 10:53:16
By script: Fig_STD_steady_state.m

WHAT IT SHOWS
  Conceptual, ANALYTIC figure (no simulation). Setting db/dt = 0 in
    db_k/dt = (1 - b_k)/tau_rec_k - b_k*r/tau_rel_k
  gives b_k(r) = 1/(1 + r*tau_rec_k/tau_rel_k); the synapse multiplies
  the timescales, so it sees prod_k b_k(r).

  Left  : prod(b) vs r, with one component b_k dashed for reference.
          Both timescales share tau_rec/tau_rel here, so b_1 and b_2
          coincide at steady state and the dashed curve is either.
  Right : prod(b)*r vs r -- what the recurrent sum actually receives,
          again with one component b_k*r dashed for reference. The
          faint green identity line y = r is the UNDEPRESSED synapse
          (b == 1); the vertical gap down to a curve is what
          depression costs at that rate.

TWO VERSIONS, differing ONLY in the right panel's y range
  Fig_STD_steady_state      : y over [0, 1] with equal pixel scaling on both
      axes (daspect [1 1 1]), so the identity line is a literal 45
      degrees and a vertical distance can be compared against a
      horizontal one by eye. Honest but small: both curves sit in the
      bottom tenth of the panel.
  Fig_STD_steady_state_zoom : y over [0, 0.12] -- just above the largest value a
      SINGLE b_k ever delivers (max b_k*r = 0.1111, at r = 1). The 1:1
      aspect ratio is DELIBERATELY DROPPED: it cannot survive a y range
      9x shorter than x's without squashing the panel to a sliver, so
      no cross-axis angle in this version means anything and the
      identity line leaves the top of the panel at r = 0.12. What it
      buys is the turnover, which is the point of the panel and is
      nearly invisible at full scale.
  The left panel is identical in both: prod(b) is a gain, already
  bounded by [0, 1], so there is nothing to zoom into.

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
  Fig_STD_steady_state_zoom.png
  Fig_STD_steady_state_zoom.svg
  Fig_STD_steady_state_zoom.fig
