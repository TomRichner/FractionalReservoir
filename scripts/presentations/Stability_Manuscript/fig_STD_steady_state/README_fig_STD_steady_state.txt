Stability_Manuscript figure: multi-timescale STD steady state
=============================================================

Generated: 27-Aug-2026 18:58:36
By:        fig_STD_steady_state.m

WHAT IT SHOWS
  Conceptual, ANALYTIC figure (no simulation). Setting db_k/dt = 0 gives
  b_k(r) = 1/(1 + r*tau_rec_k/tau_rel_k), and the synapse multiplies the
  timescales, so it sees prod_k b_k(r). Left: prod(b) against r with one
  component dashed. Right: prod(b)*r -- what the recurrent sum actually
  receives -- with the faint identity line marking the undepressed synapse.

SOURCE
  preset  celltype_pairs_Sc0p2_noise0p025_dualStd_7cond
  route   E -> E

PARAMETERS AS RUN
  tau_rec               [2 4]
  tau_rel               [0.25 0.5]
  tau_rec_over_tau_rel  [8 8]
  n_b                   2

FIGURES PRODUCED (in this folder)
  Fig_STD_steady_state.png
  Fig_STD_steady_state.svg
  Fig_STD_steady_state.fig
  Fig_STD_steady_state_zoom.png
  Fig_STD_steady_state_zoom.svg
  Fig_STD_steady_state_zoom.fig

NUMBERS WORTH QUOTING
  prod(b) at r = 0.3 is 0.08651; at r = 1 it is 0.01235 -- a 81x reduction in
  synaptic gain relative to an undepressed synapse (b = 1). The output
  prod(b)*r PEAKS at r = 0.125 with value 0.03125, and DECREASES beyond it:
  depression outruns the rate driving it, so firing harder delivers LESS to
  the network. A SINGLE timescale does not do this -- b_k*r rises
  monotonically to an asymptote of 0.125 (= tau_rel/tau_rec), which is 4.0x
  the product's peak. The turnover is specific to more than one depression
  timescale.

CONTRAST WITH SFA
  CONTRAST WITH SFA (fig_SFA_steady_state). SFA enters the dynamics as a SUM,
  so splitting c as a budget over its timescales makes the count invisible at
  steady state. STD enters as a PRODUCT, so two timescales SQUARE the
  depression however the taus are chosen -- there is no budget-split of tau
  that would make dual STD match single STD.

THE TWO VERSIONS
  TWO VERSIONS, differing ONLY in the right panel y range. The full-scale one
  uses equal pixel scaling on both axes, so the identity line is a literal 45
  degrees and a vertical distance can be compared against a horizontal one by
  eye -- honest, but both curves sit in the bottom tenth of the panel. The
  zoom DELIBERATELY drops the 1:1 aspect, which cannot survive a much shorter
  y range without squashing the panel to a sliver, so no cross-axis angle
  there means anything; what it buys is the turnover. The left panel is
  identical in both: prod(b) is a gain, already bounded by [0, 1], so there is
  nothing to zoom into.

