Stability_Manuscript figure: SFA and STD reshape the F-I curve
==============================================================

Generated: 25-Aug-2026 16:36:03
By:        fig_FI_curve.m

WHAT IT SHOWS
  Conceptual, ANALYTIC figure (no simulation): the synaptic output br = b *
  phi(x - c*a) against input x, with phi(x) = 1/(1+exp(-4*(x - S_c))) and S_c
  = 0.40.

  Left panel (SFA): b = 1, a swept 0 -> 1. Subtracting c*a inside phi shifts
  the curve HORIZONTALLY, raising the threshold.

  Right panel (STD): a = 0, b swept 1 -> 0. The b prefactor rescales the curve
  VERTICALLY, reducing gain and the saturation ceiling.

PARAMETERS AS RUN
  S_c        0.4
  c_concept  0.6
  a_levels   [0 0.25 0.5 0.75 1]
  b_levels   [1 0.75 0.5 0.25 0]

FIGURES PRODUCED (in this folder)
  Fig_FI_curve.png
  Fig_FI_curve.svg
  Fig_FI_curve.fig

READING THIS FIGURE
  c = 0.60 is ENLARGED for visibility; the model default c_E = 0.5/3 is too
  small to see. The functional form b*phi(x - c*a) is exact -- only the
  magnitude of c is chosen for clarity.

  This figure is deliberately NOT built from the paper's preset. It is a
  schematic of the form of the equations, not a picture of a particular
  network, so it keeps the logistic phi at S_c = 0.4 rather than the preset's
  piecewise nonlinearity.

