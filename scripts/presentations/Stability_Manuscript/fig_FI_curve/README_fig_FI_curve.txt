Stability_Manuscript figure: SFA & STD reshape the F-I curve
===========================================================

Generated: 28-Jul-2026 12:48:22
By script: Fig_FI_curve.m

WHAT IT SHOWS
  Conceptual, ANALYTIC figure (no SRNNModel2 simulation): plots the
  synaptic output  br = b * phi(x - c*a)  vs input x, with
  phi(x) = 1/(1+exp(-4*(x - S_c))), S_c = 0.40 (class-default logistic).

  Left panel  (SFA): b = 1, a swept 0->1. Subtracting c*a inside phi
    shifts the curve HORIZONTALLY -> raises threshold / changes bias.
  Right panel (STD): a = 0, b swept 1->0. The b prefactor rescales the
    curve VERTICALLY -> reduces gain/slope and the saturation ceiling.

CONCEPTUAL PARAMETER
  c = 0.60 is enlarged for visibility (model default c_E = 0.15/3 ~ 0.05
  is too small to see). The functional form b*phi(x - c*a) is exact.
  a levels: [0 0.25 0.5 0.75 1] ; b levels: [1 0.75 0.5 0.25 0]

FIGURES PRODUCED (in this folder)
  Fig_FI_curve.png
  Fig_FI_curve.svg
  Fig_FI_curve.fig
