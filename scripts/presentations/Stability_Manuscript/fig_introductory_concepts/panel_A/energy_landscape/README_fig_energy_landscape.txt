Stability_Manuscript figure 1 panel A (top): effective potential
================================================================

Generated: 26-Aug-2026 01:10:44
By:        fig_energy_landscape.m

WHAT IT SHOWS
  The effective potential U(x) = x^2/2 - gamma*log(cosh(x)) for the 1-D
  reduction of the network, at three gains. The curvature at the origin is 1 -
  gamma, so x = 0 is a stable minimum below gamma = 1 and an unstable maximum
  above it. Arrows show the flow direction, always downhill; the ball sits at
  x = 0 in every panel, so the story is the landscape changing underneath it.

HOW IT WAS MADE
  ANALYTIC -- no model is constructed. This is the one part of panel A that
  takes no preset, because it is the 1-D reduction rather than the network.

PARAMETERS AS RUN
  gammas  [0.9 1.6 2.5]

FIGURES PRODUCED (in this folder)
  panelA_energy_landscape.png
  panelA_energy_landscape.svg
  panelA_energy_landscape.fig

READING THIS FIGURE
  THE GAMMAS MUST MATCH fig_introductory_concepts, which draws the traces and
  eigenspectra directly below this panel at the same gains. They are an
  argument in both functions rather than a literal in each, so the master
  script can pass one value to both and they cannot drift apart.

