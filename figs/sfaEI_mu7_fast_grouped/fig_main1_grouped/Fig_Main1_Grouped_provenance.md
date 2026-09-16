# Main group 1

Run: `/Users/tom/Desktop/local_code/FractionalReservoir/data/sfaEI_mu7_fast`

- Source: `/Users/tom/Desktop/local_code/FractionalReservoir/figs/sfaEI_fast/fig_introductory_concepts/eigenspectra/panelA_eigenspectrum.fig`
- Source: `/Users/tom/Desktop/local_code/FractionalReservoir/figs/sfaEI_fast/fig_introductory_concepts/statetraces/panelA_bottom_traces.fig`

For the current mu7 figure-only configuration, the native introductory source is explicitly figs/sfaEI_fast/fig_introductory_concepts (clean source commit4406b56). It uses the same sompolinsky_pairs preset, gammas[0.9,1.6,2.5], seed0, 15 traces and [0,60]s as the mu7 intro. The intervening fig_introductory_concepts code change only combined the two rows; simulation code is unchanged.

STD endpoint: tau_rec=[2 4] s, tau_rel=[0.25 0.5] s, raw reference rate r_ref=0.25. Each b_m,ss=1/(1+r_ref*tau_rec,m/tau_rel,m)=[0.333333333333333 0.333333333333333]; total product B_min=0.111111111. Five displayed frozen total factors=[1 0.777777777777778 0.555555555555556 0.333333333333333 0.111111111111111]. Intermediate factors are illustrative, not simulated steady-state points.

The sigmoid panels are conceptual logistic curves (center0.4; enlarged total SFA shift0.6), not fits to the reference piecewise activation. STD uses a frozen product B multiplying the entire curve, not the self-consistent rate-dependent steady-state input-output relation. The two-timescale product must not be identified with a single depression variable.

Each SFA curve and shifted disk shares the exact black-to-orange color; each STD curve and shrinking disk shares the exact black-through-dark-blue-to-teal color. Disk shifts/radii are effective-connectivity intuition, not exact transformations of the full active Jacobian.

Top panels are copied native saved graphics; all scientific XData/YData, neuron selections and gains are unchanged. The10-unit scale bar and its text are omitted because the example time is arbitrary up to rescaling. Only position, row labels, 14-point fonts and trace y-axis width1.0 are changed. No simulation or Lyapunov calculation is run.


Final compact-canvas update (fb6a3fb): Fixed canvas910×665 (70% of1300×950 creation baseline), fonts14 unchanged; scientific native traces and paired analytic curves retained.
