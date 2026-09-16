# Main group 1

Run: `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\mu7revised_fast`

- Source: `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7revised_fast\fig_introductory_concepts\Fig_Intro_Concepts.fig`

For the current mu7 figure-only configuration, the native introductory source is explicitly figs/sfaEI_fast/fig_introductory_concepts (clean source commit4406b56). It uses the same sompolinsky_pairs preset, gammas[0.9,1.6,2.5], seed0, 15 traces and [0,60]s as the mu7 intro. The intervening fig_introductory_concepts code change only combined the two rows; simulation code is unchanged.

STD endpoint: tau_rec=[2 4] s, tau_rel=[0.25 0.5] s, raw reference rate r_ref=0.25. Each b_m,ss=1/(1+r_ref*tau_rec,m/tau_rel,m)=[0.333333333333333 0.333333333333333]; total product B_min=0.111111111. Five displayed frozen total factors=[1 0.777777777777778 0.555555555555556 0.333333333333333 0.111111111111111]. Intermediate factors are illustrative, not simulated steady-state points.

The sigmoid panels are conceptual logistic curves (center0.4; enlarged total SFA shift0.6), not fits to the reference piecewise activation. STD uses a frozen product B multiplying the entire curve, not the self-consistent rate-dependent steady-state input-output relation. The two-timescale product must not be identified with a single depression variable.

Each SFA curve and shifted disk shares the exact black-to-orange color; each STD curve and shrinking disk shares the exact black-through-dark-blue-to-teal color. Disk shifts/radii are effective-connectivity intuition, not exact transformations of the full active Jacobian.

Top panels are copied native saved graphics; all scientific XData/YData, neuron selections and gains are unchanged. The10-unit scale bar and its text are omitted because the example time is arbitrary up to rescaling. Only position, row labels, 14-point fonts and trace y-axis width1.0 are changed. No simulation or Lyapunov calculation is run.

