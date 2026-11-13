---
title: "Lab meeting figures: blocks A, B and C"
date: 2026-11-13
abstract: |
  Every figure of blocks A, B and C of the final-draft work so far. Sections A and C show the four conditions (no adaptation, single-timescale, chained MTS, product MTS) wherever a figure was drawn for them; the block B pilots ran the conditions each question needed, and every entry with fewer than four says which.

  **Equations of the model:** [Equations_chained_MTS.pdf](file:///C:/Users/m218089/Desktop/github_repos/FractionalReservoir/docs/EquationsParametersDocs/chained_MTS/Equations_chained_MTS.pdf) (the chained multiple-timescale model, its conditions, and in section 10 the Jacobian in block form), at `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\docs\EquationsParametersDocs\chained_MTS\Equations_chained_MTS.pdf`.
header-includes:
  - |
    ```{=latex}
    \AtBeginDocument{\let\oldtextbackslash\textbackslash
    \renewcommand{\textbackslash}{\oldtextbackslash\allowbreak}\sloppy}
    ```
---

# A. Chained vesicle-pool model and the essential run

## Essential run

### A1. Example dynamics (main group 2)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_main2_grouped/Fig_Main2_Grouped.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_main2_grouped\Fig_Main2_Grouped.png`

**Question:** How do the representative network trajectories (states, local rate and finite-time $\lambda_1$) differ across the four adaptation conditions after a uniform input step?

- Columns are no adaptation, 1TS, chained MTS and product MTS, each with the same 12 E and 12 I neurons; row G shows the local rate (low-pass filtered at 2 Hz for display only) and the finite-time $\lambda_1$.

### A2. Eigenvalues and switching (main group 3)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_main3_grouped/Fig_Main3_Grouped.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_main3_grouped\Fig_Main3_Grouped.png`

**Question:** How do the Jacobian eigenvalue occupancy, the local leading Lyapunov rates and the expansive transients differ by condition?

- Row C shows the local and accumulated leading rates (K = 50, accumulated over 1 to 40 s); row D the positive local QR growth in bits/s, which is not a measure of KS entropy.

### A3. Sensitivity of $\lambda_1$ and $\lambda_1$ against rate (main group 4)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_main4_grouped/Fig_Main4_Grouped.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_main4_grouped\Fig_Main4_Grouped.png`

**Question:** How does the median $\lambda_1$ respond to the E:I ratio, network size and the four weight blocks, and how does it relate to the firing rate in each condition?

- No adaptation swings far beyond plus or minus 1/s across every sweep; both MTS conditions stay near zero for most parameter values (chained MTS rises to about +1/s at high E to I weights and about +0.5/s at moderate E to E and I to E weights).
- In panel B, chained and product MTS sit at $\lambda_1$ near 0 across the whole rate range; no adaptation falls to about -5/s as it saturates.

### A4. $\lambda_1$ sensitivity distributions (core sweeps)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_analysis_allStd/Fig_Sensitivity_LLE_core.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_analysis_allStd\Fig_Sensitivity_LLE_core.png`

**Question:** What is the full distribution of $\lambda_1$ across the core 1-D sweeps (E:I ratio, network size, overall gain) for each condition?

- Distributions per condition and sweep level; the essential run has 4 levels x 2 reps per sweep.

### A5. $\lambda_1$ sensitivity distributions (weight blocks)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_analysis_allStd/Fig_Sensitivity_LLE_mu.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_analysis_allStd\Fig_Sensitivity_LLE_mu.png`

**Question:** What is the distribution of $\lambda_1$ across the four weight-block sweeps for each condition?

- The four mu sweeps span -75 % to +75 % of the preset value, at 4 levels x 2 reps.

### A6. Mean rate sensitivity (core sweeps)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_analysis_allStd/Fig_Sensitivity_mean_rate_core.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_analysis_allStd\Fig_Sensitivity_mean_rate_core.png`

**Question:** How does the mean firing rate depend on the core sweep parameters in each condition?

- Distributions of the population mean rate per sweep level for all four conditions (4 levels x 2 reps).

### A7. Mean rate sensitivity (weight blocks)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_analysis_allStd/Fig_Sensitivity_mean_rate_mu.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_analysis_allStd\Fig_Sensitivity_mean_rate_mu.png`

**Question:** How does the mean firing rate depend on the four weight blocks in each condition?

- Distributions of the population mean rate across the mu sweeps for all four conditions.

### A8. KS entropy sensitivity (core sweeps)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_analysis_allStd/Fig_Sensitivity_hKS_core.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_analysis_allStd\Fig_Sensitivity_hKS_core.png`

**Question:** How does the Kolmogorov-Sinai entropy h_KS vary across the core sweeps in each condition?

- h_KS distributions (bit/s) per sweep level for all four conditions.

### A9. KS entropy sensitivity (weight blocks)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_analysis_allStd/Fig_Sensitivity_hKS_mu.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_analysis_allStd\Fig_Sensitivity_hKS_mu.png`

**Question:** How does h_KS vary across the four weight-block sweeps in each condition?

- h_KS distributions per mu sweep level for all four conditions.

### A10. Kaplan-Yorke dimension sensitivity (core sweeps)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_analysis_allStd/Fig_Sensitivity_DKY_core.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_analysis_allStd\Fig_Sensitivity_DKY_core.png`

**Question:** How does the attractor dimension D_KY vary across the core sweeps in each condition?

- D_KY distributions per sweep level for all four conditions.

### A11. Kaplan-Yorke dimension sensitivity (weight blocks)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_analysis_allStd/Fig_Sensitivity_DKY_mu.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_analysis_allStd\Fig_Sensitivity_DKY_mu.png`

**Question:** How does D_KY vary across the four weight-block sweeps in each condition?

- D_KY distributions per mu sweep level for all four conditions.

### A12. $\lambda_1$ sensitivity medians

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_medians/Fig_Sensitivity_LLE_medians.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_medians\Fig_Sensitivity_LLE_medians.png`

**Question:** How does the median $\lambda_1$ (with IQR) move with each swept parameter in each condition?

- Same data as main group 4, panel A: both MTS conditions stay close to zero, while no adaptation varies widely.

### A13. Mean rate sensitivity medians

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_medians/Fig_Sensitivity_mean_rate_medians.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_medians\Fig_Sensitivity_mean_rate_medians.png`

**Question:** How does the median mean firing rate move with each swept parameter in each condition?

- Medians and IQR of the mean rate per sweep for all four conditions.

### A14. KS entropy sensitivity medians

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_medians/Fig_Sensitivity_hKS_medians.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_medians\Fig_Sensitivity_hKS_medians.png`

**Question:** How does the median h_KS move with each swept parameter in each condition?

- Medians and IQR of h_KS per sweep for all four conditions.

### A15. Kaplan-Yorke dimension sensitivity medians

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_medians/Fig_Sensitivity_DKY_medians.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_medians\Fig_Sensitivity_DKY_medians.png`

**Question:** How does the median D_KY move with each swept parameter in each condition?

- Medians and IQR of D_KY per sweep for all four conditions.

### A16. Joint parameter space: $\lambda_1$

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_param_space_allStd/Fig_ParamSpace_LLE.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_param_space_allStd\Fig_ParamSpace_LLE.png`

**Question:** Over a joint random sample of parameter space, how is $\lambda_1$ distributed in each condition?

- No adaptation splits into two modes, about 70 % at the negative clip and about 30 % at the positive clip; 1TS is spread out, with a mode near -0.4/s.
- Chained MTS has about 70 % and product MTS about 75 % of points in the bin just below zero.

### A17. Joint parameter space: mean rate

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_param_space_allStd/Fig_ParamSpace_mean_rate.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_param_space_allStd\Fig_ParamSpace_mean_rate.png`

**Question:** Over the joint parameter sample, how is the mean firing rate distributed in each condition?

- Histograms of the mean rate over the joint sample, one per condition.

### A18. Joint parameter space: KS entropy

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_param_space_allStd/Fig_ParamSpace_hKS.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_param_space_allStd\Fig_ParamSpace_hKS.png`

**Question:** Over the joint parameter sample, how is h_KS distributed in each condition?

- Histograms of h_KS over the joint sample, one per condition.

### A19. Joint parameter space: Kaplan-Yorke dimension

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_param_space_allStd/Fig_ParamSpace_DKY.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_param_space_allStd\Fig_ParamSpace_DKY.png`

**Question:** Over the joint parameter sample, how is D_KY distributed in each condition?

- Histograms of D_KY over the joint sample, one per condition.

### A20. E:I weight balance: $\lambda_1$

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_EI_weights_param_space/Fig_EI_Weights_ParamSpace_LLE.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_EI_weights_param_space\Fig_EI_Weights_ParamSpace_LLE.png`

**Question:** How does $\lambda_1$ depend on the realized E:I weight balance in each condition?

- $\lambda_1$ against the E:I weight balance across the joint sample, one panel per condition.

### A21. E:I weight balance: mean rate

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_EI_weights_param_space/Fig_EI_Weights_ParamSpace_mean_rate.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_EI_weights_param_space\Fig_EI_Weights_ParamSpace_mean_rate.png`

**Question:** How does the mean rate depend on the E:I weight balance in each condition?

- Mean rate against the E:I weight balance, one panel per condition.

### A22. E:I weight balance: KS entropy

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_EI_weights_param_space/Fig_EI_Weights_ParamSpace_hKS.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_EI_weights_param_space\Fig_EI_Weights_ParamSpace_hKS.png`

**Question:** How does h_KS depend on the E:I weight balance in each condition?

- h_KS against the E:I weight balance, one panel per condition.

### A23. E:I weight balance: Kaplan-Yorke dimension

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_EI_weights_param_space/Fig_EI_Weights_ParamSpace_DKY.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_EI_weights_param_space\Fig_EI_Weights_ParamSpace_DKY.png`

**Question:** How does D_KY depend on the E:I weight balance in each condition?

- D_KY against the E:I weight balance, one panel per condition.

### A24. $\lambda_1$ against mean firing rate

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_lle_vs_rate/Fig_LLE_vs_Rate.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_lle_vs_rate\Fig_LLE_vs_Rate.png`

**Question:** Is $\lambda_1$ tied to the firing rate, and does adaptation remove that dependence?

- Spearman rho [95 % CI], n = 144 per condition: no adaptation -0.40 [-0.56, -0.23], 1TS -0.43 [-0.57, -0.29], chained MTS -0.25 [-0.37, -0.12], product MTS +0.08 [-0.08, +0.21].
- Median $\lambda_1$ in the mid-rate band: +2.24 (none), -0.35 (1TS), -0.10 (chained), -0.10 (product); in saturated networks: -9.59 (n = 25), -1.93 (n = 10), -0.10 (n = 9), -0.09 (n = 2).

### A25. Jacobian eigenvalue heatmap

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_eig_heatmap/fig_eig_heatmap.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_eig_heatmap\fig_eig_heatmap.png`

**Question:** Where do the instantaneous Jacobian eigenvalues sit along the trajectory in each condition, relative to the finite-time $\lambda_1$?

- Finite-time $\lambda_1$ over the last 20 s: no adaptation +3.092, 1TS +0.233, chained MTS -0.105, product MTS -0.098.
- All three adapting conditions pile up eigenvalue density just left of Re = 0; the product MTS bulk is the most compact, reaching to about plus or minus 10 on the imaginary axis.

### A26. Memory capacity

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_mc_4cond/fig_memory_capacity/Fig_Memory_Capacity.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_mc_4cond\fig_memory_capacity\Fig_Memory_Capacity.png`

**Question:** How much past input can a linear readout of the synaptic output reconstruct in each condition, and for how long?

- Total memory capacity, mean [95 % CI]: no adaptation 0.105, 1TS 2.89 [1.74, 4.17], chained MTS 9.61 [7.06, 13.29], product MTS 14.04 [13.55, 14.31]; memory horizon 0 s, 2.0 s, 5.7 s and 8.7 s in the same order.
- Only 3 paired trials, so every pairwise sign-flip test has p = 0.25 (the smallest possible); chained against product MTS: difference -4.43, Cohen's d_z -1.39.

## A explorations

### A27. Single- against multiple-timescale depression (FD-014)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/prelim_fd014/4cond/Fig_STD_1TS_vs_MTS.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\prelim_fd014\4cond\Fig_STD_1TS_vs_MTS.png`

**Question:** On one E to E route, how does chained two-pool depression differ from one pool and from the product of two pools, under a sinusoidal rate, rate steps and constant rates?

- The four-condition version; it shows three traces because the no-adaptation condition has no depression. First-cycle peak output relative to 1TS: chained MTS 0.986, product MTS 0.725; by the last cycle 0.736 and 0.511, as the chain's reserve (b_2) runs down to about 0.5.
- Sustained output r F has a ceiling for one pool (0.125) and for the chain (0.056); the product has none and peaks at 0.031 (r = 0.125) before falling.

### A28. Random input steps with the top 30 Lyapunov exponents

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/examples_show/mu7chain_steps_top30_seed78/fig_example_timeseries.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\examples_show\mu7chain_steps_top30_seed78\fig_example_timeseries.png`

**Question:** On the mu7chain network (seed [7 8]), how do the four conditions respond to a new random input step every 5 s from 10 to 30 s, and how many of the top 30 Lyapunov exponents are positive?

- $\lambda_1$ (accumulated from 5 to 40 s): no adaptation +2.551 with 5 of the top 30 positive; 1TS -0.099, chained MTS -0.101, product MTS -0.093, none positive.
- Without adaptation the dendritic states span about plus or minus 10 and the rates switch between 0 and 1; with adaptation the states stay near 0 and the synaptic output below about 0.3.

### A29. The input of the steps example

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/examples_show/mu7chain_steps_top30_seed78/fig_example_stimulus.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\examples_show\mu7chain_steps_top30_seed78\fig_example_stimulus.png`

**Question:** What input does the steps example (A28) deliver, identically in all four conditions?

- Four 5-s steps from 10 to 30 s, each to 15 % of each cell type drawn anew, amplitude 0.5 randn of either sign: 66, 87, 74 and 89 neurons per step, 253 of 500 in at least one step (129 E, 124 I), largest amplitude 1.508.

### A30. Eigenvalues of the connectivity of the steps example

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/examples_show/mu7chain_steps_top30_seed78/fig_example_W_spectrum.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\examples_show\mu7chain_steps_top30_seed78\fig_example_W_spectrum.png`

**Question:** What is the spectrum of the simulated W of the mu7chain network at seed [7 8], shared by all four conditions?

- Spectral abscissa 4.880, spectral radius 5.276, bulk radius R = 4.799; the outlier eigenvalues are 0, so the bulk disc alone sets the spectrum, with a few eigenvalues just past its edge.

# B. Stimulation pilots

## B.1 P0a: the run13 network with chained MTS

### B1. run13 under three adaptation forms

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/run13_chain_mts/run13_chain_mts.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\run13_chain_mts\run13_chain_mts.png`

**Question:** Does the run13 IED-like network still produce brief, localized discharges with chained MTS adaptation, and does the chained form behave like the product form? (Three conditions: original, that is run13's single-timescale adaptation, chained MTS and product MTS; no no-adaptation.)

- Chained MTS keeps the discharges but makes them rarer: 12 events in 100 s against 30 in the original, still local (median 8 % of neurons, one module each, none global). Product MTS nearly silences them: 1 event.
- Mean E rate 0.125 (original), 0.091 (chained), 0.072 (product); the product depresses excitatory transmission most (mean depression factor 0.557 against 0.585 chained), the likely reason it loses its events.

## B.2 P0b: the bursting network with chained MTS

### B2. The bursting network under three adaptation forms

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/bursting_chain_mts/bursting_chain_mts.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\bursting_chain_mts\bursting_chain_mts.png`

**Question:** Does the small, strongly coupled bursting network still burst with chained MTS, and does the chained form behave like the product of two pools? (Three conditions: original, three SFA timescales on E and one depression pool, chained MTS and product MTS.)

- 33, 34 and 0 bursts in 100 s (original, chained, product). Chained MTS keeps the rhythm (one burst every 2.95 s against 3.0 s) with smaller bursts (70 % of neurons active at the peak against 82 %).
- Product MTS settles into a small oscillation (mean E rate peaking at 0.10): its second pool holds the synapses depressed (factor 0.67 to 0.76) and the population never ignites.

## B.3 P1a: stimulation of the run13 network

### B3. Example pair, stimulation off and on

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p1a_stim_run13/p1a_example_pairs.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p1a_stim_run13\p1a_example_pairs.png`

**Question:** For network seed 42, does a sustained uniform input of 0.20 remove run13's IED-like events under chained MTS as under the original, and does it recruit adaptation? (Original and chained MTS, each off and on.)

- The events drop in both ON columns while the SFA term steps up and the depression factor steps down; under chained MTS the SFA term keeps rising for about 30 s through its slow rungs.

### B4. Suppression and recruitment over 8 seeds

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p1a_stim_run13/p1a_summary.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p1a_stim_run13\p1a_summary.png`

**Question:** Over 8 seeds, is the suppression as large under chained MTS as under the original, and how much adaptation does it recruit? (Original and chained MTS.)

- At 0.20 the event rate fell from 15.3 to 2.55 per min (chained) and from 20.3 to 3.55 (original), lower in 8 of 8 seeds in each (sign-flip p = 0.0078); relative reduction 0.83 in both, 0.63 against 0.54 at 0.10 and 0.93 against 0.92 at 0.30.
- E SFA term +61 % (chained) and +51 % (original), depression factor -20 % and -14 %, mean E rate +62 % and +51 %: every seed moved the same way.

## B.4 P1b: stimulation of the bursting network

### B5. Original adaptation under five stimuli

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p1b_stim_bursting/p1b_original.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p1b_stim_bursting\p1b_original.png`

**Question:** Under the original adaptation, does each tonic stimulus stop the spontaneous bursts, and what replaces them? (Original only.)

- Every all-neuron step (0.1 to 0.4) turns the bursts into a flat rate (reference: mean E rate 0.081 to 0.297, adaptation term 0.041 to 0.147, depression factor 0.82 to 0.57).
- 0.4 to 20 % of the E neurons leaves a small fast oscillation (0.88 Hz, mean E rate swinging 0.14 to 0.24).

### B6. Chained MTS under five stimuli

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p1b_stim_bursting/p1b_chain_mts.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p1b_stim_bursting\p1b_chain_mts.png`

**Question:** Under chained MTS, does each tonic stimulus stop the bursts, and is adaptation recruited as they stop? (Chained MTS only.)

- Every case, partial stimulation included, ends the bursts at a flat rate (reference: mean E rate 0.060 to 0.200, adaptation term 0.030 to 0.099, depression factor 0.75 to 0.47); the grey stripes at the reference are noise-driven crossings of a rate sitting at the 0.200 threshold, not bursts.
- Graded with amplitude: mean E rate 0.200, 0.133 and 0.097 at 0.4, 0.2 and 0.1.

### B7. Before against during, every stimulus

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p1b_stim_bursting/p1b_summary.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p1b_stim_bursting\p1b_summary.png`

**Question:** Per stimulus case, how do bursts, rates, adaptation and depression change from before to during stimulation? (Original and chained MTS; product MTS at the reference stimulus only.)

- The largest population mean E rate and its peak-to-peak range collapse (chained: 0.39 to 0.21, range 0.39 to 0.01) while the mean rate and adaptation rise: suppression is steady elevated firing, not quiet, and no neuron saturates.
- Product MTS, which had no bursts, fires one onset transient (largest mean E rate 0.57) and settles at a steady 0.170.

## B.5 P2: adaptation frozen at its pre-stimulation value

### B8. The paper's network, adaptation intact and frozen

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p2_frozen_state/p2_frozen_mu7chain.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p2_frozen_state\p2_frozen_mu7chain.png`

**Question:** Can the paper's network be run with adaptation held at its value at stimulus onset, and what does the stimulus do then? (Chained MTS only; intact and frozen, each with and without the stimulus.)

- The frozen arm is exact by script: with nothing frozen the continuation reproduces the model's run with a largest difference of 0, and the frozen variables drift by 0.
- E / I rates after onset: 0.293 / 0.253 intact, 0.459 / 0.438 frozen (no stimulus 0.202 / 0.211): free adaptation absorbs about two thirds of the rate increase.

### B9. The bursting network, adaptation intact and frozen

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p2_frozen_state/p2_frozen_bursting.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p2_frozen_state\p2_frozen_bursting.png`

**Question:** Does the stimulus still stop the bursts when adaptation cannot respond? (Original condition only; intact and frozen, each with and without the stimulus.)

- Intact, the stimulus stops the bursts and leaves steady firing (E / I 0.359 / 0.223); frozen, the network runs away to saturation (0.906 / 0.984).
- Freezing alone is not neutral: frozen between two bursts, the unstimulated network falls silent (0.007 / 0.016), so a frozen arm needs its own no-stimulus control.

## B.6 P3: adaptation removed, and a fair ablated baseline

### B10. Bursting network, five ablations

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p3_ablation/bursting_ablations.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p3_ablation\bursting_ablations.png`

**Question:** What does the bursting network do, unstimulated, when adaptation is removed in each of five ways? (Chained MTS intact and ablated.)

- Bursts survive only when the slow parts alone are removed, and then twice as often (70 against 34 per 100 s); removing SFA gives steady firing (mean E rate 0.10).
- Removing depression, or both, saturates the network (mean E rate 0.83 to 0.92, 74 to 88 % of neurons pinned); the mean-matched constants saturate it too (0.77).

### B11. run13, five ablations

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p3_ablation/run13_ablations.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p3_ablation\run13_ablations.png`

**Question:** What does run13 do, unstimulated, when adaptation is removed in each of five ways? (Chained MTS intact and ablated.)

- Events survive SFA removal (11 against 12 per 100 s, similar size and duration): the only ablation with comparable aberrant activity, and it tests SFA alone.
- Removing depression saturates the network (mean E rate 1.00, 99 to 100 % pinned), and so do the mean-matched constants (0.97).

### B12. Bursting network, variants

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p3_ablation/bursting_variants.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p3_ablation\bursting_variants.png`

**Question:** Does a different matched constant, or keeping only the fast SFA share, give the bursting network a fair ablated baseline? (Chained MTS variants, and the original condition intact and ablated.)

- No. The rate-weighted depression constant (0.53) silences the network (mean E rate 0.003); constants between 0.53 and 0.75 are silent up to 0.71 and saturate at 0.75, with no bursts anywhere.

### B13. run13, variants

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p3_ablation/run13_variants.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p3_ablation\run13_variants.png`

**Question:** Does a different matched constant, or keeping only the fast SFA share, give run13 a fair ablated baseline? (Chained MTS variants, and the original condition intact and ablated.)

- No. The rate-weighted constant (0.33) gives steady firing at 0.14 with no events; constants per neuron saturate two modules (mean E rate 0.55).
- Removing the slow parts makes run13 exactly its original single-timescale condition, with more events (30 against 12).

### B14. Every ablation, both networks

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p3_ablation/summary.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p3_ablation\summary.png`

**Question:** Across every ablation and both networks, which removal leaves a baseline with comparable aberrant activity? (Chained MTS and original.)

- None that removes depression or replaces adaptation by constants: the aberrant events depend on adaptation changing in time.
- Clamping adaptation from stimulus onset (P2) is the only candidate for a fair "adaptation ablated" arm.

## B.7 P4: one measure of aberrant activity

### B15. Excess co-activation against the detectors, unstimulated

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p4_aberrant_measure/p4_recommended_unstimulated.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p4_aberrant_measure\p4_recommended_unstimulated.png`

**Question:** Does excess co-activation A, computed from the rates alone, agree with each network's own detector? (Original, chained MTS and product MTS on both networks.)

- A, in thousandths: bursting network 47.5, 37.8, 3.1 (detector 33, 34, 0 bursts); run13 1.34, 0.54, 0.10 (detector 30, 12, 1 events). Same order on both networks; every detector event coincides with co-activation above chance.
- The non-bursting product network reads 3.1 and counts as aberrant: a regular low-amplitude oscillation, since A measures co-activation, not amplitude.

### B16. The four arms of P2 on the bursting network

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p4_aberrant_measure/p4_recommended_bursting_arms.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p4_aberrant_measure\p4_recommended_bursting_arms.png`

**Question:** Is the measure fooled by a steady, a saturated or a silent network? (Original condition; P2's four arms.)

- Only the unstimulated intact arm has co-activation. The steady, the saturated (mean rate 0.93, 86 % pinned) and the silent arm (0.010) all read A = 0 and are told apart by the mean rate and the pinned fraction.

### B17. Stimulation, and the paper's network

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p4_aberrant_measure/p4_recommended_stim_and_paper.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p4_aberrant_measure\p4_recommended_stim_and_paper.png`

**Question:** Does A fall under stimulation, and what does it read on the paper's own network? (run13 original and chained MTS, with and without input; the paper's network with no adaptation and chained MTS.)

- Under an input of 0.20, A falls from 1.32 to 0.13 (run13 original) and 0.54 to 0.02 (chained MTS).
- The paper's chained-MTS network reads 1.85 (about 15 episodes per min above chance), so by the proposed threshold of 0.25 it is hyperexcitable at baseline; its no-adaptation condition is classed saturated (41 % of neurons pinned).

### B18. The ten candidate measures

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p4_aberrant_measure/p4_candidates_summary.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p4_aberrant_measure\p4_candidates_summary.png`

**Question:** How do the ten candidate measures compare with the detector counts and under stimulation? (Original, chained MTS and product MTS on each network.)

- Every other candidate fails a criterion: measures from the mean rate miss run13's events (about 8 % of neurons), synchrony reads 0.49 for the non-bursting product network and 0.61 for a stimulated arm whose bursts stopped, and co-activation against a fixed level is too steep on run13 (10 : 1 : 0).

## B.8 P5: DC against pulses of equal mean

### B19. Chained MTS at mean input 0.1

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p5_dc_vs_pulses/p5_chain_mts_mean0p1.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p5_dc_vs_pulses\p5_chain_mts_mean0p1.png`

**Question:** At mean input 0.1, do DC and pulses at 2, 10 and 50 Hz act alike on the bursts and on the adaptation they recruit? (Chained MTS only.)

- Every form stops the network's own bursts (0.33 per s before). At 2 Hz each pulse evokes exactly one locked population event (vector strength 1.00) and depresses the synapses more (factor 0.565 against 0.631 under DC); 10 and 50 Hz act as DC.

### B20. Chained MTS at mean input 0.2

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p5_dc_vs_pulses/p5_chain_mts_mean0p2.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p5_dc_vs_pulses\p5_chain_mts_mean0p2.png`

**Question:** Does the same hold at mean input 0.2? (Chained MTS only.)

- 2 Hz against DC: depression factor 0.459 against 0.560, adaptation term 0.082 against 0.066, mean E rate 0.166 against 0.133; the evoked events (mean E rate peak 0.62) are larger than the spontaneous bursts (0.39).
- The grey band under 10 Hz is the detector merging peaks that graze the threshold (0.206), not a burst.

### B21. The original condition as a reference

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p5_dc_vs_pulses/p5_original_reference.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p5_dc_vs_pulses\p5_original_reference.png`

**Question:** Does the original adaptation respond to DC and 2 Hz pulses the same way? (Original only, DC and 2 Hz at both levels.)

- Yes: both stop its bursts, 2 Hz evokes one event per pulse and recruits more (mean 0.2: depression factor 0.572 against 0.629 under DC, adaptation term 0.121 against 0.107).

### B22. Summary over forms and levels

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p5_dc_vs_pulses/p5_summary.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p5_dc_vs_pulses\p5_summary.png`

**Question:** During stimulation, how do bursts, rate, depression and adaptation depend on the form of the input? (Chained MTS and original.)

- Only 2 Hz differs from DC: at 10 and 50 Hz mean rate, depression and adaptation are within 2.5 % and 0.1 % of DC. The pulse-triggered swing of the mean E rate is 0.011 and 0.024 at 50 Hz against 0.36 and 0.62 at 2 Hz (mean 0.1 and 0.2).

## B.9 P6: hyperexcitability in the paper's parameter space

### B23. Chaos and rate across the 32 sampled networks

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p6_hyperexcitability/p6_lle_vs_rate.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p6_hyperexcitability\p6_lle_vs_rate.png`

**Question:** Of the 32 networks sampled across the paper's parameter space, how many are hyperexcitable at baseline (chaotic, high rate, or both), and does adaptation move any network across those thresholds?

- In all four conditions about 10 of 32 networks are chaotic, 11 to 12 are high-rate (mean rate above 0.40), and only 2 of 32 are both (6 %, 95 % interval 2 to 20 %). Adaptation moves almost no network across a threshold: all 10 chaotic networks stay chaotic under chained MTS.
- The stored scalars miss aberrant activity: network 761 reads as stable and low-rate under chained MTS but bursts about every 4 s, so a parameter-space experiment needs an event measure and a large sample (about 145 to 1600 networks for 50 hyperexcitable ones).

### B24. Where the classes lie in connectivity space

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p6_hyperexcitability/p6_where_in_parameter_space.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p6_hyperexcitability\p6_where_in_parameter_space.png`

**Question:** Which connectivity parameters decide whether a sampled network is chaotic or high-rate? (No adaptation and chained MTS.)

- High rate follows the excitatory fraction (all 8 networks at 80 % E, 3 or 4 of 11 at 50 %, none of 13 at 20 %) and the E:I weight balance (above about 0.5).
- Chaos follows the gain on all weights (0 of 8 at 0.25, 4 of 12 at 1, 6 of 12 at 1.75) and prefers large, inhibition-dominated networks.

### B25. Six rebuilt networks

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p6_hyperexcitability/p6_timeseries_examples.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p6_hyperexcitability\p6_timeseries_examples.png`

**Question:** What does the activity of networks in each class actually look like? (No adaptation and chained MTS.)

- The two chaotic, high-rate networks have every excitatory neuron pinned at its maximum, with the chaos in a sparsely active inhibitory population.
- Network 761, stable and low-rate by its stored scalars under chained MTS, is a regular population burster (one full-amplitude burst about every 4 s).

# C. Memory capacity across gain, bias and self-feedback

## Gain by bias

### C1. Memory capacity and $\lambda_1$ across gain and bias

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/block_C/main/four_conditions/Fig_MC_gain_bias_maps.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\block_C\main\four_conditions\Fig_MC_gain_bias_maps.png`

**Question:** Is the memory-capacity result a property of adaptation, or of one well-tuned network?

- In every condition, memory capacity is highest along a ridge just on the stable side of $\lambda_1$ = 0. Without adaptation the edge stays at gain 0.33 to 0.5 at every bias; with adaptation it moves to higher gain, and further the lower the bias is.
- At the preset bias, memory capacity stays at or above 7.5 over 0.55 octaves of gain without adaptation, against 2.0 to 2.4 octaves with adaptation. Peaks: 10.2 (none), 15.0 (1TS), 15.0 (chained MTS), 12.8 (product MTS).

### C2. Memory capacity and $\lambda_1$ across gain and bias, linear scales

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_C/linear_scales/four_conditions/Fig_MC_gain_bias_maps_linear.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_C\linear_scales\four_conditions\Fig_MC_gain_bias_maps_linear.png`

**Question:** Do the gain-by-bias maps of C1 read the same when memory capacity and $\lambda_1$ are drawn on linear color scales?

- The ridge of high memory capacity hugs the stable side of the edge in every condition. The preset network (star) sits on the edge for 1TS and chained MTS, but inside the stable region for product MTS.
- The spread over the two seeds is concentrated at cells on the edge. Product MTS's two seeds agree least (r = 0.61).

### C3. Memory capacity across gain at three E:I balances

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/block_C/main/four_conditions/Fig_MC_gain_curves.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\block_C\main\four_conditions\Fig_MC_gain_curves.png`

**Question:** How do memory capacity and $\lambda_1$ change with gain at the preset bias, and does that change when excitation or inhibition dominates?

- When excitation dominates (3 : 2), no adaptation saturates (memory capacity at most 3.1, and 0 from gain 0.7), while every adapting condition keeps memory over the whole gain range; there chained MTS (6.5 to 9.7) is ahead of 1TS (6 to 7.8). One seed only.
- When inhibition dominates (2 : 3), the three adapting conditions are level (peaks 12.2 to 12.9 at gain 0.7). The order of the conditions changes from one network to the next.

### C4. Robustness of memory capacity over the map

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/block_C/main/four_conditions/Fig_MC_robustness.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\block_C\main\four_conditions\Fig_MC_robustness.png`

**Question:** Over all gains, biases and seeds, what share of the networks keeps useful memory in each condition?

- Without adaptation, fewer than 40 % of the networks have any memory (memory capacity above about 0.5); with adaptation, 65 to 85 % do.
- Above a memory capacity of about 8 the three adapting conditions are close together. Chained MTS against 1TS: higher in both seeds at 40 % of grid points, lower in both at 22 %, median difference +0.03.

### C5. Memory capacity against $\lambda_1$, every run

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/block_C/main/four_conditions/Fig_MC_vs_lambda1.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\block_C\main\four_conditions\Fig_MC_vs_lambda1.png`

**Question:** Does memory capacity track the distance to the edge of chaos ($\lambda_1$) in every condition?

- High memory capacity sits just below $\lambda_1$ = 0 in all four conditions; past zero it falls to about 0.1 to 1. Pearson r over all runs: -0.68 (none), -0.65 (1TS), -0.70 (chained), -0.61 (product).
- With chained and product MTS, $\lambda_1$ is pinned near -0.1/s across the whole stable region (by the slowest adaptation timescale), so it says little about how close a network is to the edge: on the stable side r = 0.87 without adaptation (n = 57) against -0.20 for product MTS (n = 120). (Also drawn by the exploration script, as `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_C\mc_vs_lambda1\Fig_MC_vs_lambda1_scatter.png`.)

## Gain-by-bias explorations

### C6. Memory capacity against $\lambda_1$, marker size by bias

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_C/bias_rate_lambda1/Fig_MC_vs_lambda1_size_bias.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_C\bias_rate_lambda1\Fig_MC_vs_lambda1_size_bias.png`

**Question:** Where do the different biases fall along the relation between memory capacity and $\lambda_1$?

- Points at all biases fall on the same curve. With MTS, networks from every bias collapse onto $\lambda_1$ near -0.1 on the stable side, with memory capacity from 4 to 16.

### C7. Population mean rate against bias

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_C/bias_rate_lambda1/Fig_bias_vs_rate_mean.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_C\bias_rate_lambda1\Fig_bias_vs_rate_mean.png`

**Question:** How strongly does the bias set the network's mean firing rate in each condition?

- With adaptation, the bias sets the mean rate almost by itself (r = -0.94 for 1TS, -0.97 for chained and product MTS). Without adaptation it does so less tightly (r = -0.67), and gain moves the rate widely.

### C8. $\lambda_1$ against bias

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_C/bias_rate_lambda1/Fig_bias_vs_lambda1.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_C\bias_rate_lambda1\Fig_bias_vs_lambda1.png`

**Question:** Does bias, independently of gain, push the network toward chaos?

- Without adaptation, $\lambda_1$ does not depend on bias (r = 0.00); gain alone decides it. With adaptation there is a weak positive trend (r = 0.26 to 0.37): a lower bias keeps more high-gain networks stable.

### C9. Population mean rate against $\lambda_1$

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_C/bias_rate_lambda1/Fig_rate_mean_vs_lambda1.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_C\bias_rate_lambda1\Fig_rate_mean_vs_lambda1.png`

**Question:** Is chaos ($\lambda_1$ > 0) explained by the network's mean firing rate?

- No. Mean rate and $\lambda_1$ are essentially uncorrelated in every condition (Pearson r from -0.07 to -0.24): the low memory capacity on the chaotic side is chaos, not silence or saturation.

### C10. Population mean rate against $\lambda_1$, colored by bias

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_C/bias_rate_lambda1/Fig_colorbias_rate_mean_vs_lambda1.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_C\bias_rate_lambda1\Fig_colorbias_rate_mean_vs_lambda1.png`

**Question:** Seen by bias rather than gain, does the mean firing rate decide where the network turns unstable?

- No. Each bias level traces an arch as gain rises, and arches of every bias cross $\lambda_1 = 0$ at mean rates from about 0.1 to 0.6 without adaptation; with both MTS forms the stable runs sit at $\lambda_1 \approx -0.1$ across the whole range of mean rate (0 to 0.8).

### C11. Population mean rate against $\lambda_1$, marker size by bias

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_C/bias_rate_lambda1/Fig_rate_mean_vs_lambda1_size_bias.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_C\bias_rate_lambda1\Fig_rate_mean_vs_lambda1_size_bias.png`

**Question:** Does bias organize the relation between rate and $\lambda_1$?

- Yes, along the rate axis: negative biases (rings) sit at high mean rates and positive biases (filled markers) at low ones. Gain sets which side of $\lambda_1$ = 0 a network lands on.

### C12. Share of neurons on the nonlinearity's slope against gain

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_C/bias_rate_lambda1/Fig_gain_vs_frac_on_slope.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_C\bias_rate_lambda1\Fig_gain_vs_frac_on_slope.png`

**Question:** How does the share of neurons neither pinned nor silent change with gain, and does adaptation keep neurons in the responsive range?

- Without adaptation, 25 to 95 % of neurons are on the slope, and many networks drop below 50 % at low and middle gains. With adaptation most networks keep 60 to 100 % on the slope across gain; the exception is bias 0.6, where only about 25 % are on the slope at low gain.

## Gain by self-negative feedback (one seed, first look)

### C13. Memory capacity and $\lambda_1$ across gain and self-feedback

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/block_C/selffb/four_conditions/Fig_MC_gain_bias_maps.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\block_C\selffb\four_conditions\Fig_MC_gain_bias_maps.png`

**Question:** The bias could not move the no-adaptation edge of chaos; does self-negative feedback on the diagonal of W ($W = g,W_mathrm{drawn} - k,R_W I$, at the preset bias 0.2) move it, and does memory without adaptation then follow?

- The edge moves to higher gain as k grows: gain 0.33 to 0.5 at k = 0, 1 to 1.4 at k = 0.5, 1.4 to 2 at k = 1, and none up to gain 2 at k = 2.
- Memory without adaptation does not follow: its best at k of 0.5 or more is 6.3, below its k = 0 peak of 8.1 and far below the adapting conditions' peaks (13.4 to 15.5; chained MTS 15.5 at gain 1.4, k = 0.5).

### C14. The same maps, linear scales

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/block_C/selffb/four_conditions/Fig_MC_gain_selffb_maps_linear.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\block_C\selffb\four_conditions\Fig_MC_gain_selffb_maps_linear.png`

**Question:** Do the gain-by-self-feedback maps read the same on linear color scales?

- At k = 2 no condition turns chaotic; memory stays at 3 to 7 with adaptation, and the no-adaptation network sits at the leak floor ($\lambda_1$ about -10/s, beyond the plus or minus 5 clip) with little memory.
- With adaptation the ridge of high memory follows the edge diagonally up to gain 2 at k = 1.

### C15. Memory capacity and $\lambda_1$ across gain, no self-feedback

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/block_C/selffb/four_conditions/Fig_MC_gain_curves.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\block_C\selffb\four_conditions\Fig_MC_gain_curves.png`

**Question:** Without self-feedback (k = 0), how do memory capacity and $\lambda_1$ change with gain in this run?

- The k = 0 row reproduces the main run's seed-3001 values at bias 0.2.
- No adaptation keeps memory near 8 only at gains 0.22 and 0.33 and loses it from gain 0.5, where $\lambda_1$ turns positive; peaks 13.7 for 1TS (gain 0.7), about 12.7 for chained MTS (0.7) and 13.4 for product MTS (1).

### C16. Robustness of memory capacity over the self-feedback map

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/block_C/selffb/four_conditions/Fig_MC_robustness.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\block_C\selffb\four_conditions\Fig_MC_robustness.png`

**Question:** Over all gains and self-feedback levels, what share of the networks keeps useful memory in each condition?

- Median memory capacity: 1.03 (none), 5.50 (1TS), 5.80 (chained MTS), 6.04 (product MTS); share of the map at or above 7.77 (half the best of any condition): 0.04, 0.29, 0.29, 0.29.

### C17. Memory capacity against $\lambda_1$, every run, colored by gain

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/block_C/selffb/four_conditions/Fig_MC_vs_lambda1.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\block_C\selffb\four_conditions\Fig_MC_vs_lambda1.png`

**Question:** With self-feedback added, does memory capacity still track $\lambda_1$?

- Pearson r over all 28 runs: 0.15 (none, 19 runs clipped at $|\lambda_1| > 5$), -0.55 (1TS), -0.51 (chained), -0.45 (product). On the stable side only chained MTS tracks it (r = -0.80, n = 24; 1TS and product 0.03).

### C18. Memory capacity against $\lambda_1$, marker size by self-feedback

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/block_C/selffb/four_conditions/Fig_MC_vs_lambda1_size_selffb.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\block_C\selffb\four_conditions\Fig_MC_vs_lambda1_size_selffb.png`

**Question:** Where do the self-feedback levels fall along the relation between memory capacity and $\lambda_1$?

- With adaptation every level sits just below $\lambda_1$ = 0, and the largest self-feedback (k = 2) gives the lowest memory (3 to 7); without adaptation the large-k runs are clipped below $\lambda_1$ = -5 with memory under 4.
- Reading (lead's, FD-024): a fixed, instantaneous self-inhibition buys stability but not memory; adaptation's slow, rate-dependent feedback buys both.
