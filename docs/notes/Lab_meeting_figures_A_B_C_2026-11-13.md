---
title: "Lab meeting figures: blocks A, B and C"
date: 2026-11-13
header-includes:
  - |
    ```{=latex}
    \AtBeginDocument{\let\oldtextbackslash\textbackslash
    \renewcommand{\textbackslash}{\oldtextbackslash\allowbreak}\sloppy}
    ```
---

Every four-condition figure (no adaptation, single-timescale, chained MTS, product MTS) of blocks A, B and C of the final-draft work so far.

**Equations of the model:** [Equations_chained_MTS.pdf](file:///C:/Users/m218089/Desktop/github_repos/FractionalReservoir/docs/EquationsParametersDocs/chained_MTS/Equations_chained_MTS.pdf) (the chained multiple-timescale model, its conditions, and in section 10 the Jacobian in block form), at `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\docs\EquationsParametersDocs\chained_MTS\Equations_chained_MTS.pdf`.

# A. Chained vesicle-pool model and the essential run

## A1. Example dynamics (main group 2)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_main2_grouped/Fig_Main2_Grouped.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_main2_grouped\Fig_Main2_Grouped.png`

**Question:** How do the representative network trajectories (states, local rate and finite-time $\lambda_1$) differ across the four adaptation conditions after a uniform input step?

- Columns are no adaptation, 1TS, chained MTS and product MTS, each with the same 12 E and 12 I neurons; row G shows the local rate (low-pass filtered at 2 Hz for display only) and the finite-time $\lambda_1$.

## A2. Eigenvalues and switching (main group 3)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_main3_grouped/Fig_Main3_Grouped.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_main3_grouped\Fig_Main3_Grouped.png`

**Question:** How do the Jacobian eigenvalue occupancy, the local leading Lyapunov rates and the expansive transients differ by condition?

- Row C shows the local and accumulated leading rates (K = 50, accumulated over 1 to 40 s); row D the positive local QR growth in bits/s, which is not a measure of KS entropy.

## A3. Sensitivity of $\lambda_1$ and $\lambda_1$ against rate (main group 4)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_main4_grouped/Fig_Main4_Grouped.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_main4_grouped\Fig_Main4_Grouped.png`

**Question:** How does the median $\lambda_1$ respond to the E:I ratio, network size and the four weight blocks, and how does it relate to the firing rate in each condition?

- No adaptation swings far beyond plus or minus 1/s across every sweep; both MTS conditions stay near zero for most parameter values (chained MTS rises to about +1/s at high E to I weights and about +0.5/s at moderate E to E and I to E weights).
- In panel B, chained and product MTS sit at $\lambda_1$ near 0 across the whole rate range; no adaptation falls to about -5/s as it saturates.

## A4. $\lambda_1$ sensitivity distributions (core sweeps)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_analysis_allStd/Fig_Sensitivity_LLE_core.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_analysis_allStd\Fig_Sensitivity_LLE_core.png`

**Question:** What is the full distribution of $\lambda_1$ across the core 1-D sweeps (E:I ratio, network size, overall gain) for each condition?

- Distributions per condition and sweep level; the essential run has 4 levels x 2 reps per sweep.

## A5. $\lambda_1$ sensitivity distributions (weight blocks)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_analysis_allStd/Fig_Sensitivity_LLE_mu.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_analysis_allStd\Fig_Sensitivity_LLE_mu.png`

**Question:** What is the distribution of $\lambda_1$ across the four weight-block sweeps for each condition?

- The four mu sweeps span -75 % to +75 % of the preset value, at 4 levels x 2 reps.

## A6. Mean rate sensitivity (core sweeps)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_analysis_allStd/Fig_Sensitivity_mean_rate_core.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_analysis_allStd\Fig_Sensitivity_mean_rate_core.png`

**Question:** How does the mean firing rate depend on the core sweep parameters in each condition?

- Distributions of the population mean rate per sweep level for all four conditions (4 levels x 2 reps).

## A7. Mean rate sensitivity (weight blocks)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_analysis_allStd/Fig_Sensitivity_mean_rate_mu.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_analysis_allStd\Fig_Sensitivity_mean_rate_mu.png`

**Question:** How does the mean firing rate depend on the four weight blocks in each condition?

- Distributions of the population mean rate across the mu sweeps for all four conditions.

## A8. KS entropy sensitivity (core sweeps)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_analysis_allStd/Fig_Sensitivity_hKS_core.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_analysis_allStd\Fig_Sensitivity_hKS_core.png`

**Question:** How does the Kolmogorov-Sinai entropy h_KS vary across the core sweeps in each condition?

- h_KS distributions (bit/s) per sweep level for all four conditions.

## A9. KS entropy sensitivity (weight blocks)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_analysis_allStd/Fig_Sensitivity_hKS_mu.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_analysis_allStd\Fig_Sensitivity_hKS_mu.png`

**Question:** How does h_KS vary across the four weight-block sweeps in each condition?

- h_KS distributions per mu sweep level for all four conditions.

## A10. Kaplan-Yorke dimension sensitivity (core sweeps)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_analysis_allStd/Fig_Sensitivity_DKY_core.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_analysis_allStd\Fig_Sensitivity_DKY_core.png`

**Question:** How does the attractor dimension D_KY vary across the core sweeps in each condition?

- D_KY distributions per sweep level for all four conditions.

## A11. Kaplan-Yorke dimension sensitivity (weight blocks)

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_analysis_allStd/Fig_Sensitivity_DKY_mu.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_analysis_allStd\Fig_Sensitivity_DKY_mu.png`

**Question:** How does D_KY vary across the four weight-block sweeps in each condition?

- D_KY distributions per mu sweep level for all four conditions.

## A12. $\lambda_1$ sensitivity medians

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_medians/Fig_Sensitivity_LLE_medians.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_medians\Fig_Sensitivity_LLE_medians.png`

**Question:** How does the median $\lambda_1$ (with IQR) move with each swept parameter in each condition?

- Same data as main group 4, panel A: both MTS conditions stay close to zero, while no adaptation varies widely.

## A13. Mean rate sensitivity medians

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_medians/Fig_Sensitivity_mean_rate_medians.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_medians\Fig_Sensitivity_mean_rate_medians.png`

**Question:** How does the median mean firing rate move with each swept parameter in each condition?

- Medians and IQR of the mean rate per sweep for all four conditions.

## A14. KS entropy sensitivity medians

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_medians/Fig_Sensitivity_hKS_medians.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_medians\Fig_Sensitivity_hKS_medians.png`

**Question:** How does the median h_KS move with each swept parameter in each condition?

- Medians and IQR of h_KS per sweep for all four conditions.

## A15. Kaplan-Yorke dimension sensitivity medians

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_sensitivity_medians/Fig_Sensitivity_DKY_medians.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_sensitivity_medians\Fig_Sensitivity_DKY_medians.png`

**Question:** How does the median D_KY move with each swept parameter in each condition?

- Medians and IQR of D_KY per sweep for all four conditions.

## A16. Joint parameter space: $\lambda_1$

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_param_space_allStd/Fig_ParamSpace_LLE.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_param_space_allStd\Fig_ParamSpace_LLE.png`

**Question:** Over a joint random sample of parameter space, how is $\lambda_1$ distributed in each condition?

- No adaptation splits into two modes, about 70 % at the negative clip and about 30 % at the positive clip; 1TS is spread out, with a mode near -0.4/s.
- Chained MTS has about 70 % and product MTS about 75 % of points in the bin just below zero.

## A17. Joint parameter space: mean rate

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_param_space_allStd/Fig_ParamSpace_mean_rate.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_param_space_allStd\Fig_ParamSpace_mean_rate.png`

**Question:** Over the joint parameter sample, how is the mean firing rate distributed in each condition?

- Histograms of the mean rate over the joint sample, one per condition.

## A18. Joint parameter space: KS entropy

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_param_space_allStd/Fig_ParamSpace_hKS.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_param_space_allStd\Fig_ParamSpace_hKS.png`

**Question:** Over the joint parameter sample, how is h_KS distributed in each condition?

- Histograms of h_KS over the joint sample, one per condition.

## A19. Joint parameter space: Kaplan-Yorke dimension

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_param_space_allStd/Fig_ParamSpace_DKY.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_param_space_allStd\Fig_ParamSpace_DKY.png`

**Question:** Over the joint parameter sample, how is D_KY distributed in each condition?

- Histograms of D_KY over the joint sample, one per condition.

## A20. E:I weight balance: $\lambda_1$

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_EI_weights_param_space/Fig_EI_Weights_ParamSpace_LLE.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_EI_weights_param_space\Fig_EI_Weights_ParamSpace_LLE.png`

**Question:** How does $\lambda_1$ depend on the realized E:I weight balance in each condition?

- $\lambda_1$ against the E:I weight balance across the joint sample, one panel per condition.

## A21. E:I weight balance: mean rate

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_EI_weights_param_space/Fig_EI_Weights_ParamSpace_mean_rate.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_EI_weights_param_space\Fig_EI_Weights_ParamSpace_mean_rate.png`

**Question:** How does the mean rate depend on the E:I weight balance in each condition?

- Mean rate against the E:I weight balance, one panel per condition.

## A22. E:I weight balance: KS entropy

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_EI_weights_param_space/Fig_EI_Weights_ParamSpace_hKS.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_EI_weights_param_space\Fig_EI_Weights_ParamSpace_hKS.png`

**Question:** How does h_KS depend on the E:I weight balance in each condition?

- h_KS against the E:I weight balance, one panel per condition.

## A23. E:I weight balance: Kaplan-Yorke dimension

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_EI_weights_param_space/Fig_EI_Weights_ParamSpace_DKY.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_EI_weights_param_space\Fig_EI_Weights_ParamSpace_DKY.png`

**Question:** How does D_KY depend on the E:I weight balance in each condition?

- D_KY against the E:I weight balance, one panel per condition.

## A24. $\lambda_1$ against mean firing rate

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_lle_vs_rate/Fig_LLE_vs_Rate.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_lle_vs_rate\Fig_LLE_vs_Rate.png`

**Question:** Is $\lambda_1$ tied to the firing rate, and does adaptation remove that dependence?

- Spearman rho [95 % CI], n = 144 per condition: no adaptation -0.40 [-0.56, -0.23], 1TS -0.43 [-0.57, -0.29], chained MTS -0.25 [-0.37, -0.12], product MTS +0.08 [-0.08, +0.21].
- Median $\lambda_1$ in the mid-rate band: +2.24 (none), -0.35 (1TS), -0.10 (chained), -0.10 (product); in saturated networks: -9.59 (n = 25), -1.93 (n = 10), -0.10 (n = 9), -0.09 (n = 2).

## A25. Jacobian eigenvalue heatmap

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_4cond/fig_eig_heatmap/fig_eig_heatmap.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_4cond\fig_eig_heatmap\fig_eig_heatmap.png`

**Question:** Where do the instantaneous Jacobian eigenvalues sit along the trajectory in each condition, relative to the finite-time $\lambda_1$?

- Finite-time $\lambda_1$ over the last 20 s: no adaptation +3.092, 1TS +0.233, chained MTS -0.105, product MTS -0.098.
- All three adapting conditions pile up eigenvalue density just left of Re = 0; the product MTS bulk is the most compact, reaching to about plus or minus 10 on the imaginary axis.

## A26. Memory capacity

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/mu7chain_essential_mc_4cond/fig_memory_capacity/Fig_Memory_Capacity.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\mu7chain_essential_mc_4cond\fig_memory_capacity\Fig_Memory_Capacity.png`

**Question:** How much past input can a linear readout of the synaptic output reconstruct in each condition, and for how long?

- Total memory capacity, mean [95 % CI]: no adaptation 0.105, 1TS 2.89 [1.74, 4.17], chained MTS 9.61 [7.06, 13.29], product MTS 14.04 [13.55, 14.31]; memory horizon 0 s, 2.0 s, 5.7 s and 8.7 s in the same order.
- Only 3 paired trials, so every pairwise sign-flip test has p = 0.25 (the smallest possible); chained against product MTS: difference -4.43, Cohen's d_z -1.39.

Not included: the 1TS-versus-MTS depression figure (`Fig_STD_1TS_vs_MTS`), whose latest version shows three conditions (no no-adaptation trace).

# B. Stimulation pilots

## B1. P6: Chaos and rate across the 32 sampled networks

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_B/p6_hyperexcitability/p6_lle_vs_rate.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p6_hyperexcitability\p6_lle_vs_rate.png`

**Question:** Of the 32 networks sampled across the paper's parameter space, how many are hyperexcitable at baseline (chaotic, high rate, or both), and does adaptation move any network across those thresholds?

- In all four conditions about 10 of 32 networks are chaotic, 11 to 12 are high-rate (mean rate above 0.40), and only 2 of 32 are both (6 %, 95 % interval 2 to 20 %). Adaptation moves almost no network across a threshold: all 10 chaotic networks stay chaotic under chained MTS.
- The stored scalars miss aberrant activity: network 761 reads as stable and low-rate under chained MTS but bursts about every 4 s, so a parameter-space experiment needs an event measure and a large sample (about 145 to 1600 networks for 50 hyperexcitable ones).

Left out (fewer than four conditions; none of these pilots ran the no-adaptation condition): `run13_chain_mts` and `bursting_chain_mts` (3), P1a (2), P1b (1 to 3), P2 (1 per figure), P3 (2), P4 (3), P5 (2), and P6's parameter-space and time-series figures (2).

# C. Memory capacity across gain and bias

## C1. Memory capacity and $\lambda_1$ across gain and bias

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/block_C/main/four_conditions/Fig_MC_gain_bias_maps.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\block_C\main\four_conditions\Fig_MC_gain_bias_maps.png`

**Question:** Is the memory-capacity result a property of adaptation, or of one well-tuned network?

- In every condition, memory capacity is highest along a ridge just on the stable side of $\lambda_1$ = 0. Without adaptation the edge stays at gain 0.33 to 0.5 at every bias; with adaptation it moves to higher gain, and further the lower the bias is.
- At the preset bias, memory capacity stays at or above 7.5 over 0.55 octaves of gain without adaptation, against 2.0 to 2.4 octaves with adaptation. Peaks: 10.2 (none), 15.0 (1TS), 15.0 (chained MTS), 12.8 (product MTS).

## C2. Memory capacity across gain at three E:I balances

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/block_C/main/four_conditions/Fig_MC_gain_curves.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\block_C\main\four_conditions\Fig_MC_gain_curves.png`

**Question:** How do memory capacity and $\lambda_1$ change with gain at the preset bias, and does that change when excitation or inhibition dominates?

- When excitation dominates (3 : 2), no adaptation saturates (memory capacity at most 3.1, and 0 from gain 0.7), while every adapting condition keeps memory over the whole gain range; there chained MTS (6.5 to 9.7) is ahead of 1TS (6 to 7.8). One seed only.
- When inhibition dominates (2 : 3), the three adapting conditions are level (peaks 12.2 to 12.9 at gain 0.7). The order of the conditions changes from one network to the next.

## C3. Robustness of memory capacity over the map

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/block_C/main/four_conditions/Fig_MC_robustness.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\block_C\main\four_conditions\Fig_MC_robustness.png`

**Question:** Over all gains, biases and seeds, what share of the networks keeps useful memory in each condition?

- Without adaptation, fewer than 40 % of the networks have any memory (memory capacity above about 0.5); with adaptation, 65 to 85 % do.
- Above a memory capacity of about 8 the three adapting conditions are close together. Chained MTS against 1TS: higher in both seeds at 40 % of grid points, lower in both at 22 %, median difference +0.03.

## C4. Memory capacity against $\lambda_1$, every run

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/block_C/main/four_conditions/Fig_MC_vs_lambda1.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\block_C\main\four_conditions\Fig_MC_vs_lambda1.png`

**Question:** Does memory capacity track the distance to the edge of chaos ($\lambda_1$) in every condition?

- High memory capacity sits just below $\lambda_1$ = 0 in all four conditions; past zero it falls to about 0.1 to 1. Pearson r over all runs: -0.68 (none), -0.65 (1TS), -0.70 (chained), -0.61 (product).
- With chained and product MTS, $\lambda_1$ is pinned near -0.1/s across the whole stable region (by the slowest adaptation timescale), so it says little about how close a network is to the edge: on the stable side r = 0.87 without adaptation (n = 57) against -0.20 for product MTS (n = 120).

## C5. Memory capacity and $\lambda_1$ maps, linear scales

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_C/linear_scales/four_conditions/Fig_MC_gain_bias_maps_linear.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_C\linear_scales\four_conditions\Fig_MC_gain_bias_maps_linear.png`

**Question:** Do the gain-by-bias maps read the same when memory capacity and $\lambda_1$ are drawn on linear color scales?

- The ridge of high memory capacity hugs the stable side of the edge in every condition. The preset network (star) sits on the edge for 1TS and chained MTS, but inside the stable region for product MTS.
- The spread over the two seeds is concentrated at cells on the edge. Product MTS's two seeds agree least (r = 0.61).

## C6. Memory capacity against $\lambda_1$, exploration scatter

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_C/mc_vs_lambda1/Fig_MC_vs_lambda1_scatter.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_C\mc_vs_lambda1\Fig_MC_vs_lambda1_scatter.png`

**Question:** Within the stable and chaotic sides separately, how does memory capacity relate to $\lambda_1$ in each condition?

- Without adaptation, memory capacity rises toward the edge on the stable side (r = 0.87, n = 57); for product MTS, $\lambda_1$ is almost uninformative on the stable side (r = -0.20, n = 120).
- This exploration became C4 (same data and layout).

## C7. Memory capacity against $\lambda_1$, marker size by bias

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_C/bias_rate_lambda1/Fig_MC_vs_lambda1_size_bias.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_C\bias_rate_lambda1\Fig_MC_vs_lambda1_size_bias.png`

**Question:** Where do the different biases fall along the relation between memory capacity and $\lambda_1$?

- Points at all biases fall on the same curve. With MTS, networks from every bias collapse onto $\lambda_1$ near -0.1 on the stable side, with memory capacity from 4 to 16.

## C8. Population mean rate against bias

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_C/bias_rate_lambda1/Fig_bias_vs_rate_mean.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_C\bias_rate_lambda1\Fig_bias_vs_rate_mean.png`

**Question:** How strongly does the bias set the network's mean firing rate in each condition?

- With adaptation, the bias sets the mean rate almost by itself (r = -0.94 for 1TS, -0.97 for chained and product MTS). Without adaptation it does so less tightly (r = -0.67), and gain moves the rate widely.

## C9. $\lambda_1$ against bias

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_C/bias_rate_lambda1/Fig_bias_vs_lambda1.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_C\bias_rate_lambda1\Fig_bias_vs_lambda1.png`

**Question:** Does bias, independently of gain, push the network toward chaos?

- Without adaptation, $\lambda_1$ does not depend on bias (r = 0.00); gain alone decides it. With adaptation there is a weak positive trend (r = 0.26 to 0.37): a lower bias keeps more high-gain networks stable.

## C10. Population mean rate against $\lambda_1$

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_C/bias_rate_lambda1/Fig_rate_mean_vs_lambda1.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_C\bias_rate_lambda1\Fig_rate_mean_vs_lambda1.png`

**Question:** Is chaos ($\lambda_1$ > 0) explained by the network's mean firing rate?

- No. Mean rate and $\lambda_1$ are essentially uncorrelated in every condition (Pearson r from -0.07 to -0.24): the low memory capacity on the chaotic side is chaos, not silence or saturation.

## C11. Population mean rate against $\lambda_1$, colored by bias

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_C/bias_rate_lambda1/Fig_colorbias_rate_mean_vs_lambda1.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_C\bias_rate_lambda1\Fig_colorbias_rate_mean_vs_lambda1.png`

**Question:** Seen by bias rather than gain, does the mean firing rate decide where the network turns unstable?

- No. Each bias level traces an arch as gain rises, and arches of every bias cross $\lambda_1 = 0$ at mean rates from about 0.1 to 0.6 without adaptation; with both MTS forms the stable runs sit at $\lambda_1 \approx -0.1$ across the whole range of mean rate (0 to 0.8).

## C12. Population mean rate against $\lambda_1$, marker size by bias

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_C/bias_rate_lambda1/Fig_rate_mean_vs_lambda1_size_bias.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_C\bias_rate_lambda1\Fig_rate_mean_vs_lambda1_size_bias.png`

**Question:** Does bias organize the relation between rate and $\lambda_1$?

- Yes, along the rate axis: negative biases (rings) sit at high mean rates and positive biases (filled markers) at low ones. Gain sets which side of $\lambda_1$ = 0 a network lands on.

## C13. Share of neurons on the nonlinearity's slope against gain

![](C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/explore_C/bias_rate_lambda1/Fig_gain_vs_frac_on_slope.png){width=100%}

`C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_C\bias_rate_lambda1\Fig_gain_vs_frac_on_slope.png`

**Question:** How does the share of neurons neither pinned nor silent change with gain, and does adaptation keep neurons in the responsive range?

- Without adaptation, 25 to 95 % of neurons are on the slope, and many networks drop below 50 % at low and middle gains. With adaptation most networks keep 60 to 100 % on the slope across gain; the exception is bias 0.6, where only about 25 % are on the slope at low gain.

The self-feedback run (`figs\block_C\selffb\`) had produced no figures when this report was compiled, so it is not included.
