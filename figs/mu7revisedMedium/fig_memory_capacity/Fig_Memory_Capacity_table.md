# Memory capacity

Source: `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\mu7revisedMedium\memory_capacity\MC_sample_hold_20260916_090756_trials15_results.mat`

Preset `celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0_dualStd_3cond_mu7revisedMedium`, run mode `medium`, 15 paired trials, readout `synaptic`, horizon threshold R^2 > 0.10, 2000 bootstrap samples.

| Condition | Total MC mean [95% CI] | Horizon (s) mean [95% CI] | median MC | median horizon |
|---|---|---|---|---|
| No Adaptation | 0.106 [0.096, 0.116] | 0.000 [0.000, 0.000] | 0.108 | 0.000 |
| Single-Timescale Adaptation | 3.841 [2.473, 5.577] | 2.340 [1.680, 3.180] | 2.069 | 1.500 |
| Multiple-Timescale Adaptation | 13.664 [13.139, 14.124] | 8.500 [8.070, 8.940] | 13.697 | 8.700 |

| Pair | mean diff (total MC) | p (sign-flip) | patterns | Cohen's d_z |
|---|---|---|---|---|
| no_adaptation vs sfa1_std1 | -3.736 | 6.1e-05 | 32768 (exact) | -1.16 |
| no_adaptation vs sfa3_std2 | -13.559 | 6.1e-05 | 32768 (exact) | -13.50 |
| sfa1_std1 vs sfa3_std2 | -9.823 | 6.1e-05 | 32768 (exact) | -2.74 |
