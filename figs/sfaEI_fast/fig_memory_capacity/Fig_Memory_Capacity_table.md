# Memory capacity

Source: `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\sfaEI_fast\memory_capacity\MC_sample_hold_20260914_092045_trials5_results.mat`

Preset `celltype_pairs_sfaEI_Sc0p2sig0p1_noise0_dualStd_3cond_mu8p25`, run mode `fast`, 5 paired trials, readout `synaptic`, horizon threshold R^2 > 0.10, 2000 bootstrap samples.

| Condition | Total MC mean [95% CI] | Horizon (s) mean [95% CI] | median MC | median horizon |
|---|---|---|---|---|
| No Adaptation | 0.106 [0.096, 0.117] | 0.000 [0.000, 0.000] | 0.107 | 0.000 |
| Single-Timescale Adaptation | 0.481 [0.387, 0.570] | 0.360 [0.300, 0.480] | 0.509 | 0.300 |
| Multiple-Timescale Adaptation | 16.321 [15.260, 17.332] | 10.560 [9.420, 12.060] | 16.362 | 10.500 |

| Pair | mean diff (total MC) | p (sign-flip) | patterns | Cohen's d_z |
|---|---|---|---|---|
| no_adaptation vs sfa1_std1 | -0.375 | 0.062 | 32 (exact) | -3.19 |
| no_adaptation vs sfa3_std2 | -16.215 | 0.062 | 32 (exact) | -12.11 |
| sfa1_std1 vs sfa3_std2 | -15.840 | 0.062 | 32 (exact) | -11.13 |
