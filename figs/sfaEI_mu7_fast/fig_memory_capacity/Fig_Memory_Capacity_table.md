# Memory capacity

Source: `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\sfaEI_mu7_fast\memory_capacity\MC_sample_hold_20260915_131030_trials5_results.mat`

Preset `celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0_dualStd_3cond_mu7`, run mode `fast`, 5 paired trials, readout `synaptic`, horizon threshold R^2 > 0.10, 2000 bootstrap samples.

| Condition | Total MC mean [95% CI] | Horizon (s) mean [95% CI] | median MC | median horizon |
|---|---|---|---|---|
| No Adaptation | 0.110 [0.097, 0.122] | 0.000 [0.000, 0.000] | 0.111 | 0.000 |
| Single-Timescale Adaptation | 2.302 [1.755, 2.918] | 1.620 [1.380, 1.860] | 2.051 | 1.500 |
| Multiple-Timescale Adaptation | 14.213 [13.504, 14.671] | 8.700 [7.920, 9.240] | 14.552 | 9.000 |

| Pair | mean diff (total MC) | p (sign-flip) | patterns | Cohen's d_z |
|---|---|---|---|---|
| no_adaptation vs sfa1_std1 | -2.192 | 0.062 | 32 (exact) | -2.96 |
| no_adaptation vs sfa3_std2 | -14.103 | 0.062 | 32 (exact) | -17.93 |
| sfa1_std1 vs sfa3_std2 | -11.911 | 0.062 | 32 (exact) | -9.16 |
