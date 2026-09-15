# Memory capacity

Source: `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\sfaEI_mu5_fast\memory_capacity\MC_sample_hold_20260914_181458_trials5_results.mat`

Preset `celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0_dualStd_3cond_mu5`, run mode `fast`, 5 paired trials, readout `synaptic`, horizon threshold R^2 > 0.10, 2000 bootstrap samples.

| Condition | Total MC mean [95% CI] | Horizon (s) mean [95% CI] | median MC | median horizon |
|---|---|---|---|---|
| No Adaptation | 0.117 [0.094, 0.136] | 0.000 [0.000, 0.000] | 0.120 | 0.000 |
| Single-Timescale Adaptation | 15.813 [15.303, 16.563] | 6.660 [6.480, 6.840] | 15.406 | 6.600 |
| Multiple-Timescale Adaptation | 10.556 [10.082, 11.044] | 6.960 [6.540, 7.440] | 10.637 | 6.900 |

| Pair | mean diff (total MC) | p (sign-flip) | patterns | Cohen's d_z |
|---|---|---|---|---|
| no_adaptation vs sfa1_std1 | -15.696 | 0.062 | 32 (exact) | -18.27 |
| no_adaptation vs sfa3_std2 | -10.439 | 0.062 | 32 (exact) | -16.23 |
| sfa1_std1 vs sfa3_std2 | +5.257 | 0.062 | 32 (exact) | 5.01 |
