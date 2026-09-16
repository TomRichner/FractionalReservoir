# Memory capacity run MC_sample_hold_20260915_211044_trials5

## Source

- commit: `fb6a3fb814c094181b904e9dddc1d0ff7a570cdd` (fb6a3fb) on branch `main`
- host: R5611351 (PCWIN64), MATLAB 2026a
- preset: `celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0_dualStd_3cond_mu7revised` (SRNN_ESN_reservoir), run mode `fast`
- captured: 15-Sep-2026 21:16:05

## Protocol

- n = 500, fs = 200 Hz, integrator `sra1`, sigma_u_noise = 0
- input `sample_hold`, T_hold = 0.3 s (MC in hold units), T_wash = 10 s, T_train = 600 s, T_test = 150 s, d_max = 15 s (50 delays scored)
- readout `synaptic`, horizon threshold R^2 > 0.10
- bootstrap 2000 samples; sign-flip test exact up to N = 20, else 10000 Monte Carlo patterns
- seed bases: net 3000, stim 4000

## Trials

Attempted 5, completed 5, failed 0 (a failed trial aborts the run; see the header of write_provenance_md).

| trial | seed_net | seed_stim | complete | MC no_adaptation | MC sfa1_std1 | MC sfa3_std2 |
|---|---|---|---|---|---|---|
| 1 | 3001 | 4001 | 1 | 0.084 | 1.348 | 14.685 |
| 2 | 3002 | 4002 | 1 | 0.111 | 2.051 | 14.715 |
| 3 | 3003 | 4003 | 1 | 0.119 | 3.239 | 14.552 |
| 4 | 3004 | 4004 | 1 | 0.127 | 2.015 | 14.265 |
| 5 | 3005 | 4005 | 1 | 0.109 | 2.855 | 12.846 |

## Summary

| Condition | Total MC mean [95% CI] | Horizon (s) mean [95% CI] |
|---|---|---|
| no_adaptation | 0.110 [0.097, 0.122] | 0.000 [0.000, 0.000] |
| sfa1_std1 | 2.302 [1.755, 2.918] | 1.620 [1.380, 1.860] |
| sfa3_std2 | 14.213 [13.504, 14.671] | 8.700 [7.920, 9.240] |

| Pair | mean diff | p (sign-flip) | patterns | Cohen's d_z |
|---|---|---|---|---|
| no_adaptation vs sfa1_std1 | -2.192 | 0.0625 | 32 (exact) | -2.963 |
| no_adaptation vs sfa3_std2 | -14.103 | 0.0625 | 32 (exact) | -17.929 |
| sfa1_std1 vs sfa3_std2 | -11.911 | 0.0625 | 32 (exact) | -9.162 |

## Outputs

- `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\mu7revised_fast\memory_capacity\MC_sample_hold_20260915_211044_trials5_results.mat`
- `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\mu7revised_fast\memory_capacity\MC_sample_hold_20260915_211044_trials5_MC_Horizon.csv`
- `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\mu7revised_fast\memory_capacity\MC_sample_hold_20260915_211044_trials5_summary.txt`
- `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\mu7revised_fast\memory_capacity\git_provenance.txt`
