# Memory capacity run MC_sample_hold_20260914_092045_trials5

## Source

- commit: `ea9e756a2dc2ec494739ec358dadf4d01a76f6a1` (ea9e756) on branch `main`
- host: R5611351 (PCWIN64), MATLAB 2026a
- preset: `celltype_pairs_sfaEI_Sc0p2sig0p1_noise0_dualStd_3cond_mu8p25` (SRNN_ESN_reservoir), run mode `fast`
- captured: 14-Sep-2026 09:26:05

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
| 1 | 3001 | 4001 | 1 | 0.107 | 0.323 | 17.774 |
| 2 | 3002 | 4002 | 1 | 0.095 | 0.402 | 15.678 |
| 3 | 3003 | 4003 | 1 | 0.120 | 0.560 | 16.362 |
| 4 | 3004 | 4004 | 1 | 0.118 | 0.509 | 17.376 |
| 5 | 3005 | 4005 | 1 | 0.090 | 0.610 | 14.415 |

## Summary

| Condition | Total MC mean [95% CI] | Horizon (s) mean [95% CI] |
|---|---|---|
| no_adaptation | 0.106 [0.096, 0.117] | 0.000 [0.000, 0.000] |
| sfa1_std1 | 0.481 [0.387, 0.570] | 0.360 [0.300, 0.480] |
| sfa3_std2 | 16.321 [15.260, 17.332] | 10.560 [9.420, 12.060] |

| Pair | mean diff | p (sign-flip) | patterns | Cohen's d_z |
|---|---|---|---|---|
| no_adaptation vs sfa1_std1 | -0.375 | 0.0625 | 32 (exact) | -3.189 |
| no_adaptation vs sfa3_std2 | -16.215 | 0.0625 | 32 (exact) | -12.113 |
| sfa1_std1 vs sfa3_std2 | -15.840 | 0.0625 | 32 (exact) | -11.126 |

## Outputs

- `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\sfaEI_fast\memory_capacity\MC_sample_hold_20260914_092045_trials5_results.mat`
- `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\sfaEI_fast\memory_capacity\MC_sample_hold_20260914_092045_trials5_MC_Horizon.csv`
- `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\sfaEI_fast\memory_capacity\MC_sample_hold_20260914_092045_trials5_summary.txt`
- `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\sfaEI_fast\memory_capacity\git_provenance.txt`
