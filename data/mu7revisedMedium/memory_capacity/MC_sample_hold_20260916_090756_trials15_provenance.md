# Memory capacity run MC_sample_hold_20260916_090756_trials15

## Source

- commit: `a413e198dd2e7bef504c6802212773cda539895d` (a413e19) on branch `main`
- host: R5611351 (PCWIN64), MATLAB 2026a
- preset: `celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0_dualStd_3cond_mu7revisedMedium` (SRNN_ESN_reservoir), run mode `medium`
- captured: 16-Sep-2026 09:47:39

## Protocol

- n = 500, fs = 200 Hz, integrator `sra1`, sigma_u_noise = 0
- input `sample_hold`, T_hold = 0.3 s (MC in hold units), T_wash = 10 s, T_train = 600 s, T_test = 150 s, d_max = 15 s (50 delays scored)
- readout `synaptic`, horizon threshold R^2 > 0.10
- bootstrap 2000 samples; sign-flip test exact up to N = 20, else 10000 Monte Carlo patterns
- seed bases: net 3000, stim 4000

## Trials

Attempted 15, completed 15, failed 0 (a failed trial aborts the run; see the header of write_provenance_md).

| trial | seed_net | seed_stim | complete | MC no_adaptation | MC sfa1_std1 | MC sfa3_std2 |
|---|---|---|---|---|---|---|
| 1 | 3001 | 4001 | 1 | 0.084 | 1.348 | 14.685 |
| 2 | 3002 | 4002 | 1 | 0.111 | 2.051 | 14.715 |
| 3 | 3003 | 4003 | 1 | 0.119 | 3.239 | 14.552 |
| 4 | 3004 | 4004 | 1 | 0.127 | 2.015 | 14.265 |
| 5 | 3005 | 4005 | 1 | 0.109 | 2.855 | 12.846 |
| 6 | 3006 | 4006 | 1 | 0.082 | 1.878 | 13.697 |
| 7 | 3007 | 4007 | 1 | 0.143 | 9.294 | 13.464 |
| 8 | 3008 | 4008 | 1 | 0.087 | 12.211 | 12.686 |
| 9 | 3009 | 4009 | 1 | 0.077 | 1.972 | 11.154 |
| 10 | 3010 | 4010 | 1 | 0.097 | 6.317 | 13.308 |
| 11 | 3011 | 4011 | 1 | 0.137 | 6.095 | 14.140 |
| 12 | 3012 | 4012 | 1 | 0.108 | 1.963 | 14.707 |
| 13 | 3013 | 4013 | 1 | 0.094 | 1.806 | 12.842 |
| 14 | 3014 | 4014 | 1 | 0.112 | 2.069 | 13.321 |
| 15 | 3015 | 4015 | 1 | 0.096 | 2.510 | 14.582 |

## Summary

| Condition | Total MC mean [95% CI] | Horizon (s) mean [95% CI] |
|---|---|---|
| no_adaptation | 0.106 [0.096, 0.116] | 0.000 [0.000, 0.000] |
| sfa1_std1 | 3.841 [2.473, 5.577] | 2.340 [1.680, 3.180] |
| sfa3_std2 | 13.664 [13.139, 14.124] | 8.500 [8.070, 8.940] |

| Pair | mean diff | p (sign-flip) | patterns | Cohen's d_z |
|---|---|---|---|---|
| no_adaptation vs sfa1_std1 | -3.736 | 6.104e-05 | 32768 (exact) | -1.162 |
| no_adaptation vs sfa3_std2 | -13.559 | 6.104e-05 | 32768 (exact) | -13.496 |
| sfa1_std1 vs sfa3_std2 | -9.823 | 6.104e-05 | 32768 (exact) | -2.735 |

## Outputs

- `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\mu7revisedMedium\memory_capacity\MC_sample_hold_20260916_090756_trials15_results.mat`
- `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\mu7revisedMedium\memory_capacity\MC_sample_hold_20260916_090756_trials15_MC_Horizon.csv`
- `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\mu7revisedMedium\memory_capacity\MC_sample_hold_20260916_090756_trials15_summary.txt`
- `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\mu7revisedMedium\memory_capacity\git_provenance.txt`
