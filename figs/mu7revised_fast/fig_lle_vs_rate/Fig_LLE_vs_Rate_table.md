# lambda_1 vs mean firing rate

Run: `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\mu7revised_fast`. Sources: param_space sample + 1-D sweeps (f_E, level_of_chaos, mu_EE_relative, mu_EI_relative, mu_IE_relative, mu_II_relative, n). Spearman rank correlation with a 1000-sample bootstrap 95% CI; Pearson deliberately not used (nonlinear, saturating relation). Quiet: mean rate < 0.02; saturated: > 0.9.

| Condition | Spearman rho [95% CI] | n | n quiet | n saturated | lambda_1 median quiet | mid | saturated |
|---|---|---|---|---|---|---|---|
| no_adaptation | -0.32 [-0.44, -0.19] | 303 | 28 | 66 | +0.002 (n=28) | +1.463 (n=209) | -9.993 (n=66) |
| sfa1_std1 | -0.41 [-0.51, -0.30] | 303 | 28 | 27 | -0.506 (n=28) | -0.511 (n=248) | -1.899 (n=27) |
| sfa3_std2 | +0.09 [-0.02, +0.20] | 303 | 21 | 3 | -0.102 (n=21) | -0.103 (n=279) | -0.092 (n=3) |
