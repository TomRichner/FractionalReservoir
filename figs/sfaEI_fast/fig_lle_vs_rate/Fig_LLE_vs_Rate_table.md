# lambda_1 vs mean firing rate

Run: `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\sfaEI_fast`. Sources: param_space sample + 1-D sweeps (f_E, level_of_chaos, mu_EE_relative, mu_EI_relative, mu_IE_relative, mu_II_relative, n). Spearman rank correlation with a 1000-sample bootstrap 95% CI; Pearson deliberately not used (nonlinear, saturating relation). Quiet: mean rate < 0.02; saturated: > 0.9.

| Condition | Spearman rho [95% CI] | n | n quiet | n saturated | lambda_1 median quiet | mid | saturated |
|---|---|---|---|---|---|---|---|
| no_adaptation | -0.46 [-0.63, -0.25] | 111 | 13 | 26 | +0.556 (n=13) | +2.769 (n=72) | -6.430 (n=26) |
| sfa1_std1 | -0.47 [-0.63, -0.28] | 111 | 14 | 12 | -0.504 (n=14) | +0.068 (n=85) | -3.574 (n=12) |
| sfa3_std2 | -0.03 [-0.22, +0.14] | 111 | 12 | 1 | -0.114 (n=12) | -0.109 (n=98) | -0.102 (n=1) |
