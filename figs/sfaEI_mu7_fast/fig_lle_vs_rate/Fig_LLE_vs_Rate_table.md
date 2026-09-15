# lambda_1 vs mean firing rate

Run: `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\sfaEI_mu7_fast`. Sources: param_space sample + 1-D sweeps (f_E, level_of_chaos, mu_EE_relative, mu_EI_relative, mu_IE_relative, mu_II_relative, n). Spearman rank correlation with a 1000-sample bootstrap 95% CI; Pearson deliberately not used (nonlinear, saturating relation). Quiet: mean rate < 0.02; saturated: > 0.9.

| Condition | Spearman rho [95% CI] | n | n quiet | n saturated | lambda_1 median quiet | mid | saturated |
|---|---|---|---|---|---|---|---|
| no_adaptation | -0.41 [-0.51, -0.29] | 303 | 16 | 73 | -0.304 (n=16) | +1.615 (n=214) | -9.996 (n=73) |
| sfa1_std1 | -0.44 [-0.54, -0.33] | 303 | 17 | 31 | -0.508 (n=17) | -0.511 (n=255) | -1.835 (n=31) |
| sfa3_std2 | +0.10 [-0.00, +0.21] | 303 | 11 | 3 | -0.104 (n=11) | -0.102 (n=289) | -0.091 (n=3) |
