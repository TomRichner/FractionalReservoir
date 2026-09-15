# lambda_1 vs mean firing rate

Run: `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\sfaEI_mu5_fast`. Sources: param_space sample + 1-D sweeps (f_E, level_of_chaos, mu_EE_relative, mu_EI_relative, mu_IE_relative, mu_II_relative, n). Spearman rank correlation with a 1000-sample bootstrap 95% CI; Pearson deliberately not used (nonlinear, saturating relation). Quiet: mean rate < 0.02; saturated: > 0.9.

| Condition | Spearman rho [95% CI] | n | n quiet | n saturated | lambda_1 median quiet | mid | saturated |
|---|---|---|---|---|---|---|---|
| no_adaptation | -0.50 [-0.62, -0.37] | 212 | 24 | 53 | -2.609 (n=24) | -0.443 (n=135) | -9.996 (n=53) |
| sfa1_std1 | -0.52 [-0.63, -0.40] | 212 | 17 | 30 | -0.508 (n=17) | -0.538 (n=165) | -1.853 (n=30) |
| sfa3_std2 | +0.27 [+0.17, +0.39] | 212 | 14 | 2 | -0.103 (n=14) | -0.101 (n=196) | -0.094 (n=2) |
