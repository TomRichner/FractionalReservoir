# lambda_1 vs mean firing rate

Run: `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\mu7revisedMedium`. Sources: param_space sample + 1-D sweeps (f_E, level_of_chaos, mu_EE_relative, mu_EI_relative, mu_IE_relative, mu_II_relative, n). Spearman rank correlation with a 1000-sample bootstrap 95% CI; Pearson deliberately not used (nonlinear, saturating relation). Quiet: mean rate < 0.02; saturated: > 0.9.

| Condition | Spearman rho [95% CI] | n | n quiet | n saturated | lambda_1 median quiet | mid | saturated |
|---|---|---|---|---|---|---|---|
| no_adaptation | -0.48 [-0.53, -0.43] | 1283 | 83 | 210 | +0.272 (n=83) | +2.113 (n=990) | -6.233 (n=210) |
| sfa1_std1 | -0.41 [-0.45, -0.36] | 1283 | 70 | 77 | -0.508 (n=70) | -0.480 (n=1136) | -1.908 (n=77) |
| sfa3_std2 | +0.41 [+0.36, +0.45] | 1283 | 40 | 2 | -0.104 (n=40) | -0.101 (n=1241) | -0.087 (n=2) |
