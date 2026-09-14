# Local rate vs finite-time lambda_1

Run: `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\sfaEI_fast`. Near-default = 1-D sweeps at the level nearest the preset default (f_E, level_of_chaos, mu_EE_relative, mu_EI_relative, mu_IE_relative, mu_II_relative, n); joint = param_space sample.

| Condition | set | n | lambda_1 median [IQR] | share lambda_1 < 0 | local rate share > 0 | median frac_local_positive | median mean excursion (s) | median p95 finite 0.2 s |
|---|---|---|---|---|---|---|---|---|
| no_adaptation | near-default | 21 | +3.652 [+0.066, +4.304] | 0.24 | 0.74 | 0.90 | 0.84 | +6.98 |
| no_adaptation | joint | 27 | -8.083 [-9.996, +4.054] | 0.70 | 0.27 | 0.00 | 0.41 | -5.81 |
| sfa1_std1 | near-default | 21 | +0.614 [-0.378, +1.404] | 0.29 | 0.54 | 0.60 | 0.30 | +3.44 |
| sfa1_std1 | joint | 27 | -0.514 [-1.888, +4.597] | 0.70 | 0.28 | 0.00 | 0.86 | -0.50 |
| sfa3_std2 | near-default | 21 | -0.118 [-0.136, +0.049] | 0.71 | 0.18 | 0.02 | 0.21 | -0.04 |
| sfa3_std2 | joint | 27 | -0.108 [-0.112, +2.333] | 0.67 | 0.29 | 0.00 | 0.53 | -0.09 |
