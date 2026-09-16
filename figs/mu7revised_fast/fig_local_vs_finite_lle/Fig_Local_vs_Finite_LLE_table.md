# Local rate vs finite-time lambda_1

Run: `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\mu7revised_fast`. Near-default = 1-D sweeps at the level nearest the preset default (f_E, level_of_chaos, mu_EE_relative, mu_EI_relative, mu_IE_relative, mu_II_relative, n); joint = param_space sample.

| Condition | set | n | lambda_1 median [IQR] | share lambda_1 < 0 | local rate share > 0 | median frac_local_positive | median mean excursion (s) | median p95 finite 0.2 s |
|---|---|---|---|---|---|---|---|---|
| no_adaptation | near-default | 35 | +3.556 [+2.821, +3.915] | 0.00 | 0.93 | 0.95 | 1.19 | +5.83 |
| no_adaptation | joint | 128 | -5.559 [-8.622, -0.480] | 0.78 | 0.19 | 0.00 | 0.35 | -4.51 |
| sfa1_std1 | near-default | 35 | +0.323 [+0.027, +0.330] | 0.26 | 0.50 | 0.51 | 0.23 | +2.25 |
| sfa1_std1 | joint | 128 | -0.515 [-0.565, -0.505] | 0.85 | 0.14 | 0.00 | 0.29 | -0.50 |
| sfa3_std2 | near-default | 35 | -0.126 [-0.144, -0.123] | 1.00 | 0.00 | 0.00 | 0.05 | -0.10 |
| sfa3_std2 | joint | 128 | -0.099 [-0.101, -0.090] | 0.84 | 0.13 | 0.00 | 0.23 | -0.09 |
