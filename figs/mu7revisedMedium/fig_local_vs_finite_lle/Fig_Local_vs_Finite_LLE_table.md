# Local rate vs finite-time lambda_1

Run: `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\mu7revisedMedium`. Near-default = 1-D sweeps at the level nearest the preset default (f_E, level_of_chaos, mu_EE_relative, mu_EI_relative, mu_IE_relative, mu_II_relative, n); joint = param_space sample.

| Condition | set | n | lambda_1 median [IQR] | share lambda_1 < 0 | local rate share > 0 | median frac_local_positive | median mean excursion (s) | median p95 finite 0.2 s |
|---|---|---|---|---|---|---|---|---|
| no_adaptation | near-default | 105 | +3.036 [+2.113, +3.610] | 0.01 | 0.87 | 0.91 | 0.82 | +5.88 |
| no_adaptation | joint | 128 | -4.833 [-8.077, -0.095] | 0.77 | 0.21 | 0.00 | 0.37 | -4.11 |
| sfa1_std1 | near-default | 105 | +0.009 [-0.229, +0.153] | 0.50 | 0.47 | 0.47 | 0.22 | +2.30 |
| sfa1_std1 | joint | 128 | -0.519 [-0.572, -0.505] | 0.87 | 0.12 | 0.00 | 0.31 | -0.51 |
| sfa3_std2 | near-default | 105 | -0.102 [-0.103, -0.101] | 1.00 | 0.00 | 0.00 | 0.05 | -0.09 |
| sfa3_std2 | joint | 128 | -0.100 [-0.102, -0.095] | 0.86 | 0.12 | 0.00 | 0.21 | -0.09 |
