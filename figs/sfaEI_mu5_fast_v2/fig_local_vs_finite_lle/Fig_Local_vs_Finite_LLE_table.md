# Local rate vs finite-time lambda_1

Run: `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\sfaEI_mu5_fast`. Near-default = 1-D sweeps at the level nearest the preset default (f_E, level_of_chaos, mu_EE_relative, mu_EI_relative, mu_IE_relative, mu_II_relative, n); joint = param_space sample.

| Condition | set | n | lambda_1 median [IQR] | share lambda_1 < 0 | local rate share > 0 | median frac_local_positive | median mean excursion (s) | median p95 finite 0.2 s |
|---|---|---|---|---|---|---|---|---|
| no_adaptation | near-default | 21 | +1.853 [+0.072, +2.238] | 0.24 | 0.71 | 0.80 | 0.55 | +4.03 |
| no_adaptation | joint | 128 | -5.946 [-9.607, -2.609] | 0.81 | 0.17 | 0.00 | 0.33 | -5.20 |
| sfa1_std1 | near-default | 21 | -0.550 [-0.562, -0.545] | 1.00 | 0.01 | 0.00 | 0.10 | -0.45 |
| sfa1_std1 | joint | 128 | -0.533 [-0.786, -0.508] | 0.88 | 0.11 | 0.00 | 0.32 | -0.51 |
| sfa3_std2 | near-default | 21 | -0.111 [-0.116, -0.107] | 1.00 | 0.04 | 0.00 | 0.14 | -0.09 |
| sfa3_std2 | joint | 128 | -0.100 [-0.101, -0.095] | 0.87 | 0.09 | 0.00 | 0.24 | -0.09 |
