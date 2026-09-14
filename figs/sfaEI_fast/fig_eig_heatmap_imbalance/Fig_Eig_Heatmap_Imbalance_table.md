# Jacobian occupancy examples

Source: `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\sfaEI_fast\eig_heatmap\eig_heatmap_data.mat`

n (reference) = 500; non-reference examples at n = 250; colour = log10(1 + log10(1 + density)).

| Example | Override | Condition | lambda_1 | mean rate | B_E | median alpha(J_xx) | median omega(J_xx) | eigenvalues |
|---|---|---|---|---|---|---|---|---|
| inhibition-dominant | mu_EE_relative = 4.125 | no_adaptation | +1.9829 | 0.056 | 0.425 | +2.879 | +54.207 | 10000 |
| inhibition-dominant | mu_EE_relative = 4.125 | sfa1_std1 | -0.1892 | 0.050 | 0.425 | +4.332 | +33.965 | 40000 |
| inhibition-dominant | mu_EE_relative = 4.125 | sfa3_std2 | -0.1041 | 0.079 | 0.425 | +1.947 | +26.730 | 80000 |
| reference | none | no_adaptation | +3.9380 | 0.433 | 0.498 | +7.190 | +89.737 | 20000 |
| reference | none | sfa1_std1 | +0.6297 | 0.252 | 0.498 | +6.151 | +79.477 | 80000 |
| reference | none | sfa3_std2 | -0.1121 | 0.217 | 0.498 | +0.273 | +47.022 | 160000 |
| excitation-dominant | mu_EE_relative = 12.38 | no_adaptation | -4.2269 | 0.873 | 0.552 | -2.918 | +40.220 | 10000 |
| excitation-dominant | mu_EE_relative = 12.38 | sfa1_std1 | -0.5258 | 0.707 | 0.552 | +0.239 | +23.110 | 40000 |
| excitation-dominant | mu_EE_relative = 12.38 | sfa3_std2 | -0.1120 | 0.265 | 0.552 | +2.497 | +28.484 | 80000 |
