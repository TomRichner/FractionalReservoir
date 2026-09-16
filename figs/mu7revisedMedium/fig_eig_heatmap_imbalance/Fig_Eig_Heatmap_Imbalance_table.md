# Jacobian occupancy examples

Source: `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\mu7revisedMedium\eig_heatmap\eig_heatmap_data.mat`

n (reference) = 500; non-reference examples at n = 500; colour = log10(1 + log10(1 + density)).

| Example | Override | Condition | lambda_1 | mean rate | B_E | median alpha(J_xx) | median omega(J_xx) | eigenvalues |
|---|---|---|---|---|---|---|---|---|
| inhibition-dominant | mu_EE_relative = 3.5 | no_adaptation | +1.1656 | 0.025 | 0.426 | +2.317 | +59.794 | 30000 |
| inhibition-dominant | mu_EE_relative = 3.5 | sfa1_std1 | -0.5091 | 0.027 | 0.426 | +2.593 | +38.737 | 120000 |
| inhibition-dominant | mu_EE_relative = 3.5 | sfa3_std2 | -0.0921 | 0.040 | 0.426 | +0.983 | +28.362 | 240000 |
| reference | none | no_adaptation | +3.5179 | 0.364 | 0.498 | +5.987 | +81.082 | 75000 |
| reference | none | sfa1_std1 | +0.1591 | 0.214 | 0.498 | +5.267 | +74.504 | 300000 |
| reference | none | sfa3_std2 | -0.0940 | 0.212 | 0.498 | -0.851 | +38.574 | 600000 |
| excitation-dominant | mu_EE_relative = 10.5 | no_adaptation | -4.6448 | 0.911 | 0.553 | -4.435 | +29.787 | 30000 |
| excitation-dominant | mu_EE_relative = 10.5 | sfa1_std1 | -0.5099 | 0.756 | 0.553 | -2.024 | +14.280 | 120000 |
| excitation-dominant | mu_EE_relative = 10.5 | sfa3_std2 | -0.0224 | 0.254 | 0.553 | -1.328 | +18.642 | 240000 |
