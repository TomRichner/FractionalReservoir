# Jacobian occupancy examples

Source: `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\mu7revised_fast\eig_heatmap\eig_heatmap_data.mat`

n (reference) = 500; non-reference examples at n = 250; colour = log10(1 + log10(1 + density)).

| Example | Override | Condition | lambda_1 | mean rate | B_E | median alpha(J_xx) | median omega(J_xx) | eigenvalues |
|---|---|---|---|---|---|---|---|---|
| inhibition-dominant | mu_EE_relative = 3.5 | no_adaptation | +1.5416 | 0.051 | 0.425 | +2.173 | +48.430 | 10000 |
| inhibition-dominant | mu_EE_relative = 3.5 | sfa1_std1 | -0.5000 | 0.048 | 0.425 | +3.472 | +29.334 | 40000 |
| inhibition-dominant | mu_EE_relative = 3.5 | sfa3_std2 | -0.0942 | 0.085 | 0.425 | +0.845 | +22.074 | 80000 |
| reference | none | no_adaptation | +3.0917 | 0.397 | 0.498 | +6.495 | +81.995 | 20000 |
| reference | none | sfa1_std1 | +0.0960 | 0.216 | 0.498 | +5.360 | +75.080 | 80000 |
| reference | none | sfa3_std2 | -0.0980 | 0.213 | 0.498 | -0.978 | +38.423 | 160000 |
| excitation-dominant | mu_EE_relative = 10.5 | no_adaptation | -5.5175 | 0.870 | 0.552 | -4.972 | +30.705 | 10000 |
| excitation-dominant | mu_EE_relative = 10.5 | sfa1_std1 | -0.5201 | 0.676 | 0.552 | -0.744 | +21.155 | 40000 |
| excitation-dominant | mu_EE_relative = 10.5 | sfa3_std2 | -0.0964 | 0.252 | 0.552 | +1.411 | +23.651 | 80000 |
