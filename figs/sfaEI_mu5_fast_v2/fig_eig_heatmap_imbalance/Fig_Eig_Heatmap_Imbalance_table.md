# Jacobian occupancy examples

Source: `C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\sfaEI_mu5_fast\eig_heatmap\eig_heatmap_data.mat`

n (reference) = 500; non-reference examples at n = 250; colour = log10(1 + log10(1 + density)).

| Example | Override | Condition | lambda_1 | mean rate | B_E | median alpha(J_xx) | median omega(J_xx) | eigenvalues |
|---|---|---|---|---|---|---|---|---|
| inhibition-dominant | mu_EE_relative = 2.5 | no_adaptation | +0.6176 | 0.059 | 0.424 | +1.639 | +34.795 | 10000 |
| inhibition-dominant | mu_EE_relative = 2.5 | sfa1_std1 | -0.5180 | 0.059 | 0.424 | +1.492 | +19.337 | 40000 |
| inhibition-dominant | mu_EE_relative = 2.5 | sfa3_std2 | -0.0934 | 0.124 | 0.424 | -2.278 | +12.923 | 80000 |
| reference | none | no_adaptation | +1.7198 | 0.340 | 0.498 | +3.837 | +62.429 | 20000 |
| reference | none | sfa1_std1 | -0.5311 | 0.198 | 0.498 | +2.536 | +57.480 | 80000 |
| reference | none | sfa3_std2 | -0.0978 | 0.208 | 0.498 | -3.231 | +24.195 | 160000 |
| excitation-dominant | mu_EE_relative = 7.5 | no_adaptation | -4.5996 | 0.866 | 0.552 | -3.983 | +21.596 | 10000 |
| excitation-dominant | mu_EE_relative = 7.5 | sfa1_std1 | -0.5221 | 0.591 | 0.552 | -1.620 | +18.798 | 40000 |
| excitation-dominant | mu_EE_relative = 7.5 | sfa3_std2 | -0.0952 | 0.239 | 0.552 | -0.120 | +15.067 | 80000 |
