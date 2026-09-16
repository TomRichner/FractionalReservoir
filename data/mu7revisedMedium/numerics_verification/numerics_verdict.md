# Numerics verification verdict (celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0p025_dualStd_3cond_mu7revisedMedium, medium)

Acceptance criteria fixed on 2026-09-14, before the run. Trials: 5 reshoot, 25 LLE/QR; 3 Jacobian states per condition.

| Condition | Check | Value | Threshold | Result | What |
|---|---|---|---|---|---|
| No Adaptation | A | 2.01 | >= 1 | PASS | noise-free reshoot order slope (median over trials) |
| No Adaptation | B | 1.99 | >= 1 | PASS | noisy reshoot strong-order slope (median; ratio per halving = 2^slope) |
| No Adaptation | L | 0.296 | <= 0.05 and <= scatter 0.6 | FAIL | median |lambda_1 sra1 - lambda_1 ode45| (paired Benettin, full network) |
| No Adaptation | C | [0.00244 0.0023] | <= 0.05 and <= 0.02 | PASS | median |Benettin - QR| and |top-K - QR| lambda_1 (reduced network) |
| No Adaptation | J | [1.05e-10 7.48e-09] | <= 1e-06 and <= 1e-06 | PASS | max rel Frobenius and top-20 eigenvalue error, 0 kink rows excluded |
| Single-Timescale Adaptation | A | 2 | >= 1 | PASS | noise-free reshoot order slope (median over trials) |
| Single-Timescale Adaptation | B | 1.92 | >= 1 | PASS | noisy reshoot strong-order slope (median; ratio per halving = 2^slope) |
| Single-Timescale Adaptation | L | 0.0474 | <= 0.05 and <= scatter 0.161 | PASS | median |lambda_1 sra1 - lambda_1 ode45| (paired Benettin, full network) |
| Single-Timescale Adaptation | C | [0.00654 0.0048] | <= 0.05 and <= 0.02 | PASS | median |Benettin - QR| and |top-K - QR| lambda_1 (reduced network) |
| Single-Timescale Adaptation | J | [4.91e-11 2.93e-10] | <= 1e-06 and <= 1e-06 | PASS | max rel Frobenius and top-20 eigenvalue error, 0 kink rows excluded |
| Multiple-Timescale Adaptation | A | 2 | >= 1 | PASS | noise-free reshoot order slope (median over trials) |
| Multiple-Timescale Adaptation | B | 1.86 | >= 1 | PASS | noisy reshoot strong-order slope (median; ratio per halving = 2^slope) |
| Multiple-Timescale Adaptation | L | 5.97e-07 | <= 0.05 and <= scatter 0.00771 | PASS | median |lambda_1 sra1 - lambda_1 ode45| (paired Benettin, full network) |
| Multiple-Timescale Adaptation | C | [0.00986 0.0125] | <= 0.05 and <= 0.02 | PASS | median |Benettin - QR| and |top-K - QR| lambda_1 (reduced network) |
| Multiple-Timescale Adaptation | J | [4.97e-11 4.96e-12] | <= 1e-06 and <= 1e-06 | PASS | max rel Frobenius and top-20 eigenvalue error, 0 kink rows excluded |

## Manuscript claims

- 400-Hz SRA1 converges to the ode45 reference under step refinement (noise-free): **PASS** (PASS / PASS / PASS)
- SRA1 converges under step refinement on one Brownian path (strong order >= 1): **PASS** (PASS / PASS / PASS)
- the Benettin LLE at 400-Hz SRA1 agrees with ode45-Benettin: **FAIL** (FAIL / PASS / PASS)
- the Benettin LLE agrees with the QR spectrum's lambda_1 and with top-K: **PASS** (PASS / PASS / PASS)
- the analytic Jacobian matches central finite differences: **PASS** (PASS / PASS / PASS)

AT LEAST ONE CRITERION FAILED: see the table. Resolve in code or reflect in the Methods before production analyses.
