# Numerics verification verdict (celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25, fast)

Acceptance criteria fixed on 2026-09-14, before the run. Trials: 1 reshoot, 2 LLE/QR; 3 Jacobian states per condition.

| Condition | Check | Value | Threshold | Result | What |
|---|---|---|---|---|---|
| No Adaptation | A | 2.01 | >= 1 | PASS | noise-free reshoot order slope (median over trials) |
| No Adaptation | B | 2.01 | >= 1 | PASS | noisy reshoot strong-order slope (median; ratio per halving = 2^slope) |
| No Adaptation | L | 0.108 | <= 0.05 | FAIL | median |lambda_1 sra1 - lambda_1 ode45| (paired Benettin, full network) |
| No Adaptation | C | [0.0372 0.00608] | <= 0.05 and <= 0.02 | PASS | median |Benettin - QR| and |top-K - QR| lambda_1 (reduced network) |
| No Adaptation | J | [9.68e-11 7.48e-09] | <= 1e-06 and <= 1e-06 | PASS | max rel Frobenius and top-20 eigenvalue error, 0 kink rows excluded |
| Single-Timescale Adaptation | A | 2.01 | >= 1 | PASS | noise-free reshoot order slope (median over trials) |
| Single-Timescale Adaptation | B | 1.98 | >= 1 | PASS | noisy reshoot strong-order slope (median; ratio per halving = 2^slope) |
| Single-Timescale Adaptation | L | 0.169 | <= 0.05 | FAIL | median |lambda_1 sra1 - lambda_1 ode45| (paired Benettin, full network) |
| Single-Timescale Adaptation | C | [0.308 0.0163] | <= 0.05 and <= 0.02 | FAIL | median |Benettin - QR| and |top-K - QR| lambda_1 (reduced network) |
| Single-Timescale Adaptation | J | [5.01e-11 4.82e-10] | <= 1e-06 and <= 1e-06 | PASS | max rel Frobenius and top-20 eigenvalue error, 0 kink rows excluded |
| Multiple-Timescale Adaptation | A | 2.01 | >= 1 | PASS | noise-free reshoot order slope (median over trials) |
| Multiple-Timescale Adaptation | B | 1.84 | >= 1 | PASS | noisy reshoot strong-order slope (median; ratio per halving = 2^slope) |
| Multiple-Timescale Adaptation | L | 0.00587 | <= 0.05 | PASS | median |lambda_1 sra1 - lambda_1 ode45| (paired Benettin, full network) |
| Multiple-Timescale Adaptation | C | [0.00453 0.0132] | <= 0.05 and <= 0.02 | PASS | median |Benettin - QR| and |top-K - QR| lambda_1 (reduced network) |
| Multiple-Timescale Adaptation | J | [6.04e-11 8.11e-12] | <= 1e-06 and <= 1e-06 | PASS | max rel Frobenius and top-20 eigenvalue error, 0 kink rows excluded |

## Manuscript claims

- 400-Hz SRA1 converges to the ode45 reference under step refinement (noise-free): **PASS** (PASS / PASS / PASS)
- SRA1 converges under step refinement on one Brownian path (strong order >= 1): **PASS** (PASS / PASS / PASS)
- the Benettin LLE at 400-Hz SRA1 agrees with ode45-Benettin: **FAIL** (FAIL / FAIL / PASS)
- the Benettin LLE agrees with the QR spectrum's lambda_1 and with top-K: **FAIL** (PASS / FAIL / PASS)
- the analytic Jacobian matches central finite differences: **PASS** (PASS / PASS / PASS)

AT LEAST ONE CRITERION FAILED: see the table. Resolve in code or reflect in the Methods before production analyses.
