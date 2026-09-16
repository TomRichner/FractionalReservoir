# Numerics verification verdict (celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0p025_dualStd_3cond_mu7revised, fast)

Acceptance criteria fixed on 2026-09-14, before the run. Trials: 1 reshoot, 2 LLE/QR; 3 Jacobian states per condition.

| Condition | Check | Value | Threshold | Result | What |
|---|---|---|---|---|---|
| No Adaptation | A | 2.01 | >= 1 | PASS | noise-free reshoot order slope (median over trials) |
| No Adaptation | B | 1.98 | >= 1 | PASS | noisy reshoot strong-order slope (median; ratio per halving = 2^slope) |
| No Adaptation | L | 0.589 | <= 0.05 | FAIL | median |lambda_1 sra1 - lambda_1 ode45| (paired Benettin, full network) |
| No Adaptation | C | [0.0134 0.00906] | <= 0.05 and <= 0.02 | PASS | median |Benettin - QR| and |top-K - QR| lambda_1 (reduced network) |
| No Adaptation | J | [6.04e-11 1.4e-09] | <= 1e-06 and <= 1e-06 | PASS | max rel Frobenius and top-20 eigenvalue error, 0 kink rows excluded |
| Single-Timescale Adaptation | A | 2 | >= 1 | PASS | noise-free reshoot order slope (median over trials) |
| Single-Timescale Adaptation | B | 1.9 | >= 1 | PASS | noisy reshoot strong-order slope (median; ratio per halving = 2^slope) |
| Single-Timescale Adaptation | L | 0.243 | <= 0.05 | FAIL | median |lambda_1 sra1 - lambda_1 ode45| (paired Benettin, full network) |
| Single-Timescale Adaptation | C | [0.0631 0.00723] | <= 0.05 and <= 0.02 | FAIL | median |Benettin - QR| and |top-K - QR| lambda_1 (reduced network) |
| Single-Timescale Adaptation | J | [4.42e-11 3.38e-10] | <= 1e-06 and <= 1e-06 | PASS | max rel Frobenius and top-20 eigenvalue error, 0 kink rows excluded |
| Multiple-Timescale Adaptation | A | 2.01 | >= 1 | PASS | noise-free reshoot order slope (median over trials) |
| Multiple-Timescale Adaptation | B | 1.86 | >= 1 | PASS | noisy reshoot strong-order slope (median; ratio per halving = 2^slope) |
| Multiple-Timescale Adaptation | L | 0.000502 | <= 0.05 | PASS | median |lambda_1 sra1 - lambda_1 ode45| (paired Benettin, full network) |
| Multiple-Timescale Adaptation | C | [0.0083 0.0146] | <= 0.05 and <= 0.02 | PASS | median |Benettin - QR| and |top-K - QR| lambda_1 (reduced network) |
| Multiple-Timescale Adaptation | J | [4.85e-11 7.11e-12] | <= 1e-06 and <= 1e-06 | PASS | max rel Frobenius and top-20 eigenvalue error, 0 kink rows excluded |

## Manuscript claims

- 400-Hz SRA1 converges to the ode45 reference under step refinement (noise-free): **PASS** (PASS / PASS / PASS)
- SRA1 converges under step refinement on one Brownian path (strong order >= 1): **PASS** (PASS / PASS / PASS)
- the Benettin LLE at 400-Hz SRA1 agrees with ode45-Benettin: **FAIL** (FAIL / FAIL / PASS)
- the Benettin LLE agrees with the QR spectrum's lambda_1 and with top-K: **FAIL** (PASS / FAIL / PASS)
- the analytic Jacobian matches central finite differences: **PASS** (PASS / PASS / PASS)

AT LEAST ONE CRITERION FAILED: see the table. Resolve in code or reflect in the Methods before production analyses.
