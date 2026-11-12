# Tolerance tier: mutation checks, 2026-11-11

FD-020, Plan A step 0c. The tolerance regression tier
(`regression_tol_check`, anchor `scripts/tests/fixtures/regression_tol_anchor.mat`)
gates the rewrite of the model's arithmetic. This note records what the tier
caught when 55 small deliberate errors were put into the model code one at a
time, how far each error exceeded the tolerance that caught it, and what
that says about each tolerance.

- Code under test: commit `d87ea6e`, `src/model/SRNNCellTypePairs.m` at blob
  `2b0118d`, `src/model/SRNN_ESN_reservoir.m` at blob `4b20396`.
- Harness: `scripts/tests/regression/regression_tol_mutate.m`. The list:
  `scripts/tests/regression/regression_tol_mutations_pre_refactor.m`.
- Reference: the anchor. 14 computational threads, as at capture, so the
  unmutated code compares bit for bit (checked before and after the runs:
  `regression_tol_check('anchor', 'rtol', 0)` passes with no leaf moved).
- To repeat: `rows = regression_tol_mutate(regression_tol_mutations_pre_refactor())`,
  about 21 s per mutation, 19 minutes for the list. The mutated code is
  never committed; the harness restores the file with `git checkout` under
  an `onCleanup` and checks the blob after every mutation.

## Result

**All 55 mutations fail the tier, and every one fails by compared values.**
None was caught only because a case threw. No tolerance, override or
comparison rule was changed.

| | Mutations |
|---|---|
| Run | 55 (39 named in the brief, 16 added) |
| Caught | 55 |
| Caught by compared values | 55 |
| Caught only by a thrown case | 0 |
| Also had cases that threw | 1 (`tau_rec_I_to_E_param_1e-6`: two memory-capacity cases, 12 leaves missing, beside 329 value failures) |
| Caught by the pointwise layer | 52 |
| Caught by the stable layer | 48 |
| Caught by both | 45 |
| Caught ONLY by the pointwise layer | 7 |
| Caught ONLY by the stable layer | 3 |
| Caught ONLY by overridden leaves | 2 (both are errors in the reservoir class, not the model class) |
| Made an exact (discrete) leaf change | 12, all of them gross errors (below) |

A "factor" below is a leaf's relative difference divided by the tolerance
that applied to it: above 1 the leaf fails. "Weakest / strongest" are the
failing leaf that exceeded its tolerance least and the one that exceeded it
most, within the layer.

## Every mutation

| Mutation | From | Pointwise: leaves failing | weakest / strongest | Stable: leaves failing (of which overridden) | weakest / strongest | Overrides that failed | Exact leaf changed | Cases that threw |
|---|---|---|---|---|---|---|---|---|
| `tau_rec_1e-6` | brief | 169 | 4.3e3 / 1.8e6 | 263 (1) | 1.1 / 2.4e4 | FD-J (1.1) | no | 0 |
| `tau_rel_1e-6` | brief | 169 | 5e4 / 2.1e6 | 271 (1) | 1.5 / 1.2e4 | FD-J (3.6) | no | 0 |
| `tau_a_rung_1e-6` | brief | 137 | 2.3e4 / 1e6 | 262 (1) | 1 / 2e5 | FD-J (7.7) | no | 0 |
| `c_eff_no_K` | brief | 113 | 2.4e10 / 2e12 | 172 (7) | 1.4e4 / 4.2e9 | FD-J (2.2e6), QR (6.8e4), MC (2.1e5) | yes | 0 |
| `G_1e-6` | brief | 21 | 1.3e5 / 6.2e6 | 42 (0) | 1.4 / 2.4e3 | none | no | 0 |
| `tau_fac_1e-6` | brief | 21 | 6.5e4 / 1.6e6 | 40 (0) | 5.8 / 5.1e2 | none | no | 0 |
| `tau_dec_1e-6` | added | 21 | 4.2e4 / 1.2e6 | 42 (0) | 6.2 / 2.8e3 | none | no | 0 |
| `tau_d_1e-6` | added | 151 | 6.6e5 / 1e6 | 351 (1) | 1 / 3e5 | FD-J (10) | no | 0 |
| `route_scale_not_folded` | brief | 17 | 8.3e10 / 1.6e12 | 42 (0) | 1.1e5 / 2.2e8 | none | no | 0 |
| `tau_rec_I_to_E_param_1e-6` | added | 106 | 4.3e3 / 1e6 | 223 (1) | 1.1 / 2e5 | FD-J (1.1) | no | 2 |
| `recovery_sign` | brief | 14 | 3.9e10 / 3.6e12 | 291 (9) | 5.2e3 / 2.2e11 | FD-J (8.6e6), QR (5.3e4), MC (1e6) | yes | 0 |
| `recovery_I_to_E_1e-6` | brief | 12 | 1.9e4 / 3.4e4 | 217 (1) | 1 / 2.2e5 | FD-J (1.1) | no | 0 |
| `release_1e-6` | added | 14 | 1.5e5 / 2.1e6 | 258 (1) | 1.3 / 2.6e4 | FD-J (3.6) | no | 0 |
| `fac_buildup_sign` | added | 2 | 5.4e11 / 3.2e12 | 42 (0) | 1.1e6 / 2.8e9 | none | no | 0 |
| `sfa_drive_1e-6` | added | 14 | 2.2e5 / 3.1e6 | 261 (1) | 1 / 7.2e4 | FD-J (10) | no | 0 |
| `input_1e-6` | added | 14 | 3.5e5 / 9.2e6 | 346 (1) | 1.2 / 4.1e5 | FD-J (18) | no | 0 |
| `geomean_plain_product_dS` | brief | 2 | 6.1e11 / 1.2e12 | 42 (0) | 1.1e5 / 1.3e8 | none | no | 0 |
| `zero_floor_off_dS` | added | 2 | 8.9e11 / 1.2e12 | 42 (0) | 6.4e4 / 3.8e8 | none | no | 0 |
| `p_min_1e-6` | brief | 21 | 7.6e4 / 1.8e6 | 40 (0) | 1.3 / 8.1e2 | none | no | 0 |
| `geomean_gain_omitted_J` | brief | 11 | 9.5e9 / 1e11 | 0 (0) | none / none | none | no | 0 |
| `geomean_gain_omitted_jt` | added | 2 | 4.6e10 / 7.6e10 | 18 (0) | 5.5e4 / 3.2e7 | none | no | 0 |
| `phi_prime_1e-6` | added | 13 | 2.1e5 / 6.6e5 | 16 (0) | 4.9 / 2e3 | none | no | 0 |
| `j_drop_a_x` | brief | 81 | 9.9e10 / 6.6e11 | 4 (4) | 4.5e3 / 7.1e4 | QR (7.1e4) | yes | 0 |
| `j_drop_a_a_coupling` | brief | 81 | 1.5e10 / 2.8e11 | 4 (4) | 9.8e3 / 5.7e4 | QR (5.7e4) | yes | 0 |
| `j_drop_a_a_diagonal` | added | 81 | 1.4e11 / 6.6e11 | 4 (4) | 4.5e4 / 5.9e4 | QR (5.9e4) | yes | 0 |
| `j_drop_b_a` | brief | 81 | 1.3e10 / 1.6e11 | 4 (4) | 7.4e3 / 5.2e4 | QR (5.2e4) | yes | 0 |
| `j_drop_b_b` | brief | 81 | 1.8e11 / 5.9e11 | 4 (4) | 1.1e4 / 6e4 | QR (6e4) | yes | 0 |
| `j_drop_b_x` | brief | 81 | 8.6e10 / 3.2e11 | 4 (4) | 4e3 / 4.7e4 | QR (4.7e4) | yes | 0 |
| `j_drop_g_a` | brief | 11 | 9.7e9 / 3.3e10 | 0 (0) | none / none | none | no | 0 |
| `j_drop_g_g` | brief | 11 | 1.4e11 / 5.2e11 | 0 (0) | none / none | none | no | 0 |
| `j_drop_g_x` | brief | 11 | 6.5e10 / 2.2e11 | 0 (0) | none / none | none | no | 0 |
| `j_drop_x_a` | brief | 81 | 1.8e10 / 4.6e11 | 4 (4) | 4.1e3 / 6.9e4 | QR (6.9e4) | yes | 0 |
| `j_drop_x_b` | brief | 81 | 7.6e10 / 1e12 | 4 (4) | 4.9e3 / 5.7e4 | QR (5.7e4) | yes | 0 |
| `j_drop_x_g` | brief | 11 | 3.3e10 / 3.8e11 | 0 (0) | none / none | none | no | 0 |
| `j_drop_x_x_rate_block` | brief | 111 | 1.1e11 / 1e12 | 10 (8) | 1.1e4 / 1.8e8 | QR (6.4e4) | yes | 0 |
| `j_drop_x_x_leak` | added | 111 | 3.4e11 / 9.7e11 | 15 (7) | 4.8e4 / 6.2e9 | QR (9.2e4) | yes | 0 |
| `jt_a_x` | brief | 14 | 57 / 3e2 | 12 (0) | 1.4 / 3.4 | none | no | 0 |
| `jt_a_a_coupling` | brief | 14 | 17 / 84 | 2 (0) | 1 / 53 | none | no | 0 |
| `jt_a_a_diagonal` | added | 14 | 58 / 3e2 | 16 (0) | 1.1 / 1.3e2 | none | no | 0 |
| `jt_b_a` | brief | 14 | 15 / 63 | 14 (0) | 1 / 1.5e2 | none | no | 0 |
| `jt_b_b` | brief | 14 | 39 / 5.2e2 | 34 (0) | 1 / 1.8e2 | none | no | 0 |
| `jt_b_x` | brief | 14 | 32 / 2.5e2 | 9 (0) | 1.2 / 17 | none | no | 0 |
| `jt_g_a` | brief | 2 | 13 / 49 | 1 (0) | 1 / 1 | none | no | 0 |
| `jt_g_g` | brief | 2 | 1.1e2 / 4.2e2 | 1 (0) | 4.2 / 4.2 | none | no | 0 |
| `jt_g_x` | brief | 2 | 49 / 2e2 | 0 (0) | none / none | none | no | 0 |
| `jt_x_a` | brief | 14 | 96 / 2.6e2 | 7 (0) | 1.1 / 34 | none | no | 0 |
| `jt_x_b` | brief | 14 | 1.6e2 / 8.7e2 | 25 (0) | 1.2 / 1.6e2 | none | no | 0 |
| `jt_x_g` | brief | 2 | 1.6e2 / 2.7e2 | 0 (0) | none / none | none | no | 0 |
| `jt_x_x_rate_block` | brief | 20 | 3.7e2 / 9.7e2 | 43 (0) | 1 / 1.8e2 | none | no | 0 |
| `jt_x_x_leak` | added | 20 | 87 / 7.5e2 | 57 (0) | 1 / 1.9e2 | none | no | 0 |
| `init_b_1e-6` | brief | 0 | none / none | 269 (0) | 1.2 / 1e5 | none | no | 0 |
| `unpack_no_facilitation` | brief | 6 | 1.1e11 / 4.4e11 | 3 (0) | 2.6e8 / 4.4e8 | none | no | 0 |
| `unpack_no_zero_floor` | brief | 6 | 1.6e11 / 7.4e11 | 3 (0) | 5.5e8 / 6e8 | none | no | 0 |
| `esn_readout_rate_for_synaptic` | added | 0 | none / none | 4 (4) | 4.3e3 / 1.8e5 | MC (1.8e5) | no | 0 |
| `esn_ridge_eta_1e-3` | added | 0 | none / none | 6 (6) | 17 / 1.2e2 | MC (1.2e2) | no | 0 |

"Stable: leaves failing (of which overridden)" counts leaves under one of
the four overrides in the bracket. "Overrides that failed" gives, in
brackets, the largest factor in that group. FD-J is the finite-difference
Jacobian, QR the full-spectrum estimator, MC the memory capacity and its
R2; the fourth override, `orth_defect_max`, failed in no mutation.

### What each mutation is

"From" in the table above: `brief` is an error named in the FD-020 brief
or the issue; `added` was added while running the checks, with the reason
given here.

| Mutation | What it is |
|---|---|
| `tau_rec_1e-6` | every depression tau_rec off by 1e-6 relative |
| `tau_rel_1e-6` | every depression tau_rel off by 1e-6 relative |
| `tau_a_rung_1e-6` | the first SFA rung of every neuron off by 1e-6 relative, in the tau_a_matrix the physics reads (after the loop that fills it, so both the spread and the no-spread branch are covered) |
| `c_eff_no_K` | c_eff: the division by the number of SFA timescales removed |
| `G_1e-6` | every facilitation ceiling G off by 1e-6 relative |
| `tau_fac_1e-6` | every facilitation tau_fac off by 1e-6 relative |
| `tau_dec_1e-6` | every facilitation tau_dec off by 1e-6 relative (the brief names G and tau_fac; this is the third facilitation parameter) |
| `tau_d_1e-6` | tau_d off by 1e-6 relative in params only. The parameter view stores the OBJECT's tau_d, so this is a parameter error the view cannot see: only dS, J and JY can |
| `route_scale_not_folded` | the route scale not folded into params.W |
| `tau_rec_I_to_E_param_1e-6` | the I->E route's tau_rec alone off by 1e-6, as a PARAMETER: the case known to make the reservoir refuse the synaptic readout, kept to record what a thrown case looks like beside value failures |
| `recovery_sign` | dynamics_fast: the sign of the depression recovery term flipped |
| `recovery_I_to_E_1e-6` | dynamics_fast: the recovery term of the I->E route alone scaled by 1 + 1e-6 (one route, in the arithmetic, so no parameter differs between routes and nothing throws) |
| `release_1e-6` | dynamics_fast: the depression release term scaled by 1 + 1e-6 (the brief covers recovery; this is the other half of db/dt) |
| `fac_buildup_sign` | dynamics_fast: the sign of the facilitation build-up term flipped (the facilitation counterpart of recovery_sign) |
| `sfa_drive_1e-6` | dynamics_fast: the rate driving SFA scaled by 1 + 1e-6 (no listed mutation touches da/dt in the arithmetic) |
| `input_1e-6` | dynamics_fast: the external input scaled by 1 + 1e-6 (the term of dx/dt no route or Jacobian mutation reaches) |
| `geomean_plain_product_dS` | dynamics_fast only: the depression factor left as the plain product on a geomean route |
| `zero_floor_off_dS` | dynamics_fast only: the zero floor not applied (the brief lists this for unpack_and_compute_states; the same omission in the right-hand side) |
| `p_min_1e-6` | the zero floor's p_min off by 1e-6 relative |
| `geomean_gain_omitted_J` | compute_Jacobian_fast only: the geomean gain dD/dP omitted from the (x, b) block |
| `geomean_gain_omitted_jt` | jacobian_times only: the geomean gain dD/dP omitted (the mirror of geomean_gain_omitted_J) |
| `phi_prime_1e-6` | the logistic derivative off by 1e-6: every Jacobian block wrong in the same way in J and in jacobian_times, with the right-hand side untouched, so only the stored Jacobians can see it |
| `j_drop_a_x` | compute_Jacobian_fast: the (a, x) block dropped |
| `j_drop_a_a_coupling` | compute_Jacobian_fast: the coupling part of the (a, a) block dropped |
| `j_drop_a_a_diagonal` | compute_Jacobian_fast: the diagonal part of the (a, a) block dropped (the brief lists the coupling part only) |
| `j_drop_b_a` | compute_Jacobian_fast: the (b, a) block dropped |
| `j_drop_b_b` | compute_Jacobian_fast: the (b, b) block dropped |
| `j_drop_b_x` | compute_Jacobian_fast: the (b, x) block dropped |
| `j_drop_g_a` | compute_Jacobian_fast: the (g, a) block dropped |
| `j_drop_g_g` | compute_Jacobian_fast: the (g, g) block dropped |
| `j_drop_g_x` | compute_Jacobian_fast: the (g, x) block dropped |
| `j_drop_x_a` | compute_Jacobian_fast: the (x, a) block dropped |
| `j_drop_x_b` | compute_Jacobian_fast: the (x, b) block dropped |
| `j_drop_x_g` | compute_Jacobian_fast: the (x, g) block dropped |
| `j_drop_x_x_rate_block` | compute_Jacobian_fast: the rate_block of (x, x) dropped (the recurrent part; the leak stays) |
| `j_drop_x_x_leak` | compute_Jacobian_fast: the leak -1/tau_d of (x, x) dropped (the brief lists the rate_block only) |
| `jt_a_x` | jacobian_times: the (a, x) term scaled by 1 + 1e-9 |
| `jt_a_a_coupling` | jacobian_times: the (a, a) coupling term scaled by 1 + 1e-9 |
| `jt_a_a_diagonal` | jacobian_times: the (a, a) diagonal term scaled by 1 + 1e-9 (mirror of j_drop_a_a_diagonal) |
| `jt_b_a` | jacobian_times: the (b, a) term scaled by 1 + 1e-9 |
| `jt_b_b` | jacobian_times: the (b, b) term scaled by 1 + 1e-9 |
| `jt_b_x` | jacobian_times: the (b, x) term scaled by 1 + 1e-9 |
| `jt_g_a` | jacobian_times: the (g, a) term scaled by 1 + 1e-9 |
| `jt_g_g` | jacobian_times: the (g, g) term scaled by 1 + 1e-9 |
| `jt_g_x` | jacobian_times: the (g, x) term scaled by 1 + 1e-9 |
| `jt_x_a` | jacobian_times: the (x, a) term scaled by 1 + 1e-9 |
| `jt_x_b` | jacobian_times: the (x, b) term scaled by 1 + 1e-9 |
| `jt_x_g` | jacobian_times: the (x, g) term scaled by 1 + 1e-9 |
| `jt_x_x_rate_block` | jacobian_times: the recurrent (x, x) term scaled by 1 + 1e-9 |
| `jt_x_x_leak` | jacobian_times: the leak term of (x, x) scaled by 1 + 1e-9 (mirror of j_drop_x_x_leak) |
| `init_b_1e-6` | initialize_state: depression starting at 1 - 1e-6 instead of 1 |
| `unpack_no_facilitation` | unpack_and_compute_states: the facilitation factor dropped from the synaptic output |
| `unpack_no_zero_floor` | unpack_and_compute_states: the zero floor not applied |
| `esn_readout_rate_for_synaptic` | SRNN_ESN_reservoir: the synaptic readout returns the rate (depression left out of the readout). FD-006 rewrites this readout, and only the memory-capacity leaves, under their 1e-6 override, can see it |
| `esn_ridge_eta_1e-3` | SRNN_ESN_reservoir: the ridge regulariser off by 1e-3 relative. A PROBE of how fine an error in the reservoir's own arithmetic the 1e-6 override on MC and R2 still sees, not an error of the listed kind |

Sites that did not exist as the brief described them, and what was used:

- **"An SFA `tau_a` rung off by 1e-6 (in the `tau_a_matrix` the physics
  reads)".** `get_params` fills `params.tau_a_matrix` in two branches (the
  per-neuron draw, or the nominal row replicated). A one-branch edit would
  miss half the cases, so the mutation scales column 1 after the loop.
- **"`D` left as the plain product in `dynamics_fast` only".** `dynamics_fast`
  does not form `D`; it calls `combine_std_pair`. The mutation replaces that
  one call by `prod(b, 2)`.
- **"The matching term for each block scaled by `1 + 1e-9`" in
  `jacobian_times`.** The function has no term per block: the `(., a)` and
  `(., x)` blocks share one perturbation `Ye = Yx - c_eff * sum(Ya)`. Each
  mutation adds `1e-9` times the one part of `Ye` that belongs to its block,
  at the one site that block uses, which is the same thing to round-off.
- **"`(a, a)` coupling" and "the `rate_block` in `(x, x)`".** Both blocks
  have a second part (the diagonal `-1/tau`, the leak `-1/tau_d`). Those
  were added as mutations of their own, in the assembled Jacobian and in
  `jacobian_times`.

### Which layer caught what

**Caught only by the pointwise layer (7).** `geomean_gain_omitted_J`,
`j_drop_g_a`, `j_drop_g_g`, `j_drop_g_x`, `j_drop_x_g`, `jt_g_x`, `jt_x_g`.
The first five are errors in the ASSEMBLED Jacobian of a network with
facilitation or the geomean combine. The stable layer reaches the assembled
Jacobian only through the QR estimator and the frozen transient gain, and
both run on the paper's network, which has neither option; the options
network's stable run uses top-K, which goes through `jacobian_times`. The
pointwise layer catches all five by factors of 1e10 and more, so nothing
slips through, but **for the facilitation and geomean blocks of
`compute_Jacobian_fast` the stored pointwise `J` is the tier's only
witness.** The last two are 1e-9 scalings of `jacobian_times` terms whose
effect on the options network's exponents is below 1e-9.

**The stable layer sees the other assembled-Jacobian errors only through
QR.** The eight dropped blocks in the SFA, depression and `x` rows and
columns (`j_drop_a_x` to `j_drop_x_b`) fail four stable leaves each, all of
them under the QR override (three QR floats and the QR `sort_idx`). The two
`(x, x)` mutations also fail the frozen transient gain, under the layer's
own tolerance.

**Caught only by the stable layer (3).** `init_b_1e-6` (269 leaves under
the layer's own tolerance; `initialize_state` is not evaluated in the
pointwise layer, by design) and the two reservoir mutations.

**Caught only by overridden leaves (2).** `esn_readout_rate_for_synaptic`
and `esn_ridge_eta_1e-3`, both by the `MC` and `R2` leaves. This is by
construction: memory capacity is the only place the tier runs the
reservoir class's own arithmetic. No mutation of the MODEL class was
caught only by an overridden leaf.

**Exact leaves.** Twelve mutations changed a discrete leaf: `c_eff_no_K`,
`recovery_sign` and the ten dropped Jacobian blocks that QR sees. The
leaves were `sort_idx` (the order of the exponents) in every case, and
`transient_gain.frozen_x.info.i_peak` for the dropped leak. **No mutation
of size 1e-6 or 1e-9 flipped any discrete leaf**: the counts of positive
exponents, `K_used`, `retries`, the clip counts and the sort orders all
held. That is the margin Plan A, section 8, asked this note to record for
"discrete outputs can flip under round-off": errors nine orders of
magnitude above round-off did not move them.

**The thrown case.** `tau_rec_I_to_E_param_1e-6` makes the I->E and I->I
routes differ as parameters, so `SRNN_ESN_reservoir.assert_route_redundancy`
refuses the synaptic readout in the two adapting memory-capacity cases.
The tier reported both as ERROR with their 12 leaves missing and went on;
the same mutation fails 106 pointwise and 223 stable leaves by value. The
arithmetic version of a one-route error (`recovery_I_to_E_1e-6`) changes no
parameter, throws nothing, and fails 12 pointwise leaves (all `dS`, by
1.9e4 to 3.4e4) and 217 stable leaves.

## Detection-limit scan

Three errors repeated at five sizes. Each cell is the largest factor among
the leaves of that kind; **bold** means at least one leaf of that kind
failed. `view` is the canonical parameter view, `rate, theta` the rates and
synaptic outputs of `unpack_and_compute_states`; the first five columns
are the pointwise layer, the rest the stable layer, the last three its
overrides.

| Error | Size | view | dS | J | JY | rate, theta | trajectories | top-K | gain | Benettin | sweep | FD-J | QR | MC | leaves failing: pointwise / stable |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| tau_rec | 1e-4 | **1e8** | **1.8e8** | **3.3e7** | **2.7e7** | **7.2e7** | **3.6e4** | **3.8e5** | **1.2e4** | **1.2e4** | **2.4e6** | **1.1e2** | **2.4** | **53** | 169 / 282 |
| tau_rec | 1e-6 | **1e6** | **1.8e6** | **3.3e5** | **2.7e5** | **7.2e5** | **3.6e2** | **3.8e3** | **1.2e2** | **1.2e2** | **2.4e4** | **1.1** | 0.024 | 0.53 | 169 / 263 |
| tau_rec | 1e-8 | **1e4** | **1.8e4** | **3.3e3** | **2.7e3** | **7.2e3** | **3.6** | **38** | **1.2** | **1.2** | **2.4e2** | 0.011 | 0.0055 | 0.0054 | 169 / 181 |
| tau_rec | 1e-10 | **1e2** | **1.8e2** | **33** | **27** | **72** | 0.036 | 0.38 | 0.012 | 0.012 | **2.4** | 0.0017 | 0.0058 | 0.00087 | 161 / 1 |
| tau_rec | 1e-12 | **1** | **1.8** | 0.33 | 0.27 | 0.72 | 0.00036 | 0.0038 | 0.00012 | 0.0063 | 0.019 | 0.0017 | 0.006 | 0.00033 | 56 / 0 |
| tau_a | 1e-4 | **1e8** | **1e8** | **8.4e7** | **5.6e7** | 0 | **7.4e4** | **2.3e5** | **6.6e2** | **3.5e4** | **2e7** | **7.7e2** | **14** | **86** | 137 / 283 |
| tau_a | 1e-6 | **1e6** | **1e6** | **8.4e5** | **5.6e5** | 0 | **7.4e2** | **2.3e3** | **6.6** | **3.5e2** | **2e5** | **7.7** | 0.14 | 0.86 | 137 / 262 |
| tau_a | 1e-8 | **1e4** | **1e4** | **8.4e3** | **5.6e3** | 0 | **7.4** | **23** | 0.066 | **3.5** | **2e3** | 0.077 | 0.0086 | 0.0084 | 137 / 62 |
| tau_a | 1e-10 | **1e2** | **1e2** | **84** | **56** | 0 | 0.074 | 0.23 | 0.00066 | 0.034 | **20** | 0.00086 | 0.0059 | 0.00045 | 137 / 1 |
| tau_a | 1e-12 | **1** | **1** | 0.84 | 0.56 | 0 | 0.00074 | 0.0023 | 6.5e-06 | 0.0079 | 0.2 | 0.0017 | 0.006 | 0.00063 | 15 / 0 |
| jt_x_b | 1e-4 | 0 | 0 | 0 | **8.7e7** | 0 | 0 | **9.9e5** | **4.9e3** | 0 | **1.6e7** | 0 | 0 | 0 | 14 / 170 |
| jt_x_b | 1e-6 | 0 | 0 | 0 | **8.7e5** | 0 | 0 | **9.9e3** | **49** | 0 | **1.6e5** | 0 | 0 | 0 | 14 / 155 |
| jt_x_b | 1e-8 | 0 | 0 | 0 | **8.7e3** | 0 | 0 | **99** | 0.49 | 0 | **1.6e3** | 0 | 0 | 0 | 14 / 59 |
| jt_x_b | 1e-10 | 0 | 0 | 0 | **87** | 0 | 0 | 0.99 | 0.0049 | 0 | **16** | 0 | 0 | 0 | 14 / 1 |
| jt_x_b | 1e-12 | 0 | 0 | 0 | 0.87 | 0 | 0 | 0.0099 | 4.9e-05 | 0 | 0.16 | 0 | 0 | 0 | 0 / 0 |

The smallest size each kind of leaf still catches:

| Leaves | `tau_rec` | `tau_a` rung | `jacobian_times` `(x, b)` |
|---|---|---|---|
| Pointwise: parameter view | 1e-12 (see below) | 1e-12 (see below) | not applicable |
| Pointwise: `dS` | 1e-12 (factor 1.8) | 1e-12 (factor 1.0) | not applicable |
| Pointwise: `J`, `JY` | 1e-10 | 1e-10 | `JY`: 1e-10; missed at 1e-12 (factor 0.87) |
| Stable: trajectories | 1e-8 | 1e-8 | not applicable |
| Stable: top-K results | 1e-8 | 1e-8 | 1e-8 (0.99 at 1e-10) |
| Stable: sweep | 1e-10 (one leaf) | 1e-10 (one leaf) | 1e-10 (one leaf) |
| Stable: transient gain | 1e-8 | 1e-6 | 1e-6 |
| Stable: Benettin | 1e-8 | 1e-8 | not applicable |
| Override: finite-difference `J` (1e-7) | 1e-6 | 1e-6 | not applicable |
| Override: QR (1e-5) | 1e-4 | 1e-4 | not applicable |
| Override: MC, R2 (1e-6) | 1e-4 | 1e-4 (0.86 at 1e-6) | not applicable |

Reading it:

- **The pointwise layer resolves a relative error of about 1e-12 in a
  parameter and about 1e-11 to 1e-12 in a single term of `jacobian_times`.**
  Factors scale linearly with the size of the error all the way down, which
  is what a clean detector looks like: nothing saturates and nothing is lost
  in a floor.
- The catches at 1e-12 are at the edge and partly an accident of
  representation: `1 + 1e-12` in double precision is `1 + 1.00009e-12`, so
  the stored parameter moves by 1.00009e-12 relative and fails a 1e-12 test
  by a factor of 1.00009. Treat 1e-12 as "at the tolerance", not as caught.
- **The stable layer resolves about 1e-8 across the board** (trajectories,
  top-K, Benettin) and 1e-10 in one leaf, `stable.sweep.sfa1_std1{2}.lambda_1_drift`.
  That leaf is the difference of two nearly equal exponents, so it amplifies
  a relative error by 20 to 200 times; it is also the stable leaf whose
  round-off floor is highest (1.3e-12 under a harmless reordering).
- The three looser overrides resolve 1e-6 (finite-difference Jacobian) and
  1e-4 (QR, memory capacity). They are two to four orders of magnitude
  blunter than the leaves beside them, which is why no model-class mutation
  depended on them.

## Margins

For each tolerance: the measured round-off floor (from the header of
`regression_tol_settings.m`: a harmless reordering of five expressions,
and thread counts 1 to 14), and the weakest catch among the 55 mutations.
"Weakest catch" is, over the mutations in which that group of leaves
failed, the smallest value of the mutation's LARGEST factor in the group:
how far above the tolerance the least visible error still stood.

| Tolerance | Value | Round-off floor (as a factor of the tolerance) | Mutations it caught | Weakest catch | Smallest factor on any failing leaf | Largest factor that still passed |
|---|---|---|---|---|---|---|
| Pointwise layer | 1e-12 | 6.0e-16 on `dS` (6e-4) | 52 of 55 | 49 (`jt_g_a`, a 1e-9 error) | 13 (`jt_g_a`, `options_chaotic.JY`) | none moved and passed |
| of which `dS` | | 6.0e-16 | 19 | 3.4e4 (a 1e-6 error on one route) | | |
| of which `J` | | 2.3e-16 | 27 | 4.8e4 (a 1e-6 error on one route) | | |
| of which `JY` | | 3.4e-16 | 27 | 49 (`jt_g_a`) | | |
| Stable layer, own tolerance | 1e-9 | 2.0e-12, Benettin local rates (2e-3) | 38 of 55 | 1.03 (`jt_g_a`); 510 among the 1e-6 errors (`tau_fac_1e-6`) | 1.0 | 0.98 |
| Override: finite-difference `J` | 1e-7 | 1.7e-10 (1.7e-3) | 11 | 1.1 | 1.1 | 0.22 |
| Override: QR | 1e-5 | 6.1e-8 (6e-3) | 12 | 4.7e4 (`j_drop_b_x`) | 4.0e3 | 0.23 |
| Override: MC, R2 | 1e-6 | 8.9e-10 reordering, 6e-9 thread count (6e-3) | 4 | 120 (`esn_ridge_eta_1e-3`) | 17 | 0.94 (`input_1e-6`) |
| Override: `orth_defect_max` | 1e3 | moves by tens of percent of itself | 0 | never failed | | 4.5e-4 |

The stable layer as a whole caught 48 mutations; 38 of them through leaves
under its own tolerance and 10 only through QR or MC leaves (the eight
dropped Jacobian blocks and the two reservoir errors).

So the tolerances sit between a floor and a catch like this:

- **Pointwise, 1e-12**: 1 600 times above the worst harmless reordering,
  49 times below the least visible of the listed errors (a 1e-9
  disagreement in the smallest block of `jacobian_times`), and 3e4 times
  or more below every 1e-6 error.
- **Stable, 1e-9**: 500 times above its floor, 500 times or more below
  every 1e-6 error in the dynamics. It does NOT have a margin on the 1e-9
  `jacobian_times` errors (factors 1 to 190, two not seen at all), and it
  cannot: an error of 1e-9 in one term moves an exponent by at most about
  1e-9. Those are the pointwise layer's to catch, and it does.

## Recommendation for each tolerance

| Tolerance | Recommendation | Reason |
|---|---|---|
| `rtol.pointwise = 1e-12` | **Keep.** | Every one of the 52 mutations that reach a pointwise quantity fails it, the weakest by 49 times, and the floor is 1 600 times below. Tightening to 1e-13 would still clear the measured floor by 170 times and would catch a 1e-12 disagreement in `jacobian_times` (missed here at factor 0.87), but the floor was measured on five reorderings and the refactor of step 0d will reorder far more than five expressions; the margin is better spent there. If the refactor's own before/after comparison shows nothing above about 1e-14, tightening afterwards is cheap and safe. |
| `floor_frac = 1` | **Keep.** | The worry was that judging a small Jacobian block against the whole leaf's scale would hide an error in it. It did not: the smallest block scaled by `1 + 1e-9` (`jt_g_a`) still failed by 13 to 49 times, and every dropped block failed by 1e10 or more. A block-by-block rule would buy about one more order of magnitude on the smallest blocks and would cost a comparison rule that has to know the state layout. |
| `rtol.stable = 1e-9` | **Keep.** | 500 times above the floor; every 1e-6 error in the dynamics fails by 500 times or more; no discrete leaf flipped under any small error. Tightening to 1e-10 would leave only 50 times above a floor measured on a handful of reorderings, for a layer whose job is the end-to-end check, not resolution. Loosening is not indicated either: 1e-8 is exactly where the trajectories stop seeing a parameter error. |
| Override, finite-difference `J`, 1e-7 | **Keep.** | It is redundant by construction: the leaf changes only when `dS` changes, and `dS` is compared directly at 1e-12. It failed in 11 mutations and was never the only failure. It cannot be tightened usefully (the floor is 1.7e-10) and loosening it would lose nothing; it is kept because it costs nothing and its failure names the right function. |
| Override, QR, 1e-5 | **Keep.** | 160 times above its floor, so no room to tighten. It matters more than its looseness suggests: inside the stable layer it is the only leaf group that sees a dropped block of the assembled Jacobian, and it sees each by 4e3 times or more. |
| Override, `MC` and `R2`, 1e-6 | **Keep.** | 170 times above the thread-count floor. It is the tier's only guard on the reservoir class: a readout that returns the rate instead of the synaptic output fails by 1.8e5, a ridge regulariser off by 1e-3 by 120. A 1e-6 error in the model's dynamics moves `MC` by up to 0.94 of the tolerance and passes; that is acceptable only because the same error fails hundreds of other leaves. |
| Override, `orth_defect_max`, 1e3 | **Keep, as a sanity bound.** | No mutation moved it by more than 4.5e-4 of its tolerance. It is not a detector of model arithmetic and no mutation in this list targets the re-orthonormalisation in `lyapunov_topk`; the list was confined to the two model classes. |
| The discrete list | **Keep exact.** | See "Exact leaves" above. |

## What slipped through, and what is thin

Nothing in the list passed the tier. Three things are thin and should be
known before the refactor leans on the tier:

1. **The assembled Jacobian's facilitation and geomean blocks have one
   witness.** The stable layer never assembles `compute_Jacobian_fast` on a
   network with facilitation, dual facilitation, the geomean combine, the
   zero floor or a route scale. Five mutations there were caught by the
   stored pointwise `J` alone (by 1e10, so not marginally). The chain routes
   of FD-004 will be in the same position unless a case exercises them. A
   second witness would be a QR run or a frozen transient gain on the stable
   options network; **that needs new stored values and so a new anchor, and
   was not done.** It is a question for the lead, not a defect in the tier.
2. **The stable layer adds nothing on `jacobian_times` errors of 1e-9.**
   It catches them at factors near 1 or not at all. The contract
   "`jacobian_times` equals the assembled Jacobian times Y" rests on the
   pointwise `JY` leaves in this tier, at about 1e-11 resolution, and on
   `test_jacobian_times` outside it.
3. **Parameter errors are caught first by the parameter view**, which
   compares the stored parameter itself. That is a real check (a refactor
   that compiles a route's parameters wrongly is exactly what it is for),
   but it says nothing about the arithmetic. The evidence that the
   arithmetic is guarded is in the non-view columns: `tau_d_1e-6`, which the
   view cannot see, fails `dS`, `J` and `JY` by 1e6 each, and the scan shows
   `dS`, `J` and `JY` following a parameter error down to 1e-10 without help
   from the view.

Not tested here: errors in `lyapunov_topk`, in the integrators, and in
`ParamSpaceAnalysis2`. The tier runs all three, but the list was confined
to the model class and the reservoir class.

## How the restore was checked

After every mutation the harness compares `git hash-object <file>` with the
blob in `HEAD` and requires `git status --porcelain -- <file>` to be empty;
it raises an error otherwise, and every one of the 70 runs (55 and the 15
of the scan) reported the blob restored. After the last run
`git status --porcelain -- src/` was empty, `SRNNCellTypePairs.m` hashed to
`2b0118d42edca1ca7cc62e4227220ca9cb5a58db`, and
`regression_tol_check('anchor', 'rtol', 0)` passed with no leaf moved,
which also shows the restored code, not a cached mutated class, was what
MATLAB ran.

## For the list that follows the refactor

The texts in `regression_tol_mutations_pre_refactor` are tied to the code
as it stood; after FD-021 most will not be found, and the harness then
refuses the entry instead of running a mutation that changes nothing. The
errors to carry over are the "What it is" column above, at their new
sites. Two things to keep when rewriting it: one mutation per Jacobian
block in each of the two Jacobian routines, and at least one error that
the parameter view cannot see.

---

# After the refactor: the list re-run on the route interface (2026-11-11, FD-021)

Plan A step 0d moved product depression and today's facilitation out of
`SRNNCellTypePairs` and behind the route interface (`src/model/routes/`:
`route_model`, `route_product`). This section records the same exercise
on the refactored code: the list written again for it, run with the same
harness against the same anchor, with no tolerance, override or
comparison rule changed.

- Code under test: commit `3857fe9`. `src/model/SRNNCellTypePairs.m` at
  blob `8c94d21`, `src/model/routes/route_product.m` at `b3f5b2e`,
  `src/model/routes/route_model.m` at `dc18579`,
  `src/model/SRNN_ESN_reservoir.m` at `4b20396` (unchanged).
- The list: `scripts/tests/regression/regression_tol_mutations_post_refactor.m`,
  71 mutations (48 carried over from the first list at their new sites, 23
  new) and a 15-entry scan. Its header says what moved where.
- Reference: the anchor, 14 threads. The unmutated refactored code passes
  it, but no longer bit for bit: 240 of 1749 float leaves move (64
  pointwise, 176 stable). The leaf closest to its tolerance is
  `pointwise.options_chaotic.dS` at 4.8e-16 against 1e-12 (2100 times
  inside) in the pointwise layer; in the stable layer
  `stable.estimators.sfa3_std2.qr.local_LE_spectrum_t` at 7.6e-8 against
  its 1e-5 override (130 times inside), and among leaves at the layer's
  own 1e-9, `stable.sweep.sfa1_std1{2}.lambda_1_drift` at 1.3e-12 (770
  times inside). So a mutation's factors below sit on a background of at
  most 5e-4 (pointwise) and 1.3e-3 (stable, default tolerance).
- To repeat: `rows = regression_tol_mutate(regression_tol_mutations_post_refactor())`,
  about 22 s per mutation, 26 minutes for the list: longer than is safe
  for one MATLAB MCP call, so it was run in three parts with `save_dir`.

## Result of the re-run

**All 71 mutations fail the tier, and every one fails by compared
values.** None was caught only because a case threw.

| | Mutations |
|---|---|
| Run | 71 (48 carried, 23 new) |
| Caught | 71 |
| Caught by compared values | 71 |
| Caught only by a thrown case | 0 |
| Also had cases that threw | 1 (`tau_rec_I_to_E_param_1e-6`: two memory-capacity cases, 12 leaves missing, beside 329 value failures; as before the refactor) |
| Caught by the pointwise layer | 66 |
| Caught by the stable layer | 68 |
| Caught by both | 63 |
| Caught ONLY by the pointwise layer | 3 (`lin_A_spurious_coupling_1e-6`, `lin_Fs_g_1e-9`, `lin_Fs_b_with_g_1e-9`) |
| Caught ONLY by the stable layer | 5 (`init_b_1e-6`, `init_class_1e-6`, `unpack_b_reversed`, and the two reservoir mutations) |
| Caught ONLY by overridden leaves | 2 (the two reservoir mutations, as before) |
| Made an exact (discrete) leaf change | 16, all gross errors |

The numbers of the carried mutations are the first list's to the digits
shown wherever the error is the same error (`tau_rec_1e-6`: 169 pointwise
leaves, 4.3e3 / 1.8e6; 263 stable, 1.1 / 2.4e4, in both runs). The tier
sees the same error the same way on the refactored code.

## Every mutation of the re-run

Columns as in the first table; the stable column gives the failing
leaves and, in brackets, how many of them are under an override.

| Mutation | From | Pointwise: leaves failing | weakest / strongest | Stable: leaves failing (overridden) | weakest / strongest | Exact leaf changed | Cases that threw |
|---|---|---|---|---|---|---|---|
| `tau_rec_1e-6` | carried | 169 | 4.3e3 / 1.8e6 | 263 (1) | 1.1 / 2.4e4 | no | 0 |
| `tau_rel_1e-6` | carried | 169 | 5e4 / 2.1e6 | 271 (1) | 1.5 / 1.2e4 | no | 0 |
| `tau_a_rung_1e-6` | carried | 137 | 2.3e4 / 1e6 | 262 (1) | 1 / 2e5 | no | 0 |
| `c_eff_no_K` | carried | 113 | 2.4e10 / 2e12 | 172 (7) | 1.4e4 / 4.2e9 | yes | 0 |
| `G_1e-6` | carried | 21 | 1.3e5 / 6.2e6 | 42 (0) | 1.4 / 2.4e3 | no | 0 |
| `tau_fac_1e-6` | carried | 21 | 6.5e4 / 1.6e6 | 40 (0) | 5.8 / 5.1e2 | no | 0 |
| `tau_dec_1e-6` | carried | 21 | 4.2e4 / 1.2e6 | 42 (0) | 6.2 / 2.8e3 | no | 0 |
| `tau_d_1e-6` | carried | 151 | 6.6e5 / 1e6 | 351 (1) | 1 / 3e5 | no | 0 |
| `route_scale_not_folded` | carried | 17 | 8.3e10 / 1.6e12 | 42 (0) | 1.1e5 / 2.2e8 | no | 0 |
| `tau_rec_I_to_E_param_1e-6` | carried | 106 | 4.3e3 / 1e6 | 223 (1) | 1.1 / 2e5 | no | 2 |
| `P_tau_rest_1e-6` | new | 115 | 4.3e3 / 1.8e6 | 263 (1) | 1.1 / 2.4e4 | no | 0 |
| `P_tau_drive_1e-6` | new | 115 | 5e4 / 2.1e6 | 271 (1) | 1.5 / 1.2e4 | no | 0 |
| `P_target_1e-6` | new | 15 | 1.3e5 / 6.2e6 | 42 (0) | 1.4 / 2.4e3 | no | 0 |
| `P_one_neuron_tau_rest_1e-6` | new | 115 | 1.6e3 / 1.5e6 | 246 (0) | 1 / 8.7e4 | no | 0 |
| `recovery_sign` | carried | 14 | 3.9e10 / 3.6e12 | 298 (9) | 5.2e3 / 1.1e12 | yes | 0 |
| `recovery_I_to_E_1e-6` | carried | 12 | 1.2e5 / 1.9e5 | 218 (1) | 1 / 4e4 | no | 0 |
| `release_1e-6` | carried | 14 | 1.5e5 / 2.1e6 | 260 (1) | 1.3 / 2.6e4 | no | 0 |
| `fac_buildup_sign` | carried | 2 | 5.4e11 / 3.2e12 | 42 (0) | 1.1e6 / 2.8e9 | no | 0 |
| `sfa_drive_1e-6` | carried | 14 | 2.2e5 / 3.1e6 | 261 (1) | 1 / 7.2e4 | no | 0 |
| `input_1e-6` | carried | 14 | 3.5e5 / 9.2e6 | 346 (1) | 1.2 / 4.1e5 | no | 0 |
| `theta_1e-6` | new | 20 | 4.5e5 / 4.6e6 | 345 (2) | 1.2 / 3.9e5 | no | 0 |
| `geomean_plain_product_dS` | carried | 2 | 6.1e11 / 1.2e12 | 42 (0) | 1.1e5 / 1.3e8 | no | 0 |
| `zero_floor_off_dS` | carried | 2 | 8.9e11 / 1.2e12 | 42 (0) | 6.4e4 / 3.8e8 | no | 0 |
| `output_no_facilitation` | carried | 8 | 1.1e11 / 1.5e12 | 42 (0) | 1.1e5 / 4.4e8 | no | 0 |
| `p_min_1e-6` | carried | 21 | 7.6e4 / 1.8e6 | 40 (0) | 1.3 / 8.1e2 | no | 0 |
| `geomean_gain_omitted` | carried | 13 | 9.5e9 / 1e11 | 18 (0) | 5.5e4 / 3.2e7 | no | 0 |
| `floor_gain_omitted` | new | 13 | 4.3e10 / 3.3e11 | 18 (0) | 5.8e4 / 3.5e7 | no | 0 |
| `geomean_exponent_1e-6` | new | 17 | 2.5e4 / 5.7e5 | 40 (0) | 1.3 / 3.8e2 | no | 0 |
| `phi_prime_1e-6` | carried | 13 | 2.1e5 / 6.6e5 | 16 (0) | 4.9 / 2e3 | no | 0 |
| `lin_A_last_1e-9` | new | 95 | 20 / 4.1e2 | 33 (0) | 1 / 1.8e2 | no | 0 |
| `lin_A_first_dropped` | new | 95 | 3.9e10 / 5.9e11 | 297 (4) | 6.7e3 / 1.4e11 | yes | 0 |
| `lin_A_spurious_coupling_1e-6` | new | 45 | 3.4e4 / 9.7e4 | 0 (0) | - | no | 0 |
| `lin_B_1e-9` | new | 95 | 47 / 3.2e2 | 19 (0) | 1 / 1.7e2 | no | 0 |
| `lin_B_facilitation_dropped` | carried | 13 | 4.8e10 / 2.2e11 | 18 (0) | 1.1e5 / 5.1e8 | no | 0 |
| `lin_F_1e-9` | new | 131 | 1.1e2 / 1e3 | 47 (0) | 1 / 1.8e2 | no | 0 |
| `lin_Fs_1e-9` | new | 95 | 76 / 1e3 | 25 (0) | 1.2 / 1.6e2 | no | 0 |
| `lin_Fs_g_1e-9` | carried | 13 | 33 / 3.8e2 | 0 (0) | - | no | 0 |
| `lin_Fs_b_with_g_1e-9` | new | 13 | 51 / 6e2 | 0 (0) | - | no | 0 |
| `lin_others_not_swapped` | new | 53 | 1.2e9 / 3.4e11 | 100 (4) | 5.8e2 / 2.9e8 | yes | 0 |
| `j_drop_a_x` | carried | 81 | 9.9e10 / 6.6e11 | 4 (4) | 4.5e3 / 7.1e4 | yes | 0 |
| `j_drop_a_a_coupling` | carried | 81 | 1.5e10 / 2.8e11 | 4 (4) | 9.8e3 / 5.7e4 | yes | 0 |
| `j_drop_a_a_diagonal` | carried | 81 | 1.4e11 / 6.6e11 | 4 (4) | 4.5e4 / 5.9e4 | yes | 0 |
| `j_drop_s_a` | carried | 81 | 1.3e10 / 1.6e11 | 4 (4) | 7.4e3 / 5.2e4 | yes | 0 |
| `j_drop_s_s` | carried | 81 | 1.8e11 / 5.9e11 | 4 (4) | 1.1e4 / 6e4 | yes | 0 |
| `j_drop_s_x` | carried | 81 | 8.6e10 / 3.2e11 | 4 (4) | 4e3 / 4.7e4 | yes | 0 |
| `j_drop_x_a` | carried | 81 | 1.8e10 / 4.6e11 | 4 (4) | 4.1e3 / 6.9e4 | yes | 0 |
| `j_drop_x_s` | carried | 81 | 7.6e10 / 1e12 | 4 (4) | 4.9e3 / 5.7e4 | yes | 0 |
| `j_drop_x_x_rate_block` | carried | 111 | 1.1e11 / 1e12 | 10 (8) | 1.1e4 / 1.8e8 | yes | 0 |
| `j_drop_x_x_leak` | carried | 111 | 3.4e11 / 9.7e11 | 15 (7) | 4.8e4 / 6.2e9 | yes | 0 |
| `jt_a_x` | carried | 14 | 57 / 3e2 | 12 (0) | 1.4 / 3.4 | no | 0 |
| `jt_a_a_coupling` | carried | 14 | 17 / 84 | 2 (0) | 1 / 53 | no | 0 |
| `jt_a_a_diagonal` | carried | 14 | 58 / 3e2 | 16 (0) | 1.1 / 1.3e2 | no | 0 |
| `jt_s_a` | carried | 14 | 15 / 63 | 15 (0) | 1 / 1.5e2 | no | 0 |
| `jt_s_s` | carried | 14 | 39 / 5.2e2 | 35 (0) | 1 / 1.8e2 | no | 0 |
| `jt_s_s_last_variable` | new | 14 | 20 / 3.7e2 | 33 (0) | 1 / 1.8e2 | no | 0 |
| `jt_s_x` | carried | 14 | 32 / 2.5e2 | 9 (0) | 1.2 / 17 | no | 0 |
| `jt_x_a` | carried | 14 | 96 / 2.6e2 | 7 (0) | 1.1 / 34 | no | 0 |
| `jt_x_s` | carried | 14 | 1.6e2 / 9.9e2 | 25 (0) | 1.2 / 1.6e2 | no | 0 |
| `jt_x_s_last_variable` | new | 14 | 1.6e2 / 9.9e2 | 25 (0) | 1.2 / 1.6e2 | no | 0 |
| `jt_x_x_rate_block` | carried | 20 | 3.7e2 / 9.7e2 | 43 (0) | 1 / 1.8e2 | no | 0 |
| `jt_x_x_leak` | carried | 20 | 87 / 7.5e2 | 57 (0) | 1 / 1.9e2 | no | 0 |
| `route_idx_reversed` | new | 167 | 7.1e10 / 4.4e12 | 183 (5) | 3.9e3 / 6.5e10 | yes | 0 |
| `route_idx_blocks_swapped` | new | 17 | 1.5e11 / 5.1e12 | 19 (0) | 2.4e6 / 2.4e9 | no | 0 |
| `init_b_1e-6` | carried | 0 | - | 269 (0) | 1.2 / 1e5 | no | 0 |
| `init_class_1e-6` | new | 0 | - | 269 (0) | 1.2 / 1e5 | no | 0 |
| `registry_wrong_output` | new | 58 | 1.1e11 / 2.9e12 | 33 (4) | 4.3e3 / 3.3e9 | no | 0 |
| `registry_wrong_rhs` | new | 14 | 3.8e11 / 7e12 | 312 (9) | 6.4e4 / 1.5e11 | yes | 0 |
| `unpack_no_zero_floor` | carried | 6 | 1.6e11 / 7.4e11 | 3 (0) | 5.5e8 / 6e8 | no | 0 |
| `unpack_b_reversed` | new | 0 | - | 9 (0) | 2.1e8 / 3.4e8 | no | 0 |
| `esn_readout_rate_for_synaptic` | carried | 0 | - | 4 (4) | 4.3e3 / 1.8e5 | no | 0 |
| `esn_ridge_eta_1e-3` | carried | 0 | - | 6 (6) | 17 / 1.2e2 | no | 0 |

What each one is: the `what` field of the list.

## What the new structure changes about what the tier sees

1. **The parameter view no longer sees every parameter error.** A
   route's model carries the values the arithmetic reads per neuron in
   `P`, beside the nominal rows that `route_fields` (and so the view)
   reads. `P_tau_rest_1e-6`, `P_tau_drive_1e-6`, `P_target_1e-6` and
   `P_one_neuron_tau_rest_1e-6` put an error in `P` alone: 54 fewer
   pointwise leaves fail than for the same error in the nominal row (115
   against 169), those being the view's, and `dS`, `J`, `JY` and the
   stable layer catch them at the same factors. One neuron in each type
   is enough (1.6e3 at the weakest failing pointwise leaf).
2. **An error in the linearisation moves J and JY together.** Both
   Jacobian routines read `route_product`'s `linearize`, so the equality
   "`jacobian_times` = J * Y" cannot see an error there; before the
   refactor a slip in one routine showed up against the other. What
   catches such an error is the stored `J` and `JY` of the anchor, and
   `test_route_models` (central differences of the model's own `rhs`).
   At 1e-9 the stored values see it by factors of 20 to 1e3
   (`lin_A_last_1e-9`, `lin_B_1e-9`, `lin_F_1e-9`, `lin_Fs_1e-9`), the
   same resolution the `jt_*` mutations have. This rests on the anchor
   staying fixed: a linearisation error made BEFORE a re-capture of the
   pointwise references would be in both.
3. **Facilitation's (x, s) entries have one witness, at 1e-9.**
   `lin_Fs_g_1e-9` and `lin_Fs_b_with_g_1e-9` fail 13 pointwise leaves
   (the two options networks' `J` and `JY`) and nothing in the stable
   layer, by 33 to 600. This is item 1 of "What slipped through, and what
   is thin" above, unchanged by the refactor.
4. **`P.coupling` is trusted by `jacobian_times` and not by the assembled
   Jacobian.** `lin_A_spurious_coupling_1e-6` writes an off-diagonal
   entry of A that the model says is zero: 45 `J` leaves fail (3.4e4 to
   9.7e4) and no `JY` leaf, because `jacobian_times` applies only the
   diagonal for a model that states one. A model that understates its
   coupling therefore breaks the equality of the two routines, which
   `test_jacobian_times` checks; `test_route_models` checks the
   statement itself (A exactly zero wherever `coupling` is false).

## What the tier cannot see on the refactored code

Four things are not in the list because no case of the tier reaches
them, so a mutation there would pass:

| Not reached by the tier | What covers it |
|---|---|
| The full (non-diagonal) application of A in `jacobian_times`. No model of today has a coupled state. | `test_route_models`, check 4: forcing every route down that branch must reproduce the diagonal branch and J * Y. Checked by hand on 2026-11-11: that branch scaled by 1 + 1e-9 fails the test (7.8e-10 against 1e-14) while the tier passes. |
| `active_routes` compiling the routes again when `params.std_zero_floor` is switched after `get_params`. | `test_SRNNCellTypePairs`, `test_jacobian_times`, `test_SRNNCellTypePairs_tau_a_spread`, which switch it and check the Jacobian at both settings. |
| `route_product`'s leave-one-out products for three or more components. No network of the tier has three. | `test_route_models`, the three-component case. |
| A presynaptic type with ONE neuron. | `test_c_over_K` and `test_jacobian_times`. They caught a real defect there during the refactor: a route's state came back transposed when its index was a 1 x ns row. |

And one error is invisible everywhere until the chain exists: entry
(m, k) of A written at (k, m) in the assembled Jacobian. For a diagonal A
it changes nothing. The chain's 2 x 2 is not symmetric, so the chain's
own tests (finite differences on a chain route) have to catch it.

## Scan on the refactored code

`regression_tol_mutations_post_refactor('scan')`. Largest factor (relative
difference / tolerance) per kind of leaf; above 1 fails.

| Error | Size | View | dS | J | JY | Trajectories | top-K | Caught |
|---|---|---|---|---|---|---|---|---|
| depression `tau_rec` | 1e-4 | 1e8 | 1.8e8 | 3.3e7 | 2.7e7 | 3.6e4 | 3.8e5 | yes |
| | 1e-6 | 1e6 | 1.8e6 | 3.3e5 | 2.7e5 | 3.6e2 | 3.8e3 | yes |
| | 1e-8 | 1e4 | 1.8e4 | 3.3e3 | 2.7e3 | 3.6 | 38 | yes |
| | 1e-10 | 1e2 | 1.8e2 | 33 | 27 | 0.036 | 0.38 | yes (161 pointwise leaves, 1 stable) |
| | 1e-12 | 1.0 | 1.8 | 0.33 | 0.27 | 3.6e-4 | 3.8e-3 | yes, by 56 pointwise leaves, at the limit |
| first SFA rung | 1e-4 | 1e8 | 1e8 | 8.4e7 | 5.6e7 | 7.4e4 | 2.3e5 | yes |
| | 1e-6 | 1e6 | 1e6 | 8.4e5 | 5.6e5 | 7.4e2 | 2.3e3 | yes |
| | 1e-8 | 1e4 | 1e4 | 8.4e3 | 5.6e3 | 7.4 | 23 | yes |
| | 1e-10 | 1e2 | 1e2 | 84 | 56 | 0.074 | 0.23 | yes |
| | 1e-12 | 1.0 | 1.0 | 0.84 | 0.56 | 7.4e-4 | 2.3e-3 | yes, by 15 pointwise leaves, at the limit |
| `jacobian_times` (x, s) | 1e-4 | 0 | 5e-4 | 3e-4 | 9.9e7 | 4e-7 | 9.9e5 | yes |
| | 1e-6 | 0 | 5e-4 | 3e-4 | 9.9e5 | 4e-7 | 9.9e3 | yes |
| | 1e-8 | 0 | 5e-4 | 3e-4 | 9.9e3 | 4e-7 | 99 | yes |
| | 1e-10 | 0 | 5e-4 | 3e-4 | 99 | 4e-7 | 0.99 | yes (14 pointwise leaves, 1 stable) |
| | 1e-12 | 0 | 5e-4 | 3e-4 | 0.99 | 4e-7 | 9.8e-3 | **no** |

The detection limit: 1e-12 relative for a parameter (at the limit),
between 1e-12 and 1e-10 for a term of `jacobian_times`. In the last five
rows the `dS`, `J` and trajectory columns are not the mutation: they are
the refactor's own distance from the anchor, the background mentioned at
the top of this section.

## How the restore was checked, and two entries that were wrong

After every mutation the harness compared the file's blob with `HEAD`'s
and required `git status --porcelain -- <file>` to be empty; all 86 runs
reported the blob restored. After the last run
`git status --porcelain -- src/` was empty and
`regression_tol_check('anchor')` and `regression_tol_check('pre_0d')`
passed (the figures at the top of this section are from that run),
which also shows that the restored code, not a cached mutated class, was
what MATLAB ran.

Two entries of the first draft of the list were wrong and were corrected
before the numbers above were taken. `lin_A_last_1e-9` and
`lin_A_first_dropped` indexed `A(:, end, end)` and `A(:, 1, 1)` without
asking whether the route has any state: on a route with none the first
throws and the second would have grown the array. Twelve cases threw on
the first run of `lin_A_last_1e-9`, which is how it was seen. Both now
act only where `ns > 0`.

---

# The chain suite: mutation checks on chained depression (2026-11-12, Plan A step 5)

Chained vesicle-pool depression (`route_pool`, `std.model = 'chain'`) was
merged on 2026-11-11 with unit tests and no regression reference: the
anchor holds no chain route. Step 5 gives the chain references of its
own, before chained facilitation is built into the same file. This
section records what the tolerance tier's new **chain suite** was shown
to catch, and that the anchor still catches what it caught.

- Code under test: commit `6edf1f6`. `src/model/routes/route_pool.m` at
  blob `df97415`, `route_model.m` at `274c8c7`, `route_product.m` at
  `01746cf`, `src/model/SRNNCellTypePairs.m` at `2182216`. (The last
  three differ from the refactor's by one tidy that moves no number:
  `route_product` has a `depression` handle and `std_route_factor` is one
  call for every model.)
- The suite: `regression_tol_run(ref, 'chain')`, reference
  `scripts/tests/fixtures/regression_tol_chain.mat` (9.3 MB, 14 threads),
  run with `regression_tol_check('chain')`. 2441 leaves, 1821 of them
  floating point: 1671 pointwise (the regression network with four chain
  routes at gain 1, with and without the staircase, noise off and on;
  its stable twin; a four-cell-type MIXED network in which one
  presynaptic type has a chain, a two-component product, a one-component
  and an undepressed route), 770 stable (the chain and the mixed network
  end to end, transient gain, finite-difference Jacobian, Benettin and
  QR, an 8-job sweep, memory capacity on four identical chains). The
  tolerances are the tier's, unchanged; no override was added and no
  leaf name was added to the discrete list.
- The list: `regression_tol_mutations_chain`, 52 one-site errors, and
  its leaf groups (`regression_tol_mutations_chain('groups')`).
- To repeat: `regression_tol_mutate(regression_tol_mutations_chain(),
  'reference', 'chain', 'groups', regression_tol_mutations_chain('groups'))`,
  about 35 s per mutation with two other MATLAB sessions running; it was
  run in four parts with `save_dir`.
- Background: the unmutated code matches the chain reference **exactly**
  (`'rtol', 0`), so a factor below is all mutation.

## Result for the chain

**All 52 mutations fail the suite, and every one fails by compared
values.** No case threw in any run.

| | Mutations |
|---|---|
| Run | 52 |
| Caught | 52 |
| Caught by compared values | 52 |
| Caught only by a thrown case | 0 |
| Had a case that threw | 0 |
| Caught by the pointwise layer | 50 |
| Caught by the stable layer | 51 |
| Caught by both | 49 |
| Caught ONLY by the pointwise layer | 1 (`chain_params_std_n2_1e-6`) |
| Caught ONLY by the stable layer | 2 (`chain_init_b1_1e-6`, `chain_init_b2_1e-6`) |
| Caught ONLY by overridden leaves | 0 |
| Made an exact (discrete) leaf change | 25, all gross errors |

## Every chain mutation

Pointwise and stable: the number of failing leaves and the weakest /
strongest factor by which a failing float leaf exceeded its tolerance.
A<m><k> is entry (m, k) of A = d(ds_m)/d(s_k), 1 the releasable pool and
2 the reserve. What each one is: the `what` field of the list.

| Mutation | Pointwise: leaves failing | weakest / strongest | Stable: leaves failing | weakest / strongest | Exact leaf changed | Caught by | Cases that threw |
|---|---|---|---|---|---|---|---|
| `chain_A11_dropped` | 102 | 2.6e10 / 3.7e11 | 284 | 4.1e3 / 1.7e12 | yes | values | 0 |
| `chain_A11_sign` | 102 | 5.3e10 / 7.5e11 | 298 | 1.4e4 / 2.5e12 | yes | values | 0 |
| `chain_A11_1e-9` | 102 | 26 / 3.7e2 | 29 | 1 / 9.5e2 | no | values | 0 |
| `chain_A21_dropped` | 102 | 1.3e9 / 3.8e11 | 170 | 1.1e3 / 1.5e11 | yes | values | 0 |
| `chain_A21_sign` | 102 | 2.7e9 / 7.6e11 | 170 | 1.4e3 / 1.1e11 | yes | values | 0 |
| `chain_A21_1e-9` | 102 | 1.3 / 3.8e2 | 8 | 1 / 2.3e2 | no | values | 0 |
| `chain_A12_dropped` | 93 | 4.8e8 / 1.9e11 | 170 | 8.2e2 / 8e10 | yes | values | 0 |
| `chain_A12_sign` | 93 | 9.5e8 / 3.8e11 | 170 | 1.4e3 / 2.2e10 | yes | values | 0 |
| `chain_A12_1e-9` | 91 | 2.9 / 1.9e2 | 22 | 1.1 / 1.6e2 | no | values | 0 |
| `chain_A22_dropped` | 102 | 2.6e9 / 3.9e11 | 255 | 2.3e2 / 3.8e10 | yes | values | 0 |
| `chain_A22_sign` | 102 | 5.3e9 / 7.8e11 | 298 | 6e2 / 1.1e12 | yes | values | 0 |
| `chain_A22_1e-9` | 102 | 2.6 / 3.9e2 | 27 | 1 / 2.1e2 | no | values | 0 |
| `chain_B_dropped` | 102 | 1.1e10 / 3.3e11 | 170 | 4.1e3 / 1.8e11 | yes | values | 0 |
| `chain_B_sign` | 102 | 2.1e10 / 6.6e11 | 170 | 1.2e3 / 3.6e11 | yes | values | 0 |
| `chain_B_1e-9` | 102 | 11 / 3.3e2 | 7 | 1 / 78 | no | values | 0 |
| `chain_B_on_reserve` | 102 | 1.1e10 / 3.3e11 | 170 | 1.6e3 / 6.4e10 | yes | values | 0 |
| `chain_lin_F_b2` | 93 | 4.3e9 / 8.2e11 | 188 | 5.5e2 / 2.1e11 | yes | values | 0 |
| `chain_lin_F_1e-9` | 102 | 8.7 / 9.3e2 | 5 | 1.1 / 2.4e2 | no | values | 0 |
| `chain_Fs_on_reserve` | 102 | 7.7e10 / 1e12 | 170 | 7e2 / 5.7e11 | yes | values | 0 |
| `chain_Fs_1e-9` | 102 | 77 / 1e3 | 18 | 1 / 3.7e2 | no | values | 0 |
| `chain_Fs_as_product` | 102 | 3.7e10 / 9.5e11 | 169 | 1.9e2 / 2e11 | yes | values | 0 |
| `chain_C1_refill_dropped` | 9 | 9.7e9 / 4.3e11 | 327 | 1.9e3 / 3.4e10 | yes | values | 0 |
| `chain_C1_refill_1e-6` | 9 | 9.7e3 / 4.3e5 | 283 | 1 / 2e4 | no | values | 0 |
| `chain_C1_release_dropped` | 9 | 7.5e10 / 4.8e11 | 348 | 1.4e4 / 3e11 | yes | values | 0 |
| `chain_C1_release_1e-6` | 9 | 7.5e4 / 4.8e5 | 301 | 1 / 4.6e4 | no | values | 0 |
| `chain_C2_recovery_dropped` | 9 | 4.5e9 / 8.1e10 | 320 | 3.5e2 / 1.8e10 | yes | values | 0 |
| `chain_C2_recovery_1e-6` | 9 | 4.5e3 / 8.1e4 | 238 | 1.1 / 1.7e4 | no | values | 0 |
| `chain_C2_loss_dropped` | 9 | 4.8e9 / 8.5e11 | 326 | 1.3e2 / 7.8e10 | yes | values | 0 |
| `chain_C2_loss_1e-6` | 9 | 4.8e3 / 8.5e5 | 270 | 1 / 6.9e4 | no | values | 0 |
| `chain_rhs_F_b1b2` | 9 | 3.1e10 / 1e11 | 324 | 82 / 3.9e10 | yes | values | 0 |
| `chain_P_tau_rec1_1e-6` | 111 | 2.9e3 / 8.5e5 | 285 | 1 / 2.4e5 | no | values | 0 |
| `chain_P_tau_rec2_1e-6` | 111 | 1.2e3 / 8.1e4 | 249 | 1 / 3.7e5 | no | values | 0 |
| `chain_P_tau_rel_1e-6` | 111 | 4.3e4 / 4.8e5 | 311 | 1 / 7.1e5 | no | values | 0 |
| `chain_P_n2_1e-6` | 111 | 1.5e3 / 8.5e5 | 279 | 1 / 4.5e5 | no | values | 0 |
| `chain_P_one_neuron_n2_1e-6` | 111 | 1.4e3 / 3.8e5 | 219 | 1 / 2.2e5 | no | values | 0 |
| `chain_nominal_n2_1e-6` | 156 | 1.5e3 / 1e6 | 279 | 1 / 4.5e5 | no | values | 0 |
| `chain_params_std_n2_1e-6` | 9 | 1e6 / 1e6 | 0 | - | no | values | 0 |
| `chain_coupling_understated` | 9 | 2.9e9 / 2.9e11 | 162 | 1.1e3 / 8.7e10 | yes | values | 0 |
| `chain_A_is_diagonal_forced` | 9 | 2.9e9 / 2.9e11 | 162 | 1.1e3 / 8.7e10 | yes | values | 0 |
| `chain_init_b1_1e-6` | 0 | - | 298 | 1 / 2.2e4 | no | values | 0 |
| `chain_init_b2_1e-6` | 0 | - | 287 | 1 / 8.8e4 | no | values | 0 |
| `chain_output_b2` | 36 | 6.1e11 / 1.2e12 | 26 | 5.7e4 / 1.2e9 | no | values | 0 |
| `chain_output_b1b2` | 36 | 2.5e10 / 5.5e11 | 26 | 4.8e2 / 5.8e8 | no | values | 0 |
| `chain_depression_b2` | 36 | 8.5e11 / 8.5e11 | 16 | 4.6e8 / 7.5e8 | no | values | 0 |
| `chain_depression_b1b2` | 36 | 4.7e10 / 4.7e11 | 16 | 5e2 / 5.7e8 | no | values | 0 |
| `chain_registry_product` | 381 (72 of them missing) | 2.5e10 / 9.5e11 | 328 | 2.5e3 / 1.2e12 | yes | values | 0 |
| `chain_j_ss_transposed` | 93 | 8e9 / 3.5e11 | 8 | 2.9e3 / 4.2e4 | yes | values | 0 |
| `chain_jt_full_transposed` | 9 | 2.9e9 / 2.9e11 | 158 | 4.4e3 / 2.5e10 | yes | values | 0 |
| `chain_jt_full_1e-9` | 9 | 26 / 5.2e2 | 24 | 1 / 7.7e2 | no | values | 0 |
| `chain_jt_full_offdiag_1e-9` | 9 | 2.9 / 1.4e2 | 22 | 1.1 / 1.6e2 | no | values | 0 |
| `chain_idx_b1_b2_swapped` | 183 | 2.6e11 / 1.7e12 | 174 | 1.5e3 / 7.9e9 | yes | values | 0 |
| `chain_unpack_b_reversed` | 36 | 8.5e11 / 8.5e11 | 36 | 1.6e8 / 7.5e8 | no | values | 0 |

(`chain_registry_product`: the 72 missing leaves are the view's
`signature.n2` and `signature.std_model` and the like, which a product
route does not have. They are missing because the value is different,
not because a case threw.)

## What the chain suite sees, and through which leaves

1. **The errors the anchor could not hold are caught.** The full-A
   branch of `jacobian_times` (`chain_jt_full_*`): the nine pointwise
   `JY` leaves and, through top-K, the stable runs and the sweep; at
   1 + 1e-9 by factors of 26 to 520 in `JY`. A transposed 2 x 2 in the
   assembled Jacobian (`chain_j_ss_transposed`): all 93 stored `J`, and
   in the stable layer only the QR estimator (8 leaves, all under its
   override), because top-K and the transient gain do not assemble J.
   An understated coupling (`chain_coupling_understated`, and the same
   error on the class's side): `JY` and everything that propagates with
   `jacobian_times`, not `J`. The two errors are separable from the
   table: J alone moved, or JY alone.
2. **One entry of the 2 x 2 at 1e-9 is resolved**, by factors of 190 to
   390 at the strongest pointwise leaf, and by every stored `J` that
   holds a chain route (91 to 102 leaves).
3. **An error in the per-neuron parameters is invisible to the view and
   caught by everything else**: 111 pointwise leaves (dS, J, JY) against
   156 when the nominal row is wrong too. One neuron's n2 is enough
   (1.4e3 at the weakest).
4. **The right-hand side's terms**: each of the four terms of (C1) and
   (C2) at 1 + 1e-6 fails the nine `dS` leaves by 4.5e3 to 8.5e5 and 238
   to 301 stable leaves.
5. **Three leaves were added because nothing else read the code they
   check.** `std_route_factor` (the factor display code draws),
   `params.std_n2` and `params.std_chain` feed no arithmetic, so no
   leaf of the suite as first drafted depended on them. The suite
   stores every route's depression factor through `std_factors` at the
   pointwise states and along the stable runs, and the two params
   fields in the view. `chain_depression_b2` and `chain_depression_b1b2`
   fail 36 + 16 such leaves and no other; `chain_params_std_n2_1e-6`
   fails the nine view leaves and no other.
6. **Memory capacity sees the readout**: `chain_output_b2` fails
   `stable.mc.chain.MC` and `R2` by 6.4e4 times their override. The
   synaptic readout on a chain network works with `SRNN_ESN_reservoir`
   as it stood on 2026-11-12 (blob `4b20396`), with four identical chain
   routes.

## What the chain suite cannot see

| Not reached | What covers it |
|---|---|
| An error in `route_pool`'s `steady_state`. No case of the tier calls it. | `test_route_models` (the closed form against the right-hand side) and `test_std_chain_model`, check 2. |
| The refusals (`ChainZeroFloor`, `ChainCombine`, ...). | `test_std_chain_model`, check 8. |
| A chain route beside the zero floor: refused by design, so no case. | The same check. |
| A presynaptic type with ONE neuron on a chain route (the mixed network's smallest type has 2 in the reduced runs, 8 otherwise). | `test_route_models`; `test_jacobian_times` for the product. |
| `initial` returning different values per neuron. | Nothing yet; FD-015 (the spread) is where that becomes possible. |

## Round-off floor of the chain suite

Measured on the same day, nothing committed:

- **Thread count.** Captured at 14 threads; at 1 and at 4 threads 24 of
  1821 float leaves move, none in the pointwise layer, none of them an
  exact leaf, and none at its layer's own tolerance: all 24 are under
  existing overrides (the closest: the top-K orthonormality defect, and
  QR at 7.1e-8 against 1e-5).
- **A harmless reordering in `route_pool`** (through the harness, as two
  "mutations" that must NOT be caught): the right-hand side's two rows
  with their divisions regrouped, and the four entries of A written with
  divisions in place of the precomputed reciprocals. Neither is caught.
  The rhs probe moves 225 float leaves, the closest to its tolerance
  being `stable.sweep.chain_all{2}.lambda_gap` at 2.7e-11 against 1e-9
  (**36 times inside**); everything else is 2600 times inside or more
  (`dS` 2.8e-17, Benettin local rates 3.8e-13, the finite-difference
  Jacobian 1.1e-11 against 1e-7). The linearisation probe moves 55, the
  closest at 1.7e-4 of its tolerance.

  The one thin margin is worth knowing: `lambda_gap` is the difference
  of the two leading exponents, which are nearly equal on a stable chain
  network (a slow band), so round-off in the exponents is a larger
  fraction of the gap than of either exponent. It holds 1e-9 by 36,
  where the anchor suite's closest leaf under its own probe held it by
  500. No tolerance was changed; the lead decides whether 36 is enough.

## The anchor still catches what it caught

`regression_tol_mutations_post_refactor` (71 mutations) was run again
against the anchor on this code, where a chain model now exists beside
the product: **71 of 71 caught, all by compared values**; 66 by the
pointwise layer, 68 by the stable layer; one with cases that threw
beside its value failures (`tau_rec_I_to_E_param_1e-6`, as before).
`tau_rec_1e-6` gives the previous run's numbers to the digits shown (169
pointwise leaves, 4.3e3 / 1.8e6; 263 stable, 1.1 / 2.4e4). The unmutated
code passes the anchor with the refactor's own distances: 240 of 1749
float leaves moved, the closest `pointwise.options_chaotic.dS` at
4.8e-16 against 1e-12.

## How the restore was checked

All 125 runs (52 + 71 + the two probes) reported the mutated file's blob
restored. After the last, `git status --porcelain -- src/` was empty,
`regression_tol_check('chain', 'rtol', 0)` passed with no leaf moved,
and `run_all_tests('quick')` passed 59 of 59 with a pool of 4, the
byte-exact `test_regression_rapid` and `test_regression_rapid_chain`
among them, which also shows that the restored code, not a cached
mutated class, was what MATLAB ran.

## A chaotic chain network, and the three estimators on it

`test_std_chain_model`'s check 9 asserts only that the three estimators
agree within 0.2 on a small chain network, because every chain network
of 8 to 24 neurons tried there was stable. It takes the regression
network's size and a higher gain to get a clearly chaotic one, and not
much of either: the regression network (n = 100) with its four routes
switched to the chain (tau_rec [2 5], tau_rel 0.25, n2 2), noise off,
`rk4` at fs = 200, top-K lambda_1 over [2, 6] s:

| Gain | chain | product twin (the preset's own `sfa3_std2`) |
|---|---|---|
| 1.0 | +0.095 | -0.176 |
| 1.5 | +0.74 | -0.169 |
| 2.0 | +2.17 | +0.39 |
| 3.0 | +3.74 | +1.81 |

So at the regression preset's own gain the chain network is already
weakly chaotic where the product network is stable, and it stays the
more chaotic of the two at every gain tried. (This is one network and
one seed; it is not the manuscript's network.)

The estimators, on the same seed, noise off:

| Network | Window | fs | top-K lambda_1 | Benettin | QR, largest |
|---|---|---|---|---|---|
| chain, n = 100, gain 2 | [5, 20] s | 200 | +1.2316 | +1.2256 | not run |
| product twin | [5, 20] s | 200 | +0.1865 | +0.1814 | not run |
| chain, n = 100, gain 2 | [10, 60] s | 200 | +1.2085 | +1.2018 | not run |
| product twin | [10, 60] s | 200 | +0.4944 | +0.4658 | not run |
| chain, n = 100, gain 2 | [10, 60] s | 800 | +1.2332 | +1.2369 | not run |
| product twin | [10, 60] s | 800 | +0.5332 | +0.5331 | not run |
| chain, n = 40, gain 4 | [4, 12] s | 200 | +1.166 | +0.757 | +1.000 (822 s) |
| product twin | [4, 12] s | 200 | +0.518 | +0.531 | +0.539 (471 s) |
| chain, n = 40, gain 4 | [10, 60] s | 200 | +1.170 | +1.058 | not run |
| product twin | [10, 60] s | 200 | +0.645 | +0.626 | not run |
| chain, n = 40, gain 4 | [10, 60] s | 800 | +0.9300 | +0.9312 | not run |
| product twin | [10, 60] s | 800 | +0.4429 | +0.4426 | not run |

Three things follow. On a chaotic chain network top-K and Benettin agree
to 0.5 % at gain 2 with the usual step. Where they differ by more (10 %
on the chain and 3 % on the product at gain 4, fs = 200), the cause is
the integration step and not the chain: at fs = 800 they agree to 0.1 %
on both networks. (A different fs is a different trajectory of a chaotic
system, so the rows at 200 and at 800 are not estimates of one number.)
And the full-spectrum QR is not affordable on a chaotic network of this
kind: 14 minutes for 8 s at n = 40 (N = 320), against about a second for
the other two; over that short window, with a 3 s warm-up, the three
differ by finite-time effects (0.76 to 1.17), which the 50 s rows remove
for the two estimators that can afford them.

The chain's exact reference (`regression_rapid_chain_digest.mat`) holds
the n = 100, gain 2 network's top-K and Benettin results over [4, 12] s
(+1.2748 and +1.2742).


# The stf_chain suite: mutation checks on release-coupled facilitation (2026-11-12, Plan A step 8)

Release-coupled facilitation (`route_pool`, `stf.model = 'chain'`: a
release probability `g` that multiplies the output and the rate at which
the releasable pool is used, optionally driven by a calcium trace `h`)
was built on 2026-11-12 (FD-018). Neither the anchor nor the chain
reference holds a coupled route, so it gets a suite of its own, captured
before the per-neuron spread is built on top of it. This section records
what that suite was shown to catch, and that the chain suite still
catches what it caught on the rewritten `route_pool`.

- Code under test: commit `2328d19`. `src/model/routes/route_pool.m` at
  blob `5eecc19`, `route_model.m` at `f41e980`, `route_product.m` at
  `01746cf` (unchanged), `src/model/SRNNCellTypePairs.m` at `7167da0`.
- The suite: `regression_tol_run(ref, 'stf_chain')`, reference
  `scripts/tests/fixtures/regression_tol_stf_chain.mat` (13.0 MB, 14
  threads), run with `regression_tol_check('stf_chain')`. 2679 leaves,
  2017 of them floating point: 1869 pointwise, 810 stable. Two networks:
  `facil`, the regression network with all four routes a chain with `g`
  and `h` (tau_dec 1.5, tau_fac 0.25, G 5, tau_h 0.5), with and without
  the staircase, noise off and on, and its stable twin; and `fmixed`,
  four cell types with every variant of the pool model that has
  facilitation (no depression, a single pool or the chain; with and
  without the calcium trace) beside product routes, a chain without
  facilitation, a coupled route at G = 1 and a stiff one (G 8, tau_fac
  0.1, tau_h 0.05). Stable layer: both end to end, transient gain,
  finite-difference Jacobian, Benettin and QR on n = 6, an 8-job sweep,
  memory capacity with the synaptic readout on four identical coupled
  routes. The tolerances are the tier's, unchanged; no override was
  added and no leaf name was added to the discrete list.
- The list: `regression_tol_mutations_stf_chain`, 91 one-site errors,
  and its leaf groups (`regression_tol_mutations_stf_chain('groups')`).
- To repeat: `regression_tol_mutate(regression_tol_mutations_stf_chain(),
  'reference', 'stf_chain', 'groups',
  regression_tol_mutations_stf_chain('groups'))`, a median of 24 s per
  mutation (59 s at most) with two other MATLAB sessions running; it was
  run in four parts with `save_dir`.
- Background: the unmutated code matches the reference **exactly**
  (`'rtol', 0`; a second capture is `isequaln` to it), so a factor below
  is all mutation.

## Result for coupled facilitation

**All 91 mutations fail the suite, and every one fails by compared
values.** One of them also had cases that threw.

| | Mutations |
|---|---|
| Run | 91 |
| Caught | 91 |
| Caught by compared values | 91 |
| Caught only by a thrown case | 0 |
| Had a case that threw | 1 (`stf_registry_product`, 11 cases, beside 394 stable leaves failing by value) |
| Caught by the pointwise layer | 88 |
| Caught by the stable layer | 85 |
| Caught by both | 82 |
| Caught ONLY by the pointwise layer | 6 (`stf_B_g_1e-9`, `stf_single_A_b1_b1_1e-9`, `stf_nodep_Fs_1e-9`, `stf_nodep_lin_F_1e-9`, `stf_params_tau_h_1e-6`, `stf_params_stf_chain_false`) |
| Caught ONLY by the stable layer | 3 (`stf_init_h_1e-6`, `stf_init_h_one`, `stf_init_g_1e-6`) |
| Caught ONLY by overridden leaves | 0 |
| Made an exact (discrete, logical or text) leaf change | 48 |

## Every facilitation mutation

Pointwise and stable: the number of failing leaves and the weakest /
strongest factor by which a failing float leaf exceeded its tolerance.
The variables are named: `A_b1_g` is d(db1/dt)/dg, `B_h` is d(dh/dt)/dr.
"trace" is a route with the calcium trace, "rate" one without (the rate
drives g), "single" a single pool, "nodep" a route without depression.
What each one is: the `what` field of the list.

| Mutation | Pointwise: leaves failing | weakest / strongest | Stable: leaves failing | weakest / strongest | Exact leaf changed | Caught by | Cases that threw |
|---|---|---|---|---|---|---|---|
| `stf_A_b1_g_dropped` | 108 | 2e9 / 1.8e11 | 180 | 1.6e3 / 2.8e9 | yes | values | 0 |
| `stf_A_b1_g_sign` | 108 | 4e9 / 3.6e11 | 180 | 3.6e3 / 4.8e9 | yes | values | 0 |
| `stf_A_b1_g_1e-9` | 108 | 2 / 1.8e2 | 6 | 1 / 15 | no | values | 0 |
| `stf_A_g_g_trace_dropped` | 108 | 9.6e9 / 5.2e11 | 302 | 3.1e3 / 1.6e10 | yes | values | 0 |
| `stf_A_g_g_trace_sign` | 108 | 1.9e10 / 1e12 | 310 | 9.6e3 / 2.3e10 | yes | values | 0 |
| `stf_A_g_g_trace_1e-9` | 108 | 9.6 / 5.2e2 | 27 | 1 / 13 | no | values | 0 |
| `stf_A_g_g_rate_dropped` | 36 | 3e10 / 1.7e11 | 90 | 1.1e4 / 1.5e10 | yes | values | 0 |
| `stf_A_g_g_rate_sign` | 36 | 6.1e10 / 3.4e11 | 150 | 1.3e4 / 1.6e10 | yes | values | 0 |
| `stf_A_g_g_rate_1e-9` | 36 | 30 / 1.7e2 | 4 | 1.5 / 1.7 | no | values | 0 |
| `stf_A_g_h_dropped` | 108 | 1.6e10 / 1e12 | 180 | 9.2e3 / 2e9 | yes | values | 0 |
| `stf_A_g_h_sign` | 108 | 3.2e10 / 2e12 | 180 | 3.4e3 / 4.1e9 | yes | values | 0 |
| `stf_A_g_h_1e-9` | 108 | 16 / 1e3 | 18 | 1 / 11 | no | values | 0 |
| `stf_A_h_h_dropped` | 108 | 4.1e9 / 1e12 | 276 | 5.9e3 / 2e10 | yes | values | 0 |
| `stf_A_h_h_sign` | 108 | 8.2e9 / 2e12 | 318 | 1.7e4 / 5.7e10 | yes | values | 0 |
| `stf_A_h_h_1e-9` | 108 | 4.1 / 1e3 | 22 | 1.1 / 12 | no | values | 0 |
| `stf_B_h_dropped` | 108 | 5.7e9 / 9.9e11 | 180 | 6.2e3 / 5.5e9 | yes | values | 0 |
| `stf_B_h_sign` | 108 | 1.1e10 / 2e12 | 180 | 3.4e3 / 9.7e9 | yes | values | 0 |
| `stf_B_h_1e-9` | 108 | 5.7 / 9.9e2 | 13 | 1 / 8.4 | no | values | 0 |
| `stf_B_g_dropped` | 36 | 2.9e10 / 3.3e11 | 88 | 3.8e3 / 4.7e8 | yes | values | 0 |
| `stf_B_g_sign` | 36 | 5.7e10 / 6.5e11 | 88 | 1.9e3 / 3.7e8 | yes | values | 0 |
| `stf_B_g_1e-9` | 36 | 29 / 3.3e2 | 0 | - | no | values | 0 |
| `stf_B_b1_dropped` | 108 | 9.9e9 / 3.3e11 | 178 | 1.5e3 / 2.4e9 | yes | values | 0 |
| `stf_B_b1_sign` | 108 | 2e10 / 6.5e11 | 184 | 2.4e3 / 1.1e10 | yes | values | 0 |
| `stf_B_b1_1e-9` | 108 | 9.9 / 3.3e2 | 3 | 1.6 / 2.4 | no | values | 0 |
| `stf_A_b1_g_transposed` | 108 | 2e9 / 1.8e11 | 180 | 2.5e3 / 2.5e9 | yes | values | 0 |
| `stf_A_g_h_transposed` | 108 | 1.6e10 / 1e12 | 180 | 1.2e4 / 6.3e9 | yes | values | 0 |
| `stf_B_h_on_g` | 108 | 5.7e9 / 9.9e11 | 180 | 4.7e3 / 4.9e9 | yes | values | 0 |
| `stf_A_g_g_d_r_for_h` | 108 | 7.6e9 / 4.4e11 | 178 | 8.7e2 / 1.5e9 | yes | values | 0 |
| `stf_lin_d_r_for_h` | 108 | 1.5e10 / 1e12 | 184 | 1.8e3 / 9.2e9 | yes | values | 0 |
| `stf_lin_d_h_for_r` | 36 | 2.9e10 / 3.3e11 | 88 | 3.8e3 / 4.7e8 | yes | values | 0 |
| `stf_A_b1_b1_g_dropped` | 99 | 1.1e7 / 7.1e11 | 186 | 2.1e3 / 2.2e10 | yes | values | 0 |
| `stf_A_b1_b1_g_1e-9` | 108 | 21 / 8.7e2 | 23 | 1.1 / 12 | no | values | 0 |
| `stf_B_b1_g_dropped` | 99 | 1.5e7 / 1.5e11 | 178 | 8e2 / 1.8e9 | yes | values | 0 |
| `stf_single_A_b1_b1_dropped` | 36 | 3.7e10 / 6e11 | 105 | 2.6e3 / 1.4e11 | yes | values | 0 |
| `stf_single_A_b1_b1_1e-9` | 36 | 37 / 6e2 | 0 | - | no | values | 0 |
| `stf_single_A_b1_b1_as_chain` | 36 | 2.2e10 / 8.3e10 | 88 | 2e3 / 1.1e9 | yes | values | 0 |
| `stf_lin_F_g_dropped` | 99 | 7.8e6 / 1.9e11 | 182 | 2.6e3 / 1.5e9 | yes | values | 0 |
| `stf_lin_F_1e-9` | 108 | 23 / 8.1e2 | 1 | 1.3 / 1.3 | no | values | 0 |
| `stf_Fs_b1_dropped` | 108 | 1.3e10 / 1e12 | 179 | 3.1e3 / 3.8e9 | yes | values | 0 |
| `stf_Fs_b1_g_dropped` | 99 | 5.5e6 / 8e11 | 179 | 4.9e2 / 2.4e9 | yes | values | 0 |
| `stf_Fs_b1_1e-9` | 108 | 13 / 1e3 | 21 | 1.1 / 6.1 | no | values | 0 |
| `stf_Fs_g_dropped` | 108 | 2.9e9 / 7.8e11 | 179 | 1.6e3 / 3.1e9 | yes | values | 0 |
| `stf_Fs_g_1e-9` | 108 | 2.9 / 7.8e2 | 5 | 1 / 4.1 | no | values | 0 |
| `stf_Fs_h_as_factor` | 108 | 8.7e9 / 7.9e11 | 179 | 3e3 / 1.8e9 | yes | values | 0 |
| `stf_nodep_Fs_dropped` | 36 | 1e10 / 1.8e11 | 45 | 97 / 1.1e9 | yes | values | 0 |
| `stf_nodep_Fs_1e-9` | 36 | 10 / 1.8e2 | 0 | - | no | values | 0 |
| `stf_nodep_lin_F_1e-9` | 36 | 17 / 1.1e3 | 0 | - | no | values | 0 |
| `stf_rhs_trace_g_decay_dropped` | 9 | 6.1e10 / 1.9e11 | 418 | 1.3e2 / 1e10 | yes | values | 0 |
| `stf_rhs_trace_g_decay_1e-6` | 9 | 6.1e4 / 1.9e5 | 370 | 1.1 / 3.5e3 | no | values | 0 |
| `stf_rhs_trace_g_drive_dropped` | 9 | 1.3e11 / 3.9e11 | 400 | 7.7e3 / 1.2e10 | yes | values | 0 |
| `stf_rhs_trace_g_drive_1e-6` | 9 | 1.3e5 / 3.9e5 | 378 | 1.2 / 8e3 | no | values | 0 |
| `stf_rhs_trace_h_rate_dropped` | 9 | 4.6e10 / 5.8e11 | 404 | 7.4e3 / 1.1e10 | yes | values | 0 |
| `stf_rhs_trace_h_rate_1e-6` | 9 | 4.6e4 / 5.8e5 | 376 | 1.3 / 1.6e4 | no | values | 0 |
| `stf_rhs_trace_h_decay_dropped` | 9 | 4.6e10 / 1.2e12 | 418 | 7.9e3 / 6.9e9 | yes | values | 0 |
| `stf_rhs_trace_h_decay_1e-6` | 9 | 4.6e4 / 1.2e6 | 377 | 1.1 / 1.3e4 | no | values | 0 |
| `stf_rhs_d_r_for_h` | 9 | 3.7e11 / 1e12 | 398 | 7.3e2 / 2.1e9 | yes | values | 0 |
| `stf_rhs_rate_g_decay_dropped` | 3 | 1.3e11 / 1.5e11 | 238 | 6.2 / 6.6e8 | yes | values | 0 |
| `stf_rhs_rate_g_decay_1e-6` | 3 | 1.3e5 / 1.5e5 | 222 | 1 / 3.6e2 | no | values | 0 |
| `stf_rhs_rate_g_drive_dropped` | 3 | 4.7e11 / 5.3e11 | 239 | 5.7e3 / 1.3e9 | yes | values | 0 |
| `stf_rhs_rate_g_drive_1e-6` | 3 | 4.7e5 / 5.3e5 | 219 | 1.1 / 4.8e2 | no | values | 0 |
| `stf_rhs_release_g_dropped` | 9 | 3.2e10 / 5.8e10 | 418 | 15 / 1.3e10 | yes | values | 0 |
| `stf_rhs_release_g_1e-6` | 9 | 8.1e4 / 1.9e5 | 374 | 1.8 / 1.2e4 | no | values | 0 |
| `stf_rhs_F_g_dropped` | 9 | 3.6e10 / 4.6e11 | 400 | 3.1e2 / 3.7e9 | yes | values | 0 |
| `stf_rhs_F_1e-6` | 9 | 9.6e4 / 1.2e6 | 376 | 1.3 / 6.5e3 | no | values | 0 |
| `stf_rhs_F_g_h` | 9 | 1.2e11 / 1.2e12 | 399 | 1.1e3 / 4.7e9 | yes | values | 0 |
| `stf_rhs_nodep_F_1e-6` | 3 | 4.2e5 / 2.5e6 | 172 | 1 / 9.2e2 | no | values | 0 |
| `stf_rhs_single_refill_dropped` | 3 | 8.2e10 / 9.5e10 | 238 | 1e2 / 8.5e8 | yes | values | 0 |
| `stf_rhs_single_refill_1e-6` | 3 | 8.2e4 / 9.5e4 | 213 | 1.2 / 5.5e2 | no | values | 0 |
| `stf_P_tau_dec_1e-6` | 117 | 1.4e3 / 1.9e5 | 383 | 1.2 / 4.5e3 | no | values | 0 |
| `stf_P_tau_fac_1e-6` | 117 | 3.3e4 / 1e6 | 386 | 1.1 / 6.3e3 | no | values | 0 |
| `stf_P_G_1e-6` | 117 | 4.1e4 / 4.6e6 | 390 | 1.4 / 1.1e4 | no | values | 0 |
| `stf_P_tau_h_1e-6` | 117 | 9.6e3 / 1e6 | 377 | 1 / 1.8e4 | no | values | 0 |
| `stf_P_one_neuron_G_1e-6` | 117 | 4.1e4 / 4.5e6 | 372 | 1.1 / 1.9e3 | no | values | 0 |
| `stf_nominal_G_1e-6` | 219 | 4.1e4 / 4.6e6 | 390 | 1.4 / 1.1e4 | no | values | 0 |
| `stf_nominal_tau_h_1e-6` | 162 | 9.6e3 / 1e6 | 377 | 1 / 1.8e4 | no | values | 0 |
| `stf_params_tau_h_1e-6` | 9 | 1e6 / 1e6 | 0 | - | no | values | 0 |
| `stf_params_stf_chain_false` | 9 | - | 0 | - | yes | values | 0 |
| `stf_coupling_b1_g_understated` | 3 | 2.3e10 / 3.8e10 | 84 | 1.5e4 / 1.2e8 | no | values | 0 |
| `stf_coupling_g_h_understated` | 3 | 2.1e11 / 2.3e11 | 42 | 1.5e6 / 3.7e9 | yes | values | 0 |
| `stf_init_h_1e-6` | 0 | - | 388 | 1.2 / 2.8e4 | no | values | 0 |
| `stf_init_h_one` | 0 | - | 404 | 3.8e2 / 4e9 | yes | values | 0 |
| `stf_init_g_1e-6` | 0 | - | 378 | 1.2 / 5.7e3 | no | values | 0 |
| `stf_output_g_dropped` | 42 | 3.8e10 / 8.1e11 | 30 | 1.1e5 / 7.6e8 | no | values | 0 |
| `stf_output_g_h` | 33 | 5.5e11 / 1e12 | 24 | 2.8e5 / 1e9 | no | values | 0 |
| `stf_output_nodep_g_h` | 3 | 3.9e11 / 5.2e11 | 2 | 5.5e8 / 5.8e8 | no | values | 0 |
| `stf_facilitation_g_h` | 36 | 5.3e11 / 7.6e11 | 16 | 7.1e8 / 9e8 | no | values | 0 |
| `stf_facilitation_h` | 36 | 6.8e11 / 8.8e11 | 16 | 8.2e8 / 9.5e8 | no | values | 0 |
| `stf_registry_product` | 0 (and 1869 missing) | - | 394 (and 92 missing) | 1.3e5 / 2.5e10 | yes | values | 11 |
| `stf_route_factor_prod` | 36 | 5.3e11 / 7.6e11 | 16 | 7.1e8 / 9e8 | no | values | 0 |
| `stf_idx_g_h_swapped` | 189 | 6e11 / 3.4e13 | 188 | 2.1e3 / 9.9e9 | yes | values | 0 |
| `stf_unpack_g_reversed` | 36 | 6.8e11 / 8.8e11 | 34 | 2.4e8 / 9.5e8 | no | values | 0 |

(`stf_registry_product`: the registry hands every pool route the product
model. A coupled route then has output-only facilitation and no calcium
trace, so the state vector is shorter (1000 where the reference has 1200
in `facil`): the pointwise cases cannot be evaluated at the stored states
and throw, which is the 1869 missing leaves, while the stable networks
run end to end and fail by value, 394 leaves, plus 92 missing because an
array changed size or a signature field is gone. `stf_params_stf_chain_false`
fails nine logical leaves, `view.stf_chain`, which have no factor.)

## What the stf_chain suite sees, and through which leaves

1. **Every entry of the linearisation that facilitation adds or changes
   is resolved at 1 + 1e-9**: (b1, g), (g, g) with and without the trace,
   (g, h), (h, h), the rate columns of h, of g and of b1, the factor g in
   (b1, b1), F and Fs. By factors of 170 to 1100 at the strongest
   pointwise leaf and 2 to 37 at the weakest, in every stored `J` and
   `JY` of a network that has the entry (108 leaves, or 36 where only
   `fmixed` has the variant). The stable layer sees 11 of the 15 such
   errors too, through top-K and the transient gain, at factors of 1 to
   15: the pointwise layer is what carries a 1e-9 error here.
2. **An entry in the wrong place is caught as a gross error**: (b1, g)
   written at (g, b1), (g, h) at (h, g), the rate column of h in g's
   row, all at 1e9 and above.
3. **d = r where d = h belongs, and the reverse**: in the right-hand
   side (`stf_rhs_d_r_for_h`: 9 `dS` leaves at 3.7e11 and above, 398
   stable leaves), in (g, g) and in the rate column of g in both
   directions (`stf_A_g_g_d_r_for_h`, `stf_lin_d_r_for_h`,
   `stf_lin_d_h_for_r`).
4. **The factor g in the release term and in F are separable.** Dropped
   from the release term alone (`stf_rhs_release_g_dropped`: the
   output-only form on a pool route) the nine `dS` leaves fail at 3e10;
   dropped from F alone (`stf_rhs_F_g_dropped`) at 4e10 and above; at
   1 + 1e-6 either fails `dS` by 8e4 to 1e6 and about 375 stable leaves.
   In the linearisation the pointwise count is 99 and not 108 where an
   error only removes a factor g: each of the nine cases holds its run's
   first sample, where g is exactly 1 and that state's `J` cannot differ.
5. **Each term of the g and h equations at 1 + 1e-6** fails `dS` by 5e4
   to 1e6 and 219 to 378 stable leaves; the trace's two terms
   separately.
6. **An error in the per-neuron parameters is invisible to the view and
   caught by everything else**: 117 pointwise leaves against 162 or 219
   when the nominal row is wrong too. One neuron's G is enough (4e4 at
   the weakest pointwise leaf).
7. **The leaves that exist because nothing else reads the code they
   check** did their work: `stf_facilitation_g_h`, `stf_facilitation_h`
   and `stf_route_factor_prod` (the display factor returned as g h) fail
   36 pointwise and 16 stable `stf_factor` leaves and nothing else;
   `stf_params_tau_h_1e-6` and `stf_params_stf_chain_false` fail the
   nine view leaves and nothing else.
8. **The initial state is seen by the stable layer only**, as for the
   chain: h starting at 1e-6, or at 1 like every other route variable,
   and g at 1 + 1e-6 (378 to 404 leaves).
9. **The understated coupling is caught through the routes whose block
   it makes diagonal**: without (b1, g), the single pool with g
   (`jacobian_times` then drops the entry: 3 `JY` leaves at 2e10, 84
   stable); without (g, h), the route with g and h alone (3 at 2e11, 42
   stable). On the chain with g and h the same understatement changes
   nothing, because that block is not diagonal either way and the class
   reads only whether it is.

## What the stf_chain suite cannot see, and what is thin

| Not reached | What covers it |
|---|---|
| An error in `route_pool`'s `steady_state`. No case of the tier calls it. | `test_route_models` (the closed form against the right-hand side, a long integration, and Plan A 4.6's numbers) and `test_stf_chain_model`, check 2. |
| d = h where d = r belongs in the RIGHT-HAND SIDE of a route without a trace. It cannot be written as an edit that runs: the route has no h to read. | Its form in the linearisation is in the list; `test_route_models` compares the right-hand side with an independent loop. |
| The refusals (`ChainStfOutputOnly`, `ChainStfProductStd`, ...). | `test_route_models`, section 3c; `test_stf_chain_model`, check 6. |
| The zero floor beside a pool route without depression (allowed, not in force). | `test_stf_chain_model`, check 6. |
| `P.coupling` understated where the block stays non-diagonal. Nothing reads it there. | `test_route_models` checks the statement itself against central differences, per variant. |
| `initial` returning different values per neuron. | Nothing yet; FD-015 (the spread) is where that becomes possible. |

Thin, though caught: **the output factor of a route with g and h and no
depression** is read by five leaves only (`stf_output_nodep_g_h`: three
pointwise, two stable, all `synaptic_output` of `fmixed`'s E->SST), at
4e11, because `output` feeds the unpacked synaptic output and nothing
else and only one route of the suite is of that kind. The memory-capacity
case reads `output` on a chain with g and h, not on this variant.

## Round-off floor of the stf_chain suite

Measured on the same day, nothing committed:

- **Thread count.** Captured at 14 threads; at 1 and at 4 threads 24 of
  2017 float leaves move, none in the pointwise layer, none of them an
  exact leaf, and all 24 under existing overrides (the closest: the
  top-K orthonormality defect at 8.6e-5 of its override).
- **Harmless reorderings in `route_pool`** (through the harness, as three
  "mutations" that must NOT be caught): the g and h equations with their
  divisions regrouped, the release term as `(g r) b1 / tau_rel`, and the
  (g, h) entry as a product with the reciprocal. None is caught. The
  release probe moves 360 float leaves, the closest to its tolerance
  being `stable.estimators.facil.benettin.local_lya` at 1.3e-11 against
  1e-9 (**80 times inside**); the g and h probe 130, the closest the same
  kind of leaf at 220 times inside; the linearisation probe 70, the
  closest `JY` at 2.7e-16 against 1e-12. `dS` moves by 2e-16 at most,
  trajectories by 3.5e-16, top-K by 7e-14.

**One criterion of the other two suites does not hold here, and it is
not instability.** The leading direction's LOCAL rate is positive at
times in the stable layer (+1.6 per second in `facil`, +0.37 in
`fmixed`), where the base and chain suites' never reaches zero. Lowering
the gain does not change it (+2.1 at gain 0.1 and at 0.05): it comes
from the coupled route's own block, which is non-normal ((g, h) =
(G - g)/tau_fac is up to 16 per second beside decay rates of about 2,
and (b1, g) = -b1 r/tau_rel), so a tangent vector can grow for a while
although every exponent is negative (top-K lambda_1 -0.36 and -0.58; no
estimator finds a positive exponent; no K retry; sweep lambda_1 -0.69 to
-0.33). The round-off evidence above is what says the margin is enough:
a reordering moves a stable leaf by 1.3e-11 at most.

## The chain suite still catches what it caught

`route_pool` was rewritten from fixed two-pool lines into one code path
over the variables a route has, so the 52 entries of
`regression_tol_mutations_chain` no longer found their text. The same 52
errors were re-expressed on the new lines (names, order and meaning
unchanged) and run against the chain reference: **52 of 52 caught, all
by compared values, no case threw**, and the number of failing leaves in
each layer is the 2026-11-12 table's for every one of the 52
(`chain_A11_dropped` 102 and 284, `chain_registry_product` 309 failing
and 72 missing pointwise and 328 stable, and so on). The unmutated code
passes the chain reference at `'rtol', 0` with no leaf moved, and the
byte-exact `test_regression_rapid_chain` passes: the chain without
facilitation computes what it computed.

## How the restore was checked

All 146 runs (91 + 52 + the three probes) reported the mutated file's
blob restored. After the last, `git status --porcelain` was empty,
`regression_tol_check('stf_chain', 'rtol', 0)` and
`regression_tol_check('chain', 'rtol', 0)` passed with no leaf moved,
and `run_all_tests('quick')` passed 63 of 63 with a pool of 4, the
byte-exact `test_regression_rapid` and `test_regression_rapid_chain`
among them.
