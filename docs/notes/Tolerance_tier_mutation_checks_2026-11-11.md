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
