# Presets worth coming back to

Running notes on which `srnn_param_preset` entries produced usable dynamics, and
why. Presets themselves live in `srnn_param_preset`; this file records the
judgement about them, which the code cannot.

> ### ⚠ Every LLE number below predates the Lyapunov window fix (2026-08-13)
>
> These notes were written on the `refactorPSA` branch. The `SDE` branch then
> fixed `SRNNModel2`'s Benettin — and both QR paths — which had **ignored
> `lya_T_interval(1)` entirely** and accumulated from `t = 0` rather than over
> the requested window. Exponents measured before that fix therefore averaged in
> part of the settling transient, and **will shift**.
>
> This is not only bookkeeping. The claims below about *where* the
> `level_of_chaos` crossing sits are what commits `c2d7f22` and `45821dc` used
> to narrow both sweeps to `[0.5, 1.5]`. Re-confirm that the crossing still
> falls inside that range rather than assuming it; a range chosen from shifted
> numbers may no longer bracket it.
>
> The same branch also added `lya_warmup` (default 5 s), which biases nothing at
> that value for Benettin but leaves QR roughly 12% short of its plateau — see
> the table on `SRNNModel2.lya_warmup`. Numbers here were taken with neither
> mechanism in place.

---

## `celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5`

**Status: good.** Current best operating point (noted 2026-08-13).

An E/I network on `SRNNCellTypePairs` with short-term depression on all four
routes, spike-frequency adaptation on E only, and nothing distinguishing the two
populations except the sign of their weights.

### Resolved parameters

| | |
|---|---|
| `n` / `indegree` | 500 / 100 |
| `n_cellTypes`, `cell_type_names`, `f` | 2, `{E, I}`, `[0.5 0.5]` |
| `mu_tilde_relative` | `[5.5 -5.5; 5.5 -5.5]` |
| `sigma_tilde_relative` | `[1.5 1.5; 1.5 1.5]` |
| `level_of_chaos` | 1.0 |
| `F_tracks_network` | false, pinned at `F_ref_n = 500`, `F_ref_indegree = 100` |
| `activation` | `piecewise`, `S_a = 0.8`, `S_c = 0` (scalar; `mu_S_c`/`sigma_S_c` empty) |
| `c` | `[0.5/3, 0]` — SFA scaling, E only |
| `input_config.intrinsic_drive` | 0 |
| STD (all four routes) | `tau_rec = 2`, `tau_rel = 0.25` |
| `n_a` | `[3 0]` in the SFA conditions |

Derived: bulk radius **R = 3.833**, realized spectral radius **4.148**,
`N_sys_eqs = 2250`.

### Why it reads well

- **The transition sits inside the swept range.** `level_of_chaos` crosses zero
  LLE around 0.8 without adaptation and 1.4-1.5 with STD, so a `[0.5, 1.5]`
  sweep resolves the crossing rather than the saturated tails. (Both sweep
  scripts were narrowed to `[0.5, 1.5]` on the strength of this — see commits
  `c2d7f22`, `45821dc`.)
- **The adaptation conditions separate cleanly.** Across the `n` sweep the
  no-adaptation network runs off the ceiling by n≈300, `std_only` sits flat near
  −0.65 from n≈400 up, and `sfa_and_std` tracks zero across the whole range.
- **The connectivity blocks are not interchangeable.** `mu_IE` is the one
  monotone chaos-promoting block under both STD conditions, crossing zero around
  6-7; `mu_EE` and `mu_EI` stay flat; `mu_II` is chaotic only at its
  strongest-inhibition end. That structure is what makes the four-block sweep
  worth running at all.
- **Nothing is silent.** With `S_c = 0` the piecewise nonlinearity sits at
  `phi(0) = 0.5`, so removing the tonic drive leaves the network active rather
  than quiet — which is why `intrinsic_drive = 0` is usable here and would not
  be at a positive setpoint.

### Where it came from

Each step below is one preset, derived from the one above it by a recursive call
in `srnn_param_preset`, so only the named field differs at each step:

```
celltype_pairs_S_c_by_type              piecewise S_a 0.8, loc 1.4, c [0.5/3 0]
 └─ ..._n500                            n 300 -> 500
     └─ ..._n500_fixedF                 F pinned, F_ref_n 500
         └─ celltype_pairs_uniform_std_n500      uniform STD, S_c 0, mu [4 -4]
             └─ ..._mu5p5                        mu -> 5.5, loc 1.4 -> 1
                 └─ ..._mu5p5_nodrive            intrinsic_drive 0.1 -> 0
                     └─ ..._mu5p5_nodrive_sig1p5 sigma 1 -> 1.5
```

### Evidence

`data/param_space/sens_only_sig1p5_aug_13_26_22_09/` — sensitivity analysis
only, seven 1-D sweeps (`n`, `f_E`, `level_of_chaos`, and the four mu blocks).
Combined sheets in `replot_sensitivity_aug_13_26_22_14/figures/`:
`sensitivity_LLE_combined_network.png` and `..._combined_mu_blocks.png`.

**Caveat on that evidence.** It is a `fast` run: 4 levels x 3 reps, `T = [0 10]`,
`fs = 200`. The shapes are what to trust; the per-level histograms are built from
three samples and the magnitudes are not settled. Nothing above has been
confirmed at `medium`, and no tau-sensitivity or param-space stage was run for
this preset.

---

## Also on record

- **`celltype_pairs_uniform_std_n500`** — the `mu [4 -4]`, `loc 1.4`,
  `drive 0.1`, `sigma 1` ancestor. Confirmed at **medium** (11 levels x 15 reps),
  `data/param_space/run_all_aug_13_26_17_12/`, and the medium numbers reproduced
  the `fast2` ones closely, which is the main reason to trust `fast`-mode shapes
  in this family at all.
- **`celltype_pairs_all_std_n500`** — the heterogeneous counterpart: E→E recovers
  faster (`tau_rec` 1 vs 3), I→I depresses harder (`tau_rel` 0.15 vs 0.25), and
  the setpoints split E 0.15 / I 0.25. Under it, *no* grid point in either STD
  condition came out chaotic, against 15-26% for the uniform version. Worth
  revisiting, but the two presets differ in three ways at once (STD
  heterogeneity, setpoints, and symmetric vs asymmetric weights), so that
  contrast is not attributable to any one of them yet.

## Known caveats that apply to all of these

- **`sfa_only` and `sfa_and_std` report LLE = −0.100 far too often**, identical
  to three decimals across different parameters, `n` and `T_range`. That is a
  floor in the Lyapunov computation, not a measurement, and it has not been
  diagnosed. Treat those medians as uninformative.
- **`R` is not independent of `n`** for any of these presets. The cancellation
  that makes it so requires `mu_rel^2 = sigma_rel^2 = 1`; these run `mu_rel` at
  4-5.5. So the `n` sweep varies size *and* criticality together. (CLAUDE.md
  states the independence unconditionally, which is only right in the unity
  case.)
- **A figure export fails roughly one run in four** with
  `MATLAB:class:InvalidHandle`, undiagnosed. Since `4e04c70` it warns and
  continues instead of killing the run, so a run may be quietly missing one PNG;
  the `.fig` is normally still written and can be re-exported.
