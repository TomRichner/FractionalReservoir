# STD strength matching between one and two depression timescales (2026-09-13)

*Decision record for the paper's single- vs multiple-timescale comparison.
Equations and the authoritative statement:
`docs/EquationsParametersDocs/Equations_stability_paper.md`, § "2026-09-13: STD
strength matching". Presets: `srnn_param_preset` cases
`celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStdScaled_3cond_mu8p25` and
`..._dualStdUsage_3cond_mu8p25`. Test: `scripts/tests/test_route_scale.m`.
Steady-state curves over the occupied rates: `fig_STD_steady_state`.*

## 1. The decision

**Made by TR, 2026-09-13**, on the plan for the overnight of 2026-09-13/14.

The paper compares no adaptation, one-timescale adaptation (1 SFA + 1 STD) and
multiple-timescale adaptation (3 SFA + 2 STD) and attributes the differences to
the *temporal structure* of adaptation. SFA is normalised for that comparison:
`c_eff = c/K` makes the steady-state adaptation `c·r` whatever K is. STD was
not. With both depression timescales at ρ = τ_rel/τ_rec = 0.125 (τ_rec [2 4] s,
τ_rel [0.25 0.5] s) and the single-timescale routes at (2, 0.25), the
two-timescale steady state 1/(1 + r/ρ)² is the **square** of the one-timescale
1/(1 + r/ρ), so adding a timescale also multiplied the depression. Codex's audit
(`scripts/paper/fig_to_do_in_future.md` § 1) and
`StochasticPlasticDynamicalSystemPaper/reports/suggest_updates_main_matlab_model.md`
both flagged this as the model-defining decision to freeze before any final run.

**What is matched:** the steady-state synaptic output θ_ss(r) = r·Π_m b_m(r) of
the two-timescale routes equals that of the one-timescale routes at a reference
rate r_ref, chosen by rule as the occupied median rate (below).

**How (primary): a multiplicative scale on the two-timescale routes.** The
depression variables, θ and the ODE are untouched; the weights of the four
`sfa3_std2` routes are multiplied by

    s = θ_single(r_ref) / θ_dual(r_ref) = 1 + r_ref/ρ = 3.

Implemented as `synapse_config.<pre>.<post>.scale`, applied to `params.W` in
`SRNNCellTypePairs.get_params` only (the drawn `W` stays unscaled; dynamics,
both Jacobians and `jacobian_times` read `params.W`, so they agree by
construction; the ESN `'synaptic'` readout's route-redundancy check compares the
scale too). `no_adaptation` and `sfa1_std1` are bit-identical to the unmatched
preset's.

**Direction:** scale the dual condition *up* to the single one (s = 3 on dual,
1 on single), so the single-timescale condition stays exactly the current
network. Alternatives offered and not taken: single scaled down to dual
(s ≈ 0.33), or both adapted conditions scaled to no adaptation (s = 3 and 9).

**Control (also run):** the usage-matched preset, which needs no code change.

## 2. The alternative that was raised and overruled

The recommendation put to TR was to match by **usage**: keep τ_rec = [2 4] s and
choose τ_rel so the product equals the single factor at r_ref. With equal usage
ρ_u on both timescales,

    (1 + r_ref/ρ_u)² = 1 + r_ref/ρ   →   ρ_u = r_ref / (√(1 + r_ref/ρ) − 1) = 0.25/(√3 − 1) = 0.34151
    τ_rel = ρ_u · τ_rec = [0.68301, 1.36603] s.

Reasons it was recommended: τ_rel is Varela's release strength d_m, so this is
the tuning the product model itself offers; it changes no code; and it leaves
the **undepressed (low-rate) recurrent gain identical across the three
conditions**. Its cost: each timescale depresses less and the depression
dynamics are slower (rate 1/τ_rec + r/τ_rel per variable).

The concern with the scale, stated to TR before the choice: because b_m → 1 as
r → 0, the scaled dual condition runs at **3× the undepressed recurrent gain**
of the other two conditions. The match holds at the operating point, not at
low rates, and a change in effective W scale is itself a change the paper's
`level_of_chaos` sweep treats as a stability axis. TR reaffirmed the scale as
primary and asked for the usage variant as the control, so the morning
comparison (§ 5) shows the consequence rather than arguing it.

## 3. The numbers

**r_ref = 0.25.** Rule: the median mean firing rate of the multiple-timescale
condition at the default point of the 1-D sweeps of the medium run
`data/topk_med` (top-K, 15 reps per level), rounded to 0.05. Per sweep, the
level nearest the preset default:

| sweep            | sfa3_std2 median rate |
|------------------|-----------------------|
| `f_E`            | 0.236                 |
| `level_of_chaos` | 0.223                 |
| `mu_EE_relative` | 0.239                 |
| `mu_IE_relative` | 0.251                 |
| `n`              | 0.231                 |

(`mu_EI_relative` and `mu_II_relative` were excluded: their default is
negative and the nearest-level lookup landed on the wrong sign.) Pooled over
the 105 near-default reps:

| condition       | median | p5    | p25   | p75   | p95   |
|-----------------|--------|-------|-------|-------|-------|
| no_adaptation   | 0.483  | 0.027 | 0.274 | 0.729 | 0.915 |
| sfa1_std1       | 0.269  | 0.033 | 0.167 | 0.432 | 0.791 |
| sfa3_std2       | 0.237  | 0.050 | 0.220 | 0.249 | 0.346 |

The joint 64-point parameter-space sample, which ranges far from the default,
has medians 0.250 / 0.159 / 0.154 for the same three conditions; it was not
used because it is not the operating point.

**At r_ref = 0.25, ρ = 0.125:** one factor 1/(1 + 2) = 1/3; unmatched product
1/9; scale s = 3; usage ρ_u = 0.34151, τ_rel = [0.68301, 1.36603] s. Both
matched variants give θ_ss(0.25) = 0.25/3 = 0.0833, against 0.0278 unmatched.

| quantity at r_ref            | single (1 TS) | dual unmatched | dual scaled | dual usage |
|------------------------------|---------------|----------------|-------------|------------|
| Π b_m                        | 0.333         | 0.111          | 0.111       | 0.333      |
| θ_ss = s · r · Π b_m         | 0.0833        | 0.0278         | 0.0833      | 0.0833     |
| undepressed gain factor (r→0)| 1             | 1              | 3           | 1          |

## 4. Replacement manuscript text

The manuscript repo is not edited here; this is ready to paste into
`Manuscript5.md`.

### 4a. What it says now (Supplemental Methods, "Short-term depression", ~line 375)

> Unlike adaptation, depression is deliberately not normalized by its number of
> timescales: it enters as a product rather than a sum, each factor rests at
> one, and adding a timescale is meant to deepen depression rather than
> redistribute it, following the product-of-components description of cortical
> depression [@varelaQuantitativeDescriptionShortTerm1997].

and, earlier in the same paragraph, "so the product in {eq:theta} settles to
$(1 + 8r)^{-2}$. At $r = 1$ a single timescale reduces synaptic output ninefold
whereas two reduce it by a factor of 81, the square ({supFig:std_steady_state})."

### 4b. Replacement paragraph (scaled variant primary, usage control)

> **Short-term depression.** The same depression time constants were applied to
> all four connection routes. In the single-timescale condition the recovery
> time constant was eight times the release time constant, so at a constant
> rate $r$ the depression variable settles to $b = 1/(1 + 8r)$. Depression
> enters {eq:theta} as a product of resource variables, each resting at one,
> following the product-of-components description of cortical depression
> [@varelaQuantitativeDescriptionShortTerm1997]; with the same ratio on both
> timescales the two-timescale product would settle to $(1 + 8r)^{-2}$, the
> square of the single-timescale value. Adding a depression timescale would
> then change the strength of depression as well as its temporal structure,
> and the comparison between the adapting conditions would not isolate the
> number of timescales. We therefore matched the steady-state synaptic output
> of the two-timescale routes to that of the single-timescale routes at a
> reference rate $r_{ref} = 0.25$, the median mean firing rate of the
> multiple-timescale network at the reference parameters. In the primary
> configuration the two-timescale routes retained their release and recovery
> time constants and their weights were multiplied by a route scale $s =
> \theta^{(1)}(r_{ref})/\theta^{(2)}(r_{ref}) = 1 + 8 r_{ref} = 3$, so that the
> recurrent input to neuron $i$ is $\sum_j s_{q_j \to q_i} w_{ij} \theta_j$
> with $s = 3$ on all four routes of the multiple-timescale condition and $s =
> 1$ otherwise. This equalises synaptic output at the operating point; because
> the resource variables approach one at low rates, the undepressed recurrent
> gain of the multiple-timescale condition is three times that of the other
> conditions. As a control we also ran a usage-matched configuration in which
> the weights were unchanged and the release time constants of the two
> depression variables were lengthened to $\tau_{rel} = 0.683$ and $1.366$ s,
> giving both variables the usage ratio $\tau_{rel}/\tau_{rec} = 0.342$ for
> which $(1 + r_{ref}\tau_{rec}/\tau_{rel})^{2} = 1 + 8 r_{ref}$. The two
> configurations agree at $r_{ref}$ and differ away from it
> ({supFig:std_steady_state}). The single-timescale and no-adaptation
> conditions are identical in both. Where the text attributes a difference
> between the adapting conditions to the number of adaptation timescales, it
> refers to this strength-matched comparison; the earlier unmatched
> configuration compared timescale count and depression strength together.

### 4c. Replacement rows for {#supTable:adaptation_timescales}

| Condition                              | $\tau_a$ (s)    | $\tau_{rel}$ (s) | $\tau_{rec}$ (s) | route scale $s$ |
| -------------------------------------- | --------------- | ---------------- | ---------------- | --------------- |
| No Adaptation                          | none            | none             | none             | 1               |
| Single-Timescale Adaptation            | 0.25            | 0.25             | 2                | 1               |
| Multiple-Timescale Adaptation (primary)| 0.25, 1.581, 10 | 0.25, 0.5        | 2, 4             | 3               |
| Multiple-Timescale, usage-matched (control) | 0.25, 1.581, 10 | 0.683, 1.366 | 2, 4         | 1               |

Caption addition: "The multiple-timescale routes are strength-matched to the
single-timescale routes at $r_{ref} = 0.25$: in the primary configuration by a
route weight scale $s = 3$, in the control by the release time constants."

The parameters table rows for $\tau_{rel}$ (E to E … I to I) stay `[0.25 0.5]`
for the primary configuration; add a row "route scale $s$ (all four routes,
multiple-timescale condition) | STD strength match | `3` | –".

## 5. Results: TO FILL after the fast runs

Medians over the successful jobs of each run; the unmatched column is the
*medium* run (`data/topk_med`, 15 reps per level, T = 20 s) and the two matched
columns are *fast* runs (`data/stdScaled_fast`, `data/stdUsage_fast`; fewer
reps, shorter T), so compare direction and size, not decimals.

| condition | measure   | unmatched (topk_med, medium) | scaled (fast) | usage (fast) |
|-----------|-----------|------------------------------|---------------|--------------|
| no_adaptation | λ₁    | TO FILL | (identical physics) | (identical physics) |
| no_adaptation | mean rate | TO FILL | | |
| sfa1_std1 | λ₁        | TO FILL | (identical physics) | (identical physics) |
| sfa1_std1 | mean rate | TO FILL | | |
| sfa3_std2 | λ₁        | TO FILL | TO FILL | TO FILL |
| sfa3_std2 | mean rate | TO FILL | TO FILL | TO FILL |
| sfa3_std2 | h_KS      | TO FILL | TO FILL | TO FILL |
| sfa3_std2 | D_KY      | TO FILL | TO FILL | TO FILL |

Also to record: the new occupied median rate of `sfa3_std2` in each matched run
(if it moves far from 0.25 the match point should be revisited once), and the
`fig_STD_steady_state` matching table (θ_ss at r_ref and at the 5th/95th
occupied percentiles, the maximum relative mismatch over the occupied range for
both matchings).

## 6. Open questions for TR

1. **Direction.** Confirm scaling the dual condition up (s = 3) after seeing
   § 5. The alternative that keeps every condition at the same undepressed gain
   is the usage control; if its λ₁ ordering across the three conditions is the
   same as the scaled run's, the choice is presentational and the scale's
   3× low-rate gain is a caveat in the Methods; if the orderings differ, the
   paper's claim depends on the choice and that needs a decision.
2. **Which is primary in the manuscript.** § 4b is written with the scale as
   primary and usage as control, per the decision. Swapping them is a
   paragraph edit, not a rerun, since both are run.
3. **Does the product model stay primary?** Codex's option 2 (per-timescale
   exponents q_m with Σq_m = 1) was not taken; the Methods keep the Varela
   product and add the match. Confirm.
4. **r_ref after the rerun.** r_ref was set from the unmatched network. If the
   matched network's occupied median is not ≈ 0.25, decide whether to iterate
   once (the rule stays the same, the literal changes) or to keep 0.25 and
   report the difference.
