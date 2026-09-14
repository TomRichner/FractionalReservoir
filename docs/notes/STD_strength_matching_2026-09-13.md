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

## 2b. The third direction, built after the fast runs (2026-09-14, 03:15)

Both matchings above make the two-timescale routes WEAKER. Their fast runs
(§ 5) turned the multiple-timescale network chaotic and removed its memory
advantage, so a third preset was added that leaves the two-timescale routes
exactly as published and instead STRENGTHENS the single-timescale route to the
dual's steady state at r_ref:

    1/(1 + r_ref/rho_s) = 1/(1 + r_ref/rho)^2 = 1/9   at r_ref = 0.25, rho = 0.125
    rho_s   = r_ref / 8 = 0.03125
    tau_rel = rho_s * tau_rec = 0.03125 * 2 s = 0.0625 s     (was 0.25 s)

`celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStdSingleMatched_3cond_mu8p25`:
`no_adaptation` and `sfa3_std2` are byte-identical to the unmatched preset;
only `sfa1_std1` changes (tau_rec 2 s kept, tau_rel 0.25 -> 0.0625 s, i.e. the
usage d = 1/tau_rel quadrupled). No route scale, no model code. Consequences:
the single route's steady state is monotone in r and saturates at rho_s, so it
matches the dual product only at r_ref (the dual peaks near r = rho and falls);
its relaxation rate 1/tau_rec + r/tau_rel at r_ref rises from 1.5 to 4.5 s^-1,
so the single-timescale condition depresses both more strongly and faster than
before. The undepressed gain is unchanged in every condition.

Three matchings were therefore tried, in this order: (i) route scale s = 3 on
the dual routes (TR's choice), (ii) equal usage rho_u = 0.342 on both dual
timescales (the control), (iii) single route strengthened, dual as published.
(i) and (ii) equalise by weakening the dual; (iii) by strengthening the single.

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

## 5. Results (fast runs of all three matchings, 2026-09-14)

Near-default = the seven 1-D sweeps at the level nearest the preset default,
pooled (21 networks per condition at fast; 105 in the medium unmatched run).
Memory capacity: total MC over 15 s of delays, 5 paired trials at fast, 15 at
medium. The no-adaptation row is identical physics in every column and the
single-timescale row is identical in the scaled and usage columns; the small
differences between the unmatched and the fast columns there are the run mode
(fewer reps, shorter T), which is the calibration for reading the rest.

| measure | unmatched (`data/topk_med`, medium) | (i) scaled s = 3 (`data/stdscaled_fast`) | (ii) usage rho_u = 0.342 (`data/stdusage_fast`) | (iii) single strengthened (`data/stdsinglematched_fast`) |
|---|---|---|---|---|
| lambda_1 near default, none / single / multiple | +3.47 / +0.35 / **-0.112** | +3.65 / +0.61 / **+2.24** | +3.65 / +0.61 / **+1.01** | +3.65 / (row missing from the figure table, to check: its joint-sample median is -0.51) / **-0.118** |
| networks with lambda_1 < 0, multiple | 105 / 105 | 0 / 21 | 0 / 21 | 15 / 21 |
| occupied median rate, none / single / multiple | 0.483 / 0.269 / 0.230 | 0.429 / 0.215 / 0.300 | 0.429 / 0.215 / 0.276 | 0.429 / 0.211 / 0.236 |
| K the top-K needed for the multiple-timescale jobs | 15 | 60 (D_KY resolved 0-50%) | 15-60 (39-100%) | 15 (67-100%) |
| tau sweep, multiple-timescale, slowest tau 1 -> 30 s | -0.112 -> -0.043, all negative | +2.3 flat, all positive | +0.65 to +1.14, all positive | (figure pending; K = 15 at every level) |
| total MC, none / single / multiple | 0.103 / 0.259 / **0.590** | 0.111 / 0.243 / **0.103** | 0.111 / 0.243 / **0.146** | 0.111 / **0.385** / **0.583** |
| MC horizon (s) | 0.00 / 0.12 / 0.52 | 0 / 0.18 / 0 | 0 / 0.18 / 0 | 0 / 0.30 / 0.48 |
| single vs multiple MC | p = 1.2e-4, d_z = -2.04 | p = 0.0625, d_z = +1.66 (wrong way) | p = 0.0625, d_z = +1.24 (wrong way) | p = 0.125 (floor 0.0625), d_z = -1.05 |

Reading. (i) and (ii): at equal steady-state depression obtained by weakening
the dual routes, the multiple-timescale network is chaotic and has no fading
memory; the unmatched advantage leaned on the squared depression. (iii): with
the dual network as published and the single route brought to the same steady
state, the multiple-timescale network keeps lambda_1 ~ -0.12 and MC ~ 0.58,
while the single-timescale network improves (MC 0.24 -> 0.39, horizon 0.12 ->
0.30). Both mechanisms contribute: depression strength moves the single
condition part of the way, and two timescales at the same strength go the rest.
The claim that survives is "at matched steady-state depression, two depression
timescales extend fading memory further than one", with a smaller effect size
than the unmatched comparison showed.

The occupied median rate of the multiple-timescale condition stays within
0.23-0.30 in every variant, so r_ref = 0.25 does not need revisiting.
`fig_STD_strength_matching` (in each run's `figs/<run>/`) has theta_ss at r_ref
and at the occupied percentiles and the max |log ratio| over the occupied band:
1.11 unmatched, 0.146 scaled, 0.034 usage (the single-strengthened curve can be
drawn by passing that preset as `reference_preset`).

h_KS and D_KY medians per condition are in each run's `parameters.md` Lyapunov
quality table and the `Fig_Sensitivity_*_medians` sheets; they were not
tabulated here because the K column above already carries the message.

### 4d. Replacement paragraph if (iii) becomes the paper preset

> **Short-term depression.** Depression enters {eq:theta} as a product of
> resource variables, each resting at one, following the product-of-components
> description of cortical depression [].
> In the multiple-timescale condition the two depression variables had recovery
> time constants of 2 and 4 s and release time constants of 0.25 and 0.5 s, so
> that at a constant rate  settled to ## 3. The numbers/(1 + 8r) their product to
> 10490891 + 8r)^{-2}1 A single depression variable with the same ratio would settle
> to ## 3. The numbers/(1 + 8r) so a one-timescale condition built that way would be weaker as
> well as structurally simpler, and the comparison between the adapting
> conditions would not isolate the number of timescales. We therefore matched
> the steady-state synaptic output of the single-timescale routes to that of
> the two-timescale routes at a reference rate  = 0.25 the median mean
> firing rate of the multiple-timescale network at the reference parameters:
> the single depression variable kept its recovery time constant of 2 s and its
> release time constant was shortened to 0.0625 s, for which
> ## 3. The numbers/(1 + r_{ref}\,\tau_{rec}/\tau_{rel}) = (1 + 8 r_{ref})^{-2} = 1/91 The two
> conditions agree at  differ away from it ({supFig:std_steady_state});
> the multiple-timescale and no-adaptation conditions are unchanged. Where the
> text attributes a difference between the adapting conditions to the number of
> adaptation timescales, it refers to this strength-matched comparison.

Conditions-table row for the single-timescale condition under (iii):
1 SFA (tau_a 0.25 s) / 1 STD (tau_rel 0.0625 s, tau_rec 2 s), all four routes.

## 6. Open questions for TR (updated 2026-09-14 morning)

1. **Which matching is the paper's.** The fast runs settle the direction
   question TR was asked on 2026-09-13: (i) scale-up and (ii) usage both
   destroy the multiple-timescale result; (iii) keeps it. My recommendation is
   (iii), `..._dualStdSingleMatched_3cond_mu8p25`, run at medium (clone
   `stdUsage_med_config`). `paper_config` still points at (i) and should be
   moved with one line.
2. **What the paper claims.** Under (iii) the honest sentence is that two
   depression timescales at matched strength extend fading memory and keep
   lambda_1 near zero where one does not, AND that depression strength itself
   contributes (the single condition improved when strengthened). The
   product-of-components model stays primary; the matching is a control on the
   single-timescale condition, which is the least disruptive option the
   suggest_updates report recommended.
3. **Route scale.** Keep `synapse_config.<pre>.<post>.scale` in the model (it
   is tested and harmless at 1) or remove it now that no paper preset uses it.
4. **Effect size.** At 5 trials the single-vs-multiple MC test cannot go below
   p = 0.0625; the medium run of (iii) with 15 trials is what the manuscript
   needs, and the unmatched 15-trial p of 1.2e-4 should not be quoted for the
   matched comparison.

