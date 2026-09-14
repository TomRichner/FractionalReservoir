# Interpretation of the sfaEI_fast run (2026-09-14): where the paper's MATLAB results stand

*Run: `data/sfaEI_fast` (analyses 53 min) and `figs/sfaEI_fast` (30 of 30 entries, 16
in-paper, 11 min; `report.md` and `manifest.md` there; transcript in `log.txt`).
Preset `celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25`, the paper
network with STD **unmatched** (TR's decision of 2026-09-14), run mode `fast`:
4 levels × 3 reps per 1-D sweep, 7 τ levels × 7 reps, 27 joint-sample points, 20-s
simulations, 5 paired memory-capacity trials, 1 seed for the spectrum and transient-gain
stages. Memory capacity ran with the Wiener process OFF on the SRA1 integrator
(`cfg.mc_preset` chained to the paper preset with σ_u = 0, `cfg.mc_ode_solver = 'sra1'`).
Where a medium-resolution number exists it is quoted from `data/topk_med` (same physics,
15 reps, 11 levels, 15 MC trials with noise ON), so the reader can see which fast numbers
are already stable. Fast numbers are for direction and size; manuscript numbers come from
the medium run of this bundle (`sfaEI_fast_config` with `run_mode = 'medium'`).*

## 1. The headline numbers

| measure | none | single-timescale | multiple-timescale | source |
|---|---|---|---|---|
| λ₁ near the default point, median [IQR] (fast, 21 networks) | +3.65 [+0.07, +4.30] | +0.61 [−0.38, +1.40] | **−0.118** [−0.136, +0.049] | `fig_local_vs_finite_lle` table |
| same at medium (105 networks) | +3.47 [+1.81, +4.33] | +0.35 [−0.06, +0.86] | **−0.112** [−0.113, −0.110] | `data/topk_med` |
| networks with λ₁ < 0 at the default | 24% / 4% (fast / medium) | 29% / 31% | 71% / 100% | same |
| share of time locally expanding (leading local rate > 0) | 74% / 83% | 54% / 53% | 18% / 3% | same |
| Spearman ρ(λ₁, mean rate), 111 networks | −0.46 [−0.63, −0.25] | −0.47 [−0.63, −0.28] | **−0.03** [−0.22, +0.14] | `fig_lle_vs_rate` |
| slowest τ_a sweep 1 → 30 s (E AND I moved), λ₁ | – | – | −0.30 → −0.039, every rep negative | `fig_sfa_EOC_allStd` |
| top-K spectrum, 1 seed, noise on: n₊ / h_KS (bit/s) / D_KY | 8 / 22.6 / 18.1 | 6 / 2.4 / 14.0 | 0 / 0 / 0 | `fig_lyapunov_spectrum` |
| memory capacity, NOISE-FREE, SRA1, 5 trials: total MC | 0.106 [0.096, 0.117] | 0.48 [0.39, 0.57] | **16.3** [15.3, 17.3] | `fig_memory_capacity` |
| memory horizon (R² > 0.1) | 0 s | 0.36 s | **10.6 s** [9.4, 12.1] | same |
| memory capacity, noise ON, 15 trials (medium) | 0.103 | 0.259 | 0.590 (horizon 0.52 s) | `data/topk_med` |
| transient gain, active propagator, peak G / t_peak | 1209 (still rising at 1 s) | 30 (0.87 s) | **8.6 at 0.18 s** | `fig_transient_gain` |
| same, J_xx frozen (the conventional rate-network picture) | 47 000 | 11 000 | 39 | same |
| operating point: α(J_xx) / α(J) / λ₁ | +8.7 / +8.7 / +4.0 | +7.2 / +2.6 / +0.44 | +1.8 / +1.1 / −0.11 | same |
| imbalance examples, λ₁ at μ_EE × 0.5 / 1 / 1.5 | +1.98 / +3.94 / −4.23 | −0.19 / +0.63 / −0.53 | **−0.10 / −0.11 / −0.11** | `fig_eig_heatmap_imbalance` |
| mean rate at μ_EE × 0.5 / 1 / 1.5 | 0.06 / 0.43 / 0.87 | 0.05 / 0.25 / 0.71 | 0.08 / 0.22 / 0.27 | same |

Numerics gate (fast, 2 seeds, 10-s windows): convergence checks A and B pass in every
condition (order slopes 2.0, 2.0, 1.8); the Jacobian check J passes everywhere
(relative error 10⁻¹⁰); the finite-time LLE agreement checks L and C **pass for the
multiple-timescale network** (|Δλ₁| 0.006 between SRA1 and ode45; 0.005 Benettin vs QR)
and **fail for the chaotic and intermittent regimes** (0.11 and 0.17 between integrators;
0.31 Benettin vs QR for single-timescale). That is not a defect of the integrator: over a
10-s window two integrators follow different trajectories of a chaotic system, and the
finite-time exponent scatters by that much between trajectories. The gate is judged at
medium (25 seeds, 20-s windows) and the Methods must state the claim at the level the data
support (section 4).

## 2. What each result means

**Stability (Results subsection 1).** The multiple-timescale network sits at λ₁ = −0.11
in 100% of 105 medium networks and its top-K spectrum has no positive exponent (h_KS = 0,
D_KY = 0). It is stable, not "slightly chaotic". The stronger claim in Codex's §3
("transient expansion within long-interval stability") is **not** what this regime does:
its leading local rate is positive 3% of the time at medium (18% at fast, 1 seed per
level). The regime that shows intermittency is the **single-timescale** one: locally
expanding half the time, λ₁ straddling zero (31% of networks stable), which is exactly the
edge-of-chaos picture. So the paper should say: multiple-timescale adaptation makes the
network reliably and quietly stable; single-timescale adaptation leaves it intermittent at
the boundary; no adaptation is chaotic (λ₁ +3.5, D_KY 18, h_KS 23 bit/s). Keep the
"transient expansion" idea for the Discussion or attach it to the single-timescale
regime.

**Robustness across connectivity (subsection 2).** The medians sheet
(`fig_sensitivity_medians`, now with IQR bands and unclipped axes) shows the
multiple-timescale curve flat at ≈ 0 across E:I ratio, network size, μ_EE, μ_EI and μ_IE,
while the no-adaptation curve swings between −10 and +5 and the single-timescale curve
between −4 and +1.5. The one axis that defeats it is μ_II: past +50% all three conditions
go chaotic together (disinhibition). The imbalance heatmap sheet makes the same point in
one picture: λ₁ pinned at −0.11 at half and 1.5 times the reference E→E mean, while the
no-adaptation network goes from chaotic (+2, +3.9) to saturated-silent (−4.2). The
"less sensitive" statement can now be a defined statistic: range of the median λ₁ over
each sweep (multiple ≈ 0.3 s⁻¹ vs 5–14 for none), or the share of networks in a target
band. The medium run of this bundle gives the manuscript's values.

**Rate and stability are not interchangeable (subsection 2).** `fig_lle_vs_rate`:
Spearman ρ between λ₁ and mean rate is −0.46 for no adaptation (an inverted U: quiet
+0.6, mid +2.8, saturated −6.4), −0.47 for single, and **−0.03 for multiple** (λ₁ ≈ −0.11
at every rate from 0.02 to 0.9). Adaptation decouples stability from operating point:
quiet, mid-rate and saturated multiple-timescale networks are all equally stable. That is
the sentence for the paper; do not use Pearson.

**Long adaptation timescales (subsection 3).** With BOTH cell types' ladders swept
(`tau_levels.md` proves it), λ₁ rises monotonically from −0.30 at τ_slow = 1 s to
−0.039 at 30 s, every rep negative, tracking −1/τ_slow (−0.033 at 30 s). The slowest SFA
rung sets the exponent of the stable network. Note the change from the E-only sweep in
`data/topk_med` (−0.112 at 1 s): moving the I ladder too makes the short-τ network more
stable; the manuscript's Methods already say E and I were swept, and now the data match.

**Reservoir computation (subsection 4).** Noise-free and on SRA1, the multiple-timescale
network's fading memory is **16.3** with a 10.6-s horizon, against 0.48 / 0.36 s for
single-timescale and 0.11 / 0 for no adaptation; d_z = −11 for single vs multiple (p at
the 5-trial floor of 0.0625; medium gives 15 trials and an exact test over 2¹⁵ patterns).
The horizon of ~10 s is the 10-s SFA rung: the stable network holds the input in its
slowest adaptation state. With the paper's dendritic noise ON the same network gave MC
0.59 and a 0.52-s horizon, so the noise-free number is what the untuned recurrent
dynamics can carry and the noisy number is what survives the network's own fluctuations.
Both should be reported; the noise-free one is the headline and the noisy one the
caveat. Do not claim "requires little tuning" (nothing was tuned, but nothing measured
how much tuning would have been needed).

**Non-normal transient amplification (new, for subsection 1 or the Discussion).** A unit
perturbation of the dendritic states is amplified 8.6× within 0.18 s in the stable
multiple-timescale network before decaying (back to 1 by ~2 s at the 3-s horizon of the
earlier note). The conventional frozen-J_xx picture predicts 39× and, in the adapted
regimes, a runaway (α(J_xx) = +1.8 at a typical state); adaptation's feedback at the
operating point brings the frozen rate to +1.1 and the trajectory's own rate is −0.11.
The amplifying pattern is a ~70-neuron E/I-difference pattern (balanced amplification),
unaligned with the Lyapunov direction (|cos| ≈ 0.05); the E/I sum mode is damped, and
the isotropic (noise-average) gain is ≈ 1. Reading for the paper: adaptation is part of
the recurrent dynamics, not a correction to them; instantaneous Jacobian spectra
(the occupancy heatmaps) describe the local expansion that the trajectory then removes.
Full account: `docs/notes/Non_normal_amplification.md`,
`Transient_gain_frozen_vs_active.md`, `Transient_gain_summary_2026-09-13.md`.

**Jacobian occupancy (subsection 1 support).** The heatmaps are now a defined
comparison (three conditions, shared seeds, three E:I balances, matched λ₁ and mean rate
in every panel). The median spectral abscissa of J_xx at sampled states is +1.8 s⁻¹ in
the stable regime: the instantaneous spectrum crosses zero while λ₁ is −0.11, which is
the caption's point that occupancy is a local description and λ₁ the stability test.

## 3. The STD strength question, and what we can still say

Codex's §1 was right that the one- vs multiple-timescale comparison changes depression
STRENGTH along with timescale count (the two depression variables at the same
τ_rel/τ_rec ratio square the steady-state depression). Overnight three matchings were
built and run at fast (`docs/notes/STD_strength_matching_2026-09-13.md` §5):

| matching | multiple-timescale λ₁ near default | MC (noisy, 5 trials) |
|---|---|---|
| unmatched (this run, the paper) | −0.118 | 0.59 (medium, 15 trials) |
| dual routes scaled ×3 (equal output at r = 0.25) | +2.24, chaotic | 0.10 |
| dual usage weakened (ρ_u = 0.34) | +1.01, chaotic | 0.15 |
| single route strengthened (τ_rel 0.25 → 0.0625), dual as published | −0.118 (unchanged) | 0.58; single rises 0.24 → 0.39 |

TR's decision: the manuscript keeps the unmatched comparison. What can honestly be said:

1. The comparison is between three **adaptive architectures** as they would be built
   (no adaptation; one SFA + one STD timescale; three SFA + two STD timescales), not
   between timescale counts at fixed strength. The Methods sentence that depression is
   deliberately not normalised stays, with one added clause: the two-timescale condition
   therefore also has stronger steady-state depression, and the comparison is of the
   combined timescale-plus-strength architecture. Do not write that the difference is
   caused by timescale count alone.
2. The matched runs are a supplementary control that supports a bounded version of the
   timescale claim: with the single-timescale depression strengthened to the dual's
   steady state, the single-timescale network gains fading memory (0.24 → 0.39) but the
   multiple-timescale network still holds more (0.58) and stays stable; and weakening the
   two-timescale depression to the single's strength makes the network chaotic. So both
   strength and temporal structure contribute. That sentence can go in the Results or
   Discussion with the matched runs as supplementary figures once they are at medium.
3. The `fig_STD_steady_state` panel already in the supplement shows the squared
   depression; it should be captioned as such rather than as a matched pair.

## 4. Against the manuscript, section by section

**Results 1 (stability boundary).** Ready: `fig_example_timeseries` (composite with local
rate and accumulating λ₁), `fig_eig_heatmap` / `fig_eig_heatmap_imbalance`,
`fig_lyapunov_spectrum`, `fig_local_vs_finite_lle`. Numbers to quote from the medium run.
Text change: drop the transient-expansion framing for the multiple-timescale regime;
attach it to the single-timescale one; add the non-normal transient paragraph.

**Results 2 (connectivity robustness).** Ready: the sensitivity sheets with bands,
`fig_EI_weights_param_space`, `fig_lle_vs_rate` (new), `fig_eig_heatmap_imbalance`
(new). Text change: state "less sensitive" as a range or a share, quote Spearman not
Pearson, and add the μ_II exception.

**Results 3 (long timescales).** Ready: `fig_sfa_EOC_allStd` now unclipped, on both cell
types, with `tau_levels.md`. The image link in the manuscript still points at the old
folder name; the tag is unchanged.

**Results 4 (recurrent computation).** Ready after the medium run: `fig_memory_capacity`
redrawn from saved statistics, with provenance. Text change: the noise-free protocol
must be stated (σ_u = 0, SRA1, otherwise identical), the noisy result kept as a caveat,
the 10-s horizon tied to the 10-s SFA rung, and the "prospective" PyTorch sentence left as
is (that repository was not touched).

**Results 5 (stimulation).** Nothing changed. `fig_stim_engages_adaptation` still runs on
the separate hand-tuned bursting preset with one seed; Codex §8–§9 remain open.

**Methods / Supplemental Methods.** Now backed by saved artifacts: SRA1 vs ode45 step
refinement (A), strong-order convergence on one Brownian path (B), the analytic Jacobian
vs central finite differences (J), Benettin vs QR vs top-K (C). The LLE-agreement claims
(L, C) hold in the stable regime and must be qualified for the chaotic regimes: agreement
to ~0.1–0.3 s⁻¹ on 10-s windows, i.e. within the finite-time scatter between trajectories;
the medium run is the test. The Methods' statement that the Lyapunov exponent is
Benettin's needs updating: the sweeps use the top-K QR method (K = 15, retried to 60) and
Benettin is the cross-check. The τ sweep description (E and I) is now true.

## 5. Against `fig_to_do_in_future.md` (Codex's ranked list)

| § | item | status |
|---|---|---|
| 0 | verbose setting | done (three levels, `minimal` default) |
| 1 | freeze the 1-TS vs MTS comparison, STD normalisation | decided: unmatched stays; three matchings built, run at fast, documented; one sentence of scope in the Methods, matched runs as a supplementary control |
| 2 | numerical-validity gate | built (five checks, pre-registered thresholds, verdict file, Jacobian artifact); passes for the stable regime at fast; judged at medium |
| 3 | main stability result with local + finite-time LLE across trials | figures built; the across-trial evidence says the stronger claim belongs to the single-timescale regime |
| 4 | robustness + rate vs stability | figures built (bands, `fig_lle_vs_rate`); needs the medium rerun for the 15-network numbers |
| 5 | longest-timescale sweep on both E and I, unclipped | done; medium rerun pending |
| 6 | memory capacity audit and redo | figure and provenance done; noise-free protocol adopted; 15-trial medium run pending (30-trial production optional) |
| 7 | Jacobian occupancy as a defined comparison | done (`fig_eig_heatmap_imbalance`); medium at n = 500 pending |
| 8–9 | bursting / stimulation redesign, DC-LLE in the bursting net | not started |
| 10 | the two composite Methods/Intro figures | the time-series half of figure 2 is done (`fig_example_timeseries`); the Sompolinsky intro figure and the assembly are not |
| 11 | pulse trains | not started, by design |
| 12 | final assembly and consistency checks | after the medium run |

## 6. Against the two reports

`reality_manuscript_suggested.md` status changes: SRA1-vs-ode45 and step-halving are
now "reproducibly verified" with an acceptance record; the Jacobian finite-difference
check is automated and saved; Benettin vs QR exists as check C with the qualification
above; the longest-timescale sweep on both cell types is real; the memory-capacity audit
is done and the protocol changed to noise-free. Still prospective: PyTorch three-condition
run, the CSCS analyses, the absolute figure path, per-figure provenance (partly: each
figure root now has `manifest.md`, `command_window.log` and the tables).

`suggest_updates_main_matlab_model.md`: its recommended option 1 (keep the product model,
add a strength-matched control) is what the matched presets are; the result is that the
control changes the conclusion, which is why the paper's claim must be stated at the
architecture level. Its Jacobian-occupancy recommendations are implemented. The bursting
and PyTorch recommendations are untouched.

## 7. Message to the manuscript agent

1. Keep the three conditions and the unmatched comparison. State in the Methods that the
   multiple-timescale condition has both more timescales and stronger steady-state
   depression, and phrase every contrast as between adaptive architectures. Use the
   matched runs (when at medium) as a supplementary control with the two-sentence
   summary in section 3 above.
2. Results 1: the multiple-timescale network is stable (λ₁ = −0.11, no positive
   exponents, locally expanding 3% of the time); intermittent transient expansion is the
   single-timescale regime's property. Add the non-normal transient amplification
   paragraph (8.6× at 0.18 s; the frozen picture overstates it and predicts a runaway).
3. Results 2: quote the range of the median λ₁ per sweep, the μ_II exception, and
   Spearman ρ of −0.46 / −0.47 / −0.03; add `fig_lle_vs_rate` and
   `fig_eig_heatmap_imbalance`.
4. Results 3: unchanged claim, now on both cell types, with the −1/τ reading.
5. Results 4: noise-free SRA1 protocol; headline MC and horizon from the medium run;
   noisy result as the caveat; horizon = the slowest SFA timescale.
6. Methods: top-K QR is the estimator, Benettin the check; the numerics claims are
   backed by `numerics_verdict.md`, with the finite-time qualification for the chaotic
   regimes; the τ sweep moves E and I.
7. Every number in the text comes from the medium run of `sfaEI_fast_config` (with
   `run_mode = 'medium'`), not from this fast run.

## 8. What is pinned down, what changed, what is next

Pinned down: the sign and size of λ₁ in the three regimes; the flatness of the
multiple-timescale network across five of six connectivity axes and across E:I balance;
the −1/τ dependence on the slowest SFA rung; the decoupling of rate and stability under
multiple-timescale adaptation; the order of the memory-capacity result (three orders of
magnitude noise-free); the transient-gain picture; the numerical convergence and the
Jacobian.

Changed since the ranked list was written: the STD-matching question is answered rather
than open (and the answer is to keep the unmatched comparison with honest scope); memory
capacity is measured noise-free on SRA1; the τ sweep moves both cell types; the sweeps
store local rates and the transient-divergence scalars; every entry point logs its
transcript; the model has a per-route weight `scale` (unused by the paper preset).

Next: the medium run of this bundle (about two hours plus memory capacity), then the
numerics verdict at 25 seeds, then the figure assembly. The matched presets at medium are
optional supplementary material. The bursting/stimulation redesign (§8–9) is the one
Results subsection with no new evidence.
