# Handoff: top-K Lyapunov measures, non-normal amplification, intermittency, stimulation

*Written 2026-09-13 13:50 on R5611351, Claude Code session
`90c5825b-19c1-4546-90a4-703d562404ba`, repo `FractionalReservoir`, branch
`main`, HEAD `3a64568` ("Drop the config-level model override; the top-K
cap is 60 in analysis_run_config"). Working tree clean. Written because
the session is near its context limit; the next session should be able
to resume from this file plus the commit messages of `5a9d409..3a64568`.*

## 1. What was built in this session (2026-09-12 → 13)

All committed on `main`. Read the commit messages: they record the numbers
and the wrong turns.

### 1.1 Estimators and model features (earlier in the session)

* **`lya_method = 'topk'`** (`src/model/lyapunov/lyapunov_topk.m`): the K
  largest Lyapunov exponents by the discrete QR method on the stored
  trajectory; Heun on the fiducial grid; economy QR every `lya_dt`;
  returns `LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya,
  n_positive, h_KS(_bits), D_KY, D_KY_resolved, cond_max, orth_defect_max,
  Q_final`. Verified vs the ode45 QR, Benettin and Liouville
  (`scripts/tests/test_lyapunov_topk.m`).
* **Matrix-free Jacobian product** `SRNNCellTypePairs.jacobian_times(S, Y, params)`
  (== `compute_Jacobian_fast * Y` to 1e-16; `test_jacobian_times`). Any
  change to `compute_Jacobian_fast` or `dynamics_fast` MUST be mirrored.
  Top-10 on the 4000-state network: 139 → 8 s.
* **Per-neuron SFA ladders** `tau_a_spread` / `tau_a_seed` / read-only
  `tau_a_matrix` (log-normal at both ladder ends, geometric between;
  preset `..._tauSpread0p05_...`; `test_SRNNCellTypePairs_tau_a_spread`).
  Measured: does NOT speed Benettin or top-K warm-up
  (`docs/notes/Lyapunov_estimation_methods.md` §2.4).
* **Retry on an unresolved D_KY**: `lya_K_auto` / `lya_K_max` on the class;
  doubles K on the stored trajectory; nests exactly (`test_lya_K_auto`).
* **`SRNNCellTypePairs.lya_summary()`** — the one definition of the
  per-run scalars: `LLE, lambda_2, lambda_gap, n_positive, h_KS_bits,
  D_KY, D_KY_resolved, K_used, retries, n_positive_at_K, cond_max,
  orth_defect_max, lambda_1_drift, frac_local_positive, p95_finite_0p2s,
  mean_positive_excursion_s, lead_frac_x/_sfa/_std/_stf`, plus the
  leading direction's local-rate series. Statics `lya_summary_fields`,
  `transient_divergence(local, dt, window_s)`, `leading_vector_blocks`.

### 1.2 Sweeps, registry, figures, stages (this part of the session)

* `src/presets/analysis_run_config.m` `pack()`: every sweep and mode now
  uses `lya_method 'topk'`, `lya_K 15`, `lya_K_auto true`,
  **`lya_K_max 60`**, `lya_dt 0.05`. Benettin comparison: a config sets
  `cfg.model.lya_method = 'benettin'` after calling it (same seeds, same
  networks). NO config-level `model_overrides` layer (added and removed
  again on TR's request — `6eb7c69`/`3a64568`).
* `ParamSpaceAnalysis2.run_single_job` copies every `lya_summary` scalar
  into each result (`nan_lya_fields` for failures) plus
  `local_rate_lead` / `t_lya_lead`.
* **`src/presets/sweep_metrics.m`** — the registry of measures the figures
  plot (key, field, label, stem, ranges, zero line, `in_sheets`). Sheets
  now exist for `lle, r, hks, dky`: `fig_sensitivity_medians`
  (`Fig_Sensitivity_<stem>_medians`), `fig_sensitivity_analysis_allStd`
  (`Fig_Sensitivity_<stem>_core/_mu`), `fig_param_space_allStd`
  (`Fig_ParamSpace_<stem>`), `fig_EI_param_space` /
  `fig_EI_weights_param_space` (`Fig_EI_[Weights_]ParamSpace_<stem>`,
  shared layout in `src/figures/helpers/ei_metric_sheet.m`).
  `fig_sfa_EOC_allStd` gained `Fig_SFA_EOC_blocks` (leading-vector
  SFA/STD/x fractions vs τ_max: 80-85% SFA across 1-30 s).
* `write_run_parameters_md`: a **Lyapunov quality** table per sweep and
  condition (share D_KY resolved, share n₊ at K, K used, worst
  conditioning, median |λ₁ drift|, orthogonality defect).
* **`src/analysis/run_eig_heatmap.m`** now runs top-K and records, per
  sampled state, the numerical abscissa ω = max eig((J_xx+J_xxᵀ)/2) and
  spectral abscissa α = max Re eig(J_xx) of the **dendritic block**
  (`num_abscissa_by_cond`, `spec_abscissa_by_cond`, `J_times_by_cond`,
  `lya_by_cond`). Figure: `src/figures/fig_transient_amplification.m`.
* **`src/analysis/run_lyapunov_spectrum.m`** + `src/figures/fig_lyapunov_spectrum.m`:
  per regime, seeds, noise on AND off, top-K at K = 50/200/300 by mode,
  lya_K_max 2K; `n_override` for tests. Registered in
  `run_all_paper_analyses` (stage `lyapunov_spectrum`) and in
  `paper_config`, `single_multi_TS_independent_config`, `_med_config`,
  `topk_smoke_fast_config`, `topk_med_config`; `test_run_modes` lists it.
* Tests added: `test_lya_K_auto`, `test_sweep_metrics`,
  `test_lyapunov_spectrum_stage` (+ `test_jacobian_times`,
  `test_SRNNCellTypePairs_tau_a_spread` earlier). Baseline suite that
  passes today (2026-09-12): `test_SRNNCellTypePairs`,
  `_S_c_heterogeneity`, `test_lyapunov_topk`, `test_preset_golden`,
  `test_preset_conditions`, `test_run_modes`, `test_sensitivity_refactor`,
  `test_psa_model_class`, `test_psa_loaders`, `test_psa_validate_defaults`,
  `test_psa_saveload`, `test_figure_report`, `test_resolve_data_file`,
  `test_c_over_K`, `test_n_a_dependent`, `test_pairs_single_celltype`,
  `test_esn_route_redundancy`. Many other tests in `scripts/tests/` are
  stale (SRNNModel2-era); do not rely on them.
* Configs/runs: `scripts/paper/topk_smoke_fast_{config,run}.m` (fast smoke,
  0.67 h, 24/24 figures) and **`scripts/paper/topk_med_{config,run}.m`**
  (medium, run overnight 2026-09-12/13 into `data/topk_med` and
  `figs/topk_med`; TR "really likes" the output).

### 1.3 Medium-run headline numbers (`figs/topk_med`, n = 500, 3 seeds, K = 200, 40 s)

| regime | n₊ | h_KS (bit/s) | D_KY | λ₁ |
|---|---|---|---|---|
| no adaptation, noise on | 5 [4, 8] | 12.2 [8.7, 20.9] | 11.6 [9.6, 16.7] | +2.98 |
| no adaptation, noise off | 5 [4, 8] | 10.3 [6.9, 18.3] | 11.2 | +2.79 |
| single-timescale, noise on | 7 [4, 8] | 2.5 | 15.7 | +0.48 |
| single-timescale, noise off | 8 [8, 10] | 3.6 | 19.4 | +0.70 |
| multiple-timescale, both | 0 | 0 | 0 | −0.109 |

Transient amplification (J_xx, median over 300 states): ω = 82 / 74 / 46
s⁻¹, α = 7 / 6 / 0.3 for no / single / multiple-timescale adaptation.
Note noise RAISES h_KS in the no-adaptation regime here (opposite of
Engelken et al.'s fluctuating-input result); with 3 seeds the ranges
overlap — worth a look, not yet explained.

### 1.4 Notes written

`docs/notes/Lyapunov_estimation_methods.md` (§2.4 alignment/spread, §4.2
cost, §4.4 what the sweeps report), `Cov_Lya_vectors.md`,
`topK_reading_list.md`, `CLV_and_spectrum_reading_list.md` (papers now in
`C:\Users\m218089\Desktop\github_repos\PDF2md\StabilityPaper\*_figures\*.md`),
**`Non_normal_amplification.md`** (+ PDFs). Manuscript:
`C:\Users\m218089\Desktop\github_repos\StochasticPlasticDynamicalSystemPaper\Manuscript5.md`
(TR will rewrite it around top-K; image links to
`Fig_Sensitivity_mean_rate_medians` change case, param-space sheets are
per measure). Kreiss constant code (Mitchell 2020/2021, wrapper
`test_Kreiss/Kreiss.m`) lives in `C:\Users\m218089\Desktop\github_repos\Kreiss`;
decided NOT to use it (see §3.3).

## 2. The non-normal amplification discussion (2026-09-13)

### 2.1 What the current figure shows and why J_xx

`fig_transient_amplification`: ω(J_xx) is the largest instantaneous
growth rate of ‖δx‖ (max eigenvalue of the symmetric part), α(J_xx) the
asymptotic rate (max Re eigenvalue); ω − α per state is the "non-normal
margin". ω ≫ α in every regime because a Dale's-law W has a rank-one
mean structure whose largest singular value (√n·μ) far exceeds its
eigenvalue outliers; the direction achieving ω is the E-minus-I
difference mode: Murphy & Miller's balanced amplification, Hennequin's
non-normal amplification. Adaptation lowers ω 82 → 46 (operating point:
STD depresses W_eff, SFA moves neurons to shallower φ′).

Why the dendritic block and not the full J: the numerical abscissa is not
invariant to a diagonal rescaling of the state, and the full J mixes rows
in units 1/τ_d = 10, 1/τ_a = 0.1-4, 1/τ_rec + r/τ_rel ≈ 0.3-5 s⁻¹ with
asymmetric cross couplings. TR's point: in THIS model x, a, b are all
dimensionless and O(1), so the full-J ω is not meaningless (it came out
50-180, same order as J_xx) — it measures whole-state instantaneous
growth in a stated norm, but is not invariant, so show it with the norm
stated and do not build a claim on it. The deeper limitation of BOTH is
that ω is instantaneous (t = 0⁺): adaptation's *dynamic* negative
feedback (a, b responding over 0.25-10 s) is invisible to it. The J_xx
figure therefore shows only the static (operating-point) part of
adaptation's effect and UNDERSTATES the reduction — conservative.

### 2.2 Agreed recommendations (TR agreed to all three)

1. **Transient gain, frozen vs active — the primary measure.**
   G(t) = ‖P_x e^{Jt} P_xᵀ‖₂: perturb only dendritic states, let the whole
   system respond, measure only on dendritic states. Scale-invariant
   (x-in, x-out, one unit), finite-time (no stability requirement),
   contains adaptation's dynamic feedback. Compute twice per sampled
   state: through J_xx (adaptation FROZEN = the conventional rate-network
   picture) and through the full J (adaptation ACTIVE); the gap between
   the two curves is adaptation's dynamic contribution and the figure
   that says adaptation belongs in the conceptual model of connectivity
   and dynamics (a thesis of the paper). Report peak gain G_max and peak
   time t_peak per regime, frozen and active.
2. Keep the J_xx ω/α figure as the "conventional view" (maps onto the
   literature; shows the operating-point effect).
3. Add the full-J ω as a supplementary panel with the norm stated.

**Direction of the perturbation — no choice needed.** G(t) is the worst
case over all dendritic directions; the optimal direction is the top
right singular vector of the x-to-x propagator at that t (differs across
t). Because propagating the n dendritic unit vectors gives the whole
n × n propagator, three further readings come free: the **noise-average
gain** (Frobenius/√n = mean-square amplification of isotropic
perturbations = what the model's own additive noise on x experiences),
the gain along the **E/I difference mode** (Murphy-Miller's prediction;
compare with worst case), and along the **leading Lyapunov direction**
(x part of `Q_final(:,1)`; how much of the available amplification the
dynamics actually recruit). Also report the optimal vector's E/I
structure and participation ratio at the peak.

**Cost.** Full J at N = 4000: propagate n = 500 columns with the
matrix-free `jacobian_times` (~30 ms per call at K = 500) → about 1 min
per state for a 2 s horizon at 400 Hz; 20 states per regime ≈ 1 h total.
J_xx (500 × 500 dense): `expm` at a few t, trivial. Horizon: 2 s first
(dendritic + fast STD timescales); the 10 s SFA rung needs a longer
horizon on fewer states if the curves have not settled.

**Caveat for all of it**: tangent-space linearisation along a noisy
trajectory → distributions over states; ignores saturation, so G is an
upper envelope of real amplification.

### 2.3 Kreiss constant — considered, left out

K(A) = sup_{Re z>0} Re z ‖(zI−A)⁻¹‖ is a certified LOWER bound on the peak
transient gain (K ≤ sup_t ‖e^{At}‖ ≤ e·N·K). Left out because: it needs a
stable matrix (J_xx has α = +7, +6 in two regimes; only the
multiple-timescale regime qualifies); it has the same state-scaling
dependence as ω on the full J (an x-in/x-out K would need modifying
Mitchell's inner resolvent-norm evaluation); on the full J (N = 4000)
Mitchell's dense code is hours per state whereas G(t) is a minute; and G
is the quantity itself, not a bound. Could return as a cross-check on
J_xx in the stable regime (K ≤ G_max; compare Kreiss's optimal direction
with G's top singular vector).

## 3. Where we are, and what to add next (agreed ordering)

**What we have.** Jacobian side: ω, α of J_xx along the trajectory
(static). Trajectory side, per sweep job: fraction of time the local
Lyapunov rate is positive, p95 of the 0.2 s finite-time exponent, mean
positive-excursion length, leading-vector block fractions. Spectrum side:
n₊, h_KS, D_KY per regime. Missing: the DYNAMIC picture — what a transient
does once adaptation responds, and whether the network's own transients
are the ones non-normality predicts.

**Next 1 — Transient-gain stage (build first).** New stage
`run_transient_gain` (model on `run_numerics_verification` /
`run_lyapunov_spectrum`: name-value cfg, cost table with `:badMode`,
`P` struct, `ensure_pool`, `results/cond_names/condition_titles/settings`
schema, register in `run_all_paper_analyses`, `test_run_modes`, the
configs) + `fig_transient_gain`. Per regime, ~20 states from the same
seeds as the eig stage (`build_from_preset`, `rng_seeds [1 2]`, sample
after the warm-up), horizon 2 s: G_worst(t), G_noise(t), G_EI(t),
G_Lyap(t), frozen (J_xx via `expm`) and active (full J via
`jacobian_times` on the n-column basis `[zeros; I_n]` in the x rows,
Heun as in `lyapunov_topk`, top singular value of the x rows every ~0.02 s),
optimal-vector E/I structure at the peak. Outputs: G_max, t_peak,
per regime × frozen/active; figure: three panels (regimes), two curves
each with seed/state bands, plus a table.

**Next 2 — Intermittency: are the excursions recruited transients?**
Same machinery, different sampling: in the single-timescale regime,
sample states at local-rate excursion ONSETS and in QUIET stretches (the
stored `local_rate_lead` / `t_lya_lead` mark them). Ask: is the optimal
transient direction at onset states aligned with the leading Lyapunov
direction, and is G(t) at onset states larger than at quiet ones? If yes,
the intermittent divergence is non-normal amplification triggered by the
network's own fluctuations — a mechanistic account of the Introduction's
"transient divergence" claim. Fold into the same stage as a sampling
rule.

**Next 3 — Stimulation that engages adaptation without a large transient
(TR's idea).** The manuscript's stimulation section engages adaptation
crudely (DC) on a limit-cycle bursting network. TR's proposal: on a
network that is asymptotically stable but has positive-LLE transients
(intermittent), find stimulation directions that are LEAST aligned with
transient amplification in the short run (avoid an immediate large
transient) yet strongly engage adaptation, so that adaptation reduces
transient amplification later in time. Linear formulation: for an input
direction v on x, short-time cost = ∫₀^τ ‖P_x e^{Jt} P_xᵀ v‖² dt
(τ ≈ a few hundred ms), long-time benefit = ‖P_a e^{Jt} P_xᵀ v‖ at t of a
few seconds (engagement of SFA/STD states); both are quadratic forms in
v, so the best direction is a generalised eigenvector of two n × n
matrices built from the same propagators as Next 1. Prior: the optimum
is near the E/I SUM mode (the difference mode is what amplifies; a
uniform push engages every neuron's adaptation) — which would put the
manuscript's DC stimulation on principled ground, or reveal a better
direction. Then a NONLINEAR check: stimulate the actual model along v
vs along the worst-case direction with a matched pulse and measure the
subsequent ω, the local-rate excursion rate and discharge rate.
Separate stage; likely a separate section or follow-on paper.

**Prerequisite for 2 and 3 — the right network.** `bursting_pairs` is a
limit cycle, not an intermittent stable system. The medium sweep
(`data/topk_med/1D_sensitivity_*`, `param_space_*`) stores per job λ₁,
`frac_local_positive`, `p95_finite_0p2s`, `mean_positive_excursion_s`:
QUERY it for grid points with λ₁ < 0, a high positive-time fraction and
long excursions (load with `ParamSpaceAnalysis2.from_dir`, filter results).
Do this before designing Next 3.

## 4. If resuming after a compaction: read these first

* This file; then commit messages `git log 5a9d409..3a64568`.
* `docs/notes/Non_normal_amplification.md` (§3 J_xx rationale, §5 the four
  quantification options; item 3 there = Next 1 above).
* `src/model/lyapunov/lyapunov_topk.m` (Heun loop with `opts.jac_times`;
  copy its propagation for the transient-gain stage),
  `SRNNCellTypePairs.jacobian_times` (~line 2141), `lya_summary`,
  `make_state_layout` (layout.x = last n indices; `.a{q}`, `.b{pre,post}`).
* `src/analysis/run_eig_heatmap.m` (`sample_eigenvalues`: how states are
  sampled, J_xx = J(end-n+1:end, end-n+1:end)), `run_lyapunov_spectrum.m`
  (stage template), `src/figures/fig_transient_amplification.m`.
* `scripts/paper/run_all_paper_analyses.m` stage table (~line 135),
  `scripts/tests/test_run_modes.m` (~line 69), `scripts/paper/topk_med_config.m`
  figure registry (`add(F, name, fn, in_paper, args)`).
* Figure contract: `arguments` with `data_file/out_dir/save/visible/run_dir/preset_name`;
  `setup_paths; default_out_dir; manuscript_style; resolve_data_file`;
  return `struct('figs', fig, 'files', {{}}, 'source', data_file)` and
  `save_figure_stable` + `existing_outputs`.
* House rules: MATLAB only via MCP (long runs: launch from the MATLAB
  prompt or capture with `evalc` to a log — the MCP client aborts after
  30 min of silence); `clear classes` after classdef edits (say so);
  stage by name; **no Co-Authored-By trailer in commit messages**
  (CLAUDE.md, updated); `figs/` and `data/` gitignored; note PDFs are
  force-added.


## 5. Addendum 2026-09-13 (after compaction): Next 1 and Next 2 built

*Same session (`90c5825b-...`), same machine; commits `478052a`, `71bf5a5`
and the docs commit after them. Working tree clean at the end.*

**Built.** `SRNNCellTypePairs.transient_gain` / `leading_direction_at` /
`excursion_samples` (+ `test_transient_gain`); stage
`src/analysis/run_transient_gain.m` (registered after `lyapunov_spectrum`
in `run_all_paper_analyses`, `test_run_modes`, and the figure registries
of `paper_config`, `single_multi_TS_independent_config`, `_med_config`,
`topk_smoke_fast_config`, `topk_med_config`); figures
`fig_transient_gain` (`Fig_Transient_Gain`) and
`fig_transient_gain_excursions` (`Fig_Transient_Gain_Excursions`), each
with a `_table.md`; `test_transient_gain_stage`;
`scripts/examples/find_intermittent_stable.m`; CLAUDE.md bullet;
`Non_normal_amplification.md` §6 (the results, read that first).

**Smoke** (`data/transient_gain_smoke`, `figs/transient_gain_smoke`, fast,
n = 500, 1 seed, 11 min): active 1-s G_max 1 200 / 30 / 8.6 for no /
single / multiple-timescale adaptation; freezing the full J instead of
J_xx alone cuts the frozen gain ~80× / ~6× in the adapted regimes (the
dynamic feedback). The frozen gain in a regime with α(J_xx) > 0 grows
without bound, so its "peak" is the horizon -- compare at fixed t. The
active gain is NOT bounded by the frozen ones (above them on the small
test nets and in the stable regime at n = 500). The worst-case direction
has |cos| ≤ 0.15 with the leading Lyapunov direction; the E/I difference
mode is the best named direction, the E/I sum mode the worst (< 1 when
adapted); the noise-average gain is ~1 in the adapted regimes.

**Next 2 status.** The sampling works but 20 s gives 0-2 onsets/quiets per
regime; the medium mode (T 40 s, 2 seeds, up to 8 + 8) is the first real
contrast. `topk_med_config` now includes the stage: rerunning
`topk_med_run` (delete `data/topk_med` first, or run the stage alone into
`data/topk_med/transient_gain` with `run_transient_gain('run_mode',
'medium', 'out_dir', ...)` -- about 1 h) gives it.

**Next 3 prerequisite done.** `find_intermittent_stable()` over
`data/topk_med`: 78 / 78 / 114 intermittent-but-stable jobs per regime;
best candidates the multiple-timescale regime at `mu_IE_relative` 2.06 or
`mu_EE_relative` 13.4-15.7 (1-D sensitivity, several reps, λ₁ ≈ −0.01 to
−0.1, p95 finite-time exponent 10-13 s⁻¹), and single-timescale at n = 190.

**Not done / raise only.** Next 3 itself (the generalised eigenproblem for
stimulation directions and the nonlinear check); Kreiss; CLVs; a noise-off
replicate of the transient-gain stage; a longer horizon for the 10-s SFA
rung. Also unexplained from §1.3: noise RAISES h_KS in the no-adaptation
regime.
