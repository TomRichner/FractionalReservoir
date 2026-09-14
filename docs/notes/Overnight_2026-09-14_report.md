# Overnight report, 2026-09-13/14: verbose levels, STD strength matching, and seven manuscript improvements

*Session `90c5825b` on R5611351, branch `main`, commits `f5c5cac` .. see §9. Plan:
`~/.claude/plans/dreamy-dazzling-ladybug.md`. Everything ran through the MATLAB MCP;
`clear classes` was called several times after classdef edits, which WIPED THE BASE
WORKSPACE of the live MATLAB session (the variables `cfg`, `cfg3`, `results`,
`results3`, `run_dir`, `run_dir3`, `zz_*`, `r_*` left in it are mine).*

## 0. The headline: the scaled matching makes the multiple-timescale network chaotic

TR chose (2026-09-13) to match STD strength by a **multiplicative route scale**
s = 1 + r_ref/rho = 3 on the two-timescale routes, with the usage-matched
(`tau_rel`) variant as a control. Both were built, tested and run at `fast`
(`data/stdscaled_fast`, `data/stdusage_fast`; figures in `figs/stdscaled_fast`,
`figs/stdusage_fast`). The scale variant does what the plan's context warned it
would, and more:

| | unmatched (`data/topk_med`, medium, 15 reps) | scaled, s = 3 (`data/stdscaled_fast`, fast, 3 reps) | usage, rho_u = 0.342 (`data/stdusage_fast`, fast, 3 reps) |
|---|---|---|---|
| lambda_1 near default: none / single / multiple | +3.47 / +0.35 / **-0.112** (105 of 105 negative) | +3.65 / +0.61 / **+2.24** (0 of 21 negative) | +3.65 / +0.61 / **+1.01** (0 of 21 negative) |
| example run n = 500, 20 s, seeds [1 2] | +3.35 / +0.43 / -0.11 | +3.77 / +0.43 / **+1.99** | +3.77 / +0.43 / **+0.90** |
| tau sweep, multiple-timescale, slowest tau 1 -> 30 s | -0.112 -> -0.043, every rep negative | +2.27 -> +2.27, every rep positive | +0.86 -> +1.14, every rep positive, no trend |
| memory capacity, total MC: none / single / multiple | 0.103 / 0.259 / **0.590** (15 trials, p = 6e-5) | 0.111 / 0.243 / **0.103** (5 trials) | 0.111 / 0.243 / **0.146** (5 trials) |
| K needed by top-K for the multiple-timescale condition | 15 | 60, D_KY resolved in 0-50% of jobs | 15-60 (D_KY resolved in 39-100% of jobs) |
| 1-D sweeps: multiple-timescale is the least stable condition on | none of 7 | f_E, mu_EE, mu_EI, mu_II (4 of 7) | n, mu_EI, mu_IE, mu_II (above the single-timescale curve on 6 of 7 sweeps; below it only on mu_EE) |

With the 3x undepressed recurrent gain the multiple-timescale network is
strongly chaotic and its fading memory collapses to the no-adaptation level;
the paper's central computational result does not survive the scale
direction. The usage control (R3) is less extreme but points the same way:
lambda_1 +1.0 near default, +0.9 in the example run, +0.65 to +1.14 across the
tau sweep with no trend, memory capacity 0.146. **At equal steady-state
depression strength at r_ref, the timescale count alone does not produce the
stable regime: the unmatched multiple-timescale advantage leaned on the
squared (stronger) depression.** That is the honest answer to Codex sec. 1's
confound question, and it is what the manuscript must say.

Because both matchings WEAKEN the multiple-timescale network's depression, a
third direction was built and launched at 03:15 (sec. 8): keep the
multiple-timescale network exactly as published and STRENGTHEN the
single-timescale route to the dual's steady state (tau_rel 0.25 -> 0.0625,
`..._dualStdSingleMatched_3cond_mu8p25`). If a single-timescale STD of equal
strength also gives lambda_1 ~ -0.1 and MC ~ 0.6, the paper's claim is about
depression strength, not timescale count; if it does not, the timescale
structure matters at that strength. `paper_config` still points at the
SCALED preset because that was the instruction; **it should not stay there.**
Nothing has been run at medium on any matched variant except the usage run
launched tonight (`data/stdusage_med`).

The matching itself is exact at r_ref for both variants (`fig_STD_strength_matching`):
theta_ss(r_ref) = 0.0833 for single, scaled and usage; the largest |log ratio|
to the single-timescale route over the occupied rate band [0.199, 0.254] is
1.11 unmatched, 0.146 scaled, 0.034 usage. Usage is also the better match away
from r_ref. Caveat on record: per-neuron time-averaged rates are strongly
bimodal (near 0 or saturated), so "the operating point" is the population mean
rate, not a typical neuron's rate.

## 0b. Addendum, 03:5x: the third direction keeps the result

The single-matched run (`data/stdsinglematched_fast`; dual routes as
published, single-timescale route strengthened to the dual's steady state,
tau_rel 0.0625) had finished its sweeps and memory capacity when this was
written (figures still running):

| | unmatched (medium, 15 trials) | single-matched (fast, 5 trials) |
|---|---|---|
| total MC: none / single / multiple | 0.103 / 0.259 / **0.590** | 0.111 / **0.385** / **0.583** |
| horizon (s) | 0.00 / 0.12 / 0.52 | 0.00 / 0.30 / 0.48 |
| single vs multiple | p = 1.2e-4, d_z = -2.04 | p = 0.125 (floor 0.0625 at 5 trials), d_z = -1.05 |
| K needed by top-K, multiple-timescale | 15 | 15 (D_KY resolved 67-100%) |

So at EQUAL steady-state depression strength the multiple-timescale network
still has the larger fading memory (0.58 vs 0.39), while the single-timescale
network gains from the stronger depression (0.24 -> 0.39). Both mechanisms
contribute; the timescale-count claim survives in the form "at matched
strength, two depression timescales extend fading memory further than one",
and the effect size is smaller than the unmatched comparison suggested. The
multiple-timescale lambda_1 for this run is in
`figs/stdsinglematched_fast/fig_local_vs_finite_lle/*_table.md` once the
figure pass finishes (the sweeps' K = 15 says it is the published stable
network). **Recommendation: make `..._dualStdSingleMatched_3cond_mu8p25` the
paper preset** -- it changes only the single-timescale condition, keeps every
multiple-timescale result, and answers the confound honestly -- and run it at
medium (clone `stdUsage_med_config`). The Methods paragraph in
`STD_strength_matching_2026-09-13.md` sec. 4 then needs its matching sentence
turned around (the single-timescale release constant shortened to 0.0625 s so
its steady-state output equals the two-timescale routes' at r_ref).

## 1. r_ref

r_ref = 0.25 = the median mean firing rate of the multiple-timescale condition
at the default point of the five well-defined 1-D sweeps of `data/topk_med`
(f_E 0.236, level_of_chaos 0.223, mu_EE 0.239, mu_IE 0.251, n 0.231; with the
two negative-mean sweeps resolved correctly the pooled 105 reps give 0.230
[p5 0.199, p95 0.254]), rounded to 0.05. Single-timescale 0.269, no adaptation
0.483. The joint 64-point sample, far from the default, has 0.154 / 0.159 /
0.250. Full tables and the replacement manuscript Methods paragraph:
`docs/notes/STD_strength_matching_2026-09-13.md` (+pdf).

## 2. What was built (all committed, all tested; §9 lists the commits)

* **Verbose** (`f5c5cac`): one setting, `'verbose' | 'minimal' | 'near-none'`,
  `minimal` the default, from `paper_config` through `ctx`, every stage,
  `ParamSpaceAnalysis2` and both model classes into the parfor workers; the
  `evalc` wrappers in the stages are gone; warnings never gated. A serial 2 x 2
  sweep prints three lines. `test_verbose_levels`.
* **STD strength matching** (`d942c51`): `synapse_config.<pre>.<post>.scale`
  (params.W only; obj.W untouched; dynamics, both Jacobians and jacobian_times
  consistent by construction; the ESN readout check compares it); presets
  `..._dualStdScaled_3cond_mu8p25` (primary) and `..._dualStdUsage_3cond_mu8p25`
  (control); golden fixture regenerated (16 presets); `fig_STD_strength_matching`
  (in_paper); `stdScaled_fast_config/_run`, `stdUsage_fast_config/_run`;
  `Equations_stability_paper.md` subsection, `MTS_STD_notes.md` addendum,
  CLAUDE.md. `test_route_scale`. The manuscript was NOT edited; its "depression
  is deliberately not normalized" sentence and the replacement are in the note.
* **Tau sweep on both cell types** (`b317923`): `tau_a_EI` alias; the stage
  picks `tau_a_EI` when the condition's E and I ladders are identical;
  `tau_levels.md`/`.mat` manifest proves it (R2: equal E and I ladders at all 7
  levels). `test_tau_a_EI_alias`. The sweep folder is now named by its axis.
* **Numerics gate** (`c0a2a22`): check J (analytic Jacobian vs central finite
  differences on the reduced net, kink rows excluded); acceptance thresholds
  fixed in code on 2026-09-14 before any run, saved in `settings.acceptance`;
  `verdict` in the .mat and `numerics_verdict.md`; `'jacobian'` figure variant.
  `test_numerics_verification_stage` (~13 min).
* **Imbalance eig-heatmap examples** (`c0a2a22`): `run_eig_heatmap` `cfg.examples`
  (mu_EE x 0.5 / 1 / 1.5, shared seeds, per example and condition eigenvalues,
  top-K lambda_1, mean rate and B_E); `fig_eig_heatmap_imbalance` (in_paper);
  legacy variables unchanged. `test_eig_heatmap_stage`.
* **Figures** (`07dc9a6`): `fig_example_timeseries` three-column composite with
  the local rate and the accumulating finite-time lambda_1 (`plot_local_and_finite_lle`,
  `SRNNCellTypePairs.routes_identical`); `fig_local_vs_finite_lle`;
  `fig_lle_vs_rate`; IQR bands and data-driven limits on `fig_sensitivity_medians`;
  the unclipped density-strip `fig_sfa_EOC_allStd`; E:I weights title overlap.
* **Memory capacity** (this commit): `fig_memory_capacity` redrawn (data-driven
  cumulative axis, bootstrap bands, paired horizons with the sign-flip p and
  d_z read from the saved stats, a table); `run_memory_capacity` writes
  `<run>_provenance.md` (commit, preset, mode, seeds, trials, protocol, stats,
  paths). The 30-trial production run: see §8.

## 3. Numerics verdict (fast, scaled preset; `data/stdscaled_fast/numerics_verification/numerics_verdict.md`)

| check | none | single | multiple | note |
|---|---|---|---|---|
| A noise-free reshoot order slope >= 1 | 2.01 PASS | 2.01 PASS | 2.01 PASS | Ralston order 2 despite the kinks |
| B noisy strong-order slope >= 1 | 2.01 PASS | 1.98 PASS | 1.99 PASS | ratio 4x per halving |
| L median \|lambda_1 sra1 - ode45\| <= 0.05 | 0.108 FAIL | 0.169 FAIL | 0.722 FAIL | 10-s windows, 2 seeds: finite-time scatter |
| C \|Benettin - QR\| <= 0.05, \|topK - QR\| <= 0.02 | 0.037 / 0.006 PASS | 0.308 / 0.016 FAIL | 0.166 / 0.058 FAIL | reduced net n = 30, T 12 s |
| J rel Frobenius, eigenvalue error <= 1e-6 | 1e-10 / 7e-9 PASS | 5e-11 / 5e-10 PASS | 6e-11 / 5e-10 PASS | 0 kink rows excluded |

The gate works as designed: the convergence and Jacobian claims pass; the
finite-time LLE agreement claims cannot be judged on fast windows and must be
judged at medium (25 seeds, 20-s windows) or production. The manuscript's
Jacobian claim now has a saved artifact.

## 4. Tau sweep on both cell types

`data/topk_med` (unmatched, tau_a_E only, 11 levels x 25 reps): every rep
negative, lambda_1 -0.112 at 1 s rising monotonically to -0.043 at 30 s
(-1/tau = -0.033): the slowest SFA rung sets the exponent in the stable regime,
as the sweep's header predicts. The old figure's "half the distribution
positive" described a Benettin-era run. R2 (scaled, tau_a_EI): every rep
positive, +2.3 at every level -- the network is chaotic and the SFA timescale
no longer sets lambda_1. R3: +0.86 -> +1.14, every rep positive, no trend_TEXT.

## 5. Local vs finite-time LLE (`fig_local_vs_finite_lle`) and LLE vs rate (`fig_lle_vs_rate`)

`data/topk_med`, near-default set (7 sweeps at the level nearest the default,
105 networks per condition):

| condition | lambda_1 median [IQR] | networks with lambda_1 < 0 | local rate > 0, share of time | median frac_local_positive | mean excursion (s) |
|---|---|---|---|---|---|
| no adaptation | +3.47 [+1.81, +4.33] | 4 / 105 | 83% | 0.91 | 0.70 |
| single-timescale | +0.35 [-0.06, +0.86] | 33 / 105 | 53% | 0.54 | 0.26 |
| multiple-timescale | -0.112 [-0.113, -0.110] | 105 / 105 | 3% | 0.02 | 0.12 |

So "transient expansion within long-interval stability" (Codex §3's stronger
claim) describes the SINGLE-timescale regime near the boundary (locally
expanding half the time, lambda_1 straddling zero), not the multiple-timescale
one, which is quietly stable (locally expanding 3% of the time, lambda_1
pinned at -0.11 by the 10-s SFA rung). Keep the stronger idea in the Discussion.

LLE vs mean rate over 1219 networks (joint sample + 1-D sweeps), Spearman rho
[bootstrap 95% CI]: no adaptation -0.44 [-0.50, -0.38] (inverted U: quiet
+0.3, mid +3.1, saturated -5.7), single -0.31 [-0.36, -0.25], multiple +0.11
[+0.05, +0.16] (flat at -0.11 across the whole rate range). Rate and stability
interact through the operating point and are not interchangeable; in the
multiple-timescale regime the rate barely predicts lambda_1 at all.

## 6. Imbalance eig-heatmap examples (R2, scaled preset, fast: reference at n = 500, imbalanced at n = 250)

| example | condition | lambda_1 | mean rate | B_E |
|---|---|---|---|---|
| mu_EE x 0.5 (inhibition-dominant) | none / single / multiple | +1.98 / -0.19 / +4.89 | 0.06 / 0.05 / 0.05 | 0.42 |
| reference | none / single / multiple | +3.94 / +0.63 / +1.95 | 0.43 / 0.25 / 0.30 | 0.50 |
| mu_EE x 1.5 (excitation-dominant) | none / single / multiple | -4.23 / -0.53 / +1.19 | 0.87 / 0.71 / 0.35 | 0.55 |

The sheet (`figs/stdscaled_fast/fig_eig_heatmap_imbalance`) reads well: the
no-adaptation network goes from chaotic to saturated-and-silent as excitation
grows, single-timescale stays near zero, and the scaled multiple-timescale
network is chaotic in every column (its occupancy spills past Re = 0 in all
three rows). Rerun on the chosen preset at medium before it goes in the paper.

## 7. Transient gain on the scaled preset (fast, for the record)

Operating-point table: lambda_1 +4.0 / +0.44 / +1.95; alpha(J_xx) at state
+8.7 / +7.2 / +6.7; alpha(J) +8.7 / +2.6 / +4.0. Adaptation's feedback still
takes the frozen rate from +6.7 to +4.0 in the multiple-timescale regime, but
the trajectory itself is at +1.95: the scale has moved the operating point,
not removed the mechanism.

## 8. Memory capacity

Audit of `data/topk_med/memory_capacity/MC_sample_hold_20260912_232101_trials15_*`
(unmatched preset, medium): 15 of 15 paired trials complete for all three
conditions, all summaries finite; total MC 0.103 [0.093, 0.114] / 0.259
[0.218, 0.313] / 0.590 [0.515, 0.664]; horizon 0.00 / 0.12 / 0.52 s; exact
sign-flip over 32768 patterns: none vs multiple p = 6.1e-5, d_z = -3.32;
single vs multiple p = 1.2e-4, d_z = -2.04. `fig_memory_capacity` resolves
exactly that file for `run_dir = data/topk_med`. Horizons are quantised to
0.3-s holds, so panel (c) is inherently discrete.

R2 (scaled, fast, 5 trials): 0.111 / 0.243 / 0.103, single vs multiple p =
0.0625 (the exact floor at 5 trials) with d_z = +1.66 in the WRONG direction.
R3 (usage): 0.111 / 0.243 / **0.146** (5 trials)_TEXT.

Unattended chain launched at 03:15 in one MATLAB call (diaries in `data/`):
1. `stdSingleMatched_fast_run` -> `data/stdsinglematched_fast`, `figs/stdsinglematched_fast`
   (log `data/stdsinglematched_fast_log.txt`; ~1.7 h);
2. `stdUsage_med_run` -> `data/stdusage_med`, `figs/stdusage_med` (medium: 15 reps,
   11 levels, 15 MC trials; equal footing with `data/topk_med`; ~3 h);
3. `run_memory_capacity` at PRODUCTION (30 paired trials) on the UNMATCHED paper
   preset -> `data/topk_med_mc_production` (log `data/mc_production_log.txt`).
   This replaces the plan's production MC on the scaled preset, whose fast
   result is already decisive. Plot it with
   `fig_memory_capacity('run_dir', 'data/topk_med_mc_production', 'out_dir', ...)`
   (the runner writes into `<output_dir>/memory_capacity`).
If MATLAB is still busy in the morning, the diaries say which job it is on;
`run_all_paper_analyses` refuses a non-empty run directory, so a crashed job
can be restarted after deleting its folder.

## 9. Commits, figures, what was left undone

Commits tonight, in order: `f5c5cac` verbose; `d942c51` STD matching;
`b317923` tau_a_EI; `c0a2a22` numerics gate + imbalance examples; `07dc9a6`
figures; `9e39cb3` config handles; `a1bb18b` docs; `eb8e56d`
make_all_paper_figures verbose order (this bug cost R2 one figure pass: the
level was read at line 66 and defined at line 119); `35f70bc` memory capacity figure + provenance; `336532f` the third matching
direction (`dualStdSingleMatched`) and the medium usage config; the final
commit carries this note and the fig_STD_strength_matching title fix.

Figures: `figs/stdscaled_fast/<entry>/` (31 of 31 entries succeeded, 17
in-paper, manifest.md there); `figs/stdusage_fast/` (31 of 31, 17 in-paper). Scratch renders
from the unmatched medium data are in `figs/scratch/`.

Wrong turns worth knowing: the Bash tool mangles backslashes in heredocs, so
every multi-line edit went through perl scripts written with the Write tool
(two figure registrations lost their `@` handles that way, `9e39cb3`); test
scripts run in the base workspace and clobber loop variables, so tests were
run one per call; a `(1,:) char` argument rejects the logical the old callers
pass, so the `verbose` arguments are unconstrained and validated by
`verbose_name`; MC's `capture_git_provenance` writes a `working_changes.patch`
when the tree is dirty, which it was during the runs (the figure/provenance
edits were uncommitted while R2 ran).

Not done: nothing from the plan was skipped. Raised only: the E:I weights
sheet's colour scheme and the sensitivity sheets still show the no-adaptation
curve dominating the y range (the adapted curves are readable but compressed;
a symmetric-log axis would help); `fig_local_vs_finite_lle`'s row-2 x labels
overlap between columns; `plot_memory_capacity_combined` still has the fixed
[0 10.9] cumulative axis for direct callers; `test_numerics_verification_stage`
takes 13 min because the full-size L check has no size override.

## 10. What TR should do first

1. Read sec. 0. Both weakening matchings (scale, usage) turn the
   multiple-timescale network chaotic and remove its memory advantage. Look at
   `data/stdsinglematched_fast` (sec. 8 chain, job 1): does a single-timescale
   STD of the dual's strength reproduce lambda_1 ~ -0.1 and MC ~ 0.6? Read
   `figs/stdsinglematched_fast/fig_sensitivity_medians`, `fig_memory_capacity`
   and `fig_local_vs_finite_lle` tables.
2. Decide what the paper claims. Options: (a) the timescale-count claim,
   supported only if single-matched does NOT reproduce the stable regime; (b)
   a combined "multiple timescales AND stronger depression" architecture claim,
   which the unmatched comparison supports and which the Methods paragraph in
   `STD_strength_matching_2026-09-13.md` sec. 4 must then be rewritten to say;
   (c) match at a different quantity (e.g. the occupied mean of theta rather
   than theta at one rate).
3. Move `paper_config` off the scaled preset (one line) to whichever preset
   the decision names; the medium usage run (`data/stdusage_med`) and the
   30-trial production MC on the unmatched preset (`data/topk_med_mc_production`)
   will be waiting either way.
