# Issues

Bugs and observed problems, tracked **by issue** rather than chronologically —
`WorkLog.md` is the chronological record and this is the index into it.

Every issue records when it was identified, on which branch and commit, by which
host, and how it ended: **fixed**, **won't fix** (a deliberate decision), or
still **open**. Issues are never deleted; a resolved one keeps its entry with the
resolving commit, so the same ground is not re-covered in six months.

Newest first. Every heading carries a **plain-text status token** next to the
marker, so the backlog can be filtered without depending on emoji surviving a
shell's locale (`grep '^## 🔴'` silently matches nothing in some terminals):

| | token | meaning |
|---|---|---|
| 🔴 | `OPEN` | still outstanding |
| 🟢 | `FIXED` | resolved; entry kept with the resolving commit |
| ⚪ | `WONTFIX` | deliberate decision not to change anything |
| 🔵 | `OBSERVATION` | not a defect, but affects interpretation |

```bash
grep '^## ' Issues.md              # whole backlog, one line each
grep '^## .* OPEN ·' Issues.md     # just what is outstanding
```

---

## ⚪ WONTFIX · ISSUE-013 · Wide PNG exports drop glyphs from tick labels

| | |
|---|---|
| Identified | 2026-08-18 21:37 · `dev` @ `fc6af94` · R5611351 · Claude Code (Opus 5), session b651fa7a |
| Closed | 2026-08-19 · WONTFIX by TR |
| Area | `src/plotting/plot_saving/save_some_figs_to_folder_2.m` (the png branch) |

**WONTFIX by TR 2026-08-19.** The defect is real and is diagnosed below, but it
is a MATLAB rendering bug, it only touches PNGs (the `.svg` is clean), and the
one remedy tried cost more than the defect — see "Attempted fix, reverted". The
saver stays at a plain 600 dpi `exportgraphics`. Do not reopen or "fix" this
without asking; the diagnosis is recorded so the next person does not re-derive
it, not as an invitation to act on it.

Building `fig_sensitivity_medians`, x-tick labels in the 600-dpi PNG came out
with characters missing — `700` as `)0`, `+50%` as `+%`, later `500` as `50`,
`1000` as `10`, `+100%` as `)0%`. One label per sheet, which reads as a typo in
the tick rather than as a broken export, and is easy to ship without noticing.

**Not a layout collision and not a data problem.** Reading the axes back live,
`ax.XTickLabel` is `{'100','400','700','1000'}` and `XTickLabelRotation` is 0;
the labels are correct in the figure. The corruption is introduced by the
raster export:

- Re-exporting the *same live figure* at `'Resolution', 300` renders `700`
  perfectly.
- Re-exporting at 600 reproduces the identical corruption, byte-for-byte in the
  same place — so it is deterministic, not the intermittent
  `MATLAB:class:InvalidHandle` export failure that `save_some_figs_to_folder_2`
  already documents.
- Both corrupted labels sat at roughly the same **horizontal pixel position**
  (~x = 3750 of a 7027-px-wide raster) in different panels, which points at a
  tile seam in the high-resolution rasteriser rather than at anything about the
  glyphs.
- The vector `.svg` from the same figure is unaffected.

**First workaround, WRONG, recorded so it is not retried:** widening the figure
from 1200 to 1300 made the bad label render. It only slid a *different* label
onto the seam — the next regeneration ate `500` instead. Any "fix" that changes
the figure's size is this same coincidence. Likewise wrong: forcing
`Renderer = 'painters'` (no effect), and `print -dpng -r600` instead of
`exportgraphics` (eats the same label in the same place). It is not resolution
*per se* either — it is the resulting pixel width.

**Diagnosis.** The damage always lands near the end of a rasteriser tile, so
which label is hit moves with the image width. Measured on these figures at
2075 / 3520 / 3759 / 4001 / 5999 / 7639 / 8125 px. Bisected on this machine:
**3520 px clean, 3759 px damaged.**

**Attempted fix, REVERTED by TR 2026-08-19.** The PNG branch was made to measure
the file it had just written with `imfinfo` and, when the larger dimension
exceeded `png_max_px = 3400`, re-export at the highest resolution that fits. It
did stop the glyph loss, but the cost is the reason it was rejected: it silently
drops any wide sheet from 600 dpi to roughly 290, it does so for **every** script
in the repo that saves figures, and `png_max_px` is an empirical bound measured
on two figures on one machine rather than a documented limit. Trading the
resolution of every figure against a rendering defect in a few is not the right
bargain, and it was made without asking. The saver is back to a plain 600 dpi
`exportgraphics`.

**If it ever does need addressing**, the untried options are: rasterising the
PNG from the vector path rather than the raster one; treating the `.svg` as the
only trustworthy output for wide sheets; or reporting it to MathWorks, since
nothing about it looks like project code. None of these are to be attempted
without TR asking for them.

**Practical consequence to know about:** a tick label on a wide sheet can lose a
character or two in the PNG. If a figure looks like it has a typo in an axis
tick, check the `.svg` before believing it.

---

## 🟢 FIXED · ISSUE-012 · Default-value marker never appears on the `_allStd` sensitivity figures

| | |
|---|---|
| Identified | 2026-08-14 · `dev` @ `b880c41` · R5456622 · Claude Code (Opus 5), session 054a29ca |
| Closed | 2026-08-19 · confirmed fixed by TR |
| Area | `scripts/presentations/Stability_Manuscript/fig_sensitivity_analysis_allStd{,2}/` |

**FIXED — confirmed by TR 2026-08-19.** The marker draws. The narrative below is
the original report and the intermediate evidence; the resolution is at the
bottom.

`Fig_sensitivity_analysis_allStd.m` is supposed to draw a short reddish-gray
tick rising off the x-axis at the preset's default for each row (`:244-266`,
helper `mark_default_value` at `:650`). **It does not appear.** Every other
change in the same edit landed correctly — percent axes, `clim_frac = 0.6`
(`CLim = [0 9]` on a 15-rep run), horizontal tick labels, the dropped
`(% of default)` suffix — so the figures are usable; only the marker is missing.

Inspecting the saved `Fig_Sensitivity_LLE_core.fig`, each axes contains exactly
**one** line object: the blue median (`n=11`, spanning the full x-range,
LineWidth 3). No second, short line exists at all.

**What has been ruled out:**

- Not the in-range guard. All seven defaults were verified independently to sit
  inside their swept ranges (`f_E` 0.5 in [0.2, 0.8], gain 1 in [0.5, 1.5],
  `n` 500 in [100, 1000], `mu_*` ±5.5 in ±[2.75, 11]) for both runs, so
  `mark_default_value` should not be returning early.
- Not a failed default lookup — the run emitted no `NoDefault`,
  `PresetProbeFailed` or `NoManifest` warning, so `preset_default_values`
  apparently resolved values without complaint.
- Not the median-line restyle overwriting it: the marker is drawn *after* that
  block, deliberately, since the restyle selects on `Type 'line'`.

**Not yet checked:** whether `default_value` is genuinely non-empty at runtime
(the map is built before the sheet loop; `isKey` may be failing on the exact
key strings), whether `mark_default_value` is reached at all, and whether the
later axis-repositioning / `AddLetters2Plots` passes discard added children.

Affects `_allStd2` identically — same code, regenerated from the same source —
though `_allStd2` has not been re-run since the marker was added.

**Deferred by TR 2026-08-14**: the figures are otherwise correct and were
needed; come back to this rather than chase it inline.

**2026-08-18 21:37 · `dev` @ `fc6af94` · R5611351 · session b651fa7a — probably
already resolved, but not re-verified from the script.** Two data points found
while building `fig_sensitivity_medians`:

1. The **currently committed** `fig_sensitivity_analysis_allStd/Fig_Sensitivity_LLE_core.png`
   and `_mu.png` (regenerated 2026-08-18 19:38, after this issue was filed) DO
   show the reddish tick at the bottom of every panel.
2. `mark_default_value` and `preset_default_values`, extracted verbatim into
   `scripts/run_all_analyses/replot/`, draw the marker correctly on all six
   panels of both new median sheets.

So the helper works and the lookup resolves. What is not established is whether
the `_allStd` script itself was ever broken or whether the original report
predates a fix; `_allStd` was NOT re-run in this session (deliberately, to avoid
churning its committed figures). Close this by re-running that script and
confirming, not on the strength of these two observations alone.

**RESOLVED 2026-08-19 · TR.** `Fig_sensitivity_analysis_allStd.m` has since been
re-run several times and the reddish default tick is present on every panel of
all four sheets. TR confirms this is fixed. Whether the original report described
a real defect that was incidentally repaired, or a figure that predated the
marker being added, is not worth establishing now — the current output is
correct. `_allStd2` has still not been re-run, but it is the same code.

---

## 🟢 FIXED · ISSUE-011 · `load_results` silently returns level indices for vector parameters

| | |
|---|---|
| Identified | 2026-08-14 · `dev` @ `d58b7fe` · R5456622 |
| Area | `src/ParamSpaceAnalysis2.m` |
| Fixed | 2026-08-14 · `dev` @ `b6e5262` (+ test `b445e9b`) · R5456622 |

`param_space_summary.mat` never stores `vector_param_lookup` or
`vector_param_config` (see the save list at `:2489-2510`). After
`load_results`, both are empty, so `effective_param` (`:464`) falls through and
returns `res.config.<name>` — the **level index**.

Measured on the tau run: `effective_param(res,'tau_a_E')` returns `1` where the
answer is `5 s`. **No error is raised** — indices 1…13 quietly stand in for
5…60 s. The plotting path (`:936`) takes the same wrong branch and produces a
wrong x-axis.

**Workaround:** load `psa_object.mat` directly. `saveobj`/`loadobj` round-trip
both fields correctly — verified, it returns the real
`tau_a_E(end) = [5, 9.58, …, 55.4, 60]`.

**Blast radius — smaller than first recorded.** An earlier version of this entry
claimed `load_results` is "the path the replot scripts actually use". That is
**wrong**; checked 2026-08-14:

| caller | path |
|---|---|
| all `replot_*.m` | `psa_object.mat` only, skips dirs without it ✅ |
| `combine_runs.m` (pooling) | `psa_object.mat` + `same_config` ✅ |
| `load_and_make_unit_histograms.m` | prefers `psa_object.mat`, falls back to `load_results` |
| `Fig_2_fraction_excitatory_load_and_plot.m` | errors if `psa_object.mat` absent |

`load_results` also only *sets* six fields without clearing anything, so calling
it on an object already restored from `psa_object.mat` is harmless. The bug
therefore bites in exactly one case: **`psa_object.mat` is missing**, i.e. a
crashed or interrupted run — which is precisely when salvaging the data matters,
and precisely how it was found.

**Why the tests missed it:** `test_psa_saveload` exercises `saveobj`/`loadobj`
(the path that works) and never calls `load_results` at all.

**Shared root with ISSUE-010:** `saveobj`/`loadobj` and `load_results` have
diverged — the former persists the whole object, the latter restores a
hand-picked six fields from a summary file that stores its own different
hand-picked set. Two lists, drifting, with no test asserting they agree. Fix
both together, with a test that the two load paths yield equivalent objects.

### Resolution

Fixed by removing the second load path rather than repairing it. `run()` now
writes `psa_object.mat` itself — once before batching, once on completion — and
`ParamSpaceAnalysis2.from_dir` is the only supported way to read a run back.
The config-restore blocks in `load_results` and `consolidate` are gone;
`load_results` is now the private, results-only `load_condition_results`, and
`param_space_summary.mat` is documented as metadata, not a restore path.

`effective_param` now **errors** (`ParamSpaceAnalysis2:MissingVectorLookup`)
when a parameter is registered in `vector_param_config` but absent from
`vector_param_lookup`, instead of falling through to the scalar branch. That
guard is what makes this class of bug loud rather than silent.

Verified on the data that exposed it — `from_dir` on
`run_all_aug_13_26_23_37/FAILED_OOM_tau_sensitivity_*` returns
`effective_param(res,'tau_a_E') = [0.25 1.5478 9.5833]`, where it previously
returned `1`, and recovers 107/195 results from the partial run. All seven
completed sweeps in that directory replot cleanly through the migrated reader.

**Backward compatibility deliberately dropped:** a run with a summary but no
`psa_object.mat` no longer loads at all. `FAILED_OOM_param_space_*` from
2026-08-14 is the one such casualty on disk, and it now fails with a message
naming `consolidate()` rather than loading wrongly. Git provenance covers
re-running old data with old code.

**Test:** `scripts/tests/test_psa_loaders.m`, which sweeps a vector parameter
alongside a scalar one — the vector is the only kind that exposed this.

---

## 🟢 FIXED · ISSUE-010 · `load_results` does not restore `model_class`

| | |
|---|---|
| Identified | 2026-08-14 · `dev` @ `d58b7fe` · R5456622 |
| Area | `src/ParamSpaceAnalysis2.m:1612-1621` |
| Fixed | 2026-08-14 · `dev` @ `b6e5262` (+ test `b445e9b`) · R5456622 |

`summary_data.model_class` **is** saved (`:2495`) but `load_results` restores
only six fields and this is not one of them. A loaded `SRNNCellTypePairs` run
reports `SRNNModel2`.

Consequences: `effective_param`'s class-default fallback consults the wrong
class, and `same_config`'s documented refusal to pool runs from different model
classes cannot work on loaded runs.

Same root cause as ISSUE-011 (an incomplete `load_results`) and same test blind
spot; worth fixing together.

### Resolution

Fixed with ISSUE-011 — see that entry. The hand-maintained restore list is gone
entirely, so there is no longer a field that can be forgotten: `from_dir`
returns the saved object itself. Confirmed on last night's sweeps, which now
report `model_class = SRNNCellTypePairs` where they previously reported
`SRNNModel2`.

---

## ⚪ WONTFIX · ISSUE-009 · `medium2` tau/param_space stages died — cause NOT established

| | |
|---|---|
| Identified | 2026-08-14 · `dev` @ `ac59b42` · R5456622 |
| Closed | 2026-08-19 · WONTFIX by TR |
| Area | run infrastructure / MATLAB stability |
| Mitigated by | `c6055ea` (insurance, not a fix) |

**WONTFIX by TR 2026-08-19.** Stale. Several good production runs have completed
since — including `run_all_aug_14_26_17_25`, which the manuscript figures are
built from — so the crash never became a pattern. Closed, not parked.

**Original disposition (2026-08-14), kept for the record.** This was the
**first crash in a big run, ever**. A one-off with no reproduction is not worth
open-ended debugging, and the diagnostic cost is a whole night per attempt. It
stops blocking anything.

**Reopen at normal priority if it happens a second time** — a repeat makes it a
pattern rather than an incident, and the evidence from two failures is worth far
more than more speculation about one. If that happens, FR-006's memory
pre-flight is the instrumentation to add first.

Overnight `medium2` run: the 7 sensitivity sweeps completed 624/624 clean, then
`tau_sensitivity` returned `Out of memory` on **88 of 195** runs and
`param_space` took the MATLAB process down. Crash dumps at `05:41` (a worker —
the run continued) and `06:21` (the client).

**Two diagnoses were proposed and both are wrong**, recorded so they are not
retried:

1. *Total RAM exhausted with 14 workers* — the machine has 6 physical cores and
   the pool was 10.
2. *Heap fragmentation* — a 32-bit-era failure mode; on 64-bit the address space
   is ~128 TB and contiguity essentially never binds.

The data rules out resource exhaustion: commit limit 128.7 GB, page-file **peak
7.8 GB** over 10 days including the run, `MaxPossibleArrayBytes` 99.1 GB, no
array-size-limit preference, estimated peak ~14 GB.

**Leading hypothesis:** a worker process crashed and its jobs surfaced as memory
errors — making `Out of memory` a symptom, not the cause. Unproven.

**Note the apparent fix proves nothing.** The re-run that went 25/25 changed
**three** things at once:

| | overnight (failed) | re-run (25/25) |
|---|---|---|
| pool age | ~6 h old | fresh |
| workers | **10** | **5** |
| profile | `ten_workers_one_thread` (default) | `Processes` (built-in) |

`run_all_analyses` never creates a pool, so the overnight run got 10 workers
auto-started from `parallel.defaultProfile`. The re-run hard-coded
`parpool('Processes', 5)` — the wrong profile, and a worker count chosen from a
memory estimate that was itself built on a mis-read cap (`parcluster('Processes')`
reports 6, which is the built-in profile, not the default).

Halving the workers alone would explain the success under a memory theory; the
fresh pool alone would explain it under a pool-age theory. The observation
discriminates nothing.

**Next step:** run `tau_sensitivity` alone from a fresh MATLAB at the full
default pool, with per-job memory logging, so any failure carries numbers.

---

## ⚪ WONTFIX · ISSUE-008 · The LLE tracks `−1/max(tau_a)` when the network is stable — NOT A DEFECT

| | |
|---|---|
| Identified | 2026-08-13 · `dev` @ `ac59b42` · R5456622 |
| Reclassified | 2026-08-14 · `dev` @ `b84e898` · R5456622 — **expected consequence of the equations; this is a result, not a bug** |
| Area | interpretation |

**Resolution: won't fix, because there is nothing to fix.** This was first
written up as an artefact that "says nothing about the network". That framing
was wrong. `a` is a state variable of the closed-loop system like any other, so
its relaxation rate is a network property, and the largest Lyapunov exponent
correctly reports the slowest mode of the coupled linearisation. When the
recurrent dynamics contract faster than the adaptation, the slowest mode simply
*is* the adaptation — and the LLE saying so is the right answer, not a masked
one.

The measured values support that reading rather than a bare eigenvalue: they sit
slightly *below* −0.1 (−0.1004 … −0.1025), which is the `a`↔`x` coupling
perturbing the uncoupled −1/tau.

Nor is a near-zero exponent from slow adaptation misleading in itself: a network
with `max(tau_a) = 60 s` has LLE ≈ −0.017 because perturbations genuinely
persist for ~60 s. That is a real statement about the system's memory.

**What survives as a caution, and it is narrow:** where the LLE sits flat across
a swept parameter, it means the exponent is *insensitive* to that parameter in
that regime — the recurrent dynamics may be changing underneath while the
slowest mode does not move. If the recurrent contraction rate specifically is
what is wanted, that is a different exponent in the QR spectrum, not the LLE.
`best_presets.md`'s flat plateaus should be read that way.

**Turned into a deliberate measurement.** The tau sweep exists precisely to show
that `tau_a_E` controls the exponent in the stable regime, and its range moved
from [5, 60] to **[1, 30]** to serve that. The prediction to check is not a bare
−1/tau line but `max(−1/tau_a_E(end), STD mode ≈ −0.6)` — a knee near
`tau_a_E(end) ≈ 1.6 s`, which [1, 30] brackets and [5, 60] never reached.

No preset overrides `tau_a`, so all use the class default
`logspace(log10(0.25), log10(10), n_a)` — slowest timescale **10 s**, giving a
Jacobian eigenvalue of exactly **−0.1**.

Whenever the recurrent dynamics contract faster than that, the largest Lyapunov
exponent *is* the slowest linear adaptation mode and says nothing about the
network. Evidence: all 20 top QR exponents in `test_benettin_vs_qr` within 0.002
of −0.1; `sfa_only` / `sfa_and_std` pinned at −0.101 across unrelated sweeps.

The depression equivalent is `−(1/tau_rec + r/tau_rel)` ≈ **−0.62**, matching
`best_presets.md`'s *"`std_only` sits flat near −0.65"*.

**The tau stage is its own test** — it sweeps `tau_a_E(end)` over [5, 60], i.e.
the floor itself. Plot LLE against `−1/tau_a_E_max`. Supporting evidence already
in hand: the minimum successful LLE was **−0.2159** against a predicted floor of
−1/5 = −0.2.

**Consequence:** every flat plateau in these sweeps must be read as
"contracting, floor-limited", not as a measured exponent — including claims in
`best_presets.md`.

---

## ⚪ WONTFIX · ISSUE-007 · `best_presets.md` exponents predate the Lyapunov window fix

| | |
|---|---|
| Identified | 2026-08-13 · `dev` @ `21acf41` · R5456622 |
| Closed | 2026-08-19 · WONTFIX by TR |
| Documented by | `9777966` (caveat added; numbers not re-derived) |

**WONTFIX by TR 2026-08-19: this is not an issue.** `best_presets.md` is a
working note for choosing presets, not a results document, and the caveat added
in `9777966` is enough. Re-deriving the exponents is not worth the compute. Do
not re-open this to "refresh the numbers".

Every LLE in that file was measured before ISSUE-001, when accumulation ran from
`t = 0` rather than over `lya_T_interval`, so all include part of the settling
transient.

Scope: **the recorded exponents only.** Commits `c2d7f22` and `45821dc` also
narrowed the `level_of_chaos` sweeps to `[0.5, 1.5]`, and an earlier version of
this entry treated that range as suspect too. It is not — sweep ranges are TR's
experimental choice, and this file is for code defects and things that make
results inaccurate. Re-deriving the numbers is the whole of this issue.

---

## 🟢 FIXED · ISSUE-006 · `test_srnn_param_preset` checked 2 of 11 presets

| | |
|---|---|
| Identified / fixed | 2026-08-13 · `dev` @ `21acf41` → `a7e43a9` · R5456622 |

`preset_names` was a hand-maintained `{'default','overconnected'}`, so six
presets added on `refactorPSA` plus one added here were never checked against the
banned-key list or `validate_model_defaults`. The mode list also omitted
`fast2`.

**Fixed** by reading the preset list out of `srnn_param_preset`'s own
unknown-preset error, so a preset added without touching the test is still
covered, and threading `model_class` through the validation.

---

## 🟢 FIXED · ISSUE-005 · An unaccumulable Lyapunov window returned `LLE = 0` silently

| | |
|---|---|
| Identified / fixed | 2026-08-13 · `refactorPSA` @ `200a1cf` → `fe767ff` · R5456622 |

If no Lyapunov segment fell inside the window, the estimate came back as a bare
`0`, which reads as edge-of-chaos rather than as not-measured. Now warns.

---

## 🟢 FIXED · ISSUE-004 · `lya_method = 'qr'` could not work with a fixed-step integrator

| | |
|---|---|
| Identified / fixed | 2026-08-13 · `refactorPSA` @ `200a1cf` → `7f1e5cd` · R5456622 |

The QR variational solve was handed the model's own solver and integrates on a
**2-point span**, which `ode_rk4` rejects outright — so `'qr'` only ever worked
with `ode45`. The variational equation is deterministic and independent of how
the fiducial trajectory was produced, so it now pins `ode45` regardless.

---

## 🟢 FIXED · ISSUE-003 · Dead min/max clipping in `SRNNModel2`'s Benettin

| | |
|---|---|
| Identified / fixed | 2026-08-13 · `refactorPSA` @ `200a1cf` → `fe767ff` · R5456622 |

`get_minMaxRange_internal` returned all-NaN, so the clipping never fired. Left
in place it was a trap: had real bounds ever been filled in, clipping would have
shortened the perturbation below `d0` while the algorithm kept dividing by `d0`,
biasing every exponent downward with no visible symptom. Deleted.

---

## 🟢 FIXED · ISSUE-002 · `lya_dt` was hardcoded

| | |
|---|---|
| Identified / fixed | 2026-08-13 · `refactorPSA` @ `200a1cf` → `9c45023` · R5456622 |

Fixed at 0.02 (Benettin) / 0.1 (QR) inside both classes, so the reshooting
interval could not be tuned without editing the class. Now a property; empty
means the per-method default, so nothing moved.

---

## 🟢 FIXED · ISSUE-001 · `SRNNModel2`'s Benettin ignored `lya_T_interval(1)`

| | |
|---|---|
| Identified / fixed | 2026-08-13 · `refactorPSA` @ `200a1cf` → `fe767ff` · R5456622 |

It used the interval only to trim the last sample and gated accumulation on
`t >= 0`, so the documented "skip the first 15 s" had **no effect**.
`SRNNCellTypePairs` honoured it, so the two classes reported different windows
for the same nominal config — which mattered once `run_all_analyses` began
driving both. Both QR paths ignored it too.

Fixed alongside the new `lya_warmup` property. Measured consequence: with
adequate warmup the fixed code reproduces the pre-fix values exactly on a
stationary window (Benettin −0.1013, QR −0.1004), so the fix is
behaviour-preserving where it should be and functional where it was not.

**Downstream:** see ISSUE-007 — everything measured before this needs
re-deriving.
