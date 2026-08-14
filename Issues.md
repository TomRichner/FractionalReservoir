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

## 🔴 OPEN (LOW) · ISSUE-009 · `medium2` tau/param_space stages died — cause NOT established

| | |
|---|---|
| Identified | 2026-08-14 · `dev` @ `ac59b42` · R5456622 |
| Area | run infrastructure / MATLAB stability |
| Mitigated by | `c6055ea` (insurance, not a fix) |
| Priority | **LOW** — deliberately parked 2026-08-14 · R5456622 |

**Parked by TR, 2026-08-14: do not spend more time chasing this.** This was the
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

## 🔴 OPEN · ISSUE-007 · `best_presets.md` exponents predate the Lyapunov window fix

| | |
|---|---|
| Identified | 2026-08-13 · `dev` @ `21acf41` · R5456622 |
| Documented by | `9777966` (caveat added; numbers not yet re-derived) |

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
