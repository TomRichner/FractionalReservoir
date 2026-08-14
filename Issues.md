# Issues

Bugs and observed problems, tracked **by issue** rather than chronologically —
`WorkLog.md` is the chronological record and this is the index into it.

Every issue records when it was identified, on which branch and commit, by which
host, and how it ended: **fixed**, **won't fix** (a deliberate decision), or
still **open**. Issues are never deleted; a resolved one keeps its entry with the
resolving commit, so the same ground is not re-covered in six months.

Newest first. Status: 🔴 open · 🟢 fixed · ⚪ won't fix · 🔵 observation (not a
defect, but affects interpretation).

---

## 🔴 ISSUE-011 · `load_results` silently returns level indices for vector parameters

| | |
|---|---|
| Identified | 2026-08-14 · `dev` @ `d58b7fe` · R5456622 |
| Area | `src/ParamSpaceAnalysis2.m` |

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

**Why the tests missed it:** `test_psa_saveload` exercises `saveobj`/`loadobj`
(the path that works) and never calls `load_results` (the path the replot
scripts and post-hoc analysis actually use).

---

## 🔴 ISSUE-010 · `load_results` does not restore `model_class`

| | |
|---|---|
| Identified | 2026-08-14 · `dev` @ `d58b7fe` · R5456622 |
| Area | `src/ParamSpaceAnalysis2.m:1612-1621` |

`summary_data.model_class` **is** saved (`:2495`) but `load_results` restores
only six fields and this is not one of them. A loaded `SRNNCellTypePairs` run
reports `SRNNModel2`.

Consequences: `effective_param`'s class-default fallback consults the wrong
class, and `same_config`'s documented refusal to pool runs from different model
classes cannot work on loaded runs.

Same root cause as ISSUE-011 (an incomplete `load_results`) and same test blind
spot; worth fixing together.

---

## 🔴 ISSUE-009 · `medium2` tau/param_space stages died — cause NOT established

| | |
|---|---|
| Identified | 2026-08-14 · `dev` @ `ac59b42` · R5456622 |
| Area | run infrastructure / MATLAB stability |
| Mitigated by | `c6055ea` (insurance, not a fix) |

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

## 🔵 ISSUE-008 · The LLE has a floor at `−1/max(tau_a)`

| | |
|---|---|
| Identified | 2026-08-13 · `dev` @ `ac59b42` · R5456622 |
| Area | interpretation of every sweep, not a code defect |

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

## 🔴 ISSUE-007 · `best_presets.md` exponents predate the Lyapunov window fix

| | |
|---|---|
| Identified | 2026-08-13 · `dev` @ `21acf41` · R5456622 |
| Documented by | `9777966` (caveat added; numbers not yet re-derived) |

Every LLE in that file was measured before ISSUE-001, when accumulation ran from
`t = 0` rather than over `lya_T_interval`, so all include part of the settling
transient.

Not just bookkeeping: commits `c2d7f22` and `45821dc` narrowed both
`level_of_chaos` sweeps to `[0.5, 1.5]` **on the strength of those numbers**. A
range chosen from shifted values may no longer bracket the crossing.

---

## 🟢 ISSUE-006 · `test_srnn_param_preset` checked 2 of 11 presets

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

## 🟢 ISSUE-005 · An unaccumulable Lyapunov window returned `LLE = 0` silently

| | |
|---|---|
| Identified / fixed | 2026-08-13 · `refactorPSA` @ `200a1cf` → `fe767ff` · R5456622 |

If no Lyapunov segment fell inside the window, the estimate came back as a bare
`0`, which reads as edge-of-chaos rather than as not-measured. Now warns.

---

## 🟢 ISSUE-004 · `lya_method = 'qr'` could not work with a fixed-step integrator

| | |
|---|---|
| Identified / fixed | 2026-08-13 · `refactorPSA` @ `200a1cf` → `7f1e5cd` · R5456622 |

The QR variational solve was handed the model's own solver and integrates on a
**2-point span**, which `ode_rk4` rejects outright — so `'qr'` only ever worked
with `ode45`. The variational equation is deterministic and independent of how
the fiducial trajectory was produced, so it now pins `ode45` regardless.

---

## 🟢 ISSUE-003 · Dead min/max clipping in `SRNNModel2`'s Benettin

| | |
|---|---|
| Identified / fixed | 2026-08-13 · `refactorPSA` @ `200a1cf` → `fe767ff` · R5456622 |

`get_minMaxRange_internal` returned all-NaN, so the clipping never fired. Left
in place it was a trap: had real bounds ever been filled in, clipping would have
shortened the perturbation below `d0` while the algorithm kept dividing by `d0`,
biasing every exponent downward with no visible symptom. Deleted.

---

## 🟢 ISSUE-002 · `lya_dt` was hardcoded

| | |
|---|---|
| Identified / fixed | 2026-08-13 · `refactorPSA` @ `200a1cf` → `9c45023` · R5456622 |

Fixed at 0.02 (Benettin) / 0.1 (QR) inside both classes, so the reshooting
interval could not be tuned without editing the class. Now a property; empty
means the per-method default, so nothing moved.

---

## 🟢 ISSUE-001 · `SRNNModel2`'s Benettin ignored `lya_T_interval(1)`

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
