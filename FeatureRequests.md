# Feature Requests

Wanted changes that are **not** defects. Bugs and observed problems live in
`Issues.md`; the chronological narrative lives in `WorkLog.md`.

Each entry records when it was raised, on which branch and commit, by which
host, and how it ended: **done**, **declined** (a deliberate decision, kept so
the same idea is not re-proposed), or still **open**.

Newest first. Headings carry a plain-text status token beside the marker so the
list greps reliably (see `Issues.md` for why emoji alone is not enough):
🔴 `OPEN` · 🟢 `DONE` · ⚪ `DECLINED`.

```bash
grep '^## ' FeatureRequests.md            # whole list
grep '^## .* OPEN ·' FeatureRequests.md   # just what is outstanding
```

---

## 🟢 DONE · FR-007 · A human-readable `parameters.md` in every run directory

| | |
|---|---|
| Raised | 2026-08-18 · `dev` @ `e5f7f44` · R5611351 |
| Done | 2026-08-18 · `dev` @ `e5f7f44` · R5611351 |

A `run_all_*` directory described itself only in `.mat` form, so learning what a
run used meant starting MATLAB, loading objects, and knowing that
`resolved_defaults` deliberately excludes the grid axes and the condition
fields.

`src/write_run_parameters_md.m` writes `parameters.md` next to
`run_manifest.mat`: preset and run mode at the top, what the preset set, the
per-analysis sweep axes and timings, the four adaptation conditions, and one
unified list of every parameter in effect as run under `sfa_and_std`, each
tagged with what set it. `run_all_analyses.m` calls it after stage 3; it can be
pointed at any past directory by hand.

It reads the run's own saved artifacts, never the current source — see the
`WorkLog` entry of the same date for why that matters and what the current
preset is still used for (a drift diff).

---

## ⚪ DECLINED · FR-006 · Memory pre-flight for `ParamSpaceAnalysis2`

| | |
|---|---|
| Raised | 2026-08-14 · `dev` @ `d58b7fe` · R5456622 |
| Declined | 2026-08-20 · `dev` @ `b7fe0dc` · R5469844 — not wanted |

**Closed by TR, 2026-08-20.** Not needed. The request was written by the agent,
not by TR, and the pre-flight is not something he wants carried.

Estimate the per-worker footprint at config time and warn (or refuse) before a
long sweep starts, alongside `validate_noise_settings`. The arithmetic is
already known:

```
S_out  = nt * N_sys_eqs * 8
noise  = n * (nt-1) * 2 * 8      (zero when sigma_u_noise == 0)
u_ex   = 2 * n * nt * 8          (array + interpolant)
```

times the pool size. `medium2`'s tau stage wants ~1.3 GB/worker at `nt = 40001`;
that is the number a pre-flight would have surfaced before burning a night.

Worth doing **even though ISSUE-009 shows memory was not the cause** — the
estimate is cheap, and having it printed would itself have shortened that
diagnosis considerably.

---

## 🔴 OPEN · FR-005 · Separate the network seed from the noise seed for reps

| | |
|---|---|
| Raised | 2026-08-13 · `dev` @ `ac59b42` · R5456622 |

`noise_seed` derives from `rng_seeds(1)`, so a rep varies the network **and**
the noise realisation together. That is the right default — it is what makes
the noise path shared across the four adaptation conditions at a grid point —
but it means the within-level spread in a reps histogram mixes two sources and
neither can be attributed.

Wanted: a way to hold one fixed while varying the other, so "how much of this
spread is the network and how much is the noise" becomes answerable.

---

## ⚪ DECLINED · FR-004 · Retune `analysis_run_config` windows after the Lyapunov fix

| | |
|---|---|
| Raised | 2026-08-13 · `dev` @ `21acf41` · R5456622 |
| Declined | 2026-08-14 · R5456622 — the window is deliberate |

**Closed by TR, 2026-08-14.** Accumulating over the second half of `T_range` is
the intended design: it discards the initial-condition transient and saves
compute. The post-fix accumulation time is what was always meant to happen; the
pre-fix behaviour was the wrong one. Nothing to retune.

Every run mode sets `lya_T_interval` to exactly the second half of `T_range`.
Once ISSUE-001 was fixed and the window actually took effect, the accumulation
time **halved** in fast/fast2/medium and dropped ~40% at production, making
every LLE roughly √2 noisier for the same compute.

Deliberately left alone at the time so the fix stayed a pure bug fix. Now worth
revisiting with real post-fix output in hand — the natural question is whether
to start the window earlier (say the last three quarters) or lengthen `T_range`.

---

## 🔴 OPEN · FR-003 · Reconsider `medium2`'s sweep sizes

| | |
|---|---|
| Raised | 2026-08-13 · `dev` @ `ac59b42` · R5456622 |

`medium2` is **not** the medium/production midpoint, despite the name implying
it. Calibrated against the `fast2` run, medium alone is ~9 h and production
~80 h, so the true midpoint is a multi-day job. It was sized to fit a night and
therefore sits nearer medium (13 levels / 12 reps against medium's 11 / 15).

If overnight is no longer the binding constraint — a weekend, or the other
machine — the levels and reps are the knobs to raise.

---

## 🟢 DONE · FR-002 · Record the machine in run provenance

| | |
|---|---|
| Raised / done | 2026-08-14 · `dev` @ `3794a09` · R5456622 |

`capture_git_provenance` recorded commit and branch but not **where** the run
happened. "Which commit" alone does not distinguish two runs of the same code on
different hardware, and the machines differ in exactly the ways that mattered to
ISSUE-009 — 10 workers here against 14 on the other.

Done: `git_provenance.txt` now carries `hostname`, `platform` and the MATLAB
release, and `info` returns hostname/platform for machine-readable use.
`hostname` was chosen over a UUID file or MAC address because it needs no setup,
is stable, is human-meaningful, and reads identically on macOS, Windows and
Linux.

---

## 🟢 DONE · FR-001 · Fresh parallel pool before each analysis

| | |
|---|---|
| Raised / done | 2026-08-14 · `dev` @ `d58b7fe` → `c6055ea` · R5456622 |

`run_all_analyses` ran three multi-hour analyses back to back against one
long-lived pool. `src/restart_parpool.m` now restarts it before each stage,
using `parallel.defaultProfile` rather than a hard-coded worker count so the two
machines each get their own.

Explicitly **insurance, not a fix** for ISSUE-009 — it bounds whatever the
workers accumulate, at ~30 s per stage, without any claim about what actually
went wrong.
