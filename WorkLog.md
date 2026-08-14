# Work Log

Chronological record of work on this repo: plans executed, bugs found, bugs
fixed, and analysis runs performed. Newest entries at the **bottom**.

This exists because the work happens on more than one machine and across many
sessions, and neither git history nor the `data/` run directories capture the
*reasoning* — what was tried, what turned out to be wrong, and what is still
open.

Companion documents, indexed by topic rather than by date:

- **`Issues.md`** — bugs and observed problems, with status and resolution
- **`FeatureRequests.md`** — wanted changes that are not defects

This file is the narrative; those two are the indexes into it. `CLAUDE.md`
tells the AI to keep all three current.

## Conventions

Every entry opens with a signature block identifying **who**, **where**, and
**against what commit**:

```
### 2026-08-14 · dev @ d58b7fe · R5456622 · Claude Code (Opus 5), session 054a29ca
```

- **date** — ISO, local time
- **branch @ commit** — as of the START of the entry (`git rev-parse --short HEAD`)
- **machine** — the hostname (see below)
- **agent** — the tool and model, plus the session id, so an entry can be traced
  back to a transcript

### Identifying the machine

`hostname` is the identifier. It needs no setup, is stable, and works
identically on macOS, Windows and Linux:

```bash
hostname                       # -> R5456622
```

```matlab
h = getenv('COMPUTERNAME');                       % Windows
if isempty(h), [~, h] = system('hostname'); end    % macOS / Linux
h = strtrim(h);
```

Known machines:

The asset-tag names are the identifier as-is — they are what the institution
uses, so they are unambiguous across machines in a way a nickname would not be.

| hostname | OS | notes |
|---|---|---|
| `R5456622` | Windows 11, MATLAB R2025b | 6 physical cores, 64 GB RAM, default pool profile `ten_workers_one_thread` (10 workers) |
| _(other)_ | — | the second machine; add its hostname on first use. 14 workers. |

`src/capture_git_provenance.m` writes the same hostname (plus platform and
MATLAB release) into every run's `git_provenance.txt`, so a log entry here and a
directory under `data/` can be matched to the same machine.

### Merging this file across machines

A single chronological history is worth more than conflict-free separate files,
so this is one file **by choice**. Two machines appending will collide at the
same place in git. The rule:

> **Keep both sides and order by timestamp. Never choose between them.**

Entries are self-contained and `---`-separated precisely so that is always safe
and needs no judgement.

---

### 2026-08-13 · SDE branch created from refactorPSA @ 200a1cf · R5456622 · Claude Code (Opus 5), session 054a29ca

**Plan:** add an additive Wiener process to `dx/dt` in both model classes, with
a choice of SDE integrators, keeping Benettin's LLE estimate valid.

Design decisions taken before writing code:

- Noise enters **only `x`**, as an input-referred current. A state-independent
  diffusion is what makes Itô = Stratonovich, kills the Milstein term, leaves
  the QR variational equation untouched, and — because the fiducial and
  perturbed trajectories share one noise path — makes the additive noise
  **cancel** in Benettin's difference. So the LLE stays valid at any noise
  level, not just small ones. Noise on `b` or `a` would break all four.
- `sigma_u_noise` is input-referred (units of `u`); `sigma_x_raw` and
  `x_noise_std` are Dependent views.
- The noise cannot live in `dynamics_fast` — for additive noise it never enters
  `f` at all, only the integrator's update formula.

**Bugs found while investigating (all pre-existing):**

1. **`SRNNModel2`'s Benettin ignored `lya_T_interval(1)`.** It used the interval
   only to trim the last sample and gated accumulation on `t >= 0`, so the
   documented "skip the first 15 s" had no effect. `SRNNCellTypePairs` honoured
   it, so the two classes reported **different windows for the same config**.
   Both QR paths ignored it too.
2. **`lya_dt` was hardcoded** in both classes, so the reshooting interval could
   not be tuned without editing the class.
3. **Dead min/max clipping** in `SRNNModel2`'s Benettin.
   `get_minMaxRange_internal` returned all-NaN, so it never fired — but had real
   bounds ever been filled in, clipping would have shortened the perturbation
   below `d0` while the algorithm kept dividing by `d0`, biasing every exponent
   downward with no visible symptom.
4. **`lya_method = 'qr'` could not work with a fixed-step integrator.** The QR
   variational solve was handed the model's own solver and integrates on a
   2-point span, which `ode_rk4` rejects outright.
5. **An unaccumulable Lyapunov window returned `LLE = 0` silently**, which reads
   as edge-of-chaos rather than as not-measured.

**Fixed** — `fe767ff`, `9c45023`, `7f1e5cd`:

- Both methods now accumulate over `[max(0, lya_T_interval(1)), lya_T_interval(2)]`
  via one shared `lyapunov_sample_grid` helper per class.
- New `lya_warmup` (default 5 s): iteration begins that much earlier without
  accumulating, so the perturbation/basis can align. Clamped to `T_range(1)`
  with a warning rather than dropped.
- `lya_dt` promoted to a property, empty meaning the per-method default.
- `ode_solver` became a **name** (`'ode45' | 'ode15s' | 'rk4' | 'euler' | 'heun' | 'sra1'`);
  handles now raise `RenamedProperty`. `resolved_defaults` consequently
  contains **no function handles at all**.
- QR's variational solve pins `ode45` regardless of `ode_solver`, so
  `'qr'` + `'rk4'` works for the first time.

**Measured, and it set the `lya_warmup` default.** Too little warmup biases the
exponent badly, and QR needs more of it than Benettin:

```
warmup:   0.5s     1s      2s      5s     10s     25s
Benettin -0.1318 -0.1177 -0.1062 -0.1015 -0.1013 -0.1013
QR            -       -  -0.2197 -0.1121 -0.1004 -0.1003
```

Convergence tracks **physical time**, not renormalisation count. Default 5 s is
comfortable for Benettin (0.2% out) but **not** converged for QR (~12% out).
`test_benettin_vs_qr` therefore asks for 10 s explicitly. At that setting it
reports Benettin −0.1013 / QR −0.1004 — identical to the pre-fix values, so the
window fix is behaviour-preserving on a stationary window.

**SDE core** — `7d88bdb`, `2c90372`:

- `src/sde_fixed_step.m` implements Euler–Maruyama, stochastic Heun and Rößler
  SRA1 in one file, so the absolute-time noise indexing is written once.
- **Prefer `'sra1'`**: same two drift evaluations per step as `'heun'` but
  strong order 1.5 rather than 1.0 — measured ~85× more accurate at the same
  step (3.9e-5 vs 3.3e-3).
- The SRA1 tableau was **not** trusted from transcription. The test measures the
  convergence order, and the threshold was validated against deliberately
  broken variants: correct → 1.83, `dZ` term dropped → **1.03**, `B0(2,1)`
  mistyped as 1 → **0.90**. A typo degrades silently to order 1.0, which is
  exactly what the ≥1.4 floor catches.
- Verified end to end that `sigma_u_noise` means what the docs say: with `W = 0`
  the dendritic state is an OU process whose measured stationary variance lands
  **0.2%** (SRNNModel2) / **0.5%** (Pairs) from analytic.
- Mechanism confirmed: noise lowers the mean effective gain
  `⟨φ'⟩` 0.3227 → 0.2169 at σ = 0.6. That is how additive noise lowers the
  Lyapunov exponent.

**Example** — `0751b0a`: `example_SRNNCellTypePairs_noise` sweeps σ and shows
the network crossing from chaotic to contracting at σ ≈ 0.03, with `std(x)`
*rising* as the LLE falls (more stable ≠ quieter) and the LLE scatter
collapsing from ±0.28 to ±0.003.

---

### 2026-08-13 23:01 · dev @ 0751b0a · R5456622 · Claude Code (Opus 5), session 054a29ca

**Plan:** merge `origin/refactorPSA` (12 commits: preset chain, four mu-block
sweeps, narrowed grid ranges, `best_presets.md`) into `dev` alongside the SDE
work, then run `run_all_analyses` at `fast2` on a noisy preset.

**Merge** — `21acf41`. Conflict-free, verified in advance with `git merge-tree`.
Only two files were touched by both sides, in different hunks. Regression gate
after merging: `test_benettin_vs_qr` still −0.1013 / −0.1004.

**Structural problem found.** The repo keeps two knobs deliberately orthogonal —
a preset says *which network* (physics), `analysis_run_config` says *how much
compute*. But `ode_solver` lives in the compute table while `sigma_u_noise` is
physics, and σ > 0 **requires** a stochastic integrator. Since the sub-scripts
merge with `cfg.model` winning, a noise preset at `fast2` produced
`sigma_u_noise = 0.02` alongside `ode_solver = 'rk4'` and died in PSA's
pre-flight with no way to say otherwise. Not a one-off: all **12** table cells
hardcoded a deterministic solver.

**Fixed** — `a7e43a9`: every cell now names two integrators (deterministic
`rk4` up to medium / `ode45` at production; stochastic `sra1` everywhere) and
the **preset** picks between them. `merge_struct` precedence and the
preset ban on `ode_solver` are both untouched; a σ = 0 preset is bit-identical
to before. `tau_sensitivity/medium` changed `ode45` → `rk4` as part of
normalising the deterministic column — a deliberate behaviour change.

**Test gap closed:** `test_srnn_param_preset` checked 2 of 11 presets. It now
reads the list from `srnn_param_preset`'s own error message, so a preset added
without touching the test is still covered.

**Documented** — `9777966`: every LLE in `best_presets.md` predates the window
fix and will shift. This matters because commits `c2d7f22` / `45821dc` narrowed
both `level_of_chaos` sweeps to `[0.5, 1.5]` *on the strength of those numbers*.

**RUN — `fast2`, preset `..._sig1p5_noise0p02`, σ = 0.02, SRA1, fs = 200.**
✅ Completed in **19.68 min**, every batch 100% successful.
Output: `data/param_space/run_all_aug_13_26_23_13/`.

---

### 2026-08-13 23:37 · dev @ ac59b42 · R5456622 · Claude Code (Opus 5), session 054a29ca

**Plan:** a `medium2` run mode between medium and production at fs = 800 with
SRA1, plus a σ = 0.01 preset, run overnight.

Sizing note: fs = 800 is **free relative to medium per simulated second** —
SRA1 takes two drift evaluations per step where rk4 takes four, so `sra1@800`
and `rk4@400` both cost 1600 evaluations per second of simulation.

Calibrating against the completed `fast2` run, medium alone is already ~9 h and
production ~80 h, so a true medium/production midpoint is a **multi-day** job.
`medium2` was therefore sized to ~9 h and sits nearer medium. This is recorded
in the file so it is not later mistaken for the midpoint.

**RUN — `medium2`, preset `..._sig1p5_noise0p01`, σ = 0.01, SRA1, fs = 800.**
⚠️ **Partially failed.** Output: `data/param_space/run_all_aug_13_26_23_37/`.

| stage | result |
|---|---|
| 7 × sensitivity sweeps | ✅ **624/624 clean each** (4 368 runs, zero failures) |
| replot / combined sheets | ✅ 34 figures |
| tau_sensitivity | ❌ **88 of 195 runs: `Out of memory`** |
| param_space | ❌ crashed MATLAB after 1 of 5 batches |

Config verified on all eight completed stages: `sra1`, σ = 0.01, fs = 800 — the
mechanism worked; the failure is resource/stability, not correctness.

The two failed stage directories were renamed `FAILED_OOM_*` rather than
deleted, so nothing masquerades as good data.

---

### 2026-08-14 07:06 · dev @ d58b7fe · R5456622 · Claude Code (Opus 5), session 054a29ca

**Diagnosis of the overnight failure — TWO wrong answers before the right
question.** Recorded in full because the wrong turns are the useful part:

1. First claim: total RAM exhausted with 14 workers. **Wrong** — the machine has
   6 physical cores and the pool was 10, not 14.
2. Second claim: heap fragmentation. **Weak** — that is a 32-bit-era failure
   mode; on 64-bit the virtual address space is ~128 TB and contiguity
   essentially never binds.

What the data actually shows:

| | |
|---|---|
| commit limit | 128.7 GB (64 GB RAM + 65 GB page file) |
| page file **peak** | **7.8 GB** over 10 days of uptime, including the run |
| `MaxPossibleArrayBytes` | 99.1 GB — no array-size-limit preference set |
| estimated peak | ~1.4 GB/worker × 10 = ~14 GB |

**Nothing was close to running out.** Cause **NOT ESTABLISHED.** The leading
hypothesis is now that a *worker process crashed* — there are two crash dumps,
`05:41` (exactly when the tau stage began; the run continued afterwards, so this
was a worker not the client) and `06:21` (the client). If so, `Out of memory` is
a **symptom of the worker dying, not the reason**, and no memory setting will
help.

Also recorded: the apparent "fix" (a re-run whose first batch went 25/25)
changed **three things at once** — pool age (6 h → fresh), workers (**10 → 5**)
and cluster profile (`ten_workers_one_thread` → the built-in `Processes`). The
overnight run never created a pool explicitly, so `parfor` auto-started 10
workers from the default profile; the re-run hard-coded `parpool('Processes', 5)`
using a count derived from a mis-read cap. It validates nothing, and the
worker change was a halving rather than the trim first recorded here.

**Fixed** — `c6055ea`: `src/restart_parpool.m`, called before each of the three
analyses in `run_all_analyses`. Uses `parallel.defaultProfile` rather than a
hard-coded worker count, since the two machines differ (10 here, 14 there).
Framed as **insurance, not a diagnosis**: it bounds whatever workers accumulate
over a multi-hour run, at ~30 s per stage.

**Also fixed** — `d58b7fe`: stopped tracking `debug.log`, swept in by a blanket
`git add -A`.

---

### 2026-08-14 07:25 · dev @ 3794a09 · R5456622 · Claude Code (Opus 5), session 054a29ca

**Plan:** set up cross-session, cross-machine bookkeeping.

- `WorkLog.md` (this file) started, seeded with everything above.
- `Issues.md` and `FeatureRequests.md` created, indexing the same material by
  topic. Eleven issues and six feature requests back-filled from this session.
- `CLAUDE.md` now instructs the AI to read all three at session start and update
  them as work proceeds, with the signature format and the merge rule.
- `capture_git_provenance` extended (FR-002) to stamp `hostname`, `platform` and
  MATLAB release into every run directory, so a log entry and a `data/` folder
  can be matched to a machine. Verified: `dev @ 3794a09 on R5456622`, `PCWIN64`,
  `2025b`.

Decision recorded: a **single** chronological `WorkLog.md` rather than
`docs/worklog/<date>-<host>.md`. Per-file-per-session would never conflict, but a
single history is worth more; the cost is a trivial merge, resolved by keeping
both sides in timestamp order.

## Open items

Tracked in full in **`Issues.md`** and **`FeatureRequests.md`**; summarised here:

- 🟢 **ISSUE-011 / ISSUE-010** — **fixed 2026-08-14** (`b6e5262`, test
  `b445e9b`). `psa_object.mat` is now the one authoritative artifact and
  `from_dir` the only reader; `effective_param` errors rather than returning a
  level index.
- 🔴 **ISSUE-009** — the `medium2` tau/param_space failure. Cause **not
  established**; two diagnoses refuted.
- 🔵 **ISSUE-008** — the `−1/max(tau_a)` LLE floor. Affects how every flat
  plateau in these sweeps must be read.
- 🔴 **ISSUE-007** — `best_presets.md` numbers predate the window fix and need
  re-deriving, along with the `[0.5, 1.5]` range narrowing that rests on them.
- 🔴 **FR-006** — memory pre-flight for PSA.

## Unfinished work — and the order it has to happen in

**`medium2` `tau_sensitivity` + `param_space` re-run — two of the three blockers
are now clear:**

1. ~~**ISSUE-008 (the LLE floor)**~~ — **settled 2026-08-14: not a defect.** See
   below; the sweep range changed as a result, so the re-run should use the new
   one.
2. ~~**ISSUE-010 / ISSUE-011**~~ — **fixed 2026-08-14** (`b6e5262`). A crashed
   run is exactly the case the loader got wrong, and the last run *did* crash;
   that data now decodes correctly and a future crash stays salvageable, since
   `run()` writes the object before batching.
3. **ISSUE-009's cause is still unknown.** `restart_parpool` is insurance, not a
   fix. Relaunching without instrumentation risks burning another night and
   learning nothing again — see FR-006 for the memory pre-flight that would make
   any repeat failure legible.

Remaining order: fold ISSUE-009's diagnostic (tau alone, fresh MATLAB, full
default pool, per-job memory logging) into the `medium2` re-run.

---

### 2026-08-14 · dev @ b84e898 · R5456622 · Claude Code (Opus 5), session 054a29ca

**ISSUE-008 settled — it is a result, not a bug.** I had written up the LLE
sitting at `−1/max(tau_a)` as an artefact that "says nothing about the network".
Corrected by TR: `a` is a state variable of the closed-loop system, so its
relaxation rate *is* a network property, and the largest Lyapunov exponent
correctly reports the slowest mode of the coupled linearisation. When the
recurrent dynamics contract faster than the adaptation, the slowest mode simply
is the adaptation. Supporting detail: the measured values sit slightly *below*
−0.1 (−0.1004…−0.1025), which is the `a`↔`x` coupling perturbing the uncoupled
−1/tau — i.e. it is the coupled system's eigenvalue, not a bare parameter
readback.

What survives is narrow: a flat LLE across a swept parameter means the exponent
is *insensitive* to it in that regime, and if the recurrent contraction rate is
what is wanted specifically, that is a different exponent in the QR spectrum.

**Consequence — the tau sweep became a deliberate measurement of this.** Its
whole point is to show `tau_a_E` controlling the exponent in the stable regime.
Range changed **[5, 60] → [1, 30] s** on `tau_a_E(end)`: a factor of 30 in
−1/tau rather than 12, and crucially the fast end now reaches where STD's own
slowest mode (≈ −0.6) overtakes adaptation. The prediction to check is therefore
`max(−1/tau_a_E(end), STD mode)` — a **knee near `tau_a_E(end) ≈ 1.6 s`** that
[5, 60] never reached.

Also outstanding:

- ISSUE-007: `best_presets.md` numbers need re-deriving post-window-fix, and the
  `[0.5, 1.5]` range narrowing that rests on them re-confirming.

---

### 2026-08-14 · dev @ c72d1e8 · R5456622 · Claude Code (Opus 5), session 054a29ca

**ISSUE-010 / ISSUE-011 fixed — `psa_object.mat` is now the one authoritative
run artifact.** Commits `b6e5262` (fix) and `b445e9b` (test).

The reported symptom was that `effective_param(res,'tau_a_E')` returned `1`
where the answer was `5 s`. `tau_a_E` is a vector parameter, so it cannot be a
grid coordinate directly: `add_vector_parameter` pre-builds the concrete vectors
into `vector_param_lookup` and the grid axis carries a **level index**. The runs
were always correct — the index is stored by design. But
`param_space_summary.mat` never stored the lookup, so after `load_results` there
was nothing to decode with and `effective_param` fell through to its scalar
branch. No error. The same list drift dropped `model_class`, which is ISSUE-010.

**The fix removes the second load path rather than repairing it.** Three restore
paths each carried a hand-maintained field list, against a `saveobj` that has
long been complete. So:

- `run()` saves `psa_object.mat` itself, **twice** — once after `generate_grid`
  and before batching, once on completion. The early write is the valuable half:
  a run that dies part-way now keeps the configuration needed to interpret its
  `temp_batches/`. The four scripts that used to save it no longer do, which
  also retires the non-canonical `psa_tau_a` / `psa_tau_b` variable names.
- `ParamSpaceAnalysis2.from_dir` is the only way in. It selects the saved
  variable **by class**, not by name, so legacy `psa_tau_a` files load for free.
  Seven readers dropped their hand-rolled three-way choice; three lost a latent
  hardcoded `S.psa` bug.
- `effective_param` **fails closed** on a vector parameter it cannot decode.
  This is the highest-value single change: it converts a silent wrong number
  into a loud one, for the whole class of bug rather than this instance.

**Backward compatibility was deliberately dropped** (TR: "I'm not too worried
about being backward compatible with old runs... I have the git provenance").
A run with a summary but no object file no longer loads. Exactly one such
directory exists on disk — `FAILED_OOM_param_space_*` — and it was already the
broken case.

**Verified against the data that exposed it**, not only against fixtures:

| directory | before | after |
|---|---|---|
| `FAILED_OOM_tau_sensitivity_*` | `tau_a_E → 1` | `[0.25 1.5478 9.5833]`, 107/195 results recovered |
| the 7 completed sensitivity sweeps | `model_class → SRNNModel2` | `SRNNCellTypePairs`; all 7 replot, 14 figures |
| `FAILED_OOM_param_space_*` | loaded wrongly | errors, naming `consolidate()` |

Full suite green: `test_psa_loaders` (new, 19 checks), `test_psa_saveload`,
`test_psa_validate_defaults`, `test_psa_model_class`,
`test_sensitivity_refactor`.

**Two things I could not honestly test, and did not fake.** The crash-recovery
state is not constructible through the public API — `consolidate` deletes
`temp_batches/` on success, and `results`/`has_run` are `SetAccess=private` — so
`test_psa_loaders` asserts instead that the saved object carries the decode
table and the full grid, and the recovery itself rests on the real
`FAILED_OOM_tau_*` directory above. Separately, **file timestamps cannot prove
the early write happened**, because the final write overwrites it; that one
rests on control flow (`run()`:625 precedes `run_batched_simulation`:630).

**Incidental find:** `test_sensitivity_refactor` printed the constant
`vary_range` inside a loop over levels while labelling it the per-level lookup —
so it showed `[5 60]` three times and would have looked like a bug in the
sweep. It also printed before `run()` had built the lookup at all. Moved after
the run, pointed at `vector_param_lookup`, and given an assertion that the
levels actually differ.

**The tau floor analysis this unblocked is inconclusive on the existing data.**
Only 4 of 13 tau levels have surviving runs, because 88/195 OOM'd. That needs
the clean `medium2` re-run, which now waits only on ISSUE-009.
