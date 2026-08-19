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
- 🔴 **ISSUE-009 (LOW)** — the `medium2` crash. Cause **not established** and
  **deliberately parked**: first crash ever in a big run, no reproduction.
  Reopen only if it recurs.
- 🔵 **ISSUE-008** — the `−1/max(tau_a)` LLE floor. Affects how every flat
  plateau in these sweeps must be read.
- 🔴 **ISSUE-007** — `best_presets.md` exponents were measured pre-window-fix
  and so include settling transient. An accuracy issue about those recorded
  numbers; the sweep *ranges* are TR's experimental choice and not part of it.
- 🔴 **FR-006** — memory pre-flight for PSA.
- ⚪ **FR-004** — closed by TR: the second-half accumulation window is
  deliberate (discards the IC transient, saves compute). Not open.

## Unfinished work

**The `medium2` σ = 0.01 re-run is cancelled — TR is moving on** (2026-08-14).
Both of its substantive blockers had cleared (ISSUE-008 settled as a result, not
a defect; ISSUE-010/011 fixed), and only the unexplained crash remained, which
is now parked at low priority. The decision is not to re-run it: the salvaged
107/195 tau results and the seven clean sensitivity sweeps stand as the record
of that night, and the next big run will be a **`medium`** one instead.

**Nothing blocks the `medium` run.** No open code defects. Two things I raised
as blockers were neither, and both were my error rather than a finding:

- the `level_of_chaos` range — TR's deliberate experimental choice, not a defect
  (see the correction in the 2026-08-14 entry below)
- the Lyapunov accumulation window — the second half of `T_range` is intended:
  it discards the IC transient and saves compute

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

---

### 2026-08-14 · dev @ e3e3102 · R5456622 · Claude Code (Opus 5), session 054a29ca

**Two decisions by TR, and one measurement that came out of asking what still
blocks a `medium` run.**

**1. ISSUE-009 parked at LOW.** No further effort on the overnight crash. It was
the first crash ever in a big run, there is no reproduction, and each diagnostic
attempt costs a night. Reopen at normal priority only if it recurs — two
failures are worth more than more speculation about one.

**2. The `medium2` σ = 0.01 re-run is cancelled.** Moving on. (TR said "noise
0.001"; the mode was actually built at `sigma_u_noise = 0.01` — recording the
value that ran, in case the 0.001 was the intent rather than a slip.)

**3. Where each condition's LLE crosses zero on `[0.5, 1.5]`** — post-fix data,
last night's 13-level sweep, `celltype_pairs` at σ = 0.01, median LLE per level:

| condition | LLE at 0.500 | zero crossing |
|---|---|---|
| `sfa_only` | +0.3975 | below 0.5 (not in the swept range) |
| `no_adaptation` | −0.5958 | 0.500 → 0.583 |
| `sfa_and_std` | −0.1117 | 1.000 → 1.083 |
| `std_only` | −0.6163 | 1.333 → 1.417 |

One preset only; a different preset moves the crossings.

**This is a description of the results, not a defect — corrected after I filed
it as one.** I wrote this up as the range "failing" and as ISSUE-007 being
"worse than" recorded, on the unstated assumption that the sweep ought to
bracket every condition's crossing. TR: *"You seem to be deciding what you think
the results should be and then creating issues if it doesn't do as you assume I
want... this is an experiment, and I have chosen the ranges for a reason."*
Correct. `[0.5, 1.5]` is a deliberate choice, and a condition being chaotic
across the whole range is a finding about that condition.

**The standard for `Issues.md`, restated because I just violated it:** code
defects, bugs, and things that make results *inaccurate* — NaNs, wrong numbers,
silent fallbacks. Not experimental design, not parameter ranges, not results
that differ from what an agent guessed the user wanted. ISSUE-007's actual
content — that `best_presets.md` exponents were measured pre-window-fix and so
include settling transient — is a genuine accuracy issue and stands unchanged.
Nothing about the sweep range belongs in it.

**Note on method.** This came from reading a run already on disk, not from a new
sweep — and it was only readable *because* the loader fix landed first: the same
directory is one of the seven that previously reported `model_class =
SRNNModel2`. Worth recording as a case where the bookkeeping fix immediately
paid for itself.

**Also noted:** `grep`/`sed` patterns containing the status emoji silently
matched nothing again this session, exactly as `CLAUDE.md` warns. The
plain-text tokens (`OPEN`, `FIXED`, …) worked. The warning is doing its job;
the habit is the part that needs repeating.

### 2026-08-18 · dev @ d8886e2 · R5611351 · Claude Code (Opus 5), session 4ac7e853

**Flattened `src/srnn_param_preset.m` — removed the preset recursion chain.**

Ten of the fourteen presets were built by the function calling *itself*, up to
seven hops deep (`celltype_pairs_Sc0p2_noise0p025` → `..._sig1p5` → `..._nodrive`
→ `..._mu5p5` → `..._uniform_std_n500` → `..._n500_fixedF` → `..._n500` →
`celltype_pairs_S_c_by_type`). Every case is now self-contained: its own
`model_class`, its own complete override struct, its own `std_routes`. Lineage
survives as a `Copied from X. Changed: ...` comment block plus a "Derived
presets:" back-reference, so the history of which preset modifies which is still
readable — it just isn't executable any more.

Two latent hazards removed along the way:

- `model_class` defaulted to `'SRNNModel2'` at the top of the function and was
  overwritten only inside some cases. A new `SRNNCellTypePairs` preset that
  forgot the line would have inherited the wrong class silently and failed much
  later inside `validate_model_defaults`. There is no default now.
- Whether a derived preset inherited its `conditions` depended on whether its
  recursive call requested **two** outputs or **three**, and whether it
  `return`ed early. `all_std_n500` and `uniform_std_n500` took two and rebuilt
  conditions; the other eight took three and returned. Load-bearing, and
  invisible at a glance.

**Method (the part worth reusing).** Edited the file *in place* rather than
building a `_rf.m` beside it and renaming later — the function keeps its name
throughout, all ten callers stay untouched, and the whole refactor is one
reviewable diff instead of an add-then-swap pair. The reference is a frozen copy
of the chained implementation at `scripts/tests/srnn_param_preset_old.m`, with
its ten recursive calls repointed at itself.

**That repointing is the whole ballgame**, and it is the thing to get wrong. Had
the frozen copy kept calling `srnn_param_preset`, it would have delegated to the
very implementation under test and the differential test would have passed
unconditionally. `test_srnn_param_preset_equivalence.m` therefore greps the
frozen file for `srnn_param_preset(` and fails if it finds one.

**A wrong turn worth recording.** That guard failed on its first run — and the
only match was the *prose in the frozen file's own header* explaining the guard.
A comment is not a delegation. Fixed by stripping full-line comments before
searching, and then, because a guard that quietly stops matching looks exactly
like a guard that passes, added a positive control asserting the pattern still
fires on a synthetic `srnn_param_preset('default')` call while sparing
`srnn_param_preset_old(` and `srnn_param_preset_names(`. A vacuity guard needs
its own vacuity guard.

**Verified:** all 14 presets agree on `model_defaults` (`isequaln`), field
*order*, `model_class`, and the full `conditions` cell — plus the 1-, 2- and
3-output call forms and the `UnknownPreset` identifier and message.
`test_srnn_param_preset.m` still passes in full (every preset through
`validate_model_defaults`, the banned-name rules, the sra1 solver selection).

`scripts/tests/srnn_param_preset_old.m` is kept rather than deleted, which makes
the equivalence test a standing regression guard on all 14 presets, at the cost
of a ~440-line frozen file in `scripts/tests/`.

**Scoped to the frozen K, so new presets are free** (same session, follow-up).
As first written the test compared *every* preset both functions knew about,
which would have made adding a 15th preset an instant failure — a tax on all
future work, and the fastest way to get a test disabled rather than fixed. It
now compares only the K presets the frozen copy knows about (K read from that
file, not hardcoded), requires them to remain the **first K in order**, and
*names* anything newer as out-of-scope in the output so a green banner is never
misread as "all presets verified". The rule that buys this: **append new presets
at the bottom** of the switch and of `srnn_param_preset_names`. Inserting or
reordering fails, deliberately.

**Wrong turn, caught only by actually trying it.** I verified the change by
temporarily adding a 15th preset — and it *failed*, on a check I had not thought
of: "unknown-preset message unchanged" compared the whole error string, which
ends with the full preset list and therefore grows whenever one is appended.
Reasoning about the scoping change would not have found that; running the probe
did. The check now compares the message *preamble* exactly and separately
asserts every frozen name is still listed. Both directions are verified:
appending passes with the new preset named as uncovered, inserting at position 2
fails and prints the position-by-position divergence.

### 2026-08-18 19:35 · dev @ aea08ea · R5611351 · Claude Code (Opus 5), session 71af9711

**Repointed the three `fig_*_allStd` presentation figures at the new production
run.** `data/param_space/run_all_aug_14_26_12_14` (11 sensitivity levels, 4
param-space levels) → `run_all_aug_14_26_17_25` (25 levels throughout, reps
halved to 25 per d8886e2/a5033e6). One `data_root` line changed in each of
`Fig_param_space_allStd.m`, `Fig_sensitivity_analysis_allStd.m` and
`Fig_sfa_EOC_allStd.m`; all three were then re-run through the MATLAB MCP
session, which regenerated the .fig/.png/.svg deliverables and rewrote each
folder's `README_*.txt` and `git_provenance.txt` itself.

Nothing else needed touching — the new run has the same directory shape
(one `param_space_*`, seven `1D_sensitivity_*`, one `tau_sensitivity_*`
subfolder), so the scripts' glob-based subfolder discovery picked up the
renamed `nLevs_25` folders with no change. The "folder already exists"
warnings from `save_some_figs_to_folder_2` are expected: the figures are
saved next to the scripts, overwriting the previous run's copies.

### 2026-08-18 20:38 · dev @ e5f7f44 · R5611351 · Claude Code (Opus 5), session b2c13259

**A run directory now describes itself in prose.** New
`src/write_run_parameters_md.m` takes a `data/param_space/run_all_<dt>/` path
and writes `parameters.md` into it: preset and run mode at the top, a table of
what the preset set, the per-analysis sweep axes and timings, the four
adaptation conditions, and one unified alphabetical list of every parameter in
effect — preset, run mode and class default alike — as run in the `sfa_and_std`
condition, each tagged with what set it. Wired into `run_all_analyses.m` right
after stage 3 (in a `try`, so a failure to describe a run cannot be the last
thing a finished overnight run says), and backfilled onto
`data/param_space/run_all_aug_14_26_17_25`.

**The design point worth keeping.** It reads the *run*, never the source.
`saveobj` persists `psa.model_defaults`, and every sub-script sets that as
`merge_struct(preset_defaults, cfg.model)` — so subtracting the fields of
`cfg.model` (reconstructed by re-calling `analysis_run_config` with the recorded
`run_mode`) recovers the preset **as it was at run time**. Today's
`srnn_param_preset` is consulted only to diff against that, and any
disagreement is flagged as a ⚠ *Preset drift* block. That is what makes the tool
safe on old directories, where the named preset no longer matches the source.
`analysis_run_config` is likewise reconstructible from saved data: its only
preset dependence is `preset_is_stochastic`, which reads `sigma_u_noise` out of
`model_defaults`.

**Wrong turns worth recording.**

- `ParamSpaceAnalysis2.from_dir` looked like the obvious loader and is the wrong
  one — with no results on the object it calls `load_condition_results`, pulling
  in every 4 MB per-condition `.mat`. A plain `load` + pick-by-`isa` reads the
  47–215 KB psa and nothing else.
- `srnn_property_info` is `Static, Access = private`, so it is not reachable
  from a plain function; `class_default` and `effective_param` are the public
  way in.
- `SRNNCellTypePairs` cannot be default-constructed, so recognising "this is
  just the class default" needs a reference model built with the preset's five
  required constructor arguments. The consequence is real and is now stated in
  the file's own Caveats: scalar **aliases** of those arguments (`f_E` from `f`,
  `mu_EE_relative`/`sigma_EE_relative` from the tilde blocks) compare equal to
  that reference and therefore read as *class default* when the preset is what
  actually set them.
- `reps` is both a PSA grid axis and a real scalar property of the model
  classes. The first cut reported a `reps` axis for `param_space`, which has
  none; the axis branch is now gated on actual membership in `grid_params`, and
  the analyses with no reps axis honestly report the model property's value of 1.
- Two structs with the same fields render identically as a one-line field list,
  so the first drift table flagged `input_config` as changed while showing the
  same text in both columns. `fmt_diff`/`struct_diff_leaves` now name the moved
  leaf (`intrinsic_drive: 0 → 0.1`).

Verified on three directories through the MATLAB MCP session: the Aug-14
production run (`sra1` reported rather than the table's `ode45`, which is the
check that proves the stochastic-solver swap is being read from what the model
got and not from the run-mode table; `replot_sensitivity_*` appears only under
Caveats), the July-06 `SRNNModel2` run (pre-`resolved_defaults`, so a short
parameter list plus one caveat per analysis rather than an error), and an empty
directory (a complete file whose every section says "None recovered"). A
synthetic manifest naming the wrong preset confirmed the drift table is live
rather than vacuously passing.

---

### 2026-08-18 21:37 · dev @ fc6af94 · R5611351 · Claude Code (Opus 5), session b651fa7a

**New presentation figure: sensitivity medians collapsed across conditions**
(`scripts/presentations/Stability_Manuscript/fig_sensitivity_medians/`).

TR asked for a compact counterpart to the `_allStd` sensitivity sheets: instead
of tiling the four adaptation conditions as columns and drawing the full rep
distribution per panel (28 panels over two sheets per metric), overlay the four
conditions as **median lines only** in one axes per swept parameter. Dropping
`level_of_chaos` leaves six parameters, which fit one 2 x 3 sheet per metric.
Two sheets produced: `Fig_Sensitivity_LLE_medians` and
`Fig_sensitivity_mean_rate_medians` (png/svg/fig each), from the same source run
the `_allStd` figures use (`data/param_space/run_all_aug_14_26_17_25`). No
simulation re-run.

Decisions worth recording:

- **The replot pipeline is bypassed entirely.** `replot_sensitivity` ->
  `plot_sensitivity` -> `assemble_sensitivity_figure` exists to build the
  `imagesc` sheets; recovering median lines from its saved `.fig` files would be
  fragile. The medians are computed straight off the saved PSA objects, from the
  same `ParamSpaceAnalysis2.collect_level_values` that `plot_sensitivity` medians
  internally, so the two figures cannot disagree about what "the median" is.
  That helper was **moved from `methods (Static, Access = private)` to the public
  static block** — no behaviour change, but it means a `clear classes` is needed
  after pulling this.
- **f_E is labelled as an E:I neuron COUNT**, `100:400 / 250:250 / 400:100`, not
  the reduced ratios `1:4 / 1:1 / 4:1` the `_allStd` sheets use. TR's call: with
  `mu_EI`/`mu_IE` held fixed the sweep is really varying counts, and the reduced
  ratio hides the `n = 500` those counts come from.
- **LLE y-window kept at `[-1.75, 1.75]`** to match the `_allStd` sheets, so
  medians that leave it are clipped rather than rescaling every panel. This run's
  LLEs span roughly p1 = -10 to p99 = +3.7, so No Adaptation runs off the bottom
  on four of the six panels. TR chose clipping over autoscaling.
- **Percentiles are plumbed generally** (`pcts` / `band_pcts` at the script top,
  `prctile` rather than `median`), so adding an IQR band later is a two-line
  change; only the median is drawn now.
- **Four helpers extracted** from `Fig_sensitivity_analysis_allStd.m`'s local
  subfunctions into standalone files under `scripts/run_all_analyses/replot/`:
  `preset_default_values`, `apply_percent_axis`, `mark_default_value`,
  `save_figure_stable`. Bodies unchanged. `Fig_sensitivity_analysis_allStd.m` was
  **deliberately left untouched** so its committed figures cannot shift; its local
  copies still shadow these inside that file. Deleting them is a safe follow-up
  the next time that figure is regenerated.

Two things that went wrong and how they were diagnosed, so they are not
re-derived:

- **Tick labels came out with characters missing in the PNG** (`700` -> `)0`,
  `+50%` -> `+%`). The obvious reading — tiles too close, labels colliding — is
  **wrong**: `ax.XTickLabel` reads back correct, the `.svg` is clean, and
  re-exporting the same live figure at 300 dpi renders it perfectly while 600 dpi
  reproduces the corruption identically. Both corrupted labels sat at the same
  horizontal pixel position in different panels, i.e. a tile seam in
  `exportgraphics`' high-resolution rasteriser. Filed as ISSUE-013; worked around
  here by widening the figure 1200 -> 1300, which is a coincidence fix, not a
  repair.
- **The `1` y-tick of the mean-rate sheet was hidden under the panel letter.**
  `AddLetters2Plots` at the `_allStd` figure's `HShift/VShift = -0.03/-0.04` puts
  the letter exactly where the topmost tick label sits when the ticks are pinned
  to `[0 1]`. Now `-0.075/-0.05`.

Also added a dated note to **ISSUE-012** (default marker missing on the `_allStd`
sheets): the marker is present in the currently committed `_allStd` PNGs and the
extracted `mark_default_value` draws it correctly on all six panels of both new
sheets, so the issue looks stale — but `_allStd` was not re-run, so the status is
left OPEN rather than closed on indirect evidence.

---

### 2026-08-18 22:00 · dev @ 4810b17 · R5611351 · Claude Code (Opus 5), session b651fa7a

**Tick-label changes on both sensitivity figures, and the export bug that came
out from under them.**

TR asked for two label changes, applied to `fig_sensitivity_medians` and to
`fig_sensitivity_analysis_allStd`, both regenerated:

- Network Size ticks `100 / 500 / 1000` instead of `100 / 400 / 700 / 1000`.
  With three ticks they also fit horizontally, so the rotated labels the `n` row
  used to need are gone.
- The percent axes label their reference `+0%`, not `0%`. One line in
  `apply_percent_axis` — `sprintf('%+g%%')` already signs zero, so the fix was
  deleting the `labels(tp == 0) = {'0%'}` override. Every mu axis and the
  Synaptic Gain axis now reads as a signed departure scale end to end.

`Fig_sensitivity_analysis_allStd.m` was edited this time rather than left alone:
its four local subfunctions were deleted in favour of the copies extracted last
session into `scripts/run_all_analyses/replot/`, so the `+0%` change reaches both
figures from one place. Its "Local helpers" section is now a pointer comment.

**The export bug (ISSUE-013) is the part worth reading.** Changing the ticks
moved which label sat on the seam, so the corruption reappeared on a different
label after every regeneration. Four wrong turns, each of which looked like a fix
for one round:

1. *Widen the figure 1200 -> 1300.* Made `700` render; the next run ate `500`.
   Any figure-size change is this same coincidence.
2. *Force `Renderer = 'painters'`.* No effect. The `.svg` is clean because it
   goes through the vector path, not because of the renderer property.
3. *Use `print -dpng -r600` instead of `exportgraphics`.* Eats the same label in
   the same place, so it is not exportgraphics-specific.
4. *Assume a fixed seam at x = 4096 and cap the raster at 4000 px.* Fixed the
   middle column and broke the right-hand one. The seam tracks the end of a
   rasteriser *tile*, so it moves with the image width rather than sitting at one
   absolute pixel.

What actually settled it: the labels read back correct from `ax.XTickLabel`, the
`.svg` is clean, and bisecting the width gave **3520 px clean / 3759 px damaged**.
So `save_some_figs_to_folder_2` now measures the PNG it just wrote and re-exports
below `png_max_px = 3400` when needed, warning with the substituted dpi. Wide
sheets land at ~290 dpi instead of 600; narrow figures never take the branch and
are bit-identical to before.

A second, genuinely separate clipping problem was hiding underneath it: the
sweeps run edge to edge, so the outermost tick sits exactly ON the axis limit and
half its label overhangs the axes and is clipped. That is what turned `+100%`
into `00%` in the right-hand column, and it is not the seam — it survived every
resolution change. Both scripts now pad the x-limits by 3.5% each side
(`x_pad_frac`), which fixes it in data space instead of chasing figure margins.
`tiledlayout` `Padding`/`OuterPosition` were both tried first and neither helped;
the measured `TightInset` showed the label was never overflowing the canvas.

Verified by reading all six regenerated PNGs at native resolution: every tick
label complete on both medians sheets and all four allStd sheets.
