# Work Log

Chronological record of work on this repo: plans executed, bugs found, bugs
fixed, and analysis runs performed. Newest entries at the **bottom**.

This exists because the work happens on more than one machine and across many
sessions, and neither git history nor the `data/` run directories capture the
*reasoning* — what was tried, what turned out to be wrong, and what is still
open.

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

| hostname | OS | notes |
|---|---|---|
| `R5456622` | Windows 11 | 6 physical cores, 64 GB RAM, default pool profile `ten_workers_one_thread` |
| _(other)_ | — | the second machine; add its hostname on first use. 14 workers. |

If the asset-tag hostnames get unwieldy, a friendlier per-clone alias can be set
without touching tracked files, since `.git/config` is local to each clone:

```bash
git config --local machine.name "mayo-desktop"
git config --get machine.name
```

### A note on merge conflicts

Two machines appending to one chronological file will collide at the same place
in git. Entries are therefore self-contained and separated by `---`, so a
conflict is resolved by keeping both sides and putting them in date order —
never by choosing between them. If this becomes tiresome, the alternative is
`docs/worklog/YYYY-MM-DD-<hostname>.md`, one file per session, which cannot
conflict at all.

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
changed **two variables at once** — fresh pool *and* 6→5 workers — and was on
the wrong cluster profile. It validates nothing.

**Fixed** — `c6055ea`: `src/restart_parpool.m`, called before each of the three
analyses in `run_all_analyses`. Uses `parallel.defaultProfile` rather than a
hard-coded worker count, since the two machines differ (10 here, 14 there).
Framed as **insurance, not a diagnosis**: it bounds whatever workers accumulate
over a multi-hour run, at ~30 s per stage.

**Also fixed** — `d58b7fe`: stopped tracking `debug.log`, swept in by a blanket
`git add -A`.

## Open items

**Bugs, not yet fixed:**

- **`load_results` loses the vector-parameter lookup.**
  `param_space_summary.mat` never stores `vector_param_lookup` or
  `vector_param_config`, so after `load_results` both are empty and
  `effective_param(res, 'tau_a_E')` returns the **level index** (1…13) instead
  of the value (5…60 s) — **silently**, no error. Plotting (`:936`) takes the
  same wrong branch and produces a wrong x-axis. Workaround: load
  `psa_object.mat` directly, which round-trips correctly.
- **`load_results` does not restore `model_class`**, though
  `summary_data.model_class` is saved. A loaded `SRNNCellTypePairs` run reports
  `SRNNModel2`, which corrupts `effective_param`'s class-default fallback and
  `same_config`'s cross-class refusal.
- Both are missed by `test_psa_saveload`, which exercises `saveobj`/`loadobj`
  (the path that **works**) and never calls `load_results` (the path the replot
  scripts actually use).

**Interpretation, not a bug — but it affects every sweep:**

- **The LLE has a floor at `−1/max(tau_a)`.** No preset overrides `tau_a`, so
  all use the class default `logspace(log10(0.25), log10(10), n_a)` → slowest
  timescale 10 s → eigenvalue **−0.1** exactly. Whenever the recurrent dynamics
  contract faster than that, the largest Lyapunov exponent is the slowest
  *linear adaptation mode* and says nothing about the network. Evidence: all 20
  top QR exponents in `test_benettin_vs_qr` within 0.002 of −0.1; `sfa_only` /
  `sfa_and_std` pinned at −0.101 across unrelated sweeps. The depression
  equivalent is `−(1/tau_rec + r/tau_rel)` ≈ −0.62, matching
  `best_presets.md`'s *"`std_only` sits flat near −0.65"*.
  **The tau stage sweeps the floor itself** (`tau_a_E(end)` over [5, 60]), so it
  is its own test: plot LLE against `−1/tau_a_E_max`. Supporting evidence
  already in hand — the minimum successful LLE was **−0.2159** against a
  predicted floor of −1/5 = −0.2.

**Unfinished work:**

- `medium2` `tau_sensitivity` and `param_space` still need a clean run.
- `best_presets.md` numbers need re-deriving post-window-fix, and the
  `[0.5, 1.5]` range narrowing re-confirming.
- `dev` is **13 commits ahead and unpushed** (`git push -u origin dev`).
