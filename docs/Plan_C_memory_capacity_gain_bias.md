# Plan C: memory capacity across gain and bias, after three fixes

## Context

Block C of the final draft (`../StochasticPlasticDynamicalSystemPaper/docs/final_draft_to_do.md`,
section C) must show that the memory-capacity result comes from adaptation and
not from tuning: across a range of connectivities, not at one network. The
manuscript's present result (Figure 6, `fig:sfa_longest_timescale`) is one
network. Tom's instruction (2026-11-12, evening): do C next; he will be asleep;
the lead makes the decisions, logs them and proceeds until C is done, using the
six worker agents.

**The experiment** (Tom): scale the overall connection weights with
`level_of_chaos` (gain) and, as a second axis, the nonlinearity's bias `S_c`
(the manuscript's $a_0$), for the four conditions, and look at the distribution
of memory capacity. The question: is chained MTS adaptation the most robust,
that is, consistently higher memory capacity than single-timescale (1TS) or no
adaptation? Expected: no adaptation needs precise tuning. Shown as 2-D maps,
gain by bias, memory capacity as colour, one per condition.

**Three fixes come first** (Tom: yes to all three): (F1) the model seeds its
random numbers without naming the generator, so a pool worker and an ordinary
session build different networks; (F2) the `essential` 1-D sweeps ran 16 jobs
per condition where 8 were meant; (F3) the 1TS SFA constant must be the
neuron's fastest MTS rung.

State at the start: `main` at `691ae9d`; block A built and closed; nine block B
pilots done; quick suite and stages tier last run at `d9cd4cf`; the full tier
deferred by Tom; nothing running.

## Rules that govern the night (all already in CLAUDE.md or Tom's messages)

- The lead runs nothing in MATLAB over about 30 s. Workers run, gate and watch.
- No questions to Tom. Decisions and their reasons go under the issue in
  `docs/ISSUES.md` or in `docs/WORKLOG.md`; paper-affecting ones also into the
  to-do list, marked as the lead's.
- Agents never run `medium` or larger. Block C's runs are `essential` and
  `fast` sizes of the new stage (section 4).
- Tests for a change: unit tests, the tolerance and ultrarapid tiers, review of
  the diff; the quick suite and stages tier once before the main run. No
  mutation checks (nothing here adds arithmetic to the dynamics or Jacobian).
- A reference is replaced only deliberately, with the reason in the commit.
- Workers: one instance per name, `isolation: worktree`, brief starts with
  `git merge --ff-only main` and names the commit, pool stated on line 1, at
  most 12 pool workers on the machine. A worker's report text comes back in
  its hand-back; the lead saves reports and builds PDFs.
- Nothing is pushed. No dates or deadlines are written into the repos' docs.
- The lead's context is nearly full and will be compacted during the night.
  **After every merge the lead writes the state (what merged, what is
  running, what is next, decisions) into the top entry of `docs/WORKLOG.md`**,
  so the work can be picked up from the worklog, the issue and this plan.

## 1. The three fixes (step C0), in parallel where files allow

| Fix | Worker | Files | What |
|---|---|---|---|
| F1 + F3 | `matlab-fable-1`, pool 2 | `src/model/SRNNCellTypePairs.m` (and the other seeding sites the exploration lists), tests, references | see below |
| F2 | `matlab-opus-1`, pool 2 | the sweep drivers / `analysis_run_config`, `test_essential_mode` | the reps axis has exactly `n_reps` integer levels in every mode |

Code sites, from the exploration (line numbers at `691ae9d`):
- F1, bare `rng(seed)` calls to name the generator: `SRNNCellTypePairs.m`
  831 (noise), 1699 (W, and x0 later from the same stream), 1973 (setpoints),
  2572 (input); `SRNNModel2.m` 1387, 1545, 1635, 3242; `SRNN_ESN_reservoir.m`
  126 (input weights), 806 (stimulus); `stimulus/dc_staircase_stimulus.m` 53;
  `stimulus/stim_target_indices.m` 56; `analysis/run_stim_neurons.m` 101.
  The SFA ladder, the synapse draws and `lyapunov_topk` already name it. The
  stage functions that call `rng(0, 'twister')` before building are why the
  stages tier is already generator-independent.
- F3: `build_tau_a_matrix`, line 1907, `w = 1` for K = 1 becomes `w = 0` (the
  fast-end normal, column 1 of the draw); the property comment at line 101 and
  the check "K = 1 is the slow draw alone" in
  `test_SRNNCellTypePairs_tau_a_spread` change with it. Presets affected: the
  `regression_pairs_*` family, `regression_stable_steps_*` and the `mu7chain`
  family (all have `tau_a_spread` 0.25 and a one-timescale condition).
- F2: `psa.add_grid_parameter('reps', 1:ctx.n_reps)` in
  `run_sensitivity_analysis.m` 125 and `run_tau_sensitivity_analysis.m` 179.
  A two-element value is read as a `[min max]` range and expanded to
  `n_levels` values; three or more is an explicit list. So only `n_reps = 2`
  is wrong. Fix: pass reps so it is always an explicit list. The rep value
  never enters the seed (seed = grid index x 100), so every other mode keeps
  its jobs and seeds.

**F1, the generator.** Every place that seeds the global stream names the
generator: `rng(seed, 'twister')`. Decision: pin to twister, the generator an
ordinary session already uses, so everything built in a client (the model
tests, the tolerance and ultrarapid tiers, the examples, run13) keeps its
numbers, and only what was built on pool workers changes (the sweeps; the
`essential` run's sweep results are then not reproducible to the bit, which is
recorded). Done when: a model built in a client and one built on a pool worker
from the same seeds are bit-identical (new unit test, needs a pool of 2);
`rebuild_model` of a sweep job reproduces its stored values without setting a
generator; the ultrarapid and tolerance tiers pass unchanged.

**F3, the 1TS SFA constant.** With one SFA timescale the neuron's constant uses
the same standard normal as the fastest rung of a multi-timescale ladder on the
same seed (today it uses the slow-end normal). Done when: on one seed the 1TS
constant equals the MTS ladder's first rung, neuron by neuron (unit test, in
`test_SRNNCellTypePairs_tau_a_spread`); zero spread unchanged. This moves
every single-timescale result with `tau_a_spread > 0`, so references move:
the worker shows that the leaves that move are exactly those of
single-timescale conditions with a spread and that all others are
bit-identical, then the lead re-captures deliberately (tolerance suites and
exact tiers one at a time, reason in each commit). The tolerance anchor is
re-captured too, as a new anchor from this commit, with the old one's
agreement on all unmoved leaves stated in the commit; no tolerance is loosened.
Both equation files get the one-sentence statement of the rule.

**F2, the reps axis.** `reps` gets integer levels `1:n_reps`. Existing modes
have levels equal to reps (fast 5 and 5), so their job counts and seeds do not
change; the worker confirms that per mode.

Gate after C0 (one worker, `matlab-opus-2`, pool 6): quick suite, stages tier
with `settings.n_workers` ignored, the four tolerance suites, both ultrarapid
tiers, on merged main. Stages or sweep references that move because of F1 or
F3 are looked at once, then re-captured deliberately.

## 2. The stage (step C1): memory capacity and growth rate at every sweep point

New stage in `src/analysis/` (working name `run_mc_gain_sweep`), built by
`matlab-fable-2` while C0 is in progress on other files. FD-010 in ISSUES.

- **Inputs:** a preset (the `mu7chain` noise-free memory-capacity twin), a list
  of gain levels (`level_of_chaos`), a list of bias levels (`S_c`, the centre
  of the per-neuron setpoint draw; the spread `sigma_S_c` stays), seeds, the
  conditions (all four), a size mode.
- **Per job** (one network seed at one gain and bias): build the four
  conditions on the same seeds, so the comparison is paired exactly as
  `run_memory_capacity` pairs it (`verify_shared_build` on the first job).
  Compute memory capacity with the synaptic readout by the SAME functions
  `run_memory_capacity` uses (no second implementation of the readout), and
  $\lambda_1$ by a short separate run of the same network.
- **Stored per job and condition, scalars and one short curve:** total memory
  capacity, memory horizon, the $R^2$-against-delay curve, $\lambda_1$ with
  the memory-task input on (decision: this is the quantity memory theory
  relates to, and the one plotted) and with no input (recorded beside it),
  mean rate of E and I, the fraction of neurons pinned above 0.9 and below
  0.01 (so a low memory capacity can be read as saturated, silent or
  chaotic). Nothing else is kept.
- **Memory.** Today `SRNN_ESN_reservoir.run_reservoir_esn` integrates the
  whole 760 s in one call and returns every state at every step (152,000 rows
  by about 4,000 states, about 5 GB, then per-variable copies: the "16 GB per
  trial"), and only afterwards keeps one column per hold. There is no lean
  option. The readout needs only the synaptic output of each neuron at the
  last sample of each 0.3 s hold (2,533 columns). So the stage gets a lean
  path: integrate hold by hold (or a few holds at a time) with the same
  fixed-step integrator, keep the readout column, carry the state forward.
  It is a NEW method beside `run_reservoir_esn` (opt-in; the existing path and
  everything that calls it unchanged), and the scoring code is shared, not
  copied. 12 jobs then fit at once.
  **Check that makes it trustworthy:** with the fixed-step integrator and
  noise off, segmenting changes no step, so on a small network and on the
  n = 500 preset network the lean path must give the same $R^2$ curve and
  total memory capacity as `run_memory_capacity` to round-off; any difference
  is found and explained before a sweep is run.
- **$\lambda_1$** does not come from the memory-capacity run
  (`run_reservoir_esn` does not go through `run()`): a separate 30 s `run()`
  of the same reservoir with the top-K estimator, `lya_T_interval` set
  explicitly (its default comes from `T_range`, not from the task's input).
- **Known and kept:** the input weights `W_in` use a fixed seed
  (`rng_seed_input = 3`), so they are the same in every trial; the network
  and the input sequence vary with the seed. Stated in the report.
- **Cost, measured today:** about 2 worker-minutes per network and condition
  at n = 500 (24 min for 3 trials of 4 conditions, one at a time).
- `ParamSpaceAnalysis2` has no hook for an extra per-job quantity, so the
  stage is a plain `parfor` over a job table (with `pool_ensure`, `vprintf`,
  provenance and the run-directory conventions of the other stages).
- **Protocol:** the durations of `run_memory_capacity` unchanged (10 s
  washout, 600 s training, 150 s test, delays to 15 s): a shorter training
  window makes the readout under-determined (Tom). Noise off, `sra1`.
- **Sizes** (the stage's own table; agents run the first two only):
  `essential` = the pilot grid of section 3; `fast` = the main grid of
  section 4; `medium` and `production` defined for Tom (more seeds, finer
  grid) and not run.
- **Resumable by job file:** each finished job writes its own small file and a
  rerun skips finished jobs, so a crashed worker loses minutes, not hours.
- Tests: `test_mc_gain_sweep` (tiny network, 2 x 2 grid, 1 seed: runs, pairs,
  skips finished jobs, equals `run_memory_capacity` at the preset point);
  ultrarapid tiers unchanged.

Beside it, `matlab-opus-3` (pool 0): the two E:I presets within the `mu7chain`
bundle, chained from the noise-free twin, mean weights 8.4 / 5.6
(excitation-dominated) and 5.6 / 8.4 (inhibition-dominated) on both
postsynaptic types, sum held at 14 (the to-do list's proposal, adopted); preset
fixture regenerated deliberately; and the figure functions of section 5 drawn
first on a small synthetic result so they are ready when data arrives.

## 3. Pilot (step C2): where the no-adaptation network works

`essential` size, one Opus worker with all 12 pool workers, about half an hour.
No adaptation and chained MTS only, 2 seeds, a coarse grid: gain over a wide
log-spaced range (about 0.1 to 3, nine levels) and bias `S_c` over about -0.2
to 0.8 in seven levels (the preset's value is 0.2). Purpose: find where
$\lambda_1$ crosses zero and where memory capacity peaks without adaptation,
confirm the time per job and the memory per job, and set the main grid so that
it brackets the no-adaptation peak in both directions with room on each side.
The lead reads the pilot's table and fixes the main grid (logged).

## 4. Main runs (step C3)

Sized from the pilot's measured time per job so that the compute stays modest
; the default below assumes about 2 worker-minutes per
network-and-condition and 12 jobs at once.

| Run | Grid | Seeds | Conditions | Jobs x conditions | Estimate on 12 workers |
|---|---|---|---|---|---|
| C3a gain by bias, standard network | 9 gains x 7 biases | 3 | 4 | 756 | about 2 h 10 |
| C3b gain only, excitation-dominated (3 : 2) | 9 gains, default bias | 3 | 4 | 108 | about 20 min |
| C3c gain only, inhibition-dominated (2 : 3) | 9 gains, default bias | 3 | 4 | 108 | about 20 min |

If the pilot shows jobs are cheaper, seeds go to 5 before the grid is refined
(seeds give the distribution Tom asked for). If dearer, the bias axis drops to
5 levels before seeds drop below 3. One Opus worker launches and watches each
run (blocking call through a wrapper, marker files, foreground waits,
relaunch once on a transient failure; the stage resumes from its job files).
C3a first; C3b and C3c after it on the same pool.

## 5. Analysis, figures, report (step C4)

One Opus worker for figures, one for statistics and the report, after C3a.

- **Main figure (three conditions; a four-condition version beside it):** the
  2-D maps, gain by bias, memory capacity as colour on one shared scale, one
  panel per condition; below them the same maps for $\lambda_1$, with the
  $\lambda_1 = 0$ contour drawn on both rows.
- **Gain curves at the preset bias:** memory capacity against gain,
  $\lambda_1$ against gain, memory capacity against $\lambda_1$, per
  condition, seeds as points with the median and interquartile band; the
  3 : 2 and 2 : 3 networks as further rows.
- **Robustness statistics per condition:** the distribution of memory
  capacity over the whole gain-by-bias map (median, interquartile range, 5th
  percentile); the share of the map where memory capacity is at least half
  that condition's own best and at least half the best of any condition; the
  peak, its location and the width at half maximum along gain (the Gaussian
  fit of the agreed design, against log gain, with the measured half-maximum
  width beside it when the fit is poor); paired differences between conditions
  per grid point across seeds.
- **Report:** `docs/notes/Block_C_results_<date>.md` and PDF, in the pilots'
  format (question, method, result, caveats, figures, values), and a draft of
  the Results paragraph and figure caption for the manuscript, written to
  `../StochasticPlasticDynamicalSystemPaper/docs/` as a proposal. Decision:
  `Manuscript5.md` itself is not edited overnight; Tom reviews the proposal.
- The figure functions take `cfg.conditions`, the manuscript style and
  `save_figure_stable`, and are registered in `mu7chain_config` as a stage and
  figures that run only when selected (`cfg.stages`), so existing bundles are
  unchanged.

## 6. Order and who

| Step | What | Worker (pool) | Beside it |
|---|---|---|---|
| C0 | F1 + F3 | `matlab-fable-1` (2) | F2 `matlab-opus-1` (2); C1 stage `matlab-fable-2` (4); presets and figure scaffolding `matlab-opus-3` (0) |
| C0 gate | quick, stages, tolerance, ultrarapid on merged main; re-captures | `matlab-opus-2` (6) | |
| C1 merge | stage merged after its tests; rebased on the fixes | lead | |
| C2 | pilot | one Opus worker (12) | |
| C3a | main map | one Opus worker (12) | figure functions tested on the pilot's data |
| C3b, C3c | E:I gain sweeps | same worker (12) | statistics on C3a |
| C4 | figures, statistics, report, manuscript proposal | two Opus workers (2 each) | timeline update (Opus, no MATLAB) |
| Close | ISSUES (FD-010 done), WORKLOG, to-do checkboxes and decisions, timeline | lead | |

Merges that touch `SRNNCellTypePairs.m` go one at a time. The stage worker is
told that F1 and F3 are landing in the model class under it and rebases before
its final tests.

## 7. Decisions already taken by the lead for this plan (to be logged)

1. Generator pinned to twister, not to Threefry (keeps every client-built
   number; changes pool-built sweeps, which are rerun anyway).
2. $\lambda_1$ is measured with the memory-task input on and plotted that way;
   the no-input value is stored beside it.
3. E:I variants hold $|\mu_E| + |\mu_I|$ fixed: 8.4 / 5.6 and 5.6 / 8.4.
4. The bias axis is `S_c`, the centre of the per-neuron setpoint draw; the
   per-neuron spread is kept.
5. The Gaussian fit is against log gain, with the measured half-maximum width
   reported beside it.
6. The 2-D map is run on the standard network only; the two E:I variants get
   the gain axis at the preset bias (cost; a fuller E:I axis was already
   deferred in the to-do list).
7. `Manuscript5.md` is not edited; a proposed paragraph and caption are
   written beside the to-do list.
8. Noise stays off for block C, as the agreed design has it.

## 8. Verification

- F1: client-built and worker-built models bit-identical; a sweep job rebuilt
  in a client reproduces its stored $\lambda_1$ and rate.
- F2: `essential` 1-D sweep has 8 jobs per condition; other modes' job tables
  unchanged.
- F3: 1TS constant equals the MTS first rung per neuron on one seed; unmoved
  leaves bit-identical; moved references re-captured with reasons.
- Stage: lean memory capacity equals `run_memory_capacity` at the preset
  point; a rerun skips finished jobs; paired build verified.
- Gates on merged main before C2: quick suite, stages tier, four tolerance
  suites, both ultrarapid tiers.
- Runs: every job file present; no failed job; the pilot's and main run's logs
  free of errors; the no-adaptation peak inside the grid, not on its edge.
- Figures looked at by the worker and by the lead (PNG) before the report is
  built; the report's numbers traced to the saved table.

## 9. If something goes wrong

- A fix that cannot pass its gate by the time the stage is ready: the stage
  proceeds on the fixes that did pass, and the gap is logged in the worklog
  and stated in the report (F3 in particular changes the 1TS results, so if it
  slips, the 1TS column is marked as on the old rule).
- The lean memory path disagrees with `run_memory_capacity`: fall back to the
  existing path with three jobs at once, cut the grid as in section 4, and say
  so.
- Compute running long: finish C3a at the cost of C3b and C3c; the report is
  written from whatever has finished, stating what has not.
