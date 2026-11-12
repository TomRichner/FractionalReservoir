# Issues: code work toward the final draft

The list of discrete pieces of work in **this repo** needed for the final
draft of the manuscript. Started 2026-11-10.

## How this file is used

Three documents, three jobs:

| Document | Answers | Lives in |
|---|---|---|
| `../StochasticPlasticDynamicalSystemPaper/docs/final_draft_to_do.md` | **What and why.** The feedback, the decisions, the manuscript-level tasks (blocks A to G). | manuscript repo |
| `docs/ISSUES.md` (this file) | **What exactly, in the code.** One entry per piece of work that can be finished and checked. | this repo |
| `docs/WORKLOG.md` | **What happened, session by session.** What was done, what was verified, what is next. | this repo |

Rules:

- **One issue is one checkable piece of work.** It names the to-do block it
  serves, what "done" means, and what it depends on. If an issue cannot be
  given an acceptance test, it is not ready and stays in "To be planned".
- **Issues are created by planning with Tom**, not on an agent's own
  initiative. A defect noticed during other work is raised (in the reply and
  in the worklog) and becomes an issue only when Tom agrees. This is the
  repo's scope rule.
- **Status is one of:** `open`, `in progress`, `blocked` (say on what),
  `done` (give the commits), `dropped` (say why). Update the status in the
  same commit as the work, or in the commit that closes the session.
- **IDs are permanent.** `FD-001`, `FD-002`, ... Never renumber or reuse; a
  dropped issue keeps its number.
- **Manuscript text is not tracked here.** Edits to `Manuscript5.md` are
  tracked in `final_draft_to_do.md`. An issue here ends at "the run exists
  and the figure is on disk".
- **Decisions are not recorded here.** They go in `final_draft_to_do.md`; an
  issue links to the block that holds the decision.
- **Reasoning stays in commit messages**, as before. An issue's "Result"
  line points at commits; it does not restate them.

A lead agent orchestrates and worker agents build (TR, 2026-11-10; the
earlier arrangement, one agent doing all the building, was dropped the same
day). The lead reviews, runs the regression gates and merges. Who builds
which issue, in what order, is in `docs/Plan_A_chained_vessicles.md` for
block A.

## Order

Dependencies, not priority:

```
FD-001 test audit ─┐
FD-002 fast regression bundle ─┴─> FD-003 regression reference
        FD-017 pool module ──> FD-020 tolerance tier ──> FD-021 refactor
        FD-021 ──> FD-004 chained-pool STD ──> FD-005 chain unit tests
        FD-004 ──> FD-018 chained facilitation ──> FD-015 per-neuron spread
        FD-004 ──> FD-006 display and reservoir code
        FD-004 ──> FD-007 chain presets ──> FD-009 reruns
```

FD-003 and FD-020 come before any change to the model class: their
references have to be captured from the code as it stands.

**Decided against, 2026-11-11 (Tom): "bit-identical" as the requirement
on existing behaviour.** Issues FD-004, FD-005, FD-015 and FD-018 were
first written to require that the product model and existing presets
stay bit for bit as they were, and the plan turned that into leaving the
product code textually untouched and branching per model. Tom rejected
it: it put the convenience of a byte-exact test ahead of clean, general
code. The requirement is now that existing behaviour gives the **same
numbers to round-off, shown by the tolerance regression tier** (FD-020;
Plan A, sections 2 and 6). Byte-exact comparison remains the gate only
for changes that should move no number at all (the pool module, configs,
figures), and the exact references are re-captured once after the
refactor (FD-021). Do not reintroduce a bit-identity requirement on the
model's arithmetic.

## Issues

### FD-001 Audit the existing tests

- **Status:** done (2026-11-10)
- **Serves:** all blocks (the safety net for everything below)
- **What:** `scripts/tests/` holds about 60 files of uncertain state. Run
  every test, classify each as passing, failing because the test is stale,
  or failing because the code is wrong. Remove stale tests. Write missing
  ones for the current pipeline (`SRNNCellTypePairs`, `SRNN_ESN_reservoir`,
  `ParamSpaceAnalysis2`, presets).
- **Done when:** every remaining test passes or has a recorded reason; a
  single entry point runs them all and prints a summary; the audit's
  findings (including any real defects found, not fixed) are in the worklog.
- **Depends on:** nothing.
- **Result:** commits `ba0508c` (entry point), `a356bd4` (removals),
  `32bbf5e` (test fixes), `155b5ef` (new tests). All 58 remaining tests pass.
  `run_all_tests('quick')` is 52 tests in about 3.5 minutes;
  `run_all_tests('slow')` is 6 tests in about 18. No defect was found in
  `src/`: the four failing tests were all wrong on the test side. Five
  run-and-plot scripts that asserted nothing were removed. Three tests were
  added (`test_std_product_model`, `test_pairs_reproducibility`,
  `test_psa_pairs_end_to_end`) and `test_sample_hold_mc` was ported to the
  current reservoir class. Done by a subagent; its findings and the gaps it
  listed are in the worklog entry of 2026-11-10.

### FD-002 Fast regression bundle

- **Status:** done (2026-11-10), but slower than intended; see Result
- **Serves:** all blocks
- **What:** a preset and a config that exercise the same stages as the
  manuscript's bundle (`mu7revisedMedium`) on a network small and short
  enough to run in minutes. The existing `fast` mode takes about 67 minutes
  for the full pipeline, which is too slow to run before and after every
  change.
- **Done when:** the bundle runs end to end without error; its wall time is
  recorded; it follows the existing bundle conventions.
- **Depends on:** FD-001 (so the preset fixture and tests are in a known
  state before a preset is added).
- **Result:** commit `a24334b`. Four presets (`regression_pairs_*`,
  `single_neuron_regression`: the manuscript network at n = 100, indegree =
  20), `regression_fast_config(tag)` and `regression_fast_run(tag)`. First
  run, tag `reference`: all ten analysis stages OK in 31.7 min and 36 of 36
  figures OK in 18.3 min, so about 50 min, not minutes. The 1-D sweeps
  shrank with the network; the joint sample (9 min), memory capacity (6 min)
  and numerics (4 min) are sized by run mode, which a config cannot change.
  A smaller run mode would cut it but means editing every stage; not done.

### FD-003 Regression references and comparison tests, in three tiers

- **Status:** done (2026-11-10)
- **Serves:** all blocks
- **What:** TR, 2026-11-10: the regression has to be faster, so tiers.
  **Ultrarapid** (`regression_rapid_run`, about 35 s): the three classes
  called directly, no stage functions. **Stages** (`regression_stages_run`,
  2 min): three stage functions and four figures. **Full**
  (`regression_fast_run`, about 50 min): every stage and figure of the
  FD-002 bundle. All are reduced by the generic `regression_digest` and
  compared by `regression_compare`.
- **Done when:** each tier has a tracked reference made from a clean tree;
  its check passes on unchanged code; and two runs on the same machine show
  how reproducible it is, with any ignore patterns set from that.
- **Depends on:** FD-002.
- **Result:** tooling `7fcf533`, `10879e1`, `cae664d`. References, all made
  from a clean tree before any change to the model class:
  - ultrarapid `bef29cb`: 994 leaves, exact, no ignore patterns; two runs
    bit-identical; `test_regression_rapid` passes in 36 s (quick suite).
  - stages `58ed2ce`: 1245 leaves and 18 files; wall-clock fields ignored;
    two runs otherwise bit-identical; `test_regression_stages`, 2 min
    (slow list).
  - full `c6b71ec`: 4583 leaves, 176 run files, 183 figure files; a second
    complete analysis run was bit-identical apart from wall-clock fields,
    the commit recorded in the manifest, and the joint parameter-space
    sample; `test_regression_bundle`, about 50 min (slow list).
  The joint sample is NOT covered: `ParamSpaceAnalysis2` seeds it with
  `rng('shuffle')`, so its 128 points differ on every run. Raised, not
  changed.

### FD-004 Chained-pool depression in `SRNNCellTypePairs`

- **Status:** in progress (step 1 briefed to `matlab-fable-1`, with FD-005)
- **Serves:** block A
- **What:** add the chained vesicle-pool model as an opt-in per route beside
  the current product, which stays the default and gives the same numbers
  to round-off. Built as a second implementation of the route interface
  that FD-021 puts in place, not as a branch at each site. A route with
  one component is the 1TS condition: equivalent to the reserve held full,
  implemented with no reserve variable. Equations (C1) to (C3) in the manuscript repo's
  `adaptation_normalization.md`. Touches the route configuration, the
  dynamics, the analytic Jacobian, the matrix-free Jacobian product, the
  state unpacking, and the route signature. `std_zero_floor` is refused for
  chain routes.
- **Done when:** FD-005's tests pass; the tolerance tier (FD-020) passes
  against its anchor; and the exact references re-captured after FD-021
  pass unchanged (adding a model must not move a number on routes that do
  not use it).
- **Depends on:** FD-003, FD-020, FD-021.
- **Result:** _(pending)_

### FD-005 Unit tests for the chain

- **Status:** open
- **Serves:** block A
- **What:** analytic Jacobian against finite differences; Jacobian product
  against the assembled Jacobian; a very large reserve reproducing the single
  pool; steady state against (C4); default routes unchanged (the regression
  gates of FD-004, not a bit-identity test here); a
  one-component route identical under either model.
- **Done when:** all pass and are in the FD-001 entry point.
- **Depends on:** FD-004 (written alongside it).
- **Result:** _(pending)_

### FD-006 Display and reservoir code for the chain

- **Status:** open
- **Serves:** block A
- **What:** code outside the dynamics that recomputes the product itself:
  the reservoir class's route signature and time-series plot, the STD
  steady-state figure, the generated equation and parameter tables, the
  depression panel of the example time series. Each should ask the model for
  the route's depression factor.
- **Done when:** the paper figures build on a chain preset and show the right
  quantity; the product presets' figures are unchanged.
- **Depends on:** FD-004.
- **Result:** _(pending)_

### FD-007 Chain presets and conditions

- **Status:** open
- **Serves:** block A
- **What:** a preset bundle on the current baseline network with the MTS
  condition using the chain (reserve size 2, reserve refill 5 s), named
  `sfa3_stdc2`, in a **four-condition** preset that also keeps the product
  MTS (`sfa3_std2`) for the side-by-side; title and colour entries (the
  chain takes the product's current colour, the product becomes green);
  two figure sets per run, three-condition and four-condition, selected by
  `cfg.conditions`; `full_adaptation_condition` learns the chain's name
  (today it would silently return `sfa1_std1` for the paper's three
  conditions); the SFA-timescale sweep runs both MTS conditions (Tom);
  regenerate the preset fixture. Plan A section 4.4.
- **Done when:** the bundle runs at `fast`; preset tests pass.
- **Depends on:** FD-004.
- **Result:** _(pending)_

### FD-015 Per-neuron spread on the chain's parameters

- **Status:** open (requested and design confirmed by Tom 2026-11-10; see
  block A of the to-do list, "Per-neuron spread in the depression
  parameters"). Allowed on all three kinds of route: chain, one-component
  and multi-component product (Tom, 2026-11-11). Set in the
  preset, not in a condition (Tom, 2026-11-11). Built after chained
  facilitation, so it also covers `stf_spread`.
- **Serves:** block A
- **What:** a per-neuron log-normal spread on the chain's time constants
  (releasable-pool refill and release, reserve refill) and on the reserve
  size, sigma = 0.25 as for SFA. Values stay positive and a longer time
  constant gets a proportionally larger variance. Depression has no spread
  today: each route shares one row of time constants. Follow the SFA
  pattern (`tau_a_spread`, `tau_a_matrix`, its own seeded stream with the
  RNG state restored). Confirmed design: the four quantities drawn
  independently; one draw per presynaptic neuron, shared by its routes;
  the single-timescale condition jittered with the same releasable-pool
  draws; the reserve kept slower than the releasable pool with a minimum
  separation; **the default is no spread, so every existing preset runs
  as today** (to round-off, by the tolerance tier; one code path that
  always carries per-neuron values is preferred to a branch on whether a
  spread is set).
- **Done when:** a zero spread reproduces the non-jittered model within
  the tolerance tier; the draws are reproducible from the seed and leave the connectivity,
  setpoints and SFA ladders unchanged; the analytic Jacobian and the
  Jacobian product pass their tests with a nonzero spread; the regression
  gates pass for existing presets.
- **Depends on:** FD-004, FD-005, FD-018.
- **Result:** _(pending)_

### FD-016 A seed for the joint parameter-space sample

- **Status:** done (2026-11-11)
- **Serves:** blocks A to D (comparable reruns), and the regression
- **What:** `ParamSpaceAnalysis2` draws the joint sample with
  `rng('shuffle')` and `randperm`, so its 128 grid points differ on every
  run. Add a seed for that draw as an option. **Unset, the behaviour is
  today's**, so existing bundles and the regression references do not
  move. `run_param_space_analysis` passes it from the config; new bundles
  set it.
- **Done when:** two runs with the same seed sample the same points and
  give identical results for that stage; with the seed unset two runs
  differ as now; the seed is recorded with the run; the ultrarapid, stages
  and quick gates pass unchanged.
- **Depends on:** nothing; must land before the four-condition `fast` run
  of Plan A step 6.
- **Result:** commit `68da56f` (worker `matlab-fable-2`), merged as
  `fafd990`. `ParamSpaceAnalysis2.sample_seed`, set from a bundle's
  `cfg.param_space_seed`. On merged main: quick 54 of 54 (ultrarapid
  exact), stages tier pass, `test_psa_sample_seed` in the quick suite. The
  seed is in `psa_object.mat` always and in the summary only when set (so
  the full reference does not move). No existing config sets it; the
  chain bundle's will.

### FD-017 One place for the parallel pool: `pool_*` functions and a session allocation

- **Status:** done (2026-11-11; Plan A step 0a).
- **Serves:** every block; the multi-agent working arrangement
- **What:** a small module in `src/util/` through which all pool handling
  goes, so that no code opens the default profile or a hard-coded size.
  - The session's allocation lives in an **environment variable of the
    MATLAB process** (`setenv`), set in one step by `pool_set(N)`. It
    survives `clear all` and `clear classes`, and a launcher can set it
    before MATLAB starts.
  - `pool_set(N)` records N and opens exactly N workers on a named profile
    (`Processes`); `pool_set(0)` means no pool. `pool_ensure()` returns the
    session's pool, opening it if down. `pool_restart()` replaces
    `restart_parpool`: it reads the size of the pool that is up and reopens
    the same number. `pool_release()` closes it. `pool_status()` shows who
    holds what across the machine.
  - A machine-level registry (one small file per session: process id, a
    label such as the agent's name, worker count, start time), so agents can
    see each other. **Visibility only: the 12-worker budget is not enforced
    in code** (Tom: exceeding it bogs the CPU down but does not crash; the
    aim is to leave a couple of cores free for the CLI and remote desktop).
  - **With no allocation set, library code runs serially with a loud
    warning; `run_all_paper_analyses` and `make_all_paper_figures` refuse to
    start.**
  - The default profile stays `Just1` (Tom): it also stops legacy code from
    auto-starting a large pool. The profile name and the budget figure each
    appear once, in the module.
  - `restart_parpool`, `wait_for_parpool`, and the pool-opening lines in
    `run_lyapunov_spectrum`, `run_transient_gain` and
    `run_numerics_verification` go through the module. The six worker-agent
    prompts and `MATLAB_worker_agents.md` already describe `pool_set(N)`
    with a `parpool('Processes', N)` fallback for the time before the
    module exists (rewritten 2026-11-11); when it lands, the fallback and
    the usage note's "before the module" column are removed.
  - An option to run a session's workers at lowered priority, so the
    desktop, the CLI and remote desktop stay responsive (Tom, 2026-11-11).
    **Pinning workers to cores is out of scope for now**
    (`docs/notes/set_parpool_cpu_affinity.m` is Tom's earlier, Windows-only
    script, kept for reference).
- **Checked 2026-11-11** (a subagent, in the lead's MATLAB, R2026a): the
  variable round-trips, is private to the MATLAB process (a child process
  inherits it, an unrelated shell does not see it), and survives
  `clear functions`; `clear all` and `clear classes` were tested on pool
  workers, not on the client, and it survived there. A `persistent` is
  reset by `clear functions` unless the function is `mlock`ed. What the
  wrappers must allow for:
  - **Workers copy the client's environment once, when the pool starts.**
    A later `setenv` on the client is not seen by them. So the variable is
    client-side only, and `pool_set` must restart the pool when N changes
    while one is up.
  - It does not outlive the MATLAB process, so every new session sets it.
  - Reading a pool's size and profile and reopening the same works
    (`gcp('nocreate').NumWorkers`, `.Cluster.Profile`); delete took about
    2.5 s and reopening two workers about 5.5 s (15 s the first time).
  - From the operating system, a pool worker is a `MATLAB.exe` with
    `-parallelserver` on its command line whose parent process id is the
    client's `feature('getpid')`. So `pool_status` can count real workers
    per client without relying on the registry alone.
- **Done when:** no `parpool(` call remains outside the module; a session
  set to N opens N and restarts at N; an unset session runs a `parfor`
  serially with the warning and the two entry points refuse; `pool_status`
  lists two concurrent sessions correctly; the stages and full regression
  tiers pass unchanged (pool handling must not move a number).
- **Depends on:** nothing.
- **Result:** commit `ce771b7` (worker `matlab-fable-1`), merged as
  `ff0bd02`. `src/util/pool_{set,ensure,restart,release,status}.m` with
  `pool_config` and `pool_registry` behind them; `restart_parpool`
  deleted; the only `parpool(` call left is in `pool_ensure`. On merged
  main: quick 53 of 53; `test_pool_module` pass (slow list, about 50 s);
  two sessions listed correctly by `pool_status` (lead 10, worker 2);
  stages tier pass with `settings.n_workers` ignored; **the analysis half
  of the full tier (tag `fd017`, 38 min on 9 workers, run after FD-016
  was merged too) shows no value difference**: one row, the file
  `memory_capacity/working_changes.patch`, which is in the reference only
  because that stage of the reference run was made with uncommitted
  changes in the tree (its `git_provenance.txt` says `dirty: true`; the
  patch is empty). The figure half was not run, by agreement with Tom: the
  module changes no figure code beyond the refusal `test_pool_module`
  covers. Found on the way, and now in `CLAUDE.md`: serial runs differ
  from pool runs at round-off; `settings.n_workers` is a compared leaf
  captured at 14. Not tested: the non-Windows paths.

### FD-018 Chained facilitation, after chained depression

- **Status:** open (decided by Tom 2026-11-11; design written and reviewed,
  Plan A section 4.6; parameter values chosen)
- **Serves:** a later paper. **Facilitation does not appear in this one.**
  Built now so the model does not have to be refactored twice.
- **What:** the release-coupled, serial form: option C (section 5) of
  `docs/notes/Depression_facilitation_model_options_2026-11-11.md`, built
  as an extension of the chain's route model (`route_pool`), not as a
  third model (Tom; Plan A section 4.0). The full design is Plan A
  section 4.6, to be written by the lead and reviewed by Tom first. Release
  probability multiplies both the output and the rate at which the
  releasable pool is used; in the multiple-timescale form a calcium trace
  drives release probability. Per route, like depression, so different
  cell-type pairs can carry different short-term plasticity (several E and
  several I types are expected later). A per-neuron spread on its time
  constants (`stf_spread`), set in the preset like `std_spread`.
- **Done when:** its own unit tests pass (equations by hand, steady state,
  finite-difference Jacobian, Jacobian product); **with facilitation off,
  the chain agrees with the references saved for it before this work
  began** (Plan A step 5: within tolerance on the tolerance tier, and
  exactly on the chain's exact reference unless the lead re-captures it
  with the tolerance evidence in hand); the existing gates pass; today's
  facilitation on product routes (`single_neuron_stf`) gives the same
  numbers to round-off.
- **Depends on:** FD-004, FD-005 and the chain's own reference. **Built
  straight after the chain** (Tom, 2026-11-11), before the spread (FD-015),
  which then covers `stf_spread` too. The calcium trace is per route.
- **Result:** _(pending)_

### FD-020 Tolerance regression tier

- **Status:** done (2026-11-11). Steps 0b and 0c are merged and **Tom
  fixed the tolerances as proposed**. Acceptance test agreed with Tom in
  the execution plan; Plan A sections 2.2 and 7.3, steps 0b and 0c.
- **Serves:** block A; every later edit to the model's arithmetic
- **What:** in `scripts/tests/regression/`: `regression_tol_run` (builds
  the cases, returns raw arrays), `regression_tol_compare` (element-wise;
  per leaf the largest absolute and relative difference, reported on a
  pass too; integers, logicals and text exact),
  `make_regression_tol_reference(name)`, `regression_tol_check(name)`,
  `test_regression_tolerance`. `regression_compare` is not reused for the
  first two layers: it compares digest statistics, not elements. Three
  layers:
  - **Pointwise.** At states stored in the fixture (from the stable runs,
    from chaotic runs of `regression_pairs_*` in all three conditions with
    noise off and on, and from the "options" network of
    `regression_rapid_run`: three cell types, geomean, zero floor, STF,
    dual STF, route scale): `dynamics_fast`, `compute_Jacobian_fast`,
    `jacobian_times`, and a **canonical parameter view** (scaled `W`,
    `tau_d`, `c_eff`, `tau_a_matrix`, setpoints, `route_fields` per route).
    Proposed rtol 1e-12. The view is canonical because `params` changes
    shape in the refactor: the fixture stores physics, not the struct.
    Leaves present in the anchor must match; new leaves are reported, not
    failed.
  - **Stable, end to end.** New presets `regression_stable_*`, chained
    within the regression bundle: the regression network with the gain
    lowered until every condition is stable with a clear margin (all three
    estimators), driven by the step staircase, noise off and on; plus a
    stable variant of the options network. Trajectories, `plot_data`, all
    three Lyapunov estimators, transient gain, a 12-job sweep, memory
    capacity. Proposed rtol 1e-9; discrete outputs exact.
  - **Chaotic and bundle.** The existing rapid, stages and full digests
    through `regression_compare` with a loose `rtol`: missing, extra,
    shape, file and figure rows fail; value rows are reported as a band.
  - **Mutation checks** (step 0c, in a worktree, never merged). Each must
    fail the tier: a `tau_rec` off by 1e-6 relative; a `tau_a` rung off by
    1e-6; the `(b, a)` block dropped; the `(x, b)` block dropped; the
    `(g, x)` block dropped; one block of `jacobian_times` scaled by
    `1 + 1e-9`; the recovery term's sign flipped; `c / K` replaced by `c`;
    a change confined to one route (`I->E`); the geomean gain omitted from
    the Jacobian only; the route scale dropped from `params.W`. The list,
    adapted, is re-run on the refactored code after FD-021 and on the
    chain when its reference is captured.
- **Done when:** (1) two captures from one tree are bit-identical; (2) the
  tier passes against its own anchor at tolerance zero; (3) the anchor is
  captured from a clean tree whose `SRNNCellTypePairs.m` blob equals the
  one at the FD-003 references, with the commit and blob hash stored in
  it; (4) every mutation fails the tier at the tolerances then fixed, and
  the list, results and margins are in `docs/notes/`; (5) the rapid, quick
  and stages gates are unaffected, the preset fixture regenerated
  deliberately; (6) run time recorded, in the quick suite if under a
  minute; fixture size reported (target under about 20 MB, sparse
  Jacobians).
- **Depends on:** FD-017. Comes before any edit to the model class.
- **Result:** step 0b `85bda9b`, `222a43c`, `6eeb904`, `dfc5da2`
  (merged `4dc3d84`): the tier, two stable presets, the anchor (8.7 MB,
  2370 leaves, from model-class blob `2b0118d4`); the lead's capture on
  merged main is bit-identical to it. Step 0c `988feae` (merged
  `24e515b`): 55 mutations, all caught by compared values; note
  `docs/notes/Tolerance_tier_mutation_checks_2026-11-11.md`. Tolerances
  fixed by Tom as proposed: pointwise 1e-12, stable 1e-9, QR 1e-5, memory
  capacity 1e-6, finite-difference Jacobian 1e-7. **Left as they are, by
  Tom's decision:** the two thin spots the mutation checks found. The
  reservoir class is guarded only by memory capacity (at 1e-6), and the
  assembled Jacobian's facilitation and geomean blocks have one witness,
  the stored pointwise Jacobian. The lead proposed adding anchor cases
  for both before the refactor; Tom did not take that up and said to
  proceed to step 0d. (The lead misread "proceed" as approval and started
  a worker on it; the worker was stopped with nothing committed.) Both
  are recorded in the note's "thin spots". An anchor case can no longer
  be added from the untouched model class once step 0d merges.

### FD-021 Refactor behind one route interface

- **Status:** merged and gated (2026-11-11); set to done when the third
  exact reference, the full one, is re-captured (its run is under way).
  Interface design approved by Tom; Plan A sections 2.1 and 4.0.
- **Serves:** block A
- **What:** product depression and today's facilitation move behind the
  route interface of Plan A section 4.0, as `src/model/routes/route_product`
  resolved through `route_model`. No new model, property or behaviour.
  The public signatures of `dynamics_fast`, `compute_Jacobian_fast`,
  `jacobian_times`, `unpack_and_compute_states`, `get_params`,
  `route_fields` and `routes_identical` do not change; the state layout
  does not change; `params` keeps the fields read outside the class.
- **Done when:** the tolerance tier passes against the anchor and a fresh
  pre-edit capture; `test_route_models` (new: each model's `rhs` against a
  hand loop, `linearize` against central differences of its own `rhs`,
  `steady_state` against a long integration) passes; the quick suite
  passes unedited; the four functions contain no depression- or
  facilitation-specific arithmetic; `rhs`, the assembled Jacobian and
  `jacobian_times` (K = 15) are timed before and after at n = 100 and
  n = 500, and **the lead decides** whether a slowdown above 10 % is
  acceptable or worth fixing, recording the decision here (Tom,
  2026-11-11: this was to go to him; he wants the lead to decide); the
  adapted mutation list fails the tier on the refactored code; then the
  lead re-captures the three exact references once, with the reason in the
  commit.
- **Depends on:** FD-020.
- **Result:** worker `matlab-fable-1`, commits `d3d1e82`..`964e8b2`, merged
  `7d49430`. `src/model/routes/route_model.m` and `route_product.m`;
  `test_route_models`; `regression_time_model`;
  `regression_tol_mutations_post_refactor`. On merged main, lead's
  MATLAB: tolerance tier pass against the anchor (240 of 1749 float
  leaves moved; worst pointwise 4.8e-16 against 1e-12; QR 7.6e-8 against
  1e-5) and against the pre-edit capture; quick 55 pass with
  `test_regression_rapid` the one expected failure; band checks of the
  old rapid and stages references, no structural difference; 71 of 71
  post-refactor mutations caught by values (worker). Ultrarapid reference
  re-captured `8d7d2c9`, stages `dfc93f0` (which now ignores
  `settings.n_workers`, so it passes on a pool of any size).
  **Decisions by the lead** (reasons in the merge commit and Plan A 4.0):
  1. *Three departures from the interface as designed are accepted*
     (`params.route` a cell; `linearize` with four outputs; a model
     declares `P.coupling`). Each was made for measured speed; without
     the third `jacobian_times` was 30 % slower at n = 100.
  2. *Timing is accepted, nothing to fix.* After over before, measured by
     alternating the class with a copy of the pre-edit class: between
     0.97 and 1.05 for `dynamics_fast`, the assembled Jacobian,
     `jacobian_times` (K = 15) and build plus `get_params`, at n = 100 and
     n = 500. The same code against itself reads 0.98 to 1.05, so this is
     noise.
  3. *`active_routes` is accepted as a workaround.* Three existing tests
     switch `std_zero_floor` on the `params` struct after `get_params`,
     and the floor is now compiled into each route, so the routes are
     recompiled when the two disagree. The cleaner fix is for those tests
     to set the property before `build`; **suggested as a later tidy-up,
     not done** (the refactor had to pass the suite unedited).
  4. *The full reference is re-captured from a run compared loosely
     against the old one first* (tag `reference_fd021`), as the other two
     were.
  **Carried into step 1** (Plan A 4.0, "What step 1 inherits"): the
  full-`A` branch of `jacobian_times` is under no tier case until the
  chain's reference exists; `route_fields` cannot yet tell a chain from a
  product route; the legacy delegations are product-only.

### FD-022 Remove the remnants of resuming a `ParamSpaceAnalysis2` run

- **Status:** open (Tom, 2026-11-11: resume is not needed; remove what is
  left of it, in the future). Not scheduled; it does not block block A.
- **Serves:** housekeeping; keeps the documentation true
- **What:** `CLAUDE.md` and the class header say a sweep "can be resumed by
  loading the PSA and calling `run()` again". The FD-016 worker reported
  that this no longer works: `run()` appends a new dated folder to
  `output_dir` and calls `generate_grid` again, so it starts a new run in a
  nested directory and never finds the old `temp_batches`. Tom's reading
  is that resume existed and was not carried through an earlier refactor.
  **The lead has not verified the report**; the first step is to confirm
  it. Then remove the claim from `CLAUDE.md` and the class header, and
  remove code that exists only for resuming (for example the skip of a
  batch whose file already exists in `run_batched_simulation`), keeping
  whatever the batch files are still needed for (consolidation, memory).
- **Done when:** the report is confirmed or corrected in the commit
  message; no document or comment says a run can be resumed; no code path
  exists only to resume one; the ultrarapid, quick and stages gates pass
  unchanged (removing dead code must not move a number).
- **Depends on:** nothing. `ParamSpaceAnalysis2.m` is one file, so not
  beside other work on it.
- **Result:** _(pending)_

## To be planned

Named so they are not lost. Each needs a planning session with Tom before it
becomes an issue with an acceptance test.

- **FD-008 Frozen-reserve control** (block A). **Dropped 2026-11-10:** Tom
  settled that 1TS and MTS are not strength-matched at steady state, so the
  control it was proposed for is not needed. Was: blocked on the decision
  whether to run it.
- **FD-009 MTS reruns at `medium` or better** (block A): sweeps, memory
  capacity, stability and transient stages on the chain preset. **Agents do
  not run these**: only Tom starts them and says when they are done (Tom,
  2026-11-11). Plan A ends with one `fast` run on four conditions.
- **FD-010 Gain sweep of growth rate and memory capacity** (block C),
  including the no-adaptation pilot and the 3:2 and 2:3 networks.
- **FD-011 Stimulation: adaptation is necessary** (block B): the 2 × 2 plus
  frozen-state experiment, more than one kind of aberrant activity, the
  parameter-space version, and which network carries it.
- **FD-012 Combined stimulation figure** (block B; merges manuscript
  Figures 7 and 8).
- **FD-013 Closing figure** (block D).
- **FD-014 Paper figure of 1TS versus MTS depression** (block A), from the
  demonstration scripts.
- **FD-019 A best-practices repo for multi-agent MATLAB** (Tom's idea,
  2026-11-11; after the final-draft work). If the worker agents, the `pool_*`
  functions and the session environment variable work well here, lift
  them into a separate repo with examples and documentation so other
  projects can use the same arrangement (MATLAB MCP server, one private
  session per agent, Parallel Computing Toolbox), and test it on a Mac.
