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
        FD-003 ──> FD-004 chained-pool STD ──> FD-005 chain unit tests
        FD-004 ──> FD-006 display and reservoir code
        FD-004 ──> FD-007 chain presets ──> FD-009 MTS reruns
        FD-004, FD-005 ──> FD-015 per-neuron spread on the chain
```

FD-003 comes before any change to the model class: the reference has to be
captured from the code as it stands.

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

- **Status:** open (decision made, implementation to be planned)
- **Serves:** block A
- **What:** add the chained vesicle-pool model as an opt-in per route beside
  the current product, which stays the default and bit-identical. A route
  with one component is unchanged and is the 1TS condition; no reserve
  variable exists for it. Equations (C1) to (C3) in the manuscript repo's
  `adaptation_normalization.md`. Touches the route configuration, the
  dynamics, the analytic Jacobian, the matrix-free Jacobian product, the
  state unpacking, and the route signature. `std_zero_floor` is refused for
  chain routes.
- **Done when:** FD-005's tests pass and the FD-003 regression is unchanged
  for every existing preset.
- **Depends on:** FD-003.
- **Result:** _(pending)_

### FD-005 Unit tests for the chain

- **Status:** open
- **Serves:** block A
- **What:** analytic Jacobian against finite differences; Jacobian product
  against the assembled Jacobian; a very large reserve reproducing the single
  pool; steady state against (C4); default routes bit-identical to before; a
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
  two figure sets per run, three-condition and four-condition; regenerate
  the preset fixture.
- **Done when:** the bundle runs at `fast`; preset tests pass.
- **Depends on:** FD-004.
- **Result:** _(pending)_

### FD-015 Per-neuron spread on the chain's parameters

- **Status:** open (requested and design confirmed by Tom 2026-11-10; see
  block A of the to-do list, "Per-neuron spread in the depression
  parameters"). Allowed on all three kinds of route: chain, one-component
  and multi-component product (Tom, 2026-11-11). Still under discussion:
  whether it is set in the preset or per condition (plan, section 9).
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
  bit for bit as today**.
- **Done when:** a zero spread reproduces the non-jittered chain bit for
  bit; the draws are reproducible from the seed and leave the connectivity,
  setpoints and SFA ladders unchanged; the analytic Jacobian and the
  Jacobian product pass their tests with a nonzero spread; the regression
  bundle is unchanged for existing presets.
- **Depends on:** FD-004, FD-005.
- **Result:** _(pending)_

### FD-016 A seed for the joint parameter-space sample

- **Status:** open (decided by Tom 2026-11-11)
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
