# Work log: sessions toward the final draft

One entry per working session, **newest first**. Started 2026-11-10.

## How this file is used

- **Read the top entry before starting.** It says where the last session
  stopped and what it left unverified. Then read `docs/ISSUES.md` for what is
  open.
- **Add an entry when a session ends** (or at a natural stopping point in a
  long one). Insert it directly below this section.
- **Keep entries short.** An entry is a pointer to commits and a statement of
  state, not a narrative. Reasoning and wrong turns go in the commit
  messages, which remain the record; link to them.
- **Say what was verified and how**, and separately what was not. "Ran X, got
  Y" and "did not run Z" are both required. A later session should be able
  to tell from the entry alone what it can rely on.
- **Do not edit old entries**, except to add a dated correction line under
  one that turned out to be wrong.

Entry template:

```
## YYYY-MM-DD HH:MM (session <id>)

**Issues touched:** FD-00x (status before -> after), ...
**Commits:** <hash> <subject>; ... (repo if not this one)
**Done:** ...
**Verified:** ... (what was run, what it showed)
**Not verified / known gaps:** ...
**Raised, not acted on:** ... (defects or ideas for Tom to decide)
**Next:** ...
```

---

## 2026-11-12 01:40 (session cc2136bc-6794-4ea6-962d-9eb831682265, continued: refactor and chain merged; steps 3 and 5 running)

**Standing instructions from Tom during this stretch** (in CLAUDE.md,
"Current work"): the lead carries Plan A through to completion without
stopping to ask, deciding and recording decisions in ISSUES.md or here;
rules written in advance are guides for judgement, so a tolerance may be
adjusted after an investigation that explains it, with the reason
documented.

**Issues touched:** FD-021 -> done; FD-004 -> done; FD-005 -> done;
FD-006 -> in progress (reservoir part); FD-007 -> in progress.

**Commits (this repo):** FD-021 `d3d1e82`..`964e8b2` merged `7d49430`;
references re-captured: ultrarapid `8d7d2c9`, stages `dfc93f0`, full
`cc51085`; FD-004 / FD-005 `2b5f6c2`..`97aad1c` merged in the commit
before `f5759cd`; `f5759cd` the equations file and the issue results;
`6727ee5`, `b05a94e` Tom's two standing instructions; `da9cd86`,
`08fa37a` Plan A 4.0 as built and the lead's decisions on the refactor.

**Done:**
- Step 0d: product depression and facilitation behind the route interface
  (`src/model/routes/route_model.m`, `route_product.m`).
- The three exact references re-captured once, each after the tolerance
  evidence and a loose comparison with the old one.
- Step 1: the chain as `route_pool`, its tests, `std_route_factor`.
- `Equations_stability_paper.md` states the chain.

**Verified (lead's MATLAB, merged main):**
- Refactor: tolerance tier pass against the anchor (worst pointwise
  4.8e-16 against 1e-12) and the pre-edit capture; quick 55 pass with the
  ultrarapid tier the one expected failure; band checks of the old rapid
  and stages references structurally clean.
- Full bundle on the refactored code (40.5 min analyses, 36 of 36
  figures): against the old full reference no structural difference; 605
  of 4583 leaves moved, median 5.7e-16; the 13 rows above 1e-9 were
  looked into and are eigenvalue ordering (the sets agree to 5.5e-11) and
  a 1e-12 quantity moving in its fourth digit.
- Chain: tolerance tier pass against the anchor; quick 57 of 57 with the
  ultrarapid tier byte-exact.

**Decisions by the lead** (details under the issues): the refactor's
three departures from the designed interface, its timing (0.97 to 1.05,
noise) and the `active_routes` workaround are accepted; the chain's STD
panel and its check 9 as built; the bundle is `mu7chain`, its presets
carry `stdSpread0p25` in their names before the spread exists; the full
reference was digested from a copy of the run without one locked folder.

**Not verified / known gaps:**
- **A file is locked by a process the lead could not identify:**
  `data/regression_fast/reference_fd021/replot_sensitivity_nov_12_26_00_39/figures/sensitivity_mu_EI_relative_f_3.fig`.
  It made `fig_sensitivity_analysis_allStd` fail its cleanup once (the
  figure pass was re-run and passed). Not the lead's MATLAB (no open file
  identifiers) and not its pool (released, still locked). The folder is
  left in place; delete `data/regression_fast/reference_fd021` when the
  lock is gone. If it recurs, look at what else reads `data/` (indexer,
  antivirus, another MATLAB).
- The full reference's reproducibility on the refactored code rests on
  one run; the next full run (step 10) is the second.
- The chain has no regression reference yet (step 5, running).

**Raised, not acted on:** the cleaner fix for `active_routes` (three
tests set `std_zero_floor` on the params struct); `route_pool` and
`route_product` each carry a copy of two small helpers; seven stages save
display titles with their data, so a retitled condition changes exact
references.

**Next:** review and merge, as they report: the reservoir class
(`matlab-fable-2`), the presets and bundle (`matlab-opus-2`; then
re-capture the stages reference for the retitled product condition), the
chain references (`matlab-fable-1`; re-capture on main, require
bit-identity). Then step 8 (chained facilitation, design in Plan A 4.6)
on the model lane, with the display work (FD-006, three Opus workers,
packages in the inventory note) beside it. Then FD-015, step 10, the
`fast` four-condition run, FD-014.

---

## 2026-11-11 19:45 (session cc2136bc-6794-4ea6-962d-9eb831682265, continued: steps 0a to 0c closed, 0d running)

**Issues touched:** FD-017 -> done; FD-016 -> done; FD-020 -> done;
FD-021 open -> in progress (step 0d briefed); FD-022 created (remove the
remnants of resuming a PSA run; open, not scheduled).

**Commits (this repo):** `68da56f` FD-016 (worker `matlab-fable-2`),
merged `fafd990`; `85bda9b`, `222a43c`, `6eeb904`, `dfc5da2` the tolerance
tier and its anchor (worker `matlab-fable-1`), merged `4dc3d84`; `988feae`
the mutation checks, merged `24e515b`; `ac4ab98`, `d1f1a0f` CLAUDE.md for
the sample seed, the timer lesson and the tolerance tier; `2054f0e`
FD-022; `5864ac2`, `d87ea6e`, `aa95aca`, `429101d` issue statuses, the
last one correcting `aa95aca`.
**Commits (manuscript repo):** `f409f1b` tolerances fixed.

**Done:**
- FD-016: `ParamSpaceAnalysis2.sample_seed`, set by `cfg.param_space_seed`.
- FD-020: the tolerance tier (`regression_tol_*`, 2370 leaves, 24 s, in
  the quick suite), two stable presets, the anchor; 55 mutations all
  caught by compared values; **tolerances fixed by Tom as proposed**.
- Tom's further decisions: short labels "SFA+STD" / "SFA+prodSTD"; figure
  layout is not a question for him; resume of a PSA run is not needed.

**Verified (lead's MATLAB, merged main):**
- After FD-017 and FD-016: quick 54 of 54; stages tier pass
  (`settings.n_workers` ignored); **analysis half of the full tier, tag
  `fd017`, 38 min on 9 workers: no value difference.** Its one row is the
  file `memory_capacity/working_changes.patch`, in the reference only
  because that stage of the reference run had a dirty tree.
- After the tolerance tier: quick 55 of 55 (314 s); the lead's capture is
  `isequaln` to the worker's anchor; zero-tolerance check passes on all
  leaves. Anchor blob `8f071b7a`; model class still blob `2b0118d4`.

**Not verified / known gaps:**
- The figure half of the full tier has not been run since FD-017.
- Two thin spots of the tolerance tier, left open by Tom's decision: the
  reservoir class is guarded only by memory capacity (1e-6), and the
  assembled Jacobian's facilitation and geomean blocks have one witness
  (the stored pointwise Jacobian). After step 0d merges, an anchor case
  can no longer be added from the untouched model class.
- The harness `regression_tol_mutate` was run 70 times by the worker and
  read by the lead, not run by the lead.

**Mistakes by the lead this session, so they are not repeated:**
- The first launch of the full tier did nothing for 30 minutes (a
  TimerFcn that used `evalc` with a `catch ME` inside an anonymous
  function). The right way is in CLAUDE.md.
- Tom answered one of two questions and said "proceed"; the lead took
  that as approval of its own recommendation to add anchor cases and
  briefed a worker. Tom corrected it; the worker was stopped, its
  worktree discarded, nothing committed. "Proceed" means the next step of
  the agreed plan.

**Raised, not acted on:** the rapid digest stores a copy of `params`
inside `lya_results`, so after the refactor the rapid reference will show
`missing` / `extra` keys as well as last-bit differences; handle at the
re-capture. `assert_route_redundancy`'s error message rounds the values
it prints, so two routes that differ by 1e-6 look identical. The earlier
entry's list still stands.

**Next:** when `matlab-fable-1` reports on step 0d (FD-021): read the
whole diff against Plan A 4.0; merge; `clear classes`; tolerance tier
against the anchor; quick suite (the ultrarapid tier is expected to fail
there, and only there); the timing table, to Tom if any path is more than
10 % slower; the loose-band check of the old exact references; then
re-capture the three exact references once from a clean tree, with
`settings\.n_workers$` stored as an ignore pattern, and say why in the
commit. Then step 1 (the chain, FD-004 and FD-005).

---

## 2026-11-11 17:10 (session cc2136bc-6794-4ea6-962d-9eb831682265, Plan A execution begins; entry written mid-session)

**Issues touched:** FD-017 open -> in progress -> merged (set to done when
the full-tier gate below reports); FD-020 and FD-021 written as full
issues; FD-020 open -> in progress (step 0b); FD-016 open -> in progress;
FD-007 and FD-018 reworded.

**Commits (this repo):** `53000f5` Plan A 4.0, FD-020 / FD-021, 4.6
outline; `5cf8cd8`, `6ed8476`, `f0b0941` Plan A 4.6 in full, Tom's
answers, the facilitation parameter values; `e737f68` the step 2 inventory
(worker `matlab-opus-1`); `7f52b12`, `96f1806` Tom's figure decisions;
`ce771b7` FD-017 (worker `matlab-fable-1`), merged as `ff0bd02`; `a5c5911`
CLAUDE.md, Plan A and the six worker prompts for the pool module;
`33480de`, `f6ae523` issue statuses with gate results.
**Commits (manuscript repo):** `6426c5f`, `6f2ec58`, `a1a4913`, `f76f198`
decisions into the to-do list; `cff93d0` the worker-agents note.

**Done:**
- Tom approved the execution plan (the route interface, FD-020 and FD-021
  as issues, the sequence). His decisions, all in the to-do list: state
  order kept; route models as function files in `src/model/routes/`;
  chained facilitation extends `route_pool`; the tau sweep runs both MTS
  conditions; facilitation ceiling jittered on `G - 1`; coupled
  facilitation allowed without depression; facilitation values (G 5,
  `tau_fac` tied to `tau_rel`, `tau_dec` 1.5 s, `tau_h` 0.5 s, spread
  0.25); product MTS green `[0 0.7 0]`; both figure versions for every
  three-condition figure, the single-neuron one included; four-condition
  layouts need not be polished; short labels "SFA+STD" and "SFA+prodSTD";
  layout is not a question for Tom.
- Step 0: pre-flight and the docs above.
- Step 0a: the pool module (`src/util/pool_*.m`), every pool call site
  routed through it, `restart_parpool` deleted, both entry points refuse
  without an allocation.
- Step 2: the display and conditions inventory,
  `docs/notes/Inventory_chain_display_and_conditions_2026-11-11.md`.

**Verified:**
- Pre-flight on `53000f5`: quick 53 of 53, 297 s.
- FD-017 on merged `ff0bd02`, lead's MATLAB: quick 53 of 53 (247 s, pool
  10); `test_pool_module` pass (53 s, session left as found); stages tier
  pass on pool 12 with `settings.n_workers` ignored; a second session
  (`matlab-opus-1`, pool 2) ran `pool_status` and listed both, 12 of 12,
  no warning. The lead read the whole diff.
- `SRNNCellTypePairs.m` is still blob `2b0118d4`, the one the FD-003
  references were made from.

**Not verified / known gaps:**
- **The analysis half of the full tier is running** (tag `fd017`, pool 6,
  started 17:02; result in `data/regression_fast/fd017.result.txt`). The
  figure half is skipped for FD-017 by agreement with Tom: the module
  changes no figure code.
- The pool module's non-Windows paths are untested.
- The facilitation section 4.6 is hand-derived, unsimulated.
- FD-020 (`matlab-fable-1`, pool 3) and FD-016 (`matlab-fable-2`, pool 3)
  are running and unreviewed.

**Raised, not acted on:**
- **`settings.n_workers` is a compared leaf** in the stages and full
  references, captured on 14 workers; the lead passes an ignore pattern
  until the re-capture after step 0d stores it.
- **Serial differs from a pool run at round-off** (466 leaves of the
  stages tier, 1e-13 to 1e-15), so the stages and full tiers need a pool
  up. Predates FD-017.
- **Worktrees start from `origin/main`**, so briefs begin with
  `git merge --ff-only main`. Nothing is pushed (Tom's call).
- The full tier takes 50 min; Tom asked whether 20 is possible. Measured
  breakdown and two options (cheaper figure export in the regression
  bundle; a run mode below `fast`) were discussed, nothing decided. It
  would be its own issue.
- From the workers: `write_run_parameters_md` prints only the first
  condition's synapse config and calls it common to all; seven stages
  freeze display titles at compute time; `run_overnight_queue`'s
  pre-flight does not check the pool allocation; `fig_to_do_in_future.md`
  and two archived scripts ask for 13 workers; nine crashed jobs sit in
  the `Processes` profile's job storage and are announced on every pool
  start; `wait_for_parpool` and `replot_sensitivity` print with bare
  `fprintf`.
- Citation check owed: release probability 0.2 and 1.5 s (Mongillo et al.
  2008, from memory); `tau_h` has no source.

**Next:** read `fd017.result.txt`, set FD-017 to done. Review and merge
FD-016 and FD-020 as they report; the lead re-captures the tolerance
anchor on merged main and requires bit-identity with the worker's. Then
step 0c (mutation checks), the tolerances to Tom, step 0d.

---

## 2026-11-11 15:10 (session 1766c926-8e9d-42f1-ba3a-5851de78596b, scripts/paper archive)

**Issues touched:** none. Tom asked which `scripts/paper` files produced the
manuscript's figures and to archive the rest; this is an exception to the
deferred reorganisation, recorded in the to-do list and in `CLAUDE.md`. He
then asked for the STD-matching bundles back and for `reproduce_paper_run`
to run the manuscript's run.

**Commits:** this repo, `d477013` scripts/paper: archive the 22 config/run
pairs the manuscript's figures do not use; then "reproduce_paper_run runs
the manuscript's bundle; STD-matching bundles back out of the archive".
Manuscript repo, `9a55969` and "To-do: STD-matching bundles stay in
scripts/paper; reproduce_paper_run runs the manuscript's bundle".

**Done:** 36 files (18 `*_config.m` / `*_run.m` pairs) are in
`scripts/paper/archive/`, moved with `git mv`, contents unchanged (44 were
moved at first; the 8 files of `stdScaled_fast`, `stdUsage_fast`,
`stdUsage_med`, `stdSingleMatched_fast` went back, being live for block A).
Seventeen `.m` files stay: the manuscript bundle (`mu7revisedMedium_config`
/ `_run`, `run_all_paper_analyses`, `make_all_paper_figures`,
`grouped_figure_registry`), the regression bundle, the STD-matching
bundles, and `paper_config` with `reproduce_paper_run`. The three `.md`
notes were left where they are. `reproduce_paper_run.m` now passes
`mu7revisedMedium_config()` to both entry points (it ran `paper_config()`,
the mu8p25 network); README and `CLAUDE.md` say so.

**Verified:** in MATLAB, all 25 `*_config` functions return the same preset
names before the move, after it, and after the move back; all 53 `.m` files
resolve by name to where they now sit. `test_mu7revised_bundle` and
`test_run_modes` pass in the final state; `test_verbose_levels`,
`test_preset_golden` and `test_preset_conditions` passed after the first
move. `reproduce_paper_run` was run on this machine: it built the
mu7revisedMedium cfg (medium, `data/mu7revisedMedium`,
`figs/mu7revisedMedium`) and stopped at
`run_all_paper_analyses:TargetNotEmpty`, as intended where the manuscript's
run already exists; `git status` shows nothing written under `data/` or
`figs/`. Code Analyzer: no issues.

**Not verified / known gaps:** `reproduce_paper_run` has not been run
through (7 h at medium; Tom's to start, and it needs
`data/mu7revisedMedium` moved aside first). No regression tier was run (no
file a tier reads was moved or edited). No archived `*_run.m` was executed.
A MATLAB session started before the move needs `setup_paths` again.

**Raised, not acted on:** `paper_config` still names the mu8p25 network, so
a bare `run_all_paper_analyses()` / `make_all_paper_figures()` is not the
manuscript's run; the README and `CLAUDE.md` sections that call it "the one
file to edit" describe that default. `README_mu7revised.md` mostly
describes archived bundles.

**Next:** unchanged from the entry below.

---

## 2026-11-11 (session 7ec5eb22-1e10-4e89-8903-1803563a876a, planning, second half)

**Issues touched:** FD-020 (tolerance tier) and FD-021 (refactor) named
under "To be planned"; FD-004, FD-005, FD-015, FD-018 reworded (no
bit-identity on the model's arithmetic); FD-017 gained a note on the agent
prompts. No code changed, so no regression tier was run.

**Commits (this repo):** `81f32e5` Plan A: generalise behind a tolerance
gate; pool module, then tolerance tier, first; `cb764a1` CLAUDE.md, README:
align with the revised Plan A; `6a463b9` CLAUDE.md, README, Plan A: how
ISSUES, WORKLOG and the worker agents are used; `e891de1` Worker agents:
stop opening the default pool profile; describe the pool module; `74abae4`
Docs: no dates govern the work; 1TS has no clamped variable; four
conditions, three in the paper; `849e591` ISSUES: bit-identity on the
model's arithmetic is decided against; FD-020, FD-021 named.
**Commits (manuscript repo):** `96cf126` To-do: general code behind a
tolerance regression; build order for block A; `6190c2f` To-do: how the
issues file, the worklog and the worker agents are used; `71ec7d9` Worker
agents note: pools after AutoCreate off and Just1 as default profile;
`4ea83c8` To-do: remove the dates; 1TS without a clamped variable; four
conditions analysed, three shown; stale lines; `57fff5a` To-do: the
tolerance tier and the refactor are FD-020 and FD-021.

**Done:** Tom's decisions, written into Plan A, CLAUDE.md, README, ISSUES,
the to-do list, the six agent prompts and `MATLAB_worker_agents.md`:
- the code is general, behind one route interface; this replaces "branch,
  do not generalise";
- a tolerance regression tier in three layers, with an anchor from the
  untouched model class, mutation checks on it, and the exact references
  re-captured once after the refactor;
- build order: pool module, tolerance tier, mutation checks, refactor,
  chain, the chain's references, chained facilitation, spread;
- 1TS carries no reserve variable (removed, not clamped); the same for
  facilitation;
- four conditions in every analysis, three in the paper; chained MTS takes
  product MTS's colour and position, product MTS is green in a fourth
  column;
- no dates in the repos' documents; quality first;
- worker prompts never open the default profile (`pool_set(N)` if present,
  else `parpool('Processes', N)`).

**Verified:** after the compaction the lead re-read the plan, the issues,
the to-do list, the facilitation note, the synapse sites of
`SRNNCellTypePairs.m` (`compile_synapse_config`, `get_params`, the layout,
`dynamics_fast`, `compute_Jacobian_fast`, `jacobian_times`, `route_fields`,
`unpack_and_compute_states`), `regression_compare` and
`regression_rapid_run`. Both trees clean; the model class is untouched.
Tried to start `matlab-opus-1`: "Agent type not found".

**Not verified / known gaps:** the worker agents do not load in this
session (it was started in the manuscript repo; an exit and resume did not
load them), so a fresh session started in FractionalReservoir is needed
before step 0a. The rewritten agent prompts are untested on a worker. Not
yet written: the route-interface design, the facilitation design section,
FD-020 and FD-021 as full issues. The tolerances (1e-12 pointwise, 1e-9
stable end to end) and the 10 % slowdown limit are the lead's proposals.
Nothing pushed (this repo 6 ahead of origin before this commit, the
manuscript repo 5).

**Raised, not acted on:** (1) The route interface must be one local system
per route: chained facilitation puts g in the b1 release term and a calcium
trace drives g, so an "std block" and an "stf block" cannot stay separate.
Whether the state vector keeps today's ordering (all b, then all g) or
becomes contiguous per route is Tom's to decide. (2) `regression_compare`
already takes `rtol` / `atol`, but compares digest statistics (sum,
weighted sum, min, max), not elements; the pointwise layer should compare
element by element and report the largest difference per leaf. (3) Carried
from earlier: `build()` reseeds the global RNG; `close all` in two tests
and apparently in `run_sensitivity_analysis`; stale pointers in
`srnn_param_preset.m:119` and three example headers; `run_all_tests('slow')`
not rerun since the regression tests joined.

**Next:** plan mode with Tom: the route-interface design into plan
section 4, FD-020 and FD-021 as full issues, the facilitation design
section; Tom's go; then a fresh session in FractionalReservoir (check
`/agents`) and step 0a (FD-017).

---

## 2026-11-11 (session 7ec5eb22-1e10-4e89-8903-1803563a876a, planning; written before compaction)

**Issues touched:** FD-016, FD-017, FD-018 created; FD-019 added to "to be
planned"; FD-007, FD-009, FD-015 revised. No code changed, so no regression
tier was run.

**Commits (this repo):** `6b3d874` Plan A written; `ab02708`, `caf6744`,
`172f38b`, `ebb0b68` Plan A revised with Tom's answers; `77ae6b0`, `e806fb2`
issues; `65b7832` the facilitation options note.
**Commits (manuscript repo):** `05c4aef`, `5f8c499`, `7221748`, `d0c4e2e`,
`4563473` (to-do list, CLAUDE.md pointer, equations file, worker-agent note).

**Done:** `docs/Plan_A_chained_vessicles.md`, the implementation plan for
block A, and Tom's decisions on it, all written into the plan and the to-do
list:
- reserve refill 5 s; chained MTS named `sfa3_stdc2`;
- spread allowed on chain, one-component and product routes; set in the
  preset (`std_spread`), condition names unchanged;
- four conditions analysed, two figure sets per run (three- and
  four-condition), chained MTS takes product MTS's colour, product MTS green;
- agents never run `medium` or larger;
- the joint sample gets an optional seed (FD-016);
- all pool handling goes through a `pool_*` module with the session's
  allocation in an environment variable (FD-017), **built first**;
- chained facilitation (option C of the note) is built after chained
  depression, for a later paper (FD-018);
- working arrangement: the lead orchestrates and merges, at most two Fable
  5.2 and three Opus 5.5 medium workers build.

**Verified:** by a subagent in the lead's MATLAB: an environment variable is
private to the MATLAB process and survives `clear functions`; pool workers
copy it once at pool start and never see a later change; a pool worker is
identifiable from the OS by `-parallelserver` and its parent process id.
Checked directly: `parallel.Settings.Pool.AutoCreate` is off and
`parallel.defaultProfile` is `Just1` (one worker).

**Not verified / known gaps:** the facilitation note's derivations were not
all checked by the lead and none of its options was simulated. `clear all` /
`clear classes` on the client were not tested against the environment
variable. The worker agents are not loaded in this session (it predates
`.claude/agents/`); Tom is restarting the session.

**Raised, not acted on:** until FD-017 exists, code that opens the default
profile misbehaves (`restart_parpool` opens one worker; three stage
functions and `wait_for_parpool` error), so `test_regression_stages` and the
full bundle need a `Processes` pool opened by name first. Tom's
licence-holding one-worker pool was not running when checked.

**Next:** after the restart, the lead re-reads the plan and the code, asks
Tom any remaining questions, and waits for go. Then Plan A step 0a (FD-017).

---

## 2026-11-10 22:18 CST (session ceb409cc-a8b9-485a-83e9-993fc3383d2f)

**Issues touched:** none. Tooling for the working arrangement, at Tom's
request: can several subagents run MATLAB at once?

**Commits:** `1be5011` six MATLAB worker agents (this repo, force-added
under the gitignored `.claude/`); `5cb7e53` usage note and CLAUDE.md pointer
(manuscript repo). The reasoning and the wrong turn are in `1be5011`.

**Done:** `.claude/agents/matlab-fable-1..3` (Fable 5.2) and
`matlab-opus-1..3` (Opus 5.5), medium effort, each with a private MATLAB
through its own inline MCP server. Usage note:
`../StochasticPlasticDynamicalSystemPaper/docs/MATLAB_worker_agents.md`.
No source, test or fixture was touched.

**Verified:** from scratch space. Ordinary subagents share the lead's one
MATLAB (same pid, workspace and working folder; calls queue). Two instances
of one worker type share a MATLAB and overwrite each other. Two differently
named workers ran concurrently in separate MATLABs. The installed
`matlab-fable-1` and `matlab-opus-1`, spawned together with `pool: 2`, ran
on the right models, each opened and closed a 2-worker pool, and left no
process behind.

**Not verified / known gaps:** all six at once; desktop display mode; a
worker resumed after finishing; a job longer than the 30-minute call limit.
A session that was running before `.claude/agents/` existed must be
restarted to see the agents. With `AutoCreate` on, a worker briefed
`pool: 0` can still acquire the default pool by reaching a `parfor`.

**Raised, not acted on:** (1) turning `parallel.Settings.Pool.AutoCreate`
off machine-wide would close that gap but changes the desktop MATLAB too;
Tom to decide. (2) This repo's CLAUDE.md does not mention the worker agents
or the pool budget; the pointer went into the manuscript repo's CLAUDE.md as
asked. (3) The local cluster keeps jobs 1 to 9 "because they contain crash
dump files" (printed on every pool start). (4) The agent files hold this
workstation's paths.

**Next:** nothing owed. Regression tiers were not run: no code they cover
changed.

---

## 2026-11-10 21:40 CST (session 7ec5eb22-1e10-4e89-8903-1803563a876a, continued)

**Issues touched:** FD-001 in progress -> done; FD-002 open -> done; FD-003
open -> done; FD-008 dropped; FD-015 created.

**Commits (this repo):** FD-001 `ba0508c`, `a356bd4`, `32bbf5e`, `155b5ef`
(by a subagent), closed in `7793b6f`. FD-002 `a24334b`. FD-003 tooling
`7fcf533`, `10879e1`, `cae664d`; references `bef29cb` (ultrarapid),
`c6b71ec` (full), `58ed2ce` (stages). Issues list `5057e16`, `1cb82de`.

**Commits (manuscript repo):** chain parameters and implementation
decisions `7897083`; strength matching settled `7d083ac`; depression spread
`4756b8c`, `e7ce1df`.

**Done:**
- Test audit: 58 tests, all passing, behind `run_all_tests`. Four tests were
  broken on the test side and fixed; five scripts that asserted nothing were
  removed; three tests added.
- Regression bundle: the manuscript network at n = 100 with every stage and
  figure (`regression_fast_config`, `regression_fast_run`).
- Regression in three tiers, each with a saved reference made before any
  change to the model class: `test_regression_rapid` (36 s),
  `test_regression_stages` (2 min), `test_regression_bundle` (about 50 min).
- Decisions recorded in the to-do list: both depression models supported;
  1TS has no reserve variable; reserve size 2 and reserve refill 10 s;
  1TS and MTS not strength-matched at steady state; a per-neuron log-normal
  spread of 0.25 on the chain's time constants and reserve size, default off.

**Verified:**
- `run_all_tests('quick')`: 53 of 53 passed in 262 s, run at the end of the
  session with the regression presets and tests in place.
  `run_all_tests('slow')`: 6 of 6 passed earlier (subagent, before the
  regression tests were added).
- Each regression tier run twice on the same code: ultrarapid bit-identical;
  stages bit-identical apart from wall-clock fields; full (two complete
  analysis runs, 31.7 and 40.5 min) bit-identical apart from wall-clock
  fields, the commit in the manifest, and the joint sample.
- The full reference against the second run: no differences.
  `test_regression_rapid` and `test_regression_stages` pass against their
  references.
- The comparison catches a 1e-12 change to one element of a million-element
  array (after a hash was added; the first version did not).
- The regression preset differs from the manuscript preset in four size
  fields only; conditions are identical.

**Not verified / known gaps:**
- `run_all_tests('slow')` has not been rerun as a whole since the
  regression tests were added to it (it now includes the 50-minute full
  tier). `test_regression_stages` was run on its own and passed;
  `test_regression_bundle` has not been run through the runner.
- The full tier's figures half has run once, not twice. Figures are compared
  by file name only, in every tier.
- No tier covers the joint parameter-space sample (see below).
- The full tier takes about 50 minutes; a smaller run mode would shorten it
  but means editing every stage function.

**Raised, not acted on:**
- **The joint parameter-space sample is not reproducible.**
  `ParamSpaceAnalysis2` calls `rng('shuffle')` before `randperm`
  (`ParamSpaceAnalysis2.m`, near line 2611), so the 128 sampled grid points
  are drawn from the clock and differ on every run, including the
  manuscript's. A seed option would make them repeatable.
- `SRNNCellTypePairs.build()` leaves the caller's global random stream
  reseeded (`generate_external_input` calls `rng(seed)` without restoring).
- `run_sensitivity_analysis` appears to call `close all` (read, not run),
  and two tests call it mid-run, so running the suite closes open figures.
- A comment at `srnn_param_preset.m:119` and three `scripts/examples`
  headers point at test files the audit removed.
- Git's background repack failed twice with "Permission denied" on a pack
  file during commits. The commits themselves succeeded.

**Next:** planning with Tom. The model class has not been touched. FD-004
(chained-pool depression) is the next piece of work, then FD-005, FD-015.

## 2026-11-10 18:31 CST (session 7ec5eb22-1e10-4e89-8903-1803563a876a)

**Issues touched:** created FD-001 to FD-007 and the "to be planned" list;
FD-001 open -> in progress.

**Commits (this repo):** `db48b4c` codebase map; `abe39aa` CLAUDE.md current
work; `8275a76`, `a0485d0`, `020bce3`, `810c31f` the four demonstration
scripts in `scripts/explorations/adaptation_normalization/`; this commit
(ISSUES.md, WORKLOG.md, CLAUDE.md).

**Commits (manuscript repo):** the to-do list, the title change, the
equation build-up in Methods, `adaptation_normalization.md`, and the
chained-pool decision (`7f5dae7`).

**Done:**
- Brian's feedback turned into `final_draft_to_do.md` (blocks A to G).
- Decision on block A: chained vesicle pools for depression, 1TS as the
  chain with no reserve. Both the product and the chain will be supported in
  the model class, product as the default.
- Read the depression code paths in `SRNNCellTypePairs` to confirm the chain
  fits the existing per-route state layout.
- Set up this tracking (ISSUES.md, WORKLOG.md).

**Verified:** the four demonstration scripts run through the MATLAB MCP and
pass `check_matlab_code`; the chain's closed-form steady state matches a
long integration to about 1e-12.

**Not verified / known gaps:**
- No model code has been changed. The chain exists only in the
  demonstration scripts.
- The state of `scripts/tests/` is unknown (FD-001).
- The reserve-size and time-constant values for the chain are choices made
  for the model, not taken from checked sources.

**Raised, not acted on:**
- The chain's two conditions are nested but not matched in strength; a
  frozen-reserve control would answer "is it simply more depression". Open
  decision in the to-do list.

**Next:** FD-001 (test audit, by a subagent), then FD-002 and FD-003 (fast
regression bundle and its reference), then stop for planning with Tom.
