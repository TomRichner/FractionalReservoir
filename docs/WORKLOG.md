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
