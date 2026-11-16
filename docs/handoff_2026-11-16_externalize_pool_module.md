# Handoff: move the `pool_*` module into its own repo so several repos can share it

- **Date:** 2026-11-16
- **From:** the agent working with Tom in `cscs_stim_stability`
  (`C:\Users\m218089\Desktop\gitlab_repos\cscs_stim_stability`)
- **For:** the agent working in this repo (`FractionalReservoir`)
- **Status:** request, **deferred by Tom** ("sometime later"). Nothing in
  this repo was changed except adding this note. No issue entry was added
  to `docs/ISSUES.md`; give it a number when the work is taken up.

## Why

`cscs_stim_stability` now has its own MATLAB worker agents
(`cscs-matlab-opus-1..3`, `cscs-matlab-fable-1..3` in its
`.claude/agents/`), modelled on the six here: each has a private MATLAB
through its own MCP server. Both repos' workers run on the same 14-core,
64 GB workstation (R5611351) and share the 12-worker budget, but only this
repo has the pool module.

So the cscs workers open their pool with a bare `parpool('Processes', N)`.
`pool_status` then lists that session as `(unregistered)` and, by design,
leaves it out of the budget sum. A lead here that reads "0 of 12 allocated"
can hand out 12 workers while a cscs worker already holds 6. (On
2026-11-16 a cscs worker held exactly such an unregistered pool of 6.)

Tom wants the allocation visible, and preferably coordinated, across
repos. Messaging between Claude sessions is possible but is not a ledger;
the registry this module already keeps is the right place.

## What is here now (read 2026-11-16, not run)

Seven files in `src/util/`, 649 lines:

| File | Role |
|---|---|
| `pool_config.m` | Profile (`'Processes'`), budget (12), env-variable names, registry directory, this session's allocation |
| `pool_set.m` | Record the allocation, open exactly N workers, write the registry entry |
| `pool_ensure.m` | The session's pool, opened at its allocation if down; serial with a warning when none is set |
| `pool_restart.m`, `pool_release.m` | Restart at the same size; close and remove the entry |
| `pool_status.m` | Table of sessions: allocated (registry) and running (operating system), pruning dead entries |
| `pool_registry.m` | One JSON per session, plus the Windows process query and the priority setter |

Facts that matter for the move:

- **The registry is already machine-level**, not per repo:
  `%LOCALAPPDATA%\matlab_pool_registry\session_<pid>.json` (under `$HOME`
  elsewhere). Any MATLAB that calls the same code appears in every other
  session's `pool_status`. No new shared file or log is needed for
  visibility.
- **The running count already ignores the registry**: it counts
  `MATLAB.exe -parallelserver` children per client from the operating
  system, which is how unregistered pools show up at all.
- **The budget is reported, never enforced**, and unregistered sessions
  are not counted in it (`pool_status.m`, `pool_config.m`).
- **Tied to this repo:** the environment variables are named
  `SRNN_POOL_WORKERS`, `SRNN_POOL_LABEL`, `SRNN_POOL_LOW_PRIORITY`; the
  module calls `vprintf` (`src/util/vprintf.m`); comments cite FD-017 and
  this repo's entry points; the test is
  `scripts/tests/test_pool_module.m`, run through `run_all_tests`.
- **Callers here:** about 245 references in 37 files. Outside the module
  itself: `src/util/wait_for_parpool.m`, `src/analysis/run_all_analyses.m`,
  `run_lyapunov_spectrum.m`, `run_memory_capacity.m`,
  `run_transient_gain.m`, `run_numerics_verification.m`,
  `run_noise_switch_rollouts.m`, `run_noise_switch_stage.m`,
  `run_mc_gain_sweep.m`, `run_mc_gain_sweep_set.m`,
  `scripts/paper/run_all_paper_analyses.m`, `make_all_paper_figures.m`,
  the `mu7chain_*_run.m` scripts, several tests and examples, plus
  `CLAUDE.md`, `docs/ISSUES.md`, `docs/WORKLOG.md` and two plan documents.
  None of these needs to change if the function names stay the same.

## The request

1. **A new small repo** holding the seven `pool_*.m` files, their test and
   a README. Suggested name `matlab-pool`, beside the other repos in
   `C:\Users\m218089\Desktop\github_repos\`. Name and host (GitHub or
   GitLab) are Tom's to choose.
2. **Keep the function names and behaviour.** Callers in this repo then
   need no edit, and the byte-exact regression tiers are unaffected (the
   module only decides whether and how large a pool is; see
   `pool_ensure`'s note that pool size moves no bytes).
3. **Remove the seven files from `src/util/` in the same change** that puts
   the new repo on the path. Two copies on the path shadow each other
   silently, and which one wins depends on path order.
4. **Make it self-contained:** replace the `vprintf` calls with a private
   helper inside the new repo (or copy `vprintf` in under another name),
   so the new repo depends on nothing here.
5. **Each repo finds it through its own `setup_paths.m`**, with a clear
   error or warning when the folder is missing, as `cscs_stim_stability`'s
   `setup_paths` already does for its sibling repos.

## Points that need care

- **Worktrees.** Workers here run in git worktrees under
  `.claude/worktrees/agent-*/`. A path written relative to the repo root
  (`../matlab-pool`) points somewhere else from a worktree. `setup_paths`
  must resolve the shared repo in a way that works from both, for example
  from a fixed base folder as `cscs_stim_stability`'s `setup_paths` does
  (`base_github = 'C:\Users\m218089\Desktop\github_repos'`), with the
  other machines handled the way `figure_sources.json` handles them in the
  manuscript repo, or from an environment variable.
- **Mixed versions during the change.** A MATLAB still running the old copy
  and one running the new copy should keep seeing each other. Keep the
  registry directory and the JSON fields as they are; add fields, do not
  rename them.
- **Environment-variable names.** `SRNN_POOL_*` is a misnomer once other
  repos use the module. Either keep the names (no risk) or rename to
  something neutral and read the old names as a fallback for one release.
  A launcher that sets them before MATLAB starts must change with them.
- **`wait_for_parpool`** is host-gated to R5611351 and calls the module. It
  can stay here as a caller; say in the handoff's follow-up which way it
  went.
- **Other platforms.** The non-Windows paths (registry under `$HOME`,
  pruning by `kill -0`, priority as a no-op) are still untested, as the
  worker note already says.

## Improvements Tom was offered (decide before or after the move)

These came out of the discussion on 2026-11-16. None is required for the
move itself; Tom has not chosen among them.

1. **A `repo` field in the registry entry and the label**, so
   `pool_status` shows which repo a session belongs to
   (`cscs/cscs-matlab-opus-1`).
2. **Count unregistered pools against the budget.** The running count is
   already measured from the operating system; the total compared with the
   budget could be `max(allocated, running)` summed over sessions, so a
   pool opened by hand or by old code still reduces what is free. Tom's
   one-worker `Just1` licence-holding session would then need an explicit
   exemption.
3. **Enforce the budget in `pool_set`:** refuse, naming the current
   holders, when the request would take the machine over budget, with an
   override option. This needs a lock around the read-and-write (two
   sessions calling at once must not both pass), for example an exclusive
   lock file in the registry directory. Today's design note says exceeding
   the budget "slows the machine, it does not break a run"; enforcement is
   a change of policy, so it is Tom's call.
4. **A history log** (append-only, one line per `pool_set` /
   `pool_release` with time, pid, label, repo, count) written as a side
   effect, for looking back at who held what. It must not be read to
   decide anything: a worker's MATLAB dies when the agent ends, so
   check-ins are missed and a log goes stale, which the pid-pruned
   registry does not.

Memory is the other shared limit and is outside this module.

## What changes in the other repos afterwards

`cscs_stim_stability` (its agent can do these once the shared repo exists):

- `setup_paths.m`: add the shared repo.
- The six `.claude/agents/cscs-matlab-*.md` prompts: replace the bare
  `parpool('Processes', N)` / `delete(gcp('nocreate'))` instructions with
  `pool_set(N, 'label', ...)` and `pool_release()`.
- Code that opens a pool by itself: `manuscript_stim/ms_tica_prepare.m`
  (line 39), `stage3_xltek_tea/run_stage3.m` (107),
  `stage4_analysis/compute_channel_psds.m` (77),
  `stage5_ica/plot_component_psds.m` (81). The last three call `parpool`
  with no arguments, which opens the default `Just1` profile on this
  machine.

Documentation that describes the module's location:

- This repo: `CLAUDE.md`, the six `.claude/agents/matlab-*.md` prompts
  (they say "the project's pool module (`src/util/pool_*.m`)"),
  `docs/ISSUES.md` (FD-017).
- `StochasticPlasticDynamicalSystemPaper/docs/MATLAB_worker_agents.md`,
  the usage note for the worker agents.

## Done means

- The seven functions live in one repo, are gone from `src/util/`, and
  `which pool_set -all` shows a single copy in this repo's MATLAB, in a
  worktree's MATLAB and in `cscs_stim_stability`'s MATLAB.
- `test_pool_module` passes from its new home, and this repo's ultrarapid
  tier and one pool-dependent tier pass unchanged.
- With a pool open in a `FractionalReservoir` session and one in a
  `cscs_stim_stability` session, `pool_status` in either lists both with
  their allocations, and the total is the sum.
- The documents listed above name the new location, and Tom has decided
  on each of the four offered improvements (done, declined or later).
