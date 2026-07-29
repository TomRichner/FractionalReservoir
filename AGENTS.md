# Repository Guidelines

## Project Structure & Module Organization

This repository contains MATLAB research code for simulating and analyzing an adaptive spiking-rate reservoir. Core implementation lives in `src/`: use `SRNNModel2.m` for network dynamics and `ParamSpaceAnalysis2.m` for parameter sweeps. Supporting algorithms, connectivity, stimuli, nonlinearities, and plotting utilities are grouped in matching subdirectories. Runnable analyses and figure-generation workflows live in `scripts/`; ad-hoc verification scripts are in `scripts/tests/`. Mathematical notes, architecture descriptions, and refactor records belong in `docs/`. Generated runs go under `data/` and figures under `figs/`; both are ignored by Git.

Prefer the current class methods over similarly named standalone helpers. Treat `scripts/sine_stim/` and `scripts/paired_pulse/` as legacy references until they are ported to `SRNNModel2`.

## Build, Test, and Development Commands

There is no compilation step. **Never invoke MATLAB from the shell — `matlab -batch` and any other CLI invocation are off-limits.** All MATLAB work goes through the matlab MCP server, which drives the user's running MATLAB desktop (figures appear there, which is the point):

- `check_matlab_code` — static analysis of a `.m` file (replaces `checkcode`)
- `evaluate_matlab_code` — inline commands, e.g. bootstrapping the session
- `run_matlab_file` / `run_matlab_test_file` — run a script or test; the working folder is set to the script's own location

Typical session:

```matlab
setup_paths                                    % once, cwd at the repository root
run('scripts/tests/test_SRNN2_defaults.m')     % model smoke test
run_mode = 'fast'; run('scripts/run_all_analyses/run_all_analyses.m')   % pipeline, reduced settings
```

`setup_paths.m` lives at the repository root, so it resolves from a cold session once the cwd is the root; it recursively adds `src/` and `scripts/` and is idempotent. Call it **once per session**, plus at the start of new entry-point scripts so they can be launched cold. Small test and example scripts assume the session is already bootstrapped and carry no path code. Production sweeps can be expensive and require the Parallel Computing Toolbox.

Because the MCP session is the user's live desktop, treat its state as shared: don't `clear`/`close all` beyond what a script already does, and if a check needs a pristine path (`restoredefaultpath`), restore the previous path afterwards.

## Coding Style & Naming Conventions

Follow existing MATLAB style: four-space indentation, one primary class or function per file, `%%` section headers, and concise `%` comments explaining scientific assumptions. Use `UpperCamelCase` for classes, `snake_case` for functions and scripts, and descriptive variable names such as `master_output_dir`. Construct paths with `fullfile` and derive the project root as `fileparts(which('setup_paths'))`, not through fixed directory-depth assumptions.

## Testing Guidelines

The project has no standalone test framework or coverage threshold. Add focused scripts named `scripts/tests/test_<behavior>.m`. Keep smoke tests small, deterministic through explicit seeds, and serializable with `use_parallel = false` when practical. Run affected tests plus a fast pipeline pass; inspect generated figures when plotting changes.

## Commit & Pull Request Guidelines

Recent commits use short imperative subjects, often with a scope, for example `SRNNModel2: make centered logistic sigmoid the default activation`. Keep each commit focused and omit AI attribution trailers. Pull requests should explain the scientific or behavioral change, list validation commands, identify parameter/default changes, and include before/after figures for visualization work. Do not commit generated `.mat`, `.fig`, `.png`, `.pdf`, `.svg`, or files under output directories unless the artifact is intentionally force-added for a documented presentation workflow.
