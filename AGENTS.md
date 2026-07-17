# Repository Guidelines

## Project Structure & Module Organization

This repository contains MATLAB research code for simulating and analyzing an adaptive spiking-rate reservoir. Core implementation lives in `src/`: use `SRNNModel2.m` for network dynamics and `ParamSpaceAnalysis2.m` for parameter sweeps. Supporting algorithms, connectivity, stimuli, nonlinearities, and plotting utilities are grouped in matching subdirectories. Runnable analyses and figure-generation workflows live in `scripts/`; ad-hoc verification scripts are in `scripts/tests/`. Mathematical notes, architecture descriptions, and refactor records belong in `docs/`. Generated runs go under `data/` and figures under `figs/`; both are ignored by Git.

Prefer the current class methods over similarly named standalone helpers. Treat `scripts/sine_stim/` and `scripts/paired_pulse/` as legacy references until they are ported to `SRNNModel2`.

## Build, Test, and Development Commands

There is no compilation step. Run commands from the repository root:

```sh
matlab -batch "addpath('scripts'); setup_paths; run('scripts/tests/test_SRNN2_defaults.m')"
matlab -batch "addpath('scripts'); setup_paths; run_mode='fast'; run('scripts/run_all_analyses/run_all_analyses.m')"
matlab -batch "checkcode('src/SRNNModel2.m')"
```

The first runs a model smoke test, the second exercises the analysis pipeline with reduced settings, and the third performs MATLAB static analysis. `setup_paths()` recursively adds `src/` and `scripts/`; call it at the start of new entry-point scripts. Production sweeps can be expensive and require the Parallel Computing Toolbox.

## Coding Style & Naming Conventions

Follow existing MATLAB style: four-space indentation, one primary class or function per file, `%%` section headers, and concise `%` comments explaining scientific assumptions. Use `UpperCamelCase` for classes, `snake_case` for functions and scripts, and descriptive variable names such as `master_output_dir`. Construct paths with `fullfile` and derive the project root through `which('setup_paths')`, not fixed directory-depth assumptions.

## Testing Guidelines

The project has no standalone test framework or coverage threshold. Add focused scripts named `scripts/tests/test_<behavior>.m`. Keep smoke tests small, deterministic through explicit seeds, and serializable with `use_parallel = false` when practical. Run affected tests plus a fast pipeline pass; inspect generated figures when plotting changes.

## Commit & Pull Request Guidelines

Recent commits use short imperative subjects, often with a scope, for example `SRNNModel2: make centered logistic sigmoid the default activation`. Keep each commit focused and omit AI attribution trailers. Pull requests should explain the scientific or behavioral change, list validation commands, identify parameter/default changes, and include before/after figures for visualization work. Do not commit generated `.mat`, `.fig`, `.png`, `.pdf`, `.svg`, or files under output directories unless the artifact is intentionally force-added for a documented presentation workflow.
