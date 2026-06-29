# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

MATLAB research code for simulating and analyzing a Spiking Rate Neural Network (SRNN) reservoir with spike-frequency adaptation (SFA) and short-term synaptic depression (STD). The dynamics implemented are:

```
dx_i/dt    = (-x_i + Σ_j w_ij r_j + u_i) / τ_d
r_i        = b_i · φ(x_i - c · Σ_k a_{i,k})
da_{i,k}/dt = (-a_{i,k} + r_i) / τ_k
db_i/dt    = (1 - b_i)/τ_rec - (b_i · r_i)/τ_rel
```

Connectivity uses Random Matrix Theory (Harris 2023) tilde-notation. Outputs are largest Lyapunov exponent (LLE), full Lyapunov spectrum (QR method), firing rate statistics, and parameter-sweep histograms.

## Running things

Every entry-point script begins with `setup_paths()` (defined in `scripts/setup_paths.m`), which adds `src/` and `scripts/` recursively to the MATLAB path. Scripts in `scripts/` are intended to be run from the MATLAB cwd at the project root or from the `scripts/` directory.

Primary entry points (current, not legacy) — the orchestrator and its three sub-analyses live together in `scripts/run_all_analyses/`:

- `scripts/run_all_analyses/run_param_space_analysis2.m` — multi-dimensional grid sweep across SRNNModel2 parameters
- `scripts/run_all_analyses/run_sensitivity_analysis.m` — 1D sweeps (uses `ParamSpaceAnalysis2` with `randomize_order=false` and `reps` as a grid axis)
- `scripts/run_all_analyses/run_tau_sensitivity_analysis.m` — vector-parameter sweep over `tau_a_E` / `tau_b_E_rec`
- `scripts/run_all_analyses/run_all_analyses.m` — orchestrator that runs the three above into a single dated `data/param_space/run_all_<dt>/` directory

`scripts/setup_paths.m` stays at the `scripts/` root (every entry point depends on it, and it derives the project root from its own location). `run_all_analyses.m` derives `project_root` from `which('setup_paths')`, so it tolerates living in a subdirectory.

There is no standalone test framework; ad-hoc verification scripts are named `scripts/test_*.m` (e.g. `test_SRNN2_defaults.m`, `test_psa_saveload.m`, `test_sensitivity_refactor.m`). Run them from the MATLAB editor or via the matlab MCP `run_matlab_file` tool. The matlab MCP server is available — prefer `check_matlab_code` for static analysis and `run_matlab_file` / `evaluate_matlab_code` for execution; the user's MATLAB desktop is the visible UI for any figures.

## Architecture

The codebase has converged on **two main classes** that drive everything. Legacy predecessor classes and unused standalone duplicates were removed on the `refactor` cleanup branch; what remains is the current pipeline plus a small set of example/figure scripts that use the current classes.

### `src/SRNNModel2.m` (the model)

A `handle` class encapsulating the full SRNN simulation: parameter storage (with RMT-derived dependent properties like `alpha`, `R`, `n_E`, `N_sys_eqs`), W matrix construction, stimulus generation, ODE integration, Lyapunov computation, decimation, and plotting. Lifecycle: construct with name-value pairs → `build()` → `run()` → `plot()`.

State vector packing is `S = [a_E(:); a_I(:); b_E(:); b_I(:); x(:)]` of length `N_sys_eqs = n_E*n_a_E + n_I*n_a_I + n_E*n_b_E + n_I*n_b_I + n`. The classes assume `n_b_E, n_b_I ∈ {0, 1}`.

Many functions formerly in `src/` have been **internalized as static methods on `SRNNModel2`** (commit `6e2c58f`, 2026-02-27): `compute_Jacobian_fast`, `dynamics_fast`, `generate_external_input`, `compute_lyapunov_exponents_internal`, `benettin_algorithm_internal`, `lyapunov_spectrum_qr_internal`, `decimate_states`, `initialize_state`, `unpack_and_compute_states`, all `plot_*` and colormap helpers, and the activation functions. The standalone `src/` copies that were no longer referenced by any kept script have since been **deleted**; a few standalone helpers remain only because example/figure scripts still call them directly (e.g. `plotting/plot_firing_rate.m`, `algorithms/Jacobian/compute_Jacobian_fast.m`). When editing model behavior, **prefer the class method**, not any remaining standalone file. See `docs/refactors/internalize_functions_into_classes.md` for the original mapping.

The older single-class predecessor `src/SRNNModel.m` (default `n=100`, `indegree=20`) was removed in the cleanup; `SRNNModel2` (default `n=300`, `indegree=100`) is the only model class. Its Echo State Network subclass `src/SRNN_ESN_reservoir.m` adds memory-capacity tooling.

### `src/ParamSpaceAnalysis2.m` (the analysis driver)

A `handle` class that runs grid sweeps over `SRNNModel2`. Key behaviors not obvious from individual files:

- Configure with `add_grid_parameter(name, [min, max])`, `add_vector_parameter(...)`, `set_conditions(...)`, and `model_defaults` struct entries.
- The same network seed (W) is reused across all conditions for a given grid point so adaptation conditions can be compared fairly.
- Default conditions are 4 adaptation regimes: `no_adaptation`, `sfa_only`, `std_only`, `sfa_and_std` — toggle via `n_a_E` / `n_b_E`.
- Execution is **batched parfor with checkpoint files** (`batch_size`, `output_dir`); a run can be resumed by loading the PSA and calling `run()` again.
- `randomize_order=true` (default) gives representative early-stopping; `false` is required for sensitivity sweeps where ordered axes matter.
- Reps are added as a grid axis (`add_grid_parameter('reps', 1:n_reps)`), not as a separate property.
- Output directory defaults to `<root>/data/param_space/<folder_prefix>_<note>_<timestamp>/`. `folder_prefix` is overridden to `'1D_sensitivity'` / `'tau_sensitivity'` by the sensitivity scripts.

The older 1D-only predecessor `src/SensitivityAnalysis.m` was removed in the cleanup; sensitivity work uses `ParamSpaceAnalysis2` with one grid parameter + `reps`.

### Master-orchestrator conventions

`run_all_analyses.m` sets two variables that downstream scripts conditionally honor:

- `master_output_dir` — when set, sub-scripts write into this shared dir instead of creating their own and **skip their `clear`/`clc`/`close all`**.
- `master_save_figs` — `'save_all_figs'` / `'save_no_figs'` / `'follow_scripts_save_figs'` overrides each sub-script's local `save_figs` flag.

Preserve this pattern when adding new analysis scripts.

### `src/` layout

- `algorithms/Jacobian/` — `compute_Jacobian_fast.m`, `compute_Jacobian_at_indices.m`, `compute_J_eff.m` (called directly by example scripts; `SRNNModel2` also carries internalized equivalents). The standalone `algorithms/Lyapunov/` files and `algorithms/info/mutual_info_SISO.m` were removed (internalized / unused).
- `connectivity/` — `RMTMatrix.m` (RMT-based connectivity). The legacy `create_W_matrix.m` / `create_paired_W_matrix.m` were removed.
- `nonlinearities/` — `piecewiseSigmoid`, `tanhActivation` and their derivatives (also internalized into `SRNNModel2`). `logisticSigmoid` was removed.
- `plotting/` — colormaps, line/scatter helpers, time-series panel plots, `param_space_plots/` for post-hoc visualization, `plot_saving/save_some_figs_to_folder_2.m` (used by every script that writes figures). Unused plot/colormap duplicates were removed; the standalone `plot_*` files that example/figure scripts call directly remain.
- `SRNN_ESN_reservoir.m` — Echo State Network subclass of `SRNNModel2` (memory-capacity experiments). The exploratory `SRNN_ESN.m` / `SRNN_reservoir.m` / `SRNN_reservoir_DDE.m` RHS variants were removed.

### Scripts retained beyond the `run_all_analyses` pipeline

The `refactor` cleanup removed the legacy subtrees (`old_scripts/`, `review_paper/`, `VAR_SRNN/`, `python_piecewise/`, `reference_files/`) and the old-API comparison/run scripts. What is kept, besides the four primary entry points:

- Example/figure scripts that use the current classes: `Sompolinsky_N_1000_g_1p8.m`, `Single_vs_dual_adaptation_example.m`, `fraction_excitatory_analysis.m`, `Fig_2_fraction_excitatory_*.m`, `example_memory_capacity.m`, `looped_memory_capacity.m`, `stim_engages_adaptation/bursting_SRNN_model*.m`, and the `replot_*` figure regenerators.
- `scripts/sine_stim/` and `scripts/paired_pulse/` — kept as references but **currently non-functional**: they still use the old script-based API, and their legacy dependencies were deleted. They must be ported to `SRNNModel2` before use.

When working on the current pipeline, default to `SRNNModel2` + `ParamSpaceAnalysis2`.

## Data conventions

- Output `.mat` files, figures (`.fig`/`.png`/`.pdf`/`.svg`) and the `.venv` are all gitignored — do not commit them. The `data/`, `figs/`, and `results_DDE/` output directories are now gitignored wholesale (added in the `refactor` cleanup), so their per-run `.m`/`.txt`/`.csv` contents stay out of git too.
- Each analysis writes `psa_object.mat` plus per-condition result `.mat` files and copies the launching script into the output dir for reproducibility (via `copyfile([mfilename('fullpath') '.m'], psa.output_dir)`).
- `data/` already contains historical run dirs (`data/param_space/`, `data/sensitivity/`, `data/sensitivity_tau_timescales_*`, `data/memory_capacity/`); leave these alone unless asked.

## Conventions worth knowing

- All scripts assume current MATLAB has the Parallel Computing Toolbox; set `psa.use_parallel = false` for serial debug runs.
- Activation functions are passed as function handles in `model_defaults.activation_function` / `activation_function_derivative`; both are required for the QR Lyapunov method (which needs the Jacobian).
- The `c_E` / `c_I` adaptation scaling is conventionally `0.15/3 ≈ 0.05` in current scripts (one-third of a "total" adaptation budget split over three timescales).
- `level_of_chaos` is a multiplicative scale on the W matrix; values >1 push the network past the edge of chaos. The `R` dependent property reports the theoretical spectral radius.
- Time vectors typically start negative (e.g. `T_range = [-10, 20]`) to allow transient settling before the analysis window.

## Commit messages

Never add a `Co-Authored-By: Claude ...` trailer (or any other Claude/AI attribution) to commit messages in this repo. Plain commit messages only.
