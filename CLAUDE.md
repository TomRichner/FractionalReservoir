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

Primary entry points (current, not legacy):

- `scripts/run_param_space_analysis2.m` — multi-dimensional grid sweep across SRNNModel2 parameters
- `scripts/run_sensitivity_analysis.m` — 1D sweeps (uses `ParamSpaceAnalysis2` with `randomize_order=false` and `reps` as a grid axis)
- `scripts/run_tau_sensitivity_analysis.m` — vector-parameter sweep over `tau_a_E` / `tau_b_E_rec`
- `scripts/run_all_analyses.m` — orchestrator that runs the three above into a single dated `data/param_space/run_all_<dt>/` directory

There is no standalone test framework; ad-hoc verification scripts are named `scripts/test_*.m` (e.g. `test_SRNN2_defaults.m`, `test_psa_saveload.m`, `test_sensitivity_refactor.m`). Run them from the MATLAB editor or via the matlab MCP `run_matlab_file` tool. The matlab MCP server is available — prefer `check_matlab_code` for static analysis and `run_matlab_file` / `evaluate_matlab_code` for execution; the user's MATLAB desktop is the visible UI for any figures.

## Architecture

The codebase has converged on **two main classes** that drive everything; many older standalone files exist but are now wrappers, duplicates, or legacy.

### `src/SRNNModel2.m` (the model)

A `handle` class encapsulating the full SRNN simulation: parameter storage (with RMT-derived dependent properties like `alpha`, `R`, `n_E`, `N_sys_eqs`), W matrix construction, stimulus generation, ODE integration, Lyapunov computation, decimation, and plotting. Lifecycle: construct with name-value pairs → `build()` → `run()` → `plot()`.

State vector packing is `S = [a_E(:); a_I(:); b_E(:); b_I(:); x(:)]` of length `N_sys_eqs = n_E*n_a_E + n_I*n_a_I + n_E*n_b_E + n_I*n_b_I + n`. The classes assume `n_b_E, n_b_I ∈ {0, 1}`.

Many functions formerly in `src/` have been **internalized as static methods on `SRNNModel2`** (commit `6e2c58f`, 2026-02-27): `compute_Jacobian_fast`, `dynamics_fast`, `generate_external_input`, `compute_lyapunov_exponents_internal`, `benettin_algorithm_internal`, `lyapunov_spectrum_qr_internal`, `decimate_states`, `initialize_state`, `unpack_and_compute_states`, all `plot_*` and colormap helpers, and the activation functions. The standalone `src/` files of the same name still exist for backward compatibility but the class no longer calls them. When editing, **prefer modifying the class method**, not the standalone file. See `docs/refactors/internalize_functions_into_classes.md` for the full mapping.

`src/SRNNModel.m` is the older single-class predecessor (default `n=100`, `indegree=20`); `SRNNModel2` (default `n=300`, `indegree=100`) is the current model. Treat `SRNNModel.m` as legacy unless explicitly asked.

### `src/ParamSpaceAnalysis2.m` (the analysis driver)

A `handle` class that runs grid sweeps over `SRNNModel2`. Key behaviors not obvious from individual files:

- Configure with `add_grid_parameter(name, [min, max])`, `add_vector_parameter(...)`, `set_conditions(...)`, and `model_defaults` struct entries.
- The same network seed (W) is reused across all conditions for a given grid point so adaptation conditions can be compared fairly.
- Default conditions are 4 adaptation regimes: `no_adaptation`, `sfa_only`, `std_only`, `sfa_and_std` — toggle via `n_a_E` / `n_b_E`.
- Execution is **batched parfor with checkpoint files** (`batch_size`, `output_dir`); a run can be resumed by loading the PSA and calling `run()` again.
- `randomize_order=true` (default) gives representative early-stopping; `false` is required for sensitivity sweeps where ordered axes matter.
- Reps are added as a grid axis (`add_grid_parameter('reps', 1:n_reps)`), not as a separate property.
- Output directory defaults to `<root>/data/param_space/<folder_prefix>_<note>_<timestamp>/`. `folder_prefix` is overridden to `'1D_sensitivity'` / `'tau_sensitivity'` by the sensitivity scripts.

`src/SensitivityAnalysis.m` is the older 1D-only predecessor; new sensitivity work uses `ParamSpaceAnalysis2` with one grid parameter + `reps`.

### Master-orchestrator conventions

`run_all_analyses.m` sets two variables that downstream scripts conditionally honor:

- `master_output_dir` — when set, sub-scripts write into this shared dir instead of creating their own and **skip their `clear`/`clc`/`close all`**.
- `master_save_figs` — `'save_all_figs'` / `'save_no_figs'` / `'follow_scripts_save_figs'` overrides each sub-script's local `save_figs` flag.

Preserve this pattern when adding new analysis scripts.

### `src/` layout

- `algorithms/Lyapunov/`, `algorithms/Jacobian/`, `algorithms/info/` — numerical methods. Note that `compute_Jacobian_fast.m`, `compute_lyapunov_exponents.m`, `benettin_algorithm.m`, `lyapunov_spectrum_qr.m` are also internalized into `SRNNModel2`.
- `connectivity/` — `create_W_matrix.m`, `create_paired_W_matrix.m`, `RMTMatrix.m` (RMT-based connectivity, still external).
- `nonlinearities/` — `piecewiseSigmoid`, `logisticSigmoid`, `tanhActivation` and derivatives (also internalized).
- `generate_stimulus/` — input generators (`generate_external_input`, `generate_AM_pulse_train`, `generate_paired_pulse_input`, `generate_mackey_glass`).
- `plotting/` — colormaps, line/scatter helpers, time-series panel plots, `param_space_plots/` for post-hoc visualization, `plot_saving/save_some_figs_to_folder_2.m` (used by every script that writes figures).
- `SRNN_ESN.m`, `SRNN_ESN_reservoir.m`, `SRNN_reservoir.m`, `SRNN_reservoir_DDE.m` — alternative ODE RHS implementations and an Echo State Network variant. The DDE/ESN variants are exploratory.

### Legacy and exploratory code (likely targets for removal in a refactor)

- `scripts/old_scripts/` — predecessors of `run_*` scripts, no longer wired in.
- Single-purpose ad-hoc analyses in `scripts/`: `Sompolinsky_*.m`, `Single_vs_dual_adaptation_example.m`, `SRNN_comparisons.m`, `debug_PP.m`, `looped_memory_capacity.m`, `example_*.m`, `fraction_excitatory_analysis.m`, `Fig_2_*` (the orchestrator's stage 4 is currently commented out), `paired_pulse/`, `sine_stim/`, `review_paper/`, `python_piecewise/`, `VAR_SRNN/`.
- `reference_files/` — earlier project's source kept as a reference (`Somp_relu_centered_v5d.m`, `compute_dependent_variables.m`, `RMT.m`, …); cited by docs but not on the runtime path of the current pipeline.
- `src/SensitivityAnalysis.m` — superseded by `ParamSpaceAnalysis2`.
- `src/SRNNModel.m` — superseded by `SRNNModel2`.

When working on the current pipeline, default to `SRNNModel2` + `ParamSpaceAnalysis2`. Do not import from `reference_files/` or `old_scripts/` unless explicitly asked.

## Data conventions

- Output `.mat` files, figures (`.fig`/`.png`/`.pdf`/`.svg`) and the `.venv` are all gitignored — do not commit them.
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
