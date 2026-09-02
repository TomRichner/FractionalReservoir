# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Scope

Finish what was asked, and raise anything else rather than absorbing it into the
task. A defect you notice in the middle of another job is worth reporting in a
sentence or two; it is not a reason to go fix it, least of all in shared code
(`src/plotting/`, the model classes, `ParamSpaceAnalysis2`) or folded into a
commit about something else.

Commit messages are the record of what was done. Write them so a later session
can reconstruct the reasoning without a separate narrative file, and note what
turned out to be **wrong** as well as what worked — the wrong turns are the part
that is otherwise lost, and they are what stops the next session re-deriving a
dead end.

## What this project is

MATLAB research code for simulating and analyzing a Spiking Rate Neural Network (SRNN) reservoir with spike-frequency adaptation (SFA) and short-term synaptic depression (STD). The dynamics implemented are:

```
dx_i        = (-x_i + Σ_j w_ij · θ_j + u_i)/τ_d · dt  +  (σ_u/τ_d) · dW_i
θ_i         = r_i · Π_m b_{i,m}
r_i         = φ(x_i - (c/K) · Σ_k a_{i,k})
da_{i,k}/dt = (-a_{i,k} + r_i) / τ_{a_k}
db_{i,m}/dt = (1 - b_{i,m})/τ_rec_m - (b_{i,m} · r_i)/τ_rel_m
```

**The equations live in ONE human-written file**, `docs/EquationsParametersDocs/Equations_stability_paper.md`
(plus the *generated* `doc_equations_table/equation_table.md`, rebuilt from a live
model on every figure run). They were previously restated by hand in eleven
places and drifted apart; the block above is a summary for orientation, and that
file is authoritative. **Do not add a twelfth copy** — link instead.

Three things the block above is easy to get wrong:

**Placement of `b`.** `θ_i` is the SYNAPTIC OUTPUT — the rate after depression.
The rate `r_i` is the **pre-depression** output of the nonlinearity; depression
enters as the product `θ_j` in the recurrent sum, i.e. presynaptically and
multiplicatively. Consequently SFA and STD are both driven by the raw rate `r_i`,
not by `θ_i`. This is not cosmetic: the alternative framing `r_i = b_i·φ(...)`
would make SFA integrate `b_i r_i`, make the STD ODE depend on `b_i² r_i`, and
put a factor of `b` into the `a→x` and `a→a` Jacobian blocks. The code
(`dynamics_fast`, `compute_Jacobian_fast`) implements the form above.

**`c` is divided by `K`, the SFA timescale count** (`params.c_eff = obj.c ./ max(1, obj.n_a)`,
`SRNNCellTypePairs.m:659`). Every `a_k` relaxes to the rate, so `Σ_k a_k → K·r`
and the division makes the steady-state adaptation `c·r` exactly — independent of
how many timescales carry it. **`c` is therefore the TOTAL adaptation budget**,
and changing `K` changes the timescale *structure* without changing adaptation
*strength*. Depression needs no such factor because it enters as a PRODUCT: each
`b` rests at 1, so a second timescale squares the steady state rather than
subdividing it, which is deliberate. Note `SRNNModel2` also has `c/K` but
`src/model/jacobian/compute_J_eff.m` does **not** — see its header.

**Facilitation exists but is unused by the paper's preset.** `SRNNCellTypePairs`
supports per-route STF, contributing a further presynaptic product `Π_n g_{i,n}`
to `θ`, with `dg/dt = (1-g)/τ_dec + (G-g)·r/τ_fac`.

Connectivity uses Random Matrix Theory (Harris 2023) tilde-notation. Outputs are largest Lyapunov exponent (LLE), full Lyapunov spectrum (QR method), firing rate statistics, and parameter-sweep histograms.

## Running things

Path setup is a **once-per-session** action: run `setup_paths` with the MATLAB cwd at the project root (the function lives at `setup_paths.m` in the repo root, so it resolves with no prior `addpath`). It adds `src/` and `scripts/` recursively and is idempotent. After that, any script in the repo runs from any cwd.

Entry-point scripts (the `run_all_analyses` pipeline, the `Fig_*` presentation scripts, the runnable memory-capacity scripts) still call `setup_paths()` on their first line so they can be launched cold. Smaller scripts — everything in `scripts/tests/`, examples — assume the session is already bootstrapped and contain **no path code at all**. Do not reintroduce per-script `addpath`/`genpath` bootstrap lines.

### The two entry points

**Everything the paper needs runs from two functions in `scripts/paper/`:**

```matlab
run_dir = run_all_paper_analyses();          % all heavy compute, hours
results = make_all_paper_figures();          % every figure, minutes
```

- **`paper_config.m` is the one file to edit.** It names the preset, the run mode, and the figure registry. Change the preset there and both masters follow.
- `run_all_paper_analyses` runs the sweep pipeline (which creates the dated run directory and its `run_manifest.mat`), then memory capacity, the MC example, and the eig-heatmap sampling into subfolders of it — so one run directory holds everything the paper was built from. Each stage is wrapped, so one failure does not cost the rest of an overnight run.
- `make_all_paper_figures` resolves **one** run directory and passes it to all 17 figures, then regenerates the manuscript's equation/parameter tables. It reports success as *files actually on disk*, not as "the call returned" — `save_some_figs_to_folder_2` catches export failures and carries on, so a figure can otherwise "succeed" having written nothing.

Four figures name their own preset because they are deliberately different networks (bursting, Sompolinsky, the STF single-neuron cartoon, and memory capacity); every other figure inherits `cfg.preset_name`.

### The sweep pipeline

`src/analysis/` holds the orchestrator and its three sub-analyses. **They are all functions**, taking a context struct and returning their output directory:

- `run_all_analyses(preset_name, run_mode, ...)` — orchestrator; runs the three below into a single dated `data/param_space/run_all_<dt>/` directory and returns it. Defaults to the paper's preset at `'medium'`. This is the **sweep** pipeline only.
- `run_sensitivity_analysis(ctx)` — 1D sweeps (`ParamSpaceAnalysis2` with `randomize_order=false` and `reps` as a grid axis). Returns one dir per swept parameter.
- `run_tau_sensitivity_analysis(ctx)` — vector-parameter sweep over `tau_a_E`.
- `run_param_space_analysis(ctx)` — multi-dimensional grid sweep. (Renamed from `run_param_space_analysis2.m`; the trailing `2` came from `ParamSpaceAnalysis2` and disambiguated nothing.)
- `resolve_run_context(analysis, ...)` — builds the `ctx` all three take: preset, model class, conditions, `analysis_run_config` output, merged `model_defaults`, `integer_params`, and the per-class `f_param`.

**No script in this repo runs another script and communicates through workspace variables.** Inputs are arguments; outputs are return values. The `master_*` base-workspace protocol that used to carry settings from the orchestrator into the sub-scripts (read with `exist(...,'var')`) is gone — see "Master-orchestrator conventions" below for what replaced it and why.

`setup_paths.m` lives at the **repo root** and derives everything from its own location. It deliberately enumerates `src/` and `scripts/` rather than `genpath`-ing the root: `data/`, `figs/` and `docs/` must stay off the path, and `data/param_space/run_all_*/` holds copies of the launcher scripts that would otherwise shadow the originals. It never calls `savepath`. `run_all_analyses.m` and friends derive `project_root = fileparts(which('setup_paths'))`, so they tolerate living in a subdirectory.

There is no standalone test framework; ad-hoc verification scripts live in `scripts/tests/` as `test_*.m` (e.g. `test_SRNN2_defaults.m`, `test_psa_saveload.m`, `test_SRNNCellTypePairs.m`, `test_SRNN2_S_c_heterogeneity.m`). Run them from the MATLAB editor or via the matlab MCP `run_matlab_file` tool. Two report styles coexist: some print `PASS`/`FAIL` per check and a final banner, others `assert` and are silent until they throw.

**Run the tests that cover what you touched, including the plotting assertions.** Several tests reach into figure objects (`test_SRNNCellTypePairs` checks the actual line colours of `plot_celltypes`), so a change to a colormap or a plot helper can break a test that looks unrelated to the model maths.

**Always run MATLAB through the matlab MCP server — never `matlab -batch` or any other shell invocation.** Use `check_matlab_code` for static analysis and `run_matlab_file` / `evaluate_matlab_code` for execution. The MCP server drives the user's running MATLAB desktop, which is the visible UI for any figures; a `-batch` process would be headless, would not share that session, and is not how this project is meant to be driven. Since the session is shared and live, leave it as you found it: restore the path if a check calls `restoredefaultpath`, and don't clear the workspace beyond what the script itself does.

**Gotcha: after editing a `classdef`, run `clear classes`.** `ParamSpaceAnalysis2.srnn_property_info` caches `SRNNModel2`'s property list in a `persistent` (it is called once per result by `effective_param`, so plotting loops hit it thousands of times; uncached it costs ~0.7 ms a call, which is seconds per figure on a production-size sweep). MATLAB does **not** invalidate that cache when a classdef changes on disk, so adding or renaming an `SRNNModel2` property mid-session makes the new name report as

> `'activation' is not a property of SRNNModel2. Did you mean one of: activation_function, ...`

which reads like a bug in the calling code rather than a stale cache. `clear ParamSpaceAnalysis2` and `clear functions` do **not** clear it — only `clear classes` does. That is a deliberate exception to "don't clear the workspace" above: it wipes the user's base workspace, so say so when you use it in their live session. MATLAB needs `clear classes` after a classdef edit anyway whenever live instances exist.

## Architecture

The codebase has converged on **one analysis driver** (`ParamSpaceAnalysis2`) and **two model classes it can drive** (`SRNNModel2` and `SRNNCellTypePairs`). Legacy predecessor classes and unused standalone duplicates were removed on the `refactor` cleanup branch; what remains is the current pipeline plus a small set of example/figure scripts that use the current classes.

The two model classes are **duck-typed siblings, not a hierarchy**: they share no implementation (separate `dynamics_fast`, `compute_Jacobian_fast`, state packing, and their own copies of the nonlinearity statics). A behavioural change to one is **not** inherited by the other — check whether the change belongs in both. `ParamSpaceAnalysis2` reaches either by name, so no base class is needed.

### `src/model/SRNNModel2.m` (the model)

A `handle` class encapsulating the full SRNN simulation: parameter storage (with RMT-derived dependent properties like `alpha`, `R`, `n_E`, `N_sys_eqs`), W matrix construction, stimulus generation, ODE integration, Lyapunov computation, decimation, and plotting. Lifecycle: construct with name-value pairs → `build()` → `run()` → `plot()`.

State vector packing is `S = [a_E(:); a_I(:); b_E(:); b_I(:); x(:)]` of length `N_sys_eqs = n_E*n_a_E + n_I*n_a_I + n_E*n_b_E + n_I*n_b_I + n`. The classes assume `n_b_E, n_b_I ∈ {0, 1}`.

Beyond `plot()` it carries `plot_W` (the **scaled** W — the one actually simulated, including `level_of_chaos` — imaged with a diverging colormap so zero is white and the E/I blocks read at a glance).

Many functions formerly in `src/` have been **internalized as static methods on `SRNNModel2`** (commit `6e2c58f`, 2026-02-27): `compute_Jacobian_fast`, `dynamics_fast`, `generate_external_input`, `compute_lyapunov_exponents_internal`, `benettin_algorithm_internal`, `lyapunov_spectrum_qr_internal`, `decimate_states`, `initialize_state`, `unpack_and_compute_states`, all `plot_*` and colormap helpers, and the activation functions. The standalone `src/` copies that were no longer referenced by any kept script have since been **deleted**; a few standalone helpers remain only because example/figure scripts still call them directly (e.g. `plotting/plot_firing_rate.m`, `model/jacobian/compute_Jacobian_fast.m`). When editing model behavior, **prefer the class method**, not any remaining standalone file. See `docs/refactors/internalize_functions_into_classes.md` for the original mapping.

The older single-class predecessor `src/SRNNModel.m` (default `n=100`, `indegree=20`) was removed in the cleanup; `SRNNModel2` (default `n=300`, `indegree=100`) is the E/I model class. Its Echo State Network subclass `src/model/SRNN_ESN_reservoir.m` adds memory-capacity tooling.

### `src/model/SRNNCellTypePairs.m` (the per-cell-type model)

The generalization of `SRNNModel2` from two hardwired populations to **C named cell
types with per-route synapses**. Same lifecycle (`build()` → `run()` → `plot()`) and
the same duck-typed result contract, so `ParamSpaceAnalysis2` drives it via
`model_class = 'SRNNCellTypePairs'`.

What it buys, and what it costs:

- **Per-route STD/STF.** `synapse_config.<pre>.<post>.std` / `.stf` puts depression or
  facilitation on one route only — `synapse_config.E.E.std` is "STD on E→E
  connections and nowhere else", which `SRNNModel2` cannot express (its `n_b_E = 1`
  depresses *every* outgoing E synapse). STD fields are `tau_rec`/`tau_rel`; STF fields
  are `tau_dec`/`tau_fac`/`G`, with `dg/dt = (1−g)/tau_dec + (G−g)·r/tau_fac`, so `G`
  is the ceiling the facilitated gain approaches.
- **Required constructor arguments.** `n_cellTypes`, `cell_type_names`, `f`,
  `mu_tilde_relative`, `sigma_tilde_relative` have no defaults — the class is general
  over cell types, so there is no sensible one. Everything else defaults to
  `SRNNModel2`'s values (`n=300`, `indegree=100`, `S_c=0.40`, `tau_d=0.1`, …).
- **Per-type parameters are `1 x C` rows**: `f`, `n_a`, `c`, `mu_S_c`, `sigma_S_c`, and
  `tau_a` as a `1 x C` cell. `mu_tilde_relative` / `sigma_tilde_relative` are `C x C`
  blocks `(post, pre)`, or a `1 x C` presynaptic row broadcast down the columns.
- **Scalar sweep aliases** (Dependent *with* setters, so they count as settable) exist
  because a PSA grid axis has to be a settable property in its own right, and the
  per-type parameters are rows, blocks and cells. `mu_EE_relative`, `mu_EI_relative`, …
  and `f_E` (writing it sets `f = [f_E, 1-f_E]`) are **two-type only** and error when
  `n_cellTypes ~= 2`; `tau_a_E` aliases `tau_a{1}` for any number of types but errors
  unless `cell_type_names{1}` is `'E'`, so the name cannot lie. Add an alias rather than
  teaching the sweep scripts to special-case a shape.
- **`C = 1` works, and is for figures only.** A one-cell-type network builds, runs,
  and plots — `sompolinsky_pairs`, `single_neuron_stf` and `single_neuron_dualStd`
  are all `C = 1`, the last two at `n = 1` with `W = 0`. But the aliases above still
  refuse: `f_E` and the four `mu_*_relative` blocks are the sweep axes, and both are
  meaningless with one type (`f_E` is "the fraction excitatory"; the four blocks
  collapse into one). So the boundary is **`C = 1` for figures, `C >= 2` for sweeps**,
  and those errors are what enforce it — do not "fix" them, and do not teach the sweep
  scripts to select axes conditionally. `srnn_adaptation_conditions` takes the type
  count as its 4th argument so the four named regimes get `1 x C` `n_a` rows.
  (Until 2026-08-23 `C = 1` failed in `build_network`, which configured `RMTBlocks`
  piecemeal where `set_types` is the only supported way to change `D`. If you touch
  that code, keep it atomic — and note `scripts/tests/test_pairs_single_celltype.m`
  freezes pre-fix `W` checksums for all 14 two-type presets.)
- **`plot_data.r` is keyed by `cell_type_names`**, and `br` is called
  **`synaptic_output`** here — with STF in play the quantity is not `b·r`, so it is not
  a misnamed `br`. It nests by route (`.E.E`), which is why PSA's poolers recurse.
- **Its own plots**: `plot` (compact summary), `plot_celltypes` (one column per type,
  every neuron trace, `prod(b)` collapsed across routes), `plot_eigenvalues(times_sec)`
  (needs `store_full_state = true`), `plot_W`, `plot_W_spectrum`,
  `plot_weight_histogram`. `type_colors` names its first two colours outright — type 1
  (E) warm, type 2 (I) cool — and takes `lines()` rows 3+ for further types. Wherever
  types or routes are **overlaid in one axes**, they are drawn back-to-front via
  `draw_order` (`n:-1:1`) so type 1 lands on top; colours stay indexed by type and
  legends are built in natural order, so the reversal is layering only. Per-neuron
  accents come from `neuron_colors`, which varies *lightness* within the type's hue —
  the previous `0.5*type_color + 0.5*lines(n)` blend pulled in a second hue and turned a
  fraction of the E traces blue, which the layering cannot tolerate.
- `RMTBlocks` returns a **dense** W; this class re-sparsifies, because its state
  indexing assumes sparse connectivity.

The earlier per-type class `src/SRNNCellTypes.m` — absolute tildes, no `_relative`
port, activation handles bound by capturing `obj` — was **deleted on 2026-08-19**,
along with its tests and examples and the parity section of
`test_SRNNCellTypePairs.m` that compared the two. `SRNNCellTypePairs` is now
verified against its own finite-difference Jacobian and internal consistency, not
against an independent implementation.

### `src/analysis/ParamSpaceAnalysis2.m` (the analysis driver)

A `handle` class that runs grid sweeps over `SRNNModel2`. Key behaviors not obvious from individual files:

- Configure with `add_grid_parameter(name, [min, max])`, `add_vector_parameter(...)`, `set_conditions(...)`, and `model_defaults` struct entries.
- **`model_class`** (default `'SRNNModel2'`) selects which model class to sweep — set it to `'SRNNCellTypePairs'` to sweep per-route synapse models. The classes are **duck-typed**, not a hierarchy: PSA needs a constructor and a property list, both reachable by name via `feval` / `meta.class.fromName`, so no shared base class is required. A run records its class in `resolved_defaults`/`summary_data`, and `same_config` refuses to pool runs from different classes.
- Conditions may carry **any** fields, not just `n_a_E`/`n_b_E` — a `SRNNCellTypePairs` condition can carry a whole `synapse_config` struct, which is how "STD on E→E only" becomes a condition. `run_single_job` builds constructor args in precedence order (`model_defaults` < condition < grid) and lets last-write-win, so a swept parameter is applied *after* the defaults that give it meaning.
- `model_defaults` is validated at the top of `run()` by `validate_model_defaults()`: unknown names (typos), `Dependent` properties (`alpha`, `R`, `n_E`) and non-public ones (`W`) are a hard **error** listing every problem at once. A name a **condition sets** only **warns** (it can never take effect). A name shadowed by a **grid axis** is silent — that is normal with parameter presets, which carry a value for every axis while each sweep varies a different subset; `run()` just reports it once under `verbose`. This runs before the output directory is created, so a rejected config leaves no empty dated folder.
- The condition fields `n_a_E` / `n_a_I` / `n_b_E` / `n_b_I` may **never** be grid axes — `add_grid_parameter`, `add_vector_parameter` and `run()` all reject them. Grid parameters are applied *after* condition parameters and the constructor is last-write-wins, so a gridded `n_a_E` would silently override every condition and collapse the four adaptation regimes into one. Vary adaptation with `set_conditions()` instead.
- Only the condition fields a condition **actually sets** are treated as condition-owned. The default conditions set `n_a_E`/`n_b_E` but not `n_a_I`/`n_b_I`, so `model_defaults.n_a_I` takes effect normally. (Before this was fixed it was silently dropped and the model ran at the class default.)
- To read a parameter back off a result, use `psa.effective_param(res, name)` rather than reaching into `res.config` / `model_defaults`. It applies the same precedence `run_single_job` uses (grid → condition → `resolved_defaults` → `model_defaults` → `SRNNModel2` class default) and handles vector parameters, where `res.config` holds a level *index* rather than the value. Pass `res = []` for run-wide parameters. `ParamSpaceAnalysis2.class_default(name)` gives the bare class default. Never hand-copy a class default as a plotting fallback.
- `run()` freezes the full effective parameter set into `resolved_defaults` (saved in both `psa_object.mat` and `param_space_summary.mat`), so a run directory describes itself. `same_config` compares that — runs predating the field are refused unless called with `'allow_legacy', true`. Grid axes and condition fields are excluded from it, as are the RMT `mu_*_tilde` defaults, which depend on the grid point and are filled in `build()`.
- Histogram plots colour by `f` by default; pass `'color_by', <param>` (`'ColorBy'` in the standalone `load_and_*` plotters) to colour by whatever was actually swept.
- The same network seed (W) is reused across all conditions for a given grid point so adaptation conditions can be compared fairly.
- Default conditions are 4 adaptation regimes: `no_adaptation`, `sfa_only`, `std_only`, `sfa_and_std` — toggle via `n_a_E` / `n_b_E`.
- Execution is **batched parfor with checkpoint files** (`batch_size`, `output_dir`); a run can be resumed by loading the PSA and calling `run()` again.
- `randomize_order=true` (default) gives representative early-stopping; `false` is required for sensitivity sweeps where ordered axes matter.
- Reps are added as a grid axis (`add_grid_parameter('reps', 1:n_reps)`), not as a separate property.
- Output directory defaults to `<root>/data/param_space/<folder_prefix>_<note>_<timestamp>/`. `folder_prefix` is overridden to `'1D_sensitivity'` / `'tau_sensitivity'` by the sensitivity scripts.

The older 1D-only predecessor `src/SensitivityAnalysis.m` was removed in the cleanup; sensitivity work uses `ParamSpaceAnalysis2` with one grid parameter + `reps`.

### Parameter presets vs. run modes

Two orthogonal knobs, deliberately kept apart:

All of these live together in **`src/presets/`**, which is the point: the two
orthogonal knobs used to sit in different trees, which hid that they are a pair.

- **`src/presets/srnn_param_preset.m`** — *which network*. `[d, model_class] = srnn_param_preset(name)` returns a struct of model overrides **plus the model class they are written for** (the name is only a lookup key; the struct is what reaches the model). Assign it to `psa.model_defaults` and layer tweaks on with plain field assignment. Presets carry physics **including** swept axes like `n`/`f`/`level_of_chaos` — a grid axis overrides the preset for the sweep that varies it, while sweeps holding that axis fixed use the preset's value. They must **not** carry `n_levels`/`n_reps` (not model properties), the run_mode timing knobs, or condition fields. The second output exists because the two model classes have **disjoint** parameter vocabularies, so a preset that did not carry its class would only fail later inside `validate_model_defaults`; it defaults to `'SRNNModel2'`. `'celltype_pairs'` is the `SRNNCellTypePairs` preset (E/I with STD on E→E and I→I).
- **`src/presets/srnn_adaptation_conditions.m`** — *which adaptation regimes*. `srnn_adaptation_conditions(model_class, synapse_config, n_a_sfa, n_cell_types)` returns the four conditions (`no_adaptation`, `sfa_only`, `std_only`, `sfa_and_std`) spelled for that class: `n_a_E`/`n_b_E` counts for `SRNNModel2`, an `n_a` row plus a whole `synapse_config` struct for `SRNNCellTypePairs`. The **names** are identical across classes, which is what keeps the `condition_titles` maps inside PSA's plotters working. Use this instead of `ParamSpaceAnalysis2`'s built-in defaults, which are `SRNNModel2`-shaped and would be passed verbatim to whatever `model_class` is set. The trailing arguments all default: `synapse_config = []` means "the default E→E and I→I routes", `n_a_sfa = 3` is how many SFA timescales the SFA regimes switch on, and `n_cell_types = 2` sets the `n_a` row length (a `C = 1` preset needs 1). `srnn_param_preset` sources that count from `d.n_cellTypes`, **not** `numel(d.f)` — the `'default'` preset has no `f` field at all.
- **`src/presets/analysis_run_config.m`** — *how much compute*. `analysis_run_config(analysis, run_mode, preset_defaults)` returns `n_levels`, `n_reps`, and a `cfg.model` struct of timing settings (`ode_solver`, `fs`, `T_range`, and `lya_T_interval` only when it should be set).

  **The two knobs are no longer fully orthogonal, and this is where they meet.** Every cell names *two* integrators — deterministic (`'rk4'` for fast/fast2/medium, `'ode45'` for production) and stochastic (`'sra1'` everywhere) — and the **preset** picks between them: `sigma_u_noise > 0` selects the stochastic one. That is why the third argument exists, and why the three sub-scripts resolve `preset_defaults` *before* calling this. Selecting here is what leaves everything else intact: `merge_struct` precedence is unchanged (`cfg.model` still wins), `ode_solver` stays banned from presets, and a σ = 0 preset is bit-identical to before the mechanism existed. `cfg.sde_solver` / `cfg.is_stochastic` live at the top level of `cfg`, never inside `cfg.model` — only `cfg.model` is merged into `model_defaults`, and neither is a model property. Note a stochastic `production` run is fixed-step and so does **not** carry `ode45`'s adaptive tolerance: there is no adaptive SDE solver, because step-size control is meaningless once increments are tied to the step.

`resolve_run_context` combines them with `merge_struct(preset_defaults, cfg.model)` — **preset first**, so `run_mode` keeps final say over its own knobs, and so a whole-struct assignment cannot clobber them. That call lives in exactly one place now; it used to be copy-pasted into each of the three sub-scripts. `run_all_analyses` takes `preset_name` as its first argument and records it, with the model class, in `run_manifest`.

### Master-orchestrator conventions

`run_all_analyses` passes one **`ctx` struct** to each sub-analysis, built by
`resolve_run_context(analysis, 'preset_name', …, 'run_mode', …, 'output_dir', …, 'save_figs', …)`.
Its fields:

- `output_dir` — the shared run directory. Empty means "let `ParamSpaceAnalysis2` create its own dated folder", which is what a standalone run wants.

**The invariant that generalizes: an analysis function defaults its output into `data/`, never next to its own `.m`.** Standalone runs are supported and useful — that is how one stage gets regenerated without a 3-hour sweep — but where they *land* matters. `run_all_analyses`, `run_memory_capacity`, `run_dc_lle_analysis` and `ParamSpaceAnalysis2` all obey this. `run_eig_heatmap` and `run_memory_capacity_example` did not (they defaulted to `this_dir`, being inside figure folders), and the result was a real bug: each dropped a `.mat` beside its `.m`, the matching figure read only that path, and on 2026-08-26 two figures plotted four-day-old data while the other sixteen used the run they were handed — with every provenance file claiming the same commit. **The smell is `this_dir` in an *output* path**; `this_dir` for locating a sibling *asset* is fine. On the reading side, resolve data with `_common/resolve_data_file.m` (run directory first, then `data/`, then error naming everywhere looked) rather than a single hardcoded path.
- `save_figs` — a **logical**. (The old `master_save_figs` had a third value, `'follow_scripts_save_figs'`, meaning "let each sub-script use its own local flag". There are no local flags any more, so it had nothing left to mean.)
- `model_class` + `conditions` + `preset_defaults` — all three from **one** `srnn_param_preset` call, which is what stops a Pairs preset being swept with `SRNNModel2`-shaped conditions.
- `integer_params` — `{'n','indegree'}`, not the class default (which lists `SRNNModel2`'s adaptation counts, meaningless on Pairs).
- `f_param` — `'f'` on `SRNNModel2`, `'f_E'` on `SRNNCellTypePairs`. The fraction-excitatory axis is a scalar property on one and a scalar *alias* onto a `1 x C` row on the other. Note `plot`'s `color_by` must be given explicitly for Pairs, because its default `'f'` is a row there and breaks the histogram colouring.
- `cfg` / `model_defaults` — `analysis_run_config` output and `merge_struct(preset_defaults, cfg.model)`.

**Why this replaced the `master_*` variables.** The sub-scripts used to read
`master_output_dir` / `master_save_figs` / `master_model_overrides` /
`master_model_class` / `master_conditions` / `run_mode` / `save_figs` out of
whatever workspace called them, via `exist(...,'var')` — 37 such sites. You could
not tell what a sub-script needed without grepping; a variable left behind by one
run silently applied to the next (`run_overnight_queue.m` then existed *only* to
scrub them, as its own header admitted — it has since been rewritten as a proper
function taking a queue of `{preset, run_mode}` jobs, and is live); and the
sub-scripts had to skip their own
`clear`/`clc` when `master_output_dir` was set, leaking "am I being orchestrated?"
into their cleanup logic. All three problems are properties of shared mutable
scope, and all three vanish with arguments.

Preserve this pattern when adding new analysis functions: take a `ctx` (or plain
named arguments), return your output directory, and read nothing from the caller.

### `src/` layout

- `model/` — the three classes (`SRNNModel2`, `SRNNCellTypePairs`, `SRNN_ESN_reservoir`) and everything that defines what a network *is*:
  - `stimulus/` — `dc_staircase_stimulus.m`, `generate_paired_pulse_input.m`. These were `src/stimulus/` and `src/generate_stimulus/`, two folders holding one file each for the same job, merged 2026-09-02.
  - `integrators/` — `ode_rk4.m`, `sde_fixed_step.m`.
  - `jacobian/` — `compute_Jacobian_fast.m`, `compute_Jacobian_at_indices.m`, `compute_J_eff.m` (called directly by example scripts; `SRNNModel2` also carries internalized equivalents). The standalone `algorithms/Lyapunov/` files and `algorithms/info/mutual_info_SISO.m` were removed (internalized / unused).
  - `connectivity/` — **`RMTBlocks.m`** is the generator all three model classes now use. It extends Harris (2023) to full D×D block statistics indexed **(postsynaptic, presynaptic)**, so `mu_EE` can differ from `mu_IE`, and it supplies the bulk radius `R` and the outlier eigenvalues `lambda_O`. `RMTMatrix.m` and `RMTCellTypes.m` remain but are no longer used by the model classes. Note `RMTBlocks` returns a **dense** matrix; the CellTypes classes re-sparsify after assembly.
  - `nonlinearities/` — `piecewiseSigmoid`, `tanhActivation` and their derivatives (also internalized into `SRNNModel2`; the standalone `logisticSigmoid.m` was removed, but `logisticSigmoid`/`logisticSigmoidDerivative` live as static methods). Select one with `activation` (see "Conventions" below), not by building a handle. The **default** is `'logistic'` = `logisticSigmoid(x, S_c)` = `1/(1+exp(-4*(x-c)))` (centred at `S_c = 0.4`, range (0,1), unit slope at centre) — chosen for more robust near-edge-of-chaos stability; `'piecewise'` (`S_a = 0.9`) and `'tanh'` are the alternatives. Each model class carries its **own copy** of these statics, so a fix to one is not a fix to the others.
- `presets/` — the two orthogonal knobs, together: `srnn_param_preset` (*which network*) and `analysis_run_config` (*how much compute*), plus `srnn_adaptation_conditions`, `srnn_sfa_timescales`, `srnn_condition_titles` and `resolve_run_context`. See "Parameter presets vs. run modes" above.
- `plotting/` — colormaps, line/scatter helpers, time-series panel plots, `param_space_plots/` for post-hoc visualization, `plot_saving/save_some_figs_to_folder_2.m` (used by every script that writes figures). Unused plot/colormap duplicates were removed; the standalone `plot_*` files that example/figure scripts call directly remain.
- `model/SRNN_ESN_reservoir.m` — Echo State Network subclass of `SRNNModel2` (memory-capacity experiments). The exploratory `SRNN_ESN.m` / `SRNN_reservoir.m` / `SRNN_reservoir_DDE.m` RHS variants were removed.

### Scripts retained beyond the `run_all_analyses` pipeline

The `refactor` cleanup removed the legacy subtrees (`old_scripts/`, `review_paper/`, `VAR_SRNN/`, `python_piecewise/`, `reference_files/`) and the old-API comparison/run scripts, and reorganized `scripts/` into topic subdirectories. Current layout:

- `setup_paths.m` — shared bootstrap, **at the repo root**, not under `scripts/` (self-locating; resolvable from a cold session with cwd at the root).
- `scripts/paper/` — the two entry points and `paper_config.m`. Start here.
- `src/figures/` — the `fig_*.m` functions, **flat**. They were one folder per figure under `scripts/presentations/Stability_Manuscript/`, which made sense while each folder also held that figure's outputs; once Half A moved output to `figs/`, they were 17 directories holding one file each (and `fig_energy_landscape.m` sat four levels deep inside another figure's tree). Also here: `plot_memory_capacity`, `replot_memory_capacity`, `plot_memory_capacity_combined`.
  - `helpers/` — was `_common/`: `manuscript_style` (fonts, condition palette — returns values, never sets global state), `with_manuscript_defaults`, `build_from_preset`, `resolve_run_dir`, `resolve_data_file`, `write_manuscript_tables`, `fig_doc_tables`, `write_figure_manifest`, `figure_settings`, `default_out_dir`, `existing_outputs`, `panel_letters`, `sort_axes_left_to_right`. Joined by `save_figure_stable`, `apply_percent_axis`, `preset_default_values` and `mark_default_value`, which had been stranded in the replot folder — `save_figure_stable` alone is called by 18 figures and its own header records that it ended up there only because that is where it was extracted.
  - `replot/` — the `replot_*` regenerators, `assemble_sensitivity_figure`, `combine_runs`, `combine_two_runs`.

  **Output goes to `figs/`, not into these folders.** `cfg.fig_root` (in `paper_config`) names it: set — the default `figs/paper` — means a stable directory overwritten every run; empty means `make_all_paper_figures` auto-names `figs/figures_<dt>/` so a one-off cannot clobber the paper's set. Each registry entry gets `<fig_root>/<entry name>/`, and the root carries one `manifest.md` (run_dir, preset, run mode, per-entry results) plus `git_provenance.txt`. **`figs/` is gitignored on purpose** — a regeneration is ~200 MB of `.fig`/`.svg`, and one figure already accounts for 555 MB of this repo's history.

  Per-figure `README_*.txt` and per-figure `git_provenance.txt` are **gone** (2026-09-01) — 17 and 16 copies respectively, replaced by the one manifest. Do not reintroduce them.

  The per-entry subdirectory is load-bearing, not tidiness: `save_figure_stable` deletes `<out_dir>/<fig_tag>*` before saving, and `Fig_single_neuron_SFA_STD` is a strict *prefix* of `Fig_single_neuron_SFA_STD_STF`. Sharing a directory, the first tag's save destroys the second's outputs.

  **Every figure has the same signature**, so the master can loop over them:

  ```matlab
  function out = fig_whatever(cfg)   % cfg.run_dir/out_dir/preset_name/save/visible
  % out.figs / out.files / out.source
  ```

  **No figure calls `close all force`.** It is correct standalone — `replot_*` saves *all* open figures — but in a batch it destroys the previous entry's output before it can be verified. `make_all_paper_figures` closes each entry's returned handles once that entry has been checked, which also stops earlier figures polluting the `replot_*` prep folders.

  **Never `set(0, 'Default…')` without restoring.** Root defaults are process-global; a plotter that sets `DefaultTextInterpreter = 'none'` and does not restore it breaks the `\lambda_1` and `\mu_{EE}` labels in every figure drawn afterwards *in that session*. This actually happened. Use `with_graphics_defaults` (`src/plotting/`), which returns an `onCleanup`, or better, pass style per-object from `manuscript_style`.
- `src/analysis/` — `ParamSpaceAnalysis2` plus everything that *runs* a sweep: the orchestrator `run_all_analyses` and its three sub-analyses (all functions taking `ctx`), `run_memory_capacity`, `run_eig_heatmap`, `run_memory_capacity_example`, and the standalone tools `run_dc_lle_analysis`, `check_sensitivity_sim`, `run_overnight_queue`. `mu_block_from_preset` lives here too — it decides the mu sweep ranges and is shared by the 1-D and grid sweeps so they cannot drift apart. `analysis_run_config` and `resolve_run_context` are **not** here; they moved to `src/presets/`, being the other half of the preset pair.

  `run_eig_heatmap` and `run_memory_capacity_example` used to sit inside the figure folders whose data they produce, and defaulted to writing their `.mat` next to their own `.m` — which is exactly how two figures came to read four-day-old data (see the invariant above). They are analyses; they belong here.
**Memory capacity** was `scripts/memory_capacity/` and is now split by layer like everything else: `run_memory_capacity(cfg)` (the paired-trial ensemble, a function with its own `mc_run_config` local cost table) → `src/analysis/`; `plot_memory_capacity` and `replot_memory_capacity` → `src/figures/`; `example_memory_capacity.m` → `scripts/examples/`; `test_sample_hold_mc.m` → `scripts/tests/`. `SRNN_ESN_reservoir` subclasses `SRNNModel2`, so memory capacity **cannot** run on `SRNNCellTypePairs` — it is the one part of the paper on a different model class, behind the `mc_esn` preset. Porting the ESN readout is tracked follow-up work.
- `scripts/tests/` — verification scripts (`test_SRNN2_defaults.m`, `test_psa_saveload.m`, `test_sensitivity_refactor.m`) plus the standalone example/comparison scripts `Sompolinsky_N_1000_g_1p8.m` and `Single_vs_dual_adaptation_example.m`.
**Deleted 2026-09-02**, recoverable from git history:

- `scripts/EI_balance/` — swept `f` over 0.4–0.6 under 4 conditions on `SRNNModel2`. Superseded twice over: `run_sensitivity_analysis` sweeps `f_E` over 0.2–0.8 with reps and `run_param_space_analysis` carries it as a grid axis, both under 7 conditions on `SRNNCellTypePairs`, and `fig_EI_param_space` / `fig_EI_weights_param_space` plot the result. Note the replacement is two-part — those figures plot, the pipeline sweeps.
- `scripts/sine_stim/`, `scripts/paired_pulse/` — **broken, not merely stale**: they called `SRNN_reservoir` and `SRNN_reservoir_DDE`, removed in the earlier cleanup.
- `scripts/presentations/{OCNS, STF, Angel_ESN_example}/` — one stray PNG each; the note from the last moved to `docs/notes/ENS_stretch.md`.
- The **602 MB** of figure outputs under `Stability_Manuscript/` (124 tracked files), now that `figs/` holds them. `doc_equations_table/` went with them — it is generated into `figs/paper/doc_tables/`. The hand-drawn reservoir diagram was *not* output and moved to `docs/diagrams/reservoir_diagram/`.

**Path convention:** scripts do not assume they sit at any particular depth under `scripts/`. Derive the project root as `project_root = fileparts(which('setup_paths'))` (depth-independent) rather than with a fixed-depth `fileparts(mfilename)` chain. Scripts that need their *own* folder (to locate a sibling data file or write a figure next to themselves) still use `this_dir = fileparts(mfilename('fullpath'))` — that is a data path, not a path bootstrap, and is fine to keep.

When working on the current pipeline, default to `SRNNModel2` + `ParamSpaceAnalysis2`.

## Data conventions

- Output `.mat` files, figures (`.fig`/`.png`/`.pdf`/`.svg`) and the `.venv` are all gitignored — do not commit them. The `data/`, `figs/`, and `results_DDE/` output directories are now gitignored wholesale (added in the `refactor` cleanup), so their per-run `.m`/`.txt`/`.csv` contents stay out of git too.
- Each analysis writes `psa_object.mat` plus per-condition result `.mat` files and copies the launching script into the output dir for reproducibility (via `copyfile([mfilename('fullpath') '.m'], psa.output_dir)`).
- `data/` already contains historical run dirs (`data/param_space/`, `data/sensitivity/`, `data/sensitivity_tau_timescales_*`, `data/memory_capacity/`); leave these alone unless asked.

## Conventions worth knowing

- All scripts assume current MATLAB has the Parallel Computing Toolbox; set `psa.use_parallel = false` for serial debug runs.
- **The integrator is chosen by name**, exactly like `activation`: `ode_solver` is `'ode45'` (default) | `'ode15s'` | `'rk4'` | `'euler'` | `'heun'` | `'sra1'`. Passing a function handle raises `SRNNModel:RenamedProperty` naming the string that replaces it (`@ode_rk4` → `'rk4'`). Storing a name rather than a handle is what lets the choice sit in `resolved_defaults` as comparable data, be carried by a preset, and be validated against a known list. Both model classes route their trajectory through a protected `integrate()` helper — `SRNN_ESN_reservoir` does **not** go through `run()` (it has its own `run_reservoir_esn`), so `integrate()` is what keeps it from silently missing integrator work.
- **Additive Wiener noise on `x`** — the model can be run as an SDE: `dx_i = (…)/tau_d dt + sigma_u_noise/tau_d dW_i`. `sigma_u_noise` is **input-referred** (same units as `u`, so comparable to `intrinsic_drive` and the stimulus amplitude); the Dependent read-only `sigma_x_raw = sigma_u_noise/tau_d` is what the integrator multiplies `dW` by, and `x_noise_std = sigma_u_noise/sqrt(2·tau_d)` is the nominal stationary std of `x`. Noise enters **only** `x`, which keeps the diffusion constant — that is what makes Itô = Stratonovich, kills the Milstein term, leaves the QR variational equation untouched, and makes the noise cancel in Benettin's trajectory difference (both trajectories share one path), so the LLE stays measurable at any noise level. `sigma_u_noise > 0` **requires** `'euler'`, `'heun'` or `'sra1'`; `'rk4'` stays deterministic so σ = 0 work is bit-identical to earlier runs. `ParamSpaceAnalysis2.run()` pre-flights a swept σ against the integrator before creating the output directory.
- **Prefer `'sra1'`** (Rößler SRA1): it costs the same two drift evaluations per step as `'heun'` but is strong order 1.5 rather than 1.0 — measured ~85× more accurate at the same step size. `src/model/integrators/sde_fixed_step.m` implements all three schemes in one file so the absolute-time noise indexing (`i0 = round((tspan(1) − t0)·fs) + 1`) is written once; that indexing is what lets Benettin re-integrate a segment against exactly the increments the original run used. The noise tensor is generated at the top of `run()` from `noise_seed` (empty → `rng_seeds(1) + 224737`), kept alive through `compute_lyapunov`, then cleared — it is ~96 MB at `n=300, T=50 s, fs=400` and is never saved because it is regenerable.
- **The QR Lyapunov method always integrates its variational equation with `ode45`**, whatever `ode_solver` says. It uses a 2-point span, which the fixed-step solvers reject (`ode_rk4.m:21`), so QR previously only worked with `ode45` at all; the variational equation is deterministic and independent of how the fiducial trajectory was produced. Benettin, by contrast, deliberately uses the **same** integrator as the trajectory, so discretisation error is common to the fiducial and perturbed runs and cancels in their difference.
- **Both Lyapunov methods honour `lya_T_interval`**: they accumulate over `[max(0, lya_T_interval(1)), lya_T_interval(2)]`. (Before the `SDE` branch, `SRNNModel2`'s Benettin and both QR paths ignored `lya_T_interval(1)` and accumulated from `t = 0`.) `lya_warmup` (default 5 s) is how long they iterate *before* that window without accumulating, so the perturbation or QR basis can align first; a warmup reaching before `T_range(1)` is clamped to the available data with a warning. **Too little warmup silently biases the exponent** — on the `test_benettin_vs_qr` network, 2 s puts Benettin 5% and QR 120% off their plateau, and both need ~10 s. QR needs more warmup than Benettin because convergence tracks physical time, not renormalisation count. `lya_dt` is likewise a property now; empty means the per-method default (0.02 Benettin, 0.1 QR).
- The nonlinearity is chosen **by name**, not by passing handles: `activation` is `'logistic'` (default) | `'piecewise'` | `'tanh'`, parameterised by `S_a` (piecewise only) and `S_c` (piecewise and logistic; `'tanh'` uses neither). `activation_function` / `activation_function_derivative` are **Dependent, read-only** — read them freely (the QR Lyapunov method needs the derivative), but set `activation`. Assigning a handle raises `SRNNModel:RenamedProperty`. For a nonlinearity outside the three, set `activation_custom = {fn, dfn}`, which overrides the name. Keeping this as data is what stops `S_a`/`S_c` from silently disagreeing with the function actually in use, and lets a preset express the nonlinearity without handles.
- **The setpoint `S_c` can be per-neuron.** Set a mean and/or a standard deviation per population and `build()` draws `S_c_i = mu + sigma·randn` into the read-only `S_c_vec` (`n x 1`); leave them alone and every neuron shares the scalar `S_c`, which is the bit-identical default. The knobs follow each class's own convention: `mu_S_c_E` / `sigma_S_c_E` / `mu_S_c_I` / `sigma_S_c_I` (scalars, hence usable as PSA grid axes) on `SRNNModel2`, and `1 x C` `mu_S_c` / `sigma_S_c` rows on `SRNNCellTypePairs`. An empty `mu` falls back to `S_c`, so `sigma` alone means "spread around the shared centre". `S_c_seed` pins the draw; left empty it is derived from `rng_seeds(1)` (offset well away from the stream that builds `W`), so a reps sweep varies the setpoints too. The draw saves and restores the RNG state, so `W`, the stimulus and `x0` are unaffected either way. Two consequences: `'tanh'` and `activation_custom` have no centre to vary and **error** rather than silently ignoring the request, and in heterogeneous mode the `activation_function` handle is **only valid for length-`n` input** (`c` lines up with `x` elementwise) — evaluate φ on a plotting grid only for a homogeneous model.
- The `c_E` / `c_I` adaptation scaling is conventionally `0.15/3 ≈ 0.05` in current scripts (one-third of a "total" adaptation budget split over three timescales).
- `level_of_chaos` is a multiplicative scale on the W matrix; values >1 push the network past the edge of chaos. The `R` dependent property reports the theoretical spectral radius.
- The five RMT connectivity parameters are set as **multipliers of** `F = default_val = 1/sqrt(n·α(2−α))`: `mu_E_tilde_relative` (default 3), `mu_I_tilde_relative` (−4), `sigma_E_tilde_relative` (1), `sigma_I_tilde_relative` (1), `E_W_relative` (0). The absolute `mu_E_tilde` / `sigma_E_tilde` / `E_W` are **Dependent, read-only** — read them freely, but set the `_relative` ones. (Setting an absolute name raises `SRNNModel:RenamedProperty` naming the replacement.) Storing the multiplier is what lets these be swept — `F` depends on `n` and `indegree`, so no constant absolute value is right at more than one grid point — and frozen into `resolved_defaults`. Note `_relative` ≠ `_rel`: `tau_b_E_rel` / `tau_b_I_rel` are the STD **release** constants.
- `F_tracks_network` (default `true`) computes `F` from the current `n`/`indegree`, which makes `R` *exactly* independent of `n` (the `n·α` cancels in `get.R`). Set it `false` to pin `F` to `(F_ref_n, F_ref_indegree)` — the weight distribution then stays fixed while `R` varies with network size. Freezing `F` does **not** freeze the network: `build()` still passes the real `alpha` to `RMTBlocks`, so connectivity tracks the grid point. The choice is recorded in `resolved_defaults`, and `same_config` refuses to pool runs that used different conventions.
- Time vectors typically start negative (e.g. `T_range = [-10, 20]`) to allow transient settling before the analysis window.

## Commit messages

Never add a `Co-Authored-By: Claude ...` trailer (or any other Claude/AI attribution) to commit messages in this repo. Plain commit messages only.
