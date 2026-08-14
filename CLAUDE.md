# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Session bookkeeping — read and update these

Three tracked documents at the repo root carry context across sessions and
across the two machines this work happens on. **Read them at the start of a
session** (they say what is already known, already tried, and already decided)
and **update them as you go** — not only at the end, since a session can be
interrupted.

- **`WorkLog.md`** — chronological narrative: plans executed, bugs found, bugs
  fixed, runs performed. Add an entry when you finish a meaningful unit of work
  (a plan stage, a bug diagnosis, an analysis run), not once per commit.
- **`Issues.md`** — bugs and observed problems, indexed by issue. Add one when
  you *find* a problem, even if you fix it immediately, and update its status
  when it is resolved or deliberately declined. Never delete an entry.
- **`FeatureRequests.md`** — wanted changes that are not defects, same
  lifecycle.

Every entry is signed with **date · branch @ commit · hostname · agent/model,
session id**, for example:

```
### 2026-08-14 07:06 · dev @ d58b7fe · R5456622 · Claude Code (Opus 5), session 054a29ca
```

Get the pieces with `git rev-parse --short HEAD`, `git branch --show-current`,
and `hostname` (identical on macOS, Windows and Linux). The session id is the
UUID in your scratchpad path. `src/capture_git_provenance.m` stamps the same
hostname onto every run directory, so a `WorkLog` entry and a `data/` folder can
be matched up.

**Record what turned out to be wrong, not just what worked.** The wrong turns
are the part that is otherwise lost, and they are what stops the next session
re-deriving a dead end — see ISSUE-009, where two plausible diagnoses were both
refuted by the data.

### Reading them without flooding your context

These files grow monotonically and are never pruned, so reading all three in
full will eventually cost more context than the work itself. Escalate only as
far as you need:

**1. The index is one grep.** Every issue and request is a `## ` heading
carrying its status marker, ID and title, so the whole backlog is a dozen-odd
lines:

```bash
grep '^## ' Issues.md FeatureRequests.md          # full index, status included
grep '^## .* OPEN ·' Issues.md FeatureRequests.md # just what is outstanding
grep -n '^### ' WorkLog.md | tail -5              # the last few entries
```

Filter on the plain-text token (`OPEN` / `FIXED` / `WONTFIX` / `DONE` /
`DECLINED`), **not** on the emoji — `grep '^## 🔴'` silently matches nothing in
some terminals depending on locale, which is a bad failure mode for a check
whose whole job is to tell you what is outstanding.

That is usually enough to know whether the thing you are about to do is already
known, already tried, or already decided against. Read a single entry in full
only when one looks relevant.

**2. Delegate the synthesis when you need the narrative.** For a genuinely
prioritised picture — what is blocked and on what, which open items touch the
area you are about to change, what was tried and refuted — use an **Explore
agent** rather than reading the files into your own context:

> Read `WorkLog.md`, `Issues.md` and `FeatureRequests.md` at the repo root.
> Return a compact prioritised summary: open issues and requests by status, what
> is currently blocked and on what, the last few work-log entries, and anything
> bearing on `<the task at hand>`. Quote IDs and headline facts, not whole
> entries. Do not modify anything.

Give it the current task so it prioritises for relevance rather than recency.
The point is that the agent absorbs the bulk and returns the conclusion.

**3. Read in full only when the narrative itself is the subject** — e.g.
reconstructing how a decision was reached, or auditing whether the record is
accurate.

The same escalation applies to writing: appending an entry needs no prior read
of the whole file.

`WorkLog.md` is a single chronological file by deliberate choice: when two
machines both append, resolve the conflict by **keeping both entries and
ordering by timestamp**, never by choosing between them. Entries are
self-contained so that is always safe.

## What this project is

MATLAB research code for simulating and analyzing a Spiking Rate Neural Network (SRNN) reservoir with spike-frequency adaptation (SFA) and short-term synaptic depression (STD). The dynamics implemented are:

```
dx_i/dt    = (-x_i + Σ_j w_ij · b_j r_j + u_i) / τ_d
r_i        = φ(x_i - c · Σ_k a_{i,k})
da_{i,k}/dt = (-a_{i,k} + r_i) / τ_k
db_i/dt    = (1 - b_i)/τ_rec - (b_i · r_i)/τ_rel
```

Note the placement of `b`. The rate `r_i` is the **pre-depression** output of the
nonlinearity; depression enters as the product `b_j r_j` in the recurrent sum, i.e.
presynaptically and multiplicatively. Consequently SFA and STD are both driven by
the raw rate `r_i`, not by `b_i r_i`. This is not cosmetic: the alternative framing
`r_i = b_i·φ(...)` would make SFA integrate `b_i r_i`, make the STD ODE depend on
`b_i² r_i`, and put a factor of `b` into the `a→x` and `a→a` Jacobian blocks. The
code (`dynamics_fast`, `compute_Jacobian_fast`) implements the form above.

Connectivity uses Random Matrix Theory (Harris 2023) tilde-notation. Outputs are largest Lyapunov exponent (LLE), full Lyapunov spectrum (QR method), firing rate statistics, and parameter-sweep histograms.

## Running things

Path setup is a **once-per-session** action: run `setup_paths` with the MATLAB cwd at the project root (the function lives at `setup_paths.m` in the repo root, so it resolves with no prior `addpath`). It adds `src/` and `scripts/` recursively and is idempotent. After that, any script in the repo runs from any cwd.

Entry-point scripts (the `run_all_analyses` pipeline, the `Fig_*` presentation scripts, the runnable memory-capacity scripts) still call `setup_paths()` on their first line so they can be launched cold. Smaller scripts — everything in `scripts/tests/`, examples — assume the session is already bootstrapped and contain **no path code at all**. Do not reintroduce per-script `addpath`/`genpath` bootstrap lines.

Primary entry points (current, not legacy) — the orchestrator and its three sub-analyses live together in `scripts/run_all_analyses/`:

- `scripts/run_all_analyses/run_param_space_analysis2.m` — multi-dimensional grid sweep across SRNNModel2 parameters
- `scripts/run_all_analyses/run_sensitivity_analysis.m` — 1D sweeps (uses `ParamSpaceAnalysis2` with `randomize_order=false` and `reps` as a grid axis)
- `scripts/run_all_analyses/run_tau_sensitivity_analysis.m` — vector-parameter sweep over `tau_a_E` / `tau_b_E_rec`
- `scripts/run_all_analyses/run_all_analyses.m` — orchestrator that runs the three above into a single dated `data/param_space/run_all_<dt>/` directory

`setup_paths.m` lives at the **repo root** and derives everything from its own location. It deliberately enumerates `src/` and `scripts/` rather than `genpath`-ing the root: `data/`, `figs/` and `docs/` must stay off the path, and `data/param_space/run_all_*/` holds copies of the launcher scripts that would otherwise shadow the originals. It never calls `savepath`. `run_all_analyses.m` and friends derive `project_root = fileparts(which('setup_paths'))`, so they tolerate living in a subdirectory.

There is no standalone test framework; ad-hoc verification scripts live in `scripts/tests/` as `test_*.m` (e.g. `test_SRNN2_defaults.m`, `test_psa_saveload.m`, `test_SRNNCellTypePairs.m`, `test_SRNN2_S_c_heterogeneity.m`). Run them from the MATLAB editor or via the matlab MCP `run_matlab_file` tool. Two report styles coexist: some print `PASS`/`FAIL` per check and a final banner, others `assert` and are silent until they throw.

**Run the tests that cover what you touched, including the plotting assertions.** Several tests reach into figure objects (`test_SRNNCellTypePairs` checks the actual line colours of `plot_celltypes`), so a change to a colormap or a plot helper can break a test that looks unrelated to the model maths.

**Always run MATLAB through the matlab MCP server — never `matlab -batch` or any other shell invocation.** Use `check_matlab_code` for static analysis and `run_matlab_file` / `evaluate_matlab_code` for execution. The MCP server drives the user's running MATLAB desktop, which is the visible UI for any figures; a `-batch` process would be headless, would not share that session, and is not how this project is meant to be driven. Since the session is shared and live, leave it as you found it: restore the path if a check calls `restoredefaultpath`, and don't clear the workspace beyond what the script itself does.

**Gotcha: after editing a `classdef`, run `clear classes`.** `ParamSpaceAnalysis2.srnn_property_info` caches `SRNNModel2`'s property list in a `persistent` (it is called once per result by `effective_param`, so plotting loops hit it thousands of times; uncached it costs ~0.7 ms a call, which is seconds per figure on a production-size sweep). MATLAB does **not** invalidate that cache when a classdef changes on disk, so adding or renaming an `SRNNModel2` property mid-session makes the new name report as

> `'activation' is not a property of SRNNModel2. Did you mean one of: activation_function, ...`

which reads like a bug in the calling code rather than a stale cache. `clear ParamSpaceAnalysis2` and `clear functions` do **not** clear it — only `clear classes` does. That is a deliberate exception to "don't clear the workspace" above: it wipes the user's base workspace, so say so when you use it in their live session. MATLAB needs `clear classes` after a classdef edit anyway whenever live instances exist.

## Architecture

The codebase has converged on **one analysis driver** (`ParamSpaceAnalysis2`) and **two model classes it can drive** (`SRNNModel2` and `SRNNCellTypePairs`), plus `SRNNCellTypes`, which survives only as the parity reference. Legacy predecessor classes and unused standalone duplicates were removed on the `refactor` cleanup branch; what remains is the current pipeline plus a small set of example/figure scripts that use the current classes.

The two model classes are **duck-typed siblings, not a hierarchy**: they share no implementation (separate `dynamics_fast`, `compute_Jacobian_fast`, state packing, and their own copies of the nonlinearity statics). A behavioural change to one is **not** inherited by the other — check whether the change belongs in both. `ParamSpaceAnalysis2` reaches either by name, so no base class is needed.

### `src/SRNNModel2.m` (the model)

A `handle` class encapsulating the full SRNN simulation: parameter storage (with RMT-derived dependent properties like `alpha`, `R`, `n_E`, `N_sys_eqs`), W matrix construction, stimulus generation, ODE integration, Lyapunov computation, decimation, and plotting. Lifecycle: construct with name-value pairs → `build()` → `run()` → `plot()`.

State vector packing is `S = [a_E(:); a_I(:); b_E(:); b_I(:); x(:)]` of length `N_sys_eqs = n_E*n_a_E + n_I*n_a_I + n_E*n_b_E + n_I*n_b_I + n`. The classes assume `n_b_E, n_b_I ∈ {0, 1}`.

Beyond `plot()` it carries `plot_W` (the **scaled** W — the one actually simulated, including `level_of_chaos` — imaged with a diverging colormap so zero is white and the E/I blocks read at a glance).

Many functions formerly in `src/` have been **internalized as static methods on `SRNNModel2`** (commit `6e2c58f`, 2026-02-27): `compute_Jacobian_fast`, `dynamics_fast`, `generate_external_input`, `compute_lyapunov_exponents_internal`, `benettin_algorithm_internal`, `lyapunov_spectrum_qr_internal`, `decimate_states`, `initialize_state`, `unpack_and_compute_states`, all `plot_*` and colormap helpers, and the activation functions. The standalone `src/` copies that were no longer referenced by any kept script have since been **deleted**; a few standalone helpers remain only because example/figure scripts still call them directly (e.g. `plotting/plot_firing_rate.m`, `algorithms/Jacobian/compute_Jacobian_fast.m`). When editing model behavior, **prefer the class method**, not any remaining standalone file. See `docs/refactors/internalize_functions_into_classes.md` for the original mapping.

The older single-class predecessor `src/SRNNModel.m` (default `n=100`, `indegree=20`) was removed in the cleanup; `SRNNModel2` (default `n=300`, `indegree=100`) is the E/I model class. Its Echo State Network subclass `src/SRNN_ESN_reservoir.m` adds memory-capacity tooling.

### `src/SRNNCellTypePairs.m` (the per-cell-type model)

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

`src/SRNNCellTypes.m` is the earlier per-type class. It still takes **absolute** tildes
(no `_relative` port) and still binds its activation handles by capturing `obj`. It is
kept only for `test_SRNNCellTypes_parity_lyapunov.m`; do not build new work on it.

### `src/ParamSpaceAnalysis2.m` (the analysis driver)

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

- **`src/srnn_param_preset.m`** — *which network*. `[d, model_class] = srnn_param_preset(name)` returns a struct of model overrides **plus the model class they are written for** (the name is only a lookup key; the struct is what reaches the model). Assign it to `psa.model_defaults` and layer tweaks on with plain field assignment. Presets carry physics **including** swept axes like `n`/`f`/`level_of_chaos` — a grid axis overrides the preset for the sweep that varies it, while sweeps holding that axis fixed use the preset's value. They must **not** carry `n_levels`/`n_reps` (not model properties), the run_mode timing knobs, or condition fields. The second output exists because the two model classes have **disjoint** parameter vocabularies, so a preset that did not carry its class would only fail later inside `validate_model_defaults`; it defaults to `'SRNNModel2'`. `'celltype_pairs'` is the `SRNNCellTypePairs` preset (E/I with STD on E→E and I→I).
- **`src/srnn_adaptation_conditions.m`** — *which adaptation regimes*. `srnn_adaptation_conditions(model_class)` returns the four conditions (`no_adaptation`, `sfa_only`, `std_only`, `sfa_and_std`) spelled for that class: `n_a_E`/`n_b_E` counts for `SRNNModel2`, an `n_a` row plus a whole `synapse_config` struct for `SRNNCellTypePairs`. The **names** are identical across classes, which is what keeps the `condition_titles` maps inside PSA's plotters working. Use this instead of `ParamSpaceAnalysis2`'s built-in defaults, which are `SRNNModel2`-shaped and would be passed verbatim to whatever `model_class` is set.
- **`scripts/run_all_analyses/analysis_run_config.m`** — *how much compute*. `analysis_run_config(analysis, run_mode, preset_defaults)` returns `n_levels`, `n_reps`, and a `cfg.model` struct of timing settings (`ode_solver`, `fs`, `T_range`, and `lya_T_interval` only when it should be set).

  **The two knobs are no longer fully orthogonal, and this is where they meet.** Every cell names *two* integrators — deterministic (`'rk4'` for fast/fast2/medium, `'ode45'` for production) and stochastic (`'sra1'` everywhere) — and the **preset** picks between them: `sigma_u_noise > 0` selects the stochastic one. That is why the third argument exists, and why the three sub-scripts resolve `preset_defaults` *before* calling this. Selecting here is what leaves everything else intact: `merge_struct` precedence is unchanged (`cfg.model` still wins), `ode_solver` stays banned from presets, and a σ = 0 preset is bit-identical to before the mechanism existed. `cfg.sde_solver` / `cfg.is_stochastic` live at the top level of `cfg`, never inside `cfg.model` — only `cfg.model` is merged into `model_defaults`, and neither is a model property. Note a stochastic `production` run is fixed-step and so does **not** carry `ode45`'s adaptive tolerance: there is no adaptive SDE solver, because step-size control is meaningless once increments are tied to the step.

Sub-scripts combine them with `merge_struct(preset_defaults, cfg.model)` — **preset first**, so `run_mode` keeps final say over its own knobs, and so a whole-struct assignment cannot clobber them. `run_all_analyses.m` picks the preset in one line (`preset_name`) and records it, with the model class, in `run_manifest`.

### Master-orchestrator conventions

`run_all_analyses.m` sets two variables that downstream scripts conditionally honor:

- `master_output_dir` — when set, sub-scripts write into this shared dir instead of creating their own and **skip their `clear`/`clc`/`close all`**.
- `master_save_figs` — `'save_all_figs'` / `'save_no_figs'` / `'follow_scripts_save_figs'` overrides each sub-script's local `save_figs` flag.
- `master_model_class` — the preset's model class. Each sub-script reads it into a local `model_class` (defaulting to `'SRNNModel2'` when run standalone), assigns `psa.model_class`, resets `integer_params` to `{'n','indegree'}` (the class default lists `SRNNModel2`'s adaptation counts), and takes its conditions from `srnn_adaptation_conditions(model_class)`. Two axes are renamed per class: the fraction-excitatory sweep is `f` on `SRNNModel2` and `f_E` on `SRNNCellTypePairs`, and `plot`'s `color_by` must be given explicitly for Pairs because its default `'f'` is a row there and breaks the histogram colouring.

Preserve this pattern when adding new analysis scripts.

### `src/` layout

- `algorithms/Jacobian/` — `compute_Jacobian_fast.m`, `compute_Jacobian_at_indices.m`, `compute_J_eff.m` (called directly by example scripts; `SRNNModel2` also carries internalized equivalents). The standalone `algorithms/Lyapunov/` files and `algorithms/info/mutual_info_SISO.m` were removed (internalized / unused).
- `connectivity/` — **`RMTBlocks.m`** is the generator all three model classes now use. It extends Harris (2023) to full D×D block statistics indexed **(postsynaptic, presynaptic)**, so `mu_EE` can differ from `mu_IE`, and it supplies the bulk radius `R` and the outlier eigenvalues `lambda_O`. `RMTMatrix.m` and `RMTCellTypes.m` remain but are no longer used by the model classes. Note `RMTBlocks` returns a **dense** matrix; the CellTypes classes re-sparsify after assembly.
- `nonlinearities/` — `piecewiseSigmoid`, `tanhActivation` and their derivatives (also internalized into `SRNNModel2`; the standalone `logisticSigmoid.m` was removed, but `logisticSigmoid`/`logisticSigmoidDerivative` live as static methods). Select one with `activation` (see "Conventions" below), not by building a handle. The **default** is `'logistic'` = `logisticSigmoid(x, S_c)` = `1/(1+exp(-4*(x-c)))` (centred at `S_c = 0.4`, range (0,1), unit slope at centre) — chosen for more robust near-edge-of-chaos stability; `'piecewise'` (`S_a = 0.9`) and `'tanh'` are the alternatives. Each model class carries its **own copy** of these statics, so a fix to one is not a fix to the others.
- `plotting/` — colormaps, line/scatter helpers, time-series panel plots, `param_space_plots/` for post-hoc visualization, `plot_saving/save_some_figs_to_folder_2.m` (used by every script that writes figures). Unused plot/colormap duplicates were removed; the standalone `plot_*` files that example/figure scripts call directly remain.
- `SRNN_ESN_reservoir.m` — Echo State Network subclass of `SRNNModel2` (memory-capacity experiments). The exploratory `SRNN_ESN.m` / `SRNN_reservoir.m` / `SRNN_reservoir_DDE.m` RHS variants were removed.

### Scripts retained beyond the `run_all_analyses` pipeline

The `refactor` cleanup removed the legacy subtrees (`old_scripts/`, `review_paper/`, `VAR_SRNN/`, `python_piecewise/`, `reference_files/`) and the old-API comparison/run scripts, and reorganized `scripts/` into topic subdirectories. Current layout:

- `setup_paths.m` — shared bootstrap, **at the repo root**, not under `scripts/` (self-locating; resolvable from a cold session with cwd at the root).
- `scripts/run_all_analyses/` — the orchestrator + its three sub-analyses, with `replot/` (the `replot_*` figure regenerators + `assemble_sensitivity_figure.m`) nested inside. Also `analysis_run_config.m`, the single per-script table of `run_mode` settings that replaced the duplicated `switch run_mode` blocks.
- `scripts/EI_balance/` — fraction-excitatory analyses: `fraction_excitatory_analysis.m`, `Fig_2_fraction_excitatory_analysis.m`, `Fig_2_fraction_excitatory_load_and_plot.m`.
- `scripts/memory_capacity/` — `example_memory_capacity.m`, `looped_memory_capacity.m` (Echo State Network experiments).
- `scripts/presentations/Stability_Manuscript/fig_stim_engages_adaptation/` —
  `bursting_SRNN_model*.m`.
- `scripts/tests/` — verification scripts (`test_SRNN2_defaults.m`, `test_psa_saveload.m`, `test_sensitivity_refactor.m`) plus the standalone example/comparison scripts `Sompolinsky_N_1000_g_1p8.m` and `Single_vs_dual_adaptation_example.m`.
- `scripts/sine_stim/` and `scripts/paired_pulse/` — kept as references but **currently non-functional**: they still use the old script-based API, and their legacy dependencies were deleted. They must be ported to `SRNNModel2` before use.

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
- **Prefer `'sra1'`** (Rößler SRA1): it costs the same two drift evaluations per step as `'heun'` but is strong order 1.5 rather than 1.0 — measured ~85× more accurate at the same step size. `src/sde_fixed_step.m` implements all three schemes in one file so the absolute-time noise indexing (`i0 = round((tspan(1) − t0)·fs) + 1`) is written once; that indexing is what lets Benettin re-integrate a segment against exactly the increments the original run used. The noise tensor is generated at the top of `run()` from `noise_seed` (empty → `rng_seeds(1) + 224737`), kept alive through `compute_lyapunov`, then cleared — it is ~96 MB at `n=300, T=50 s, fs=400` and is never saved because it is regenerable.
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
