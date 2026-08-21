# refactorRunAll — investigation notes

Branch: `refactorRunAll`. Written 2026-08-20.

Goal (from the request): **two entry points**, working backwards from the paper
figures.

1. `run_all_analyses` — one click, all heavy compute, into one dated run directory.
2. `make_all_figures` — one click, all final paper figures, from saved data only.

Plus: one **main preset** with per-figure overrides where a figure genuinely needs a
different network; plotting scripts that **read the preset off the run** rather than
hardcoding it; documentation parameter tables **generated** from the preset; and
scripts that are really functions turned into functions.

---

## 1. The paper figure inventory

The manuscript currently pulls these 13 PNGs (paths given on the Mac, mapped here to
the repo). "allStd substitution" applied per the request: the `_allStd` sheet wins
where one exists.

| # | Paper figure | Producing script | Kind | Data source | Source on this machine? |
|---|---|---|---|---|---|
| 1 | `fig_introductory_concepts/panel_A/statetraces/panelA_bottom_traces_figure_1.png` | `fig_introductory_concepts/panel_A/panelA_bottom_traces.m` | **computes + plots** (fast, ~1 min) | none — simulates inline | n/a |
| 2 | `fig_introductory_concepts/panel_A/eigenspectra/panelA_eigenspectrum_figure_2.png` | **same script** (`panelA_bottom_traces.m`, fig 2) | **computes + plots** | none | n/a |
| 3 | `fig_example_timeseries/fig_example_timeseries_figure_1.png` | `fig_example_timeseries/fig_example_timeseries.m` | **computes + plots** (seconds) | none | n/a |
| 4 | `fig_FI_curve/Fig_FI_curve.png` | `fig_FI_curve/Fig_FI_curve.m` | **analytic**, no sim | none | n/a |
| 5 | `fig_sensitivity_analysis_allStd/Fig_Sensitivity_LLE_{core,mu}.png` | `fig_sensitivity_analysis_allStd/Fig_sensitivity_analysis_allStd.m` | **replot** | `data/param_space/run_all_aug_14_26_17_25/1D_sensitivity_*` | ✅ present + **git-tracked** |
| 6 | `fig_EI_param_space/Fig_EI_ParamSpace.png` | `fig_EI_param_space/Fig_EI_param_space.m` | **replot** | `data/param_space/run_all_jul_06_26_22_00/param_space_*` | ✅ present + git-tracked |
| 6′ | `fig_param_space_allStd/Fig_ParamSpace_allStd.png` | `fig_param_space_allStd/Fig_param_space_allStd.m` | **replot** | `run_all_aug_14_26_17_25/param_space_*` | ✅ |
| 7 | `fig_sfa_EOC_allStd/Fig_SFA_EOC_allStd.png` | `fig_sfa_EOC_allStd/Fig_sfa_EOC_allStd.m` | **replot** | `run_all_aug_14_26_17_25/tau_sensitivity_*` | ✅ |
| 8 | `fig_memory_capacity/Fig_Memory_Capacity_figure_3.png` | `fig_memory_capacity/Fig_memory_capacity.m` | **replot** | `data/memory_capacity/paper_ready/MC_sample_hold_20260722_154245_trials30_results.mat` | ❌ **`data/memory_capacity/` does not exist here** |
| 9 | `fig_memory_capacity_example/Fig_MC_Example.png` | `fig_memory_capacity_example/Fig_memory_capacity_example.m` | **replot** | `fig_memory_capacity_example/mc_example_data.mat` (written by sibling `compute_memory_capacity_example.m`) | ❌ **missing** (gitignored `.mat`) |
| 10 | `fig_stim_engages_adaptation/bursting_piecewise_seed42_42_figure_{1,2}.png` | `fig_stim_engages_adaptation/**bursting_SRNN_example.m**` | **computes + plots** (~min) | none | n/a |
| 11 | `fig_adaptation_methods/panel_A/sfa_std_single_neuron_example_figure_5.png` | `fig_adaptation_methods/panel_A/test_single_neuron_adaptation.m` | **computes + plots** (seconds) | none | n/a |

Figures present in the tree but **not** in the paper list — decide keep/drop:
`fig_eig_heatmap/` (needs `eig_heatmap_data.mat`, ❌ missing), `fig_SFA_steady_state/`,
`fig_STD_steady_state/` (+ `_zoom`), `fig_sensitivity_medians/`,
`fig_introductory_concepts/.../energy_landscape/`,
`fig_adaptation_methods/.../sfa_std_stf_single_neuron_example_*` (an STF variant with
no script left that names it), `fig_reservoir_diagram/` (hand-drawn SVG, no script),
`fig_equations/` + `doc_equations_table/` (markdown/LaTeX, hand-maintained).

### Which of the several similar scripts is the live one

Resolved by evidence, not by name:

- **`fig_stim_engages_adaptation`** has three near-identical 400–570-line drivers.
  `bursting_SRNN_example.m` and `bursting_SRNN_model_good_ex_piecewise.m` write to the
  **same filename** `bursting_piecewise_seed42_42_*`. Commits `bd256ec` ("fixed tau_a_E
  to match default") and `85e2fd2` ("fixed tau_b_E_rel to be default") each touched
  **`bursting_SRNN_example.m` and the three committed PNGs together** → `bursting_SRNN_example.m`
  is the live one. The other two (`..._good_ex.m` → `bursting_seed19_20_*`, ReLU-era;
  `..._good_ex_piecewise.m`) are stale ancestors and should be deleted.
- The non-`allStd` siblings (`fig_param_space/`, `fig_sensitivity_analysis/`,
  `fig_sfa_EOC/`, and the `_allStd2` triplet) were **already deleted** in history, so
  there is no ambiguity left there.
- `fig_EI_param_space` is *not* a superseded sibling of `fig_param_space_allStd` — it
  is a different plot (per-network patches coloured by `f`, 2×5 with an E:I colorbar)
  built from a **different, older run** on **`SRNNModel2`**. → open question Q1.

---

## 2. What produces the data

Two independent compute pipelines, and they do not share a model class.

### 2a. `scripts/run_all_analyses/` — the sweep pipeline (`SRNNCellTypePairs`)

`run_all_analyses.m` (233 lines) is the orchestrator. It sets base-workspace
`master_*` variables and then **`run`s three sibling scripts** that read those
variables back with `exist(...,'var')`:

```
run_all_analyses.m
  ├─ master_output_dir / master_save_figs / master_model_overrides
  │  master_model_class / master_conditions / run_mode        (base workspace)
  ├─ run_sensitivity_analysis        (7 × 1-D PSA: n, f_E, level_of_chaos, mu_{EE,EI,IE,II})
  │    └─ replot_sensitivity + assemble_sensitivity_figure ×2 (network / mu_blocks sheets)
  ├─ run_tau_sensitivity_analysis    (vector sweep over tau_a_E)
  ├─ run_param_space_analysis2       (3-D grid n × f_E × level_of_chaos)
  └─ write_run_parameters_md          → parameters.md
```

This part is **already in decent shape**: `analysis_run_config` centralises the
cost/fidelity table, `srnn_param_preset` centralises the physics,
`srnn_adaptation_conditions` centralises the four regimes, `capture_git_provenance`
and `write_run_parameters_md` make a run self-describing. The problems are structural,
not conceptual — see §4.

Not wired into the orchestrator: `run_dc_lle_analysis.m` (commented out at step 4),
`scripts/EI_balance/Fig_2_fraction_excitatory_analysis.m` (commented out),
`run_overnight_queue.m` (a wrapper that runs the orchestrator N times),
`check_sensitivity_sim.m` (a debugging one-off).

### 2b. `scripts/memory_capacity/` — the MC pipeline (`SRNN_ESN_reservoir < SRNNModel2`)

- `looped_memory_capacity.m` (407 lines) — 30 paired trials × 4 conditions → the
  `paper_ready/*_results.mat` behind figure 8.
- `compute_memory_capacity_example.m` (in the manuscript folder, 150 lines) →
  `mc_example_data.mat` behind figure 9.
- `example_memory_capacity.m`, `test_sample_hold_mc.m` — exploratory, not paper.
- `plot_memory_capacity.m` / `replot_memory_capacity.m` — already functions. Good.

**`SRNN_ESN_reservoir` subclasses `SRNNModel2`, not `SRNNCellTypePairs`.** So it speaks
`n_a_E`/`n_b_E`/`mu_E_tilde_relative`, and **cannot consume the main preset**
`celltype_pairs_Sc0p2_noise0p025_dualStd` at all. This is the single biggest obstacle
to "one preset for everything" → open question Q2.

---

## 3. Preset / model-class consistency across the paper — the real problem

| Figure | Model class | Nonlinearity | Where its parameters come from |
|---|---|---|---|
| 1–2 intro traces + eigenspectra | `SRNNModel2` | `tanh` | hardcoded in script (Sompolinsky reproduction — *deliberately* different) |
| 3 example timeseries | `SRNNModel2` | `logistic`, class defaults (n=300, f=0.5) | hardcoded |
| 4 F-I curve | none (analytic) | `logistic`, `S_c = 0.4` | hardcoded |
| 5 sensitivity allStd | `SRNNCellTypePairs` | `piecewise`, `S_a=0.8`, `S_c=0.2` | **preset `celltype_pairs_Sc0p2_noise0p025`** (aug_14 run) |
| 6 EI param space | `SRNNModel2` | `logistic` era | jul_06 run, **pre-preset** |
| 6′ param space allStd | `SRNNCellTypePairs` | as row 5 | preset `celltype_pairs_Sc0p2_noise0p025` |
| 7 SFA EOC allStd | `SRNNCellTypePairs` | as row 5 | preset `celltype_pairs_Sc0p2_noise0p025` |
| 8 memory capacity | `SRNN_ESN_reservoir` | `logistic`, `S_c=0.35`, n=300, f=0.6, `loc=2.0` | hardcoded |
| 9 MC example | `SRNN_ESN_reservoir` | as row 8 (`loc=2.0`) | hardcoded |
| 10 bursting | `SRNNModel2` | `piecewise`, `S_c=0.5`, n=50, indegree=10 | hardcoded |
| 11 single-neuron adaptation | `SRNNModel2` | `piecewise`, `S_a=1.0`, `S_c=0.35`, n=1 | hardcoded |
| (extra) SFA/STD steady state | analytic | — | **reads `celltype_pairs_Sc0p2_noise0p025_dualStd` from `srnn_param_preset`** ← the pattern to copy |

**Findings:**

- **Only 3 of 13 figures currently sit on the target preset, and it is the wrong one.**
  Figures 5/6′/7 come from `run_all_aug_14_26_17_25`, whose `parameters.md` records
  `preset = celltype_pairs_Sc0p2_noise0p025`, `run_mode = production`. The target is
  `celltype_pairs_Sc0p2_noise0p025_dualStd`. → **the sweep pipeline must be re-run**
  before those three figures are final. The only `_dualStd` run on disk is
  `run_all_aug_18_26_21_41`, which is `run_mode = medium` (11 levels × 15 reps).
- Figures 3, 8, 9, 10, 11 hardcode `SRNNModel2` parameters that **contradict** the
  swept network: figure 3 shows `logistic`/n=300/`loc=1` while every sweep figure is
  `piecewise`/n=500/`S_c=0.2`/noisy. A reader comparing "example timeseries" to
  "parameter space" is looking at two different models.
- `Fig_STD_steady_state.m` already does the right thing —
  `[~,~,conditions] = srnn_param_preset(preset_name)` and pulls `tau_rec`/`tau_rel`
  out of the returned conditions. **Generalise this.**
- `preset_default_values(data_root, params)` already does the *other* right thing —
  reads `run_manifest.mat` from the run directory, reconstructs the preset, builds a
  model, and reads defaults off the object so class accessors resolve the aliases. Two
  figure scripts use it. **This is the "plotting scripts automatically grab the correct
  preset" mechanism; it just isn't used everywhere.**

---

## 4. What is wrong with the current structure

### 4a. Scripts that should be functions

| File | Lines | Why |
|---|---|---|
| `run_sensitivity_analysis.m` | 234 | Communicates with the orchestrator through **base-workspace globals read with `exist()`**. A stale `master_*` silently changes a run (`run_overnight_queue.m` exists *only* to clear them between entries — a workaround for the design). Should be `run_sensitivity_analysis(cfg)`. |
| `run_tau_sensitivity_analysis.m` | 203 | same |
| `run_param_space_analysis2.m` | 213 | same |
| `run_all_analyses.m` | 233 | Should be `run_all_analyses(preset_name, run_mode, opts)` returning the run dir. Then `run_overnight_queue` becomes a 10-line loop and needs no state hygiene. |
| `run_dc_lle_analysis.m` | 253 | same, if kept |
| `looped_memory_capacity.m` | 407 | Config block at the top, no arguments; not callable from a master script. |
| `compute_memory_capacity_example.m` | 150 | same |
| every `Fig_*.m` in `Stability_Manuscript/` | 30–569 | All hardcode `data_root` / `out_dir` as script-locals. To be driven by one `make_all_figures`, each needs to be `fig_x(data_root, out_dir, opts)`. |
| `bursting_SRNN_example.m` | 562 | `clearvars -except rng_seeds` at line ~27 — a script that *clears the caller's workspace*. Cannot be called from anything. |
| `replot_all_analyses.m` | 62 | Hardcodes `data_root = run_all_mar_02_26_17_12` — **a run that no longer exists**. Dead as written. |
| `combine_two_runs.m` | 38 | Config-only wrapper around `combine_runs`. |

### 4b. Duplication to collapse

- **README generation.** 10 figure scripts contain 18–91 `fprintf(fid, ...)` calls each
  (~350 lines total) hand-writing `README_*.txt`. Should be one
  `write_figure_readme(out_dir, info)` helper fed a struct.
- **Style constants.** `tick_fs = 14`, `label_fs = 15.4` appear in 9 scripts; the
  Okabe-Ito condition palette in 3; condition display names in ~5. Should be one
  `manuscript_style()` / `condition_style()`.
- **Local helper copies.** `sort_axes_left_to_right` (2 copies), `save_figure_stable` /
  `save_fig_stable` (4 copies, 2 spellings), `dc_staircase_stimulus` (3 copies as a
  local function *and* one in `src/stimulus/`), `population_chi2` (2 copies).
  `Fig_sensitivity_analysis_allStd.m` still carries local copies of
  `preset_default_values` / `apply_percent_axis` / `mark_default_value` /
  `save_figure_stable` even though all four now exist standalone in
  `scripts/run_all_analyses/replot/` (its own header admits this).
- **The three `run_*_analysis` scripts** share ~60 lines of identical preamble
  (save_figs resolution, `run_mode` default, preset resolution, `model_class`
  resolution, conditions resolution, `merge_struct`). One `resolve_run_context()`.

### 4c. Data fragility

- `data/memory_capacity/` **is not on this machine and is not in git.** Figures 8 and 9
  cannot be regenerated here without re-running the MC pipeline
  (`looped_memory_capacity` = 30 trials × 4 conditions × 750 s of simulated time —
  hours).
- `mc_example_data.mat` and `eig_heatmap_data.mat` are gitignored intermediates that
  live *next to their figure scripts*, so they are invisible to any run-directory
  convention.
- Six `run_all_*` directories are force-committed past `.gitignore`
  (`jul_06_26_22_00`, `aug_14_26_00_48`, `aug_14_26_12_04`, `aug_14_26_12_14`,
  `aug_14_26_17_25`, `aug_18_26_21_41` — 744 tracked files). Only two of those are
  actually referenced by a figure. The rest is dead weight in the repo.
- `data/param_space/` holds **26** run directories; the figures reference 2.

### 4d. Smaller issues

- Two spellings of the same helper (`save_figure_stable` vs `save_fig_stable`) with
  different behaviour.
- `test_single_neuron_adaptation.m` is named `test_*` but is a **paper figure**
  generator sitting in the manuscript tree — it will read as a unit test.
- `run_param_space_analysis2.m` carries a stale `2` from the class rename; its `note`
  is still `'test_refactor'`, which is why every param-space output folder on disk is
  called `param_space_test_refactor_nLevs_*`.
- The `doc_equations_table/` and `fig_equations/` markdown tables are **hand-written and
  already stale**: `equation_table.md` documents a *single* STD timescale, the
  *logistic* activation, and `c_E = 0.5/3` for `SRNNModel2` — the target preset is dual
  STD, piecewise, `SRNNCellTypePairs`. `adaptation_conditions.md` likewise lists
  `n_a_E/n_a_I/n_b_E/n_b_I` counts that only exist on `SRNNModel2`.
  `write_run_parameters_md.m` (1070 lines) **already generates exactly this content**
  from a run — the manuscript tables should be emitted by it (or a sibling) rather than
  maintained by hand.

---

## 5. Proposed target structure

```
scripts/
  paper/                              <-- NEW: the two entry points live here
    run_all_paper_analyses.m          [1] one click, all heavy compute
    make_all_paper_figures.m          [2] one click, all final figures
    paper_config.m                    the single place naming the preset,
                                      run_mode, and per-figure overrides
    quick_check.m                     [3] fast-mode variant of [1] for dial-in

  analyses/                           (was run_all_analyses/, now all functions)
    run_sensitivity_analysis.m        function(cfg) -> dirs
    run_tau_sensitivity_analysis.m    function(cfg) -> dir
    run_param_space_analysis.m        function(cfg) -> dir
    run_memory_capacity.m             function(cfg) -> mat path   (was looped_memory_capacity)
    run_memory_capacity_example.m     function(cfg) -> mat path
    analysis_run_config.m             (unchanged)
    resolve_run_context.m             NEW: the shared preamble

  presentations/
    Stability_Manuscript/             <-- STAYS PUT (Q6): the manuscript on the Mac
      fig_*/fig_*.m                   references these paths. Contents refactored,
                                      each now a function(run_dir, out_dir, opts).
      _common/                        NEW, inside the manuscript tree so it stays
        manuscript_style.m            self-contained
        write_figure_readme.m         replaces ~350 lines of fprintf(fid,...)
        sort_axes_left_to_right.m
        save_figure_stable.m          one spelling
        resolve_run_dir.m             "newest run whose manifest says preset X"
      doc_tables/
        write_manuscript_tables.m     NEW: parameter + condition tables from the preset
```

Key design decisions:

- **`paper_config.m` is the one file the user edits.** It returns
  `{preset_name, run_mode, figure_overrides}` where `figure_overrides` is a map from
  figure name → preset name, so the bursting figure and the MC figures can name their
  own preset explicitly instead of hardcoding parameters. Everything else reads it.
- **Figures never guess a run.** `make_all_paper_figures` resolves *one* run directory
  (newest whose `run_manifest.mat` names the configured preset, or an explicit path)
  and passes it down. `resolve_run_dir` errors loudly if none matches, rather than
  silently plotting last month's preset.
- **Every figure gets its preset from the run, not from source.**
  `preset_default_values`'s approach generalised: read `run_manifest.mat`, rebuild the
  model, read values off the object. For the *computing* figures (1,2,3,10,11) the
  preset is applied directly at construction.
- **`write_run_parameters_md` gains a manuscript-table mode** so
  `fig_equations/parameter_table.md` and `doc_equations_table/*.md` are generated
  artifacts, not hand-maintained prose.
- **The `master_*` base-workspace protocol dies.** Functions take a `cfg` struct.
  `run_overnight_queue` collapses into a loop.

### What `run_all_paper_analyses` would contain

| Stage | Cost at `production` | Notes |
|---|---|---|
| 7 × 1-D sensitivity | ~20 h (25 lev × 25 reps × 4 cond) | dominant |
| tau sensitivity | ~few h | |
| param space (3-D grid) | ~1 h | 5³ = 125 pts, no reps |
| memory capacity (30 trials) | hours | different model class |
| MC example (4 conditions) | ~30 min | |
| eig heatmap | ~30 min | only if kept |

(Timings extrapolated from `analysis_run_config`'s own notes, which were measured on
the aug-13 preset; re-measure.)

### What `make_all_paper_figures` would contain

All 13 figures. The five that *simulate* (intro traces, eigenspectra, example
timeseries, bursting, single-neuron adaptation) are each ≲ a few minutes, so they
belong in the light script — that matches the request ("light-weight computing and
figure making").

---

## 6. Decisions (answered 2026-08-20)

| # | Decision |
|---|---|
| Q1 | **Both.** `fig_param_space_allStd` is the paper figure. `fig_EI_param_space` is **kept and rebuilt from the same `_dualStd` run** — i.e. the `f`-coloured 2×5 layout, but sourced from the new `SRNNCellTypePairs` run instead of the old `jul_06` `SRNNModel2` one. |
| Q2 | **Own `SRNNModel2` preset now, port to `SRNNCellTypePairs` later.** MC keeps its own named preset so it is at least preset-driven and appears in the generated tables; the port is tracked follow-up work. |
| Q4 | **Keep every figure folder.** Specifically: `fig_SFA_steady_state` / `fig_STD_steady_state` read their τ from the preset and, for the single-timescale comparison curve, use **`tau(1)` and `c = value/1`** (not `value/3`); `fig_sensitivity_medians` → `_dualStd`; `fig_eig_heatmap` → `_dualStd`; `energy_landscape` kept. STF variant → see below. |
| Q5 | **`fast` for smoke tests only.** The refactor is validated at `fast`; the user runs `production` overnight personally. So the pipeline must be re-runnable unattended with one edited line. |
| Q3 | **Re-run MC here at fast settings.** Add `run_memory_capacity(cfg)` with a fast mode (few trials, short `T_train`) so the plumbing is verifiable now; the user runs the full 30-trial version overnight alongside the production sweep. |
| Q6 | **`scripts/presentations/Stability_Manuscript/` stays where it is.** Only its contents are refactored, so no manuscript image path breaks on the Mac. `scripts/paper/` holds the new entry points and calls into it. |
| Q7 | **All of them move to `SRNNCellTypePairs`** — see the table below. |
| STF | **Parked, outside this refactor.** Written up in `UserNotes.md` ("Rebuild the single-neuron STF methods figure"); the three orphan files are left untouched. |

### Q7 in detail — every figure onto `SRNNCellTypePairs`

The whole paper moves to one model class. Three figures get their own **new** presets
because they are deliberately different networks; the rest take the main preset.

| Figure | Preset | Notes |
|---|---|---|
| `fig_example_timeseries` (fig 3) | `celltype_pairs_Sc0p2_noise0p025_dualStd` | straight port |
| `fig_adaptation_methods` single-neuron (fig 11) | `celltype_pairs_Sc0p2_noise0p025_dualStd` | **rebuild on `SRNNCellTypePairs`**. SFA and STD only — **no STF column**. `c` **matches the preset** (`0.5/3`), not the exaggerated `c_E = 1.0` in use today. **No new preset.** |
| `fig_stim_engages_adaptation` bursting (fig 10) | **NEW preset**, named for bursting | port to `SRNNCellTypePairs`, preset **equivalent to the current figure** (`n=50`, `indegree=10`, `f=0.7`, piecewise `S_a`/`S_c=0.5`, the DC-staircase `input_config`, …) |
| `fig_introductory_concepts` panel A (figs 1–2) | **NEW preset** | reproduce the Sompolinsky network **on `SRNNCellTypePairs`**, preset matching the figure's current settings (`N=200`, `tanh`, `tau_d=1`, fully connected, zero-mean random connectivity, `x0_std=1`) |
| `fig_eig_heatmap` | `celltype_pairs_Sc0p2_noise0p025_dualStd` | port from `SRNNModel2` (Q4) |
| `fig_sensitivity_medians` | `celltype_pairs_Sc0p2_noise0p025_dualStd` | already replot-only; re-point at the `_dualStd` run (Q4) |
| `fig_SFA_steady_state`, `fig_STD_steady_state` | `celltype_pairs_Sc0p2_noise0p025_dualStd` | analytic; read τ from the preset. Single-timescale comparison curve uses **`tau(1)`** and **`c = value/1`** (Q4) |
| `fig_memory_capacity`, `fig_memory_capacity_example` | **NEW `SRNNModel2` preset** | `SRNN_ESN_reservoir` blocks the port; own preset now, port later (Q2) |

### Feasibility of the two `SRNNCellTypePairs` ports — **verified in MATLAB, not assumed**

Probed on 2026-08-21 (scripts in the session scratchpad). **Both ports are feasible
today with no change to shared code.** Findings:

**1. `SRNNCellTypePairs` cannot build a one-cell-type model.** This is a real latent
defect, and it is what would have blocked the obvious approach to both ports.

```
SRNNCellTypePairs(... 'n_cellTypes',1, 'f',1, 'mu_tilde_relative',0, ...)
  -> object holds n_cellTypes=1, numel(f)=1, mu_tilde [1 1]   (all correct)
  -> build() errors:  RMTBlocks:InconsistentTypes
     "numel(f) = 2 but mu_tilde is [1 1]. Use set_types(...) to change the number
      of cell types."
```

Cause, located exactly: `SRNNCellTypePairs.build_W` (`src/SRNNCellTypePairs.m:1090-1094`)
assigns the generator **piecemeal** —

```matlab
rmt = RMTBlocks(obj.n);   % defaults to D = 2  (RMTBlocks.m:155, f = [0.5 0.5])
rmt.f           = obj.f;  % scalar 1 -> the D=2 setter expands it to [1 0]
rmt.mu_tilde    = obj.mu_tilde;      % 1x1
rmt.sigma_tilde = obj.sigma_tilde;   % 1x1
```

`RMTBlocks.set_types`'s own docstring says piecemeal assignment "works fine when D is
unchanged" and that `set_types` "is the only way to CHANGE D". Confirmed directly:
`r.f = 1` leaves `numel(r.f) == 2` with `f = [1 0]`, whereas
`r.set_types(1, 0, 1)` builds a D=1 network correctly.

**Reported, not fixed** — this is shared model-class code, and per CLAUDE.md's Scope
rule a defect found mid-task gets reported rather than folded into another commit. The
one-line change (`rmt.set_types(obj.f, obj.mu_tilde, obj.sigma_tilde)`) is a candidate
for its own commit if TR wants D=1 supported; **neither port needs it.**

**2. Sompolinsky panel A — works today as D=2 with identical zero-mean blocks.**

```matlab
SRNNCellTypePairs('n_cellTypes',2, 'cell_type_names',{'A','B'}, 'f',[0.5 0.5], ...
    'mu_tilde_relative', zeros(2), 'sigma_tilde_relative', ones(2), ...
    'n',200, 'indegree',200, 'level_of_chaos',1.6, 'activation','tanh', 'tau_d',1)
```

Result: builds; `nnz(W) = 40000` of 40000 (dense survives the re-sparsify step);
**fraction of negative weights = 0.500** (Dale-free, both signs — the thing that makes
the reproduction expressible at all); `R` theory 1.600 vs realized spectral radius
1.692 (ordinary finite-size overshoot at N=200). Type *names* are arbitrary — `A`/`B`
build fine, so the two populations need not be mislabelled E/I. Zero-mean blocks make
the two types statistically identical, which is exactly a Sompolinsky network.

**3. Single-neuron adaptation figure — floor is n=2, not n=1.** `n >= n_cellTypes` is
enforced, and `indegree = 0` is rejected (`0 < indegree <= n`). But this builds:

```matlab
SRNNCellTypePairs('n_cellTypes',2, 'cell_type_names',{'E','I'}, 'f',[0.5 0.5], ...
    'mu_tilde_relative', zeros(2), 'sigma_tilde_relative', zeros(2), ...
    'n',2, 'indegree',1, 'activation','piecewise', 'S_a',0.8, 'S_c',0.2)
```

→ `n_per_type = [1 1]`, `W = [0 0; 0 0]`, `N_sys_eqs = 2`. So the figure becomes a
**two-neuron network with identically zero weights, plotting the E neuron only**.
Functionally identical to the `n=1, W=0` model it replaces.

**4. `tanh` + heterogeneous `S_c` errors as documented** —
`"Per-type setpoints (mu_S_c / sigma_S_c) require a nonlinearity with a centre; 'tanh'
has none."` So the Sompolinsky preset must leave `mu_S_c`/`sigma_S_c` empty. Confirmed.

**`c` on the single-neuron figure — settled, build it as-is.** The figure currently uses
`c_E = 1.0` *specifically* so the rate decay is visible; the preset's `0.5/3 ≈ 0.167`
split over 3 timescales is ~6× weaker, so the adaptation will read more subtly. TR has
confirmed **6× weaker is fine** and will revise later if wanted. So: take `c` from the
preset, do **not** exaggerate it, and do not "helpfully" scale it back up if the
resulting panel looks flat — that is the honest value and the decision is already made.

### Q4 follow-up — the `sfa_std_stf_*` orphan

`sfa_std_stf_single_neuron_example_{figure_1.png, figure_1.svg, f_1.fig}` were produced
by `scripts/presentations/Stability_Manuscript/fig_adaptation_methods/test_single_neuron_stf.m`,
added in `390d86a`, last touched in `60c2992`, **deleted in the `refactor` cleanup**.
Not present on any branch tip (local or remote).

It showed **seven** columns (No adaptation, SFA, STD, **STF**, SFA+STD, STD+STF,
SFA+STD+STF) against the surviving figure's four.

It cannot be restored as-is. It bypassed the model classes entirely — hand-built an
`n=1, K=1` params struct and called `SRNNModelCellTypes.dynamics_fast_ct` /
`unpack_states_ct` with `SRNNModelBase.piecewiseSigmoid`. **All of those classes are
deleted.** And its facilitation equations are superseded:

| | deleted script | current `SRNNCellTypePairs` |
|---|---|---|
| facilitation state | `p`, rest `p0 = 0.35` | `g`, ceiling `G` |
| STD depletion | `db/dt = (1−b)/τ_rec − (p·b·r)/τ_rel` — **coupled to `p`** | `db/dt = (1−b)/τ_rec − (b·r)/τ_rel` — independent |
| synaptic gain | `eff = (p/p0)·b` | `g·b`, `dg/dt = (1−g)/τ_dec + (G−g)·r/τ_fac` |

**Recommendation: delete the three orphan files.** Rebuilding the panel means
re-deriving STF on `SRNNCellTypePairs`'s facilitation — a new figure, not a
restoration — and the target preset has **no STF on any route**, so nothing in the
paper's model facilitates. Awaiting confirmation.

## 7. Still open

Nothing. See §6 for the decisions; Q8 below is settled.

### Q8 (settled) — the `f`-coloured param-space rebuild needs fixing, not just re-pointing

`load_and_make_unit_histograms` colours by `f`. On `SRNNCellTypePairs` `f` is a `1 x C`
**row**, not a scalar, so the colour axis has nothing scalar to sort on — CLAUDE.md's
note that `color_by` "must be given explicitly for Pairs because its default `'f'` is a
row there and breaks the histogram colouring". **So `fig_EI_param_space` is not a
re-point job: the colouring path has to be fixed.**

The fix is to colour by **`f_E`**, and — verified in `src/SRNNCellTypePairs.m:337-344` —
**`f_E` is simply `f(1)`**:

```matlab
function v = get.f_E(obj)
    obj.assert_two_types();          % errors unless n_cellTypes == 2
    if isempty(obj.f); v = []; else; v = obj.f(1); end
end
```

So on any two-type preset — which every `celltype_pairs_*` preset is — `f_E` is exactly
the fraction excitatory the old `SRNNModel2` scalar `f` meant. Same quantity, same
range, so the existing `blue_gray_red_colormap` ramp and the E:I ratio colorbar ticks
(`1:3, 1:2, 2:3, 1:1, 3:2, 2:1, 3:1`) carry over unchanged; only the **name read off
the result** changes.

Two things to watch when doing it:

- Read it through `psa.effective_param(res, 'f_E')`, not `res.config.f` — CLAUDE.md is
  explicit that `effective_param(res,'f')` on a Pairs run returns the **class default**
  `[0.5 0.5]` rather than the swept value.
- `get.f_E` calls `assert_two_types()`, so it **errors** on a run with any other number
  of cell types. That is the right failure (loud, not silent miscolouring), but the
  plumbing should pass `color_by` in rather than hardcoding `'f_E'`, so a future
  three-type run can name its own axis.

---

## 8. Work plan (order of operations)

Each phase leaves the tree runnable; nothing is half-migrated across a phase boundary.

**Phase 1 — turn the pipeline into functions.** `run_sensitivity_analysis`,
`run_tau_sensitivity_analysis`, `run_param_space_analysis` become
`f(cfg) -> output_dir`; `run_all_analyses` becomes
`run_all_analyses(preset_name, run_mode, opts) -> run_dir`. Delete the `master_*`
base-workspace protocol and the state-scrubbing in `run_overnight_queue`. Extract the
shared ~60-line preamble into `resolve_run_context`. **Verify with a `fast` smoke run**
and confirm it produces the same directory layout as `run_all_aug_18_26_21_41`.

**Phase 2 — `scripts/paper/` entry points.** `paper_config.m` (the one file to edit:
main preset, run_mode, per-figure preset overrides), `run_all_paper_analyses.m`,
`make_all_paper_figures.m`, `quick_check.m`. Fold the MC pipeline in as
`run_memory_capacity(cfg)` with its own named `SRNNModel2` preset (Q2).

**Phase 3 — figure helpers.** `_common/`: `manuscript_style`, `write_figure_readme`,
`sort_axes_left_to_right`, one `save_figure_stable`, `resolve_run_dir`. Collapse the
~350 lines of `fprintf(fid,...)`, the 9 copies of `tick_fs/label_fs`, the 4 copies of
the save helper, and the 3 copies of `dc_staircase_stimulus`.

**Phase 4 — figures become functions and read their preset.** Each `Fig_*.m` →
`f(run_dir, out_dir, opts)`, staying inside `Stability_Manuscript/` (Q6). Replot
figures resolve the preset from the run's `run_manifest.mat` (generalise
`preset_default_values`); computing figures take the preset at construction. Then
per Q4/Q7:

- *Replot-only, just re-point:* `fig_sensitivity_analysis_allStd`,
  `fig_param_space_allStd`, `fig_sfa_EOC_allStd`, `fig_sensitivity_medians` → the
  `_dualStd` run.
- *Analytic, read τ from the preset:* `fig_SFA_steady_state`, `fig_STD_steady_state`
  — single-timescale curve uses `tau(1)` and `c = value/1`.
- *Rebuild from `SRNNModel2` onto `SRNNCellTypePairs`:* `fig_example_timeseries`,
  `fig_adaptation_methods` single-neuron (SFA+STD only, `c` from the preset),
  `fig_eig_heatmap`.
- *Rebuild onto `SRNNCellTypePairs` behind a **new** preset:*
  `fig_stim_engages_adaptation` (bursting preset),
  `fig_introductory_concepts` panel A (Sompolinsky preset).
- *Fix the colouring path, then rebuild:* `fig_EI_param_space` — plumb `color_by`
  through `load_and_make_unit_histograms` and read the value via
  `effective_param(res, 'f_E')`. `f_E` is `f(1)`, so the colormap and E:I ticks are
  unchanged. This one is **not** a re-point; the Pairs `f`-row breaks the existing
  colouring outright (Q8).

**Phase 5 — generated documentation.** Extend `write_run_parameters_md` (or a sibling)
to emit `fig_equations/parameter_table.md` and `doc_equations_table/*.md` from the
preset + conditions, replacing the hand-written tables that currently document a
single-STD, logistic, `SRNNModel2` network the paper no longer uses.

**Phase 6 — deletions.** `bursting_SRNN_model_good_ex.m` and
`bursting_SRNN_model_good_ex_piecewise.m` (stale ancestors, the latter writing a *live*
filename); the four local helper copies inside `Fig_sensitivity_analysis_allStd.m`;
`replot_all_analyses.m`'s dead hardcoded path. **The `sfa_std_stf_*` orphans are NOT
deleted** — parked per TR, see `UserNotes.md`.

## 9. Log of what was checked (so a later session need not redo it)

- All 21 `.m` under `Stability_Manuscript/` read for `data_root` / preset / save-name.
- All 10 `README_*.txt` read — they are the reliable record of each figure's source run.
- `git log --diff-filter=D` confirms the non-`allStd` figure siblings were already
  deleted; no live duplicates there.
- `git show --stat` on `bd256ec` / `85e2fd2` proves `bursting_SRNN_example.m` (not
  `..._good_ex_piecewise.m`) produced the committed bursting PNGs, despite the two
  writing the same filename.
- `parameters.md` in `run_all_aug_14_26_17_25` and `run_all_aug_18_26_21_41` read —
  they record preset + run_mode, which is what identifies a run.
- `git ls-files data` → 6 force-committed run dirs; `run_all_aug_14_26_17_25` **is**
  tracked, so the current paper figures are reproducible from git.
- `data/memory_capacity`, `mc_example_data.mat`, `eig_heatmap_data.mat` all confirmed
  absent.
