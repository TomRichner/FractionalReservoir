# Codebase map: from the manuscript's figures back to the code (2026-11-10)

A high-level map of where the code that produces the paper lives, traced
backwards from the figure files `Manuscript5.md` embeds (in the sibling repo
`StochasticPlasticDynamicalSystemPaper`) to the scripts and classes that made
them. Written at the start of the final-draft work (to-do list:
`StochasticPlasticDynamicalSystemPaper/docs/final_draft_to_do.md`), as
orientation for the new experiments and for the reorganisation deferred until
after the 2026-11-16 deadline.

This is a snapshot. `CLAUDE.md` is the maintained reference for conventions;
this note only says which parts of the repo the current manuscript actually
exercises.

## 1. The main script

The manuscript's `figs/` link points at `figs/mu7revisedMedium`
(`figure_sources.json` in the manuscript repo). That figure set is produced by
one script, `scripts/paper/mu7revisedMedium_run.m`:

```matlab
setup_paths();
cfg     = mu7revisedMedium_config();      % preset, run mode, figure registry
run_dir = run_all_paper_analyses(cfg);    % all compute -> data/mu7revisedMedium
results = make_all_paper_figures(cfg);    % all figures -> figs/mu7revisedMedium
```

- `mu7revisedMedium_config.m` names the network preset, the run mode
  (`medium`), the side presets a few stages need, and the figure registry.
- `run_all_paper_analyses.m` and `make_all_paper_figures.m` (both in
  `scripts/paper/`) are the two entry points. Every other `*_config.m` /
  `*_run.m` pair in that folder is the same pattern for a different experiment
  bundle.

Presets used by this bundle (all in `src/presets/srnn_param_preset.m`):

| Role | Preset |
|---|---|
| The network (sweeps and most stages) | `celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0p025_dualStd_3cond_mu7revisedMedium` |
| Memory capacity (noise-free twin) | `..._noStim_noise0_dualStd_3cond_mu7revisedMedium` |
| Local Lyapunov (step-input twin) | `..._steps5s_noise0p025_dualStd_3cond_mu7revisedMedium` |
| Bursting example | `bursting_pairs` |
| Introductory concepts | `sompolinsky_pairs` |
| Single-neuron facilitation cartoon | `single_neuron_stf` |

The network preset is n = 500, in-degree 100, half the neurons excitatory,
mean weights `mu_tilde_relative = [7 -7; 7 -7]`, three conditions
(`no_adaptation`, `sfa1_std1`, `sfa3_std2`).

## 2. The analysis stages

`run_all_paper_analyses` runs the sweep pipeline first, then a list of further
stages into subfolders of the same run directory. Each stage is wrapped, so
one failure does not stop the rest.

| Stage (in `src/analysis/`) | What it does | Class it drives |
|---|---|---|
| `run_all_analyses`, which calls `run_sensitivity_analysis`, `run_tau_sensitivity_analysis`, `run_param_space_analysis` | 1-D sweeps, the adaptation-timescale sweep, the joint connectivity sample | `ParamSpaceAnalysis2` sweeping `SRNNCellTypePairs` |
| `run_memory_capacity` | paired-trial memory capacity ensemble | `SRNN_ESN_reservoir` |
| `run_memory_capacity_example` | one example reconstruction | `SRNN_ESN_reservoir` |
| `run_eig_heatmap` | Jacobian spectra along trajectories, imbalance examples | `SRNNCellTypePairs` |
| `run_numerics_verification` | solver order, Lyapunov-method agreement, Jacobian check | `SRNNNumericsProbe` |
| `run_lyapunov_spectrum` | spectrum shape per condition | `SRNNCellTypePairs` |
| `run_transient_gain` | gain with adaptation frozen vs active | `SRNNCellTypePairs` |
| `run_local_lyapunov` | local exponents under step input | `SRNNCellTypePairs` |
| `run_dc_lle_analysis` | LLE against DC input | `SRNNCellTypePairs` |
| `run_paper_illustrations` (when `cfg.illustrations`) | saved example dynamics for the illustrative panels | `SRNNCellTypePairs` |

The stages other than the sweeps build their model through
`src/figures/helpers/build_from_preset.m`.

## 3. The custom classes

| Class | Where | Lines | Role |
|---|---|---|---|
| `SRNNCellTypePairs` | `src/model/` | 3,836 | The network model: C cell types, per-route synapses, dynamics, Jacobian, Lyapunov methods, plots. Everything in the paper runs on it. |
| `SRNN_ESN_reservoir` | `src/model/` | 1,214 | Subclass of `SRNNCellTypePairs` adding the reservoir readout and training used for memory capacity. |
| `ParamSpaceAnalysis2` | `src/analysis/` | 2,849 | The sweep driver. Builds models by class name, runs batched parallel jobs with checkpoints, stores and plots results. |
| `RMTBlocks` | `src/model/connectivity/` | 1,096 | Generates the connectivity matrix; called by the model. |
| `SRNNNumericsProbe` | `src/analysis/` | 129 | Small subclass of `SRNNCellTypePairs` for the numerics stage. |
| `IEDExplorationNetwork` | `scripts/explorations/IED_model/` | 20 | Small subclass of `SRNNCellTypePairs` for the IED-like network. |
| `SRNNModel2` | `src/model/` | 4,097 | The older two-population model. Per `CLAUDE.md`, no paper preset names it; still swept and tested. |
| `RMTMatrix` | `src/model/connectivity/` | 195 | Older connectivity generator. Per `CLAUDE.md`, kept as a frozen test reference only. |

`SRNN2TestAccess` and `SRNNPairsTestAccess` in `scripts/tests/` are test
helpers.

## 4. The manuscript's figures

`make_all_paper_figures` loops over the registry in the config and calls each
`fig_*` function in `src/figures/` (48 `.m` files, flat) on the one run
directory.

### Main figures

Figures 1 to 6 are the "grouped" entries appended to the registry by
`scripts/paper/grouped_figure_registry.m`. Each `fig_mainN_grouped.m` is a
thin wrapper over `src/figures/fig_grouped_main.m`, which composes a group and
runs no simulations.

| Manuscript figure | Label | Entry | Composed from |
|---|---|---|---|
| 1 | `fig:intro` | `fig_main1_grouped` | `fig_main1_concepts`, using the `fig_introductory_concepts` output |
| 2 | `fig:adaptation_network` | `fig_main2_grouped` | `fig_main2_dynamics`, from the run's saved illustration stage |
| 3 | `fig:eig_heatmap` | `fig_main3_grouped` | `fig_main3_stability`, from the run directory |
| 4 | `fig:sensitivity_lle_medians` | `fig_main4_grouped` | `fig_sensitivity_medians` + `fig_lle_vs_rate` (sweep data) |
| 5 | `fig:transient_amplification` | `fig_main5_grouped` | `fig_transient_amplification` + `transient_gain_data.mat` |
| 6 | `fig:sfa_longest_timescale` | `fig_main6_grouped` | `fig_sfa_EOC_allStd` (tau sweep) + `fig_memory_capacity` |
| 7 | `fig:stim_clinical` | none here | `Fig_Tica_Main.png` from the separate `cscs_stim_stability` clinical pipeline, via the `Stim-Analysis-Figs` repo |
| 8 | `fig:stim_model` | none in the registry | `scripts/explorations/IED_model/`: `ied_stim_experiment.m` then `fig_ied_stim.m`, output in `figs/IED_model/stim_20261005` |

`fig_main7_grouped` and `fig_main8_grouped` exist (a placeholder learning
figure and a PSD figure) but the manuscript does not embed them; its Figures 7
and 8 are the two stimulation figures above.

### Supplementary figures

Individual registry entries, embedded directly:
`fig_single_neuron_revised`, `fig_SFA_steady_state`, `fig_STD_steady_state`,
`fig_eig_heatmap_imbalance`, `fig_local_vs_finite_lle`,
`fig_sensitivity_analysis_allStd` (LLE and mean-rate sheets),
`fig_sensitivity_medians` (mean rate), `fig_EI_weights_param_space`,
`fig_numerics_solver`, `fig_numerics_lya_method`, `fig_numerics_jacobian`,
`fig_memory_capacity_example`, and `fig_stim_engages_adaptation` (the
`bursting_pairs` network: time series and synchrony).

Registered and generated but not embedded in the manuscript, among others:
`fig_param_space_allStd`, `fig_eig_heatmap`, `fig_dc_lle`,
`fig_lyapunov_spectrum`, `fig_transient_gain`,
`fig_transient_gain_excursions`, `fig_local_lyapunov`, `fig_FI_curve`.

## 5. Stimulation work outside the two entry points

Three separate tracks touch stimulation, and none goes through
`run_all_paper_analyses`:

- **IED-like network** (`scripts/explorations/IED_model/`): its own config
  (`ied_base_config`, the frozen "run13" network), builder
  (`ied_build_network`, `IEDExplorationNetwork`), experiment
  (`ied_stim_experiment`) and figure (`fig_ied_stim`). This is manuscript
  Figure 8. Record in `docs/IED_model/`.
- **Bursting network** (`bursting_pairs` preset, `fig_stim_engages_adaptation`):
  inside the registry, so it is regenerated with the paper figures, but it
  simulates in the figure function and has no analysis stage.
- **Sine / paired-pulse stimulation** (`scripts/stim_figures/`,
  `run_stim_evoked`, `run_stim_neurons`, `fig_sine_paired_stim`,
  `fig_sine_stim_neurons`): not embedded in the manuscript.

## 6. What this means for the final-draft experiments

- **Gain sweep of LLE and memory capacity** (to-do block C). The sweeps
  (`ParamSpaceAnalysis2`) and memory capacity (`SRNN_ESN_reservoir`) are
  separate stages today: memory capacity is measured at the preset network
  only. The new experiment needs both measures for each (gain, seed,
  condition), so it has to bring the two together.
- **Stimulation causality and its parameter-space version** (block B). The
  stimulation result is the one main result outside the main pipeline, on a
  network (run13) that has its own construction path.
- **1TS versus MTS normalization** (block A). The matching variants live as
  presets plus one model option (`std.combine = 'geomean'` on
  `SRNNCellTypePairs`); see `STD_strength_matching_2026-09-13.md`.
