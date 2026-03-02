# Refactor: Port Fig_2 Scripts & Internalize PSA Plotting

**Date:** 2026-03-02

## Summary

Ported `Fig_2_fraction_excitatory_analysis.m` and `Fig_2_fraction_excitatory_load_and_plot.m` from `ConnectivityAdaptation` to `FractionalReservoir`. Internalized two external plotting functions as `ParamSpaceAnalysis2` methods.

## What Changed

### New Scripts (in `scripts/`)

| File                                        | Description                                                                              |
| ------------------------------------------- | ---------------------------------------------------------------------------------------- |
| `Fig_2_fraction_excitatory_analysis.m`      | Runs the f-sweep PSA (f=[0.4,0.6], 4 conditions, local LLE + unit histograms + beeswarm) |
| `Fig_2_fraction_excitatory_load_and_plot.m` | Loads saved PSA results and regenerates all plots + Wilcoxon stats                       |

### Changes from ConnectivityAdaptation originals

- `ParamSpaceAnalysis` → `ParamSpaceAnalysis2`
- `@(x) piecewiseSigmoid(...)` → `@(x) SRNNModel2.piecewiseSigmoid(...)` (static method)
- `datestr(now)` → `string(datetime('now'))`
- `setup_paths()` uncommented for standalone use

### Internalized Methods on ParamSpaceAnalysis2

| Method                              | Replaces                             | What it does                                                               |
| ----------------------------------- | ------------------------------------ | -------------------------------------------------------------------------- |
| `plot_unit_histograms(varargin)`    | `load_and_make_unit_histograms()`    | Unit histograms colored by f-value, one figure per metric                  |
| `plot_lle_by_stim_period(varargin)` | `load_and_plot_lle_by_stim_period()` | Paired beeswarm of mean local LLE per stim step with significance brackets |

Both methods operate directly on `obj` (no file loading). The external functions still exist for backwards compatibility.

### Dependencies

All plotting utilities already existed in FractionalReservoir — no files copied from ConnectivityAdaptation:
- `unit_histogram_patch`, `paired_beeswarm`, `blue_gray_red_colormap`, `concatenate_figs`, `AddLetters2Plots`, `save_some_figs_to_folder_2`

## Verification

Tested via MATLAB MCP with 3 levels × 2 reps:
- ✅ 36/36 simulations successful
- ✅ `plot_unit_histograms` generated 2 figures (br, lle) with f-value coloring
- ✅ `plot_lle_by_stim_period` generated figure with p-values and brackets
- ✅ `concatenate_figs` combined into 3×4 layout
- ✅ Statistical output: Wilcoxon signed-rank p-values and Cohen's d for all 4 conditions

## Files Modified

- `src/ParamSpaceAnalysis2.m` — Added `plot_unit_histograms` and `plot_lle_by_stim_period` methods (~570 lines)

## Files Created

- `scripts/Fig_2_fraction_excitatory_analysis.m`
- `scripts/Fig_2_fraction_excitatory_load_and_plot.m`
