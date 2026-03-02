# Refactor: Grouped Output Directories & folder_prefix

**Date:** 2026-03-02

## Problem

All PSA runs created directories starting with `param_space_...`, making it hard to distinguish between tau sensitivity, 1D sensitivity, and full parameter space analyses. When run via `run_all_analyses.m`, outputs were scattered across separate top-level folders with no grouping.

## Changes

### `ParamSpaceAnalysis2.m`
- Added `folder_prefix` property (default `'param_space'`)
- `run()` now uses `folder_prefix` instead of hardcoded `'param_space'` when constructing the output folder name
- Added to `saveobj`/`loadobj` for persistence

### `run_all_analyses.m`
- Creates a shared `master_output_dir` (`data/param_space/run_all_<datestr>/`)
- Sub-scripts detect this variable in the workspace and nest their outputs inside it

### Analysis scripts

| Script                                 | `folder_prefix`   | Example folder                                           |
| -------------------------------------- | ----------------- | -------------------------------------------------------- |
| `run_tau_sensitivity_analysis.m`       | `tau_sensitivity` | `tau_sensitivity_tau_a_E_max_nLevs_25_...`               |
| `run_sensitivity_analysis.m`           | `1D_sensitivity`  | `1D_sensitivity_sensitivity_level_of_chaos_nLevs_21_...` |
| `run_param_space_analysis2.m`          | `param_space`     | `param_space_test_refactor_nLevs_3_...`                  |
| `Fig_2_fraction_excitatory_analysis.m` | `fig2`            | `fig2_frac_exc_nLevs_5_...`                              |

## Standalone vs Master

Each script checks `if exist('master_output_dir', 'var')` and sets `psa.output_dir` accordingly. When run standalone, outputs go to `data/param_space/` as before.
