# Refactor: SensitivityAnalysis → ParamSpaceAnalysis2 (1D Mode)

**Date:** 2026-03-02

## Summary

Replaced the standalone `SensitivityAnalysis` class with `ParamSpaceAnalysis2` operating in 1D mode. Sensitivity analysis is now a special case of parameter space analysis: a 1D grid with `reps` as a dummy second dimension.

## What Changed

### ParamSpaceAnalysis2 Extensions

| Feature                         | Description                                                                                                                 |
| ------------------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| `randomize_order` property      | Default `true`. Set `false` for sensitivity analysis to run levels sequentially.                                            |
| `add_vector_parameter()` method | Supports sweeping vector-valued parameters (e.g., `tau_a_E`) where one end of the vector varies across levels.              |
| `vector_param_config` property  | Stores config for each vector param: `vary_element`, `fixed_value`, `vary_range`, `n_elements`, `spacing`, `level_spacing`. |
| `vector_param_lookup` property  | Pre-generated cell array of vectors, indexed by grid position.                                                              |
| `plot_sensitivity()` method     | Generates imagesc heatmaps for 1D sweeps — one subplot per condition, with median overlay and λ=0 line.                     |

### `add_vector_parameter` API

```matlab
psa.add_vector_parameter('tau_a_E', ...
    'vary_element', 'last', ...     % 'first' or 'last'
    'fixed_value', 0.25, ...        % value of non-varying end
    'vary_range', [5, 60], ...      % [min, max] of varying end
    'n_elements', 3, ...            % length of each vector
    'spacing', 'log', ...           % 'linear' or 'log' within vector
    'level_spacing', 'linear')      % 'linear' or 'log' across levels
```

Integer params: if param is in `integer_params`, vary values and vector elements are rounded.

### Run Script Changes

| Script                           | Change                                                                                                                                                                 |
| -------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `run_sensitivity_analysis.m`     | Replaced `SensitivityAnalysis` with per-parameter `ParamSpaceAnalysis2` instances. Replaced `EI_imbalance` with `f`. No `model_defaults` — uses `SRNNModel2` defaults. |
| `run_tau_sensitivity_analysis.m` | Replaced manual parfor with `ParamSpaceAnalysis2`. Uses `add_vector_parameter` for `tau_a_E` and `add_grid_parameter` for `tau_b_E_rec`.                               |
| `run_all_analyses.m`             | Removed paired-pulse MI. Updated numbering (1/3). Changed to `run_param_space_analysis2`.                                                                              |

### Key Design Decision

Each parameter gets its own `ParamSpaceAnalysis2` instance (not one big multi-param sweep). This matches the old `SensitivityAnalysis` behavior where parameters are swept independently, keeping the heatmap interpretation simple: one figure per parameter × condition.

## Files Modified

- `src/ParamSpaceAnalysis2.m` — Extended with 5 new features
- `scripts/run_sensitivity_analysis.m` — Rewritten
- `scripts/run_tau_sensitivity_analysis.m` — Rewritten
- `scripts/run_all_analyses.m` — Updated

## Verification

Tested with 3 levels × 3 reps:
- Scalar parameter (level_of_chaos): ✅ all 36 sims successful, heatmap generated
- Vector parameter (tau_a_E): ✅ all 9 sims successful, heatmap with `tau_a_E(last)` x-axis
