# Debug: High Chaos in FractionalResevoir vs Low Chaos in ConnectivityAdaptation

## Problem
- `example_memory_capacity.m` (FractionalResevoir): LLE ≈ 10 for all conditions
- `Simple_network_with_dual_adaptation.m` (ConnectivityAdaptation): LLE ≈ 0 for STD/SFA+STD

## Parameter Comparison Table

| Parameter | `example_memory_capacity.m` | `Simple_network_...m` | Effect on LLE | Status |
|-----------|-----------------------------|-----------------------|---------------|--------|
| **Network** |
| `n` | 200 | 300 (default) | Lower n → potentially higher chaos | |
| `indegree` | inherits SRNNModel2 default=100 | 100 (default) | Same | ✓ |
| `f` | 0.5 (default) | 0.5 (default) | Same | ✓ |
| `level_of_chaos` | **1.8** | **1.0** (default) | ⚠️ **1.8x higher W scaling** | |
| **Dynamics** |
| `tau_d` | **0.025** | **0.1** (default) | ⚠️ **4x faster dynamics** | |
| `S_a` (activation) | 0.9 (default) | 0.9 (default) | Same | ✓ |
| `S_c` (activation) | 0.35 (default) | 0.35 (default) | Same | ✓ |
| **Adaptation** (SFA) |
| `n_a_E` | 3 (for SFA conditions) | 3 | Same | ✓ |
| `tau_a_E` | [0.1, 1.0, 10] | logspace(0.25,10,3)≈[0.25,1.58,10] | Slightly different timescales | |
| `c_E` | **0.15/3** | **0.15/3** (default) | Same | ✓ FIXED |
| **STD** |
| `n_b_E` | 1 (for STD conditions) | 1 | Same | ✓ |
| `tau_b_E_rec` | 1.0 | 1 (default) | Same | ✓ |
| `tau_b_E_rel` | 0.25 | 0.25 (default) | Same | ✓ |
| **Simulation** |
| `fs` | 400 (SRNNModel2 default) | 400 (default) | Same | ✓ |
| `T_range` | determined by ESN protocol | [0, 50] (default) | Different | |
| **Lyapunov** |
| `lya_method` | computed internally by ESN | 'benettin' | Same method | ✓ |
| `lya_T_interval` | [15, T_end] (skips first 15s) | [15, 50] (skips first 15s) | Same | ✓ FIXED |
| **Input** |
| `input_type` | 'bandlimited' (ESN) | step pulses (default) | Different stimulus | |
| `u_f_cutoff` | 5 Hz | N/A | Different stimulus | |
| `positive_only` | false (default) | true | Different stimulus polarity | |

## Key Differences Causing High LLE

### 1. **`level_of_chaos = 1.8`** (most important)
This scales the weight matrix W by 1.8x, directly increasing spectral radius and instability.

### 2. **`tau_d = 0.025`** (very important)
4x faster dendritic dynamics means 4x faster accumulation of chaos.
Jacobian eigenvalues scale as ~1/tau_d, so LLE scales proportionally.

## Changes Made

1. ✓ **Fixed `c_E = 0.15/3`** in `example_memory_capacity.m` (was 0.5)
2. ✓ **Fixed `lya_T_interval`** to skip first 15 seconds (was 2s) in SRNNModel2
3. ✓ **Internalized Lyapunov functions** in SRNNModel2 to avoid path conflicts

## Remaining Differences

- `level_of_chaos = 1.8` vs 1.0
- `tau_d = 0.025` vs 0.1
- `n = 200` vs 300
- Input type: bandlimited vs step pulses

## Files for Reference
- ConnectivityAdaptation: `StabilityAnalysis/scripts/Simple_network_with_dual_adaptation.m`
- FractionalResevoir: `scripts/example_memory_capacity.m`
- SRNNModel defaults: `src/SRNNModel.m` / `src/SRNNModel2.m`
