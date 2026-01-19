# Parameter Comparison: SRNN Example vs Parameter Space Analysis

This document compares the parameters used in the **example script** (`SRNN_comparisons_for_review_v3.m` → `full_SRNN_run_v3.m`) versus the **parameter space analysis** (`run_param_space_analysis.m` → `ParamSpaceAnalysis` → `SRNNModel`).

---

## Grid Parameters

These parameters are varied across a range in `run_param_space_analysis.m`. The example uses fixed values.

| Parameter | Example Value | Grid Range (PSA) | In Range? | Notes |
|-----------|--------------|------------------|-----------|-------|
| `level_of_chaos` | 1.7 | [0.5, 2.5] | ✅ Yes | Abscissa scaling factor |
| `n` | 100 | [50, 250] | ✅ Yes | Number of neurons |
| `f` | 0.5 | [0.2, 0.8] | ✅ Yes | Fraction of E neurons |
| `EI_imbalance` | 1 | [0.5, 2] | ✅ Yes | Example uses `imbalance = 1` which is equivalent |
| `mu_E` | 1 | [0.2, 2] | ✅ Yes | Mean excitatory weight |

**Summary**: All example values fall within the grid parameter ranges. ✅

---

## Adaptation Condition Parameters

Both scripts test 4 conditions with matching settings:

| Condition | `n_a_E` | `n_b_E` | Example? | PSA? |
|-----------|---------|---------|----------|------|
| No Adaptation | 0 | 0 | ✅ Run 1 | ✅ |
| SFA Only | 3 | 0 | ❌ | ✅ |
| STD Only | 0 | 1 | ❌ | ✅ |
| SFA + STD | 3 | 1 | ✅ Run 4 | ✅ |

**Note**: Example only runs 2 of the 4 conditions (no adaptation and SFA+STD).

---

## Non-Grid Parameters

### Network Architecture

| Parameter | Example Value | SRNNModel Default | PSA Setting | Match? |
|-----------|--------------|-------------------|-------------|--------|
| `indegree` | 20 | 20 | — | ✅ |
| `G_stdev` | `1/sqrt(indegree)` | computed same way | — | ✅ |
| `d` (density) | `indegree/n` | computed same way | — | ✅ |

### Spike-Frequency Adaptation (SFA)

| Parameter | Example Value | SRNNModel Default | PSA Setting | Match? |
|-----------|--------------|-------------------|-------------|--------|
| `n_a_E` | 0 or 3 (condition) | 0 | per condition | ✅ |
| `n_a_I` | 0 | 0 | — | ✅ |
| `tau_a_E` | `logspace(log10(0.5), log10(10), n_a_E)` | `logspace(log10(0.25), log10(10), n_a_E)` | — | ⚠️ **Mismatch** |
| `tau_a_I` | `logspace(log10(0.5), log10(10), n_a_I)` | `logspace(log10(0.25), log10(10), n_a_I)` | — | ⚠️ **Mismatch** |
| `c_E` | 0.15/3 ≈ **0.05** | 0.1/3 ≈ 0.033 | `0.5` (typo: `pas.`) | ⚠️ **Mismatch** |
| `c_I` | 0.1 | 0.1 | — | ✅ |

### Short-Term Depression (STD)

| Parameter | Example Value | SRNNModel Default | PSA Setting | Match? |
|-----------|--------------|-------------------|-------------|--------|
| `n_b_E` | 0 or 1 (condition) | 0 | per condition | ✅ |
| `n_b_I` | 0 | 0 | — | ✅ |
| `tau_b_E_rel` | 0.25 | 0.25 | — | ✅ |
| `tau_b_I_rel` | 0.25 | 0.25 | — | ✅ |
| `tau_b_E_rec` | **2** | **1** | — | ⚠️ **Mismatch** |
| `tau_b_I_rec` | **2** | **1** | — | ⚠️ **Mismatch** |

### Dynamics

| Parameter | Example Value | SRNNModel Default | PSA Setting | Match? |
|-----------|--------------|-------------------|-------------|--------|
| `tau_d` | 0.1 | 0.1 | — | ✅ |
| `S_a` (activation param) | 0.9 | 0.9 | — | ✅ |
| `S_c` (activation center) | **0.3** | **0.35** | — | ⚠️ **Mismatch** |
| `activation_function` | `piecewiseSigmoid(x, 0.9, 0.3)` | `piecewiseSigmoid(x, 0.9, 0.35)` | — | ⚠️ **Mismatch** |

### Simulation Settings

| Parameter | Example Value | SRNNModel Default | PSA Setting | Match? |
|-----------|--------------|-------------------|-------------|--------|
| `fs` | 200 | 200 | 200 | ✅ |
| `T_range` | **[-20, 40]** | [0, 50] | **[0, 30]** | ⚠️ **Mismatch** |
| `T_plot` | [0, 40] | T_range | — | Different |
| `ode_solver` | @ode45 | @ode45 | — | ✅ |
| `RelTol` | 1e-7 | 1e-7 | — | ✅ |
| `AbsTol` | 1e-7 | 1e-7 | — | ✅ |
| `lya_method` | 'benettin' | 'benettin' | 'benettin' | ✅ |

### External Input Configuration

| Parameter | Example Value | SRNNModel Default | Match? |
|-----------|--------------|-------------------|--------|
| `u_ex_scale` | 1.5 | 1.0 | ⚠️ **Mismatch** |
| `input_config.n_steps` | 3 | 3 | ✅ |
| `input_config.step_density` | 0.2 | 0.2 | ✅ |
| `input_config.amp` | 0.5 | 0.5 | ✅ |
| `input_config.no_stim_pattern` | `[true, false, true]` | `[true, false, true]` | ✅ |
| `input_config.intrinsic_drive` | 0 | 0 | ✅ |

### Random Seeds

| Parameter | Example Value | SRNNModel Default | PSA Setting | Notes |
|-----------|--------------|-------------------|-------------|-------|
| `rng_seeds` | [105, 25] | [1, 2] | `[config_idx, config_idx+1]` | Different seed strategies |

---

## Structural Differences

### 1. Implementation Approach
| Aspect | Example | Parameter Space Analysis |
|--------|---------|-------------------------|
| Architecture | Procedural function | Object-oriented (`SRNNModel` class) |
| Parameter passing | Hardcoded in function | Configurable via constructor args |
| Network generation | In `full_SRNN_run_v3.m` | In `SRNNModel.build()` |

### 2. Simulation Strategy
| Aspect | Example | Parameter Space Analysis |
|--------|---------|-------------------------|
| Runs per execution | 2 conditions | All 4 conditions |
| Network reuse | New W per run | Same W across all 4 conditions |
| Output | Figures + workspace | Metrics (LLE, mean_rate) only |
| Plotting | Full time series + eigenvalues | Histograms of metrics |

### 3. Typo in run_param_space_analysis.m
Line 57 has: `pas.model_defaults.c_E = 0.5;`  
Should be: `psa.model_defaults.c_E = 0.5;`

This typo means `c_E = 0.5` is **not actually applied** to the PSA runs.

---

## Summary of Mismatches

| Parameter | Example | PSA/SRNNModel | Impact |
|-----------|---------|---------------|--------|
| `tau_a_E` lower bound | 0.5 s | 0.25 s | Faster adaptation possible in PSA |
| `tau_b_E_rec` | 2 s | 1 s | Faster STD recovery in PSA |
| `tau_b_I_rec` | 2 s | 1 s | Faster STD recovery in PSA |
| `c_E` | 0.05 | 0.033 (default) | Slightly stronger SFA in example |
| `S_c` | 0.3 | 0.35 | Different activation threshold |
| `u_ex_scale` | 1.5 | 1.0 | 50% stronger input in example |
| `T_range` | [-20, 40] | [0, 30] | Example has settling time; different duration |

---

## Recommendations

1. **Update `SRNNModel` defaults** to match example values if reproducing example results is desired:
   - `tau_b_E_rec = 2`, `tau_b_I_rec = 2`
   - `S_c = 0.3`
   - `tau_a_E` lower bound = 0.5

2. **Fix typo** in `run_param_space_analysis.m` line 57: `pas.` → `psa.`

3. **Consider adding** `u_ex_scale` to grid parameters or model defaults in PSA if external input strength matters.

4. **Document** whether the T_range difference (settling time) affects Lyapunov estimates.
