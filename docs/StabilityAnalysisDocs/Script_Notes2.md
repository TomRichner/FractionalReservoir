# Script Notes: Current Scripts

This document explains how to run the current scripts in the FractionalReservoir project.

> **Note:** Many scripts in the `scripts/` directory are from earlier development and reference the original (non-v2) classes. Only the scripts documented here are up to date with `SRNNModel2` and `ParamSpaceAnalysis2`. The others should not be used without refactoring first.

---

## Prerequisites

Once per MATLAB session, with the cwd at the repository root, call [setup_paths.m](../../setup_paths.m):

```matlab
setup_paths();
```

This function:
- Locates the repository root as its own folder
- Adds `src/` and `scripts/` (and all subdirectories) to the MATLAB path
- Errors if the `src/` directory is not found

> **Note:** Entry-point scripts call `setup_paths()` themselves, so they can be launched cold. Smaller test and example scripts carry no path code and assume the session has already been bootstrapped.

---

## Script 1: Parameter Space Analysis

**File:** [run_param_space_analysis2.m](../../scripts/run_param_space_analysis2.m)

**Purpose:** Multi-dimensional parameter space exploration using `ParamSpaceAnalysis2`, comparing all four adaptation conditions on networks with the same connectivity matrix W.

### Configuration

The script configures the analysis in several steps:

**1. Create the PSA object:**
```matlab
psa = ParamSpaceAnalysis2(...
    'n_levels', 3, ...
    'batch_size', 25, ...
    'note', 'test_refactor', ...
    'verbose', true);
```

**2. Add grid parameters** (all combinations tested):
```matlab
psa.add_grid_parameter('E_W', [-0.2, 0.2] ./ sqrt(100));
psa.add_grid_parameter('f', [2/3, 3/2]);
```
Additional parameters can be uncommented (`tau_d`, `u_ex_scale`, `c_E`).

**3. Set model defaults** (constant across all runs):
```matlab
psa.model_defaults.n = 300;
psa.model_defaults.T_range = [-10, 20];
psa.model_defaults.fs = 200;
psa.model_defaults.c_E = 0.15/3;
psa.model_defaults.lya_method = 'benettin';
psa.model_defaults.level_of_chaos = 1.5;
psa.model_defaults.activation_function = @logisticSigmoid;
psa.model_defaults.activation_function_derivative = @logisticSigmoidDerivative;
% ... etc
```

**4. (Optional) Set custom conditions:**
```matlab
psa.set_conditions({
    struct('name', 'no_adaptation', 'n_a_E', 0, 'n_b_E', 0),
    struct('name', 'sfa_and_std',   'n_a_E', 3, 'n_b_E', 1)
});
```
By default, all four conditions are used (no_adaptation, sfa_only, std_only, sfa_and_std).

### Execution

```matlab
psa.run();
```

This will:
1. Generate all parameter combinations
2. Randomize execution order
3. Run batched `parfor` with checkpoint files (resume-capable)
4. Consolidate results into per-condition MAT files
5. Copy the script into the output directory for reproducibility

### Outputs

**Data:**
```
data/param_space/param_space_<note>_nLevs_<N>_<timestamp>/
|-- param_space_summary.mat
|-- psa_object.mat
|-- figures/
|   +-- LLE_distribution.{png,fig}
|   +-- mean_rate_distribution.{png,fig}
|-- no_adaptation/
|   +-- param_space_results_no_adaptation.mat
|-- sfa_only/
|   +-- param_space_results_sfa_only.mat
|-- std_only/
|   +-- param_space_results_std_only.mat
+-- sfa_and_std/
    +-- param_space_results_sfa_and_std.mat
```

**Plots:** Histograms of LLE and mean firing rate distributions for each condition.

### Accessing Results Programmatically

```matlab
% After running:
results = psa.results.sfa_only;
LLEs = cellfun(@(r) r.LLE, results(~cellfun(@isempty, results)));
histogram(LLEs);
```

### Loading Previous Results

```matlab
psa_loaded = ParamSpaceAnalysis2();
psa_loaded.load_results('/path/to/param_space_...');
psa_loaded.plot('metric', 'LLE');
```

---

## Script 2: Memory Capacity Example

**File:** [example_memory_capacity.m](../../scripts/example_memory_capacity.m)

**Purpose:** Demonstrates how to use `SRNN_ESN_reservoir` to measure memory capacity under three adaptation conditions: baseline (no adaptation), SFA only, and SFA + STD.

### Configuration

The script uses a **base config + condition-specific overrides** pattern to avoid duplicating configuration:

**Base configuration** (shared across all conditions):
```matlab
base_args = { ...
    'n', 300, ...
    'fs', 200, ...
    'level_of_chaos', 1.0, ...
    'rng_seeds', [42, 43], ...
    'tau_d', 0.1, ...
    'S_c', 0.4, 'S_a', 0.9, ...
    'c_E', 0.15/3, ...
    'tau_a_E', [0.1, 1.0, 10], ...
    'tau_b_E_rec', 1.0, 'tau_b_E_rel', 0.25, ...
    'input_type', 'white', ...    % Options: 'white', 'bandlimited', 'one_over_f'
    'T_wash', T_wash, 'T_train', T_train, 'T_test', T_test, 'd_max', d_max};
```

**Condition-specific overrides** (only what differs):
```matlab
condition_args = { ...
    {'n_a_E', 0, 'n_b_E', 0}, ...   % Baseline
    {'n_a_E', 3, 'n_b_E', 0}, ...   % SFA only
    {'n_a_E', 3, 'n_b_E', 1}, ...   % SFA + STD
};
```

### MC Protocol Parameters

| Parameter     | Default    | Description                                                   |
| ------------- | ---------- | ------------------------------------------------------------- |
| `T_wash_sec`  | 20 s       | Washout duration                                              |
| `T_train_sec` | 50 s       | Training duration                                             |
| `T_test_sec`  | 50 s       | Test duration                                                 |
| `d_max`       | 3×fs (600) | Maximum delay (samples)                                       |
| `input_type`  | `'white'`  | Input spectrum: `'white'`, `'bandlimited'`, or `'one_over_f'` |

### Execution

```matlab
% Build all conditions with shared W
for i = 1:n_cond
    esn{i} = SRNN_ESN_reservoir(base_args{:}, condition_args{i}{:});
    esn{i}.build();
end

% Verify shared structure across conditions
verify_shared_build(esn, {'n_a_E', 'n_b_E'}, {'W', 'W_in', 'u_scalar', 'u_ex', 't_ex'});

% Run memory capacity for each condition
for i = 1:n_cond
    [MC(i), R2{i}, results{i}] = esn{i}.run_memory_capacity();
end
```

The call to `verify_shared_build()` asserts that all conditions share identical network structure (W, W_in, stimulus) and that only the intended parameters (n_a_E, n_b_E) differ.

### Outputs

**Plots:**
1. **3-panel comparison figure:**
   - R² vs delay for all conditions (bar chart)
   - Cumulative MC vs delay (line plot)
   - Total MC comparison (bar chart)
2. **Per-condition time series** showing target/prediction overlays at selected delays

**Saved data:**
```
data/memory_capacity/memory_capacity_results.mat
```

Contains `MC`, `R2`, per-condition results, and condition names.

---

## Handling Interrupted PSA Runs

If `run_param_space_analysis2.m` is interrupted before completion, `ParamSpaceAnalysis2` can consolidate partial data from the completed batch checkpoint files:

```matlab
psa = ParamSpaceAnalysis2();
psa.output_dir = '/path/to/interrupted/run';
psa.consolidate();  % Merge completed batch files, load config from summary
```

> **Note:** This consolidates data from completed batches but does **not** resume computation. To obtain complete results, you must restart the analysis from scratch. Batch checkpoint files are stored in `temp_batches/` within the output directory and are cleaned up after successful consolidation.

---

## See Also

- [Stability_Analysis_Code_Structure.md](./Stability_Analysis_Code_Structure.md) – Class documentation
- [cell_type_pair_equations.md](../EquationsParametersDocs/cell_type_pair_equations.md) – Parameter reference
