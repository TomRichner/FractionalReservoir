# Stability Analysis Code Structure

This document describes the object-oriented architecture of the FractionalReservoir codebase, focusing on three core classes that implement simulations of recurrent neural networks with spike-frequency adaptation (**SFA**) and short-term synaptic depression (**STD**).

> **Note:** This document describes the `...2.m` versions of the classes (`SRNNModel2`, `ParamSpaceAnalysis2`) which are the current implementations in this repository. These were forked from the original `ConnectivityAdaptation/StabilityAnalysis` codebase and have diverged: several standalone utility functions have been internalized as class methods, and new properties and methods have been added.

---

## System Equations

The SRNN dynamics are governed by the following system of differential equations (see [system_equations.md](../EquationsParametersDocs/system_equations.md)):

$$
\frac{dx_i}{dt} = \frac{-x_i + u_i + \sum_{j=1}^{N} w_{ij}\, b_j r_{j}}{\tau_d}
$$

$$
r_i = \phi\left( x_i - a_{0_i} - c \sum_{k=1}^{K} a_{ik} \right)
$$

$$
\frac{da_{ik}}{dt} = \frac{-a_{ik} + r_i}{\tau_{a_k}}
$$

$$
\frac{db_i}{dt} = \frac{1-b_i}{\tau_{rec}} - \frac{b_i\, r_i}{\tau_{rel}}
$$

| Symbol                   | Description                                        |
| ------------------------ | -------------------------------------------------- |
| $x_i$                    | Dendritic potential of neuron $i$                  |
| $r_i$                    | Firing rate (output of activation function $\phi$) |
| $a_{ik}$                 | SFA variable for neuron $i$ at timescale $k$       |
| $b_i$                    | STD variable (synaptic resources available)        |
| $w_{ij}$                 | Connection weight from neuron $j$ to neuron $i$    |
| $u_i$                    | External input to neuron $i$                       |
| $\tau_d$                 | Dendritic time constant                            |
| $\tau_{a_k}$             | SFA time constant for timescale $k$                |
| $\tau_{rec}, \tau_{rel}$ | STD recovery and release time constants            |
| $c$                      | SFA coupling strength                              |

For a complete parameter reference, see [parameter_table.md](../EquationsParametersDocs/parameter_table.md).

---

## Class: SRNNModel2

**File:** [SRNNModel2.m](../../src/SRNNModel2.m)

`SRNNModel2` is the main simulation class. It encapsulates network parameters, ODE integration, Lyapunov exponent computation, and visualization. It was forked from `SRNNModel.m` and has internalized several previously standalone functions as static/instance methods (Lyapunov algorithms, state decimation, state initialization, eigenvalue plotting).

### Usage Pattern

```matlab
model = SRNNModel2('n', 300, 'n_a_E', 3, 'n_b_E', 1, 'level_of_chaos', 1.0);
model.build();   % Initialize W, stimulus, initial conditions
model.run();     % Integrate ODEs (and compute Lyapunov exponents if lya_method != 'none')
model.plot();    % Generate time series plots
```

### Key Properties

| Category                 | Properties                                                                                  | Description                                                                                            |
| ------------------------ | ------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------ |
| **Network Architecture** | `n`, `f`, `indegree`                                                                        | Total neurons, E fraction, expected in-degree                                                          |
| **RMT Parameters**       | `mu_E_tilde`, `mu_I_tilde`, `sigma_E_tilde`, `sigma_I_tilde`                                | Tilde-notation stats (Harris 2023)                                                                     |
| **Mean Offset**          | `E_W`                                                                                       | Added to both `mu_E_tilde` and `mu_I_tilde` during W construction (default: 0)                         |
| **ZRS**                  | `zrs_mode`                                                                                  | Zero row sum mode: `'none'`, `'ZRS'`, `'SZRS'`, `'Partial_SZRS'`                                       |
| **Scaling**              | `level_of_chaos`, `rescale_by_abscissa`                                                     | W scaling factor, and optional 1/abscissa₀ normalization                                               |
| **SFA**                  | `n_a_E`, `n_a_I`, `tau_a_E`, `tau_a_I`, `c_E`, `c_I`                                        | Number of SFA timescales, time constant vectors (log-spaced defaults), coupling constants              |
| **STD**                  | `n_b_E`, `n_b_I`, `tau_b_E_rec`, `tau_b_E_rel`, etc.                                        | Enable/disable STD (0 or 1), and time constants if enabled                                             |
| **Dynamics**             | `tau_d`                                                                                     | Dendritic integration time constant                                                                    |
| **Nonlinearity**         | `activation_function`, `activation_function_derivative`, `S_a`, `S_c`                       | Firing rate nonlinearity handle, its derivative, and shape parameters (a = slope fraction, c = center) |
| **Input**                | `input_config`, `u_ex_scale`                                                                | Struct with stimulus parameters (n_steps, step_density, amp, etc.), and external input scaling factor  |
| **Simulation**           | `fs`, `T_range`, `T_plot`, `ode_solver`, `ode_opts`                                         | Sampling freq, time interval, plot interval, solver                                                    |
| **Storage**              | `store_full_state`, `store_decimated_state`, `plot_freq`, `plot_deci`                       | Control what state data is retained after run (default: store decimated only)                          |
| **Lyapunov**             | `lya_method`, `lya_T_interval`, `filter_local_lya`, `lya_filter_order`, `lya_filter_cutoff` | Method ('benettin', 'qr', 'none'), time interval, and optional lowpass filtering of local LLE          |
| **RNG**                  | `rng_seeds`                                                                                 | `[network_seed, stimulus_seed]` for reproducibility                                                    |

### Dependent Properties (Computed)

| Property                 | Description                                     |
| ------------------------ | ----------------------------------------------- |
| `alpha`                  | Sparsity = `indegree/n`                         |
| `default_val`            | Normalization factor F = 1/√(N·α·(2−α))         |
| `mu_se`, `mu_si`         | Sparse E/I means                                |
| `sigma_se`, `sigma_si`   | Sparse E/I std devs (includes `E_W` offset)     |
| `R`                      | Theoretical spectral radius (Harris 2023 Eq 18) |
| `n_E`, `n_I`             | Number of E/I neurons                           |
| `E_indices`, `I_indices` | Index vectors for E/I populations               |
| `N_sys_eqs`              | Total state dimension                           |

### Core Methods

| Method                          | Description                                                                                                           |
| ------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| `build()`                       | Delegates to `build_network()`, `build_stimulus()`, `finalize_build()` — creates W, generates stimulus, caches params |
| `run()`                         | Integrates ODEs, computes Lyapunov exponents, decimates state for plotting                                            |
| `plot()`                        | Generates 6-panel time series figure (input, x, r, a, b, Lyapunov)                                                    |
| `plot_eigenvalues(J_times_sec)` | Plot Jacobian eigenvalues at specified time points                                                                    |
| `plot_W_spectrum()`             | 2-panel figure: eigenvalues of (−I+W) and the LTI Jacobian (−I+W)/τ_d                                                 |
| `compute_lyapunov()`            | Standalone Lyapunov computation (can be called separately from `run()`)                                               |
| `filter_lyapunov()`             | Apply Butterworth lowpass to local Lyapunov time series                                                               |
| `get_params()`                  | Returns a struct for compatibility with standalone functions                                                          |
| `clear_results()`               | Free memory by clearing stored state data                                                                             |
| `reset()`                       | Clear built state to allow rebuilding with new parameters                                                             |
| `dynamics(t, S)`                | Convenience wrapper around static `dynamics_fast`                                                                     |

### Protected Build Sub-Methods

`build()` delegates to three protected sub-methods that subclasses (e.g., `SRNN_ESN_reservoir`) can override:

| Method             | Description                                                                                                 |
| ------------------ | ----------------------------------------------------------------------------------------------------------- |
| `build_network()`  | Sets RNG seed, fills in default tilde/tau parameters, creates and scales W via `RMTMatrix`                  |
| `build_stimulus()` | Generates external stimulus using `generate_external_input`, builds `griddedInterpolant`, initializes state |
| `finalize_build()` | Validates configuration, caches params struct                                                               |

### State Vector Organization

The state vector `S` is packed as:

```
S = [a_E(:); a_I(:); b_E(:); b_I(:); x(:)]
```

where:
- `a_E`: E adaptation variables (n_E × n_a_E, flattened)
- `a_I`: I adaptation variables (n_I × n_a_I, flattened)
- `b_E`: E STD variables (n_E × n_b_E, flattened)
- `b_I`: I STD variables (n_I × n_b_I, flattened)
- `x`: Dendritic states for all neurons (n × 1)

### Equation Implementation in `dynamics_fast()`

The static method `dynamics_fast(t, S, params)` implements the ODEs. Key steps:

1. **Interpolate external input**: `u = params.u_interpolant(t)'`
2. **Unpack state**: Extract `a_E`, `a_I`, `b_E`, `b_I`, `x` from `S`
3. **Compute effective potential**: `x_eff = x - c_E * sum(a_E, 2)` for E neurons (and similarly for I)
4. **Apply STD**: `b = ones(n,1)` with `b(E_indices) = b_E` if enabled
5. **Compute firing rate**: `r = activation_fn(x_eff)`
6. **Compute derivatives**:
   - `dx_dt = (-x + W * (b .* r) + u) / tau_d`
   - `da_E_dt = (r(E_indices) - a_E) ./ tau_a_E`
   - `db_E_dt = (1 - b_E) / tau_b_E_rec - (b_E .* r(E_indices)) / tau_b_E_rel`
7. **Pack derivatives**: `dS_dt = [da_E_dt(:); da_I_dt(:); db_E_dt(:); db_I_dt(:); dx_dt]`

### Internalized Methods

The following functions exist both as standalone `.m` files (in `src/algorithms/`) and as internal static/instance methods within `SRNNModel2`. The internalized versions are used by `SRNNModel2` to avoid path conflicts with the original codebase:

| Internal Method                                    | Standalone Equivalent                              | Description                            |
| -------------------------------------------------- | -------------------------------------------------- | -------------------------------------- |
| `compute_lyapunov_exponents_internal` (static)     | `algorithms/Lyapunov/compute_lyapunov_exponents.m` | Main Lyapunov dispatcher               |
| `benettin_algorithm_internal` (static)             | `algorithms/Lyapunov/benettin_algorithm.m`         | Benettin's LLE algorithm               |
| `lyapunov_spectrum_qr_internal` (static)           | `algorithms/Lyapunov/lyapunov_spectrum_qr.m`       | QR full spectrum method                |
| `get_minMaxRange_internal` (static)                | `algorithms/Jacobian/get_minMaxRange.m`            | State bounds for Benettin              |
| `compute_kaplan_yorke_dimension_internal` (static) | —                                                  | Kaplan-Yorke dimension                 |
| `initialize_state` (instance)                      | `initialize_state.m`                               | State vector initialization            |
| `decimate_states` (instance)                       | `decimate_states.m`                                | State decimation for plotting          |
| `plot_eigenvalues_helper` (static)                 | `plotting/plot_eigenvalues.m`                      | Eigenvalue scatter plot                |
| `filter_lyapunov` (instance)                       | —                                                  | Local Lyapunov lowpass filtering (new) |

### Lyapunov Methods

The `lya_method` property controls how Lyapunov exponents are computed:

| Method                 | Output        | Description                                                                                                                      | Time Complexity            | Memory                      |
| ---------------------- | ------------- | -------------------------------------------------------------------------------------------------------------------------------- | -------------------------- | --------------------------- |
| `'benettin'` (default) | LLE only      | Computes local and finite-time largest Lyapunov exponent using a shadow trace (single perturbation vector rescaled periodically) | O(T · K_steps · N_states)  | O(N_states)                 |
| `'qr'`                 | Full spectrum | Computes all N_states Lyapunov exponents via QR decomposition of the tangent space; integrates N_states² variational equations   | O(T · K_steps · N_states³) | O(N_states² + T · N_states) |
| `'none'`               | —             | Skips Lyapunov computation entirely                                                                                              | —                          | —                           |

where:
- **N_states** = system dimension (`N_sys_eqs`), not to be confused with network size N
- **T** = number of rescaling intervals = `(T_end - T_start) / lya_dt_interval`
- **K_steps** = average ODE substeps per rescaling interval (depends on solver tolerances and dynamics stiffness)

> **Recommendation:** Use `'benettin'` for large networks (N_states > 100) or parameter sweeps. Use `'qr'` only when the full spectrum is needed and N_states is small (the N_states³ scaling makes it prohibitive for N_states > ~200).

---

## Class: RMTMatrix

**File:** [RMTMatrix.m](../../src/RMTMatrix.m)

`RMTMatrix` constructs sparse connectivity matrices following Random Matrix Theory (Harris et al., 2023). It enforces Dale's law by separating excitatory and inhibitory populations with distinct mean weights.

### Usage Pattern

```matlab
rmt = RMTMatrix(N);
rmt.alpha = indegree / N;     % Sparsity
rmt.f = 0.5;                   % E fraction
rmt.mu_tilde_e = 3 * D;        % E mean (tilde notation)
rmt.mu_tilde_i = -4 * D;       % I mean (tilde notation)
rmt.sigma_tilde_e = D;         % E std dev
rmt.sigma_tilde_i = D;         % I std dev
W = rmt.W;                     % Access triggers construction
```

### Key Properties

| Property                         | Description                                                      |
| -------------------------------- | ---------------------------------------------------------------- |
| `N`                              | Network size                                                     |
| `alpha`                          | Connection probability (sparsity)                                |
| `f`                              | Fraction of excitatory neurons                                   |
| `mu_tilde_e`, `mu_tilde_i`       | Normalized population means                                      |
| `sigma_tilde_e`, `sigma_tilde_i` | Normalized population std devs                                   |
| `zrs_mode`                       | Zero row sum mode: `'none'`, `'ZRS'`, `'SZRS'`, `'Partial_SZRS'` |

### Dependent Properties (Computed)

| Property   | Formula                                           | Description                         |
| ---------- | ------------------------------------------------- | ----------------------------------- |
| `mu_se`    | $\alpha \cdot \tilde{\mu}_E$                      | Sparse excitatory mean              |
| `mu_si`    | $\alpha \cdot \tilde{\mu}_I$                      | Sparse inhibitory mean              |
| `R`        | $\sqrt{N(f \sigma_{se}^2 + (1-f) \sigma_{si}^2)}$ | Theoretical spectral radius (Eq 18) |
| `lambda_O` | $N(f \mu_{se} + (1-f) \mu_{si})$                  | Outlier eigenvalue (Eq 17)          |

### W Matrix Construction

The weight matrix W is constructed as:

```matlab
W_dense = (A * D) + M;    % A: Gaussian random, D: variance diagonal, M: low-rank mean
W = S .* W_dense;         % S: sparse binary mask
```

where:
- `A` is an N × N Gaussian random matrix (mean 0, variance 1)
- `D = diag([sigma_tilde_e, ..., sigma_tilde_e, sigma_tilde_i, ..., sigma_tilde_i])` encodes population variance
- `M = u * v'` is the rank-1 mean structure with `v = [mu_tilde_e, ..., mu_tilde_e, mu_tilde_i, ..., mu_tilde_i]'`
- `S` is a Bernoulli (0 or 1) sparsity mask with connection probability `alpha`

---

## Class: ParamSpaceAnalysis2

**File:** [ParamSpaceAnalysis2.m](../../src/ParamSpaceAnalysis2.m)

`ParamSpaceAnalysis2` performs multi-dimensional gridded parameter space analysis of stability (and other metrics). The user can choose any combination of model parameters to vary. By default it simulates all parameter combinations across four conditions: no adaptation, SFA only, STD only, and SFA + STD. This enables a systematic comparison of network behavior under different parameter ranges and conditions. The parallel computing toolbox is highly recommended but not required.

### Usage Pattern

```matlab
psa = ParamSpaceAnalysis2('n_levels', 5, 'note', 'f_sweep');
psa.add_grid_parameter('f', [0.4, 0.6]);       % Sweep E fraction (range mode)
psa.add_grid_parameter('reps', [1:10]);         % 10 repetitions per point (explicit mode)
psa.model_defaults.n = 300;                    % Set constant parameters
psa.run();                                     % Execute analysis
psa.plot('metric', 'LLE');                     % Visualize results
```

### `add_grid_parameter` Behavior

The `add_grid_parameter(param_name, values)` method supports two modes:

| Mode         | Input                   | Behavior                                                                             |
| ------------ | ----------------------- | ------------------------------------------------------------------------------------ |
| **Range**    | 1×2 vector `[min, max]` | Evenly divides the range into `n_levels` values using `linspace(min, max, n_levels)` |
| **Explicit** | Vector with 3+ elements | Uses the exact values provided (ignores `n_levels` for this parameter)               |

**Examples:**
```matlab
psa.add_grid_parameter('f', [0.4, 0.6]);       % Range mode: 5 levels → [0.4, 0.45, 0.5, 0.55, 0.6]
psa.add_grid_parameter('reps', [1, 2, 3, 4, 5]); % Explicit mode: uses [1, 2, 3, 4, 5] directly
```

> **Note:** Each call to `add_grid_parameter` adds a dimension to the parameter space. The total number of simulations grows **multiplicatively**:
> 
> `Total = n_conditions × n_levels_1 × n_levels_2 × ...`
> 
> For example, 3 grid parameters with 5 levels each and 4 conditions yields:
> `4 × 5^3 = 500` total simulations.

### Key Properties

| Category       | Property                                          | Description                                                                       |
| -------------- | ------------------------------------------------- | --------------------------------------------------------------------------------- |
| **Grid**       | `grid_params`, `param_ranges`, `explicit_vectors` | Parameter names, ranges for range mode, explicit vectors for explicit mode        |
| **Grid**       | `n_levels`                                        | Number of levels per parameter (range mode only)                                  |
| **Grid**       | `integer_params`                                  | Parameter names to round to integers (default: `{'n', 'indegree', 'n_a_E', ...}`) |
| **Conditions** | `conditions`                                      | Cell array of condition structs (`name`, `n_a_E`, `n_b_E`)                        |
| **Model**      | `model_defaults`                                  | Struct of default `SRNNModel2` properties applied to every simulation             |
| **Execution**  | `batch_size`, `use_parallel`                      | Configs per batch for checkpointing; toggle `parfor` vs `for`                     |
| **Output**     | `output_dir`, `note`                              | Base directory for saving results; optional note for folder naming                |
| **Storage**    | `store_local_lya`, `store_local_lya_dt`           | Whether to store decimated local Lyapunov time series and at what resolution      |
| **Verbosity**  | `verbose`                                         | Print progress during execution                                                   |

### Core Methods

| Method                            | Description                                                             |
| --------------------------------- | ----------------------------------------------------------------------- |
| `add_grid_parameter(name, range)` | Add a parameter to the grid (range or explicit mode)                    |
| `remove_grid_parameter(name)`     | Remove a parameter from the grid                                        |
| `set_conditions(cell_array)`      | Set custom conditions (replaces default four)                           |
| `run()`                           | Generate grid, randomize order, run batched parfor, consolidate results |
| `plot('metric', 'LLE')`           | Generate histogram plots of a metric across parameter space             |
| `load_results(results_dir)`       | Load results from a previous run                                        |
| `consolidate()`                   | Standalone recovery: merge batch results after an interrupted run       |
| `saveobj()` / `loadobj()`         | MATLAB serialization support for `save()`/`load()`                      |

### Key Features

- **Multi-dimensional grid**: All combinations of specified parameters are tested
- **Four adaptation conditions** (default):
  | Condition       | n_a_E | n_b_E | Description  |
  | --------------- | ----- | ----- | ------------ |
  | `no_adaptation` | 0     | 0     | Baseline     |
  | `sfa_only`      | 3     | 0     | SFA enabled  |
  | `std_only`      | 0     | 1     | STD enabled  |
  | `sfa_and_std`   | 3     | 1     | Both enabled |
- **Same W across conditions**: Each grid point uses identical connectivity for fair comparison
- **Batched parallel execution**: Uses `parfor` with checkpointing for resume capability
- **Randomized order**: Allows representative early stopping
- **Serialization**: `saveobj`/`loadobj` enable saving and restoring `ParamSpaceAnalysis2` objects

### Output Structure

Results are saved to `data/param_space/<timestamped_folder>/`:

```
param_space_<note>_nLevs_<N>_<timestamp>/
|-- param_space_summary.mat      # Grid configuration
|-- psa_object.mat               # Serialized PSA object
|-- figures/
|   +-- LLE_distribution.png
|   +-- mean_rate_distribution.png
|-- no_adaptation/
|   +-- param_space_results_no_adaptation.mat
|-- sfa_only/
|   +-- param_space_results_sfa_only.mat
|-- std_only/
|   +-- param_space_results_std_only.mat
+-- sfa_and_std/
    +-- param_space_results_sfa_and_std.mat
```

Each result struct contains: `LLE`, `mean_rate`, `mean_synaptic_output`, `config`, `config_idx`, `network_seed`, `success`, `run_duration`, and optionally `local_lya` and `t_lya` time series.

---

## Class: SRNN_ESN_reservoir

**File:** [SRNN_ESN_reservoir.m](../../src/SRNN_ESN_reservoir.m)

`SRNN_ESN_reservoir` extends `SRNNModel2` to provide Echo State Network (ESN) functionality with memory capacity measurement. It overrides `build_stimulus()` to generate ESN-specific inputs (scalar random sequences with configurable spectral properties) and adds read-out training.

### Usage Pattern

```matlab
esn = SRNN_ESN_reservoir('n', 300, 'n_a_E', 3, 'n_b_E', 1, ...
    'T_wash', 4000, 'T_train', 10000, 'T_test', 10000, 'd_max', 600);
esn.build();
[MC, R2_d, results] = esn.run_memory_capacity();
esn.plot_esn_timeseries(delays_to_plot);
```

### Key Additional Properties

| Category        | Property                               | Description                                                                   |
| --------------- | -------------------------------------- | ----------------------------------------------------------------------------- |
| **Input**       | `f_in`, `sigma_in`, `rng_seed_input`   | Fraction of neurons receiving input, weight scale, seed                       |
| **Input Type**  | `input_type`, `u_f_cutoff`, `u_alpha`  | `'white'`, `'bandlimited'`, or `'one_over_f'`; cutoff freq; spectral exponent |
| **MC Protocol** | `T_wash`, `T_train`, `T_test`, `d_max` | Washout/train/test durations (samples), max delay                             |

### Core Methods

| Method                        | Description                                                            |
| ----------------------------- | ---------------------------------------------------------------------- |
| `run_memory_capacity()`       | Full MC protocol: drive, train readouts, compute R²_d, return total MC |
| `run_reservoir_esn()`         | Run reservoir ODE integration only (no readout training)               |
| `plot_memory_capacity()`      | R² vs delay and cumulative MC plots                                    |
| `plot_esn_timeseries(delays)` | Time series with target/prediction overlays for selected delays        |

---

## Supporting Functions

The `src/` directory contains additional standalone utilities organized into subdirectories:

| Category           | Location               | Files                                                                                                                | Description                                                            |
| ------------------ | ---------------------- | -------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------- |
| **Lyapunov**       | `algorithms/Lyapunov/` | `compute_lyapunov_exponents.m`, `lyapunov_spectrum_qr.m`, `benettin_algorithm.m`                                     | Standalone LLE and full spectrum computation                           |
| **Jacobian**       | `algorithms/Jacobian/` | `compute_Jacobian_fast.m`, `compute_Jacobian_at_indices.m`, `compute_J_eff.m`, `get_minMaxRange.m`                   | Sparse Jacobian construction                                           |
| **Stimulus**       | `generate_stimulus/`   | `generate_external_input.m`, `generate_paired_pulse_input.m`, `generate_AM_pulse_train.m`, `generate_mackey_glass.m` | Input generation utilities                                             |
| **Nonlinearities** | `nonlinearities/`      | `piecewiseSigmoid.m`, `logisticSigmoid.m`, `tanhActivation.m` (+ derivatives)                                        | Activation functions and their derivatives                             |
| **Plotting**       | `plotting/`            | `plot_SRNN_tseries.m`, `plot_SRNN_combined_tseries.m`, colormaps, etc.                                               | Visualization utilities                                                |
| **Verification**   | `src/`                 | `verify_shared_build.m`                                                                                              | Asserts that multiple ESN conditions share identical W, W_in, stimulus |
| **State**          | `src/`                 | `unpack_and_compute_states.m`                                                                                        | Unpack state vector into named fields                                  |

> **Note:** Many of the standalone Lyapunov, Jacobian, and state management functions also exist as internalized methods within `SRNNModel2` (see [Internalized Methods](#internalized-methods) above). The class uses its internal copies; the standalone versions are retained for backward compatibility and use by other code.

---

## See Also

- [parameter_table.md](../EquationsParametersDocs/parameter_table.md) – Complete parameter reference
- [system_equations.md](../EquationsParametersDocs/system_equations.md) – Mathematical model
- [Script_Notes2.md](./Script_Notes2.md) – How to run the current scripts
