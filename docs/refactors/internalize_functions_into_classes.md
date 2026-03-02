# Refactor: Internalize External Functions into SRNNModel2 and SRNN_ESN_reservoir

**Date:** 2026-02-27  
**Branch:** `speech_recognition`  
**Commits:** `6e2c58f` (SRNNModel2), `49945cb` (SRNN_ESN_reservoir)

## Summary

Moved standalone external functions into `SRNNModel2.m` and `SRNN_ESN_reservoir.m` as static methods to make the classes self-contained. After this refactor, building and running these classes no longer requires `src/nonlinearities/`, `src/generate_stimulus/`, `src/algorithms/Jacobian/`, `src/plotting/`, or `src/verify_shared_build.m` on the MATLAB path. The only external dependency remaining is `RMTMatrix.m` (a separate class for random matrix theory connectivity generation) and `setup_paths.m`.

## Motivation

- **Reduce path fragility** — classes previously failed silently or with cryptic errors when `setup_paths.m` was not called or when specific subdirectories were missing from the path.
- **Make classes standalone** — `SRNNModel2` and `SRNN_ESN_reservoir` can now be used by copying only the class files and `RMTMatrix.m`, without needing the full `src/` tree.
- **Enable subclass reuse** — subclasses that override plotting or stimulus generation can call `SRNNModel2.methodName(...)` directly instead of depending on path-resolved standalone functions.

## Functions Internalized into SRNNModel2

### Activation Functions (`methods (Static)`)

| Method                                | Original file                                     | Notes                              |
| ------------------------------------- | ------------------------------------------------- | ---------------------------------- |
| `piecewiseSigmoid(x, a, c)`           | `src/nonlinearities/piecewiseSigmoid.m`           | Piecewise linear/quadratic sigmoid |
| `piecewiseSigmoidDerivative(x, a, c)` | `src/nonlinearities/piecewiseSigmoidDerivative.m` | First derivative                   |
| `logisticSigmoid(x)`                  | `src/nonlinearities/logisticSigmoid.m`            | Standard 4x-slope logistic sigmoid |
| `logisticSigmoidDerivative(x)`        | `src/nonlinearities/logisticSigmoidDerivative.m`  | First derivative                   |
| `tanhActivation(x)`                   | `src/nonlinearities/tanhActivation.m`             | Wrapper around `tanh`              |
| `tanhActivationDerivative(x)`         | `src/nonlinearities/tanhActivationDerivative.m`   | `1 - tanh(x)^2`                    |

**Access:** Public static. These are referenced by function handle in `set_defaults()`:
```matlab
obj.activation_function = @(x) SRNNModel2.piecewiseSigmoid(x, obj.S_a, obj.S_c);
```

### Stimulus Generation (`methods (Static, Access = protected)`)

| Method                                                           | Original file                                     | Notes                       |
| ---------------------------------------------------------------- | ------------------------------------------------- | --------------------------- |
| `generate_external_input(params, T, fs, rng_seed, input_config)` | `src/generate_stimulus/generate_external_input.m` | Sparse random step stimulus |

**Access:** Protected static. Called by `build_stimulus()` in the parent class. Supports custom generator function handles via `input_config.generator`.

### Jacobian Computation (`methods (Static)`)

| Method                                                | Original file                                           | Notes                               |
| ----------------------------------------------------- | ------------------------------------------------------- | ----------------------------------- |
| `compute_Jacobian_fast(S, params)`                    | `src/algorithms/Jacobian/compute_Jacobian_fast.m`       | Sparse/vectorized Jacobian assembly |
| `compute_Jacobian_at_indices(S_out, J_times, params)` | `src/algorithms/Jacobian/compute_Jacobian_at_indices.m` | Multi-timepoint Jacobian wrapper    |

**Access:** Public static. `compute_Jacobian_fast` is referenced by the QR Lyapunov method via:
```matlab
jacobian_wrapper = @(tt, S, p) SRNNModel2.compute_Jacobian_fast(S, p);
```

### Helper (`methods (Static, Access = private)`)

| Method                                         | Original file                 | Notes                                                        |
| ---------------------------------------------- | ----------------------------- | ------------------------------------------------------------ |
| `safe_get_param(params, field, default_value)` | New (no external predecessor) | Extracts field with default; used by `compute_Jacobian_fast` |

### Plotting Functions (`methods (Static, Access = protected)`)

| Method                                                 | Original file                             | Notes                               |
| ------------------------------------------------------ | ----------------------------------------- | ----------------------------------- |
| `plot_SRNN_tseries(...)`                               | `src/plotting/plot_SRNN_tseries.m`        | Main multi-panel time series figure |
| `plot_external_input(t, u)`                            | `src/plotting/plot_external_input.m`      | E/I input visualization             |
| `plot_dendritic_state(t, x, plot_mean)`                | `src/plotting/plot_dendritic_state.m`     | Dendritic potential dynamics        |
| `plot_firing_rate(t, r)`                               | `src/plotting/plot_firing_rate.m`         | Neural activity patterns            |
| `plot_synaptic_output(t, br)`                          | `src/plotting/plot_synaptic_output.m`     | `b .* r` product                    |
| `plot_adaptation(t, a, params)`                        | `src/plotting/plot_adaptation.m`          | Adaptation variable time courses    |
| `plot_std_variable(t, b, params)`                      | `src/plotting/plot_std_variable.m`        | Synaptic depression dynamics        |
| `plot_lyapunov(lya_results, Lya_method, plot_options)` | `src/plotting/plot_lyapunov.m`            | Lyapunov exponent visualization     |
| `plot_lines_with_colormap(t, data, cmap, ...)`         | `src/plotting/plot_lines_with_colormap.m` | Shared line-plotting utility        |
| `inhibitory_colormap(n_colors)`                        | `src/plotting/inhibitory_colormap.m`      | Reds/magentas colormap              |
| `excitatory_colormap(n_colors)`                        | `src/plotting/excitatory_colormap.m`      | Blues/greens colormap               |

**Access:** Protected static. Called by `plot()` and `decimate_and_unpack()` within the class hierarchy.

### Enhanced Plotting: `plot_eigenvalues_helper`

`plot_eigenvalues_helper` was already a static protected method before this refactor (from an earlier internalization). In this commit, its implementation was upgraded from a simple 2-tier coloring (inside/outside circle) to a **3-tier outlier coloring** scheme:

| Tier         | Condition                              | Visual                 |
| ------------ | -------------------------------------- | ---------------------- |
| Interior     | `distance ≤ R`                         | Black unfilled circles |
| Near outlier | `R < distance ≤ outlier_threshold × R` | Black X markers        |
| Far outlier  | `distance > outlier_threshold × R`     | Green filled circles   |

The `outlier_threshold` parameter (default: 1.04) is read from `circle_params.outlier_threshold`.

### Lyapunov Functions (previously internalized)

These were internalized from `ConnectivityAdaptation` in an earlier refactor and are listed here for completeness. They were **not changed** in this commit but are part of the same static method block:

| Method                                            | Original file                                                                 | Access           |
| ------------------------------------------------- | ----------------------------------------------------------------------------- | ---------------- |
| `compute_lyapunov_exponents_internal(...)`        | `ConnectivityAdaptation/src/algorithms/Lyapunov/compute_lyapunov_exponents.m` | Static protected |
| `benettin_algorithm_internal(...)`                | `ConnectivityAdaptation/src/algorithms/Lyapunov/benettin_algorithm.m`         | Static protected |
| `lyapunov_spectrum_qr_internal(...)`              | `ConnectivityAdaptation/src/algorithms/Lyapunov/lyapunov_spectrum_qr.m`       | Static protected |
| `variational_eqs_ode_internal(...)`               | (inline helper for QR method)                                                 | Static protected |
| `get_minMaxRange_internal(params)`                | (inline helper for Benettin clamping)                                         | Static protected |
| `compute_kaplan_yorke_dimension_internal(lambda)` | (inline helper)                                                               | Static protected |

---

## Functions Internalized into SRNN_ESN_reservoir

| Method                                                                     | Original file               | Access        | Notes                              |
| -------------------------------------------------------------------------- | --------------------------- | ------------- | ---------------------------------- |
| `verify_shared_build(esn_array, expected_to_differ, also_check_protected)` | `src/verify_shared_build.m` | Static public | Comparison-experiment verification |

### Calling convention change

Before:
```matlab
verify_shared_build(esn_array, expected_to_differ, also_check_protected)
```

After:
```matlab
SRNN_ESN_reservoir.verify_shared_build(esn_array, expected_to_differ, also_check_protected)
```

---

## What Was NOT Internalized

| Item                                 | Reason                                                                                                              |
| ------------------------------------ | ------------------------------------------------------------------------------------------------------------------- |
| `RMTMatrix.m`                        | Separate class with its own property/method structure. Used by `build_network()`. Shared with other projects.       |
| `setup_paths.m`                      | Path management utility; still needed for `RMTMatrix` and other shared tools in `src/`.                             |
| `src/connectivity/create_W_matrix.m` | Legacy function not used by SRNNModel2 (uses `RMTMatrix` instead).                                                  |
| `src/unpack_and_compute_states.m`    | Already internalized as `unpack_and_compute_states()` instance method in the earlier build-decomposition refactor.  |
| `src/initialize_state.m`             | Already internalized as `initialize_state()` protected instance method in the earlier build-decomposition refactor. |

---

## Impact on External Scripts

Scripts that previously called standalone functions like `piecewiseSigmoid(x, a, c)` or `compute_Jacobian_fast(S, params)` directly can now call them via the class:

```matlab
% Before (requires src/nonlinearities/ on path)
y = piecewiseSigmoid(x, a, c);

% After (no path dependency beyond the class file)
y = SRNNModel2.piecewiseSigmoid(x, a, c);
```

For protected methods (e.g., plotting, stimulus generation), external scripts cannot call them directly — they must go through the class's public API (`plot()`, `build()`, etc.). This is intentional: these methods are implementation details.

`setup_paths.m` is still needed for `RMTMatrix.m` and any other shared utilities not yet internalized.
