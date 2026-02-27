# Refactor Summary: SRNNModel2 and SRNN_ESN_reservoir

This document describes the refactor applied to `SRNNModel2.m` and `SRNN_ESN_reservoir.m`.
Use this as a reference when applying the same patterns to other SRNNModel2 subclasses.

## 1. Build Sub-Method Decomposition (SRNNModel2)

### What changed

`SRNNModel2.build()` was a single monolithic method. It is now decomposed into three protected sub-methods:

```matlab
function build(obj)
    obj.build_network();     % W generation
    obj.build_stimulus();    % Stimulus + interpolant + S0
    obj.finalize_build();    % Validation + caching
end
```

### Protected sub-methods

| Method | Access | Purpose | Override? |
|--------|--------|---------|-----------|
| `build_network(obj)` | protected | Seeds RNG, fills defaults (tilde params, tau_a), creates W via RMTMatrix | Override for custom connectivity |
| `build_stimulus(obj)` | protected | Generates stimulus, builds griddedInterpolant, initializes S0 | Override for custom stimulus (ESN does this) |
| `finalize_build(obj)` | protected | Calls validate(), caches params struct, sets is_built=true | Typically not overridden |

### How to override in a new subclass

A new subclass only needs to override the sub-method(s) that differ from the parent. Example:

```matlab
classdef MySubclass < SRNNModel2
    methods (Access = protected)
        function build_stimulus(obj)
            % Custom stimulus generation
            % Must set: obj.t_ex, obj.u_ex, obj.u_interpolant, obj.S0
            ...
        end
    end
end
```

The subclass does NOT need to override `build()` itself. The parent's `build()` dispatches to the sub-methods, and MATLAB's method resolution picks the subclass override.

### Removed method

`compute_derived_params()` was removed. Its job (computing n_E, n_I, E_indices, I_indices, N_sys_eqs) is now done by Dependent property getters that compute on access.

---

## 2. Property Changes (SRNNModel2)

### Made Dependent (removed from stored properties)

| Property | Computation | Rationale |
|----------|-------------|-----------|
| `n_E` | `round(f * n)` | Eliminates stale-value risk; trivially computable |
| `n_I` | `n - n_E` | Same |
| `E_indices` | `1:n_E` | Same |
| `I_indices` | `(n_E+1):n` | Same |
| `N_sys_eqs` | `n + n_E*n_a_E + n_I*n_a_I + n_E*n_b_E + n_I*n_b_I` | Same |

**Performance note**: These are accessed in the ODE loop via `cached_params` (a plain struct), not via property getters. Making them Dependent only affects pre-build access and the `get_params()` call. Zero ODE loop overhead.

### Removed

| Property | Reason |
|----------|--------|
| `row_zero_W` | Never referenced in SRNNModel2. Legacy from old `create_W_matrix.m`. |

### Annotated

| Property | Change |
|----------|--------|
| `reps` | Comment changed to "reserved for future use; typically unused" |

### Fixed

| Property | Issue | Fix |
|----------|-------|-----|
| `default_val` getter | Contained `fprintf(...)` — property getters should be side-effect-free | Removed the fprintf |

---

## 3. Property Changes (SRNN_ESN_reservoir)

### Moved to `SetAccess = protected`

| Property | Rationale |
|----------|-----------|
| `W_in` | Build-output. Should not be modified after build. Now protected like `W`. |

### Made Dependent

| Property | Computation | Rationale |
|----------|-------------|-----------|
| `dt_sample` | `1 / fs` | Trivially computable; was stale if fs changed after construction |

### Added

| Property | Access | Purpose |
|----------|--------|---------|
| `u_scalar` | `SetAccess = protected` | Scalar input sequence generated during build. Used by readout training and for pre-run verification. |

---

## 4. Stimulus Generation Moved to Build Time (SRNN_ESN_reservoir)

### Before

Stimulus was generated inside `run_memory_capacity()`:
1. `build()` called parent's build (which generated a wasted SRNNModel2 default stimulus)
2. `build()` then generated W_in
3. `run_memory_capacity()` generated u_scalar, then called `run_reservoir_esn(u_scalar)`
4. `run_reservoir_esn()` called `generate_esn_stimulus(u_scalar)` which overwrote t_ex, u_ex, u_interpolant, S0

### After

Stimulus is generated during `build()` via the `build_stimulus()` override:
1. `build()` calls `build_network()` (inherited — creates W)
2. `build()` calls `build_stimulus()` (overridden by ESN):
   - Generates W_in
   - Seeds RNG with `rng_seeds(2)`, generates u_scalar
   - Computes u_ex = W_in * u_scalar'
   - Creates piecewise-constant griddedInterpolant ('previous', 'nearest')
   - Initializes S0
3. `build()` calls `finalize_build()` (inherited — validates, caches)
4. `run_memory_capacity()` uses pre-built u_scalar directly; no stimulus generation
5. `run_reservoir_esn()` takes no arguments; uses obj.t_ex, obj.u_interpolant, obj.S0

### Removed methods

| Method | Reason |
|--------|--------|
| `generate_esn_stimulus(obj, u_scalar)` | Logic moved into `build_stimulus()` |

### Changed method signatures

| Method | Before | After |
|--------|--------|-------|
| `run_reservoir_esn` | `run_reservoir_esn(obj, u_scalar)` | `run_reservoir_esn(obj)` |

---

## 5. RNG Seeding Architecture

Three independent RNG streams, each explicitly seeded before use:

| Random object | Seed | Set in | Method |
|---------------|------|--------|--------|
| W (connectivity) | `rng_seeds(1)` | `build_network()` | `RMTMatrix` constructor |
| W_in (input weights) | `rng_seed_input` | `build_stimulus()` (ESN) | `generate_input_weights()` |
| u_scalar (input signal) | `rng_seeds(2)` | `build_stimulus()` (ESN) | Input generation |

Each stream resets the global RNG before use. No state leaks between them. Two ESN objects with identical user-config properties will produce identical W, W_in, and u_scalar regardless of build order.

---

## 6. Interpolation Differences

| Class | Method | Extrapolation | Rationale |
|-------|--------|---------------|-----------|
| SRNNModel2 | `'linear'` | `'none'` (NaN) | Smooth interpolation for step inputs |
| SRNN_ESN_reservoir | `'previous'` | `'nearest'` | Zero-order hold for discrete ESN input |

---

## 7. Verification Helper: `verify_shared_build.m`

External function for comparison experiments. Placed in `src/`.

```matlab
verify_shared_build(esn_array, expected_to_differ, also_check_protected)
```

- Uses `metaclass` introspection to auto-discover all properties
- Skips Dependent properties (auto-checked via their inputs)
- Skips run-output properties (S0, cached_params, mc_results, etc.)
- Checks all public non-Dependent properties match across objects
- Additionally checks specified protected properties (W, W_in, u_scalar, etc.)
- Asserts that `expected_to_differ` properties actually DO differ (catches misconfigured experiments)

---

## 8. Checklist for Applying to Another SRNNModel2 Subclass

When creating a new class that inherits from SRNNModel2:

1. **Inherit from SRNNModel2**, not SRNNModel (the old class)
2. **Do NOT override `build()`** — override `build_stimulus()` and/or `build_network()` instead
3. **In `build_stimulus()` override**, you must set:
   - `obj.t_ex` — time vector
   - `obj.u_ex` — stimulus matrix (n x T)
   - `obj.u_interpolant` — griddedInterpolant for ODE solver
   - `obj.S0` — initial state vector (use `initialize_state(obj.get_params())`)
4. **Use independent RNG seeds** — call `rng(seed)` before each random generation
5. **Store build-output properties** with `SetAccess = protected`
6. **Make trivially-computable properties Dependent** to avoid stale values
7. **Use `verify_shared_build()`** in scripts comparing multiple conditions
