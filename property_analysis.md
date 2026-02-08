# Property Analysis: SRNNModel2 and SRNN_ESN_reservoir

Reference document for the planned refactor of the class hierarchy.
Analyzes every property, its current access level, lifecycle, and recommended changes.

## Table of Contents

- [Terminology](#terminology)
- [SRNNModel2 Properties](#srnnmodel2-properties)
  - [Network Architecture (public)](#network-architecture-public)
  - [Spike-Frequency Adaptation (public)](#spike-frequency-adaptation-public)
  - [Short-Term Depression (public)](#short-term-depression-public)
  - [Dynamics (public)](#dynamics-public)
  - [Simulation Settings (public)](#simulation-settings-public)
  - [Input Configuration (public)](#input-configuration-public)
  - [Lyapunov Settings (public)](#lyapunov-settings-public)
  - [Storage Options (public)](#storage-options-public)
  - [RMT Dependent Properties](#rmt-dependent-properties)
  - [Computed / Build-Output Properties (SetAccess = protected)](#computed--build-output-properties)
  - [Results Properties (SetAccess = protected)](#results-properties)
- [SRNN_ESN_reservoir Properties](#srnn_esn_reservoir-properties)
  - [ESN Input Properties (public)](#esn-input-properties-public)
  - [Memory Capacity Protocol Properties (public)](#memory-capacity-protocol-properties-public)
  - [Memory Capacity Results (SetAccess = private)](#memory-capacity-results)
- [Refactor Recommendations Summary](#refactor-recommendations-summary)
- [Stimulus Generation Notes](#stimulus-generation-notes)
- [Build Sub-Method Decomposition Notes](#build-sub-method-decomposition-notes)
- [Verification Function Design Notes](#verification-function-design-notes)

---

## Terminology

| Term | Meaning |
|------|---------|
| **User-config** | Set by the user before `build()`. Defines the experiment. |
| **Build-output** | Computed during `build()`. Should not be modified after build. |
| **Run-output** | Computed during `run()` / `run_memory_capacity()`. |
| **Defaulted-in-build** | User-config property that gets a default value during `build()` if the user left it empty. |
| **Derived** | Computable from other properties (candidate for `Dependent`). |
| **Legacy** | Property that exists but is not used in the current code. |

---

## SRNNModel2 Properties

### Network Architecture (public)

| Property | Default | Lifecycle | Used in build? | Recommendation | Notes |
|----------|---------|-----------|----------------|----------------|-------|
| `n` | 300 | User-config | Yes (W, derived params) | **SetAccess = immutable after build** | Core dimension. Changing after build invalidates everything. |
| `f` | 0.5 | User-config | Yes (n_E/n_I split) | **SetAccess = immutable after build** | Changing after build invalidates W, indices, params. |
| `indegree` | 100 | User-config | Yes (alpha → sparsity) | **SetAccess = immutable after build** | |
| `mu_E_tilde` | `[]` | Defaulted-in-build | Yes (RMTMatrix) | Keep public. Note: written during `build()` if empty — set to `3*F`. After build, treat as immutable. | Dual lifecycle: user can set before build, or build fills default. |
| `mu_I_tilde` | `[]` | Defaulted-in-build | Yes (RMTMatrix) | Same as `mu_E_tilde` | Default: `-4*F` |
| `sigma_E_tilde` | `[]` | Defaulted-in-build | Yes (RMTMatrix) | Same as `mu_E_tilde` | Default: `F` |
| `sigma_I_tilde` | `[]` | Defaulted-in-build | Yes (RMTMatrix) | Same as `mu_E_tilde` | Default: `F` |
| `E_W` | 0 | User-config | Yes (RMTMatrix offset) | Keep public | Mean offset applied to both E and I tilde means. |
| `zrs_mode` | `'none'` | User-config | Yes (RMTMatrix) | Keep public | Row zero-sum mode. |
| `level_of_chaos` | 1.0 | User-config | Yes (W scaling) | Keep public | Scaling factor for W. |
| `rescale_by_abscissa` | `false` | User-config | Yes (W scaling) | Keep public | Whether to normalize by spectral abscissa. |
| `row_zero_W` | `true` | **Legacy** | **No** | **Remove.** Not referenced anywhere in SRNNModel2 code. Only used by the old `create_W_matrix.m` (not called by SRNNModel2). | Dead code. |

### Spike-Frequency Adaptation (public)

| Property | Default | Lifecycle | Used in build? | Recommendation | Notes |
|----------|---------|-----------|----------------|----------------|-------|
| `n_a_E` | 0 | User-config | Yes (N_sys_eqs, tau_a defaults) | Keep public. **This is the primary experimental variable for the MC comparison.** | |
| `n_a_I` | 0 | User-config | Yes (same) | Keep public | |
| `tau_a_E` | `[]` | Defaulted-in-build | Yes (params) | Keep public. Written in `build()` if empty and `n_a_E > 0` (logspace defaults). | Default: `logspace(log10(0.25), log10(10), n_a_E)` |
| `tau_a_I` | `[]` | Defaulted-in-build | Yes (params) | Same as `tau_a_E` | |
| `c_E` | 0.15/3 | User-config | No (only in dynamics at runtime) | Keep public | Adaptation coupling strength. |
| `c_I` | 0.15/3 | User-config | No (only in dynamics at runtime) | Keep public | |

### Short-Term Depression (public)

| Property | Default | Lifecycle | Used in build? | Recommendation | Notes |
|----------|---------|-----------|----------------|----------------|-------|
| `n_b_E` | 0 | User-config | Yes (N_sys_eqs) | Keep public. **Primary experimental variable for MC comparison.** | |
| `n_b_I` | 0 | User-config | Yes (N_sys_eqs) | Keep public | |
| `tau_b_E_rec` | 1 | User-config | No (only in dynamics) | Keep public | STD recovery time constant. |
| `tau_b_E_rel` | 0.25 | User-config | No (only in dynamics) | Keep public | STD release time constant. |
| `tau_b_I_rec` | 1 | User-config | No (only in dynamics) | Keep public | |
| `tau_b_I_rel` | 0.25 | User-config | No (only in dynamics) | Keep public | |

### Dynamics (public)

| Property | Default | Lifecycle | Used in build? | Recommendation | Notes |
|----------|---------|-----------|----------------|----------------|-------|
| `tau_d` | 0.1 | User-config | No (runtime only) | Keep public | Dendritic time constant. |
| `activation_function` | (set in `set_defaults`) | Derived | No (runtime only) | **Make Dependent on `S_a` and `S_c`.** Currently set once in constructor via `set_defaults()`, but if user changes `S_a` or `S_c` after construction, the function handle is stale. Making it Dependent ensures consistency. | Currently: `@(x) piecewiseSigmoid(x, obj.S_a, obj.S_c)` — captures `obj`, so the closure does update... but only because it captures `obj` by reference. This is fragile. |
| `activation_function_derivative` | (set in `set_defaults`) | Derived | No (runtime only) | **Make Dependent on `S_a` and `S_c`.** Same issue as `activation_function`. | |
| `S_a` | 0.9 | User-config | No | Keep public | Activation function shape parameter. |
| `S_c` | 0.35 | User-config | No | Keep public | Activation function center parameter. |

**Note on `activation_function` / `activation_function_derivative`**: These function handles are created in `set_defaults()` as closures that capture `obj`. Because `obj` is a handle, the closures *do* see updated `S_a`/`S_c` values. However, these closures are cached into `cached_params.activation_function` during `build()`, and `cached_params` is a plain struct — so if `S_a` or `S_c` is changed after build, the cached version is stale. Making them Dependent would not help with `cached_params` caching. The real fix is: don't change `S_a`/`S_c` after build.

### Simulation Settings (public)

| Property | Default | Lifecycle | Used in build? | Recommendation | Notes |
|----------|---------|-----------|----------------|----------------|-------|
| `fs` | 400 | User-config | Yes (dt, t_ex, stimulus generation) | Keep public. Immutable after build. | Sampling frequency. |
| `T_range` | `[0, 50]` | User-config | Yes (stimulus generation, Lyapunov) | Keep public | Simulation time interval. Not used by ESN (ESN uses T_wash+T_train+T_test). |
| `T_plot` | `[]` | User-config | No (plotting only) | Keep public | Plotting time interval, defaults to `T_range`. |
| `ode_solver` | `@ode45` | User-config | No (used in run) | Keep public | |
| `ode_opts` | `[]` | Mixed | **Problematic.** Lazy-initialized during `run()` or `run_reservoir_esn()` if empty. Has different defaults in SRNNModel2 (`RelTol=1e-9`) vs SRNN_ESN_reservoir (`RelTol=1e-5`). | Keep public. The lazy init is fine — but the divergent defaults between parent/child should be documented. Consider setting ESN defaults in ESN's `build()`. | |

### Input Configuration (public)

| Property | Default | Lifecycle | Used in build? | Recommendation | Notes |
|----------|---------|-----------|----------------|----------------|-------|
| `input_config` | (struct, set in `set_defaults`) | User-config | Yes (stimulus generation) | Keep public. **Not used by ESN** — the ESN has its own input system (`input_type`, `u_f_cutoff`, etc.). | The ESN calls `build@SRNNModel2` which calls `generate_stimulus()` using `input_config`, but that stimulus is then overwritten. This is the wasted work identified in the investigation. |
| `u_ex_scale` | 1.0 | User-config | Yes (stimulus scaling) | Keep public. **Not used by ESN.** | |
| `rng_seeds` | `[1, 2]` | User-config | Yes (W generation, stimulus) | Keep public. **Critical for reproducibility.** `rng_seeds(1)` → W. `rng_seeds(2)` → stimulus. | |
| `reps` | 1 | **Legacy** | **No** | **Remove or document.** Comment says "used by ParamSpaceAnalysis" but `obj.reps` is never referenced in any file in the project. If ParamSpaceAnalysis accesses it via the params struct, it would use `get_params()`, but `reps` is not in `get_params()`. | Dead code unless used externally. |

### Lyapunov Settings (public)

| Property | Default | Lifecycle | Used in build? | Recommendation | Notes |
|----------|---------|-----------|----------------|----------------|-------|
| `lya_method` | `'benettin'` | User-config | No (run-time) | Keep public | |
| `lya_T_interval` | `[]` | Mixed | No (set lazily during `compute_lyapunov()` if empty; also written by ESN during `run_memory_capacity()`) | Keep public. Lazy default is OK. | |
| `filter_local_lya` | `false` | User-config | No | Keep public | |
| `lya_filter_order` | 2 | User-config | No | Keep public | |
| `lya_filter_cutoff` | 0.25 | User-config | No | Keep public | |

### Storage Options (public)

| Property | Default | Lifecycle | Used in build? | Recommendation | Notes |
|----------|---------|-----------|----------------|----------------|-------|
| `store_full_state` | `false` | User-config | No (run-time) | Keep public | |
| `store_decimated_state` | `true` | User-config | No (run-time) | Keep public | |
| `plot_deci` | `[]` | Derived | No (computed in constructor from `fs/plot_freq` if empty) | **Make Dependent on `fs` and `plot_freq`**, or keep current lazy-init pattern. | Currently set once in constructor. If `fs` changes after construction, `plot_deci` is stale. |
| `plot_freq` | 10 | User-config | No | Keep public | |

### RMT Dependent Properties

All correctly marked as `Dependent`. No storage — computed on access.

| Property | Depends on | Recommendation | Notes |
|----------|------------|----------------|-------|
| `alpha` | `indegree`, `n` | Keep Dependent | `indegree/n` |
| `default_val` | `n`, `alpha` | Keep Dependent. **Remove `fprintf` from getter** — property getters should be side-effect-free. The `fprintf` fires on every access, including during build. | Currently prints "default_val (F) computed..." every time it's accessed. |
| `mu_se` | `alpha`, `mu_E_tilde`, `E_W` | Keep Dependent | |
| `mu_si` | `alpha`, `mu_I_tilde`, `E_W` | Keep Dependent | |
| `sigma_se` | `alpha`, `mu_E_tilde`, `E_W`, `sigma_E_tilde` | Keep Dependent | |
| `sigma_si` | `alpha`, `mu_I_tilde`, `E_W`, `sigma_I_tilde` | Keep Dependent | |
| `R` | `n`, `f`, `sigma_se`, `sigma_si`, `level_of_chaos` | Keep Dependent | Theoretical spectral radius. |

### Computed / Build-Output Properties

Currently `SetAccess = protected` (readable by anyone, writable by class and subclasses).

| Property | Set during | Recommendation | Notes |
|----------|------------|----------------|-------|
| `W` | `build()` | **Should be immutable after build.** Consider `SetAccess = protected` + don't expose a setter. Current access is fine structurally — no external code modifies W. | 300×300 matrix. Core random structure. |
| `n_E` | `build()` | **Make Dependent on `n` and `f`.** It's just `round(f * n)`. No reason to store it. | Currently computed in `compute_derived_params()`. |
| `n_I` | `build()` | **Make Dependent.** `n - n_E`. | |
| `E_indices` | `build()` | **Make Dependent.** `1:n_E`. | |
| `I_indices` | `build()` | **Make Dependent.** `(n_E+1):n`. | |
| `N_sys_eqs` | `build()` | **Make Dependent.** `n + n_E*n_a_E + n_I*n_a_I + n_E*n_b_E + n_I*n_b_I`. Fully determined by user-config properties. | Currently computed in `compute_derived_params()`. |
| `is_built` | `build()` | Keep as build-output flag. `SetAccess = protected`. | |
| `t_ex` | `build()` (SRNNModel2) or `generate_esn_stimulus()` (ESN) | Build-output. Should be immutable after build. | Overwritten by ESN during `run_memory_capacity()` — this is the lifecycle issue to fix. |
| `u_ex` | `build()` (SRNNModel2) or `generate_esn_stimulus()` (ESN) | Build-output. Same lifecycle issue as `t_ex`. | |
| `u_interpolant` | `build()` (SRNNModel2) or `generate_esn_stimulus()` (ESN) | Build-output. Same lifecycle issue. | SRNNModel2 uses `'linear'` interpolation; ESN uses `'previous'` (zero-order hold). |
| `S0` | `build()` + overwritten in ESN `generate_esn_stimulus()` | Build-output. Should be set once during build. | |
| `cached_params` | `build()` | Build-output. **Large struct duplicating many properties.** Exists for performance (avoids OOP property access in tight ODE loop). | Updated by ESN in `generate_esn_stimulus()` (`cached_params.u_interpolant`). |

### Results Properties

Currently `SetAccess = protected`.

| Property | Set during | Recommendation | Notes |
|----------|------------|----------------|-------|
| `t_out` | `run()` / `run_reservoir_esn()` | Run-output. Keep protected. | |
| `S_out` | `run()` / `run_reservoir_esn()` | Run-output. May be cleared after run to save memory. | |
| `plot_data` | `run()` (decimate_and_unpack) | Run-output. | |
| `lya_results` | `run()` / `run_memory_capacity()` | Run-output. | |
| `has_run` | `run()` | Run-output flag. | |

---

## SRNN_ESN_reservoir Properties

### ESN Input Properties (public)

| Property | Default | Lifecycle | Used in build? | Recommendation | Notes |
|----------|---------|-----------|----------------|----------------|-------|
| `f_in` | 0.1 | User-config | Yes (`generate_input_weights`) | Keep public. Immutable after build. | Fraction of neurons receiving input. |
| `sigma_in` | 0.5 | User-config | Yes (`generate_input_weights`) | Keep public. Immutable after build. | Input weight scaling. |
| `W_in` | `[]` | **Build-output** | Set during `build()` | **Change to `SetAccess = protected` (or private).** Currently public — user could accidentally overwrite after build. It's set by `generate_input_weights()` and should be treated as immutable after build, like `W`. | 300×1 vector. Core random structure. |
| `rng_seed_input` | 3 | User-config | Yes (`generate_input_weights`) | Keep public | Separate RNG seed for W_in generation. |
| `input_type` | `'white'` | User-config | Will be used in build (after refactor) | Keep public | `'white'`, `'bandlimited'`, or `'one_over_f'`. |
| `u_f_cutoff` | `[]` | User-config | Will be used in build (after refactor) | Keep public | Cutoff frequency for bandlimited input. If empty, defaults to `1/(2*pi*tau_d)`. |
| `u_alpha` | 1 | User-config | Will be used in build (after refactor) | Keep public | Spectral exponent for 1/f noise. |
| `u_scale` | 1 | User-config | Will be used in build (after refactor) | Keep public | Stimulus amplitude scaling. |
| `u_offset` | 0 | User-config | Will be used in build (after refactor) | Keep public | Stimulus DC offset. |

### Memory Capacity Protocol Properties (public)

| Property | Default | Lifecycle | Used in build? | Recommendation | Notes |
|----------|---------|-----------|----------------|----------------|-------|
| `T_wash` | 1000 | User-config | Will be used in build (after refactor, for T_total) | Keep public | Washout samples. |
| `T_train` | 5000 | User-config | Will be used in build (after refactor, for T_total) | Keep public | Training samples. |
| `T_test` | 5000 | User-config | Will be used in build (after refactor, for T_total) | Keep public | Test samples. |
| `d_max` | 70 | User-config | No (run-time readout training) | Keep public | Maximum delay for MC measurement. |
| `eta` | 1e-7 | User-config | No (run-time readout training) | Keep public | Ridge regression regularization. |
| `dt_sample` | (computed) | **Derived** | No | **Make Dependent on `fs`.** Currently set once in constructor as `1/fs`. If `fs` changes after construction, `dt_sample` is stale. Trivial computation — no reason to store. | `1 / obj.fs` |

### Memory Capacity Results

| Property | Access | Lifecycle | Recommendation | Notes |
|----------|--------|-----------|----------------|-------|
| `mc_results` | `SetAccess = private` | Run-output | Keep private. | Struct containing all MC results. |

### Missing Property (to be added in refactor)

| Property | Recommended access | Lifecycle | Notes |
|----------|-------------------|-----------|-------|
| `u_scalar` | `SetAccess = protected` | Build-output (after refactor) | The scalar input sequence. Currently generated as a local variable inside `run_memory_capacity()`. Should be stored for: (1) pre-run verification, (2) access during readout training. |

---

## Refactor Recommendations Summary

### Properties to make Dependent

These are currently stored but are trivially computable from other properties:

| Property | Class | Computation | Benefit |
|----------|-------|-------------|---------|
| `n_E` | SRNNModel2 | `round(f * n)` | Eliminates stale-value risk |
| `n_I` | SRNNModel2 | `n - n_E` | |
| `E_indices` | SRNNModel2 | `1:n_E` | |
| `I_indices` | SRNNModel2 | `(n_E+1):n` | |
| `N_sys_eqs` | SRNNModel2 | `n + n_E*n_a_E + n_I*n_a_I + n_E*n_b_E + n_I*n_b_I` | |
| `dt_sample` | ESN | `1 / fs` | |

**Performance note**: `n_E`, `n_I`, `E_indices`, `I_indices`, and `N_sys_eqs` are accessed frequently in the ODE loop via `cached_params`, not via the property getters. Making them Dependent only affects pre-build access and the `get_params()` call during build. The ODE loop uses `params.n_E`, etc., which are plain struct fields — no overhead.

### Properties to remove

| Property | Class | Reason |
|----------|-------|--------|
| `row_zero_W` | SRNNModel2 | Never referenced in SRNNModel2. Only exists in old `create_W_matrix.m` (not called). |
| `reps` | SRNNModel2 | Never referenced by `obj.reps` anywhere. Dead code. |

### Properties to restrict access

| Property | Class | Current | Recommended | Reason |
|----------|-------|---------|-------------|--------|
| `W_in` | ESN | public | `SetAccess = protected` | Build-output. Should not be modified after build. |

### Property getters to fix

| Property | Class | Issue |
|----------|-------|-------|
| `default_val` | SRNNModel2 | Getter contains `fprintf(...)`. Property getters should be side-effect-free. Remove the print statement. |

### New property to add

| Property | Class | Access | Notes |
|----------|-------|--------|-------|
| `u_scalar` | ESN | `SetAccess = protected` | Scalar input sequence, generated during `build()` (after refactor). |

---

## Stimulus Generation Notes

### Current flow (pre-refactor)

```
SRNNModel2.build()
├─ rng(rng_seeds(1))
├─ Create W                                ← uses rng_seeds(1)
├─ generate_stimulus()                     ← calls rng(rng_seeds(2)) internally
│  └─ generate_external_input(...)         → obj.u_ex, obj.t_ex
├─ griddedInterpolant('linear', 'none')    → obj.u_interpolant
├─ initialize_state()                      → obj.S0  [uses randn]
└─ cache params

SRNN_ESN_reservoir.build()
├─ build@SRNNModel2(obj)                   ← above flow (u_ex/S0 are WASTED for ESN)
└─ generate_input_weights()                ← rng(rng_seed_input) → obj.W_in

SRNN_ESN_reservoir.run_memory_capacity()
├─ rng(rng_seeds(2))
├─ Generate u_scalar                       ← stimulus created HERE (late)
├─ run_reservoir_esn(u_scalar)
│  └─ generate_esn_stimulus(u_scalar)
│     ├─ obj.u_ex = W_in * u_scalar'       ← OVERWRITES parent's u_ex
│     ├─ griddedInterpolant('previous')    ← OVERWRITES parent's interpolant
│     └─ obj.S0 = initialize_state(...)    ← OVERWRITES parent's S0
└─ [readout training, R² computation]
```

### Planned flow (post-refactor)

```
SRNNModel2.build()
├─ build_network()         [PROTECTED, overridable]
│  ├─ rng(rng_seeds(1))
│  ├─ compute defaults (tilde params, tau_a)
│  └─ Create W via RMTMatrix
├─ build_stimulus()        [PROTECTED, overridable]  ← ESN overrides this
│  ├─ generate_stimulus()
│  ├─ griddedInterpolant('linear', 'none')
│  └─ initialize_state() → S0
└─ finalize_build()        [PROTECTED, overridable]
   ├─ validate()
   └─ cache params, is_built = true

SRNN_ESN_reservoir.build()
├─ build@SRNNModel2(obj)
│  ├─ build_network()                      ← inherited, creates W
│  ├─ build_stimulus()                     ← OVERRIDDEN by ESN
│  │  ├─ generate_input_weights()          → W_in
│  │  ├─ rng(rng_seeds(2))
│  │  ├─ Generate u_scalar                 → obj.u_scalar  (NEW property)
│  │  ├─ obj.u_ex = W_in * u_scalar'
│  │  ├─ griddedInterpolant('previous', 'nearest')
│  │  └─ initialize_state() → S0
│  └─ finalize_build()                     ← inherited
└─ (done — fully built, ready for verification and run)

SRNN_ESN_reservoir.run_memory_capacity()
├─ [no stimulus generation — already built]
├─ run_reservoir_esn()                     ← uses obj.t_ex, obj.u_interpolant, obj.S0
│  └─ ode_solver(rhs, obj.t_ex, obj.S0, opts)
└─ [readout training uses obj.u_scalar, R² computation]
```

### Interpolation method differences

| Class | Method | Extrapolation | Rationale |
|-------|--------|---------------|-----------|
| SRNNModel2 | `'linear'` | `'none'` (NaN) | Smooth interpolation between grid points. NaN catches solver overshoot beyond `t_ex`. |
| SRNN_ESN_reservoir | `'previous'` | `'nearest'` | Zero-order hold (sample-and-hold). Correct for ESN discrete-time input semantics. `'nearest'` avoids NaN at boundaries. |

### Why griddedInterpolant exists

The ODE solver (e.g., ode45) takes adaptive internal steps that don't align with the stimulus grid.
At each internal sub-step, the RHS function calls `params.u_interpolant(t)` to get the input value.
The `griddedInterpolant` is created once (during build) and provides O(log n) lookup thereafter.
Without it, each RHS evaluation would need to manually search for the correct time bin — far slower.

### RNG seeding architecture (post-refactor)

Three independent RNG streams, each with its own explicit `rng(seed)` call:

| Random object | Seed | Set where | Consumed by |
|--------------|------|-----------|-------------|
| W (connectivity) | `rng_seeds(1)` | `build_network()` | `RMTMatrix` constructor (`randn`, `rand`) |
| W_in (input weights) | `rng_seed_input` | `build_stimulus()` (ESN override) | `generate_input_weights()` (`randperm`, `rand`) |
| u_scalar (input signal) | `rng_seeds(2)` | `build_stimulus()` (ESN override) | Input generation (`rand`, `randn` depending on type) |

Each stream is independently seeded — no RNG state leaks between them.
If two ESN objects share all user-config properties, they will produce identical W, W_in, and u_scalar regardless of build order.

---

## Build Sub-Method Decomposition Notes

The proposed `SRNNModel2.build()` decomposition:

```matlab
function build(obj)
    obj.build_network();     % W generation — override for custom connectivity
    obj.build_stimulus();    % Stimulus + interpolant + S0 — ESN overrides this
    obj.finalize_build();    % Validation + caching
end
```

These should be `Access = protected` so subclasses can override them but external code cannot call them directly.

### What goes in each sub-method

**`build_network()` (protected)**
- `rng(rng_seeds(1))`
- `compute_derived_params()` (or removed if those become Dependent)
- Set tilde defaults if empty
- Set tau_a defaults if empty
- Create W via RMTMatrix
- Scale W (level_of_chaos, optional abscissa rescaling)
- Print spectral radius info

**`build_stimulus()` (protected, overridden by ESN)**
- SRNNModel2 default: `generate_stimulus()` + `griddedInterpolant('linear')` + `initialize_state()`
- ESN override: `generate_input_weights()` + `generate_esn_stimulus_at_build()` (u_scalar → u_ex → interpolant → S0)

**`finalize_build()` (protected)**
- `validate()`
- `cached_params = get_params()`
- `is_built = true`

### Existing subclass impact

Only one subclass exists: `SRNN_ESN_reservoir`. No other classes inherit from `SRNNModel2`. The decomposition is safe.

---

## Verification Function Design Notes

### Goal

After building multiple ESN objects for a comparison experiment, verify that all properties that *should* be shared are in fact identical, while allowing a user-specified set of properties to differ.

### Approach: property-category-aware comparison

Rather than comparing all public properties and subtracting exceptions, categorize properties by their role:

1. **User-config properties that define the shared experiment** — these MUST match (n, fs, rng_seeds, tau_d, S_a, S_c, etc.)
2. **User-config properties that intentionally differ** — the user specifies these (n_a_E, n_b_E, etc.)
3. **Build-output properties that derive from shared config** — these MUST match (W, W_in, u_ex, u_scalar, u_interpolant)
4. **Build-output properties that derive from differing config** — these WILL differ (N_sys_eqs, S0, cached_params)
5. **Dependent properties** — skip (they're computed from other properties)
6. **Run-output properties** — skip (not yet populated)

### Implementation sketch

```matlab
function verify_shared_build(esn_array, allowed_to_differ)
    % VERIFY_SHARED_BUILD Verify ESN objects share structure after build
    %
    % Inputs:
    %   esn_array        - Cell array of built SRNN_ESN_reservoir objects
    %   allowed_to_differ - Cell array of property names expected to differ
    %                       e.g., {'n_a_E', 'n_b_E', 'tau_a_E', 'tau_b_E_rec', ...}

    ref = esn_array{1};
    mc = metaclass(ref);

    for i = 2:numel(esn_array)
        obj = esn_array{i};
        for p = 1:numel(mc.PropertyList)
            prop = mc.PropertyList(p);

            % Skip: Dependent, private, protected-get, Constant, run-output
            if prop.Dependent, continue; end
            if ~strcmp(prop.GetAccess, 'public'), continue; end

            name = prop.Name;

            % Skip properties expected to differ
            if ismember(name, allowed_to_differ), continue; end

            % Skip known run-output and derived-from-differing properties
            % (N_sys_eqs, S0, cached_params, mc_results, etc.)
            % These could be identified by a class-level attribute or a hardcoded skip list

            % Compare
            val_ref = ref.(name);
            val_obj = obj.(name);

            if ~isequaln(val_ref, val_obj)
                error('verify_shared_build:Mismatch', ...
                    'Property ''%s'' differs between condition 1 and %d.', name, i);
            end
        end
    end
    fprintf('Verified: all %d conditions share identical base configuration.\n', numel(esn_array));
end
```

### Key design decisions

1. **Use `metaclass` introspection** — dynamically discover all properties rather than hardcoding. This ensures new properties are automatically checked.

2. **Skip Dependent properties** — they're computed from other properties; if the inputs match, the outputs match.

3. **The `allowed_to_differ` list should be small** — typically just `{'n_a_E', 'n_b_E'}` and their associated params (`tau_a_E`, `c_E`, `tau_b_E_rec`, `tau_b_E_rel`).

4. **Properties derived from the differing config** (N_sys_eqs, S0, cached_params) must also be skipped. If `n_E`, `n_I`, `E_indices`, `I_indices`, `N_sys_eqs` are made Dependent (as recommended above), they are automatically skipped. `S0` and `cached_params` have `SetAccess = protected` but `GetAccess = public` — they would be checked. Since their values legitimately differ (different state dimensions), they need to be in a "derived-from-differing" skip list.

5. **Recommendation**: Rather than a hardcoded skip list for derived properties, compute the skip list dynamically. If a property is `SetAccess = protected` (i.e., build-output), check whether it *could* depend on the `allowed_to_differ` properties. In practice, the simplest approach is:
   - Compare all public-get, non-Dependent properties
   - For mismatches: check if the property is in `allowed_to_differ` — if so, skip
   - For remaining mismatches: report as error
   - For build-output mismatches (`S0`, `cached_params`): warn but don't error, since these are *expected* to differ when adaptation config differs

### Properties that will legitimately differ when n_a_E / n_b_E differ

Even with all other config identical, changing `n_a_E` or `n_b_E` causes these to differ:

| Property | Why it differs | Current access |
|----------|---------------|----------------|
| `n_a_E`, `n_b_E` | Intentionally different | public |
| `tau_a_E` | Defaulted differently (empty vs logspace) | public |
| `tau_a_I`, `tau_b_E_rec`, `tau_b_E_rel`, `c_E`, `c_I` | May or may not differ depending on experiment design | public |
| `N_sys_eqs` | Different state dimension | protected (Dependent after refactor) |
| `S0` | Different length (different N_sys_eqs) | protected |
| `cached_params` | Contains all of the above | protected |
| `ode_opts` | Lazy-initialized, may differ if Jacobian wrapper captures params | public (set lazily) |

If `N_sys_eqs` is made Dependent, it's automatically skipped. `S0` and `cached_params` have protected SetAccess but public GetAccess — the verification function should identify these as build-output and handle them specially.
