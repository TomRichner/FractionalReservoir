# Design: cell-typed SRNN with STF (shared base + new class)

Design spec for the `cannonical-circuit` worktree: generalize the E/I `SRNNModel2` into a
**cell-type-resolved** rate reservoir (Pyr, Pvalb, Sst, Vip) and add **short-term facilitation
(STF)** alongside the existing SFA + STD, parameterized from the Allen `aisynphys` dataset
(see `docs/Litt/Campagnola_2022_data_summary.md`).

## Decisions (locked)

1. **Architecture: new standalone class + shared abstract base.** Extract the type-agnostic
   machinery of `SRNNModel2` into an abstract base class; both the existing E/I model and the
   new cell-typed model inherit it.
2. **State/RHS representation: per-neuron type-index gather.** Each neuron carries a type
   label; per-type parameters are gathered into length-`n` vectors and the RHS vectorizes over
   all `n` neurons uniformly. E/I is the K=2 special case. This replaces the current
   contiguous-block (`a_E; a_I; b_E; b_I; x`) layout in the *new* class.

## Current state (what we're generalizing)

`SRNNModel2` (2792 lines, `handle`) is E/I-binary in three layers: parameters (`*_E`/`*_I`
pairs), state layout (`S = [a_E; a_I; b_E; b_I; x]`), and connectivity (single-E/single-I RMT
via `RMTMatrix`). STF does not exist yet (only STD `b`). Hierarchy today:
`handle → SRNNModel2 → SRNN_ESN_reservoir`. The ESN subclass relies entirely on the inherited
numeric core.

**The dispatch constraint:** `run()` (line ~295) and its Jacobian (line ~287) call the kernels
**statically by hardcoded class name** — `SRNNModel2.dynamics_fast(...)` — and MATLAB static
methods do not dispatch polymorphically. The instance wrapper `dynamics(obj,t,S)` (line ~720)
also hardcodes `SRNNModel2.dynamics_fast`. So the base's `run()` must route through an
**overridable instance seam**, not a static class-name call.

## Target class hierarchy

```
handle
└── SRNNModelBase (abstract)         % type-agnostic machinery + abstract kernel hooks
    ├── SRNNModel2                   % E/I kernels (unchanged behavior)  [existing pipeline + ESN]
    │   └── SRNN_ESN_reservoir       % unchanged
    └── SRNNModelCellTypes           % K-type kernels + STF + block connectivity  [NEW]
```

### `SRNNModelBase` — what moves up (type-agnostic)

Owns everything that operates on a generic state vector `S(t)` given a kernel:
- `build()` template → `build_network` / `build_stimulus` / `finalize_build` hooks.
- `run()` (ODE orchestration), `compute_lyapunov()` / `filter_lyapunov()` and the
  Benettin/QR/Kaplan-Yorke static algorithms (they need only `S` + a Jacobian handle).
- Decimation, plotting infrastructure, storage/results props, stimulus generation
  (`generate_external_input`), RNG handling, `ParamSpaceAnalysis2` integration surface.
- Shared props: `n`, `tau_d`, `fs`, `T_range`, activation fn, sim/Lyapunov/storage settings, `W`.

**Abstract kernel hooks** (each subclass implements): `build_network`, `initialize_state`,
`get_params`, `validate`, `unpack_and_compute_states`, and the two numeric kernels
`dynamics(obj,t,S,params)` and `compute_jacobian(obj,S,params)`.

### Dispatch seam (the enabling refactor)

In `SRNNModelBase.run()` / `compute_lyapunov()`, replace hardcoded static calls with an
overridable instance seam:
```matlab
rhs = @(t, S) obj.dynamics(t, S, params);            % was SRNNModel2.dynamics_fast(...)
jac = @(t, S) obj.compute_jacobian(S, params);        % was SRNNModel2.compute_Jacobian_fast(...)
```
`SRNNModel2.dynamics` calls its existing static `dynamics_fast` internally, so **E/I numeric
results stay bit-identical.** This is pure indirection.

### `SRNNModel2` — E/I subclass (behavior-preserving)

Change `classdef SRNNModel2 < handle` → `< SRNNModelBase`; move the type-agnostic methods up to
the base; keep the E/I kernels (`dynamics_fast`, `compute_Jacobian_fast`, block
pack/unpack, `get_params`, `validate`) as its implementation of the hooks. **No change to E/I
outputs** — guarded by a regression test (below).

### `SRNNModelCellTypes` — the new class

Implements the hooks with the per-neuron representation and STF.

## Per-neuron type-index representation

- **Type assignment:** `type_of(i) ∈ {1..K}` for each neuron `i`, plus per-type counts. E→I
  sign still applies (a type is excitatory or inhibitory).
- **Per-type param tables:** for each mechanism, a per-type row → gathered to length-`n`
  vectors once at build time (e.g. `tau_a(i,:)`, `c(i)`, `tau_b_rec(i,:)`, `tau_b_rel(i)`,
  `tau_f(i)`, `f_amount(i)`). The RHS then vectorizes over all `n` neurons — no E/I branching.
- **State families (per neuron):**
  `a` (SFA, `n × n_a`), `b` (STD, `n × n_b`), **`u` (STF, `n × n_u`) — NEW**, and `x` (`n × 1`).
  Layout: `S = [a(:); b(:); u(:); x(:)]` with per-neuron rows. K types + STF are "just more
  per-neuron families," which is why this is simpler than K contiguous blocks.

## New mechanism: short-term facilitation (STF)

Add a facilitation state `u` mirroring how `b` (STD) works, driving the synaptic efficacy
multiplicatively. Candidate rate-form Tsodyks–Markram (exact form is a design item to finalize):

```
du_i/dt = (U0_i - u_i)/tau_f_i  +  f_amount_i * (1 - u_i) * r_i     % facilitation grows with rate
r_i     = phi(x_i - c_i * sum_k a_{i,k})                            % (SFA as today)
synaptic efficacy_i = g(u_i) * prod_m b_{i,m}                        % STF * STD, multiplicative
dx_i/dt = (-x_i + sum_j W_ij * efficacy_j * r_j + u^ext_i)/tau_d
```

Parameterization from `aisynphys`: `tau_f_i` ← `synapse_model.ml_facilitation_tau`,
`f_amount_i` ← `ml_facilitation_amount`; STD `tau_b_rec/rel` ← `ml_depression_*`. Alignment
rule to preserve when indexing STP params: **excitatory** synapse dynamics key on the
**post**-synaptic type, **inhibitory** on the **pre**-synaptic type (Campagnola §5). Open
question: whether STF/STD live on the presynaptic neuron (per-neuron `u_i`/`b_i`, as written)
or per (pre,post) synapse-type — start per-presynaptic-neuron; revisit if the pre/post
alignment demands per-edge state.

## Connectivity: cell-type block construction

Replace single-E/single-I RMT with a **block-structured W**: for each (pre-type, post-type),
sample connections at the Campagnola connection probability and scale nonzero weights by the
per-element strength. Extraction recipe + table→column map in
`docs/Litt/Campagnola_2022_data_summary.md §4c`. Start with 4 subclasses pooled across layers
(4×4 matrices). Keep `level_of_chaos` as an overall scale so edge-of-chaos sweeps still work.

## Reproducibility guard (must pass before/after the base extraction)

Before refactoring, save a reference run (LLE + Lyapunov spectrum + a state checkpoint) from a
fixed-seed `SRNNModel2` config. After extracting the base + rerouting the seam, assert the E/I
outputs are **bit-identical** (or within ODE tolerance). Add as `scripts/tests/test_base_refactor.m`.

## Phased implementation plan

1. **Base extraction (behavior-preserving):** create `SRNNModelBase`, move type-agnostic
   methods up, add the dispatch seam, make `SRNNModel2 < SRNNModelBase`. Verify ESN + a sample
   analysis still run and the regression test passes. *No new physics yet.*
2. **Data extraction:** pull the 4×4 connectivity/strength/STP/release-model matrices +
   per-subclass `adaptation_index` from `aisynphys` → CSV/`.mat` (see summary doc §4c/§5).
3. **New class scaffold:** `SRNNModelCellTypes` with per-neuron type-index state, K-type params,
   block `build_network`. Reuse SFA + STD kernels generalized to per-neuron gather; validate
   against `SRNNModel2` for a degenerate K=2 E/I config.
4. **Add STF:** the `u` state family + synaptic-efficacy coupling + Jacobian terms; finalize the
   functional form; unit-check facilitation vs. a known paired-pulse response.
5. **Parameterize from data:** wire the extracted matrices/indices; SFA via Option A
   (index → per-type strength `c`, fixed `τ_a`; design keeps `τ_a` per-type for future Option B).

## Open design items

- Exact STF functional form and whether efficacy uses `u` or `u·b` normalized.
- Per-presynaptic-neuron vs. per-(pre,post)-edge STP state (start per-neuron).
- Layer resolution later (16–17 classes) vs. 4 pooled subclasses now.
- Jacobian for STF (needed for the QR Lyapunov method — both kernels must stay analytically
  consistent).
