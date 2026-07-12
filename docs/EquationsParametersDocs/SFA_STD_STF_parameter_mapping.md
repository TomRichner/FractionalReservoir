# Mapping Campagnola 2022 data to SRNN parameters (SFA / STD / STF)

How the cell-type-resolved reservoir `SRNNModelCellTypes` gets its parameters from the
Allen Institute synaptic-physiology dataset (Campagnola et al. 2022). This documents the
**data → parameter** path only; for the dynamical equations see
[`system_equations_with_stf.md`](system_equations_with_stf.md). Regenerate every figure below with
`param_visualizations.m` (co-located in this folder).

## 1. Overview

The model has **four cell types** — Pyr, Pvalb, Sst, Vip — in that fixed order. Every synaptic
data matrix is **4×4 with rows = presynaptic type, columns = postsynaptic type**; intrinsic
(single-cell) properties are 4×1 by type. Excitation/inhibition follows **Dale's law by
presynaptic type**: the Pyr row of the PSP amplitude is positive, the Pvalb/Sst/Vip rows negative,
so a neuron's outgoing sign is set by its type.

**Provenance.** Matrices are extracted from the aisynphys `small` database (mouse cortex, pooled
across layers) by `scripts/aisynphys/extract_campagnola_matrices.py`, committed as CSVs under
`src/connectivity/campagnola/`, and read by `src/connectivity/load_campagnola_matrices.m`. Elements
with < 2 sampled pairs are `NaN`. See `src/connectivity/campagnola/PROVENANCE.md` for the DB version
and queries. In the current dataset **all 4×4 elements used are populated** (no `NaN` fallbacks fire).

**Where the parameters live in code.** `SRNNModelCellTypes.load_parameter_tables` copies the
Campagnola matrices into the model's K×K tables (with NaN fallbacks + range clamps);
`SRNNModelCellTypes.get_params` then derives the runtime parameters and gathers them per
`(presynaptic neuron j, postsynaptic type T)`, i.e. `param(j,T) = table(type(j), T)` — the full
`(pre, post)` element, so both of Campagnola's alignment rules (excitatory dynamics track the
postsynaptic type; inhibitory mostly the presynaptic type) fall out with no special-casing.

## 2. The mapping at a glance

| Campagnola field (source) | Model parameter | Expression | Rationale |
|---|---|---|---|
| `conn_prob_adj` (P→Q) | Bernoulli connection mask | edge present w.p. `conn_prob_adj(P,Q)` | distance-adjusted connection probability sets block sparsity |
| `psp_amplitude` (P→Q, signed, V) | synaptic weight `w` | `w ~ psp_amplitude(P,Q)·(1 + w_cv·ξ)`, then global `level_of_chaos` scaling | signed resting PSP sets strength **and** E/I sign (Dale by pre-type) |
| `sfa_adaptation_index` (per type) | SFA strength `c` | `c_θ = c_gain · adaptation_index_θ` (`c_gain = 0.7`) | relative adaptation strength per type; single fixed `τ_a` |
| `ml_depression_tau` (P→Q, s) | STD recovery `τ_rec` | `τ_rec = clamp(ml_depression_tau, 0.05, 5)` | release-model depression time constant used directly |
| `ml_depression_amount` (P→Q) | STD release `τ_rel` | `τ_rel = max(τ_rel_ref·(1 − clamp(amount,0,0.95)), 0.05)`, `τ_rel_ref = 0.5` | **provisional heuristic**: more depression → shorter `τ_rel` → faster/deeper depletion |
| `ml_release_prob` (P→Q) | baseline release prob `p0` | `p0 = clamp(ml_release_prob, 0.05, 0.95)` | resting release probability = STF floor + depression scale |
| `ml_facilitation_tau` (P→Q, s) | facilitation decay `τ_f` | `τ_f = clamp(ml_facilitation_tau, 0.05, 5)` | facilitation (residual-calcium) decay time |
| `ml_facilitation_amount` (P→Q) | facilitation rate `κ` | `κ = max(ml_facilitation_amount, 0)` | facilitation strength, used directly as the rate `κ` |

`ξ ~ N(0,1)` per edge; `w_cv` (default 1.0) sets per-edge weight heterogeneity. The Campagnola time
constants are in **seconds** and are used directly; because the model rate `r = φ(·) ∈ (0,1)` is
dimensionless, `τ_rel` is a modeling timescale (default reference `0.5 s`), not a measured quantity.

## 3. Campagnola source data used

![Campagnola source data used by the model](figures/param_map_campagnola_inputs.png)

Connection probability, signed PSP amplitude (Pyr row red/excitatory, interneuron rows
blue/inhibitory), release probability, the SFA index per type, and the four release-model STF/STD
metrics. Note `ml_depression_amount` is **signed** — negative entries (e.g. Vip→Pvalb) indicate net
facilitation in the release model and are clamped to 0 before the `τ_rel` heuristic.

## 4. Connectivity and strength → W

For a presynaptic neuron `j` of type `P` and a postsynaptic neuron `i` of type `Q`, an edge exists
with probability `conn_prob_adj(P,Q)` and, if present, has weight
`w_ij ~ psp_amplitude(P,Q)·(1 + w_cv·ξ)`. The whole `W` is then scaled to the edge of chaos
(`rescale_by_abscissa` normalizes the spectral abscissa to 1, then multiplies by `level_of_chaos`),
so absolute PSP units drop out and `level_of_chaos` sets the operating point. The **expected**
(sparsity-weighted) block weight `conn_prob_adj ∘ psp_amplitude` is the top-left panel of the derived
figure in §7; individual realized synapses carry `~psp_amplitude`.

## 5. Spike-frequency adaptation (SFA)

SFA is an **intrinsic, per-neuron** property (not per synapse). Each neuron gets a single adaptation
current `a` with a fixed timescale `τ_a` (default `1 s`, a modeling choice — **Option A**: the
Campagnola *adaptation index* is mapped to a relative *strength*, not to a per-type `τ_a`), and
strength `c_θ = c_gain · adaptation_index_θ` entering `x_eff = x − c·a`.

| Type | `adaptation_index` (median) | `c = 0.7 · index` | n cells | SFA in model |
|---|---|---|---|---|
| **Pyr**   | 0.0687 | 0.048 | 1539 | strongest |
| **Pvalb** | 0.0019 | 0.001 |  586 | ≈ none (fast-spiking) |
| **Sst**   | 0.0571 | 0.040 |  708 | moderate |
| **Vip**   | 0.0552 | 0.039 |  451 | moderate |

![SFA per cell type](figures/param_map_sfa.png)

Pvalb (fast-spiking) shows essentially no adaptation; Pyr adapts most, with Sst/Vip intermediate —
consistent with the intrinsic-physiology literature.

## 6. Short-term depression and facilitation (STD / STF) equations

The synapse is a rate Tsodyks–Markram pair (per `(pre neuron j, post type T)`):

$$
\frac{db_{j,T}}{dt} = \frac{1 - b_{j,T}}{\tau^{\mathrm{rec}}} - \frac{p_{j,T}\,b_{j,T}\,r_j}{\tau^{\mathrm{rel}}},
\qquad
\frac{dp_{j,T}}{dt} = \frac{p^0 - p_{j,T}}{\tau^{f}} + \kappa\,(1 - p_{j,T})\,r_j,
$$

with effective synaptic gain `g = (p/p0)·b` (rest: `b=1, p=p0 ⇒ g=1`). Every parameter above
(`τ_rec, τ_rel, p0, τ_f, κ`) is a `(pre, post)` matrix from §2.

## 7. Derived model STD/STF parameters

![Derived model synaptic parameters](figures/param_map_model_std_stf.png)

- **`p0`** (release prob) — highest for Pvalb→Pvalb/Sst (≈0.44–0.51, strongly depressing, high initial release), low for Pyr→Sst (floored at 0.05).
- **`κ`** (facilitation rate) — strongest for Vip→Vip, Pvalb→Vip, Vip→Pyr; the Vip row facilitates most.
- **`τ_f`** — mostly fast (floored at 0.05 s) except a few slow facilitating edges (Sst→Vip ≈ 0.53 s).
- **`τ_rec`** — 0.05–0.24 s; the Pvalb row recovers slowest (τ_rec ≈ 0.24 s).
- **`τ_rel`** — derived from depression amount; shorter (deeper depletion) where `ml_depression_amount` is large (e.g. Sst→Pyr), pinned to `τ_rel_ref = 0.5 s` where the amount is ≤ 0 (Vip row).

## 8. Clamps and fallbacks

Applied in `load_parameter_tables` (to the tables) and `get_params` (to the derived values):

| Quantity | Range clamp | NaN fallback (under-sampled) |
|---|---|---|
| `conn_prob_adj` | — | 0.05 |
| `psp_amplitude` | — | signed type-median |
| `ml_depression_tau` → `τ_rec` | [0.05, 5] s | 1.0 s |
| `ml_depression_amount` | [0, 0.95] (negatives → 0) | 0.3 |
| `ml_release_prob` → `p0` | [0.05, 0.95] | 0.2 |
| `ml_facilitation_tau` → `τ_f` | [0.05, 5] s | 1.0 s |
| `ml_facilitation_amount` → `κ` | ≥ 0 | 0.2 |
| `adaptation_index` | — | 0.03 |

In the current data, only the **clamps** fire (e.g. `p0` Pyr→Sst floored to 0.05; several `τ_rec`/`τ_f`
floored to 0.05 s; negative depression amounts set to 0) — no NaN fallbacks are needed.

## 9. Simplifications and open items

- **Fixed `τ_a`.** SFA uses one timescale, mapped from the adaptation *index* as a strength only;
  per-type `τ_a` (from ISI kinetics) is deferred (Option B).
- **`τ_rel` heuristic is provisional.** [`system_equations_with_stf.md`](system_equations_with_stf.md)
  argues this `ml_depression_amount → τ_rel` mapping should be *dropped* (depression already flows
  from `p0` and `τ_rec` once depletion is coupled to `p`). The code still uses the heuristic; that is
  the one mapping flagged for revision.
- **Pooled across layers, median statistics**; single STD/STF timescale each; `κ` uses the
  facilitation *amount* directly as a rate (the amount and its accumulation time — residual calcium
  buildup — enter only as this ratio because `r` is normalized).

## 10. Regenerating the figures

```matlab
setup_paths();            % or let param_visualizations add src/ to the path
param_visualizations();   % writes PNGs into docs/EquationsParametersDocs/figures/
```
