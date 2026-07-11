# Campagnola et al. 2022 — Data Types & Figure Map

Reference summary for parameterizing a cell-type-resolved rate reservoir (Pyr, Pvalb, Sst,
Vip) with SFA + STD + **STF**, extending `src/SRNNModel2.m`. This catalogs *what data the
paper provides* and *which figure / supplementary item holds each kind*, and flags which
quantitative matrices are locked in figures (so they need hand-digitization or an Allen API
pull) versus available as text.

**Source:** Campagnola, Seeman, Chartrand, et al., "Local connectivity and synaptic dynamics
in mouse and human neocortex," *Science* **375**, eabj5861 (2022).
DOI: [10.1126/science.abj5861](https://doi.org/10.1126/science.abj5861).
OCR of the PDF is in
`Campagnola et al. - 2022 - Local connectivity and synaptic dynamics in mouse and human neocortex.md`.

**Data portal (open):** <https://portal.brain-map.org/explore/connectivity/synaptic-physiology>
Python client: `aisynphys` (Allen Institute synaptic physiology database).

---

## 1. Scope of the dataset

- **Technique:** up to 8-cell simultaneous whole-cell patch-clamp (multipatch) in acute
  slices; mostly current-clamp with a subset in voltage-clamp. Intralaminar connections only
  (cells < ~200 µm apart).
- **Counts:** 1931 experiments. **23,620 potential connections probed** (mouse 20,949; human
  2671); **1731 connected by chemical synapses** (overall 7.3%). Mouse: 1526 connections
  (7.3%); human: 205 (7.7%).
- **Regions:** mouse primary visual cortex (VISp, 1715 experiments, adult ~46 d) and human
  frontotemporal neurosurgical tissue (216 experiments).
- **Probing protocol for dynamics:** stimulus trains at **10–200 Hz**; recovery tested at
  delays of **125, 250, 500, 1000, 2000, 4000 ms**.

---

## 2. Cell subclasses & layers

The four subclasses we care about:

| Subclass | Class | Notes |
|---|---|---|
| **Pyr** (excitatory) | E | further split by layer/projection (below) |
| **Pvalb** | I | sampled across all targeted layers; fast, strong, depressing, low-variability |
| **Sst** | I | sampled across all layers; slow kinetics; high-variability; targets others broadly, avoids self (mostly) |
| **Vip** | I | sampled across all layers; disinhibitory (targets Sst) |

Excitatory subclasses are **layer / projection specific** (inhibitory subclasses span all
layers):

- **L4:** Nr5a1, Rorb
- **L5 ET** (extratelencephalic): Sim1 / mscRE4-FlpO
- **L5 IT** (intratelencephalic): Tlx3
- **L6 CT** (corticothalamic): Ntsr1
- **Human L2/3** further resolved into L2, L3a, L3b pyramidal subclasses (depth gradient).

Layers covered: **L2/3, L4, L5, L6** (mouse); supragranular focus (L2, L3a, L3b) for human.

---

## 3. Data types × figures (core map)

"Format" = how it appears relative to the OCR: **text** (numbers recoverable from the .md),
**image** (locked in a figure raster — the OCR did not capture it), or **model/DB** (per-
connection parameters best pulled from the portal/`aisynphys`).

| Data type | What it is | Figure(s) / Table(s) | Format |
|---|---|---|---|
| **Connection-probability matrix (mouse)** | `p_max` per pre→post subclass element, layer-resolved (unified-model adjusted) | **Fig 2B**; distance fits Fig 2C–F; layer detail **fig S1A** | **image — needs digitize/API** |
| Distance dependence of connectivity | Gaussian `p_max` + lateral spread σ per class; slice/depth/detection corrections | Fig 2A, C–F; fig S1B; **fig S2** (A–E corrections) | mostly image (some σ, p_max in text) |
| Reciprocal / bidirectional rate | reciprocity normalized to random expectation (3–5× enriched) | Fig 2A bottom row | image |
| Electrical (gap-junction) connectivity | connection prob, coupling coefficient, junctional conductance, input resistance | **fig S4A–D** | image (key nS values in text) |
| **Synaptic strength & kinetics matrices** | PSP **latency, rise time, decay tau, resting-state amplitude, 90th-pct amplitude**, per E–I subclass element | **Fig 3A–E** | **image — needs digitize/API** |
| **Short-term plasticity (STP) matrices** | STP induction ratio (facilitation ↔ depression), recovery @ 250 ms, resting-state aCV (variability) | **Fig 4B–D** | **image — needs digitize/API** |
| STP vs frequency & recovery curves | train-induced STP at 10/20/50/100 Hz; recovery at 125–4000 ms | Fig 4E–F | image |
| Stochastic release / STP generative model | per-connection **quantal params: release probability, # release sites, facilitation τ, depression mechanism** (vesicle-depletion vs release-independent), for 1196 connections (1035 mouse, 161 human) | text + **Fig 6**; figs S6–S7; **table S2** (feature rankings) | **model/DB** (portal/`aisynphys`) |
| Dimensionality-reduction organization | sPCA → UMAP of connection model params; strength vs variability axes; E vs I clustering | Fig 6A–D | image (qualitative) |
| Human intralaminar data | conn prob, kinetics, STP, L2/3 depth gradients, disynaptic inhibition | **Fig 5A–G**; fig S5 | image |
| Canonical circuit diagrams | Pvalb/Sst/Vip motifs; L2/3 vs L5 differences; facilitating-vs-depressing "dynamic modes" | **Fig 7A–D** | schematic (qualitative) |
| Intrinsic ephys & morphology | input resistance, firing/spiking, dendrite type, axon length, biocytin fills | Fig 1B–C; fig S3 (Sst/Martinotti); UMAP fig S3A | image (a few values in text) |
| **Spike-frequency adaptation (SFA)** & intrinsic firing | ISI-based adaptation index, F–I slope, rheobase, firing rate, membrane τ — per cell, from 1-s current steps | Fig 1C (protocol); underlies UMAP in Fig 5E–F | **DB `intrinsic` table** (see §5) |
| Transgenic lines / sampling | driver/reporter lines, ages, cell/experiment counts | Fig 1A; **table S1** | table S1 **not in OCR** |
| STP model validation | recorded vs simulated strength/STP/variability correlations | figs S6–S7 | image |

**The three matrices we most need for the model** are all image-only in the OCR:
**Fig 2B** (connection probability), **Fig 3** (strength + kinetics), **Fig 4B–D** (STP:
facilitation/depression, recovery, variability).

---

## 4. What we can get quantitatively, and how

### 4a. Directly usable numbers already in the text/OCR

These are representative point values quoted in Results (not the full matrices). Useful for
sanity-checks and seeding, not for a complete parameterization. Values are median [IQR] or
[95% CI] as reported.

**Connection probability (`p_max`, mouse):**
- Class-level peaks: E→E ≈ 5%; I→I ≈ 12%; E→I ≈ 12%; I→E ≈ 15%.
- Pyr→Pvalb: L2/3 = 0.79 [0.57, 0.99] vs L5 IT = 0.22 [0.12, 0.35].
- E→Vip = 0.38 [0.21, 0.60] and Vip→E = 0.11 [0.03, 0.26] in L2/3; rare/absent deeper.
- Human recurrent L4 E→E ≈ 0.01 [0, 0.04] (1/145) vs mouse 0.22 [0.16, 0.28] (44/452).
- Lateral spread σ: inhibitory chemical ≈ 127 µm; electrical ≈ 74 µm; human ≈ 130 µm;
  mouse within-class ≈ 125 µm.

**Strength & kinetics:**
- Inhibitory: median latency 1.07 ms, slow kinetics, strong PSPs (Pvalb submillisecond,
  e.g. Pvalb→L6 E = 0.94 [0.86, 0.99] ms).
- Excitatory: median latency 1.49 ms; E→I rise 2.75 ms vs E→E rise 3.88 ms.
- Sst slow kinetics: Sst→Vip rise 7.02 ms, decay 50.81 ms; Sst→L5 IT rise 7.46 ms, decay 58.59 ms.
- E→Pvalb resting amp 0.41 [0.22, 0.79] mV vs E→Sst 0.08 [0.05, 0.25] mV; L2/3 Pyr→Vip
  0.45 [0.26, 0.84] mV; biggest single PSP 9.51 mV (L2/3 Pyr→Sst).
- L4 Pyr→Pvalb rise 1.45 ms, decay 6.92 ms; L2/3 Pyr→Sst rise 3.91 ms, decay 19.41 ms.
- Amplitude class contrast: E→I 90th-pct 0.73 mV vs E→E 0.34 mV.

**STP (the sign/direction is the key qualitative fact):**
- **Facilitating:** E→Sst (strongly, high resting variability then reliable when induced);
  **Sst→Vip** (highest facilitation in dataset, 0.28 [0.0, 0.44]); Vip→Sst and recurrent Vip
  (weakly). These are the connections that motivate adding **STF**.
- **Depressing:** recurrent E→E (esp. L5 ET→L5 ET; deepens with frequency); E→Pvalb (mostly);
  Pvalb→anything (strong depression, still depressed at 125 ms).
- **Pseudolinear:** subset of L2/3 E→Pvalb.
- Variability (log aCV): Pvalb resting ≈ −1.01 → induced ≈ −0.27; E→Sst resting ≈ +1.66
  (high); Sst↔Vip resting ≈ −0.03.

### 4b. Locked in figures — need digitization

The **full subclass × subclass matrices** (the thing we actually need to build W and its
per-type dynamics) are figure rasters, not in the OCR:
- **Fig 2B** — connection-probability matrix.
- **Fig 3A–E** — latency / rise / decay / resting amp / 90th-pct amp matrices.
- **Fig 4B–D** — STP induction, recovery, resting aCV matrices.
- **fig S1A** — layer-resolved connectivity detail.

### 4c. Best obtained via the Allen portal / `aisynphys` API — **recommended**

The underlying per-connection database is open and exposes, programmatically, everything the
figure matrices summarize — connection probability with corrections, strength, kinetics, STP
metrics, **and the stochastic-release model parameters** (release probability, number of
release sites, facilitation/depression amount + time constants). This is strictly better than
pixel-digitizing figures: it gives per-connection distributions, not one median per element,
and lets us pool cells into whatever cell classes the model needs.

- Portal: <https://portal.brain-map.org/explore/connectivity/synaptic-physiology>
- Repo studied locally: `/Users/richner.thomas/Desktop/local_code/aisynphys`
  (branch `origin/figure-notebooks` holds the exact 2022 manuscript notebooks; `doc/` holds
  tutorial reproductions using identical DB fields/helpers).

**Access is a plain public SQLite download — no auth / API key.**
`SynphysDatabase.load_current('small')` fetches a manifest from GitHub, downloads the newest
`small` sqlite (~73 MB) into `~/ai_synphys_cache/database/`, and opens it read-only. Sizes are
`small` / `medium` / `full`; larger tiers only add deferred raw-waveform/array columns —
**every scalar we need is in `small`.**

**Planned granularity for the model:** 4 subclasses (Pyr, Pvalb, Sst, Vip) **pooled across
layers** → clean 4×4 matrices. (The DB also supports the layer-resolved 16–17 classes if we
later want L2/3-vs-L5 detail.)

**Deps needed** (not yet in `.venv`): `sqlalchemy`, `pandas`, `neuroanalysis` (git package),
plus `scikit-learn` + `umap` only for the Fig 6 sPCA/UMAP reduction. `numpy`/`scipy`/
`matplotlib` already present.

**General extraction recipe** (from `doc/short_term_plasticity.ipynb`):

```python
import numpy as np
from aisynphys.database import SynphysDatabase
from aisynphys.cell_class import CellClass, classify_pair_dataframe

db = SynphysDatabase.load_current('small')
pairs = db.pair_query(experiment_type='standard multipatch', species='mouse',
                      synapse=True, preload=['cell', 'synapse']).dataframe()

# 4 subclasses pooled across layers
cell_classes = {
    'pyr':  CellClass(dendrite_type='spiny'),
    'pv':   CellClass(cre_type='pvalb'),
    'sst':  CellClass(cre_type='sst'),
    'vip':  CellClass(cre_type='vip'),
}
classify_pair_dataframe(cell_classes.values(), pairs)   # adds pre_class / post_class cols
M_stp = pairs.pivot_table('dynamics.stp_induction_50hz', 'pre_class', 'post_class', aggfunc=np.mean)
```

`aisynphys/ui/notebook.py::cell_class_matrix(pre_classes, post_classes, metric, db, ...)` is
the generic per-metric pre×post matrix builder. Connection probability uses the separate
`measure_connectivity(pair_groups, sigma=100e-6, dist_measure='lateral_distance')` path from
`aisynphys.connectivity` (see `connectivity.ipynb`).

**Table → column map for each matrix** (see also §3):

| Matrix | Table.column(s) | Tutorial notebook |
|---|---|---|
| Connection probability | `measure_connectivity()` → `connection_probability`, `adjusted_connectivity` | `connectivity.ipynb` |
| Strength | `synapse.psp_amplitude`, `dynamics.pulse_amp_90th_percentile` | `synaptic_strength.ipynb` |
| Kinetics | `synapse.latency`, `synapse.psc_rise_time`, `synapse.psc_decay_tau` | `synaptic_kinetics.ipynb` |
| STP | `dynamics.stp_induction_50hz`, `dynamics.stp_recovery_250ms`, `dynamics.variability_resting_state` | `short_term_plasticity.ipynb` |
| Release model (STF/STD) | `synapse_model.ml_base_release_probability`, `ml_n_release_sites`, `ml_facilitation_amount`, `ml_facilitation_tau`, `ml_depression_amount`, `ml_depression_tau`, `ml_mini_amplitude` | `tutorial.ipynb`, `aisynphys/synapse_types.py` |

**Key for the STF extension:** the `synapse_model.ml_*` columns give per-connection
*facilitation amount + facilitation τ* and *depression amount + depression τ* directly — the
exact quantities to parameterize the new STF state variable and existing STD, indexed by
(pre-type, post-type). Alignment rule to preserve: excitatory dynamics index by **post**-type,
inhibitory dynamics by **pre**-type (§5).

---

## 5. Spike-frequency adaptation (SFA) — intrinsic data

**Yes, SFA data is available**, but it lives in a different place from the synaptic matrices:
it is a **per-cell intrinsic property** measured from the 1-second current-step recordings
(Fig 1C), stored in the DB **`intrinsic`** table (`aisynphys/database/schema/intrinsic.py`,
accessor `Cell.intrinsic`). It is *not* connection-resolved — one row per cell, keyed by
`cell_id`. Because every cell carries `cre_type` / `cell_class`, `adaptation_index` can be
pooled by the same 4 subclasses (Pyr, Pvalb, Sst, Vip) as the synaptic matrices.

**SFA-relevant `intrinsic` columns:**

| Column | Meaning | Use |
|---|---|---|
| **`adaptation_index`** | ratio of consecutive ISIs, averaged across sweeps | **primary SFA measure** — degree of firing slowdown during a sustained step |
| **`isi_adapt_ratio`** | ISI on 5th vs 1st spike | complementary SFA measure |
| `isi_cv` | CV of the ISI distribution | firing regularity |
| `fi_slope`, `rheobase` | F–I slope, threshold current | excitability / gain |
| `firing_rate_rheo`, `firing_rate_40pa` | mean firing rate at rheobase and +40 pA | firing dynamics |
| `latency_rheo`, `latency_40pa` | first-spike latency at those currents | onset dynamics |
| `tau`, `input_resistance`, `input_resistance_ss`, `sag` | membrane τ, Rin, H-current sag | passive / subthreshold |
| `upstroke_adapt_ratio`, `downstroke_adapt_ratio`, `width_adapt_ratio`, `threshold_v_adapt_ratio` | spike-shape adaptation, 5th vs 1st spike | secondary adaptation signatures |

**Access recipe (per-cell, pooled by subclass):**

```python
from aisynphys.database import SynphysDatabase
db = SynphysDatabase.load_current('small')
q = db.query(db.Cell.cre_type, db.Cell.cell_class,
             db.Intrinsic.adaptation_index, db.Intrinsic.isi_adapt_ratio,
             db.Intrinsic.fi_slope, db.Intrinsic.tau) \
      .join(db.Intrinsic, db.Intrinsic.cell_id == db.Cell.id)
df = q.dataframe()          # then groupby cre_type / cell_class → per-subclass SFA
```

Alternatively, `pair_query(...)` decorates the query with `pre_intrinsic` / `post_intrinsic`
aliases if you want intrinsic features alongside a connection query.

### Mapping the index to model parameters: strength vs. timescale

`adaptation_index` / `isi_adapt_ratio` are **phenomenological, ISI-based** measures — they say
*which cell types adapt more, and by how much*, but they do **not** directly give the SFA
**time constant** `τ_a` that `SRNNModel2`'s `a` variables represent. A single scalar cannot
independently set both the adaptation strength (`c`) and `τ_a`, because the index conflates
them:

- At steady input, `da/dt = 0 ⟹ a = r`, so the adapted rate solves `r = b·φ(x − c·n_a·r)`.
  **The steady-state amount of adaptation depends only on `c`, not on `τ_a`.** `τ_a` only sets
  *how fast* the network reaches that steady state.
- Over a long train (≫ τ_a) the index is dominated by `c`; over the early window the index is
  sensitive to both. So the index most naturally constrains **strength**, not timescale.

**Option A — index → strength `c` (CHOSEN, to get started).** Hold `τ_a` fixed at a
physiological value (the model's existing SFA τ), and set the per-subclass adaptation strength
(`c_E` / `c_I`, or a new per-subclass factor) from the relative `adaptation_index` across
Pyr/Pvalb/Sst/Vip. Clean, monotonic, one-line mapping; matches what the index actually
measures. (Expected ordering: Pvalb weak, Pyr/Sst stronger.)

**Option B — index → per-type `τ_a` (KEEP IN MIND for later; design for it now).** Possible via
a one-time calibration: fix `c` and the firing rate, sweep `τ_a`, simulate the *model's own*
adaptation index under the paper's 1-s-step protocol to get `index(τ_a)`, then invert. Caveats:
must verify monotonicity over the τ_a range, and the result is entangled with the fixed `c` /
firing rate. A *true* per-type τ_a is better fit from raw ISI-vs-time trajectories in the
**`full`** DB tier (deferred array columns), not `small`.

**Design implication:** parameterize SFA so `c` (or `c_E`/`c_I`) **and** `τ_a` are both
**per-subclass fields** even though Option A only varies `c` at first — that way switching to
Option B later means populating an existing `τ_a` vector, not reworking the state/param layout.

Also note: intrinsic characterization covered a **subset** of cells, so per-subclass SFA
coverage will be sparser than the synaptic matrices — check `n` per group before relying on it.

---

## 6. Relevance to the `SRNNModel2` extension (orientation only)

How these data map onto the planned additions (design is a later task):

| Data | Model use |
|---|---|
| Connection-probability matrix (Fig 2B / DB) | cell-type-**blocked W construction** — per (pre-type, post-type) connection probability, replacing the single RMT indegree |
| Strength matrix (Fig 3D–E / DB) | per-type **weight scaling** of nonzero W entries |
| STP: facilitating elements (E→Sst, Sst→Vip, Vip→Sst/Vip) | drives the **new STF state variable** (`b`-like facilitation) alongside existing STD |
| STP: depressing elements (E→E, Pvalb→*, E→Pvalb) | maps to existing **STD** mechanism |
| SFA (`intrinsic.adaptation_index`, §5) | SFA mechanism already present in `SRNNModel2` (the `a` variables); use per-subclass `adaptation_index` to scale SFA magnitude (`c_E`/`c_I`) per cell type |
| Stochastic-release variability (Fig 6, table S2) | optional noise term per connection type (variability is the dataset's strongest subclass discriminator) |

Key biological facts to preserve in the model:
- **Excitatory synapse dynamics align with the _postsynaptic_ subclass**; **inhibitory
  synapse dynamics align with the _presynaptic_ subclass.** (This dictates how STP params are
  indexed by pre/post type.)
- The **Sst→Vip** facilitating connection (often overlooked vs the Vip→Sst disinhibitory
  pathway) is a headline motivation for adding STF.
- "Dynamic modes" (Fig 7D): depressing Pyr↔Pvalb interactions dominate at activity onset;
  facilitating Pyr/Sst/Vip interactions dominate during sustained activity.

---

## 7. Open follow-ups

- **DECIDED:** use the `aisynphys` API pull, not figure digitization (§4c). Granularity =
  4 subclasses pooled across layers.
- **Next (deferred):** set up the env (add `sqlalchemy`, `pandas`, `neuroanalysis` to `.venv`),
  download the `small` DB, and extract the connectivity / strength / kinetics / STP /
  release-model matrices to CSV + `.mat` for use by `SRNNModel2`.
- **SFA (§5):** extract per-subclass `intrinsic.adaptation_index` (+ `isi_adapt_ratio`,
  `fi_slope`, `tau`) pooled over the 4 subclasses. **Decided: Option A** (index → strength `c`,
  fixed `τ_a`) to start; **design SFA params so `c` and `τ_a` are both per-subclass fields** so
  Option B (per-type `τ_a` via calibration / raw-ISI fit) is a drop-in later.
- `table S1` (transgenic lines) and figs S1–S14 are in the Supplementary Materials
  (`science.org/doi/10.1126/science.abj5861`), not in the main-text OCR — retrieve if needed.
- Design the cell-typed W construction and the STF state variable in `SRNNModel2`.
