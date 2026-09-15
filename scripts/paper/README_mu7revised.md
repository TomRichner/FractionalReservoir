# Grouped mu7 figures and the next analysis bundle

## Current data: figures only

From the MATLAB desktop with this repository as the current folder:

```matlab
setup_paths();
results = make_all_paper_figures(sfaEI_mu7_grouped_figures_config());
```

Equivalently, run `sfaEI_mu7_grouped_figures_run.m`. This config calls no analysis or inline-simulation figure functions. It reads `data/sfaEI_mu7_fast` and writes `figs/sfaEI_mu7_fast_grouped`, with one PNG, editable MATLAB FIG, and provenance note per group, plus the master manifest and report. Existing source-run data and figure files are not modified. Explicit repository-relative paths work when the data are copied from another computer; no archived Windows path is used as a fallback.

| Group | Current source and limitations |
|---|---|
| 1 | Archived introductory PNG; native analytic F-I axes and conceptual disk schematic. The disk is an approximation, not the full adaptive Jacobian. |
| 2 | Archived baseline PNG, visibly labeled midpoint-step pending. No step has been fabricated. Future runs prefer saved representative dynamics. |
| 3 | Native saved-data eigenvalue occupancy and simplified switching-input traces. Input is eight fixed actual neurons. Local positive QR sum is `sum(max(local_rates,0),2)/log(2)`; percentage uses the same samples and strict positivity over the saved second-half window. Broader entropy interpretation remains deferred. |
| 4 | Native median/IQR LLE sensitivity and rate–LLE panels; inherited fixed display clipping remains. |
| 5 | Native non-normal margin and **only active-direction gain**, all five curves. Different colors, dashes, markers; E/I labels use “Excite E / inhibit I” and “Excite E / excite I”. Noise is RMS gain. All three condition panels share logarithmic gain limits. The current saved horizon remains **1 s**. Margin and gain sample different stages. |
| 6 | Native timescale LLE, leading-vector fractions, and three memory panels. Both timescale axes are linear; source statistics are unchanged. |
| 7 | Explicit sibling manuscript PyTorch image when present, visibly provisional. This old experiment does not establish the intended three-condition learning claim. |
| 8 | Archived model PSD with explicitly unavailable human SOZ panel on this laptop. Uniform DC model input is distinct from clinical 2-Hz stimulation. |

The archived illustrative sources contain PNGs only, not saved MATLAB FIGs or raw illustrative trajectories. Thus groups 1/2/7/8 include raster panels; saving them inside a FIG does not make their data editable. Groups 3–6 are native data graphics.

## Future computation: not executed locally

`mu7revised_config` is independent of prior configs. Its baseline preset is an explicit copy of the mu7 physical parameters and conditions, under a new name. Its noise-free memory, switching-input and single-neuron derivatives chain only within this bundle.

```matlab
setup_paths();
cfg = mu7revised_config();             % FAST default
% cfg = mu7revised_config('medium');   % longer runs / more trials
run_dir = run_all_paper_analyses(cfg);
results = make_all_paper_figures(cfg);
```

Alternatively, edit `run_mode` in `mu7revised_run.m` and press Run. Outputs are separate: `data/mu7revised_fast` / `figs/mu7revised_fast`, or their `_medium` counterparts. The analysis master refuses a nonempty data target. The legacy source registry is retained and grouped figures are appended after it. External panel paths can be supplied in `grouped_figure_registry`; absent human/learning inputs remain explicit pending panels.

New analysis support:

- `run_paper_illustrations`: reference network conditions on paired seeds `[1 2]`; uniform positive input amplitude 0.5, onset 15 s midway through displayed `[0,30]` s; simulation starts at -15 s for settling. Save 4 fixed neurons per type, raw input/state/rate/synaptic output/SFA feedback/STD product. No LLE estimation in this illustration.
- Supplementary mechanism stage: main-preset-derived single E neuron, `n=1`, zero recurrent weights, noise/heterogeneity off. Columns are no adaptation / **1TS SFA only** / **1TS STD only**. The old `single_neuron_dualStd` figure used three SFA and two STD timescales and four columns; it is replaced in the revised registry by the saved-data figure.
- Transient gain: `cfg.transient_gain_horizon_s=5`, `cfg.transient_gain_duration_s=40` (60 in production). These are explicit knobs. Sample-start bounds are validated before any computation, requiring `T/2 + 5 < T - horizon`. The saved horizon is plotted directly; current 1-s data are never stretched.

The new illustrative and longer-horizon analyses have **not** run on this laptop. Full numerical/statistical results and solver verdicts still require the stronger-computer run. The existing fast numerical verdict remains mixed.

## Verification

MATLAB R2026a through MCP: all eight current grouped exports succeeded. Visual checks covered panel composition, all five direction curves, and timescale axes. `test_mu7revised_bundle` passed without simulations: baseline equivalence, memory noise-only difference, independent fast/medium paths, one-neuron zero weights and timescale counts, midpoint stimulus values (also checked in the built n=500 model with negative-time burn-in), and too-short-horizon rejection. Saved group 5 was also checked for exactly five curves in each of three panels and identical log limits. Code Analyzer found no new functional issues; the transient stage retains its pre-existing parfor broadcast warning. Full new-stage numerical execution is intentionally deferred.
