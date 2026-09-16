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
| 1 | Native archived Sompolinsky axes from the explicitly selected `sfaEI_fast` introductory FIGs, with unchanged trace/eigenvalue data. Four analytic lower panels: SFA sigmoid / shifted disks / STD sigmoid / shrinking disks. Exact paired colors, 14-point fonts, axis widths1.0. Disks are conceptual approximations. |
| 2 | Current baseline trace pixels retained under calibrated native axes; six rows A–F, titles No Adaptation/STA/MTA, requested limits and one 10-s scale bar. Original trace colors and 5-s finite-LLE start remain. Future runs use numeric saved trajectories and exact per-neuron E/I palettes. |
| 3 | Native saved-data eigenvalue occupancy and simplified switching-input traces. Input is eight fixed actual neurons. Local positive QR sum is `sum(max(local_rates,0),2)/log(2)`; percentage uses the same samples and strict positivity over the saved second-half window. Broader entropy interpretation remains deferred. |
| 4 | Native median/IQR LLE sensitivity and rate–LLE panels; inherited fixed display clipping remains. |
| 5 | Native non-normal margin and **only active-direction gain**, all five curves. Different colors, dashes, markers; E/I labels use “Excite E / inhibit I” and “Excite E / excite I”. Noise is RMS gain. All three condition panels share logarithmic gain limits. The current saved horizon remains **1 s**. Margin and gain sample different stages. |
| 6 | Native timescale LLE, leading-vector fractions, and three memory panels. Both timescale axes are linear; source statistics are unchanged. |
| 7 | Explicit sibling manuscript PyTorch image when present, visibly provisional. This old experiment does not establish the intended three-condition learning claim. |
| 8 | Archived model PSD with explicitly unavailable human SOZ panel on this laptop. Uniform DC model input is distinct from clinical 2-Hz stimulation. |

The archived illustrative sources contain PNGs only, not saved MATLAB FIGs or raw illustrative trajectories. Groups 2/7/8 include raster panels; saving them inside a FIG does not make their data editable. Groups 1 and 3–6 are native graphics. For group1 the figure-only config explicitly names the older native Sompolinsky archive, whose simulation settings match the mu7 intro; the intervening source change only combined its two rows. Those two FIG files must also be copied to a second computer to reproduce this exact native replot.

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

- `run_paper_illustrations`: reference network conditions on paired seeds `[1 2]`; uniform positive input amplitude 0.5, onset 15 s midway through displayed `[0,30]` s; simulation starts at -15 s for settling. Save 4 fixed neurons per type, raw input/state/rate/synaptic output/SFA feedback/STD product. The network illustration now saves top-K K = 15 leading local/finite rates, with 0.05-s segments, accumulation from t = 0 and 10-s negative-time alignment. The first finite value is displayed at the first completed-segment endpoint t = 0.05 s. Single-neuron mechanisms still skip LLE estimation.
- Supplementary mechanism stage: main-preset-derived single E neuron, `n=1`, zero recurrent weights, noise/heterogeneity off. Columns are no adaptation / **1TS SFA only** / **1TS STD only**. The old `single_neuron_dualStd` figure used three SFA and two STD timescales and four columns; it is replaced in the revised registry by the saved-data figure.
- Transient gain: `cfg.transient_gain_horizon_s=5`, `cfg.transient_gain_duration_s=40` (60 in production). These are explicit knobs. Sample-start bounds are validated before any computation, requiring `T/2 + 5 < T - horizon`. The saved horizon is plotted directly; current 1-s data are never stretched.

The new illustrative and longer-horizon analyses have **not** run on this laptop. Full numerical/statistical results and solver verdicts still require the stronger-computer run. The existing fast numerical verdict remains mixed.

## Verification

MATLAB R2026a through MCP: all eight current grouped exports succeeded. Visual checks covered panel composition, all five direction curves, and timescale axes. `test_mu7revised_bundle` passed without simulations: baseline equivalence, memory noise-only difference, independent fast/medium paths, one-neuron zero weights and timescale counts, midpoint stimulus values (also checked in the built n=500 model with negative-time burn-in), and too-short-horizon rejection. Saved group 5 was also checked for exactly five curves in each of three panels and identical log limits. Code Analyzer found no new functional issues; the transient stage retains its pre-existing parfor broadcast warning. Full new-stage numerical execution is intentionally deferred.

## Figure 1 mechanism revision

The lower panels use paired black-to-orange SFA curves/disks and blue-to-teal STD curves/disks. STD uses the actual main-preset E→E MTS recovery/release ratios at representative raw rate0.25: each steady depression factor is1/3, and their product is1/9. The five frozen total attenuation levels span1 to1/9, so the strongest illustrated sigmoid remains nonzero. These are conceptual frozen-factor curves, not a self-consistent steady-state F-I relation or an exact transformation of the active Jacobian. See the generated Figure1 provenance for parameters and the authoritative `docs/EquationsParametersDocs/Equations_stability_paper.md` for dynamics. No simulations were run for this revision.

## Figure 2 formatting and current-source limits

The current mu7 archive contains only a PNG for this example. Its exact plotted trace pixels are calibrated into native axes; no other run or stage supplies replacement trajectories. Existing neutral axis/text pixels and colored legend strokes are removed, while remaining trajectory artwork is retained. Pixels overwritten by the original transparent legend glyphs cannot be recovered. The fixed source checksum prevents applying that calibration to a different image. Current finite-LLE curves still begin at 5 s and no midpoint step has been added to the existing data.

New labels/axes use 14-point fonts, with 20-point column titles, one E/I legend, row letters A–F and gray column dividers. The 10-second scale bar occupies the lower-left lambda panel at x = 5…15, y = -4.8. Per-neuron colors cannot be exactly reassigned from flattened antialiased raster overlaps; the current E/I shades are retained. The future numeric plotter calls the standalone `excitatory_colormap` and `inhibitory_colormap` helpers (the same palettes used privately by `SRNNModel2`), preserves neuron ordering across rows and draws E over I.

For the next `mu7revised` analysis only, `illustration_lya_start_s=0` and `illustration_lya_warmup_s=10`: alignment over [-10,0] s occurs within the simulation's [-15,30] s range. Finite rates are stored at segment endpoints (first 0.05 s). This support was statically checked and model construction validated, but no new illustrative simulations or estimates were produced locally. Input trajectories remain saved; the six main rows are x, raw rate, synaptic output, SFA feedback, STD product andlambda.

## Figure 3 saved-data revision

Four native rows now show eigenvalue density, 24 actual fixed input-neuron traces (12 E and 12 I), accumulated leading exponent, and local positive QR-rate sum. Condition headings use the standard black/gold/blue colors and 16-point type; all other fonts are 14. Figure 2 headings now use those condition colors too. Rows 3–4 share [0,60] s axes. The current source is **K=30**, accumulation **[30,60] s**, with saved segment-start timestamps 30…59.95 s. The figure explicitly states this; requested K=100 and early accumulation cannot be recovered by replotting.

The future bundle explicitly sets `local_lyapunov_K=100` even in fast mode, `local_lyapunov_accumulation_start_s=1`, `local_lyapunov_warmup_s=15`, and `local_lyapunov_simulation_start_s=-15`. This provides the full alignment [-14,1] s before accumulation [1,60] s. Saved QR times retain the stage's segment-start convention (first completed [1,1.05] segment is indexed at 1). Legacy callers retain K by mode and [T/2,T] accumulation. No new switching simulation was run locally. Local positive sum remains a provisional local h_KS proxy, not an established invariant KS entropy.
