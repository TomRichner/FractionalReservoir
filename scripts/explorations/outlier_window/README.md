# outlier_window: what the 200 s window of outlier_spectra buys and costs

Question (Tom, 2026-11-23): the stage `outlier_spectra` gives `reference` and `EEx1p5`
a 200 s window (0.05 s frozen-gain step, 0.1 s full-Jacobian step) in the larger run
modes, for `fig_outlier_dwell` only; `fastmedium` spent 388 min on it on 12 workers.
Keep, shorten or drop the long window?

Data: a copy of `mu7chain_AC_fastmedium/outlier_spectra` (10 seeds) in
`data/explorations/outlier_window/fastmedium_copy/`, and of `mu7chain_AC_fast`'s stage
(2 seeds, 20 s) in `fast_copy/`. No stage or figure function is changed.

| script | does | writes |
|---|---|---|
| `ow_extract.m` | reads each run once; keeps the two series, the pass timings (`R.seconds`, `R.full.seconds`) and the window-wide figure quantities on the first 20, 50, 100, 200 s | `data/explorations/outlier_window/ow_series.mat` |
| `ow_draw_dwell.m` | `fig_outlier_dwell` unchanged, from both copies | `figs/explorations/outlier_window/dwell_200s_fastmedium/`, `dwell_20s_fast/` |
| `ow_analysis.m` | truncates each 200 s series to 20, 50, 100 s before `ospec_excursions`, `censored_survival`; the comparison figure; the cost split and wall-time estimates (list scheduling of rescaled run times in the stage's order) | `ow_truncation_table.md`, `ow_prefix_table.md`, `ow_cost_table.md`; `figs/.../ow_window_comparison.png` |
| `ow_time_passes.m` | times `ospec_snapshot` and `ospec_full_spectrum` at saved states, one thread | `ow_time_passes.md` |

## Findings (2026-11-23)

- **Full Jacobian** (solid curves): the reference with adaptation dwells a median
  0.8 s (single timescale) and 0.3 s (multiple), 1 to 5 % of excursions censored; the
  quartiles are the same at 20, 50, 100 and 200 s. Without adaptation the reference is
  beyond zero the whole time at every length (censored line at 1). E-to-E x 1.5 has no
  excursion with adaptation; without, 2 one-sample excursions in one seed at 200 s.
- **Frozen gains** (dashed): only two cells depend on the length: reference without
  adaptation (median > 20 s at 20 s, about 130 s at 200 s, half still censored) and
  reference with multiple timescales (pooled median 12.7 s at 20 s, 2.3 s from 100 s on).
- **Cost**: the 60 long runs take 89 % of the worker time; the 20 multiple-timescale long
  runs alone about 2 800 of 4 000 worker-min, nearly all in the full-Jacobian pass
  (3 to 4 s a state at 4 000 state variables, against 0.13 to 0.15 s for the frozen
  spectrum). Estimated stage wall time: 392 min at 200 s, 208 at 100, 133 at 50, 99 at
  20 s (extrapolated; the first reproduces the logged 388).
- Window-wide quantities of `fig_outlier_dose` and `fig_outlier_plane` (detachment,
  E fraction on the slope, E/I gain) change by far less than the spread across seeds
  between the first 20 s and 200 s.
