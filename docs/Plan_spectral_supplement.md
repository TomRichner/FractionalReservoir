# Plan: the spectral figures of 2026-11-20 as stages of the full run

Tom, 2026-11-20, after reading the November 20 figure report: "I'd like you to
make the changes (e.g. with more seeds) discussed so far, and integrate these
figures into the full run. Use fable 5.2 medium agents and opus 5.5 medium
agents to do it." And: "don't do a full run on this computer. I will do the
full run on a pc with more ram later."

So this plan turns the four explorations of 2026-11-19 and 20 (outlier
eigenvalues, frequency response, eigenvalue density, other kinds of outlier)
into registered opt-in stages and figure functions of the `mu7chain` bundles,
with the changes Tom asked for, built and tested on the MacBook at tiny and
`essential` size only. The full run is Tom's, on the workstation.

The lead decides what the plan leaves open and writes the decision under the
issue (FD-040 to FD-046) or in the worklog; Tom can reverse it afterwards.

## 1. What Tom decided (the record is in the manuscript's `docs/final_draft_to_do.md` and in FD-034 to FD-039)

| Report figure | Exploration output | Decision |
|---|---|---|
| 1 | Figure 3, four rows, switching run, K = 100, new timing | The paper's Figure 3 ("we will use the switching version") |
| 2 | `outliers_complex_plane_v3_alt` (E above, I below the curve) | "should go in the paper, or at least the supplement" |
| 3, 4 | ribbon insets | seen; not chosen |
| 5 | `outliers_summary` | redrawn against the theoretical radius; not chosen for the paper |
| 6 | `outliers_dose` | supplement, "with more seeds" |
| 7 | `outliers_slow_variables` | supplement |
| 8 | `outliers_timeseries` | not chosen |
| 9 | `outliers_dwell` | "could become a good figure": more seeds, medians, a much longer run |
| 10, 11 | `freqresp_network`, `freqresp_plane` | "I really like figures 10 and 11" |
| 12 | `freqresp_one_vs_multiple` | the phase-lead panel liked; left and middle panels not |
| 13 | `freqresp_single_neuron` | supplement |
| 14 to 28 | eigenvalue density (9), other outliers (6) | "much of the remaining figs below in the supplement too" |

Standing choices for these figures: complex-plane plots through
`complex_plane_axes` (axes through the origin, 1:1, one scale across compared
panels, a break in the real axis for a far eigenvalue, ticks outside the bulk
only); the bulk is the unscaled theoretical radius; eigenvalues beyond 103 %
of it are marked, no ringed markers for the rightmost one; no adaptation in
dark red `[0.75 0.10 0.15]` in every figure from the next full run on
(manuscript FR-004); E and I in the established E and I colors; no
explanatory text, counts or stacked titles on panels (sidecar and caption).

## 2. Lead's decisions where Tom left it open

- **Which of 14 to 28 are built:** all of them except the two ribbon variants,
  as supplementary candidates (not `in_paper`); choosing among them costs
  nothing once they are made by the run. The E-to-E x 1.5 copies of the
  density figures become a `network` option of the same functions.
- **Figure 12:** one new figure holding its phase-lead panel alone.
- **Figures 5 and 8:** not carried into the bundle; their numbers are in the
  sidecars of the others.
- **Outliers are judged against each state's own theoretical radius**, not
  the drawn median circle (the drawn circle inflates the single-timescale
  count in the reference network, 15.7 against 5.8 per state).
- **The far-left real eigenvalue** (E-to-E x 1.5 with adaptation, -59 to
  -42 1/s) is named in the caption and sidecar; no left-hand break.
- **Seeds and windows by run level** (section 4), sized so that the stages
  add a few hours, not a day, to the full run; the workers time them and the
  lead fixes the numbers.
- **Opt-in:** the stages are named by the `mu7chain` bundles (through one
  helper), not by other bundles; the regression bundle is unchanged.

## 3. Stages, figures, issues

Each stage: a plan function (no compute), a run function, a row in
`paper_opt_in_stages`, as FD-032 (`lyapunov_spectrum_n80`). Each figure: the
standard signature `out = fig_x(cfg)`, reads one run directory or errors
(`resolve_data_file`), `_table.md` sidecar, registry entry.

| Issue | Stage | Figures | Worker |
|---|---|---|---|
| FD-040 | `outlier_spectra`: reference and E-to-E x {1.1, 1.25, 1.5, 1.75}, three conditions, several seeds; frozen-gain spectra, gains, per-neuron states, null calibration and dense full Jacobian at calibration states; rightmost-eigenvalue series | (data only) | `matlab-fable-2` |
| FD-041 | `frequency_response`, from FD-040's calibration states | `fig_freqresp_network`, `_plane`, `_phase`, `_single_neuron` | `matlab-fable-3` |
| FD-042 | `eigenvalue_density`, from FD-040's spectra and gains | `fig_eigdens_radial`, `_gains`, `_measures`, `_plane`, `_theory_check` | `matlab-opus-1` |
| FD-043 | `outlier_structures`: the fifteen networks, the intended mode followed | `fig_ostruct_summary`, `_plane`, `_rotation`, `_h1` | `matlab-opus-2` |
| FD-044 | (none) no adaptation in dark red in `manuscript_style` | every registered figure | `matlab-opus-3` |
| FD-045 | (reads FD-040) | `fig_outlier_plane` (FD-036), `fig_outlier_dose` (FD-037), `fig_outlier_slow_variables` (FD-038), `fig_outlier_dwell` (FD-039) | after FD-040 merges |
| FD-046 | wiring: the four stages and all figures into the `mu7chain` bundles; the switching Figure 3 as the bundle default (Plan_Figure3_update, step 6); one `essential` run end to end, serial; cost table for Tom's full run | | after the rest |

Names: exploration code stays where it is until FD-046 (a job is still
running from it); stage code is COPIED into `src/analysis/` and
`src/figures/` under new names, so nothing on the MATLAB path is shadowed:
prefixes `ospec_` (FD-040), `fresp_` (FD-041), `eigdens_` (FD-042),
`ostruct_` (FD-043).

## 4. Sizes by run level (to be confirmed by timing)

| | tests | `essential` | `fastmedium` and above |
|---|---|---|---|
| network size | `n_override` (40 to 80) | 500 | 500 |
| seeds, outlier_spectra | 1 | 2 | 10 |
| measured window after 60 s settling | 2 s | 20 s | 200 s for reference and E-to-E x 1.5 (dwell); 20 s for the dose levels |
| calibration states | 3 | 12 | 12 |
| full-Jacobian series | off | reference and E-to-E x 1.5 | reference and E-to-E x 1.5 |
| seeds, outlier_structures | 1 | 1 | 5 |

## 5. The data contract of `outlier_spectra` (FD-040 owns it; the others read it through one loader of their own)

`<run_dir>/outlier_spectra/` holds one `.mat` per run and an index.

- `ospec_index.mat`: table `runs` with `network`, `condition`, `seed1`,
  `seed2`, `file`, `has_full_series`, `window_s`.
- Per run, struct `R`: everything `outlier_run` saves today (the series, the
  calibration-state spectra `ev_M`, null spectra, `R_eff`, `kappa`, measures
  inputs, dense full-Jacobian spectra and the continuation), plus
  `R.preset_name`, `R.net_args` (the name-value overrides that build the
  network), `R.rng_seeds`, and at the calibration states: `R.cal.t`,
  `R.cal.S` (state vectors, states x N), `R.cal.gain` (n x states, per
  neuron), `R.cal.u` (n x states, the argument of the nonlinearity relative
  to the neuron's own setpoint), `R.cal.rate`, `R.cal.cell_type`.
- `ospec_rebuild_model(R)` returns the model object of a run (for Jacobian
  blocks), without integrating.

FD-041 and FD-042 keep their reads of this in ONE function each, so a renamed
field is a one-file change at the merge.

## 6. Order, branches, gates

Wave 1, in parallel, each worker in its own git worktree and branch
(`../FractionalReservoir_wt_fd04N`, branch `fd04N-...`), because a worker is
still running in the main checkout: FD-040, FD-041, FD-042, FD-043, FD-044.
No wave-1 worker edits `paper_opt_in_stages.m`, `mu7chain_config.m`,
`mu7chain_AC_config.m`, `run_all_paper_analyses.m` or
`make_all_paper_figures.m`: each delivers its functions, its tests and the
registry rows as text in its report. Wave 2: FD-045, then FD-046 wires
everything on main.

Gates per branch: its own tests at tiny size; `checkcode` clean; the covering
existing tests named in its brief; no change to any existing saved output
with the stage off. The exact regression tiers are owed on the workstation
(with FD-032 to FD-034). The lead reviews each diff before merging.

Pool on the MacBook (budget 6): the ten-seed dose run holds 3; FD-040 may
hold 2; the others 0.

## 7. What Tom runs afterwards, on the workstation

The full run as in manuscript FR-003, with the bundle's new stages; the cost
of the new stages at `fastmedium` is in FD-046's table. Then the exact
regression tiers, `figure_sources.json` and `set_figs.py --sync`, and the
caption changes of Figure 3.

## 8. Open, for Tom

- Figure 2: main text or supplement; a figure of its own or the
  responsiveness panel of the closing figure.
- Figure 3: K = 100 (today) against the 150 recorded on 2026-11-17; the 40
  sampled states at 1, 2, ..., 40 s or at half-second offsets; this stage
  runs at 400 Hz while the uniform rerun is planned at 500 Hz.
- Which of the supplementary candidates stay.
