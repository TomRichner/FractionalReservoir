# Reproducing the steps / top-30 example

`example_mu7chain_steps_top30.m` in this folder draws the mu7chain network
under four adaptation conditions for 40 s: 10 s with no stimulus, a new
random step every 5 s from 10 to 30 s, then 10 s with no stimulus. The
bottom row of the time-series figure shows the top 30 Lyapunov exponents.
This note says how to get the same three figures and the same numbers
again.

## The commit

| | |
|---|---|
| Repository | `FractionalReservoir` (`git@github.com:TomRichner/FractionalReservoir.git`) |
| Branch | `main` |
| Commit | `dbca1ce2ab82709d733ec8cde9f889be17612f96` (2026-11-12 14:56 -06:00) |
| MATLAB | R2026a Update 4, Windows 64-bit, 14 computational threads |

That is the commit the example was run at. When this note was written it
existed only in the local repository: `origin/main` was at `0bc35ca`, which
has neither this script nor the two options of `fig_example_timeseries` it
uses (`lambda_K`, `lambda_link`). To reproduce on another machine the
commit has to be pushed first.

## Steps

1. Check out the commit (or any later one; see "If the numbers differ").

   ```
   git checkout dbca1ce
   ```

2. In MATLAB, with the current folder at the repository root:

   ```matlab
   setup_paths
   example_mu7chain_steps_top30
   ```

   `setup_paths` is needed once per session; the script has no path code
   of its own. No parallel pool is needed. It takes about three and a half
   minutes (214 s on the machine above) and opens three figure windows:
   the time series with the top-30 exponents, the stimulus, and the
   eigenvalues of W.

3. Compare what it prints with "What a correct run gives" below.

Nothing is written to disk. To save the figures as well, set
`save_figs = true` in the script's settings block **and give `out_dir` a
folder of your choice**. With `out_dir` left empty the figure function
saves beside its own source file in `src/figures/`, which is not where
output belongs.

## What defines the example

Everything is in the settings block at the top of the script.

| Setting | Value | Meaning |
|---|---|---|
| `preset_name` | `celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain` | The mu7chain bundle's network: n = 500, half E and half I, `mu_tilde_relative` 7 on every route, per-neuron setpoints, SFA ladders and depression parameters (`tau_a_spread` and `std_spread` 0.25), Wiener noise 0.025. Four conditions: no adaptation, single-timescale, chained multiple-timescale, product multiple-timescale. |
| `rng_seeds` | `[7 8]` | One seed pair fixes W, the setpoints, the SFA and depression draws, the stimulus and the noise path. |
| `step_len` | 5 s | A new random step every 5 s. |
| `stim_on` | `[0 0 1 1 1 1 0 0]` | Eight steps of 5 s: off for 10 s, on for 20 s, off for 10 s. |
| `stim_frac` | 0.15 | Share of each cell type chosen anew for every step. |
| `stim_amp` | 0.5 | A chosen neuron's amplitude for a step is `0.5 * randn`, either sign. |
| `lambda_K` | 30 | The bottom row draws the top 30 exponents. |
| `lambda_link` | `[2 3 4]` | Columns 2 to 4 of that row share a y axis; column 1 keeps its own. |

Fixed inside `fig_example_timeseries`, not set by the script: sampling at
400 Hz, the top-K Lyapunov estimator with `lya_dt = 0.05`, a 5 s warm-up,
and accumulation from 5 s to the end. So each exponent is one average over
the rest, the steps and the rest together.

## What a correct run gives

Printed by the script:

| Condition | lambda_1 | Positive among the top 30 |
|---|---|---|
| `no_adaptation` | +2.5514 | 5 |
| `sfa1_std1` (single-timescale) | -0.0991 | 0 |
| `sfa3_stdc2` (chained multiple-timescale) | -0.1015 | 0 |
| `sfa3_std2` (product multiple-timescale) | -0.0933 | 0 |

Stimulus: zero before 10 s and after 30 s. The four steps stimulate 66
(34 E, 32 I), 87 (41 E, 46 I), 74 (37 E, 37 I) and 89 (46 E, 43 I)
neurons; 253 of 500 are stimulated in at least one step; the largest
amplitude is 1.508.

Connectivity: spectral abscissa of W 4.8802, spectral radius 5.2761, bulk
radius R = 4.799.

Bottom row of the time-series figure: columns 2 to 4 on a shared scale of
about -0.6 to +0.3, column 1 on its own of about -7.5 to +3.7. The count of
positive exponents is printed above each panel.

## If the numbers differ

On the same machine, MATLAB release and commit, a rerun gives these
numbers exactly. Elsewhere:

- **The three adapting conditions** should agree to several digits. Their
  networks are stable here, so round-off differences between machines or
  thread counts do not grow.
- **The no-adaptation condition is chaotic.** A last-bit difference in
  arithmetic (another MATLAB release, processor or thread count) grows
  along the trajectory, so its traces and its lambda_1 can differ in the
  second or third digit while the picture stays the same. That is not a
  failure to reproduce.
- **The stimulus counts and the W numbers** do not depend on the
  trajectory. If they differ, the random streams or the preset changed.
- **A later commit** gives the same numbers as long as the model, the
  preset and `fig_example_timeseries` have not changed what they compute.
  If the preset is retuned (the connectivity was under review when this
  was written) the numbers change and this note describes the old network;
  the preset's case block in `src/presets/srnn_param_preset.m` says what it
  holds.

## Related examples

Same network and seed, shorter runs with a single step from 5 to 10 s:
`scripts/examples/mu7chain_step_example` (amplitudes of either sign) and
`scripts/examples/mu7chain_step_signed_example` (positive to E, negative
to I).
