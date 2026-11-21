# Solver and rate audit, 2026-11-17

Decision (Tom, 2026-11-17): every simulation in every paper analysis uses SRA1
at fs = 500 Hz (fixed step 2 ms). Without noise SRA1 is Ralston's second-order
Runge-Kutta, so the Methods state one scheme. The one exception is
`run_numerics_verification`, which compares solvers and rates on purpose.

Branch `solver500` (worker matlab-opus-4). Read from the code paths and from
the saved `settings` of the runs the manuscript's figures come from
(`data/mu7chain_AC_fastmedium`, 2026-11-15; `data/mu7chain_gainbias_fastmedium`,
2026-11-17; `data/IED_model/stim_blocks_20261117`).

## What the current figures used, and what the bundle now hands every stage

| Paper figure | Stage (or figure that simulates) | Solver, fs used by the current figure | Noise | After this branch (bundle run) |
|---|---|---|---|---|
| Fig 1 | `fig_introductory_concepts` (simulates while drawing, `sompolinsky_pairs`) | rk4, 200 Hz (build_from_preset at sigma 0; figure default fs) | off | sra1, 500 (registry passes `fs`, `ode_solver`) |
| Fig 1 | `fig_FI_curve`, the other concept panels | no simulation | - | - |
| Fig 2 | `paper_illustrations` / representative_dynamics | sra1, 400 | on (0.025) | sra1, 500 |
| Fig 3 | `local_lyapunov` (steps5s preset) | sra1, 400 | on | sra1, 500 |
| Fig 3 | `eig_heatmap` (reference example) | sra1, 200 (mode table) | on | sra1, 500 |
| Fig 4 | `sweeps`: seven 1-D sensitivity sweeps, joint param_space | sra1, 200 (mode table) | on | sra1, 500 |
| Fig 5 | `eig_heatmap` (transient amplification) | sra1, 200 | on | sra1, 500 |
| Fig 5 | `transient_gain` | sra1, 400 | on | sra1, 500 |
| Fig 6 | `sweeps`: tau sweep | sra1, 400 (mode table) | on | sra1, 500 |
| Fig 6 | `memory_capacity` (mc preset, lean) | sra1 (cfg.mc_ode_solver), 200 (mc_protocol) | off | sra1, 500 |
| Fig 7 | block C `mc_gain_sweep`, grid `selfcoupling` (gainbias run) | sra1, 200 (mc_protocol) | off | sra1, 500 (protocol.fs) |
| Fig 7 | block C `mc_gain_sweep`, grid `noise1d_wiener` | sra1, 200 | Wiener on (swept) | sra1, 500 |
| (Fig 7 data, not drawn) | block C `noise1d_release` | sra1, 200 | release noise on | sra1, 500 |
| Fig 8 | `ied_stim_blocks` (run `stim_blocks_20261117`) | sra1, 400 | on (0.025) | sra1, 500 already in `ied_stim_config` 'production' (no change needed) |
| Supp single neuron | `paper_illustrations` / single_neuron (`single_neuron_mu7chain`) | rk4, 400 | off | sra1, 500 |
| Supp SFA/STD steady state | analytic / built model only, no run | - | - | - |
| Supp eig heatmap imbalance | `eig_heatmap` (two imbalanced examples) | sra1, 200 | on | sra1, 500 |
| Supp local vs finite LLE, EI weights, sensitivity panels | `sweeps` | sra1, 200 | on | sra1, 500 |
| Supp numerics (three figures) | `numerics` | ode45 and sra1 over an fs ladder, by design | off/on | unchanged (the exception) |
| Supp MC example | `mc_example` | sra1, 200 | off | sra1, 500 |
| Supp bursting (two figures) | `fig_stim_engages_adaptation` (simulates while drawing, `bursting_pairs`) | rk4, 400 | off (the input carries a stimulus noise probe) | sra1, 500 |
| Supp IED state | `ied_stim_blocks` | sra1, 400 | on | sra1, 500 |
| not in paper | `lyapunov_spectrum` noise_on / noise_off | sra1 / rk4 (hard-coded), 400 | on / off | sra1 / sra1, 500 |
| not in paper | `dc_lle` | sra1, 400 | on | sra1, 500 |
| not in paper | `noise_switch` rollouts | sra1, 200 (store 200, example 100) | on | sra1, 500 (store 250, example 125) |
| not in paper | `fig_noise_levels_timeseries` (block C, simulates while drawing) | analysis_run_config's sde solver and fs at the run mode | on | unchanged (not plumbed) |

`SRNNCellTypePairs`' own default (`ode_solver = 'ode45'`) is never reached by a
paper stage: every stage builds through `build_from_preset` (sra1 with noise,
rk4 without), `analysis_run_config` (mode table, sra1 when noisy), or names
its integrator (memory capacity, block C, dc_lle, noise switch, IED).

## How it is enforced

- `src/util/stage_solver.m`, the twin of `stage_fs`: '' keeps the stage's own
  choice, a name is validated (`check_ode_solver`).
- `mu7chain_config`: `cfg.ode_solver = 'sra1'` beside `cfg.fs = 500`;
  `cfg.mc_ode_solver` is set from it; `cfg.mc_lean = true` (at 500 Hz the full
  memory-capacity path would hold about 2.5x its 200-Hz state series).
  `mu7chain_AC_config` and `mu7chain_gainbias_config` inherit both.
  `mu7chain_AC_config`: `cfg.noise_switch.store_fs = 250`, `example_fs = 125`.
- `run_all_paper_analyses` hands `'fs'` and `'ode_solver'` to the sweeps
  (`run_all_analyses` -> `resolve_run_context`), eig_heatmap,
  lyapunov_spectrum, transient_gain, local_lyapunov, dc_lle, memory_capacity
  and mc_example (new `'fs'` option), and paper_illustrations reads
  `cfg.ode_solver`. The opt-in stages read the bundle cfg in their plans:
  `mc_gain_sweep_stage_plan` (solver: `mc_ode_solver`, else `ode_solver`; rate:
  `protocol.fs`, which is in every job's signature, so a 500-Hz grid never
  resumes into 200-Hz job files) and `noise_switch_stage_plan` (`fs`,
  `ode_solver`, `store_fs`, `example_fs`; `run_noise_switch_rollouts` has new
  `fs` and `ode_solver` options). `run_manifest` records `ode_solver`.
- The two figures that simulate while drawing take `ode_solver`; the
  `mu7chain_config` registry passes `fs` and `ode_solver` to both.
- Stand-alone calls are unchanged: every new option defaults to the stage's
  previous choice.
- Test: `scripts/tests/test_bundle_fs_plumbing.m` (dry runs of
  `run_all_paper_analyses` for three bundles at 'fast' and 'production', every
  stage stubbed and checked for 500 and 'sra1'; then the real eig_heatmap,
  local_lyapunov, paper_illustrations and the two figures against a
  `build_from_preset` stub, checked to hand both to the model).

## Pilot at 500 Hz (2026-11-17, 6-worker pool, machine otherwise idle)

| Pilot | 200 Hz | 500 Hz | Ratio | Peak worker working set |
|---|---|---|---|---|
| Block A memory capacity, one lean trial, one condition (sfa3_stdc2), 760 s | about 2.0 min (from the fastmedium log: 15 trials x 4 conditions in 15.7 min, two rounds on 12 workers) | 3.96 min (3 trials in parallel) | about 2.0 | 0.97 GB |
| Block C job, gain 1 x self-coupling 1, one seed, one condition | 113 s | 273 s (run side by side) | 2.42 | 1.48 GB |
| local_lyapunov, one condition, K = 50, T = 40 s, top-K + Benettin | 1.9 min (400 Hz, the fastmedium log: 7.6 min for 4 conditions) | 3.13 min | 1.65 | client 5.9 GB peak |

## Estimates for the overnight runs (13-worker pool, three conditions)

Base: `data/mu7chain_AC_fastmedium/command_window.log` (2026-11-15, four
conditions, 12 workers, 200 Hz in most stages, 400 Hz in some) and
`data/mu7chain_gainbias_fastmedium` (2026-11-17, three conditions, 12 workers,
200 Hz). Scaling: x2.4 from 200 Hz, x1.25 to x1.65 from 400 Hz, x0.75 for three
conditions instead of four, x12/13 for the parallel stages.

| Stage | logged | estimate at 500 Hz |
|---|---|---|
| sweeps | 83 min | about 2.2 h |
| memory_capacity | 16 min | about 25 min |
| mc_example | 4 min | about 7 min |
| eig_heatmap | 17 min | about 30 min |
| numerics | 21 min | about 15 min (own rates) |
| lyapunov_spectrum | 8 min | about 10 min |
| transient_gain | 77 min | 1.2 to 1.6 h |
| local_lyapunov | 8 min | about 10 min |
| dc_lle | 28 min | about 30 min |
| paper_illustrations | 3 min | about 4 min |
| **block A total** | | **about 5.5 to 6 h** |
| mc_gain_sweep (selfcoupling + two 1-D noise grids, 715 jobs) | 443 min (gainbias run) | **about 16.5 h** |
| noise_switch | 1.4 min | about 3 min |
| figures | | under 1 h |
| **mu7chain_AC_run('fastmedium')** | | **about 22 to 23 h** |
| `ied_stim_run('production')` (96 blocks of 310 s, 108 worker-min at 400 Hz) | | about 15 to 25 min plus figures |

Memory: the measured peaks are at most 1.5 GB a worker and 6 GB in the client;
13 workers at 1.5 GB is about 20 GB of 63 GB. The stages that keep the full
state (eig_heatmap, transient_gain, lyapunov_spectrum) were not measured at
500 Hz; their trajectories are 1.25x (from 400) to 2.5x (from 200) their
previous size, and the 12-worker runs at those rates fitted.
