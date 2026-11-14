# Plan: release noise (a gamma process) in the model, and its figures

## Context

The model's only noise is a Wiener term on the dendritic state. Tom wants a
second, biophysically motivated source: spontaneous synaptic release, a
non-negative **gamma process** in the synaptic output that reaches the
network through the weights, depletes the vesicle pool, and is seen by the
memory readout. Sections 11 and 12 of
`docs/EquationsParametersDocs/chained_MTS/Equations_chained_MTS.md` worked
out a Poisson and a gamma form on paper; planning showed a simpler form,
the one Tom first wrote ($\theta = r g b_1 + b_1 \cdot$ noise), and he chose
it. **Existing analyses keep the gamma process off** and must not move by a
bit. It is turned on for new memory-capacity maps and new time-series
figures. Target classes: `src/model/SRNNCellTypePairs.m` and
`src/model/SRNN_ESN_reservoir.m`. `SRNNModel2` is defunct and gets no
support (Tom); it is left as it is.

## Decisions

Tom's (2026-11-14):
1. **Noise law: amplitude scales with the pool.** Release of neuron $j$ in a
   step is $F_j\,(r_j\,dt + dR_j)$, with $F_j$ the route's output factor
   ($g\,b_1$ for the chain, the product of pools for the product form, 1 for
   an undepressed route) and $dR_j$ a gamma increment that **does not depend
   on the state**.
2. **The noise axis is $\sigma_\gamma$ at a fixed small leak $r_s$.**
3. **Readout for the new maps: release through a $\tau_d$ filter**, sampled
   at the end of each hold.
4. **A matching Wiener map is made too**: (gain, $\sigma_u$).
5. Top-K on the mean-dynamics Jacobian, unchanged; Benettin as the pathwise
   check. The integrator name stays `sra1`; the release kick is added on top.
6. Nothing an agent runs takes more than about 10 minutes (the quick suite,
   about 13 minutes, is allowed), except **one** essential run at the very
   end. No second essential run.
7. The Fable worker writes all the difficult numerical code itself; Opus
   workers do figures, configurations and sweep plumbing.

The lead's, to be logged under the new issue (FD-027):
- Symbols and names: $r_s$ = `release_rate`, $\sigma_\gamma$ =
  `release_sigma`, both `1 x C` by presynaptic cell type, default 0 (off).
- Values for the new figures: $r_s = 0.01$; $\sigma_\gamma$ from 0 to 0.014,
  where about 0.007 gives the same noise variance in $x$ at gain 1 as the
  paper's $\sigma_u = 0.025$ (the Fable worker measures the constant on the
  network and the axis is set from the measured value). The Wiener map's
  axis is $\sigma_u$ from 0 to 0.05, so equal positions on the two axes mean
  about equal noise in $x$ at gain 1.
- The three new maps use the filtered readout, the Wiener one included, so
  they are comparable with each other; the existing maps keep today's
  readout.
- Release noise is supported on the lean memory path only.

## 1. The model as it will be implemented

For presynaptic neuron $j$, with $dR_j \sim \mathrm{Gamma}(\text{shape } r_s^2 dt / \sigma_\gamma^2,\ \text{scale } \sigma_\gamma^2 / r_s)$,
mean $r_s\,dt$, variance $\sigma_\gamma^2\,dt$, independent across neurons
and steps, one draw per neuron shared by its routes:

- effective release drive: $r_j\,dt + dR_j$ wherever a route's equations use
  the rate to drive release (the output and the depletion of its pools);
- synaptic output: $d\Theta_j = F_j\,(r_j\,dt + dR_j)$, never negative;
- dendrites: $\tau_d\,dx_i = (-x_i + u_i)\,dt + \sum_j w_{ij}\,d\Theta_j + \sigma_u\,dW_i$;
- SFA sees $r_j$ only.

Mean dynamics: the present model with $r \to r + r_s$ at the synapse (a
steady leak). With $r_s = 0$ the mechanism is off and the code path is the
present one.

## 2. Code design

**Split each increment into its mean and a zero-mean remainder.**
- **The mean goes into the drift.** The route supplies the *release drive*
  $r + r_s$ and uses it for its pools' depletion only (not for the
  facilitation or calcium equations, which follow the rate). The class forms
  the synaptic output itself ($\theta = F \cdot r$, in `dynamics_fast`,
  `compute_Jacobian_fast`, `jacobian_times` and
  `unpack_and_compute_states`), so those four functions change by one
  expression each: the route's release drive in place of the rate. No
  branch is added to them, and with $r_s = 0$ the drive is the rate to the
  bit. The Jacobian is then that of the mean dynamics and top-K runs
  unchanged.
- **The remainder $dR_j - r_s\,dt$ is a kick**, applied to the state after
  the whole deterministic step (never inside a stage): a new route-interface
  function `[s_new, dTheta] = kick(s, dRc, P)`. Pools are updated by the
  exact solution of the kick part, so a pool cannot go negative, and
  `dTheta` is the amount released along the kick, not $F \cdot dRc$:
  $-\tau_{rel}\, b_1\, \mathrm{expm1}(-g\, dRc / \tau_{rel})$ for the chain
  and the single pool (which is $\tau_{rel}$ times the pool's loss);
  $-F\, \tau_e\, \mathrm{expm1}(-dRc / \tau_e)$ with
  $1/\tau_e = \sum_m 1/\tau_{rel_m}$ for the product form; $dRc$ for a
  route with no pool. The class adds $W\, d\Theta / \tau_d$ to $x$. (The
  exact update makes the mean spontaneous release smaller than the linear
  law by about $g\, \sigma_\gamma^2 / (2 r_s \tau_{rel})$, under 1 % at the
  paper's $\tau_{rel}$; tested against the closed form and stated.)
- **The realized release is recovered, not stored.** The amount a kick
  released is a function of the state after the step and the increment
  (for the chain $\tau_{rel}\, b_1\, (e^{g\, dRc / \tau_{rel}} - 1)$), so
  the readout filter is computed from the states the lean path already
  holds per segment and the increments. No passive state variable is added
  (Tom's rule: no state that means nothing to the dynamics; a filter state
  per neuron would add $n$ exponents at $-1/\tau_d$).

| File | Change |
|---|---|
| `src/model/routes/route_pool.m`, `route_product.m`, `route_model.m` | `release_rate` in the compiled `P`; the release drive in `rhs`, `linearize`, `steady_state` (depletion terms only); new `kick` and its inverse for the released amount; a pass-through route for a presynaptic type with no synapse dynamics when release is on; `std_zero_floor` refused with release noise (a kick can take the product below its floor) |
| `src/model/SRNNCellTypePairs.m` | the four functions use the route's release drive where they form $\theta$; properties `release_rate`, `release_sigma` (scalar aliases for sweep axes), validation, `get_params` / `build_noise` / `integrate` plumbing; both added to `route_fields` |
| `src/model/integrators/sde_fixed_step.m`, `check_noise_settings`, `resolve_solver` | an opt-in hook for a state-dependent kick after the existing final Wiener update in `euler`, `heun`, `sra1`; **noise is "on" when either source is on** (today every switch tests `sigma_u_noise` alone, so a release-only run would silently get the leak and no kicks); release noise on a deterministic solver is an error |
| `src/model/integrators/` (new file beside `wiener_increments.m`) | the gamma increments: a deterministic function of (seed, absolute block of steps) through `RandStream` substreams, so a segmented run and a single-pass run get identical numbers; honors `noise_seed_schedule`; $\sigma_\gamma = 0$ is the exact leak with no draw |
| `src/presets/analysis_run_config.m`, the memory stage's solver choice | "stochastic" also when release noise is on |
| `src/model/SRNN_ESN_reservoir.m` | lean path: a windowed increment tensor per segment; the filtered-release readout (a new `readout_signal` name) computed from the segment's states and increments, carried across segments; the non-lean path refuses release noise |
| `src/model/lyapunov/lyapunov_topk.m`, `SRNNModel2.m`, `ParamSpaceAnalysis2.m` | **not changed** |

Reused as they are: `wiener_increments` and the absolute-time indexing of
`sde_fixed_step`; `noise_seed_schedule`; the route registry and the
per-neuron draw convention (`route_model.m`, "THE SPREAD");
`assert_route_redundancy` (draws are per presynaptic neuron, so a neuron's
routes stay identical); `mc_protocol`; the lean path's segment loop.

## 3. Analyses and figures

All in the A + C bundle (`scripts/paper/mu7chain_AC_config.m`), in every
run mode, block C's sizes (`mc_gain_sweep_grid`: essential 3 x 3 x 1
network, fast 7 x 7 x 3, medium 13 x 13 x 10).

**Three new memory-capacity grids** (maps with linear color scales,
robustness, statistics, as the existing grids):
1. gain x $\sigma_\gamma$ (Wiener off, $r_s = 0.01$);
2. self-coupling x $\sigma_\gamma$ at gain 1;
3. gain x $\sigma_u$ (release off).

Grid 2 has no gain axis, and the sweep's first axis is hard-wired to
`level_of_chaos` (job arguments, job-file names, table column, readers,
figures, statistics). It is generalized the way the second axis was
(`mc_sweep_second_axis`): the first axis is any scalar model property, gain
by default with today's names, so existing results and job files still read
and resume.

**Two time-series figures** (new figure function, columns the four
conditions, three bands of rows for three noise levels, each band showing
$x$, $r$, the synaptic output and the depression factor for a sample of
neurons; same network and seeds in every panel):
- Wiener: $\sigma_u$ = 0, 0.025, 0.05;
- gamma: $\sigma_\gamma$ = 0 (leak only), 0.007, 0.014 at $r_s = 0.01$,
  Wiener off.

Essential-run cost: block C grows from 24 to 51 jobs, so the stage is about
50 minutes and the whole essential run about 1 h 30. A fast run grows by
about 6 hours (Tom's to start).

## 4. Tests and references

- **Off means identical:** both ultrarapid tiers, the four tolerance suites,
  the stages tier and the quick suite pass with no reference touched.
- **Unit tests** (new, each seconds): moments of the increments (tolerance
  set from their large kurtosis) and of the realized release against the
  closed form; segmented = single pass; a schedule is bit-identical before
  its switch; release and pools never negative; kick(kick(s, a), b) =
  kick(s, a + b); the recovered release equals the kick's; for one
  releasing pool, release equals $\tau_{rel}$ times the pool's loss; a
  leak-only run equals the mean dynamics; **a release-only run has kicks**;
  release noise on `rk4` errors; finite-difference Jacobian and
  `jacobian_times` with the leak; top-K against Benettin at the paper's
  level and at four times it; strong order 1 in the step (coarse increments
  made by summing fine ones); two routes of one type stay `isequal`; the
  lean readout mode; the first-axis generalization reads old job files.
- **A new exact ultrarapid tier** for release-noise cases (four conditions
  on the regression network, Benettin, top-K, lean memory capacity), with
  its own reference, captured once, deliberately, when the numerics are
  final; and a **tolerance suite** for the same cases.
- **A mutation list of 14** (new arithmetic in the dynamics and the
  Jacobian): the leak's sign in the pool's release; the leak missing from
  $\theta$ in `dynamics_fast`; missing from one Jacobian function only; the
  leak wrongly in the facilitation equation; on the product's facilitation
  columns; $r_s + 10^{-6}$ in `linearize` only; `steady_state` without the
  leak; `dTheta = F * dRc`; the kick's sign; one $\tau_{rel}$ for every
  product pool; the kick into $x$ without $1/\tau_d$; an uncentered
  increment; the increment of step $j + 1$; the shape without $dt$.
- No run over about 10 minutes; the quick suite once per merge to main.

## 5. Order of work

| Step | What | Worker |
|---|---|---|
| 0 | Issue FD-027, this plan saved as `docs/Plan_release_noise.md`, decisions logged | lead |
| N1 | Route models, integrator hook, gamma increments, model properties; unit tests; off = identical | `matlab-fable-1` |
| N2 | Reservoir lean path and the filtered readout; tests | `matlab-fable-1` |
| N3 | Lyapunov checks (top-K against Benettin), step-size check, the new exact tier, tolerance suite and mutation list; equations: `Equations_stability_paper.md` and a new self-contained section 13 in the chained-MTS file ("as implemented"), with a note in sections 11 and 12 | `matlab-fable-1` |
| P1 (beside N1) | First-axis generalization of the sweep, old results still read | `matlab-opus-1` |
| P2 (beside N1) | The time-series figure function, Wiener version first | `matlab-opus-3` |
| P3 (after N2) | The three grids in the bundle, their maps and statistics; the gamma time-series figure; registry and bundle tests | `matlab-opus-1`, `matlab-opus-3` |
| G | Each branch reviewed and merged by the lead, one at a time; quick suite and stages tier on main; push | lead + a gate worker |
| E | **One** essential run of the A + C bundle with everything; the lead reviews it and reports | `matlab-opus-2` |
| C | `CLAUDE.md` (release noise; `SRNNModel2` defunct), issue, worklog, figure plan, timeline | lead |

If a stage of the essential run fails, the fix is tested alone and only that
stage and its figures are rerun into the same run directory (as was done for
the fast run); the essential run is not repeated.

## 6. Verification

- Unit tests and the new tier as in section 4; the quick suite (about 13
  minutes) and the stages tier on merged main.
- Numerics, reported with numbers by the Fable worker: mean and variance of
  the realized release against $r_s$ and $\sigma_\gamma$; $\lambda_1$ from
  top-K against Benettin at three noise levels on the regression network;
  the measured Wiener-equivalent $\sigma_\gamma$; cost of the release noise
  per memory run.
- The essential run: every stage and figure present in the manifest; the
  lead looks at the three new maps and both time-series figures; the
  $\sigma_\gamma = 0$ and $\sigma_u = 0$ columns compared with the existing
  self-coupling and bias maps at matching points.

## 7. Risks

- **Top-K's bias under release noise.** The mean-dynamics tangent leaves
  out the kick's derivative; its first-order effect is zero and the
  second-order shift of $\lambda_1$ grows as $\sigma_\gamma^2$ times the
  depression share of the leading vector (a rough estimate: a few
  hundredths per second at the Wiener-equivalent level on the paper's
  network, several times that at the top of the axis). The Fable worker
  measures it against Benettin at three levels in N3. If it is not small
  against the features of the maps, the maps' $\lambda_1$ under release
  noise is computed by Benettin (the sweep's Lyapunov run is 30 s), and that
  decision is logged.
- **SRA1's order** holds for the Wiener part only; the kick is first order.
  Checked by the step-size test; if 200 Hz is not enough the new grids run
  at 400 Hz and the cost is restated.
- **A defect found while planning, to be reported and not absorbed:**
  `SRNN_ESN_reservoir.run_reservoir_esn` clears the noise increments before
  `run_memory_capacity` computes Lyapunov exponents, so Benettin on a noisy
  reservoir re-integrates without the noise. The sweep does not use that
  path (its $\lambda_1$ comes from a separate `run()`), so it is noted under
  the issue for Tom.
- **Lumpy increments** (shape per step about 0.01): rare large kicks; the
  exact pool update keeps pools positive; the tests include the tail.
- **Memory** of Wiener increments on the lean path (about 1.2 GB per job at
  n = 500) for the Wiener map: measured in N2; if 12 jobs do not fit, the
  Wiener increments are generated per segment for that path too.
- **The leak changes the mean network** slightly even at $\sigma_\gamma = 0$;
  the report states the size of that shift.
