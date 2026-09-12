---
title: "Verification of the numerical methods: SRA1 vs ode45, with and without noise, and Benettin vs QR"
subtitle: "Spiking-rate network with SFA and STD, preset `celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25`"
date: 2026-09-11
geometry: margin=2cm
fontsize: 10pt
colorlinks: true
---

# Executive summary

**Run.** `numerics_verification_trials_run`, 2026-09-11, 51 min on 13 workers.
Five network seeds for the integrator-precision experiments, twenty-five for
the two Lyapunov comparisons, in each of the preset's three adaptation
regimes: no adaptation (strongly chaotic, LLE about +3.6), single-timescale
adaptation (weakly chaotic and intermittent, LLE about +0.6), and
multiple-timescale adaptation (stable, LLE about -0.12).

**SRA1 vs ode45 without noise.** At the paper's 400 Hz, SRA1 reset every
0.02 s to a 1e-10 ode45 reference accumulates a per-state error of 3.8e-5
(chaotic), 2.6e-6 (single-timescale) and 7.3e-8 (multiple-timescale), i.e.
relative errors of 1e-5, 5e-6 and 3e-7. On every one of the 15 seed x regime
combinations the error falls with step size at slope 2.00-2.01, the exact
order of SRA1's drift. The piecewise activation costs nothing.

**SRA1 with the paper's noise.** On a shared Brownian path, SRA1 at 400 Hz
differs from SRA1 at 12800 Hz by 4.1e-5, 1.2e-5 and 4.8e-6 per state in the
three regimes, converging at slopes 1.98-2.00, 1.92-1.97 and 1.82-1.88. The
drift error dominates at sigma_u = 0.025; the stochastic term is just
visible where the drift error is smallest. The integrator behaves the same
with the noise on as off.

**Does the LLE depend on the integrator?** No. Paired over 25 seeds, Benettin
with SRA1 minus Benettin with ode45 is -0.03 +- 0.45 (chaotic), +0.05 +- 0.17
(single-timescale) and 0.0000 +- 0.0000 (stable). Standard errors of those
means are 0.09, 0.03 and 0.00; the signs split 14/11, 15/10 and 15/10. In the
two chaotic regimes the per-seed scatter is large, because after a few
Lyapunov times the two integrators follow different trajectories on the same
attractor and a 10 s finite-time exponent depends on the trajectory; the
ensemble means agree to 0.03 and 0.05, and the seed-to-seed spread (sd 0.62
and 0.17) is what the paper's reps average over. In the stable regime the two
integrators give the same exponent to four decimal places on every seed.

**Benettin vs QR.** On a 40-neuron network with the preset's physics, paired
over 25 seeds, QR's largest exponent minus Benettin's LLE is +0.010 +- 0.065
(median |difference| 0.005; the exponents span -10 to +2 and correlate at
0.9999), +0.015 +- 0.058 (median 0.031; span -0.65 to +0.64; correlation
0.985) and +0.0001 +- 0.020 (median 0.005). The two independent estimators
agree across a hundred-fold range of exponents. The one systematic feature
is in the stable regime: QR's lambda_1 is -0.117 +- 0.004 on every seed
(the slowest adaptation recovery), while Benettin scatters +- 0.02 with two
low outliers, because Benettin's perturbation has not always aligned with the
slowest direction within the 10 s warm-up. That is a known finite-time
property of Benettin's method in strongly stable systems, and it biases the
estimate toward more negative values, never toward instability.

**Bottom line for the manuscript.** The 400 Hz SRA1 integrator is verified
in every regime with and without noise; its Lyapunov exponent is
integrator-independent in the ensemble mean; and the Benettin estimate is
corroborated by the full-spectrum QR method. Single-realisation exponents in
the chaotic regimes carry finite-time scatter (sd 0.6 and 0.17 over 10 s),
which is a property of the dynamics, not of the numerics; report medians
over reps, as the sweeps do.

# 1. What was verified, and why

The manuscript's simulations use a fixed-step stochastic integrator
(Roessler SRA1 at 400 Hz) on a 500-neuron rate network with 4000 state
variables, and estimate the largest Lyapunov exponent (LLE) with Benettin's
reshooting method. Two questions had to be answered before those numbers
could be trusted:

1. **Is SRA1 at 400 Hz precise enough**, both without noise (where it can be
   compared with an adaptive reference) and with the noise the paper uses
   (where no adaptive reference exists)? And does the LLE depend on which
   integrator produced it?
2. **Do the two Lyapunov estimators in the code base agree** on the class the
   paper uses (`SRNNCellTypePairs`)? Benettin's method and the QR
   full-spectrum method are independent algorithms; they had only been
   cross-checked on the older `SRNNModel2` class.

Both questions were asked in each of the preset's three adaptation regimes,
which span the dynamical range the paper is about: no adaptation (strongly
chaotic, LLE around +3 to +4), single-timescale adaptation (weakly chaotic
and intermittent, LLE around +0.3 to +0.9) and multiple-timescale adaptation
(stable, LLE around -0.1 to -0.4).

# 2. Methods

## 2.1 Why trajectories cannot simply be compared

Two of the three regimes have a positive LLE. Any difference between two
integrations, however small, grows as exp(lambda t), so a pointwise
comparison of whole trajectories measures the chaos, not the integrator. A
1e-10 adaptive reference is not exempt: it is a true trajectory of the system
only locally, and after a few Lyapunov times it too is just one trajectory
among many on the same attractor.

## 2.2 Reshooting: measuring integrator error without chaotic amplification

Precision is therefore measured the way Benettin measures the LLE. A
reference trajectory is computed once. At each of many restart times the
tested integrator is reset exactly to the reference state and integrated for
one short segment; the distance from the reference at the segment's end is
the discretisation error accumulated over that segment, with no chaotic
amplification because every segment starts from the truth. Two segment
lengths are recorded: 0.02 s, which is Benettin's renormalisation interval
(so this is the error over one Benettin segment, the quantity that matters
for the LLE), and two steps, the shortest span the fixed-step integrator
accepts.

The error is the root-mean-square over state components, reported per state
family (dendritic x, adaptation a, depression b) and relative to that
family's RMS value over the same window. Per condition and trial the
reshoot uses up to 800 restarts for each segment length.

## 2.3 The references

**Noise-free (sub-experiment A).** Noise off. Reference: MATLAB ode45 with
RelTol = AbsTol = 1e-10, output on a 1600 Hz grid, on the full 500-neuron
network. Tested: SRA1 at 400, 800 and 1600 Hz. With the noise off SRA1's
drift stages reduce to Ralston's second-order Runge-Kutta, so the segment
error should fall 4x per halving of the step (slope 2 on a log-log plot).
The preset's piecewise activation has a derivative discontinuity at its
breakpoints, which an adaptive solver steps around and a fixed-step scheme
cannot; a measured slope between 1 and 2 would have indicated
kink-dominated error.

**With noise (sub-experiment B).** Noise on at the preset's sigma_u = 0.025.
ode45 cannot integrate the SDE, so the reference is SRA1 itself at 12800 Hz,
eight times the finest tested rate. The coarser runs consume the SAME
Brownian path, rebuilt by exact aggregation (`coarsen_noise`): increments
are summed, and the second stochastic integral (the area under the path
within a step, which SRA1 uses) is carried with its base correction.
Decimating the path would have produced a different path with the wrong
variance. The difference between runs on one path is the strong error; SRA1
is strong order 1.5 for additive noise, so a slope of 1.5 is the theoretical
floor, but where the drift error dominates the slope stays at 2. Eight-fold
separation from the reference is the minimum the integrator test found
necessary: closer, and the reference's own error correlates with the tested
run's and flatters the slope.

**LLE, integrator dependence (sub-experiment L).** Noise off, 400 Hz, a 20 s
run with Benettin accumulating over the last 10 s. Benettin with ode45 at
1e-10 for both the fiducial and the perturbed trajectory, versus Benettin
with SRA1 for both. Because Benettin uses one integrator for both
trajectories most of the discretisation error is common to the two and
cancels in their difference, so the exponents should agree more closely
than the raw trajectory error suggests.

**Benettin vs QR (sub-experiment C).** The QR method integrates an N x N
variational system per segment and is not feasible at 4000 states. It was
run on a reduced network built from the same preset physics: n = 40,
in-degree 8 (the preset's connection density), `F_tracks_network = true` so
the theoretical spectral radius matches the full network, noise off, ode45,
T = [-10, 10] s with 10 s of warm-up and 10 s of accumulation from t = 0.
Both methods ran on the same fiducial trajectory. The largest QR exponent
should equal Benettin's LLE.

## 2.4 Trials

Every sub-experiment was repeated on five networks per condition, with
`rng_seeds = [k, k+1]` for trial k, which draws the weight matrix, the
stimulus amplitudes, the initial state and the per-neuron setpoints. Within
a trial all sub-experiments share the network. The trial dimension exists
because the earlier single-seed runs (2026-09-10) found that the 10 s
finite-time LLE of the intermittent single-timescale regime scatters by
about +-0.2 from one trajectory to the next with either integrator, so a
single seed cannot separate integrator bias from finite-time scatter. Paired
per-trial values can: bias would show as a consistent sign across trials,
scatter as differences of either sign.

## 2.5 A methodological trap that was found and closed

The stimulus is a three-step pattern (off / on / off). Its step edges fall
at T/3 and 2T/3 rounded to the nearest sample, and the linear interpolant
ramps over one sample, so runs at different sampling rates see slightly
different input around each edge. A first version of the analysis let the
reshoot window straddle an edge, and one edge-crossing segment dominated the
RMS: the error fell 3x from 400 to 800 Hz and then 23x to 1600 Hz, the one
rate sharing the reference grid. All reshoot windows are confined to the
stimulus-on middle third, where the input is identical on every grid.

## 2.6 Implementation

* `src/analysis/run_numerics_verification.m`: the stage; the four
  sub-experiments per trial and condition, with a cost table per run mode.
* `src/analysis/SRNNNumericsProbe.m`: a subclass of `SRNNCellTypePairs`
  that only widens access (inject a Brownian path, keep the path a run
  consumed, integrate a segment from a supplied state). It adds no physics.
* `src/model/integrators/coarsen_noise.m`: exact aggregation of a Brownian
  path (increment and second integral) onto a coarser grid.
* `src/figures/fig_numerics_verification.m`: the three figures.
* `scripts/paper/numerics_verification_trials_{config,run}.m`: this run.
* Tests: `test_numerics_probe` (injection, bit-identical segment
  re-integration from a mid-run state with noise, coarsening identities),
  `test_sde_integrators` (strong orders at large sigma), `test_run_modes`.

## 2.7 Settings of this run

| setting | value |
|---|---|
| preset | `celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25` |
| network | n = 500, in-degree 100, 4000 state variables, piecewise activation, per-neuron setpoint |
| noise | sigma_u = 0.025 (input-referred), additive on x |
| run mode | medium, 5 trials |
| A: noise-free reshoot | ode45 1e-10 at 1600 Hz; SRA1 at 400/800/1600 Hz; T = 12 s, window 4.3-7.99 s |
| B: noisy reshoot | SRA1 12800 Hz reference; SRA1 at 400/800/1600 Hz on the same path; T = 4.5 s, window 1.8-2.99 s |
| L: LLE | Benettin, 400 Hz, T = 20 s, accumulate 10-20 s, lya_dt 0.02 s, d0 = 1e-3 |
| C: Benettin vs QR | n = 40, in-degree 8, F_tracks_network, ode45 1e-9, T = [-10, 10] s, warm-up 10 s |
| restarts per segment length | up to 800 |

# 3. Acceptance criteria

Stated from theory and from the earlier single-seed runs (2026-09-10), before
this ensemble was examined:

| check | criterion | basis |
|---|---|---|
| noise-free convergence | fitted slope of segment error vs step within 1.8-2.2 on every seed and regime | SRA1's drift is second-order RK |
| noisy convergence | slope >= 1.5 on every seed and regime | strong order 1.5 is the floor |
| noise-free error at 400 Hz | relative error per 0.02 s segment < 1e-4 in every regime | the Benettin perturbation is 1e-3 in norm; per-segment error must sit well below it |
| integrator dependence of the LLE | paired mean (SRA1 - ode45) within two standard errors of zero in every regime | no systematic bias |
| Benettin vs QR | paired median |difference| < 0.05 in every regime, and no regime where the sign is consistent across seeds | independent estimators of one quantity |

# 4. Results

Figures: `fig_numerics_solver/Fig_numerics_solver.png` (seed 1 traces, error
per state family, convergence), `fig_numerics_lya_method/Fig_numerics_lya_method.png`
(seed 1 local-exponent overlay and QR spectra), `fig_numerics_ensemble/Fig_numerics_ensemble.png`
(paired per-seed comparisons). Data: `data/numerics_verification_trials/numerics_verification/numerics_verification_data.mat`.

## 4.1 Noise-free precision of SRA1 (sub-experiment A), five seeds per regime

RMS error per state component over one 0.02 s segment, SRA1 reset to the
ode45 (1e-10) reference at every restart; mean over seeds, with the range.

| regime | 400 Hz | 800 Hz (seed 1) | 1600 Hz (seed 1) | slope, all seeds | state RMS (x) | relative at 400 Hz |
|---|---|---|---|---|---|---|
| no adaptation | 3.8e-5 (1.8-5.7e-5) | 8.8e-6 | 2.2e-6 | 2.01, 2.01, 2.01, 2.01, 2.01 | 3.6 | 1.1e-5 |
| single-timescale | 2.6e-6 (2.3-3.3e-6) | 6.1e-7 | 1.5e-7 | 2.00, 2.01, 2.00, 2.01, 2.01 | 0.56 | 4.6e-6 |
| multiple-timescale | 7.3e-8 (4.5-8.7e-8) | 2.2e-8 | 5.4e-9 | 2.00, 2.01, 2.01, 2.00, 2.00 | 0.29 | 2.5e-7 |

Per state family on seed 1 at 400 Hz (x / a / b): 3.5e-5 / - / - ;
4.8e-6 / 9.8e-7 / 3.5e-7 ; 2.2e-7 / 6.0e-8 / 2.5e-8. The dendritic state
carries the largest error and depression the smallest, in proportion to
their time constants. The two-step (local) error is 1.0e-5, 6.5e-7 and
2.2e-8. **Criterion met** on all 15 combinations.

## 4.2 With noise, shared Brownian path (sub-experiment B), five seeds per regime

| regime | 400 Hz, mean (range) | 800 Hz (seed 1) | 1600 Hz (seed 1) | slope, all seeds |
|---|---|---|---|---|
| no adaptation | 4.1e-5 (2.9-5.3e-5) | 1.2e-5 | 3.0e-6 | 2.00, 1.99, 1.98, 1.99, 2.00 |
| single-timescale | 1.2e-5 (1.1-1.3e-5) | 3.0e-6 | 8.0e-7 | 1.95, 1.97, 1.96, 1.92, 1.95 |
| multiple-timescale | 4.8e-6 (4.5-5.0e-6) | 1.2e-6 | 3.4e-7 | 1.87, 1.88, 1.86, 1.84, 1.82 |

The noisy error at 400 Hz is 1.1x, 4.6x and 66x the noise-free error in the
three regimes: where the drift error is large (chaotic) the noise adds
nothing measurable, and where it is tiny (stable) the stochastic term shows
through, pulling the slope from 2 toward the strong-order floor of 1.5
without reaching it. **Criterion met** on all 15 combinations; the noisy
slopes are remarkably consistent across seeds (within 0.06).

## 4.3 Integrator dependence of the LLE (sub-experiment L), 25 seeds per regime

Benettin over the last 10 s of a 20 s run at 400 Hz, noise off. Each seed is
one network; both integrators start from its initial state.

| regime | ode45 1e-10, mean +- sd (range) | SRA1 400 Hz, mean +- sd (range) | paired SRA1 - ode45, mean +- sd | s.e.m. | signs +/- | corr |
|---|---|---|---|---|---|---|
| no adaptation | +3.62 +- 0.62 (2.16 to 4.82) | +3.59 +- 0.72 (1.77 to 4.62) | -0.033 +- 0.448 | 0.090 | 14 / 11 | 0.79 |
| single-timescale | +0.58 +- 0.17 (0.32 to 0.96) | +0.63 +- 0.18 (0.33 to 1.03) | +0.046 +- 0.170 | 0.034 | 15 / 10 | 0.54 |
| multiple-timescale | -0.124 +- 0.016 (-0.146 to -0.070) | -0.124 +- 0.016 (-0.146 to -0.071) | -0.0000 +- 0.0000 | 0.0000 | 15 / 10 | 1.000 |

The paired means are 0.4, 1.4 and 0 standard errors from zero. **Criterion
met** in all three regimes. Two further observations:

* In the stable regime the two integrators agree on every seed to the
  fourth decimal. A decaying perturbation aligns with the slowest direction
  and the finite-time exponent is then a property of the network, not of
  the trajectory.
* In the two chaotic regimes the paired scatter (sd 0.45 and 0.17) is of
  the same size as the seed-to-seed scatter (sd 0.62 and 0.17). After a few
  Lyapunov times the two integrators are on different trajectories of the
  same attractor, so a 10 s finite-time exponent from each is an
  independent draw. The 2026-09-10 controls established that changing the
  initial condition by 1% produces the same scatter with a single
  integrator, and that refining SRA1 to 1600 Hz does not move its value
  toward ode45's. The correlation of 0.79 in the chaotic regime shows that
  the network identity still explains most of the variance there.

## 4.4 Benettin vs QR on the reduced network (sub-experiment C), 25 seeds per regime

n = 40, in-degree 8, spectral radius matched to the full network, noise off,
ode45, 10 s of accumulation after 10 s of warm-up.

| regime | Benettin, mean +- sd (range) | QR lambda_1, mean +- sd (range) | paired QR - Benettin, mean +- sd | median abs | max abs | corr |
|---|---|---|---|---|---|---|
| no adaptation | -5.46 +- 4.42 (-10.0 to +1.99) | -5.45 +- 4.40 (-10.0 to +1.99) | +0.010 +- 0.065 | 0.005 | 0.25 | 0.9999 |
| single-timescale | -0.17 +- 0.31 (-0.65 to +0.60) | -0.15 +- 0.29 (-0.50 to +0.64) | +0.015 +- 0.058 | 0.031 | 0.15 | 0.985 |
| multiple-timescale | -0.117 +- 0.020 (-0.200 to -0.098) | -0.117 +- 0.004 (-0.129 to -0.110) | +0.0001 +- 0.020 | 0.005 | 0.08 | 0.02 |

**Criterion met** in all three regimes. What the ensemble adds to the
single-seed picture:

* The 40-neuron no-adaptation network is heterogeneous across seeds: some
  seeds are dead fixed points (LLE -10 = -1/tau_d), some sit at the edge
  of chaos, two are weakly chaotic (+1.0, +2.0). Benettin and QR track each
  other across that whole range; the largest disagreement (0.25, seed 22)
  is on a seed with LLE -6.5, i.e. 4%.
* In the single-timescale regime several seeds land exactly on -0.500,
  which is -1/(2 s), the recovery rate of the single depression timescale:
  those seeds have settled to a fixed point whose slowest direction is
  depression recovery. The others are bursting or weakly chaotic, and the
  two methods agree to a median of 0.03.
* In the multiple-timescale regime QR's lambda_1 is -0.117 +- 0.004 on
  every seed, which is the slowest adaptation recovery (tau_a = 10 s gives
  -0.100; with the offset and depression it lands at -0.117). Benettin
  gives the same median but scatters +- 0.02 with two low outliers (-0.153,
  -0.200). In a strongly stable system the Benettin perturbation decays
  along whichever direction it started in until the slowest one dominates,
  and 10 s of warm-up is not always enough for that alignment; the
  finite-time estimate is then biased *negative*. The correlation of 0.02
  reflects that QR's values have no spread to correlate with. This is the
  one systematic method effect in the study, and it can only make a stable
  network look more stable, never an unstable one look stable. The
  full-network stable-regime LLE of -0.12 in 4.3 is the same number, so the
  500-neuron and 40-neuron networks agree on what sets stability there.

# 5. Interpretation

## 5.1 The integrator

SRA1 at 400 Hz is second-order accurate on this network in every regime, on
every seed, with the noise on or off, and the piecewise activation's
breakpoints do not degrade the order. The per-segment relative error runs
from 1e-5 in the chaotic regime to 3e-7 in the stable one. What that error
means for a trajectory depends on the regime:

* In the contracting regimes it does not accumulate: it saturates at about
  the per-segment error times the number of segments in one relaxation
  time. At an LLE of -0.12 that is about 400 segments, so an SRA1
  trajectory tracks the true one to a few 1e-5 per state indefinitely.
* In the chaotic regimes no trajectory can be tracked by any integrator for
  more than a few Lyapunov times: a relative error of 1e-5 growing at 3.6
  per second reaches order one in about 3 s, and a 1e-10 reference buys
  about three more. Row 1 of the solver figure shows exactly that. The
  paper's stability claims must therefore rest on statistics, and the
  statistic it uses is the LLE.

## 5.2 The Lyapunov exponent

The LLE is integrator-independent in the ensemble mean, in every regime,
to within its standard error over 25 seeds. What the ensemble makes
unmistakable is the size of the finite-time scatter in the chaotic regimes:
sd 0.6 in the strongly chaotic regime and 0.17 in the intermittent one over
a 10 s window, with either integrator, and with the paired difference
between integrators as large as the seed-to-seed difference. Two
consequences:

* A single 10 s Benettin run in either chaotic regime is one draw from a
  distribution with that width. The number to quote is the median over
  reps, which the sweeps compute; a single-realisation exponent should be
  read with +- 0.6 or +- 0.2 in mind.
* The stable regime has essentially no finite-time scatter (sd 0.016
  across seeds, 0.0000 between integrators). Its exponent is set by the
  slowest adaptation recovery, and the same value emerges on 500 and on 40
  neurons.

## 5.3 Benettin vs QR

The two independent estimators agree on `SRNNCellTypePairs` to a median of
0.005-0.03 across exponents from -10 to +2, and the residual disagreement
in the stable regime is a known finite-time property of Benettin's method
(incomplete alignment within the warm-up) that biases toward more negative
values. This corroborates the class's Benettin implementation and its
analytic-Jacobian QR implementation against each other, as was already
established on `SRNNModel2`.

## 5.4 Interpretation of the ensemble figure

Row 1 of `Fig_numerics_ensemble`: points scatter about the identity line
with no offset in the chaotic regimes and lie on it in the stable one. Row
2: points lie on the identity line across three orders of magnitude in the
no-adaptation column; the stable column's horizontal band is QR's tight
-0.117 against Benettin's scatter, as discussed. Row 3: the reshoot error
varies by a factor of three across seeds in the chaotic regime (the
reference trajectories differ) and by 30% or less in the adapted regimes,
while the slopes do not vary at all.

# 6. Remaining concerns

1. **Benettin in strongly stable regimes is biased negative when the
   warm-up is shorter than the alignment time.** Two of 25 stable-regime
   seeds on the reduced network show it (by 0.04 and 0.08). The full
   network's stable regime showed no such scatter across integrators, but
   its alignment time has not been measured. Where a stable-regime LLE is
   quoted to better than 0.05, either lengthen the warm-up or report QR on
   a reduced network alongside.
2. **The finite-time scatter was measured noise-free.** With the paper's
   noise on, the LLE is a property of the noisy dynamics and the
   trajectory-to-trajectory scatter could be narrower or wider. The reps
   spread at fixed parameters in the existing sensitivity sweeps is that
   measurement; it should be read off before the manuscript states a
   scatter.
3. **The reduced network is not the paper's network.** Forty neurons with
   in-degree eight reproduce the three regimes qualitatively and show that
   the two Lyapunov methods agree on them, but individual reduced networks
   range from dead fixed points to weak chaos under the "no adaptation"
   label. The 500-neuron exponents come from Benettin only.
4. **The stochastic increment handling is verified only where it is
   visible.** At sigma_u = 0.025 the drift error dominates, so the noisy
   slopes near 2 verify the drift on this network and only glimpse the
   strong-order term. The strong order itself (1.73 measured) is verified
   separately at large sigma in `test_sde_integrators`.
5. **Seed means network only when the generator is pinned.** This stage
   pins `twister` in every worker job; the rest of the pipeline does not,
   so the sweeps' "seed k" (built on workers, `threefry`) is a different
   network from a figure's "seed k" (built on the client). Recorded as a
   known bug in `UserNotes.md`; statistically harmless, but it forbids
   claiming that a figure shows one of the sweep's reps.
6. **The analytic Jacobian's finite-difference check is not preserved
   here.** It exists in `test_SRNNCellTypePairs` as a console result;
   priority 2 of the ranked list asks for it to be saved with its
   configuration and tolerance.

# 7. Reproduction

```matlab
setup_paths();
wait_for_parpool(13);                      % holds a PCT licence seat; polls if none is free
numerics_verification_trials_run           % ~50 min on 13 workers, then three figures
```

Earlier single-seed runs and the initial-condition / lya_dt / window
controls:

* `figs/numerics_verification_test/interpretation_report_fast.md`
* `figs/numerics_verification_test_med/interpretation_report_medium.md`
* `scripts/examples/lle_finite_time_variability.m`
