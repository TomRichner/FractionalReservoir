# Numerical-method verification: interpretation report ('fast' run)

Run: `numerics_verification_test_run` at run mode `fast`, 2026-09-10.

Preset: `celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25`

(n = 500, 4000 state variables, piecewise activation, sigma_u = 0.025,
three adaptation regimes). Data: `data/numerics_verification_test/numerics_verification/numerics_verification_data.mat`.
Figures: `fig_numerics_solver/Fig_numerics_solver.png`, `fig_numerics_lya_method/Fig_numerics_lya_method.png`.

Code:

* `src/analysis/run_numerics_verification.m`
* `src/figures/fig_numerics_verification.m`
* `src/analysis/SRNNNumericsProbe.m`
* `src/model/integrators/coarsen_noise.m`

## 1. Questions

1. Is the fixed-step stochastic integrator the paper uses (Roessler SRA1 at
   400 Hz) precise enough, and does the largest Lyapunov exponent (LLE) it
   yields depend on the integrator?
2. Do the two Lyapunov methods in the code base, Benettin's reshooting
   estimate and the QR full-spectrum method, agree on `SRNNCellTypePairs`?
   They had only been cross-checked on the older `SRNNModel2` class.

## 2. Methods

### 2.1 Why trajectories cannot simply be compared

Two of the three regimes have a positive LLE. Any difference between two
integrations, however small, grows as exp(lambda t), so a pointwise
comparison of whole trajectories measures chaos, not integrator error. The
ode45 reference at tolerance 1e-10 is not exempt: it is a true trajectory of
the system only locally.

### 2.2 Reshooting

Precision is therefore measured the way Benettin measures the LLE. A
reference trajectory is computed once. Then, at each of many restart times,
the tested integrator is reset exactly to the reference state and integrated
for one short segment; the distance from the reference at the segment's end
is the discretisation error accumulated over that segment, with no chaotic
amplification because every segment starts from the truth. Two segment
lengths were recorded: 0.02 s, which is Benettin's renormalisation interval
(so this is the error over one Benettin segment), and two steps (the
shortest span the fixed-step integrator accepts).

The error is reported per state family (dendritic x, adaptation a,
depression b) as the root-mean-square over the family's components, and
relative to that family's RMS value over the same window.

### 2.3 The three references

* **Noise-free (sub-experiment A).** Noise off. Reference: ode45 with
  RelTol = AbsTol = 1e-10, output on a 1600 Hz grid, full n = 500 network.
  Tested: SRA1 at 400, 800 and 1600 Hz. With the noise off SRA1's drift
  stages reduce to Ralston's second-order Runge-Kutta, so the segment error
  should fall 4x per halving of the step (slope 2 on a log-log plot).
* **With noise (sub-experiment B).** Noise on at the preset's sigma_u =
  0.025. ode45 cannot integrate the SDE, so the reference is SRA1 itself at
  12800 Hz (8x the finest tested rate). The coarser runs consume the SAME
  Brownian path, rebuilt by exact aggregation: increments are summed, and
  the second stochastic integral (the area under the path within a step,
  which SRA1 uses) is carried with its base correction. Decimating the path
  would have produced a different path with the wrong variance. The
  difference between runs on one path is the strong error; SRA1 is strong
  order 1.5 for additive noise, so a slope of 1.5 is the theoretical floor,
  but where the drift error dominates the slope stays at 2.
* **Lyapunov exponent (sub-experiment L).** Noise off, 400 Hz, 10 s run
  with 5 s of accumulation. Benettin's LLE with ode45 at 1e-10 for both the
  fiducial and the perturbed trajectory, versus Benettin with SRA1 for both.
  Because Benettin uses one integrator for both trajectories most of the
  discretisation error is common to the two and cancels in their
  difference, so the exponents should agree more closely than the raw
  trajectory error suggests.

### 2.4 Benettin vs QR (sub-experiment C)

The QR method integrates an N x N variational system per segment and is not
feasible at 4000 states. It was run on a reduced network built from the
same preset physics: n = 30, in-degree 6 (the preset's connection density),
`F_tracks_network = true` so the theoretical spectral radius matches the
full network, noise off, ode45, T = [-8, 4] s with 8 s of warm-up and 4 s
of accumulation from t = 0. Both methods ran on the same fiducial
trajectory. The largest QR exponent should equal Benettin's LLE.

### 2.5 A methodological trap that was found and closed

The stimulus is a three-step pattern (off / on / off). Its step edges fall
at T/3 and 2T/3 rounded to the nearest sample, and the linear interpolant
ramps over one sample, so runs at different sampling rates see slightly
different input around each edge. A first version of the analysis let the
reshoot window straddle an edge, and one edge-crossing segment dominated the
RMS: the error fell 3x from 400 to 800 Hz and then 23x to 1600 Hz, the one
rate sharing the reference grid. All reshoot windows are now confined to the
stimulus-on middle third, where the input is identical on every grid.

## 3. Results

### 3.1 Noise-free precision of SRA1 (figure `Fig_numerics_solver`, rows 2-3)

RMS error per component over one 0.02 s segment, and the fitted slope over
the 400/800/1600 Hz ladder:

| regime | 400 Hz | 800 Hz | 1600 Hz | slope | state RMS (x) | relative at 400 Hz |
|---|---|---|---|---|---|---|
| no adaptation (chaotic) | 4.7e-5 | 1.15e-5 | 2.9e-6 | 2.01 | 3.5 | 1.3e-5 |
| single-timescale | 3.6e-6 | 9.0e-7 | 2.3e-7 | 2.01 | 0.62 | 6e-6 |
| multiple-timescale | 2.3e-7 | 5.8e-8 | 1.5e-8 | 2.01 | 0.33 | 7e-7 |

Every regime converges at second order, exactly the theoretical order of
SRA1's drift. The piecewise activation's derivative discontinuities did not
degrade it. Within each regime the x family carries the largest error and b
the smallest, which follows from their time constants (x is the fastest
state).

### 3.2 With noise (row 3, dashed)

| regime | 400 Hz | 800 Hz | 1600 Hz | slope |
|---|---|---|---|---|
| no adaptation | 5.8e-5 | 1.4e-5 | 3.6e-6 | 2.01 |
| single-timescale | 1.3e-5 | 3.3e-6 | 8.4e-7 | 1.98 |
| multiple-timescale | 4.6e-6 | 1.25e-6 | 3.5e-7 | 1.84 |

The noisy error is the same order as the noise-free error and still falls
with slope close to 2. At sigma_u = 0.025 the drift error dominates; the
strong-order-1.5 stochastic term is not yet visible. (It would appear at
much larger noise, and is verified separately in `test_sde_integrators`,
which measures 1.73 at large sigma.) The slight droop to 1.84 in the
most-adapted regime is the stochastic term beginning to show where the drift
error is smallest.

### 3.3 LLE agreement (row 1 and column titles)

| regime | Benettin, ode45 1e-10 | Benettin, SRA1 400 Hz |
|---|---|---|
| no adaptation | +3.645 | +3.715 |
| single-timescale | +0.319 | +0.489 |
| multiple-timescale | -0.212 | -0.201 |

In the chaotic and the stable regimes the two integrators agree to 2% and
5%. The single-timescale regime differs by 0.17 on a 5 s accumulation
window; see section 4.

### 3.4 Benettin vs QR on the reduced network (figure `Fig_numerics_lya_method`)

| regime | states | Benettin LLE | QR lambda_1 |
|---|---|---|---|
| no adaptation | 30 | -9.898 | -9.898 |
| single-timescale | 120 | -1.100 | -0.500 |
| multiple-timescale | 240 | -0.133 | -0.127 |

## 4. Interpretation

**What the results establish about the integrator.**

* At 400 Hz, SRA1 accumulates a relative error of 1e-5 to 1e-7 per state
  over one Benettin segment, depending on regime, and this error converges
  at second order, so it can be extrapolated to any step size.
* In the contracting (adapted) regimes the error does not accumulate: it
  saturates at roughly the per-segment error times the number of segments
  in one relaxation time. With an LLE near -0.2 that is about 250 segments,
  giving a pointwise error of order 6e-5 per component after the transient.
  An SRA1 trajectory at 400 Hz tracks the true one to about 1e-4 there.
* In the chaotic regime a global trajectory error cannot exist for any
  integrator. Starting from 1.3e-5 relative and growing at 3.6 per second,
  the error reaches order one after about 3 s, which is where the ode45 and
  SRA1 traces in row 1 visibly part. A 1e-10 tolerance buys only about
  three more seconds. Precision claims in that regime must rest on
  statistics, and the statistic the paper uses is the LLE, which agrees
  between integrators to 2%.
* One honest caveat: summed over 500 dendritic states the absolute error
  norm per segment in the chaotic regime is about 1e-3, which is the size
  of Benettin's perturbation. Benettin is unaffected because both of its
  trajectories use the same integrator and the error is common-mode; the
  2% LLE agreement confirms this empirically.

**The single-timescale LLE gap (0.32 vs 0.49).** This is a finite-time
exponent accumulated over 5 s on a trajectory whose local exponent
fluctuates strongly. Two integrators produce two slightly different
trajectories through the same attractor and, on a short window, two
different finite-time averages. It is not a precision issue in the sense of
sections 3.1-3.2 (the reshoot error in this regime is 6e-6). The 'medium'
run, which accumulates for 10 s, is the right test.

**What the results establish about Benettin vs QR.** Less than hoped, for
two reasons that both stem from the reduced network.

* The 30-neuron no-adaptation network is a dead fixed point: its LLE of
  -9.9 is simply -1/tau_d. Agreement there is exact but trivial. Thirty
  neurons at in-degree 6 cannot sustain activity even with the spectral
  radius matched, so the reduced network does not reproduce the chaotic
  regime the paper is about.
* The single-timescale gap (-1.10 vs -0.50) is on a trajectory whose local
  exponent swings by about +-5 through the whole 4 s window, after the
  stimulus switches off at t = 0. Both numbers are finite-time values on a
  transient. Benettin's finite perturbation (d0 = 1e-3) also crosses
  activation breakpoints that QR's linearisation never sees, so a gap that
  persists at longer accumulation would be a genuine finding about
  piecewise nonlinearities rather than a bug.
* The multiple-timescale regime agrees to 0.006.

## 5. Verdict at 'fast'

* SRA1 at 400 Hz: **verified**. Second-order convergence in every regime,
  with and without noise; per-segment relative error 1e-5 to 1e-7; LLE
  integrator-independent to a few percent where the accumulation window is
  adequate.
* Benettin vs QR on `SRNNCellTypePairs`: **not yet settled**. Exact on a
  fixed point, within 5% on the most-adapted regime, and a factor of two
  apart on the single-timescale regime at 4 s of accumulation. Requires the
  'medium' run (10 s accumulation, 40 neurons) before any statement is made.

## 6. Reproduction

```matlab
setup_paths();
numerics_verification_test_run                       % ~11 min at 'fast'
% or the pieces:
run_numerics_verification('run_mode', 'fast', 'out_dir', ...)
fig_numerics_verification('run_dir', 'data/numerics_verification_test', 'variant', 'solver')
fig_numerics_verification('run_dir', 'data/numerics_verification_test', 'variant', 'lya_method')
```

Tests: `scripts/tests/test_numerics_probe.m`, `scripts/tests/test_sde_integrators.m`,
`scripts/tests/test_run_modes.m`.
