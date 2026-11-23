# S2 panel B: why the QR local rate and Benettin's local rate differ

2026-11-23, branch `p3-s2b`. Supplementary Figure S2, panel B
(`fig_numerics_verification`, variant `lya_method`, data of sub-experiment C of
`run_numerics_verification`): in the multiple-timescale column the "QR local
lambda_1" swings by about +-3 s^-1 while Benettin's local rate stays near zero.

## Answer

Both traces are correct, but they are the growth rates of **two different
vectors**, and in that run neither vector has reached the leading direction.
The local rate is a property of the vector, and the two methods start from
different ones: Benettin from a random vector, the QR method from `eye(N)`,
whose first column is the first state variable (an adaptation variable of
neuron 1). The first QR column is exactly a single vector evolved alone, so
"QR local lambda_1" is Benettin's method started from e_1. The two agree only
after both vectors have aligned with the leading direction, which takes of
order 1 / (lambda_1 - lambda_2). In the multiple-timescale network the leading
exponents are the slow adaptive rates and are nearly equal (gap 0.005 to 0.07
s^-1, against 8 s of warm-up), so they have not aligned. It is not the segment
length, not the solver step and not precision.

A second fact about that particular network (30 neurons, seeds [1 2], no
noise): its activity is a decaying transient. It settles on a fixed point
(|dS/dt| falls to 1e-8 by 200 s), so the swings of panel B are the rates along
a transient, there is no stationary local rate to estimate, and the true
exponents are the eigenvalues of the Jacobian at the fixed point.

## What each trace is (from the code)

- Run: reduced network, `n_small` = 30 neurons (240 states with chained
  multiple-timescale adaptation, 120 with single), `sigma_u_noise = 0`, ode45
  (1e-9), `T_range = [-8, 4]`, both methods iterate from -8 s and accumulate
  over [0, 4]. No smoothing in the stage or in the figure.
- Benettin: one random vector, finite perturbation d0 = 1e-3, re-integrated
  with ode45 and renormalised every `lya_dt` = 0.02 s (600 points).
- QR: N vectors (all 240, not 20; 20 is the top-K method), ode45 on the
  variational equation, QR every 0.1 s (120 points). The columns are sorted by
  their time-averaged exponent before saving; here sorted column 1 is raw
  column 1 in both conditions.

## Evidence (30 neurons, seeds [1 2] unless said)

`s2b_01_window` (the stage's window):

- Re-derivation (RK4 on the variational equation) reproduces the saved QR
  trace to 0.001 s^-1.
- Six random vectors, each evolved alone, multiple-timescale: the swing (std
  over t >= -4) of their local rates is 0.15, 1.07, 0.97, 1.05, 1.56, 1.09
  s^-1; the QR first vector 1.15; the saved Benettin trace 0.27. Pairwise
  |cos| between the six vectors at t = 4: 0.11 to 0.80 (median 0.45). Single-
  timescale: 1.46 to 1.90, QR 1.91, |cos| 0.99 to 1.00.
- The class's Benettin (d0 = 1e-3, ode45, 20 ms) started from e_1 reproduces
  the QR first-column trace: rms difference 0.010 s^-1 (std 1.18 against 1.18).
  Finite against infinitesimal perturbation from the same random vector: 0.005.
- Segment 5, 10, 20, 40, 100 ms: swing of the random vector 0.157, 0.157,
  0.157, 0.157, 0.154; of the QR first vector 1.175, 1.175, 1.175, 1.173,
  1.154. (100 ms only clips the sharpest peaks a little.)
- Half the step: the 0.1 s rates move by at most 0.004 s^-1. ode45 at 1e-10 on
  the variational equation of e_1: within 0.003 s^-1 of the QR column.
- The random vectors are not inside the span of the leading QR vectors
  (projection on the first 20: 0.74 to 0.97), so the mean of the first k QR
  rates does not reproduce Benettin's trace (rms 0.64 for k = 2, 0.47 for
  k = 3 on the saved data). The near-degenerate group is wide: 15 exponents
  within 0.05 s^-1 of lambda_1.

`s2b_02_long` (the same networks to 200 s):

- Multiple-timescale: fixed point; exponents = Re eig J = -0.067, -0.074,
  -0.079, -0.081, -0.093 ... (QR over [100, 200]: -0.070, -0.075, -0.077,
  -0.084). Gap lambda_1 - lambda_2 = 0.005 to 0.007 s^-1. At 200 s the six
  random vectors still have |cos| 0.36 to 0.98 with the QR first vector.
- Single-timescale: sustained oscillation, lambda_1 = -0.002, lambda_2 =
  -0.135; all vectors within |cos| 0.99 by t = 4.

`s2b_03_seeds` (six seed pairs at 30 neurons, two at 60, both conditions):

- The disagreement depends on the network, not only on the condition. It is
  large (rms of random minus QR rate over [-4, 4] above 0.4 s^-1) in the
  multiple-timescale condition for 30 neurons seeds [1 2] 0.74, [4 5] 0.50,
  [5 6] 0.95 and for both 60-neuron networks (1.27, 1.29), and small (below
  0.15) for seeds [2 3], [3 4], [6 7], whose activity has already died in the
  window (both traces flat). The vectors are unaligned in 7 of the 8
  multiple-timescale networks (median pairwise |cos| 0.12 to 0.45).
- Single-timescale is not immune: seeds [3 4] 0.70, 60 neurons [2 3] 0.89,
  [1 2] 0.48 (early part of the window).

`s2b_04_warmup` (sustained activity, run to 120 s):

- Multiple-timescale, 60 neurons (lambda_1 near 0, gap 0.063 and 0.074): the
  disagreement decays as exp(-gap t): median rms 1.0 at 0 to 8 s, 0.1 to 0.25
  at 24 to 32 s, below 0.05 from 40 s, 0.001 at 90 s. A warm-up of about 50 s
  is enough there.
- Multiple-timescale, 30 neurons, seeds [4 5] (lambda_1 = +0.021, lambda_2 =
  -0.005, gap 0.026): still 0.5 to 1 at 100 s. Here the two leading exponents
  are close enough that no affordable warm-up separates them.
- Single-timescale, 60 neurons, seeds [2 3] (two exponents near zero, gap
  0.002): levels off at 0.13.

## Tom's hypotheses

- Smaller segment: ruled out (no change from 5 to 100 ms).
- Precision or integration: ruled out (finite against linear 0.005 to 0.01;
  half step 0.004; ode45 1e-10 reference 0.003 s^-1).
- Leading exponents nearly equal: confirmed, with one correction. Benettin's
  vector is not switching between two directions; neither method's vector has
  converged, each is a different mixture of a wide group of slow directions,
  and the local rate of a mixture depends on the mixture.
- Depends on the network: yes (see seeds).

## Not determined

- The stage's own Benettin vector cannot be reproduced exactly (its random
  draw depends on the generator state after the run), so the saved Benettin
  trace is matched in kind (one of six random vectors has the same small
  swing, 0.15), not sample for sample.
- With the input noise of the paper (this check runs without noise) the
  network would not settle on a fixed point; whether the leading exponent then
  separates from the slow group at 30 neurons was not tested.
- Sizes above 60 neurons were not run.

## Files

- `s2b_build.m` the stage's reduced network, full state kept.
- `s2b_tangent.m` QR vectors and independent single vectors on the
  variational equation (RK4 on the stored grid), unsorted columns.
- `s2b_benettin_fd.m` the class's finite-perturbation Benettin with a chosen
  initial vector.
- `s2b_01_window.m`, `s2b_02_long.m`, `s2b_03_seeds.m`, `s2b_04_warmup.m` the
  experiments (about 1, 2.5, 1 and 1 minute, no pool); results in
  `data/explorations/s2b_estimators/`.
- `s2b_05_figs.m` draws `figs/explorations/s2b_estimators/S2B_1_window_vectors.png`,
  `S2B_2_long_run.png`, `S2B_3_seeds_warmup.png`.

Run from this folder after `setup_paths` (the folder is not on the path; `cd`
into it). The stage's saved data is read from
`data/mu7chain_AC_fast/numerics_verification/numerics_verification_data.mat`
(a copy of the `fast` run's file).
