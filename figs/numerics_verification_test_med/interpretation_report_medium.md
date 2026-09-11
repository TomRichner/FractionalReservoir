# Numerical-method verification: interpretation report ('medium' run)

Run: `numerics_verification_test_med_run` at run mode `medium`, 2026-09-10, 17 min.
Preset: `celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25`
(n = 500, 4000 state variables, piecewise activation, sigma_u = 0.025,
three adaptation regimes). Data: `data/numerics_verification_test_med/numerics_verification/numerics_verification_data.mat`.
Figures: `fig_numerics_solver/Fig_numerics_solver.png`, `fig_numerics_lya_method/Fig_numerics_lya_method.png`.

This report supersedes `figs/numerics_verification_test/interpretation_report_fast.md`,
which describes the methods in full (reshooting, the three references, the
shared-Brownian-path construction, and the stimulus-edge trap). Only what
changed at 'medium' is restated here.

What 'medium' changed relative to 'fast':

| setting | fast | medium |
|---|---|---|
| noise-free reshoot window (stimulus-on third of T) | 2.3-3.99 s of 6 s | 4.3-7.99 s of 12 s |
| noisy reshoot window | 1.3-1.99 s of 3 s | 1.8-2.99 s of 4.5 s |
| LLE run / accumulation | 10 s / last 5 s | 20 s / last 10 s |
| reduced network for QR | n = 30, [-8, 4] s | n = 40, [-10, 10] s, 10 s accumulation |
| restarts per segment length (cap) | 400 | 800 |

Two follow-up checks were run by hand after the stage and are reported in
section 3; they are one-line `build_from_preset` calls, not part of the stage.

## 1. Results

### 1.1 Noise-free precision of SRA1 (`Fig_numerics_solver`, rows 2-3)

RMS error per state component over one 0.02 s Benettin segment, SRA1 reset
to the ode45 (1e-10) reference at every restart:

| regime | 400 Hz | 800 Hz | 1600 Hz | slope | relative at 400 Hz (x) |
|---|---|---|---|---|---|
| no adaptation (chaotic) | 3.5e-5 | 8.8e-6 | 2.2e-6 | 2.01 | 1.0e-5 |
| single-timescale | 2.5e-6 | 6.1e-7 | 1.5e-7 | 2.00 | 8.6e-6 |
| multiple-timescale | 8.7e-8 | 2.2e-8 | 5.4e-9 | 2.00 | 7.5e-7 |

Identical to 'fast' within 30%, with twice the restarts and a window that
sits later in the run. Second-order convergence in every regime, in every
state family.

### 1.2 With noise, shared Brownian path (row 3, dashed)

| regime | 400 Hz | 800 Hz | 1600 Hz | slope |
|---|---|---|---|---|
| no adaptation | 4.8e-5 | 1.2e-5 | 3.0e-6 | 2.00 |
| single-timescale | 1.2e-5 | 3.0e-6 | 8.0e-7 | 1.95 |
| multiple-timescale | 4.6e-6 | 1.2e-6 | 3.4e-7 | 1.87 |

Same picture as 'fast': at sigma_u = 0.025 the drift error dominates and the
slope stays near 2. The mild droop toward 1.5 in the adapted regimes is the
stochastic term becoming visible where the drift error is smallest.

### 1.3 LLE, ode45 vs SRA1, full network, 10 s accumulation (row 1)

| regime | Benettin, ode45 1e-10 | Benettin, SRA1 400 Hz | fast run |
|---|---|---|---|
| no adaptation | +4.015 | +3.853 | +3.645 / +3.715 |
| single-timescale | +0.317 | +0.578 | +0.319 / +0.489 |
| multiple-timescale | -0.1193 | -0.1193 | -0.212 / -0.201 |

Chaotic regime: 4% apart. Stable regime: identical to four digits. The
single-timescale regime is again apart, now by 0.26; section 3 resolves it.

### 1.4 Benettin vs QR, reduced network, n = 40 (`Fig_numerics_lya_method`)

| regime | states | Benettin LLE | QR lambda_1 | difference |
|---|---|---|---|---|
| no adaptation | 40 | +0.0576 | +0.0558 | 0.002 |
| single-timescale | 160 | -0.4007 | -0.4010 | 0.0003 |
| multiple-timescale | 320 | -0.1143 | -0.1209 | 0.007 |

All three agree. The two points that were open at 'fast' are closed:

* The 40-neuron no-adaptation network is no longer a dead fixed point (the
  30-neuron one had LLE -9.9 = -1/tau_d). At n = 40 it sits at the edge of
  chaos with an exponent near zero and a local exponent that swings between
  -15 and +15 (row 1, left), which is the demanding case for both methods,
  and they agree to 0.002.
* The single-timescale gap at 'fast' (-1.10 vs -0.50, over 4 s on a
  transient) has vanished with 10 s of accumulation: -0.4007 vs -0.4010.
  The local-exponent panel shows why the short window was unreliable: the
  trajectory is a periodic burst train whose local exponent spikes to +5
  every ~0.7 s and dips to -6 between bursts, so a 4 s average is dominated
  by how many bursts happen to fall inside it.

The QR spectra (row 2) are plausible in every case: a single near-zero
leading exponent, a plateau near -1/tau_d = -10 for the 40-state network,
and staircases at the adaptation and depression rates for the others.

## 2. Interpretation

**The integrator.** Everything the 'fast' run said is confirmed with twice
the data and a later window. At the paper's 400 Hz, SRA1 makes a relative
error of 1e-5 (chaotic), 9e-6 (single-timescale) and 8e-7
(multiple-timescale) per state over one Benettin segment, converging at
second order with or without the paper's noise. In the contracting regimes
this saturates at a few 1e-5 per component and stays there; in the chaotic
regime no integrator can track a trajectory for more than a few Lyapunov
times, and the statistic that survives, the LLE, agrees between SRA1 and a
1e-10 adaptive reference to 4% there and to four digits in the stable
regime.

**Benettin vs QR on `SRNNCellTypePairs` is now verified.** With adequate
accumulation the two independent methods agree to better than 0.01 in all
three regimes, on a network that spans a near-zero exponent, a periodic
burst train and a stable state. The class's Benettin implementation and its
QR implementation (which uses the analytic Jacobian) therefore corroborate
each other, as they already did on `SRNNModel2`.

## 3. The single-timescale LLE gap between integrators is finite-time variability, not bias

The one number that did not settle at 'medium' was the full-network
single-timescale LLE: ode45 +0.317 vs SRA1 +0.578. Two hand checks
discriminate between "SRA1 is biased" and "these are two different
trajectories with different finite-time exponents".

**Refining the step does not move SRA1 toward ode45.** Same seed, same
window, noise off:

| integrator | LLE |
|---|---|
| SRA1 400 Hz | +0.578 |
| SRA1 800 Hz | +0.711 |
| SRA1 1600 Hz | +0.691 |
| ode45 tol 1e-10 | +0.317 |
| ode45 tol 1e-6 | +0.328 |

A step-size bias would shrink monotonically with the step; it does not.
Note also that the three ode45 values (0.319 at 'fast', 0.317, 0.328 at
loose tolerance) are one trajectory: their mutual error is small enough that
they have not separated within 20 s at this exponent. Each SRA1 run, whose
per-segment error is 1e-6 rather than 1e-10, is a distinct trajectory on the
same attractor by the time the window opens.

**Across network seeds the spread is far larger than the gap, and its sign
flips.** SRA1 at 400 Hz vs ode45 at 1e-10, 20 s runs, last 10 s accumulated:

| seeds | SRA1 | ode45 |
|---|---|---|
| [1 2] | +0.578 | +0.317 |
| [2 3] | +0.988 | +0.827 |
| [3 4] | +0.331 | +0.610 |
| [4 5] | +0.543 | +0.430 |

Perturbing only the initial condition (same seed, `x0_std` 0.1001 instead
of 0.1) gives +0.565, so nearby trajectories give nearby values and the
sensitivity is to which trajectory, not to numerical noise.

The finite-time LLE of this regime over 10 s therefore varies by about 0.3
to 1.0 from realisation to realisation, and the ode45/SRA1 difference on any
one seed lies well inside that spread, with either integrator coming out
higher depending on the seed. This is the signature of an intermittent
regime (the local exponent in row 1 of the LLE figure alternates between
bursts of growth and quiet decay). It is a property of the dynamics, and it
is exactly what the paper's reps sweeps and their medians are for; it is not
a defect of the integrator.

## 4. Verdict

* **SRA1 at 400 Hz: verified** (second order, 1e-5 to 1e-7 relative error
  per Benettin segment, LLE integrator-independent wherever the finite-time
  exponent is well defined).
* **Benettin vs QR on `SRNNCellTypePairs`: verified** (agreement to < 0.01
  in all three regimes on the 40-neuron network at 10 s accumulation).
* **A caution for the manuscript, not a concern about the numerics:** in
  the single-timescale regime the 10 s finite-time LLE varies by a factor
  of three across seeds and trajectories. Single-realisation LLE values in
  that regime should not be quoted to better than about +-0.3; medians over
  reps should.

## 5. What to use in the manuscript

The 'medium' figures in this folder are the ones to cite; the 'fast' set is
a plumbing check with shorter windows. If a supplementary statement is
wanted, section 2 gives it in two sentences, and section 3's two tables
justify reporting medians over reps for the single-timescale regime.

## 6. Reproduction

```matlab
setup_paths();
numerics_verification_test_med_run          % ~17 min at 'medium'
```

The follow-up checks in section 3 are `build_from_preset(preset, 'sfa1_std1',
'sigma_u_noise', 0, 'ode_solver', <'sra1'|'ode45'>, 'fs', <fs>, 'T_range', [0 20],
'lya_method', 'benettin', 'lya_T_interval', [10 20], 'rng_seeds', <seeds>)` followed by
`run()`; for the tolerance variant set `ode_opts` before `run()`.
