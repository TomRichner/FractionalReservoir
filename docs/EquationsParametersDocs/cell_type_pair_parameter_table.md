# SRNN Pair-Specific Synaptic Dynamics

`SRNNCellTypePairs` models deterministic STD and STF separately for each
presynaptic-to-postsynaptic cell-type pair. Let (q=p(j)) be the type of
presynaptic neuron (j), and let (s=p(i)) be the type of postsynaptic neuron
(i). The recurrent dynamics are

$$
dx_i=\frac{-x_i+u_i+\sum_jw_{ij}
h_j^{q\to s}r_j}{\tau_d}\,dt\;+\;\frac{\sigma_u}{\tau_d}\,dW_i,
$$

$$
h_j^{q\to s}=\left(\prod_{m=1}^{M_{qs}}b_{jm}^{q\to s}\right)
\left(\prod_{k=1}^{L_{qs}}g_{jk}^{q\to s}\right).
$$

For every enabled STD state,

$$
\frac{db_{jm}^{q\to s}}{dt}=
\frac{1-b_{jm}^{q\to s}}{\tau_{rec,qsm}}-
\frac{b_{jm}^{q\to s}r_j}{\tau_{rel,qsm}}.
$$

For every enabled STF state,

$$
\frac{dg_{jk}^{q\to s}}{dt}=
\frac{1-g_{jk}^{q\to s}}{\tau_{dec,qsk}}+
\frac{(G_{qsk}-g_{jk}^{q\to s})r_j}{\tau_{fac,qsk}}.
$$

STD and STF states are independent filters of the same presynaptic rate. Their
products multiply one another in the route readout. Missing mechanisms have an
empty product equal to one.

## Additive Wiener noise on `x`

The (dW_i) are increments of (n) independent standard Wiener processes, one per
neuron, so the model is an Itô SDE whenever `sigma_u_noise` (\sigma_u) is
nonzero, and collapses to a deterministic ODE at (\sigma_u=0).

Noise enters **only** (x). The diffusion coefficient is therefore constant,
which is what makes Itô and Stratonovich agree, kills the Milstein correction
term, leaves the variational equation used by the QR Lyapunov method untouched,
and cancels in Benettin's trajectory difference (both trajectories are
integrated against the same noise path), so the LLE stays measurable at any
noise level.

| Property | Role |
|---|---|
| `sigma_u_noise` (\sigma_u) | Input-referred noise amplitude, in the units of (u), so it is directly comparable to `intrinsic_drive` and the stimulus amplitude. Settable; default `0`. |
| `sigma_x_raw` (\sigma_u/\tau_d) | Coefficient the integrator multiplies (dW_i) by. Dependent, read-only. |
| `x_noise_std` (\sigma_u/\sqrt{2\tau_d}) | Nominal stationary std of (x) from noise alone ((u=0), (W=0)). Dependent, read-only. |
| `noise_seed` | Pins the increment draw; empty derives it from `rng_seeds(1)`. |

`sigma_u_noise > 0` requires a stochastic integrator — `ode_solver` must be
`'euler'`, `'heun'` or `'sra1'` ("sra1" preferred: same two drift evaluations
per step as Heun, strong order 1.5 instead of 1.0). The adaptive solvers cannot
step an SDE, and `'rk4'` is deliberately kept deterministic so (\sigma_u=0) work
stays bit-identical to earlier runs.

The increments are pre-generated for the whole run by `build_noise()` and
indexed by absolute time, so a re-integrated segment sees exactly the increments
the original run used. They are held only for the duration of `run()` and
`compute_lyapunov()`, then cleared — regenerable from `noise_seed`, and large
enough (~96 MB at (n=300), (T=50) s, (f_s=400)) not to be worth saving.

## Named route configuration

Routes use presynaptic name first and postsynaptic name second:

```matlab
synapse_config.E.PV.std = struct( ...
    'tau_rec', [0.2 1], ...
    'tau_rel', 0.25);

synapse_config.E.SST.stf = struct( ...
    'tau_dec', [0.5 2], ...
    'tau_fac', [0.2 0.4], ...
    'G', [1.5 2.5]);
```

| Mechanism | Required fields | Timescale count |
|---|---|---:|
| `std` | `tau_rec`, `tau_rel` | `numel(tau_rec)` |
| `stf` | `tau_dec`, `tau_fac`, `G` | `numel(tau_dec)` |

`tau_rel`, `tau_fac`, and `G` may be scalars, which are broadcast across the
route, or vectors matching its timescale count. Omitted routes have neither STD
nor STF. A route can contain `std`, `stf`, both fields, or neither.

## State layout and readouts

For each configured route, every neuron of the presynaptic type receives the
route's states. Pair arrays use `(pre, post)` ordering:

```matlab
params.state_layout.b{pre, post}
params.state_layout.g{pre, post}
model.plot_data.b.E.PV
model.plot_data.g.E.PV
model.plot_data.synaptic_output.E.PV
```

The total state count is

$$
N_{sys}=N+\sum_qn_qK_q+\sum_{q,s}n_q(M_{qs}+L_{qs}).
$$

States are not pruned when an individual neuron happens to have no connection
to a configured target type. `dead_state_count` reports the realized number of
such states after `build()`.

Connectivity is uniform in `alpha = indegree/n`, while `mu_tilde` and
`sigma_tilde` depend only on presynaptic type.
