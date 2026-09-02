# Memory capacity: the protocol as implemented

> **This describes the code, not a plan for it.** The previous version of this
> file was a pre-implementation sketch written before `SRNN_ESN_reservoir`
> existed. It described a white `U(0,1)` input, `f_in = 0.3`, delays counted in
> integrator samples, and a model with a single `b` per neuron — none of which
> is what the code does. Its markdown was also mangled by an underscore-eating
> paste (`\dot{x}*i`, `\sum*{j=1}`), so most of its equations did not render.
> Rewritten 2026-09-02 against the implementation.

The implementation is `SRNN_ESN_reservoir` (`src/model/`), driven by
`run_memory_capacity` (the paired-trial ensemble) and
`run_memory_capacity_example` (one trial, kept per-delay reconstructions), both
in `src/analysis/`.

## What is measured

Drive the reservoir with a scalar input, and for each delay $d$ train a
**separate linear readout** to reconstruct $u(t-d)$ from the reservoir state at
time $t$. Score each on held-out data, and sum:

$$
R^2_d = \frac{\mathrm{cov}\big(u(t-d),\, y_d(t)\big)^2}
              {\mathrm{var}\big(u(t-d)\big)\,\mathrm{var}\big(y_d(t)\big)},
\qquad
\mathrm{MC} = \sum_{d=1}^{d_{\max}} R^2_d .
$$

Note this is the **squared correlation**, not $1 - SS_{res}/SS_{tot}$. The two
coincide for an unbiased least-squares fit on the training set but not on held-out
data, where the squared correlation is insensitive to a constant offset or scale
error in the prediction. `SRNN_ESN_reservoir.compute_R2` implements it, returns 0
when either variance underflows, and clamps to $[0,1]$.

The model itself is stated once, in
[`../EquationsParametersDocs/Equations_stability_paper.md`](../EquationsParametersDocs/Equations_stability_paper.md)
— not restated here. Two features of it matter for the readout and are covered
below: adaptation is normalized by $K$, and depression is **per synaptic route**.

## The input

`input_type` selects one of four; `run_memory_capacity` uses **`'sample_hold'`**
in every run mode.

| type | what it is |
|---|---|
| `sample_hold` | i.i.d. values held for `T_hold` seconds — a staircase |
| `white` | i.i.d. every sample |
| `bandlimited` | white, low-passed at `u_f_cutoff` |
| `one_over_f` | $1/f^{\alpha}$ noise, exponent `u_alpha` |

Values are `u_offset + u_scale * (rand - 0.5)`, i.e. **uniform on $[-0.5, 0.5]$**
at the defaults — centred on zero, not the $U(0,1)$ of the original sketch.

**Sample-and-hold is the load-bearing choice.** A low-pass reservoir cannot
follow white input, and white input's own autocorrelation inflates $R^2$ at
short delays. Holding each value for `T_hold` puts the input's energy where the
network can respond and keeps successive values i.i.d. *at the hold rate*, which
is the canonical setting for a memory-capacity measurement.

**Consequently MC is measured in hold units, not samples.** With
`hold_len = round(T_hold * fs)`, `run_memory_capacity` decimates the state to one
sample per hold (taking the last sample of each), and delay $d$ means $d$ holds
— `delay_scale = hold_len / fs` seconds. At the standard `T_hold = 0.3`,
`fs = 200`, one delay step is 0.3 s and 60 integrator samples. Delay *counts*
from this protocol are therefore not comparable with any measured per-sample.

Input weights: `f_in` of the neurons (default 0.1) are chosen by `randperm` and
given weights `sigma_in * (rand - 0.5)`; the rest get zero. `rng_seed_input`
pins the draw, so `W_in` is shared across conditions.

## The reservoir state that is read out

`readout_signal` chooses what the linear readout sees. Both stack every cell
type in **type order**, so row $i$ is neuron $i$ (`SRNN_ESN_reservoir.stack_by_type`).

- **`'rate'`** — the firing rate $r = \phi(x_{\text{eff}})$. Adaptation is baked
  into $x_{\text{eff}}$, but depression is not visible.
- **`'synaptic'`** (the default) — $r \prod_m b_m \prod_n g_n$, which exposes the
  depression state to the readout. Identical to `'rate'` in conditions without
  STD.

### Why `'synaptic'` needs a check on `SRNNCellTypePairs`

Depression is **per route**: $\theta_j^{q \to s} = r_j \prod b^{q\to s} \prod g^{q\to s}$.
A presynaptic neuron therefore transmits a *different* signal to each target
type, and "the synaptic output of neuron $j$" need not exist at all.

`assert_route_redundancy` requires, **per presynaptic type**, that every one of
that type's outgoing routes carries identical STD *and* STF settings; the readout
then takes route 1. This is exact rather than approximate — since

$$\frac{db_{jm}}{dt} = \frac{1-b_{jm}}{\tau_{rec_m}} - \frac{b_{jm}\,r_j}{\tau_{rel_m}}$$

depends only on the route's constants and the presynaptic neuron's own rate,
from a common $b = 1$, identical routes produce **bit-identical** trajectories.
Reading route 1 is not a shortcut; it is the same array.

The check is *per presynaptic type* rather than global because neuron $j$ of type
$q$ only has routes $q \to 1 \ldots q \to C$; what other types do cannot reach
it. So "E→\* depressed, I→\* not" is well defined and accepted, while STD on
E→E alone is refused with `SRNN_ESN_reservoir:AmbiguousSynapticOutput`. Use
`readout_signal = 'rate'` for a network with genuinely per-route synapses.

`scripts/tests/test_esn_route_redundancy.m` covers both directions and asserts
the bit-identity numerically.

## The readouts

Ridge regression, closed form, one readout per delay:

$$w^{(d)} = \big(X^\top X + \eta I\big)^{-1} X^\top y^{(d)}$$

with `eta = 1e-7`. Trained on the training segment and scored on the test
segment; the first $d$ samples of each segment are lost to the shift.

**The readout must be over-determined to mean anything.** With $n$ features and
$N_{\text{train}} = T_{\text{train}} / T_{\text{hold}}$ hold-samples, $N_{\text{train}} > n$
is required. This is exactly why `'fast'` is documented as producing meaningless
numbers: at $N_{\text{train}} = 200$ against $n = 300$ it is under-determined,
and only the ridge term keeps it finite.

## Run modes

`mc_run_config` inside `run_memory_capacity` is the cost table — the MC analogue
of `analysis_run_config`, kept separate because MC has no grid and its cost is
trials × training duration.

| mode | trials | T_train | N_train | usable? |
|---|---|---|---|---|
| `fast` | 4 | 60 s | 200 | **no** — under-determined; plumbing only |
| `medium` | 15 | 300 s | 1000 | yes |
| `production` | 30 | 600 s | 2000 | the paper's |

Constant in every mode: `input_type = 'sample_hold'`, `T_hold = 0.3`, `fs = 200`,
`readout_signal = 'synaptic'`, `eta = 1e-7`, `T_wash = 10 s`.

**Integrator.** `mc_run_config` names both a deterministic (`det_solver = 'rk4'`)
and a stochastic (`sde_solver = 'sra1'`) scheme and selects between them from the
preset's `sigma_u_noise`, the same rule `analysis_run_config` applies to the
sweeps. Both entry points also take an explicit `'ode_solver'` override: `'sra1'`
is legal at $\sigma = 0$ (`sde_fixed_step` treats an absent noise tensor as zero),
and naming it is how an MC run uses the same integrator as the rest of the
analyses with noise off.

## Conditions, and why the trials are paired

The four regimes come from `srnn_adaptation_conditions` via the preset — the same
`no_adaptation` / `sfa_only` / `std_only` / `sfa_and_std` the sweeps use, under
the same snake_case keys. They are saved as those keys;
`src/figures/helpers/mc_display_names.m` maps them to display text through
`srnn_condition_titles` at plot time.

On `SRNNCellTypePairs` a condition carries `tau_a` (a `1 x C` cell) and a
`synapse_config` struct. It does **not** carry `n_a`, which is Dependent on
`tau_a` and read-only — the reverse of `SRNNModel2`, where the count was set and
the timescales auto-filled.

**Every condition at a given trial shares one network, one `W_in` and one input
sequence**, so a difference in MC is attributable to adaptation and nothing else.
`verify_shared_build` asserts this structurally on the first seed pair rather
than trusting it: it checks that everything except `tau_a` and `synapse_config`
matches across the four conditions, including the protected `W`, `W_in`,
`u_scalar`, `u_ex` and `t_ex`. Trials differ only by seed
(`seed_net_base + k`, `seed_stim_base + k`), which is what makes the `parfor`
order-independent.

## Statistics

Per trial, per condition: `MC`, the memory horizon `H` (the largest delay whose
$R^2$ exceeds `R2_threshold_for_horizon = 0.10`), and the full $R^2_d$ curve.
Across trials: bootstrap 95% CIs on the means (`n_boot`), and **paired
sign-flip permutation tests** between conditions (`n_perm`) — paired because the
conditions share a network, which is the whole point of the design.

## Model class

`SRNN_ESN_reservoir` subclasses **`SRNNCellTypePairs`** (re-parented 2026-09-02,
commit `7966fc7`; it previously subclassed `SRNNModel2`). Memory capacity is no
longer the one part of the paper on a different model class.

The preset is `mc_pairs_dualStd`: `n = 300`, `f = [0.6 0.4]` (off perfect balance
so no-adaptation is not accidentally favoured), `level_of_chaos = 2.0`,
`S_c = 0.35`, a total SFA budget of 0.5 over three timescales, dual-timescale STD
on **all four routes**, and `sigma_u_noise = 0`.

Two of those are deliberate and worth knowing:

- **All four routes**, rather than the paper preset's own arrangement or
  `mc_esn`'s "every outgoing E synapse". Identical STD on all four is what makes
  the `'synaptic'` readout well defined here.
- **Noise off**, for now. The mechanism to turn it on exists and would select
  `sra1` automatically, but the protocol has not been re-validated with noise.

## Reading a result

`run_memory_capacity` returns the path of a `.mat` holding `results_all`, with
`conditions` (snake_case keys), `MC_trials` / `H_trials` (`n_trials x n_cond`),
`R2_trials` (`n_trials x n_cond x d_max_eff`), a `summary` with means and CIs,
and a `settings` struct recording the network, the protocol and the solver
actually used. A `.txt` summary is written alongside it.

`plot_memory_capacity` and `plot_memory_capacity_combined` read that file;
`replot_memory_capacity` regenerates from a saved run without recomputing.

## See also

- [`../EquationsParametersDocs/Equations_stability_paper.md`](../EquationsParametersDocs/Equations_stability_paper.md) — the model
- `src/model/SRNN_ESN_reservoir.m` — the protocol
- `src/analysis/run_memory_capacity.m` — the ensemble, `mc_run_config`
- `scripts/tests/test_esn_route_redundancy.m`, `scripts/tests/test_sample_hold_mc.m`
