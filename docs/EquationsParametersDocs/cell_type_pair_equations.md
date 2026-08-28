# SRNNCellTypePairs — equations and implementation reference

The **per-route** form of the model, as `SRNNCellTypePairs` implements it, plus
the pieces shared with every class: the activation function and the connectivity
scaling factor.

For the model stated generally — one E population, one I population, no route
indices — see [`Equations_stability_paper.md`](Equations_stability_paper.md).
That file is the paper's statement; this one is what the code actually does when
depression and facilitation are configured per presynaptic-to-postsynaptic pair.

---

## System equations

Let $q = p(j)$ be the type of presynaptic neuron $j$ and $s = p(i)$ the type of
postsynaptic neuron $i$. The recurrent dynamics are

$$
dx_i=\frac{-x_i+u_i+\sum_jw_{ij}\,
\theta_j^{q\to s}}{\tau_d}\,dt\;+\;\frac{\sigma_u}{\tau_d}\,dW_i,
$$

where the rate is the **pre-depression** output of the nonlinearity,

$$
r_i = \phi\!\left( x_i - a_{0_i} - \frac{c_s}{K_s} \sum_{k=1}^{K_s} a_{ik} \right),
$$

and each spike-frequency-adaptation state relaxes towards that rate,

$$
\frac{da_{ik}}{dt} = \frac{-a_{ik} + r_i}{\tau_{a,sk}}, \qquad k = 1, \dots, K_s .
$$

$a_{0_i}$ is the **setpoint** of the nonlinearity — the property `S_c`, drawn
per neuron when `mu_S_c` / `sigma_S_c` are set — so $\phi$ here is the
zero-centred function; see the activation section below.

$c_s$, $K_s$ and $\tau_{a,sk}$ are **per cell type**, not per route: adaptation is
a property of the neuron, so it is indexed by the neuron's own type $s$, whereas
depression and facilitation belong to the synapse and are indexed by the route
$q \to s$. In the class these are the `1 x C` rows `c`, `n_a` and the `1 x C`
cell array `tau_a`.

**The $1/K_s$ is not decoration.** Every $a_{ik}$ relaxes to the rate, so
$\sum_k a_{ik} \to K_s r_i$ whatever the individual $\tau_{a,sk}$. Dividing by
$K_s$ makes the steady-state adaptation exactly $c_s r_i$ — independent of how
many timescales carry it — so $c_s$ is the **total adaptation budget** and
changing $K_s$ alters the timescale *structure* without altering adaptation
*strength*. Implemented as `params.c_eff = obj.c ./ max(1, obj.n_a)`
(`SRNNCellTypePairs.m:659`); the `max` guards $K_s = 0$, where the term is absent
anyway.

Note also that $r_i$ drives the adaptation and the synaptic states, while
$\theta_j^{q\to s}$ — the synaptic output — is what the recurrent sum transmits.
Depression is presynaptic and multiplicative; it does **not** feed back into
$a$, $b$ or $g$.

### Synaptic output

$$
\theta_j^{q\to s}=r_j\left(\prod_{m=1}^{M_{qs}}b_{jm}^{q\to s}\right)
\left(\prod_{k=1}^{L_{qs}}g_{jk}^{q\to s}\right).
$$

**This is the same $\theta$ as in
[`Equations_stability_paper.md`](Equations_stability_paper.md)**, carrying a
route superscript. The paper writes $\theta_j = r_j \prod_m b_{jm}$; that is
this expression in the special case of a single route with no facilitation,
where the $g$ product is empty and equal to one. Reading the two documents
together, $\theta$ always means *the rate after the synapse has acted on it* —
the quantity a presynaptic neuron actually delivers, as opposed to the rate
$r_j$ it fires at.

The superscript is **directed**: $\theta_j^{q\to s}$ is what neuron $j$ of type
$q$ delivers *to* a postsynaptic neuron of type $s$. One presynaptic neuron
therefore has as many synaptic outputs as there are configured target types, and
they differ whenever the routes carry different depression or facilitation. That
is exactly what `SRNNCellTypePairs` buys over a single per-neuron $\theta$, and
why `plot_data.synaptic_output` is nested by route (`.E.PV`) rather than being
one trace per neuron.

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

Unlike adaptation, depression is **not** normalised by its timescale count,
because it enters as a product rather than a sum: each $b$ rests at 1, so a
second timescale squares the steady-state depression rather than subdividing it.
That is deliberate — see
[`MTS_STD/MTS_STD_product_form_assessment.md`](MTS_STD/MTS_STD_product_form_assessment.md).

## Additive Wiener noise on `x`

The $dW_i$ are increments of $n$ independent standard Wiener processes, one per
neuron, so the model is an Itô SDE whenever `sigma_u_noise` ($\sigma_u$) is
nonzero, and collapses to a deterministic ODE at $\sigma_u = 0$.

Noise enters **only** $x$. The diffusion coefficient is therefore constant,
which is what makes Itô and Stratonovich agree, kills the Milstein correction
term, leaves the variational equation used by the QR Lyapunov method untouched,
and cancels in Benettin's trajectory difference (both trajectories are
integrated against the same noise path), so the LLE stays measurable at any
noise level.

| Property | Role |
|---|---|
| `sigma_u_noise` ($\sigma_u$) | Input-referred noise amplitude, in the units of $u$, so it is directly comparable to `intrinsic_drive` and the stimulus amplitude. Settable; default `0`. |
| `sigma_x_raw` ($\sigma_u/\tau_d$) | Coefficient the integrator multiplies $dW_i$ by. Dependent, read-only. |
| `x_noise_std` ($\sigma_u/\sqrt{2\tau_d}$) | Nominal stationary std of $x$ from noise alone ($u=0$, $W=0$). Dependent, read-only. |
| `noise_seed` | Pins the increment draw; empty derives it from `rng_seeds(1)`. |

`sigma_u_noise > 0` requires a stochastic integrator — `ode_solver` must be
`'euler'`, `'heun'` or `'sra1'` (`'sra1'` preferred: same two drift evaluations
per step as Heun, strong order 1.5 instead of 1.0). The adaptive solvers cannot
step an SDE, and `'rk4'` is deliberately kept deterministic so $\sigma_u = 0$ work
stays bit-identical to earlier runs.

The increments are pre-generated for the whole run by `build_noise()` and
indexed by absolute time, so a re-integrated segment sees exactly the increments
the original run used. They are held only for the duration of `run()` and
`compute_lyapunov()`, then cleared — regenerable from `noise_seed`, and large
enough (~96 MB at $n=300$, $T=50$ s, $f_s=400$) not to be worth saving.

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

Connectivity is uniform in `alpha = indegree/n`. `mu_tilde` and `sigma_tilde`
are $C \times C$ blocks indexed **(postsynaptic, presynaptic)**, so E→E may
differ from E→I; a `1 x C` row is accepted and broadcast down the columns as a
*presynaptic* row, which is the form the paper's preset uses.

---

## Definition of the activation function $\phi(x)$

Selected by name — `activation` is `'logistic'` (the class default) `'piecewise'`
or `'tanh'`. The paper's preset uses `'piecewise'`, defined here.
Implementation: `src/nonlinearities/piecewiseSigmoid.m`.

**$\phi$ is defined here centred on zero**, with $\phi(0) = \tfrac{1}{2}$. The
setpoint does not appear below: it is subtracted in the rate equation as
$a_{0_i}$, which keeps the translation separate from the shape. The benefit is
that $\phi$ can then be swapped for any standard sigmoid — the plain logistic,
$\tanh$ — without each one needing its own shift parameter, and the setpoint
means the same thing whichever is chosen.

The code does it the other way round, passing the centre into the nonlinearity
(`piecewiseSigmoid(x, S_a, S_c)`, where the third argument is the horizontal
shift). The two agree because a translation commutes with everything else here:
$\phi_{S_c}(z) = \phi_0(z - S_c)$.

The activation function $\phi(x)$ is a **hard sigmoid with rounded corners** — a
piecewise function composed of five regions that smoothly transitions between
saturation at 0 and 1.

**Domain:** $(-\infty, \infty)$
**Range:** $[0, 1]$

**Auxiliary quantities:**

Let the half-width of the linear segment be $q = \frac{q_\phi}{2}$, and define the scaling constant for the quadratic regions:

$$k = \frac{1}{2(1 - 2q)}$$

The four breakpoints dividing the five regions are:

$$x_1 = q - 1, \quad x_2 = -q, \quad x_3 = q, \quad x_4 = 1 - q$$

**Piecewise definition:**

$$
\phi(x) =
\begin{cases}
0 & \text{if } x < x_1 \quad \text{(left saturation)} \\[4pt]
k(x - x_1)^2 & \text{if } x_1 \le x < x_2 \quad \text{(left quadratic)} \\[4pt]
x + \tfrac{1}{2} & \text{if } x_2 \le x \le x_3 \quad \text{(linear segment)} \\[4pt]
1 - k(x - x_4)^2 & \text{if } x_3 < x \le x_4 \quad \text{(right quadratic)} \\[4pt]
1 & \text{if } x > x_4 \quad \text{(right saturation)}
\end{cases}
$$

**Parameter interpretation:**

- **$q_\phi \in [0, 1]$** (`S_a`): Controls the width of the central linear region. When $q_\phi = 1$, the function reduces to a standard hard sigmoid (pure piecewise linear with corners). When $q_\phi = 0$, the function becomes purely quadratic transitions with no linear segment.

There is no centre parameter here — see the note above. The setpoint is
$a_{0_i}$ in the rate equation, held in the property `S_c`, and may be drawn
per neuron via `mu_S_c` / `sigma_S_c` into the read-only `S_c_vec`.

> **Note:** The quadratic segments ensure $C^1$ continuity (smooth corners) at the transitions between the linear and saturated regions, avoiding the discontinuous derivatives of a standard hard sigmoid.

---

## Derivation of the default scaling factor $F$

The five RMT connectivity parameters are set as multipliers of $F$ — see
`mu_E_tilde_relative` and friends — so $F$ is the unit in which weights are
expressed. It is derived from Harris (2023), Equations 16 and 18.

**Equation 16** — Sparse variance for population $k \in \{e, i\}$:

$$\sigma_{sk}^2 = \alpha(1-\alpha)\tilde{\mu}_k^2 + \alpha \tilde{\sigma}_k^2$$

**Equation 18** — Spectral radius:

$$\mathcal{R} = \sqrt{N[f \sigma_{se}^2 + (1-f)\sigma_{si}^2]}$$

**Derivation:** To find the scaling factor $F$ that yields a unit spectral radius ($\mathcal{R} = 1$) when all normalized parameters are equal, set $|\tilde{\mu}_E| = |\tilde{\mu}_I| = \tilde{\sigma}_E = \tilde{\sigma}_I = F$.

Substituting into Eq. 16:

$$\sigma_{se}^2 = \sigma_{si}^2 = \alpha(1-\alpha)F^2 + \alpha F^2 = F^2[\alpha - \alpha^2 + \alpha] = F^2 \cdot \alpha(2-\alpha)$$

Substituting into Eq. 18, noting that $f + (1-f) = 1$:

$$\mathcal{R} = \sqrt{N \cdot F^2 \cdot \alpha(2-\alpha)} = 1$$

Solving for $F$:

$$F = \frac{1}{\sqrt{N\alpha(2-\alpha)}}$$

`F_tracks_network` (default `true`) recomputes $F$ from the current $n$ and
`indegree`, which makes $\mathcal{R}$ exactly independent of network size — the
$n\alpha$ cancels in `get.R`. Setting it `false` pins $F$ to a reference
$(n, \text{indegree})$, so the weight distribution stays fixed while
$\mathcal{R}$ varies with size.

> **Note on Zero Row Sum (ZRS):** Harris (2023) describes a Zero Row Sum condition (ZRS/SZRS) that controls "local" eigenvalue outliers escaping the spectral disc. In these simulations, we deliberately did not apply the ZRS condition in order to test whether adaptation mechanisms (SFA and STD) could fulfill a similar stabilizing role—effectively examining if adaptation can substitute for ZRS in constraining network dynamics.

Further RMT background: [`../RandomMatrixTheoryDocs/RMT_notes.md`](../RandomMatrixTheoryDocs/RMT_notes.md).

---

## Where the values are

This file states **forms**, not the numbers any particular run used. For those,
read the generated
`scripts/presentations/Stability_Manuscript/doc_equations_table/equation_table.md`,
rebuilt from a live model object on every `make_all_paper_figures` run, and its
sibling `adaptation_conditions.md`.
