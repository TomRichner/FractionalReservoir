# SRNN model equations

> **This is the single human-written statement of the model.** The equations were
> previously copied by hand into eight files across five folders, which drifted
> apart — each copy carried a different subset of the corrections made during the
> `c/K` refactor. Everything else now links here.
>
> The equations also appear in **generated** form, as `equation_table.md`, which
> `make_all_paper_figures` writes into its figure root — `figs/paper/doc_tables/`
> by default. It is built from a live model object on every run and therefore
> cannot drift from the code, and it reports the values a given preset actually
> runs at. This file states the model; that one states a particular instance of
> it.
>
> That copy is **not committed** — `figs/` is gitignored, so it exists once you
> have run the pipeline. It used to be tracked, which meant a clone could hold a
> table describing a preset the code no longer used: exactly the drift the
> generated table exists to prevent.
>
> The implementation is `SRNNCellTypePairs.dynamics_fast` (`src/model/SRNNCellTypePairs.m`).
> Since 2026-11-11 the depression and facilitation equations of a route live in
> a **route model** under `src/model/routes/`: `route_product` (the product form
> below, and the facilitation further down) and `route_pool` (the chained
> vesicle pools, in their own section). `dynamics_fast` forms the rate, the SFA
> equations and the dendritic equation, and asks each route's model for its state
> derivative and its output factor.

## The model

$$
\begin{aligned}
dx_i &= \frac{-x_i + u_i + \sum_{j=1}^{N} w_{ij}\, \theta_j}{\tau_d}\, dt \;+\; \frac{\sigma_u}{\tau_d}\, dW_i \\[8pt]
\theta_i &= r_{i} \prod_{m=1}^{M} b_{im} \\[8pt]
r_i &= \phi\left( x_i - a_{0_i} - \frac{c}{K} \sum_{k=1}^{K} a_{ik} \right) \\[8pt]
\frac{da_{ik}}{dt} &= \frac{-a_{ik} + r_i}{\tau_{a,ik}}, \qquad k = 1, \dots, K \\[8pt]
\frac{db_{im}}{dt} &= \frac{1-b_{im}}{\tau_{rec_m}} - \frac{b_{im}\, r_i}{\tau_{rel_m}}, \qquad m = 1, \dots, M
\end{aligned}
$$

$\theta_i$ is the **synaptic output**: the rate after depression, and the
quantity the recurrent sum actually transmits.

Note the placement of the depression factor. The rate $r_i$ is the
**pre-depression** output of the nonlinearity, and depression enters
presynaptically as the product $\theta_j$ in the recurrent sum. Both SFA and STD
are therefore driven by the raw rate $r_i$, not by $\theta_i$. This is not
cosmetic: the alternative framing $r_i = b_i\,\phi(\cdot)$ would make SFA
integrate $b_i r_i$, make the STD equation depend on $b_i^2 r_i$, and put a
factor of $b$ into the $a \to x$ and $a \to a$ Jacobian blocks.

**The SFA timescales may be per neuron.** By default every neuron of a cell
type shares the type's ladder, $\tau_{a,ik} = \tau_{a_k}$. With a nonzero
`tau_a_spread` $\sigma_q$ (one dimensionless number per type, zero in every
preset unless its name says otherwise) `build()` draws each neuron its own
ladder from the nominal one:

$$
\log \tau_{a,ik} = \log \tau_{a_k} + (1 - w_k)\, \delta_i^{\mathrm{fast}} + w_k\, \delta_i^{\mathrm{slow}},
\qquad \delta_i^{\mathrm{fast}}, \delta_i^{\mathrm{slow}} \sim \mathcal{N}(0, \sigma_q^2),
\qquad w_k = \frac{\log \tau_{a_k} - \log \tau_{a_1}}{\log \tau_{a_K} - \log \tau_{a_1}} .
$$

The two ends of the ladder are jittered log-normally (median-preserving, so
$\tau$ stays positive and the linear spread is proportional to $\tau$) and the
interior rungs follow in log space at the nominal's own position; for a
log-spaced ladder this is "jitter the endpoints and re-run logspace". It exists
to break the degeneracy of the slow Lyapunov band (every neuron otherwise
contributes an exponent within a few percent of $-1/\tau_{a_K}$); note that
with a spread the leading exponent of a stable network is $-1/\tau$ of the
*slowest drawn neuron*. The draw is seeded, saved and restored around, and
recorded in the read-only `tau_a_matrix`. See `SRNNCellTypePairs.tau_a_spread`.

**Adaptation is normalized by $K$, depression is not.** Each $a_{ik}$ relaxes to
the rate, so $a_{ik} \to r_i$ for every timescale whatever its $\tau_{a,ik}$, and
$\sum_k a_{ik} \to K r_i$. Dividing by $K$ therefore makes the steady-state
adaptation $c\, r_i$ exactly — independent of how many timescales carry it — so
$c$ is the **total adaptation budget** and changing $K$ changes the timescale
*structure* without also changing adaptation *strength*.

Depression needs no such factor because it enters as a **product** rather than a
sum: each $b_{im}$ rests at 1, so adding a timescale multiplies rather than
subdividing. With $M = 2$ sharing a common $\tau_{rec}/\tau_{rel}$ ratio the
steady state is the square of the single-timescale value, which is deliberate.
See [`MTS_STD/MTS_STD_notes.md`](MTS_STD/MTS_STD_notes.md) for why the product
form is defensible, what it costs, and why there is deliberately **no**
depression budget to conserve. (It supersedes the earlier
`MTS_STD_product_form_assessment.md` in the same folder, which recommended
splitting a budget across timescales — that recommendation was rejected.)

### 2026-09-13: STD strength matching between one and two timescales

**The problem.** The paragraph above is still how the *unmatched* preset
(`celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25`) is
built, and it confounds two things. With the steady state of one depression
variable $b_m(r) = 1/(1 + r/\rho_m)$, $\rho_m = \tau_{rel,m}/\tau_{rec,m}$, and
both timescales at $\rho = 0.125$ ($\tau_{rec} = [2, 4]$ s, $\tau_{rel} =
[0.25, 0.5]$ s; the single-timescale routes at $(2, 0.25)$), the
two-timescale product is the **square** of the one-timescale factor at every
rate. Going from one timescale to two therefore changes the *number of recovery
timescales* and the *steady-state strength* of depression at once, so the
single-vs-multiple comparison could not attribute a difference to timescale
count alone. SFA does not have this problem because of $c/K$.

**What is matched.** The steady-state synaptic output $\theta_{ss}(r) = r
\prod_m b_m(r)$ of a two-timescale route is made equal to that of a
one-timescale route at one reference rate,

$$ r_{ref} = 0.25, $$

the median mean firing rate of the multiple-timescale condition at the default
point of the five well-defined 1-D sweeps of the medium run `data/topk_med`
(per sweep 0.223–0.251, pooled over 105 reps 0.237 with 5th–95th percentiles
0.050–0.346), rounded to 0.05. At $r_{ref}$ one factor is $1/(1 +
0.25/0.125) = 1/3$ and the unmatched product is $1/9$. Exact matching over
every rate is impossible for an unweighted product with a different number of
factors; the match is at $r_{ref}$ and the mismatch elsewhere is what
`fig_STD_steady_state` draws over the occupied rates.

**Variant 1, route scale (primary; TR's decision, 2026-09-13).** Preset
`..._dualStdScaled_3cond_mu8p25`. The depression variables and $\theta_j$ are
unchanged; the two-timescale routes carry a per-route weight scale $s$, so the
recurrent sum a postsynaptic neuron receives is

$$ \sum_j s_{q_j \to q_i}\, w_{ij}\, \theta_j, \qquad s = \frac{\theta^{(1)}_{ss}(r_{ref})}{\theta^{(2)}_{ss}(r_{ref})} = 1 + \frac{r_{ref}}{\rho} = 3 $$

on all four routes of the `sfa3_std2` condition and $s = 1$ everywhere else. In
the code this is `synapse_config.<pre>.<post>.scale`, folded into `params.W`
inside `SRNNCellTypePairs.get_params` (the drawn `W` is untouched, so the
shared-build check and `plot_W` still describe the drawn matrix; the
dynamics, both Jacobians and `jacobian_times` all read `params.W`, so they
stay consistent by construction). **Consequence on record:** because
$b_m \to 1$ as $r \to 0$, the *undepressed* recurrent gain of the
two-timescale condition is $3\times$ that of the other two conditions; the
match holds at the operating point, not at low rates. TR chose this with that
consequence stated; the second variant is the control for it.

**Variant 2, usage matching (control).** Preset
`..._dualStdUsage_3cond_mu8p25`. No scale; $\tau_{rec} = [2, 4]$ s kept, and
the usage $\rho_u = \tau_{rel}/\tau_{rec}$ set equal on both timescales so
that $(1 + r_{ref}/\rho_u)^2 = 1 + r_{ref}/\rho$:

$$ \rho_u = \frac{r_{ref}}{\sqrt{1 + r_{ref}/\rho} - 1} = \frac{0.25}{\sqrt{3} - 1} = 0.34151, \qquad \tau_{rel} = \rho_u\, \tau_{rec} = [0.68301,\ 1.36603]\ \text{s}. $$

The low-rate gain is unchanged, each timescale depresses less, and the
depression dynamics are slower (rate $1/\tau_{rec} + r/\tau_{rel}$ per
variable). $\tau_{rel}$ is Varela's $d_m$, so this is the tuning the product
model itself offers.

In both variants `no_adaptation` and `sfa1_std1` are identical to the unmatched
preset's; only `sfa3_std2` differs. The full rationale, the numbers and the
replacement manuscript text are in
[`../notes/STD_strength_matching_2026-09-13.md`](../notes/STD_strength_matching_2026-09-13.md);
`test_route_scale` checks the presets against the formulas above.

**Option, geometric-mean combine (2026-09-24; opt-in, no preset uses it).**
Ported from train-srnn's geometric-mean STD matching. A route may set
`synapse_config.<pre>.<post>.std.combine = 'geomean'` (default `'product'`,
bit-identical to the form above), which replaces the route's depression
factor by

$$ \theta_j = r_j \Big(\prod_{m=1}^{M} b_{jm}\Big)^{1/M} \prod_n g_{jn}. $$

When every pair shares one $\rho$ the steady state is $1/(1 + r/\rho)$ at
**every** constant rate, not only at $r_{ref}$: it is the depression analogue
of $c/K$, and leaves the low-rate gain and the $b_m$ ODEs unchanged. With
`std_zero_floor` the floor is applied to the combined factor (its $p_{min}$ is
combined the same way). `SRNNCellTypePairs.combine_std_pair` is the one
definition; the dynamics, `plot_data.synaptic_output`, both Jacobians and
`jacobian_times` all go through it (`test_std_geomean`).

## Chained vesicle-pool depression (opt-in per route)

Added 2026-11-12 (FD-004). **This is the model the final draft's
multiple-timescale condition uses** (condition `sfa3_stdc2`); the product form
above stays the default for a route and is kept as a fourth condition for
comparison. A route asks for it with

```matlab
synapse_config.<pre>.<post>.std = struct('model', 'chain', ...
    'tau_rec', [2 5], ...   % [releasable-pool refill, reserve refill], s
    'tau_rel', 0.25, ...    % ONE release time: only the releasable pool releases
    'n2', 2);               % reserve capacity / releasable-pool capacity
```

and then carries two variables per presynaptic neuron, a **releasable pool**
$b_{i1}$ and a **reserve** $b_{i2}$, in series:

$$
\begin{aligned}
\frac{db_{i1}}{dt} &= \frac{b_{i2}\,(1 - b_{i1})}{\tau_{rec_1}} - \frac{b_{i1}\, r_i}{\tau_{rel}} \\[8pt]
\frac{db_{i2}}{dt} &= \frac{1 - b_{i2}}{\tau_{rec_2}} - \frac{1}{n_2}\,\frac{b_{i2}\,(1 - b_{i1})}{\tau_{rec_1}} \\[8pt]
\theta_i &= r_i\, b_{i1}
\end{aligned}
$$

Only the releasable pool is released from, and only it sets the synaptic
output. It is refilled from the reserve at a rate proportional to how full the
reserve is; the reserve is refilled toward full with its own, slower time
constant and is drawn down by exactly what it hands on, divided by its relative
capacity $n_2$. Both start at 1 and both stay in $[0, 1]$. As everywhere in this
model, the drive is the raw rate $r_i$.

**The single-timescale condition is this model with the reserve removed, not
clamped.** With $b_{i2} \equiv 1$ the first equation is the one-component
depression equation of the previous section, so a route with one depression
component *is* the chain without a reserve, and it carries one state variable,
not two. No route keeps an unused or clamped variable: it would enlarge the
state and add a zero eigenvalue and a zero Lyapunov exponent that mean nothing.
So 1TS and MTS share the releasable-pool equation and differ by one mechanism.
The reserve starts full, so the two coincide at the onset of activity and come
apart as sustained activity draws the reserve down. They are **not** matched
in steady-state depth, by decision (to-do list, block A).

**Steady state at constant rate $r$.** With
$k_1 = r\,\tau_{rec_1}/\tau_{rel}$, $k_2 = r\,\tau_{rec_2}/(n_2\,\tau_{rel})$ and
$s = 1 + k_1 + k_2$,

$$
b_1 = \frac{2}{s + \sqrt{s^2 - 4k_2}}, \qquad b_2 = 1 - k_2\, b_1 .
$$

This reduces to the single pool's $1/(1 + k_1)$ as $n_2 \to \infty$, and the
sustained output $r\, b_1$ rises with rate to a ceiling of
$\big(\tau_{rec_1}/\tau_{rel} + \tau_{rec_2}/(n_2\,\tau_{rel})\big)^{-1}$, where the
product form's output peaks and then falls. With the chosen values
($\tau_{rec} = [2, 5]$ s, $\tau_{rel} = 0.25$ s, $n_2 = 2$): $b_1$ is 0.420,
0.200, 0.106, 0.054 at $r$ = 0.1, 0.25, 0.5, 1 (single pool: 0.556, 0.333,
0.200, 0.111; product with $[2, 4]$ and $[0.25, 0.5]$ s: 0.309, 0.111, 0.040,
0.012), and the ceiling is $1/18$. `route_pool`'s `steady_state` is the one
place this closed form is coded (in a form that avoids cancellation when the
reserve is nearly empty); figures and tables call it.

`std.combine`, the zero floor (`std_zero_floor`) and facilitation are refused
on a chain route; the route `scale` is allowed. The depression factor of any
route, product or chain, is `SRNNCellTypePairs.std_route_factor`.

## Facilitation (optional)

`SRNNCellTypePairs` also supports short-term facilitation, per route. The paper's
preset carries depression only, so the equations above omit it — but the
mechanism is implemented, not missing. With $n_f$ facilitation timescales the
synaptic output gains a second product,

$$
\theta_i = r_i \left( \prod_{m=1}^{M} b_{im} \right) \left( \prod_{n=1}^{n_f} g_{in} \right)
$$

and each facilitation variable follows

$$
\frac{dg_{in}}{dt} = \frac{1-g_{in}}{\tau_{dec_n}} + \frac{(G - g_{in})\, r_i}{\tau_{fac_n}}, \qquad n = 1, \dots, n_f
$$

where $G$ is the ceiling the facilitated gain approaches. Like depression,
facilitation is driven by the raw rate $r_i$ and enters presynaptically; each
$g_{in}$ rests at 1, so an absent mechanism contributes an empty product of one.

That is the **output-only** form (`stf.model` absent or `'product'`): $g$
multiplies what the synapse transmits and nothing else, so a facilitated synapse
transmits more without using its resources faster.

### Release-coupled facilitation on the vesicle pools (opt-in per route)

Added 2026-11-12 (FD-018). **Not used in the current paper**; built so the pool
model is complete. On a route served by the pool model, facilitation is a
**release probability** $g_i$ relative to rest, which multiplies both the
synaptic output and the rate at which the releasable pool is used:

$$
\begin{aligned}
\frac{db_{i1}}{dt} &= \frac{b_{i2}\,(1 - b_{i1})}{\tau_{rec_1}} - \frac{g_i\, b_{i1}\, r_i}{\tau_{rel}} \\[8pt]
\frac{dg_i}{dt} &= \frac{1 - g_i}{\tau_{dec}} + \frac{(G - g_i)\, d_i}{\tau_{fac}} \\[8pt]
\frac{dh_i}{dt} &= \frac{r_i - h_i}{\tau_h} \\[8pt]
\theta_i &= r_i\, g_i\, b_{i1}
\end{aligned}
$$

with the reserve equation unchanged. The drive $d_i$ is the **calcium trace**
$h_i$ when the route has one (`stf.tau_h` given), and the rate $r_i$ when it does
not; in that single-variable form $h$ is not carried at all. So facilitation
nests as depression does: the one-variable form is the equation of the previous
paragraph plus the coupling, and the two-variable form adds one variable that
changes one term. Likewise a route may have no reserve ($b_{i2}$ dropped, the
refill as if it were 1) or no depression at all (then $\theta_i = r_i g_i$).
$g$ starts at 1 and stays in $[1, G]$; $h$ starts at 0 and stays in $[0, 1]$.

```matlab
synapse_config.<pre>.<post>.stf = struct('model', 'chain', ...
    'tau_dec', 1.5, 'tau_fac', 0.25, 'G', 5, ...   % G >= 1: the ceiling of g
    'tau_h', 0.5);                                 % omit for the one-variable form
```

**Steady state at constant rate.** $h = r$ and
$g = (1 + G\kappa r)/(1 + \kappa r)$ with $\kappa = \tau_{dec}/\tau_{fac}$, with or
without the calcium trace; the pools are the chain's steady state evaluated at
the **effective rate** $g\,r$. So at steady state coupled facilitation stretches
the rate axis and nothing else: it moves the sustained output toward the
depression ceiling at lower rates and does not raise the ceiling. The calcium
trace changes the time course only (it delays $g$, and removes the overshoot of
the output after a rate step). A coupled `stf` requires the route's depression to
be absent, a single pool or the chain; it is refused beside a product of two or
more pools, `std.combine` or the zero floor, and a chain route refuses the
output-only form.

Facilitation is configured per presynaptic-to-postsynaptic route rather than per
neuron — `synapse_config.<pre>.<post>.stf`, with fields `tau_dec`, `tau_fac` and
`G`. See [`cell_type_pair_equations.md`](cell_type_pair_equations.md) for the
per-route form, in which
$b$ and $g$ carry route superscripts.

## Per-neuron synapse parameters (optional)

The depression and facilitation constants above are written per route. With
`std_spread` or `stf_spread` (one dimensionless number $\sigma_q$ per
presynaptic cell type $q$, default 0) each presynaptic neuron $j$ of that type
has its own values, the counterpart of the per-neuron SFA ladders:

$$
\tau_{j} = \tau \, e^{\sigma_q z_{j}}, \qquad
n_{2,j} = n_2 \, e^{\sigma_q z_{j}'}, \qquad
G_j = 1 + (G - 1)\, e^{\sigma_q z_{j}''}, \qquad z \sim \mathcal{N}(0, 1),
$$

with an independent standard normal for each quantity (every refill and release
time of every pool, the reserve size, and for facilitation $\tau_{dec}$,
$\tau_{fac}$, $\tau_h$ and the excess of the ceiling over 1). The nominal value
is the median, the values stay positive, and $G_j \ge 1$. In every equation of
the three sections above the constant is then read per presynaptic neuron; the
equations are otherwise unchanged, and with $\sigma_q = 0$ every neuron has the
nominal value exactly.

**A neuron's draws are shared.** One row of normals is drawn per neuron and
mechanism, and the column a quantity reads is fixed and the same for every
route model. So a neuron has the same values on all its routes, and on one
network seed it has the same releasable-pool refill and release time whether
its routes carry a single pool, the chain or a product of pools, and the same
normal behind its second refill time in the last two. That is what keeps the
single-timescale and multiple-timescale conditions comparable neuron by neuron.

**The reserve stays the slower pool.** On a chain route whose nominal reserve
refill is slower than the releasable pool's, a neuron whose drawn
$\log(\tau_{rec2,j} / \tau_{rec1,j})$ falls below a quarter of the nominal
$\log(\tau_{rec2} / \tau_{rec1})$ has $\tau_{rec2,j}$ raised to that floor;
$\tau_{rec1,j}$ is never moved, since it must equal the single pool's. With
$\tau_{rec} = [2, 5]$ s and $\sigma = 0.25$ this moves about 2.6 % of neurons,
and the expected reserve refill is 5.17 s (5.16 s from the log-normal alone).
The components of a product route are jittered independently and have no
ordering to keep.

The draws are `std_draws` and `stf_draws` on `SRNNCellTypePairs`, made in
`build()` from their own seeded streams; the realised values are in
`get_params().route{pre, post}.P`.

## Notes

**$a_{0_i}$ is the setpoint of the nonlinearity**, i.e. the property `S_c`, and
$\phi$ above is the **zero-centred** function. The code writes it the other way
round — `dynamics_fast` forms `x_eff = x - c_eff .* sum(a, 2)` and then applies a
$\phi$ that is itself centred at `S_c`, via `piecewiseSigmoid(x, S_a, S_c)`.
Since $\phi_{S_c}(z) = \phi_0(z - S_c)$ the two are identical; pulling the
setpoint out is what makes it visible as the swept parameter it is rather than
hidden inside the symbol $\phi$. It carries a neuron index because the setpoint
may be drawn per neuron: set `mu_S_c` / `sigma_S_c` and `build()` fills `S_c_vec`.

> A revision of this file dated 2026-08-27 deleted this term, claiming no model
> class implemented it and that it was "a symbol the model never needed". That
> was wrong. The mistake was checking `src/` for the *name* `a_0`, finding only
> `abscissa_0`, and concluding the *quantity* was absent — when the code calls it
> `S_c`. Restored 2026-08-28.

**Noise enters only $x$.** That is what keeps the diffusion constant, which in
turn makes Itô and Stratonovich coincide, kills the Milstein term, leaves the QR
variational equation untouched, and makes the noise cancel in Benettin's
trajectory difference — so the largest Lyapunov exponent stays measurable at any
noise level. $\sigma_u$ is **input-referred**, in the units of $u$, so it is
directly comparable to `intrinsic_drive` and the stimulus amplitude.
