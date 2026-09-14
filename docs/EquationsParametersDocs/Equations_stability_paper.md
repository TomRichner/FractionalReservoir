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

Facilitation is configured per presynaptic-to-postsynaptic route rather than per
neuron — `synapse_config.<pre>.<post>.stf`, with fields `tau_dec`, `tau_fac` and
`G`. See [`cell_type_pair_equations.md`](cell_type_pair_equations.md) for the
per-route form, in which
$b$ and $g$ carry route superscripts.

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
