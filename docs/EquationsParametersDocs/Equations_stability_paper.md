# SRNN model equations

> **This is the single human-written statement of the model.** The equations were
> previously copied by hand into eight files across five folders, which drifted
> apart — each copy carried a different subset of the corrections made during the
> `c/K` refactor. Everything else now links here.
>
> The one other place the equations appear is
> `scripts/presentations/Stability_Manuscript/doc_equations_table/equation_table.md`,
> which is **generated** by `write_manuscript_tables` from a live model object on
> every `make_all_paper_figures` run and therefore cannot drift from the code. It
> also reports the values a given preset actually runs at. This file states the
> model; that one states a particular instance of it.
>
> The implementation is `SRNNCellTypePairs.dynamics_fast` (`src/SRNNCellTypePairs.m`).

## The model

$$
\begin{aligned}
dx_i &= \frac{-x_i + u_i + \sum_{j=1}^{N} w_{ij}\, \theta_j}{\tau_d}\, dt \;+\; \frac{\sigma_u}{\tau_d}\, dW_i \\[8pt]
\theta_i &= r_{i} \prod_{m=1}^{M} b_{im} \\[8pt]
r_i &= \phi\left( x_i - \frac{c}{K} \sum_{k=1}^{K} a_{ik} \right) \\[8pt]
\frac{da_{ik}}{dt} &= \frac{-a_{ik} + r_i}{\tau_{a_k}}, \qquad k = 1, \dots, K \\[8pt]
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

**Adaptation is normalized by $K$, depression is not.** Each $a_{ik}$ relaxes to
the rate, so $a_{ik} \to r_i$ for every timescale whatever its $\tau_{a_k}$, and
$\sum_k a_{ik} \to K r_i$. Dividing by $K$ therefore makes the steady-state
adaptation $c\, r_i$ exactly — independent of how many timescales carry it — so
$c$ is the **total adaptation budget** and changing $K$ changes the timescale
*structure* without also changing adaptation *strength*.

Depression needs no such factor because it enters as a **product** rather than a
sum: each $b_{im}$ rests at 1, so adding a timescale multiplies rather than
subdividing. With $M = 2$ sharing a common $\tau_{rec}/\tau_{rel}$ ratio the
steady state is the square of the single-timescale value, which is deliberate.
See `MTS_STD/MTS_STD_product_form_assessment.md` for why the product form is
defensible and what it costs.

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

## Notes on terms that are *not* in the model

**There is no $a_0$ offset.** Older derivations wrote the argument of $\phi$ as
$x_i - a_{0_i} - c\sum_k a_{ik}$, with $a_{0_i}$ standing for a threshold shift
or preferred input. **No model class implements it**:
`SRNNCellTypePairs.dynamics_fast` forms `x_eff = x - c_eff .* sum(a, 2)`, and
there is no $a_0$ anywhere in `src/`. It is also redundant by construction — a
constant subtracted inside $\phi$ is indistinguishable from a constant added to
$u_i$, which is exactly what `input_config.intrinsic_drive` supplies. So this is
not a missing feature to implement; it is a symbol the model never needed.

**Noise enters only $x$.** That is what keeps the diffusion constant, which in
turn makes Itô and Stratonovich coincide, kills the Milstein term, leaves the QR
variational equation untouched, and makes the noise cancel in Benettin's
trajectory difference — so the largest Lyapunov exponent stays measurable at any
noise level. $\sigma_u$ is **input-referred**, in the units of $u$, so it is
directly comparable to `intrinsic_drive` and the stimulus amplitude.
