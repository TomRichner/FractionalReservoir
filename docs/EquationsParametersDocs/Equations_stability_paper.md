$$
\begin{aligned}
dx_i &= \frac{-x_i + u_i + \sum_{j=1}^{N} w_{ij}\, \theta_j}{\tau_d}\, dt \;+\; \frac{\sigma_u}{\tau_d}\, dW_i \\[8pt]
\theta_i &= r_{i} \prod_{m=1}^{M} b_{im} \\[8pt]
r_i &= \phi\left( x_i - \frac{c}{K} \sum_{k=1}^{K} a_{ik} \right) \\[8pt]
\frac{da_{ik}}{dt} &= \frac{-a_{ik} + r_i}{\tau_{a_k}}, \qquad k = 1, \dots, K \\[8pt]
\frac{db_{im}}{dt} &= \frac{1-b_{im}}{\tau_{rec_m}} - \frac{b_{im}\, r_i}{\tau_{rel_m}}, \qquad m = 1, \dots, M
\end{aligned}
$$

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

**No $a_{0}$ offset.** Earlier notes in this folder — `J_eff_notes.md` and
`MTS_STD/system_equations_MTS_SFA_MTS_STD.md` — write the argument of $\phi$ as
$x_i - a_{0_i} - c\sum_k a_{ik}$, with $a_{0_i}$ an offset standing for a
threshold shift or preferred input. **Neither model class implements it**:
`SRNNCellTypePairs.dynamics_fast` forms `x_eff = x - c_eff .* sum(a, 2)` and
there is no $a_0$ anywhere in `src/`. A constant offset would in any case be
indistinguishable from a constant term in $u_i$, which is what
`input_config.intrinsic_drive` provides. The term was carried here by copying and
has been removed; the older notes still show it, and also predate the $c/K$
normalisation, so read them as history rather than as the current model.

**Facilitation is omitted, not absent.** `SRNNCellTypePairs` supports per-route
STF, which contributes a further presynaptic factor $\prod_m g_{im}$ alongside
$\prod_m b_{im}$. The paper's preset carries depression only, so $\theta$ is
written above without it.
