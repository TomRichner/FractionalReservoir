$$
\begin{aligned}
dx_i &= \frac{-x_i + u_i + \sum_{j=1}^{N} w_{ij}\, \theta_j}{\tau_d}\, dt \;+\; \frac{\sigma_u}{\tau_d}\, dW_i \\[8pt]
\theta_i &= r_{i} \prod_{m=1}^{M} b_{im} \\[8pt]
r_i &= \phi\left( x_i - a_{0_i} - \frac{c}{K} \sum_{k=1}^{K} a_{ik} \right) \\[8pt]
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
