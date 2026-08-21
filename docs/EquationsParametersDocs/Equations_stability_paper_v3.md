$$
\begin{aligned}
dx_i &= \frac{-x_i + u_i + \sum_{j=1}^{N} w_{ij}\, \phi\left( x_j - a_{0_j} - c \sum_{k=1}^{K} a_{jk} \right) \prod_{m=1}^{M} b_{jm}}{\tau_d}\, dt \;+\; \frac{\sigma_u}{\tau_d}\, dW_i \\[8pt]
\frac{da_{ik}}{dt} &= \frac{-a_{ik} + \phi\left( x_i - a_{0_i} - c \sum_{k'=1}^{K} a_{ik'} \right)}{\tau_{a_k}}, \qquad k = 1, \dots, K \\[8pt]
\frac{db_{im}}{dt} &= \frac{1-b_{im}}{\tau_{rec_m}} - \frac{b_{im}\, \phi\left( x_i - a_{0_i} - c \sum_{k=1}^{K} a_{ik} \right)}{\tau_{rel_m}}, \qquad m = 1, \dots, M
\end{aligned}
$$
