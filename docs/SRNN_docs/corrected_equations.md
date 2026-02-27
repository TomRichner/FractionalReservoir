$$
\begin{aligned}
\dot{x}_i &= \frac{-x_i + u_i + \sum_{j=1}^{J} w_{ij}\, b_j r_{j}}{\tau_d}\\[4pt]
r_i &= \phi\!\left(
            x_i - a_{0_i}
            - c \sum_{k=1}^{K} a_{ik}
        \right)\\[4pt]
\dot{a}_{ik} &= \frac{-a_{ik} + r_i}{\tau_{k}}\\[4pt]
\dot{b}_i &= \frac{1-b_i}{\tau_{rec}}
            - \frac{b_i\, r_i}{\tau_{rel}}
\end{aligned}
$$