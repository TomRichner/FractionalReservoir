\begin{aligned}
\dot{x}_i(t) &= \frac{-x_i(t) + u_i(t) + \sum_{j=1}^{J} w_{ij}\, r_{j}(t - d)}{\tau_d}\\[4pt]
r_i(t) &= b_i(t)\,\phi\!\left(
            x_i(t) - a_{0_i}
            - c \sum_{k=1}^{K} a_{ik}(t)
        \right)\\[4pt]
\dot{a}_{ik}(t) &= \frac{-a_{ik}(t) + r_i(t)}{\tau_{k}}\\[4pt]
\dot{b}_i(t) &= \frac{1-b_i(t)}{\tau_{rec}}
            - \frac{b_i(t)\, r_i(t)}{\tau_{rel}}
\end{aligned}
