$$
\frac{dx_i}{dt} = \frac{-x_i + u_i + \sum_{j=1}^{J} w_{ij}\, b_j r_{j}}{\tau_d}
$$

$$
r_i = \phi\left( x_i - a_{0_i} - c \sum_{k=1}^{K} a_{ik} \right)
$$

$$
\frac{da_{ik}}{dt} = \frac{-a_{ik} + r_i}{\tau_{k}}
$$

$$
\frac{db_i}{dt} = \frac{1-b_i}{\tau_{rec}} - \frac{b_i\, r_i}{\tau_{rel}}
$$