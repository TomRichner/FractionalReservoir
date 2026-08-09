# SRNN System Equations — Multi-Timescale SFA + Multi-Timescale STD

Model equations with multi-timescale spike-frequency adaptation ($K$ timescales)
and multi-timescale short-term synaptic depression ($M$ timescales). Facilitation
is omitted.

## System Equations

$$
\frac{dx_i}{dt} = \frac{-x_i + u_i + \sum_{j=1}^{N} w_{ij}\, \left(\prod_{m=1}^{M} b_{jm}\right) r_{j}}{\tau_d}
$$

$$
r_i = \phi\left( x_i - a_{0_i} - c \sum_{k=1}^{K} a_{ik} \right)
$$

$$
\frac{da_{ik}}{dt} = \frac{-a_{ik} + r_i}{\tau_{a_k}}, \qquad k = 1, \dots, K
$$

$$
\frac{db_{im}}{dt} = \frac{1-b_{im}}{\tau_{rec_m}} - \frac{b_{im}\, r_i}{\tau_{rel}}, \qquad m = 1, \dots, M
$$

---

**Abbreviations:**
- **SRNN**: Stable Recurrent Nonlinear Network
- **SFA**: Spike frequency adaptation
- **STD**: Short-term synaptic depression
