## Table 2: Adaptation Conditions

| Condition     | $n_{a_E}$ | $n_{a_I}$ | $n_{b_E}$ | $n_{b_I}$ | Description                        |
| :------------ | :-------- | :-------- | :-------- | :-------- | :--------------------------------- |
| No Adaptation | 0         | 0         | 0         | 0         | Baseline                           |
| SFA Only      | 3         | 0         | 0         | 0         | Spike-frequency adaptation enabled |
| STD Only      | 0         | 0         | 1         | 0         | Short-term depression enabled      |
| SFA + STD     | 3         | 0         | 1         | 0         | Both mechanisms enabled            |

**Effect on parameters:**

- When $n_{a_E} = 0$: No SFA variables $a_{ik}$; the $c_E \sum_k a_{ik}$ term is zero.
- When $n_{a_E} = 3$: Three SFA timescales with $\tau_{a_k} \in \{0.25, 1.58, 10\}$ s (logspaced) and coupling $c_E = \frac{0.5}{3}$.
- When $n_{b_E} = 0$: No STD variable $b_i$; synaptic output equals $r_i$ (equivalent to $b_i = 1$).
- When $n_{b_E} = 1$: STD enabled with a single timescale, recovery constant $\tau_{rec} = 1$ s and release constant $\tau_{rel} = \frac{1}{4}$ s; the resource enters the recurrent term as the product $b_i r_i$.

Inhibitory neurons have no adaptation mechanisms ($n_{a_I} = n_{b_I} = 0$).

> **Implementation note:** When any of $n_{a_E}$, $n_{a_I}$, $n_{b_E}$, or $n_{b_I}$ is set to zero, the corresponding state variables ($a_{ik}$ or $b_i$) are excluded from the system state vector and the Jacobian matrix. This prevents spurious zero eigenvalues that would otherwise arise from including disabled adaptation dynamics.
