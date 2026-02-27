# Delay Equations Draft — What the Code Actually Implements

This document is a line-by-line translation of `SRNN_reservoir_DDE.m` and its caller `SRNN_ESN.m` into mathematical equations. The goal is to faithfully document what is **actually numerically integrated**, not what was intended.

---

## 1. The Full System of Delay Differential Equations

The system actually being integrated is:

$$
\boxed{
\begin{aligned}
\dot{x}_i(t) &= \frac{1}{\tau_d}\!\left[-x_i(t) + u_i(t) + \sum_{j \in \mathcal{E}} w_{ij}\, b_j(t)\, r_j(t) + \sum_{j \in \mathcal{I}} w_{ij}\, b_j(t-\tau_{\text{syn}})\, r_j(t-\tau_{\text{syn}})\right] \\[6pt]
r_i(t) &= \phi\!\left(x_i(t) - c_P \sum_{k=1}^{K_P} a_{ik}(t)\right), \quad P \in \{E, I\} \\[6pt]
\dot{a}_{ik}(t) &= \frac{r_i(t) - a_{ik}(t)}{\tau_{a_k}}, \quad k = 1, \ldots, K_P \\[6pt]
\dot{b}_i(t) &= \frac{1-b_i(t)}{\tau_{\text{rec}}} - \frac{b_i(t)\, r_i(t)}{\tau_{\text{rel}}}
\end{aligned}
}
$$

where the delayed firing rate and depression in the $\dot{x}_i$ equation are defined as:

$$
r_j(t - \tau_{\text{syn}}) = \phi\!\left(x_j(t-\tau_{\text{syn}}) - c_P \sum_{k=1}^{K_P} a_{jk}(t-\tau_{\text{syn}})\right)
$$

$$
b_j(t - \tau_{\text{syn}}) \text{ is read directly from the delayed state vector}
$$

Here $P$ refers to the population of the presynaptic neuron $j$. Since the delayed sum runs over $j \in \mathcal{I}$, we have $P = I$ in that context. Note that when $K_I = 0$, the adaptation sum vanishes and $r_j = \phi(x_j)$. When $n_{b_I} = 0$, $b_j = 1$ and the depression factor drops out.

---

## 2. How the Weight Matrix is Split (`SRNN_ESN.m`, lines 137--152)

When a nonzero `lags` parameter is provided, the full recurrent weight matrix $W$ is decomposed into two components:

$$
W^{(\text{inst})}_{ij} = \begin{cases} w_{ij} & \text{if } j \in \mathcal{E} \\ 0 & \text{if } j \in \mathcal{I} \end{cases}
\qquad
W^{(\tau)}_{ij} = \begin{cases} 0 & \text{if } j \in \mathcal{E} \\ w_{ij} & \text{if } j \in \mathcal{I} \end{cases}
$$

where $\mathcal{E}$ and $\mathcal{I}$ are the sets of excitatory and inhibitory neuron indices respectively. These are stored as `W_components = {W_inst, W_delayed}`.

**In words:** Excitatory presynaptic connections are instantaneous. Inhibitory presynaptic connections are delayed by $\tau_{\text{syn}}$.

---

## 3. Firing Rate Computation (`compute_rates_and_unpack`, lines 121--202)

Given any state vector $S$ (whether current or delayed), the function unpacks the state and computes:

**Firing rate:**

$$
r_i(S) = \phi\!\left(x_i - c_P \sum_{k=1}^{K_P} a_{ik}\right)
$$

where $P \in \{E, I\}$ denotes the population that neuron $i$ belongs to, $c_P$ is the SFA coupling strength for population $P$ (i.e., $c_E$ or $c_I$), and $K_P$ is the number of SFA timescales for population $P$ (i.e., $K_E$ or $K_I$). Here $x_i$ is the dendritic state and $a_{ik}$ are the adaptation variables, so the argument to $\phi$ is the adaptation-adjusted dendritic state.

> **Implementation note:** When $K_P = 0$, the adaptation variables $a_{ik}$ are excluded from the state vector entirely and the sum vanishes. In practice, $K_I = 0$ is common (no SFA for inhibitory neurons), so that only excitatory neurons carry adaptation state.

**Depression vector:**

$$
b_i(S) = \begin{cases} b_{E,i} & \text{if } i \in \mathcal{E} \text{ and } n_{b_E} = 1 \\ b_{I,i} & \text{if } i \in \mathcal{I} \text{ and } n_{b_I} = 1 \\ 1 & \text{otherwise} \end{cases}
$$

> **Implementation note:** When $n_{b_I} = 0$ (common), the $b_I$ variables are excluded from the state vector and $b_i = 1$ for all $i \in \mathcal{I}$ (no synaptic depression for inhibitory neurons). Similarly, when $n_{b_E} = 0$, $b_E$ is excluded and $b_i = 1$ for all $i \in \mathcal{E}$.

> **Key observation:** When called on the *delayed* state $S(t - \tau_{\text{syn}})$, this function returns $r_j(S(t-\tau_{\text{syn}}))$ and $b_j(S(t-\tau_{\text{syn}}))$, i.e., the firing rate and depression are both evaluated at the delayed time. The firing rate at the delayed time incorporates the delayed adaptation: $r_j = \phi(x_j(t-\tau_{\text{syn}}) - c \sum_k a_{jk}(t-\tau_{\text{syn}}))$.

---

## 4. The Recurrent Input (`SRNN_reservoir_DDE.m`, lines 39--58)

The total recurrent input to neuron $i$ is built in two parts:

**Instantaneous part** (line 44):

$$
I_i^{(\text{inst})}(t) = \sum_{j \in \mathcal{E}} w_{ij}\, b_j(t)\, r_j(t)
$$

**Delayed part** (lines 49--57):

$$
I_i^{(\tau)}(t) = \sum_{j \in \mathcal{I}} w_{ij}\, b_j(t - \tau_{\text{syn}})\, r_j(t - \tau_{\text{syn}})
$$

where both $b_j$ and $r_j$ in the delayed term are evaluated from the **full delayed state** $S(t-\tau_{\text{syn}})$.

**Total recurrent input:**

$$
I_i^{\text{rec}}(t) = I_i^{(\text{inst})}(t) + I_i^{(\tau)}(t)
$$

---

## 5. What Is Delayed and What Is Not

| Quantity | Delayed? | Evidence |
|----------|----------|----------|
| $x_j$ in recurrent sum (I neurons) | **Yes** -- $x_j(t-\tau_{\text{syn}})$ | `compute_rates_and_unpack(S_delayed_k, ...)` unpacks $x$ from the delayed state (line 174) |
| $a_{jk}$ in recurrent sum (I neurons) | **Yes** -- $a_{jk}(t-\tau_{\text{syn}})$ | Delayed state includes adaptation; it affects $r_j$ via the adaptation-adjusted dendritic state (lines 180--184) |
| $b_j$ in recurrent sum (I neurons) | **Yes** -- $b_j(t-\tau_{\text{syn}})$ | `b_vec_delayed` is extracted from delayed state (lines 188--194) |
| $r_j$ in recurrent sum (I neurons) | **Yes** (derived from above) | $r_j(t-\tau_{\text{syn}}) = \phi(x_j(t-\tau_{\text{syn}}) - c_P \sum_k a_{jk}(t-\tau_{\text{syn}}))$ |
| $x_j$ in recurrent sum (E neurons) | No -- uses $x_j(t)$ | `compute_rates_and_unpack(S, ...)` on current state (line 37) |
| $\dot{a}_{ik}$ equation | **No** -- uses $r_i(t)$ | line 90: `r_current(E_indices)` |
| $\dot{b}_i$ equation | **No** -- uses $r_i(t)$, $b_i(t)$ | lines 105--106: `r_current`, `b_current` |

**Summary:** The delay applies to the **entire presynaptic state of inhibitory neurons** when computing their contribution to the recurrent sum. All local derivative computations ($\dot{a}$, $\dot{b}$) use only the current (non-delayed) state.

---

## 6. How `dde23` is Invoked (`SRNN_ESN.m`, lines 547--563)

```matlab
sol = dde23(@(t, y, Z) SRNN_reservoir_DDE(t, y, Z, t_ex, u_ex, params), ...
            obj.lags, history, [t_span(1), t_span(end)], options_dde);
```

- **`y`** = current state $S(t)$
- **`Z`** = `S_delay`, a matrix where column $k$ = $S(t - \tau_k)$ for `lags(k)`. Since `lags` is a scalar in practice (e.g., `lags = 0.03`), there is one column: $S(t - \tau_{\text{syn}})$.
- **`history`** = `obj.S` (a constant vector), meaning for $t < 0$, the history is the initial state.

---

## 7. Comparison with the Non-Delayed ODE (`SRNN_reservoir.m`)

The non-delayed version (used when `lags` is empty) integrates:

$$
\dot{x}_i(t) = \frac{1}{\tau_d}\!\left[-x_i(t) + u_i(t) + \sum_{j=1}^{N} w_{ij}\, b_j(t)\, r_j(t)\right]
$$

The only difference is:
1. The full weight matrix $W$ is used (no E/I split)
2. All quantities are evaluated at time $t$ (no delays)

---

## 8. Notation Crosswalk: Code to Math

| Code Variable | Math Symbol | Description |
|--------------|-------------|-------------|
| `params.lags` | $\tau_{\text{syn}}$ | Synaptic delay (scalar, in seconds) |
| `W_components{1}` | $W^{(\text{inst})}$ | Weight matrix, E columns only |
| `W_components{2}` | $W^{(\tau)}$ | Weight matrix, I columns only |
| `r_current` | $r_i(t)$ | Current firing rate |
| `r_delayed` | $r_j(t - \tau_{\text{syn}})$ | Delayed firing rate |
| `b_vec_current` | $b_j(t)$ | Current STD variable |
| `b_vec_delayed` | $b_j(t - \tau_{\text{syn}})$ | Delayed STD variable |
| `input_recurrent` | $I_i^{\text{rec}}(t)$ | Total recurrent input |
| `x_current` | $x_i(t)$ | Current dendritic state |
| `S_delay(:,k)` | $S(t - \tau_k)$ | Full delayed state vector |
| `params.tau_d` | $\tau_d$ | Dendritic time constant |
| `params.tau_a_E` | $\tau_{a_k}$ (for E) | SFA time constants (vector) |
| `params.tau_a_I` | $\tau_{a_k}$ (for I) | SFA time constants (vector) |
| `params.tau_b_E_rec` | $\tau_{\text{rec}}$ (for E) | STD recovery time constant |
| `params.tau_b_E_rel` | $\tau_{\text{rel}}$ (for E) | STD release time constant |
| `params.tau_b_I_rec` | $\tau_{\text{rec}}$ (for I) | STD recovery time constant |
| `params.tau_b_I_rel` | $\tau_{\text{rel}}$ (for I) | STD release time constant |
| `params.c_E` | $c_E$ | SFA coupling strength (E) |
| `params.c_I` | $c_I$ | SFA coupling strength (I) |
| `params.n` | $N$ | Total number of neurons |
| `params.n_E` | $N_E$ | Number of excitatory neurons |
| `params.n_I` | $N_I$ | Number of inhibitory neurons |
| `params.n_a_E` | $K_E$ | Number of SFA timescales for E neurons; when 0, $a_E$ is excluded from state |
| `params.n_a_I` | $K_I$ | Number of SFA timescales for I neurons; often 0, excluding $a_I$ from state |
| `params.n_b_E` | $n_{b_E}$ | STD flag for E (0 or 1); when 0, $b_E$ is excluded from state |
| `params.n_b_I` | $n_{b_I}$ | STD flag for I (0 or 1); often 0, excluding $b_I$ from state |
| `params.activation_function` | $\phi(\cdot)$ | Activation function |

---

## 9. Discrepancies and Observations

1. **The `corrected_equations.md` describes a non-delayed system** with the recurrent sum $\sum_j w_{ij} b_j r_j$ over all neurons. The DDE version splits this sum by excitatory vs. inhibitory presynaptic identity.

2. **The user asked about delaying $x_j(t-\tau_{\text{syn}})$** specifically. The code delays the **entire** state $S(t-\tau_{\text{syn}})$, not just $x$. This means $a_{jk}$ and $b_j$ are also evaluated at the delayed time for the purpose of computing the inhibitory recurrent contribution. The delayed $a_{jk}$ feeds into $r_j$ via the adaptation-adjusted dendritic state; the delayed $b_j$ multiplies $r_j$ directly.

3. **Only inhibitory presynaptic connections carry the delay.** This is a modeling choice hardcoded in `SRNN_ESN.m` (line 138 comment: `"E connections are instant, I connections are delayed"`). The `SRNN_reservoir_DDE.m` function itself is more general -- it accepts arbitrary `W_components` and `lags` -- but `SRNN_ESN.m` specifically wires it as E-instant / I-delayed.

4. **Single delay value.** In practice `lags` is a scalar (e.g., 0.03 s = 30 ms), giving a single delay. The code structure supports multiple lags (the for-loop on line 49), but the `SRNN_ESN.m` constructor only ever creates two components: `{W_inst, W_delayed}`.

5. **The non-delayed ODE (`SRNN_reservoir.m`) has a `noise` term** (6th argument) that is absent in `SRNN_reservoir_DDE.m`.
