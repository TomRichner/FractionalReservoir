# SRNN Model Parameter Table

This document is the **parameter reference**: state variables, derived
quantities, the activation function, and the derivation of the default scaling
factor $F$.

> **Equations live in
> [`Equations_stability_paper.md`](Equations_stability_paper.md).** This file
> used to restate them, and the copy fell behind — it showed an $a_{0_i}$
> offset no model class implements, an un-normalised $c$ rather than $c/K$, a
> scalar $\tau_{rel}$ where the model indexes it per timescale, and no noise
> term. Its one unique contribution, the facilitation equations, has been moved
> into the canonical doc, where facilitation is now stated as an optional
> mechanism.
>
> **The VALUES in Table 1 are illustrative, not authoritative.** They describe an
> earlier operating point and several are now wrong for the paper's preset — it
> says $N = 300$ where the preset runs 500, $\tau_{a_k} = \{0.1, 1, 10\}$ s where
> the model uses $\{0.25, 1.581, 10\}$, $c_E = 1/12$ "per timescale" where $c$ is
> now the total adaptation budget of $0.5$, and a single shared $\tau_{rel}$
> where the model indexes it per timescale. What this table is *for* is the
> **meaning, units and role** of each parameter, and the reasoning behind the
> choices — none of which the generated tables carry.
>
> For the values a **particular preset** actually runs at, read the generated
> `scripts/presentations/Stability_Manuscript/doc_equations_table/equation_table.md`,
> rebuilt from a live model object on every `make_all_paper_figures` run. That
> file is the record of what a sweep used; this one explains what the symbols
> mean.

**Abbreviations:**
- **SRNN**: Stable Recurrent Nonlinear Network
- **SFA**: Spike frequency adaptation
- **STD**: Short-term synaptic depression
- **STF**: Short-term synaptic facilitation

## Table 1: Model Parameters

| Symbol | Name | Value | Units | Description |
|--------|------|-------|-------|-------------|
| **State Variables** |||||
| $x_i$ | Membrane potential | $(-\infty, \infty)$ | arbitrary | Dendritic potential of neuron $i$ |
| $a_{ik}$ | Adaptation variable | $[0, 1]$ | — | SFA state for neuron $i$, timescale $k$ |
| $b_{im}$ | Synaptic resource | $[0, 1]$ | — | STD variable for neuron $i$, timescale $m$ (available resources) |
| **Dependent Variables** |||||
| $r_i$ | Firing rate | $[0, 1]$ | — | Instantaneous firing rate, $r_i = \phi(\cdot)$ |
| $\left(\prod_m b_{im}\right) r_i$ | Synaptic output | $[0, 1]$ | — | Effective synaptic output with multi-timescale STD |
| $u_i$ | External input | $[0, \infty)$ | arbitrary | External drive to neuron $i$ |
| **Connection Weight parameters** |||||
| $W$ | Connection matrix | — | — | $N \times N$ sparse weight matrix |
| $N$ | Network size | 300 | — | Total number of neurons |
| $f$ | Fraction excitatory | $\frac{1}{2}$ or (0.4–0.6) | — | Fraction of excitatory neurons, remainder are inhibitory; systematically varied from 0.4 to 0.6 to produce Figure 2G |
| $S$ | Sparsity mask | — | — | Binary mask, $S_{ij} \sim \text{Bernoulli}(\alpha)$ |
| $\alpha$ | Connection probability | $\frac{1}{3}$ | — | $\alpha = \text{indegree}/N = 100/300$ |
| $F$ | Default scaling factor | $\frac{1}{\sqrt{N\alpha(2-\alpha)}}$ | — | Scaling factor which yields $R=1$ if $\tilde{\mu}_E, \tilde{\mu}_I, \tilde{\sigma}_E, \text{ and } \tilde{\sigma}_I$ are equal (see derivation below). |
| $\tilde{\mu}_E$ | Mean excitatory weight | $3 F$ | — | Normalized mean of non-zero E weights |
| $\tilde{\mu}_I$ | Mean inhibitory weight | $-4 F$ | — | Normalized mean of non-zero I weights. If $f=\frac{1}{2}$, then inhibition exceeds excitation ($-4 F$ vs. $+3 F$), creating a negative global outlier eigenvalue|
| $\tilde{\sigma}_E$ | Std dev excitatory | $F$ | — | Normalized std dev of non-zero E weights |
| $\tilde{\sigma}_I$ | Std dev inhibitory | $F$ | — | Normalized std dev of non-zero I weights |
| **Time Constants** |||||
| $\tau_d$ | Dendritic time constant | 100 | ms | Membrane integration time constant |
| $\tau_{a_k}$ | SFA time constants | $[0.1, 1, 10]$ | s | Logspaced, $K=3$ timescales |
| $\tau_{rec_m}$ | STD recovery timescales | vector, $M$ timescales | s | Synaptic vesicle recovery time constant, one per STD timescale $m$ |
| $\tau_{rel}$ | STD release | $\frac{1}{2}$ | s | Synaptic vesicle release time constant, shared across all $M$ STD timescales |
| $M$ | Number of STD timescales | — | — | Count of STD timescales, $M = n_{b_E}$ (analogous to $K = n_{a_E}$ for SFA) |
| **Adaptation Strength** |||||
| $c_E$ | SFA coupling | $\frac{1}{12}$ | — | SFA strength per timescale |
| **Activation Function** |||||
| $\phi$ | Piecewise sigmoid | $(-\infty,\infty) \to [0,1]$ | — | Hard sigmoid with rounded corners (see definition below) |
| $q_\phi$ | Linear fraction | 0.9 | — | Fraction of the range [0, 1] that is linear (slope 1) |
| $a_0$ | Sigmoid center | 0.4 | — | Horizontal shift; $\phi(a_0) = \frac{1}{2}$ |
| **Stimulus Configuration** |||||
| $n_{steps}$ | Number of periods | 3 | — | Total number of no-stim and stim periods |
| $\rho_E$ | E input density | 0.15 | — | Fraction of E neurons receiving external input |
| $\rho_I$ | I input density | 0 | — | Fraction of I neurons receiving external input |
| $A$ | Input amplitude | $\frac{1}{2}$ | — | Amplitude scale for Gaussian random step input |
| **Simulation Settings** |||||
| $f_s$ | Sampling frequency | 400 | Hz | ODE solver sampling rate |
| $T$ | Simulation interval | $[-15, 45]$ | s | Start and end time |
| **ODE Integration** |||||
| — | Solver | RK45 | — | MATLAB ode45 (Dormand-Prince) |
| — | Relative tolerance | $10^{-9}$ | — | ODE solver RelTol |
| — | Absolute tolerance | $10^{-9}$ | — | ODE solver AbsTol |
| — | Maximum step | 2.5 | ms | ODE solver MaxStep ($1/f_s$) |
| **Lyapunov Settings** |||||
| — | LLE method | — | — | Benettin rescaling shadow trace method |
| $\Delta t_{lya}$ | Rescaling period | 20 | ms | Shadow trace rescaling interval |
| $d_0$ | Perturbation norm | $10^{-3}$ | — | Initial and rescaled perturbation magnitude |
| $f_{corner}$ | LLE filter corner | 0.25 | Hz | corner frequency for lowpass filter of local LLE for plotting, 4th order bidirectional Butterworth |

---

## Definition of the Activation Function $\phi(x)$

![Hard sigmoid with rounded corners](hard_sigmoid_rounded_corners.png)

The activation function $\phi(x)$ is a **hard sigmoid with rounded corners**—a piecewise function composed of five regions that smoothly transitions between saturation at 0 and 1.

**Domain:** $(-\infty, \infty)$  
**Range:** $[0, 1]$

**Auxiliary quantities:**

Let the half-width of the linear segment be $q = \frac{q_\phi}{2}$, and define the scaling constant for the quadratic regions:

$$k = \frac{1}{2(1 - 2q)}$$

The four breakpoints dividing the five regions are:

$$x_1 = a_0 + q - 1, \quad x_2 = a_0 - q, \quad x_3 = a_0 + q, \quad x_4 = a_0 + 1 - q$$

**Piecewise definition:**

$$
\phi(x) = 
\begin{cases}
0 & \text{if } x < x_1 \quad \text{(left saturation)} \\[4pt]
k(x - x_1)^2 & \text{if } x_1 \le x < x_2 \quad \text{(left quadratic)} \\[4pt]
(x - a_0) + \tfrac{1}{2} & \text{if } x_2 \le x \le x_3 \quad \text{(linear segment)} \\[4pt]
1 - k(x - x_4)^2 & \text{if } x_3 < x \le x_4 \quad \text{(right quadratic)} \\[4pt]
1 & \text{if } x > x_4 \quad \text{(right saturation)}
\end{cases}
$$

**Parameter interpretation:**

- **$q_\phi \in [0, 1]$**: Controls the width of the central linear region. When $q_\phi = 1$, the function reduces to a standard hard sigmoid (pure piecewise linear with corners). When $q_\phi = 0$, the function becomes purely quadratic transitions with no linear segment.

- **$a_0$**: The horizontal center of the sigmoid. At $x = a_0$, the function outputs $\phi(a_0) = \frac{1}{2}$. This parameter shifts the entire function along the $x$-axis, allowing the resting point of the firing rate to be tuned.

> **Note:** The quadratic segments ensure $C^1$ continuity (smooth corners) at the transitions between the linear and saturated regions, avoiding the discontinuous derivatives of a standard hard sigmoid.

---

## Derivation of the Default Scaling Factor $F$

The scaling factor $F$ is derived from Harris (2023), Equations 16 and 18.

**Equation 16** — Sparse variance for population $k \in \{e, i\}$:

$$\sigma_{sk}^2 = \alpha(1-\alpha)\tilde{\mu}_k^2 + \alpha \tilde{\sigma}_k^2$$

**Equation 18** — Spectral radius:

$$\mathcal{R} = \sqrt{N[f \sigma_{se}^2 + (1-f)\sigma_{si}^2]}$$

**Derivation:** To find the scaling factor $F$ that yields a unit spectral radius ($\mathcal{R} = 1$) when all normalized parameters are equal, set $|\tilde{\mu}_E| = |\tilde{\mu}_I| = \tilde{\sigma}_E = \tilde{\sigma}_I = F$.

Substituting into Eq. 16:

$$\sigma_{se}^2 = \sigma_{si}^2 = \alpha(1-\alpha)F^2 + \alpha F^2 = F^2[\alpha - \alpha^2 + \alpha] = F^2 \cdot \alpha(2-\alpha)$$

Substituting into Eq. 18, noting that $f + (1-f) = 1$:

$$\mathcal{R} = \sqrt{N \cdot F^2 \cdot \alpha(2-\alpha)} = 1$$

Solving for $F$:

$$F = \frac{1}{\sqrt{N\alpha(2-\alpha)}}$$

---

> **Note on Zero Row Sum (ZRS):** Harris (2023) describes a Zero Row Sum condition (ZRS/SZRS) that controls "local" eigenvalue outliers escaping the spectral disc. In these simulations, we deliberately did not apply the ZRS condition in order to test whether adaptation mechanisms (SFA and STD) could fulfill a similar stabilizing role—effectively examining if adaptation can substitute for ZRS in constraining network dynamics.

## Adaptation conditions

Stated in the generated
`scripts/presentations/Stability_Manuscript/doc_equations_table/adaptation_conditions.md`,
which is rebuilt from the preset's own condition list on every
`make_all_paper_figures` run.

A hand-written table lived here and had gone badly out of date: it described
**four** conditions where the paper now runs **seven**, spelled them in
`SRNNModel2`'s $n_{a_E}$ / $n_{b_E}$ vocabulary while the paper runs
`SRNNCellTypePairs`, gave $\tau_{a_k} \in \{0.1, 1, 10\}$ s where the model uses
$\{0.25, 1.581, 10\}$, quoted $c_E = 1/12$ where $c$ is now the total adaptation
budget of $0.5$, and claimed inhibitory neurons carry no adaptation while the
current preset depresses all four routes including I to I. Generating it is the
only way it stays true.
