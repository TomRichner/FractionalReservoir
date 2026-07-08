# SRNN Model Parameter Table

This document describes the equations, state variables, derived quantities, and
parameters used by the `run_all_analyses.m` pipeline (sensitivity, tau
sensitivity, parameter-space, and DC Lyapunov analyses). It uses a **single STD
timescale** and the **logistic activation function**, matching the current
`SRNNModel2` defaults and the master overrides set in `run_all_analyses.m`
(`c_E = 0.5/3`).

## System Equations

$$
\frac{dx_i}{dt} = \frac{-x_i + u_i + \sum_{j=1}^{N} w_{ij}\, b_j\, r_{j}}{\tau_d}
$$

$$
r_i = \phi\left( x_i - c \sum_{k=1}^{K} a_{ik} \right)
$$

$$
\frac{da_{ik}}{dt} = \frac{-a_{ik} + r_i}{\tau_{a_k}}
$$

$$
\frac{db_i}{dt} = \frac{1-b_i}{\tau_{rec}} - \frac{b_i\, r_i}{\tau_{rel}}
$$

The sparse connection weights are drawn element-wise, scaled by the synaptic gain $g$:

$$
w_{ij} = g\; S_{ij}\left( \tilde{\sigma}_j\, \xi_{ij} + \tilde{\mu}_j \right),
\qquad \xi_{ij} \sim \mathcal{N}(0,1),
\qquad S_{ij} \sim \text{Bernoulli}(\alpha)
$$

where $(\tilde{\mu}_j, \tilde{\sigma}_j) = (\tilde{\mu}_E, \tilde{\sigma}_E)$ if the presynaptic neuron $j$ is excitatory and $(\tilde{\mu}_I, \tilde{\sigma}_I)$ if inhibitory.

---

**Abbreviations:**
- **SRNN**: Stable Recurrent Nonlinear Network
- **SFA**: Spike frequency adaptation
- **STD**: Short-term synaptic depression

## Table 1: Model Parameters

\begin{longtable}{@{}>{\centering\arraybackslash}p{0.9in} >{\raggedright\arraybackslash}p{2.0in} >{\raggedright\arraybackslash}p{2.4in} >{\raggedright\arraybackslash}p{8.85in}@{}}
\toprule
\textbf{Symbol} & \textbf{Name} & \textbf{Value} & \textbf{Description} \\
\midrule
\endhead
\addlinespace[0.6em]
\multicolumn{2}{@{}l}{\hspace{0.75em}\textbf{State Variables}} & & \\
$x_i$ & Membrane potential & $(-\infty, \infty)$ & Dendritic potential of neuron $i$ \\
$a_{ik}$ & Adaptation variable & $[0, 1]$ & SFA state for neuron $i$, timescale $k$ \\
$b_i$ & Synaptic resource & $[0, 1]$ & STD variable for neuron $i$ (available resources) \\
\addlinespace[0.6em]
\multicolumn{2}{@{}l}{\hspace{0.75em}\textbf{Dependent Variables}} & & \\
$r_i$ & Firing rate & $(0, 1)$ & Instantaneous firing rate, $r_i = \phi(\cdot)$ \\
$b_i r_i$ & Synaptic output & $(0, 1)$ & Effective synaptic output after STD \\
$u_i$ & External input & arbitrary & External drive to neuron $i$ \\
\addlinespace[0.6em]
\multicolumn{2}{@{}l}{\hspace{0.75em}\textbf{Nonlinearity}} & & \\
$\phi$ & Logistic sigmoid & $\phi(x) = \dfrac{1}{1 + \exp(-4(x - S_c))}$ & Logistic function with unit slope at the center, range $(0,1)$ \\
$S_c$ & Sigmoid center & 0.35 & Horizontal center; $\phi(S_c) = \frac{1}{2}$ and slope $\phi'(S_c) = 1$ \\
\addlinespace[0.6em]
\multicolumn{2}{@{}l}{\hspace{0.75em}\textbf{Connection Weight parameters}} & & \\
$w_{ij}$ & Connection weights & — & $N \times N$ sparse weight matrix \\
$N$ & Network size & 300 & Total number of neurons \\
$f$ & Fraction excitatory & $\frac{1}{2}$ & Fraction of excitatory neurons, remainder are inhibitory \\
$S$ & Sparsity mask & — & Binary mask, $S_{ij} \sim \text{Bernoulli}(\alpha)$ \\
$\alpha$ & Connection probability & $\frac{1}{3}$ & $\alpha = \text{indegree}/N = 100/300$ \\
$F$ & Default scaling factor & $\frac{1}{\sqrt{N\alpha(2-\alpha)}}$ & Scaling factor which yields $R=1$ if $\tilde{\mu}_E, \tilde{\mu}_I, \tilde{\sigma}_E, \text{ and } \tilde{\sigma}_I$ are equal. \\
$\tilde{\mu}_E$ & Mean excitatory weight & $3 F$ & Normalized mean of non-zero E weights \\
$\tilde{\mu}_I$ & Mean inhibitory weight & $-4 F$ & Normalized mean of non-zero I weights. \\
$\tilde{\sigma}_E$ & Std dev excitatory & $F$ & Normalized std dev of non-zero E weights \\
$\tilde{\sigma}_I$ & Std dev inhibitory & $F$ & Normalized std dev of non-zero I weights \\
$g$ & Synaptic gain & 1 & factor multiplying the connection weights $w_{ij}$ \\
\addlinespace[0.6em]
\multicolumn{2}{@{}l}{\hspace{0.75em}\textbf{Time Constants (s)}} & & \\
$\tau_d$ & Dendritic time constant & 0.1 & Membrane integration time constant \\
$\tau_{a_k}$ & SFA time constants & $\{0.25, 1.58, 10\}$ & Logspaced from 0.25 to 10, $K=3$ timescales \\
$\tau_{rec}$ & STD recovery & 1 & Synaptic vesicle recovery time constant \\
$\tau_{rel}$ & STD release & $\frac{1}{4}$ & Synaptic vesicle release time constant \\
\addlinespace[0.6em]
\multicolumn{2}{@{}l}{\hspace{0.75em}\textbf{Adaptation Strength}} & & \\
$c_E$ & SFA coupling & $\frac{1}{6}$ & SFA strength per timescale for E neurons \\
\bottomrule
\end{longtable}

---

## Definition of the Activation Function $\phi(x)$

The activation function is a **logistic sigmoid** with unit slope at its center,
centered at $S_c$:

$$
\phi(x) = \frac{1}{1 + \exp\!\left(-4\,(x - S_c)\right)}
$$

**Domain:** $(-\infty, \infty)$  
**Range:** $(0, 1)$

**Derivative:**

$$
\phi'(x) = 4\,\phi(x)\,\bigl(1 - \phi(x)\bigr)
$$

**Parameter interpretation:**

- **$S_c$**: The horizontal center of the sigmoid. At $x = S_c$, the function
  outputs $\phi(S_c) = \frac{1}{2}$ and has unit slope ($\phi'(S_c) = 1$; the
  factor 4 in the exponent sets that slope). A positive $S_c$ shifts the
  inflection to positive $x$, so at the resting operating point (small $x$) the
  network sits on the lower, low-gain part of the curve.

> **Note:** The logistic sigmoid is smooth ($C^\infty$), so its derivative is
> available in closed form for the Jacobian used by the QR Lyapunov method.

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

> The adaptation conditions (Table 2) are documented separately in `adaptation_conditions.md`.
