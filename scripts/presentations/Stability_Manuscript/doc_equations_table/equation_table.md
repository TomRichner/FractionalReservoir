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
\multicolumn{2}{@{}l}{\hspace{0.75em}\textbf{Connection Weight parameters}} & & \\
$w_{ij}$ & Connection weights & — & $N \times N$ sparse weight matrix \\
$N$ & Network size & 300 & Total number of neurons \\
$f$ & Fraction excitatory & $\frac{1}{2}$ & Fraction of excitatory neurons, remainder are inhibitory \\
$S$ & Sparsity mask & — & Binary mask, $S_{ij} \sim \text{Bernoulli}(\alpha)$ \\
$\alpha$ & Connection probability & $\frac{1}{3}$ & $\alpha = \text{indegree}/N = 100/300$ \\
$F$ & Default scaling factor & $\frac{1}{\sqrt{N\alpha(2-\alpha)}}$ & Scaling factor which yields $R=1$ if $\tilde{\mu}_E, \tilde{\mu}_I, \tilde{\sigma}_E, \text{ and } \tilde{\sigma}_I$ are equal (see derivation below). \\
$\tilde{\mu}_E$ & Mean excitatory weight & $3 F$ & Normalized mean of non-zero E weights \\
$\tilde{\mu}_I$ & Mean inhibitory weight & $-4 F$ & Normalized mean of non-zero I weights. Since $f=\frac{1}{2}$, inhibition exceeds excitation ($-4 F$ vs. $+3 F$), creating a negative global outlier eigenvalue \\
$\tilde{\sigma}_E$ & Std dev excitatory & $F$ & Normalized std dev of non-zero E weights \\
$\tilde{\sigma}_I$ & Std dev inhibitory & $F$ & Normalized std dev of non-zero I weights \\
$g$ & Level of chaos & swept & Multiplicative scaling on $w_{ij}$ (\texttt{level\_of\_chaos}, e.g.\ swept over $[0.5, 3]$); the reported spectral radius is $R = g$ at the default normalization \\
\addlinespace[0.6em]
\multicolumn{2}{@{}l}{\hspace{0.75em}\textbf{Time Constants (s)}} & & \\
$\tau_d$ & Dendritic time constant & 0.1 & Membrane integration time constant \\
$\tau_{a_k}$ & SFA time constants & $\{0.25, 1.58, 10\}$ & Logspaced from 0.25 to 10, $K=3$ timescales \\
$\tau_{rec}$ & STD recovery & 1 & Synaptic vesicle recovery time constant \\
$\tau_{rel}$ & STD release & $\frac{1}{4}$ & Synaptic vesicle release time constant \\
\addlinespace[0.6em]
\multicolumn{2}{@{}l}{\hspace{0.75em}\textbf{Adaptation Strength}} & & \\
$c_E$ & SFA coupling & $\frac{1}{6}$ & SFA strength per timescale for E neurons (\texttt{master\_c\_E} in \texttt{run\_all\_analyses.m}) \\
\addlinespace[0.6em]
\multicolumn{2}{@{}l}{\hspace{0.75em}\textbf{Activation Function}} & & \\
$\phi$ & Logistic sigmoid & $\phi(x) = \dfrac{1}{1 + \exp(-4(x - S_c))}$ & Logistic function with unit slope at the center, range $(0,1)$ \\
$S_c$ & Sigmoid center & 0.35 & Horizontal center; $\phi(S_c) = \frac{1}{2}$ and slope $\phi'(S_c) = 1$ \\
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

## Table 2: Adaptation Conditions

The four conditions are defined by enabling/disabling SFA and STD via the number of timescales ($K = n_{a_E}$ for SFA, $n_{b_E}$ for STD):

| Condition | $n_{a_E}$ (K) | $n_{a_I}$ | $n_{b_E}$ | $n_{b_I}$ | Description |
|-----------|---------------|-----------|-----------|-----------|-------------|
| No Adaptation | 0 | 0 | 0 | 0 | Baseline |
| SFA Only | 3 | 0 | 0 | 0 | Spike-frequency adaptation enabled |
| STD Only | 0 | 0 | 1 | 0 | Short-term depression enabled |
| SFA + STD | 3 | 0 | 1 | 0 | Both mechanisms enabled |

**Effect on parameters:**

- When $n_{a_E} = 0$: No SFA variables $a_{ik}$; the $c_E \sum_k a_{ik}$ term is zero.
- When $n_{a_E} = 3$: Three SFA timescales with $\tau_{a_k} \in \{0.25, 1.58, 10\}$ s (logspaced) and coupling $c_E = \frac{0.5}{3}$.
- When $n_{b_E} = 0$: No STD variable $b_i$; synaptic output equals $r_i$ (equivalent to $b_i = 1$).
- When $n_{b_E} = 1$: STD enabled with a single timescale, recovery constant $\tau_{rec} = 1$ s and release constant $\tau_{rel} = \frac{1}{4}$ s; the resource enters the recurrent term as the product $b_i r_i$.

Inhibitory neurons have no adaptation mechanisms ($n_{a_I} = n_{b_I} = 0$).

> **Implementation note:** When any of $n_{a_E}$, $n_{a_I}$, $n_{b_E}$, or $n_{b_I}$ is set to zero, the corresponding state variables ($a_{ik}$ or $b_i$) are excluded from the system state vector and the Jacobian matrix. This prevents spurious zero eigenvalues that would otherwise arise from including disabled adaptation dynamics.
