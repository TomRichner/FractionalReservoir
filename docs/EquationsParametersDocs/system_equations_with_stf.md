# SRNN system equations with short-term facilitation (STF)

Cell-type-resolved rate reservoir (types Pyr, Pvalb, Sst, Vip). This extends the base
SRNN (spike-frequency adaptation + short-term depression) in two ways: short-term
plasticity is made **per-edge** — indexed by (presynaptic type, postsynaptic type) — and a
new **short-term facilitation** variable $p$ (a dynamic release probability) is added, with
the depression and facilitation coupled as in a Tsodyks–Markram synapse.

**Indices and symbols.** Neurons $i, j \in \{1,\dots,N\}$. The cell type of neuron $i$ is
$\theta_i \in \{1,\dots,K\}$ ($K=4$). SFA timescales $\ell \in \{1,\dots,n_a\}$. Synapse-dynamics
parameters carry a $(\text{pre-type},\ \text{post-type})$ superscript, e.g.
$\tau^{\mathrm{rec}}_{\theta_j,\, q}$ is the STD recovery constant for a synapse from a type-$\theta_j$
neuron onto a type-$q$ neuron. For the synapse $j \to i$ the post-type is $q = \theta_i$. The
external input to neuron $i$ is $u_i$ (reserved for stimulus; the STF variable is $p$, not $u$).
The firing rate $r_i = \phi(\cdot) \in (0,1)$ is **dimensionless** (normalized activity, not Hz).

## Dendritic state

$$
\frac{dx_i}{dt} = \frac{-x_i + u_i + \sum_{j=1}^{N} w_{ij}\, g_{j,\,\theta_i}\, r_j}{\tau_d}
$$

## Firing rate (with spike-frequency adaptation)

$$
r_i = \phi\!\left( x_i - a_{0_i} - c_{\theta_i} \sum_{\ell=1}^{n_a} a_{i\ell} \right)
$$

## Effective synaptic gain (depression $\times$ facilitation)

Gain of the synapse from presynaptic neuron $j$ onto a postsynaptic neuron of type $q$
(equivalently, the actual release $p\,b$ normalized by the resting release $p^0$):

$$
g_{j,q} = \frac{p_{j,q}}{p^{0}_{\theta_j,\, q}}\; b_{j,q}
$$

At rest $b_{j,q}=1$ and $p_{j,q}=p^{0}_{\theta_j, q}$, so $g_{j,q}=1$ (baseline synaptic weight
preserved). With facilitation disabled ($p_{j,q}\equiv p^{0}$), $g_{j,q}=b_{j,q}$ — pure depression.

## Spike-frequency adaptation — per neuron

$$
\frac{da_{i\ell}}{dt} = \frac{-a_{i\ell} + r_i}{\tau^{a}_{\ell}}
$$

## Short-term depression — per (presynaptic neuron $j$, post-type $q$)

Resources deplete in proportion to the actual release $p\,b$, so a facilitated synapse (larger
$p$) drains faster:

$$
\frac{db_{j,q}}{dt} = \frac{1 - b_{j,q}}{\tau^{\mathrm{rec}}_{\theta_j,\, q}} - \frac{p_{j,q}\, b_{j,q}\, r_j}{\tau^{\mathrm{rel}}_{\theta_j,\, q}}
$$

$\tau^{\mathrm{rel}}$ is the vesicle-release / depletion timescale (needed because $r$ is
dimensionless); default $\tau^{\mathrm{rel}}=1$.

## Short-term facilitation — per (presynaptic neuron $j$, post-type $q$) — **new**

$$
\frac{dp_{j,q}}{dt} = \frac{p^{0}_{\theta_j,\, q} - p_{j,q}}{\tau^{f}_{\theta_j,\, q}} + \kappa_{\theta_j,\, q}\,(1 - p_{j,q})\, r_j
$$

$\tau^{f}$ is the facilitation *decay* time (relaxation of $p$ back to $p^0$; residual-Ca$^{2+}$
clearance). $\kappa$ is the facilitation *rate* $[1/\text{time}]$ — how fast activity drives $p$
toward 1 (residual-Ca$^{2+}$ accumulation). $(1-p)$ makes the buildup saturate at the ceiling $p=1$.

## Notes

- **Depression–facilitation coupling.** Because release $= p\,b$, raising $p$ accelerates
  depletion — the Tsodyks–Markram coupling. This makes the pair a consistent rate-TM synapse and
  lets **the earlier `ml_depression_amount`→$\tau^{\mathrm{rel}}$ heuristic be dropped**: depression
  strength now flows from the release probability $p^0$ and $\tau^{\mathrm{rec}}$.
- **Two independent drive rates.** $1/\tau^{\mathrm{rel}}$ (vesicle release / pool depletion) and
  $\kappa$ (residual-Ca$^{2+}$ accumulation) are *different molecular processes* and need not be equal.
  Because $r$ is normalized, $\kappa$ folds the facilitation amount and its accumulation timescale
  into one effective rate ($\kappa = f_{\mathrm{amount}}/\tau^{\mathrm{fac}}$; they appear only as
  this ratio). $\kappa$ replaces the Tsodyks–Markram $f$ to avoid clashing with the base E/I
  model's $f$ (fraction excitatory).
- **Post-type-structured drive.** The synapse $j \to i$ uses the plasticity state indexed by the
  postsynaptic type $\theta_i$: $g_{j,\theta_i} = (p_{j,\theta_i}/p^{0}_{\theta_j,\theta_i})\,b_{j,\theta_i}$.
  With every parameter indexed by the full $(\theta_j, q)$ pair, excitatory synaptic dynamics vary
  with the postsynaptic type and inhibitory dynamics mainly with the presynaptic type — the
  alignment reported in Campagnola et al. 2022 — with no special casing.
- **State dimension.** $x$: $N$; $a$: $N \times n_a$; $b$ and $p$: $N \times K$ each (single-timescale
  STD/STF). $N_{\mathrm{sys}} = N\,(1 + n_a + 2K)$.
- **Reduction.** Freezing facilitation ($p \equiv p^{0}$) gives
  $\dot b = (1-b)/\tau^{\mathrm{rec}} - p^{0} b\, r/\tau^{\mathrm{rel}}$ — the base STD with depletion
  scaled by the resting release probability $p^{0}$; setting $p^{0}=1$ recovers the original
  $-b\,r/\tau^{\mathrm{rel}}$, and with $K=1$ the full base SFA+STD model. The offset $a_{0_i}$ is a
  per-neuron threshold (in the current implementation, realized through the center of $\phi$).
- **Parameter sources (Campagnola 2022).** $w_{ij}$ and connection sparsity ← connection probability
  and signed PSP amplitude; $c_{\theta}$ ← per-type adaptation index; $\tau^{\mathrm{rec}}$ ←
  $\text{ml\_depression\_tau}$; $p^{0}$ ← $\text{ml\_base\_release\_probability}$ (sets both the
  depression strength and the facilitation floor); $\tau^{f}$ ← $\text{ml\_facilitation\_tau}$;
  $\kappa$ ← $\text{ml\_facilitation\_amount}$ (as a rate); $\tau^{\mathrm{rel}}=1$ (modeling).
