# SRNN system equations with short-term facilitation (STF)

Cell-type-resolved rate reservoir (types Pyr, Pvalb, Sst, Vip). This extends the base
SRNN (spike-frequency adaptation + short-term depression) in two ways: short-term
plasticity is made **per-edge** — indexed by (presynaptic type, postsynaptic type) — and a
new **short-term facilitation** variable $p$ (a dynamic release probability) is added.

**Indices and symbols.** Neurons $i, j \in \{1,\dots,N\}$. The cell type of neuron $i$ is
$\theta_i \in \{1,\dots,K\}$ ($K=4$). SFA timescales $\ell \in \{1,\dots,n_a\}$. Synapse-dynamics
parameters carry a $(\text{pre-type},\ \text{post-type})$ superscript, e.g.
$\tau^{\mathrm{rec}}_{\theta_j,\, q}$ is the STD recovery constant for a synapse from a type-$\theta_j$
neuron onto a type-$q$ neuron. For the synapse $j \to i$ the post-type is $q = \theta_i$. The
external input to neuron $i$ is $u_i$ (reserved for stimulus; the STF variable is $p$, not $u$).

## Dendritic state

$$
\frac{dx_i}{dt} = \frac{-x_i + u_i + \sum_{j=1}^{N} w_{ij}\, g_{j,\,\theta_i}\, r_j}{\tau_d}
$$

## Firing rate (with spike-frequency adaptation)

$$
r_i = \phi\!\left( x_i - a_{0_i} - c_{\theta_i} \sum_{\ell=1}^{n_a} a_{i\ell} \right)
$$

## Effective synaptic gain (depression $\times$ facilitation)

Gain of the synapse from presynaptic neuron $j$ onto a postsynaptic neuron of type $q$:

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

$$
\frac{db_{j,q}}{dt} = \frac{1 - b_{j,q}}{\tau^{\mathrm{rec}}_{\theta_j,\, q}} - \frac{b_{j,q}\, r_j}{\tau^{\mathrm{rel}}_{\theta_j,\, q}}
$$

## Short-term facilitation — per (presynaptic neuron $j$, post-type $q$) — **new**

$$
\frac{dp_{j,q}}{dt} = \frac{p^{0}_{\theta_j,\, q} - p_{j,q}}{\tau^{f}_{\theta_j,\, q}} + f_{\theta_j,\, q}\,(1 - p_{j,q})\, r_j
$$

## Notes

- **Post-type-structured drive.** The synapse $j \to i$ uses the plasticity state indexed by the
  postsynaptic type $\theta_i$: $g_{j,\theta_i} = (p_{j,\theta_i}/p^{0}_{\theta_j,\theta_i})\,b_{j,\theta_i}$.
  Because every parameter is indexed by the full $(\theta_j, q)$ pair, excitatory synaptic dynamics
  vary with the postsynaptic type and inhibitory dynamics mainly with the presynaptic type — the
  alignment reported in Campagnola et al. 2022 — with no special casing.
- **Facilitation variable.** $p_{j,q} \in [p^{0}, 1]$ is a dynamic release probability that relaxes to
  its baseline $p^{0}_{\theta_j, q}$ with time constant $\tau^{f}_{\theta_j, q}$ and is driven upward by
  presynaptic activity $r_j$ at rate $f_{\theta_j, q}$. ($p$ replaces the Tsodyks–Markram $u$, since
  $u$ denotes the external input here.)
- **State dimension.** $x$: $N$; $a$: $N \times n_a$; $b$ and $p$: $N \times K$ each (single-timescale
  STD/STF). $N_{\mathrm{sys}} = N\,(1 + n_a + 2K)$.
- **Reduction.** With $K=1$ and $p \equiv p^{0}$, $g_j = b_j$ and the dendritic drive collapses to
  $\sum_j w_{ij}\, b_j\, r_j$ — the base SFA+STD model. The offset $a_{0_i}$ is a per-neuron threshold
  (in the current implementation it is realized through the center of the activation $\phi$).
- **Parameter sources (Campagnola 2022).** $w_{ij}$ and connection sparsity ← connection probability
  and signed PSP amplitude; $c_{\theta}$ ← per-type adaptation index; $\tau^{\mathrm{rec}}$ ←
  $\text{ml\_depression\_tau}$, $\tau^{\mathrm{rel}}$ ← (heuristic of) depression amount;
  $p^{0}, \tau^{f}, f$ ← release probability, facilitation $\tau$, facilitation amount.
