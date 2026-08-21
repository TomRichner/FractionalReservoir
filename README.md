# FractionalReservoir

MATLAB research code for simulating and analyzing a Spiking Rate Neural Network
(SRNN) reservoir with spike-frequency adaptation (SFA) and short-term synaptic
depression (STD).

## Primary class

**`src/SRNNCellTypePairs.m` is the primary model class for this repo.** It
generalizes the network to `C` named cell types with per-route synapses
(`synapse_config.<pre>.<post>`), and is what new work should be built on.
`src/SRNNModel2.m` remains as the two-population (E/I) sibling that much of the
existing sweep and figure code was written against. Sweeps are driven by
`src/ParamSpaceAnalysis2.m` with `psa.model_class = 'SRNNCellTypePairs'`.

## Equations

The equations used in the **current manuscript draft** are
[`docs/EquationsParametersDocs/Equations_stability_paper_v2.md`](docs/EquationsParametersDocs/Equations_stability_paper_v2.md)
(rendered as `_v2.pdf` beside it):

$$
\begin{aligned}
dx_i &= \frac{-x_i + u_i + \sum_{j=1}^{N} w_{ij}\, r_j \prod_{m=1}^{M} b_{jm}}{\tau_d}\, dt \;+\; \frac{\sigma_u}{\tau_d}\, dW_i \\[8pt]
r_i &= \phi\left( x_i - a_{0_i} - c \sum_{k=1}^{K} a_{ik} \right) \\[8pt]
\frac{da_{ik}}{dt} &= \frac{-a_{ik} + r_i}{\tau_{a_k}}, \qquad k = 1, \dots, K \\[8pt]
\frac{db_{im}}{dt} &= \frac{1-b_{im}}{\tau_{rec_m}} - \frac{b_{im}\, r_i}{\tau_{rel_m}}, \qquad m = 1, \dots, M
\end{aligned}
$$

Note the placement of `b`: the rate `r_i` is the **pre-depression** output of the
nonlinearity, and depression enters presynaptically and multiplicatively as the
product `r_j ∏_m b_{jm}` in the recurrent sum. Both SFA and STD are therefore
driven by the raw rate `r_i`.

Short-term **facilitation (STF)** is implemented in `SRNNCellTypePairs`
(`synapse_config.<pre>.<post>.stf`, with `dg/dt = (1−g)/τ_dec + (G−g)·r/τ_fac`)
but is **not used in the current paper**, so it does not appear above.

Two other renderings of the same system sit in the same directory:
`Equations_stability_paper.md` names the product `θ_j`, and `_v3.md` additionally
inlines `φ`. They are equivalent — **v2 is the version to cite and to update.**

## Getting started

Run `setup_paths` once per MATLAB session with the working directory at the repo
root; after that any script in the repo runs from any cwd. See `CLAUDE.md` for
the full architecture notes, entry points, and conventions.
