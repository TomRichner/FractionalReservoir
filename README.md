# FractionalReservoir

MATLAB research code for simulating and analysing a **spiking-rate neural network
(SRNN) reservoir** with spike-frequency adaptation (SFA), short-term synaptic
depression (STD) and optional facilitation (STF), built on excitatory/inhibitory
connectivity drawn from random matrix theory.

The question it exists to answer: **can adaptation mechanisms stabilise a
recurrent network that would otherwise be chaotic**, and how does the *structure*
of that adaptation — how many timescales, on which synaptic routes — change the
answer? Outputs are the largest Lyapunov exponent, the full Lyapunov spectrum,
firing-rate statistics, memory capacity, and parameter-sweep distributions.

## The model

$$
\begin{aligned}
dx_i &= \frac{-x_i + u_i + \sum_{j=1}^{N} w_{ij}\, \theta_j}{\tau_d}\, dt \;+\; \frac{\sigma_u}{\tau_d}\, dW_i \\[8pt]
\theta_i &= r_{i} \prod_{m=1}^{M} b_{im} \\[8pt]
r_i &= \phi\left( x_i - \frac{c}{K} \sum_{k=1}^{K} a_{ik} \right) \\[8pt]
\frac{da_{ik}}{dt} &= \frac{-a_{ik} + r_i}{\tau_{a_k}}, \qquad k = 1, \dots, K \\[8pt]
\frac{db_{im}}{dt} &= \frac{1-b_{im}}{\tau_{rec_m}} - \frac{b_{im}\, r_i}{\tau_{rel_m}}, \qquad m = 1, \dots, M
\end{aligned}
$$

Three details that are easy to get wrong:

- **$\theta_i$ is the synaptic output** — the rate *after* depression. $r_i$ is
  the **pre-depression** output of the nonlinearity, and depression enters
  presynaptically as a product inside the recurrent sum. Both SFA and STD are
  therefore driven by $r_i$, not by $\theta_i$.
- **Adaptation is normalised by $K$, depression is not.** Every $a_{ik}$ relaxes
  to the rate, so $\sum_k a_{ik} \to K r_i$ and dividing by $K$ makes $c$ the
  *total* adaptation budget, independent of how many timescales carry it.
  Depression multiplies rather than sums, so it needs no such factor.
- **Noise enters only $x$**, which keeps the diffusion constant and leaves the
  Lyapunov machinery measurable at any noise level.

Full statement, including the optional facilitation equations:
[`docs/EquationsParametersDocs/Equations_stability_paper.md`](docs/EquationsParametersDocs/Equations_stability_paper.md).

## Running it

Once per session, with the working directory at the repo root:

```matlab
setup_paths()
```

Everything the paper needs then runs from two entry points in `scripts/paper/`:

```matlab
run_dir = run_all_paper_analyses();   % all heavy compute, hours
results = make_all_paper_figures();   % every figure, minutes
```

Called with no argument, both take their settings from `paper_config()`. To do
both in one action, open and run
[`scripts/paper/reproduce_paper_run.m`](scripts/paper/reproduce_paper_run.m),
which is those two lines and a summary.

**[`scripts/paper/paper_config.m`](scripts/paper/paper_config.m) is the one file
to edit** — it names the preset (*which network*) and the run mode (*how much
compute*), and everything else follows it. It defaults to `'medium'`: about
3 hours of compute, and figures you can actually read. The paper's final run
uses `'production'`, which is a deliberate edit rather than the default.

## Layout

| path | what |
|---|---|
| `src/model/` | the three model classes, connectivity, nonlinearities, integrators, Jacobians |
| `src/presets/` | *which network* (`srnn_param_preset`) and *how much compute* (`analysis_run_config`) |
| `src/analysis/` | `ParamSpaceAnalysis2` and everything that runs a sweep |
| `src/figures/` | the `fig_*.m` manuscript figures, their helpers and replot tools |
| `src/plotting/` | drawing primitives the model classes themselves call |
| `src/util/` | cross-cutting helpers |
| `scripts/paper/` | the two entry points and their config — **start here** |
| `scripts/examples/` | exploratory scripts |
| `scripts/tests/` | verification scripts, run individually |
| `figs/paper/` | where the figures are written (gitignored, regenerable) |
| `data/` | run directories (gitignored) |
| `docs/` | equations, parameter reference, design notes |

`scripts/` holds only what a human invokes; everything called by something else
lives in `src/`, grouped by layer.

The two model classes are **duck-typed siblings, not a hierarchy**:
`SRNNModel2` (two hardwired E/I populations) and `SRNNCellTypePairs` (C named
cell types with per-route synapses, which is what the paper uses). They share no
implementation, so a change to one is not a change to the other.

`CLAUDE.md` carries the working conventions in more detail.
