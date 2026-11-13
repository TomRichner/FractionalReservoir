# Chained multiple-timescale adaptation: the model equations

## The model at a glance

For a first reading: the model with one set of constants for the whole network,
without cell types, routes or per-neuron constants. The equations tagged (S1) to
(S10) are simplified forms of the numbered equations of
[section 3](#3-the-full-system), and each line below names the equation it
simplifies. A state variable carries the index of its neuron ($i$, or $j$ for
the sender in the recurrent sum); a constant carries none. The symbols are
defined in [section 2](#2-notation).

**Dendritic state** (simplifies (1)): each neuron integrates its input, the
synaptic outputs it receives, and noise.

$$
\tau_d\, dx_i = \left( -x_i + u_i + \sum_{j=1}^{N} w_{ij}\, \theta_j \right) dt + \sigma_u\, dW_i \tag{S1}
$$

**Firing rate** (simplifies (2)): the output of the nonlinearity, lowered by
adaptation, before the synapse acts on it.

$$
r_i = \phi\!\left( x_i - a_0 - \frac{c}{K} \sum_{k=1}^{K} a_{ik} \right) \tag{S2}
$$

**Spike-frequency adaptation** (simplifies (3)): $K$ low-pass copies of the
neuron's own rate, one per timescale.

$$
\frac{da_{ik}}{dt} = \frac{-a_{ik} + r_i}{\tau_{a_k}}, \qquad k = 1, \dots, K \tag{S3}
$$

**Synaptic output** (simplifies (4)): the rate scaled by the releasable pool.

$$
\theta_i = r_i\, b_{i1} \tag{S4}
$$

**Releasable pool** (simplifies (5)): emptied by release, refilled from the
reserve.

$$
\frac{db_{i1}}{dt} = \frac{b_{i2} \left( 1 - b_{i1} \right)}{\tau_{rec_1}} - \frac{b_{i1}\, r_i}{\tau_{rel}} \tag{S5}
$$

**Reserve** (simplifies (6)): refills slowly toward full, and loses what it
hands to the releasable pool, divided by its relative capacity $n_2$.

$$
\frac{db_{i2}}{dt} = \frac{1 - b_{i2}}{\tau_{rec_2}} - \frac{1}{n_2}\, \frac{b_{i2} \left( 1 - b_{i1} \right)}{\tau_{rec_1}} \tag{S6}
$$

(S1) to (S6) are the multiple-timescale model of the paper's network.

**Optional: release-coupled facilitation. Not used by the paper's network.** A
synapse may also carry a release probability $g_i$ (relative to rest) and a
calcium trace $h_i$.

**Release probability** (simplifies (7)): decays to its resting value of one,
and is driven toward its ceiling $G$ by the calcium trace.

$$
\frac{dg_i}{dt} = \frac{1 - g_i}{\tau_{dec}} + \frac{\left( G - g_i \right) h_i}{\tau_{fac}} \tag{S7}
$$

**Calcium trace** (simplifies (8)): a low-pass copy of the rate.

$$
\frac{dh_i}{dt} = \frac{r_i - h_i}{\tau_h} \tag{S8}
$$

**Synaptic output and releasable pool with facilitation** (simplify (4) and
(5)): $g_i$ multiplies the output and the release term, so (S9) takes the place
of (S4) and (S10) the place of (S5). (S6) is unchanged.

$$
\theta_i = r_i\, g_i\, b_{i1} \tag{S9}
$$

$$
\frac{db_{i1}}{dt} = \frac{b_{i2} \left( 1 - b_{i1} \right)}{\tau_{rec_1}} - \frac{g_i\, b_{i1}\, r_i}{\tau_{rel}} \tag{S10}
$$

**What the simplified form leaves out, and where it is.**

- **Routes.** Depression and facilitation are per route, from a presynaptic
  cell type to a postsynaptic one, and $c$, $K$ and the $\tau_{a_k}$ are per
  cell type: [section 3](#3-the-full-system).
- **Per-neuron constants.** Each neuron draws its own time constants, reserve
  capacity, ceiling and setpoint about the values written here, which are the
  nominal ones: [section 6](#6-per-neuron-parameters).
- **The single-timescale condition**, and the other forms a synapse can take
  with a variable removed: [section 4](#4-reductions).

The simplified equations are the full system with every route alike and every
spread zero. With every route alike the route superscripts drop (the shorthand
of section 3); with every spread zero ($\sigma_{\tau_a} = \sigma_{std} =
\sigma_{stf} = \sigma_{S_c} = 0$) each neuron's constant is the nominal one and
the neuron index drops from it. (S1) to (S3) and (S6) to (S8) are then (1) to
(3) and (6) to (8), (S9) and (S10) are (4) and (5), and (S4) and (S5) are (4)
and (5) on a synapse without facilitation, which is (11). Steady states are in
[section 5](#5-steady-states-at-a-constant-rate) and the values of the
constants in [section 8](#8-parameter-values).

## 1. What this file is

The section above, "The model at a glance", is the simplified statement for a
first reading. Sections 2 to 9 are the full statement, and
[section 10](#10-the-jacobian-of-the-simplified-model) derives the Jacobian of
the simplified model.

This file states the network model with multiple-timescale (MTS) adaptation in
its **chained** form: spike-frequency adaptation with $K$ timescales, short-term
depression as a releasable vesicle pool refilled from a reserve, and
release-coupled facilitation. Depression and facilitation are per route
(presynaptic cell type to postsynaptic cell type), and every time constant is
per neuron. It is written in the notation of `Manuscript5.md` (Methods, "Network
model", and the Supplemental Methods), extended only where that notation cannot
carry something.

[`../Equations_stability_paper.md`](../Equations_stability_paper.md) remains
the statement of the model that the code is checked against, and
[`../cell_type_pair_equations.md`](../cell_type_pair_equations.md) the per-route
form of the product model. This file must agree with both in substance. Where
they differ, they are right and this file is to be corrected. The Jacobian of
the full system is not restated here; section 9 says where it is in the code,
and section 10 derives it for the simplified model.

## 2. Notation

In the column "Origin", "M" marks a symbol used as in the manuscript; "M\*" a
symbol of the manuscript that is extended here (a neuron index or a route label
is added to it, or its meaning changes with the chain); and
"new" a symbol this file adds. Rates, pool fractions and the dendritic state are dimensionless; the
firing rate is normalised to $[0, 1]$.

**Indices**

| Symbol | Meaning | Origin |
|:---------------|:-----------------------------------------------------------|:---------|
| $i$, $j$ | Neuron, $1, \dots, N$. In the recurrent sum $i$ is the postsynaptic neuron and $j$ the presynaptic one. | M |
| $q(i)$ | Cell type of neuron $i$ (E or I in the paper). | M |
| $q \to p$ | Route: from presynaptic cell type $q$ to postsynaptic cell type $p$. Written as a superscript on the synaptic variables and constants. | new |
| $k$ | SFA timescale, $1, \dots, K$. | M |
| $m$ | Pool, $1, \dots, M$. In the chain $m = 1$ is the releasable pool and $m = 2$ the reserve. In the manuscript $m$ counts independent depression variables. | M\* |

**State variables**

| Symbol | Meaning | Rest value | Code | Origin |
|:-----------|:-------------------------------------------------|:------|:----------------|:------|
| $x_i$ | Dendritic state | drawn | `x` | M |
| $a_{ik}$ | SFA variable $k$ of neuron $i$ | 0 | `a` | M |
| $b_{i1}^{q \to p}$ | Releasable pool at the terminals of neuron $i$ on route $q \to p$, as a fraction of its capacity | 1 | `b`, column 1 | M\* |
| $b_{i2}^{q \to p}$ | Reserve that refills the releasable pool, as a fraction of its capacity | 1 | `b`, column 2 | M\* |
| $g_i^{q \to p}$ | Release probability relative to rest (facilitation) | 1 | `g`, column 1 | new |
| $h_i^{q \to p}$ | Calcium trace, a low-pass copy of the rate | 0 | `g`, column 2 | new |

**Derived quantities**

| Symbol | Meaning | Origin |
|:-----------|:-----------------------------------------------------------|:------|
| $r_i$ | Firing rate: the output of the nonlinearity, **before** depression and facilitation | M |
| $\theta_i^{q \to p}$ | Synaptic output of neuron $i$ on route $q \to p$: what the recurrent sum transmits | M\* |
| $r_{\mathrm{eff}}$ | Effective rate $g\, r$ seen by the pools at steady state (section 5) | new |

**Parameters.** A constant written with a neuron index is that neuron's own
value; the same symbol without the neuron index ($\tau_{rec_m}$, $\tau_{a_k}$,
$n_2$, $G$, ...) is the **nominal** value, which the configuration states and
which is the median over neurons (section 6).

| Symbol | Meaning | Units | Code | Origin |
|:-----------|:---------------------------------------------|:-----|:-------------------|:------|
| $\tau_d$ | Dendritic time constant | s | `tau_d` | M |
| $w_{ij}$ | Connection weight onto $i$ from $j$ | -- | `W` | M |
| $u_i$ | External input | -- | `input_config` | M |
| $\sigma_u$ | Input-referred noise amplitude | units of $u$ | `sigma_u_noise` | M |
| $\phi$ | Zero-centred sigmoid, range $[0, 1]$ | -- | `activation`, `S_a` | M |
| $a_{0_i}$ | Setpoint of neuron $i$ | -- | `S_c`, `sigma_S_c`; realised `S_c_vec` | M |
| $c$ | Total SFA strength (per cell type in the code) | -- | `c` | M |
| $K$ | Number of SFA timescales (per cell type in the code) | -- | `n_a` $=$ `numel(tau_a{q})` | M |
| $\tau_{a,ik}$ | SFA time constant $k$ of neuron $i$ | s | `tau_a`; realised `tau_a_matrix` | M |
| $\tau_{rec,i1}^{q \to p}$ | Refill time of the releasable pool | s | `std.tau_rec(1)` | M\* |
| $\tau_{rec,i2}^{q \to p}$ | Refill time of the reserve | s | `std.tau_rec(2)` | M\* |
| $\tau_{rel,i}^{q \to p}$ | Release time of the releasable pool. One per route: only pool 1 releases, so it carries no pool index | s | `std.tau_rel` | M\* |
| $n_{2,i}^{q \to p}$ | Capacity of the reserve over the capacity of the releasable pool | -- | `std.n2` | new |
| $\tau_{dec,i}^{q \to p}$ | Decay time of the release probability to rest | s | `stf.tau_dec` | new |
| $\tau_{fac,i}^{q \to p}$ | Build-up time of the release probability | s | `stf.tau_fac` | new |
| $G_i^{q \to p}$ | Ceiling of the release probability, $G_i \ge 1$ | -- | `stf.G` | new |
| $\tau_{h,i}^{q \to p}$ | Time constant of the calcium trace | s | `stf.tau_h` | new |
| $\sigma_{\tau_a}$ | Log-normal spread of the SFA ladders, per cell type | -- | `tau_a_spread` | M |
| $\sigma_{std}$ | Log-normal spread of the depression constants, per presynaptic cell type | -- | `std_spread` | new |
| $\sigma_{stf}$ | Log-normal spread of the facilitation constants, per presynaptic cell type | -- | `stf_spread` | new |
| $\sigma_{S_c}$ | Standard deviation of the setpoints, per cell type | -- | `sigma_S_c` | M |

Of the extended symbols, $b_{im}$ and $\theta_i$ gain the route label, and
$b_{i2}$ a new meaning (a reserve that is not released from, where the manuscript
has a second independent pool); $\tau_{rec}$ and $\tau_{rel}$ gain the neuron
index and the route label. The symbols $g$, $h$, $G$, $\tau_{dec}$, $\tau_{fac}$, $\tau_h$ and $n_2$ are
those of `Equations_stability_paper.md`. In the manuscript $g$ is also the
overall synaptic gain in the connectivity section (`level_of_chaos` in the
code); the two never meet in one equation here, and the gain is inside $w_{ij}$.

## 3. The full system

Let neuron $i$ have cell type $q = q(i)$. Its own state is the dendritic state
$x_i$ and $K$ adaptation variables. On each route $q \to p$ that leaves its cell
type, its terminals carry up to four variables: the releasable pool, the
reserve, the release probability and the calcium trace.

**Dendritic state.** Each neuron integrates the synaptic outputs it receives,
each taken on the route from the sender's cell type to its own:

$$
\tau_d\, dx_i = \left( -x_i + u_i + \sum_{j=1}^{N} w_{ij}\, \theta_j^{q(j) \to q(i)} \right) dt + \sigma_u\, dW_i \tag{1}
$$

**Firing rate.** The rate is the output of the nonlinearity, before the synapse
acts on it:

$$
r_i = \phi\!\left( x_i - a_{0_i} - \frac{c}{K} \sum_{k=1}^{K} a_{ik} \right) \tag{2}
$$

**Spike-frequency adaptation.** Each adaptation variable is a low-pass copy of
the neuron's own rate, with that neuron's own time constant:

$$
\frac{da_{ik}}{dt} = \frac{-a_{ik} + r_i}{\tau_{a,ik}}, \qquad k = 1, \dots, K \tag{3}
$$

**Synaptic output.** On route $q \to p$ the neuron transmits its rate scaled by
the release probability and by the releasable pool:

$$
\theta_i^{q \to p} = r_i\; g_i^{q \to p}\; b_{i1}^{q \to p} \tag{4}
$$

**Releasable pool.** It empties by release, at a rate proportional to the firing
rate and to the release probability, and it is refilled from the reserve in
proportion to what the reserve holds and to its own empty fraction:

$$
\frac{db_{i1}^{q \to p}}{dt} = \frac{b_{i2}^{q \to p} \left( 1 - b_{i1}^{q \to p} \right)}{\tau_{rec,i1}^{q \to p}} - \frac{g_i^{q \to p}\, b_{i1}^{q \to p}\, r_i}{\tau_{rel,i}^{q \to p}} \tag{5}
$$

**Reserve.** It refills toward full with its own, slower time constant, and it
loses exactly what it hands to the releasable pool, divided by its relative
capacity (the same vesicles are a smaller fraction of a larger pool):

$$
\frac{db_{i2}^{q \to p}}{dt} = \frac{1 - b_{i2}^{q \to p}}{\tau_{rec,i2}^{q \to p}} - \frac{1}{n_{2,i}^{q \to p}}\, \frac{b_{i2}^{q \to p} \left( 1 - b_{i1}^{q \to p} \right)}{\tau_{rec,i1}^{q \to p}} \tag{6}
$$

**Release probability.** It decays to its resting value of one and is driven
toward its ceiling by the calcium trace:

$$
\frac{dg_i^{q \to p}}{dt} = \frac{1 - g_i^{q \to p}}{\tau_{dec,i}^{q \to p}} + \frac{\left( G_i^{q \to p} - g_i^{q \to p} \right) h_i^{q \to p}}{\tau_{fac,i}^{q \to p}} \tag{7}
$$

**Calcium trace.** It follows the firing rate:

$$
\frac{dh_i^{q \to p}}{dt} = \frac{r_i - h_i^{q \to p}}{\tau_{h,i}^{q \to p}} \tag{8}
$$

Six things to read off these equations.

1. **The synaptic variables belong to the presynaptic neuron.** Equations
   (4) to (8) are written for the neuron $i$ that fires. In the recurrent sum of
   (1) the same variables appear with the index $j$ of the sender, on the route
   $q(j) \to q(i)$. A neuron of type $q$ therefore has one copy of (4) to (8)
   for each cell type it projects to, and all its synapses onto neurons of one
   type share one copy.
2. **The rate drives everything.** $r_i$ is the output of the nonlinearity
   before depression. SFA (3), release (5) and the calcium trace (8) are driven
   by $r_i$, never by $\theta_i$. Depression and facilitation act on what a
   neuron transmits, not on how fast it fires.
3. **Only the releasable pool releases.** The reserve appears in the synaptic
   output only through $b_{i1}$; (4) holds $b_{i1}$, not $b_{i1} b_{i2}$.
   Likewise the calcium trace acts only through $g_i$.
4. **Facilitation is coupled to release.** $g_i$ multiplies the output in (4)
   and the release term in (5): a synapse that releases more per spike also
   uses its releasable pool faster.
5. **The SFA sum is normalised by $K$.** Every $a_{ik}$ relaxes to $r_i$
   whatever its $\tau_{a,ik}$, so the steady-state subtraction in (2) is
   $c\, r_i$ for any $K$ and for any spread of the time constants. $c$ is the
   total adaptation strength; $K$ and the $\tau_{a,ik}$ set only its time
   course.
6. **Noise enters only $x_i$**, with a constant coefficient. The integrator
   multiplies each Wiener increment by $\sigma_u / \tau_d$.

**Ranges and starting values.** For $r_i \in [0, 1]$ the variables stay in
$b_{i1}, b_{i2} \in [0, 1]$, $g_i \in [1, G_i]$ and $h_i \in [0, 1]$. A run
starts with full pools ($b_{i1} = b_{i2} = 1$), release probability at rest
($g_i = 1$), no calcium ($h_i = 0$), no adaptation ($a_{ik} = 0$) and a random
$x_i$.

**Shorthand when a neuron's routes are alike.** If all routes leaving a cell
type carry the same configuration, a neuron's variables on those routes obey
the same equation with the same constants from the same starting value and are
identical. The route superscript can then be dropped: $b_{im}$, $g_i$, $h_i$,
$\theta_i$, $\tau_{rec,im}$, $\tau_{rel,i}$, $n_{2,i}$, and (1) reads
$\sum_j w_{ij}\, \theta_j$ as in the manuscript. This is the case in the paper's
network (section 8), and sections 4 to 7 use the shorthand.

**Nonlinearity and setpoint.** $\phi$ is the manuscript's zero-centred piecewise
sigmoid and $a_{0_i}$ the setpoint. The code holds the setpoint inside the
function (`piecewiseSigmoid(x, S_a, S_c)`, centred at `S_c`); the two are the
same function. $c$ and $K$ are one value per cell type in the code, $c_{q(i)}$
and $K_{q(i)}$; the paper uses the same values for both types.

**Route weight.** A route may also carry a weight multiplier (`scale`, default
one, and one on every route of the paper's network), which multiplies $w_{ij}$
for the pairs $(i, j)$ of that route.

## 4. Reductions

A route carries only the variables it has. An absent variable is **removed**
from the state and from the equations; it is never kept and clamped, because a
clamped variable would enlarge the state and add a zero eigenvalue and a zero
Lyapunov exponent that mean nothing.

**No reserve: the single pool.** Drop $b_{i2}$ and equation (6), and refill the
releasable pool as if the reserve were full:

$$
\frac{db_{i1}}{dt} = \frac{1 - b_{i1}}{\tau_{rec,i1}} - \frac{g_i\, b_{i1}\, r_i}{\tau_{rel,i}}, \qquad \theta_i = r_i\, g_i\, b_{i1} \tag{9}
$$

Without facilitation ($g_i$ absent) this is the manuscript's depression equation
with $M = 1$: the single-timescale condition. It is also the limit
$n_{2,i} \to \infty$ of the chain. **The single-timescale and the chained
conditions therefore share the releasable-pool equation and differ by one
mechanism, the reserve.** On one network seed they also share, neuron by neuron,
the constants $\tau_{rec,i1}$ and $\tau_{rel,i}$ (section 6). Since the reserve
starts full, the two conditions coincide at the onset of activity and separate
as sustained firing draws the reserve down.

**No calcium trace: rate-driven release probability.** Drop $h_i$ and equation
(8); the rate takes its place in (7):

$$
\frac{dg_i}{dt} = \frac{1 - g_i}{\tau_{dec,i}} + \frac{\left( G_i - g_i \right) r_i}{\tau_{fac,i}} \tag{10}
$$

**No facilitation: the chain alone.** Drop $g_i$, $h_i$ and equations (7) and
(8), and put one in place of $g_i$ in (4) and (5):

$$
\frac{db_{i1}}{dt} = \frac{b_{i2} \left( 1 - b_{i1} \right)}{\tau_{rec,i1}} - \frac{b_{i1}\, r_i}{\tau_{rel,i}}, \qquad \theta_i = r_i\, b_{i1} \tag{11}
$$

with the reserve equation (6) unchanged. **This is the multiple-timescale
depression of the paper's network.**

**No depression.** Drop both pools; the route is a facilitating synapse,
$\theta_i = r_i\, g_i$, with (7) and (8) or with (10). With neither depression
nor facilitation $\theta_i = r_i$.

**One SFA timescale.** $K = 1$ leaves one adaptation variable and
$r_i = \phi(x_i - a_{0_i} - c\, a_{i1})$. With $K = 0$ the subtraction is
omitted.

Depression (none, single pool, chain) and facilitation (none, $g$, $g$ with $h$)
combine freely, which gives the nine variants of a route; the paper uses three
of them (none; single pool; chain), all without facilitation.

## 5. Steady states at a constant rate

Hold the rate of one neuron at $r$. The index $i$ is dropped in this section;
every constant is that neuron's own.

**Adaptation.** $a_k = r$ for every $k$, so the subtraction in (2) is $c\, r$.

**Facilitation.** $h = r$, and with or without the calcium trace

$$
g = \frac{1 + G \kappa\, r}{1 + \kappa\, r}, \qquad \kappa = \frac{\tau_{dec}}{\tau_{fac}} \tag{12}
$$

The pools then see the **effective rate** $r_{\mathrm{eff}} = g\, r$ (without
facilitation, $r_{\mathrm{eff}} = r$).

**Pools.** With

$$
k_1 = \frac{r_{\mathrm{eff}}\, \tau_{rec,1}}{\tau_{rel}}, \qquad k_2 = \frac{r_{\mathrm{eff}}\, \tau_{rec,2}}{n_2\, \tau_{rel}}, \qquad k_0 = 1 + k_1 + k_2
$$

equation (6) gives $b_2 = 1 - k_2 b_1$ and (5) then gives
$k_2 b_1^2 - k_0 b_1 + 1 = 0$, whose root in $[0, 1]$ is

$$
b_1 = \frac{2}{k_0 + \sqrt{k_0^2 - 4 k_2}}, \qquad b_2 = 1 - k_2\, b_1 \tag{13}
$$

Without a reserve $k_2 = 0$ and $b_1 = 1 / (1 + k_1)$, the single pool.

**Sustained synaptic output.** $\theta = r\, g\, b_1 = r_{\mathrm{eff}}\, b_1(r_{\mathrm{eff}})$.
It rises monotonically with the rate toward a ceiling set by how fast the two
pools can be refilled,

$$
\theta_{\max} = \frac{\tau_{rel}}{\tau_{rec,1} + \tau_{rec,2} / n_2} \tag{14}
$$

(for the single pool, $\tau_{rel} / \tau_{rec,1}$), and the steady-state reserve
never falls below $\tau_{rec,1} / (\tau_{rec,1} + \tau_{rec,2} / n_2)$. The
product of two independent pools behaves differently: its sustained output
rises, peaks and then falls with the rate (for the constants below, a peak of
0.031 at $r = 0.125$).

**What release-coupled facilitation does at steady state.** It replaces $r$ by
$g\, r$ in (13) and nothing else: it stretches the rate axis. The sustained
output approaches the ceiling (14) at lower rates, and the ceiling does not
move. The calcium trace changes only the time course: it delays $g$.

**Values** at the nominal constants of the paper's network
($\tau_{rec} = [2, 5]$ s, $\tau_{rel} = 0.25$ s, $n_2 = 2$; single pool
$\tau_{rec} = 2$ s, $\tau_{rel} = 0.25$ s; product of two pools
$\tau_{rec} = [2, 4]$ s, $\tau_{rel} = [0.25, 0.5]$ s). The depression factor is
what multiplies the rate in $\theta$.

| Rate $r$ | Single pool $b_1$ | Single pool $\theta$ | Chain $b_1$ | Chain $b_2$ | Chain $\theta$ | Product $b_1 b_2$ | Product $\theta$ |
|-----:|------:|------:|------:|------:|------:|------:|------:|
| 0.05 | 0.714 | 0.036 | 0.631 | 0.684 | 0.032 | 0.510 | 0.026 |
| 0.10 | 0.556 | 0.056 | 0.420 | 0.580 | 0.042 | 0.309 | 0.031 |
| 0.25 | 0.333 | 0.083 | 0.200 | 0.500 | 0.050 | 0.111 | 0.028 |
| 0.50 | 0.200 | 0.100 | 0.106 | 0.472 | 0.053 | 0.040 | 0.020 |
| 1.00 | 0.111 | 0.111 | 0.054 | 0.458 | 0.054 | 0.012 | 0.012 |
| limit, large $r$ | 0 | 0.125 | 0 | 0.444 | 0.0556 | 0 | 0 |

The last row is the limit of the formulas as the rate grows (the model's rate
does not exceed one); the chain's ceiling is $1/18$. The single-timescale and chained conditions are
**not** matched in steady-state depth, by decision; they share the onset.

With release-coupled facilitation at the values of section 8
($G = 5$, $\tau_{dec} = 1.5$ s, $\tau_{fac} = 0.25$ s; not used in the paper):

| Rate $r$ | $g$ | Chain $\theta$ | Single pool $\theta$ |
|-----:|------:|------:|------:|
| 0.05 | 1.92 | 0.042 | 0.054 |
| 0.10 | 2.50 | 0.050 | 0.083 |
| 0.25 | 3.40 | 0.054 | 0.109 |
| 0.50 | 4.00 | 0.055 | 0.118 |
| 1.00 | 4.43 | 0.055 | 0.122 |

## 6. Per-neuron parameters

Every time constant in section 3 carries a neuron index because each neuron
draws its own, once, when the network is built. $z$ and $\xi$ below are
independent standard normal variables, one per neuron and per quantity.

**Depression and facilitation.** With $\sigma_{std}$ and $\sigma_{stf}$ the
spreads of the neuron's (presynaptic) cell type,

$$
\begin{aligned}
\tau_{rec,i1}^{q \to p} &= \tau_{rec_1}^{q \to p}\, e^{\sigma_{std}\, z_i^{rec1}}, &
\tau_{rel,i}^{q \to p} &= \tau_{rel}^{q \to p}\, e^{\sigma_{std}\, z_i^{rel}}, \\
\tau_{rec,i2}^{q \to p} &= \tau_{rec_2}^{q \to p}\, e^{\sigma_{std}\, z_i^{rec2}}, &
n_{2,i}^{q \to p} &= n_2^{q \to p}\, e^{\sigma_{std}\, z_i^{n}}, \\
\tau_{dec,i}^{q \to p} &= \tau_{dec}^{q \to p}\, e^{\sigma_{stf}\, z_i^{dec}}, &
\tau_{fac,i}^{q \to p} &= \tau_{fac}^{q \to p}\, e^{\sigma_{stf}\, z_i^{fac}}, \\
\tau_{h,i}^{q \to p} &= \tau_{h}^{q \to p}\, e^{\sigma_{stf}\, z_i^{h}}, &
G_i^{q \to p} &= 1 + \left( G^{q \to p} - 1 \right) e^{\sigma_{stf}\, z_i^{G}} .
\end{aligned} \tag{15}
$$

Each value is log-normal about the nominal one, which is its **median**; it is
positive, and its linear spread is proportional to its size. The ceiling is
jittered on its excess over one, so $G_i \ge 1$ for every neuron. With a spread
of zero every neuron has the nominal value exactly.

**A neuron's draws are shared by its routes and across conditions.** The normals
in (15) carry a neuron index and no route label: a neuron has one $z_i^{rec1}$,
one $z_i^{rel}$ and so on, and every route it has reads them. So a neuron's
routes with the same nominal values have the same constants (which is what
allows the shorthand of section 3). And on one network seed a neuron has the
same $\tau_{rec,i1}$ and $\tau_{rel,i}$ whether its routes carry the single
pool, the chain, or the product of two pools (for the product these are its
first component's constants), and the same normal $z_i^{rec2}$ behind its second
refill time in the chain and in the product. That keeps the conditions
comparable neuron by neuron.

**The reserve stays the slower pool.** On a chain route whose nominal reserve
refill is the slower one, each neuron keeps at least a quarter of the nominal
logarithmic gap between its two refill times. If

$$
\log \frac{\tau_{rec,i2}}{\tau_{rec,i1}} < \frac{1}{4} \log \frac{\tau_{rec_2}}{\tau_{rec_1}}
\qquad \text{then} \qquad
\tau_{rec,i2} \leftarrow \tau_{rec,i1} \left( \frac{\tau_{rec_2}}{\tau_{rec_1}} \right)^{1/4} \tag{16}
$$

Only the reserve's refill time is moved; $\tau_{rec,i1}$ never is, because it
must equal the single pool's. With $\tau_{rec} = [2, 5]$ s and
$\sigma_{std} = 0.25$ the floor is a ratio of 1.257 and the rule moves 2.6 % of
neurons. The components of a product route and the facilitation constants have
no ordering to keep.

**SFA ladders.** With $\sigma_{\tau_a}$ the spread of the neuron's cell type and
$\tau_{a_1} < \dots < \tau_{a_K}$ the nominal ladder, the two ends of the ladder
are drawn and the interior rungs follow at the nominal ladder's own position in
log space:

$$
\log \tau_{a,ik} = \log \tau_{a_k} + (1 - \omega_k)\, \sigma_{\tau_a} z_i^{fast} + \omega_k\, \sigma_{\tau_a} z_i^{slow},
\qquad \omega_k = \frac{\log \tau_{a_k} - \log \tau_{a_1}}{\log \tau_{a_K} - \log \tau_{a_1}} \tag{17}
$$

For $K = 1$, $\omega_1 = 0$: the single time constant takes the fast draw alone,
$\tau_{a,i1} = \tau_{a_1} e^{\sigma_{\tau_a} z_i^{fast}}$, so that, as for
depression, a neuron's single SFA time constant in the single-timescale
condition is the fastest rung of its ladder in the multiple-timescale condition
on the same seed (the two have the same nominal value, and the equality holds
for every neuron whose ladder ends the floor below did not move). The ordering is kept
by the same quarter-of-the-nominal-span floor, applied here **symmetrically**:
where the two end offsets bring the ends closer than a quarter of
$\log (\tau_{a_K} / \tau_{a_1})$, both are moved apart about their log-midpoint
until the gap equals the floor. On the paper's ladder (ratio 40) this moves
essentially no neuron.

**Setpoints.** The setpoint is Gaussian, not log-normal:
$a_{0_i} = S_c + \sigma_{S_c}\, \xi_i$, with $S_c$ the centre.

**Expected values.** The nominal value is the median, not the mean. For a
log-normal constant $\mathrm{E}[\tau_i] = \tau\, e^{\sigma^2 / 2}$, which is
3.2 % above the nominal at $\sigma = 0.25$; the ceiling has
$\mathrm{E}[G_i] = 1 + (G - 1)\, e^{\sigma_{stf}^2 / 2}$; an interior SFA rung
has $\mathrm{E}[\tau_{a,ik}] = \tau_{a_k} \exp\!\big( \sigma_{\tau_a}^2 [(1 - \omega_k)^2 + \omega_k^2] / 2 \big)$,
1.6 % above the nominal for the middle rung of three. For the paper's chain the
expected reserve refill is 5.17 s after rule (16) (5.16 s from the log-normal
alone). The expected setpoint is $S_c$.

**Seeds.** The four draws (setpoints, SFA ladders, depression, facilitation)
come from separate random streams derived from the network seed, so a sweep over
network seeds varies them with the weights, and turning one spread on or off
changes no other draw.

## 7. How this differs from the product form

The manuscript's present equations give each neuron $M$ **independent** pools,
each released from and each refilled toward full, and multiply them:

$$
\frac{db_{im}}{dt} = \frac{1 - b_{im}}{\tau_{rec_m}} - \frac{b_{im}\, r_i}{\tau_{rel_m}}, \quad m = 1, \dots, M, \qquad \theta_i = r_i \prod_{m=1}^{M} b_{im} .
$$

In the chain the two pools are **in series**. Three things change, and (1), (2)
and (3) do not:

- the first pool's refill is gated by the second, $(1 - b_{i1}) \to b_{i2} (1 - b_{i1})$, in (5);
- the second pool is not released from; its loss term $b_{i2}\, r_i / \tau_{rel_2}$ is replaced by the flux it hands to the first pool, as in (6):

$$
\frac{db_{i2}}{dt} = \frac{1 - b_{i2}}{\tau_{rec,i2}} - \frac{1}{n_{2,i}}\, \frac{b_{i2} \left( 1 - b_{i1} \right)}{\tau_{rec,i1}}
$$

- the synaptic output holds the first pool alone, $\theta_i = r_i\, b_{i1}$, in place of $r_i\, b_{i1} b_{i2}$.

So there is one release time in place of two, and one new constant, $n_2$. With
$M = 1$ the two forms are the same equation. The consequence for the comparison
of conditions: in the product form the second component squares the steady-state
depression of the first (at equal ratios $\tau_{rec} / \tau_{rel}$), so one
versus two timescales changed depth at every rate and from the first moment; in
the chain the second pool only limits the refill of the first, the two
conditions start identical, and the sustained output saturates at (14) where the
product's falls toward zero. The manuscript also writes $\tau_{rec_m}$ and
$\tau_{rel_m}$ without a neuron index or a route; here both are per neuron and
per route.

## 8. Parameter values

**The paper's network** is the preset of the bundle `mu7chain_config`. Its
name, one word broken here in two for the page, is
`celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_` followed by
`noStim_noise0p025_chainStd_4cond_mu7chain`. Time constants are nominal values,
the medians over neurons.

| Symbol | Meaning | Value | Units |
|:-----------|:------------------------------------------|:---------------------|:-----|
| $N$ | Network size | 500 | neurons |
| | Expected in-degree; connection probability $\alpha$ | 100; 0.2 | |
| | Fraction of neurons per cell type (E, I) | 0.5, 0.5 | |
| $\tau_d$ | Dendritic time constant | 0.1 | s |
| $\tilde{\mu}$ | Weight means (post, pre), multiples of $F$ | $[7, -7;\ 7, -7]$ | |
| $\tilde{\sigma}$ | Weight standard deviations, multiples of $F$ | 1.5 on every block | |
| | Overall synaptic gain | 1 | |
| $\phi$ | Piecewise sigmoid, $S_a$ | 0.8 | |
| $a_{0_i}$ | Setpoint: centre $S_c$, spread $\sigma_{S_c}$ | 0.2, 0.1 | |
| $c$ | Total SFA strength, E and I | 0.5 | |
| $\tau_{a_k}$ | SFA ladder, $K = 3$, E and I | 0.25, 1.581, 10 | s |
| $\tau_{a_1}$ | SFA, $K = 1$ (single-timescale condition) | 0.25 | s |
| $\sigma_{\tau_a}$ | Spread of the SFA ladders, E and I | 0.25 | |
| $\tau_{rec_1}$ | Refill of the releasable pool | 2 | s |
| $\tau_{rel}$ | Release time of the releasable pool | 0.25 | s |
| $\tau_{rec_2}$ | Refill of the reserve (chain) | 5 | s |
| $n_2$ | Reserve capacity over releasable capacity (chain) | 2 | |
| $\sigma_{std}$ | Spread of the depression constants, E and I | 0.25 | |
| $\sigma_u$ | Noise amplitude | 0.025 | units of $u$ |
| $u_i$ | External input | none | |

**Routes.** In each adapting condition all four routes (E to E, E to I, I to E,
I to I) carry the same depression, so the shorthand of section 3 applies. No
route carries facilitation, and every route weight multiplier is one.

**Conditions.** Every analysis runs four; the paper shows the first three.

| Condition | Code name | SFA $\tau_{a_k}$ (s) | Depression | $\tau_{rec}$ (s) | $\tau_{rel}$ (s) | $n_2$ |
|:------------------|:---------------|:-------------|:---------------|:------|:---------|:---|
| No adaptation | `no_adaptation` | none | none | | | |
| Single timescale | `sfa1_std1` | 0.25 | single pool, (9) | 2 | 0.25 | |
| Multiple timescale, chained | `sfa3_stdc2` | 0.25, 1.581, 10 | chain, (11) and (6) | 2, 5 | 0.25 | 2 |
| Multiple timescale, product | `sfa3_std2` | 0.25, 1.581, 10 | product of two pools | 2, 4 | 0.25, 0.5 | |

**Facilitation: not used by the paper's preset.** The values chosen for the
release-coupled form (FD-018 in `docs/ISSUES.md`), used by its tests:

| Symbol | Meaning | Value | Units |
|:-----------|:------------------------------------------|:------|:-----|
| $G$ | Ceiling of the release probability ($1/U$ at a resting release probability $U = 0.2$) | 5 | |
| $\tau_{fac}$ | Build-up time, set equal to $\tau_{rel}$ | 0.25 | s |
| $\tau_{dec}$ | Decay time | 1.5 | s |
| $\tau_h$ | Calcium trace | 0.5 | s |
| $\sigma_{stf}$ | Spread of the facilitation constants | 0.25 | |

## 9. Where it is in the code

Paths are from the repository root.

| Equation | File and function |
|:----------------------|:------------------------------------------------------------------|
| (1), (2), (3): dendritic state, rate, SFA; the recurrent sum over routes | `src/model/SRNNCellTypePairs.m`, `dynamics_fast` |
| $c / K$ | `SRNNCellTypePairs.get_params` (`params.c_eff`) |
| Noise term of (1) | `src/model/integrators/sde_fixed_step.m` |
| (4) to (11): a route's state derivative and output factor $g\, b_1$ | `src/model/routes/route_pool.m`, `rhs` and `output` |
| Single pool without facilitation, (9) with $g$ absent; the product form of section 7 | `src/model/routes/route_product.m`, `rhs` |
| (12), (13), (14): steady state | `route_pool.m`, `steady_state` |
| (15), (16): per-neuron synapse constants | `route_pool.m`, `compile`; the normals are drawn in `SRNNCellTypePairs.build_synapse_spread` (`std_draws`, `stf_draws`); the column each quantity reads is fixed in `src/model/routes/route_model.m`, "THE SPREAD" |
| (17): per-neuron SFA ladders | `SRNNCellTypePairs.build_tau_a_matrix` (`tau_a_matrix`) |
| Setpoints | `SRNNCellTypePairs.build_setpoints` (`S_c_vec`) |
| Which model a route gets | `SRNNCellTypePairs.compile_route`; registry `route_model.m` |
| Jacobian (the full system's is not restated here; section 10 derives the simplified model's) | a route's block: `route_pool.m`, `linearize`, and the table in its header; assembled in `SRNNCellTypePairs.compute_Jacobian_fast`; matrix-free product in `jacobian_times` |
| A route's factors for display | `SRNNCellTypePairs.std_route_factor` ($b_1$), `stf_route_factor` ($g$) |
| Parameter values | `src/presets/srnn_param_preset.m`, the case block named in section 8 |

**Tests that pin it** (in `scripts/tests/`):

- `test_route_models`: every route variant behind the one interface: the
  linearisation against finite differences of the equations, the steady state,
  the spread.
- `test_std_chain_model`: the chain in a network: `dynamics_fast` against a
  per-neuron loop over the equations, the steady state, the limit
  $n_2 \to \infty$, the Jacobians, the ranges.
- `test_stf_chain_model`: the same for release-coupled facilitation and the
  calcium trace.
- `test_std_spread`: the per-neuron draws, the constants shared across
  conditions, rule (16).
- `test_SRNNCellTypePairs_tau_a_spread`: the SFA ladders.
- `test_jacobian_times`: the matrix-free Jacobian product.

The tolerance-tier suites `chain` and `stf_chain` hold reference values for
networks with these routes.

**Checked on 2026-11-12.** Equations (1) to (17) were coded independently, one
neuron at a time, exactly as written here, and compared with the code: a
route's derivative and output against `route_pool` for seven variants with and
without spreads (largest difference $6 \times 10^{-17}$), the steady state
against `steady_state` ($4 \times 10^{-16}$), the drawn constants against the
compiled ones ($2 \times 10^{-16}$ relative), and the full right-hand side of
two 40-neuron networks with three SFA timescales against `dynamics_fast`
($7 \times 10^{-15}$ relative).

## 10. The Jacobian of the simplified model

This section linearizes the simplified model, (S1) to (S6), and writes its
Jacobian in block form, in the manner of the effective-connectivity equation of
the previous paper (Richner et al. 2026, its Equation 3),
$J_{\mathrm{eff}} = \frac{1}{\tau_d} \left( -I + W B\, \Phi'(h) \right)$. That
equation is the block of the Jacobian that acts on the dendritic states. Here
the adaptation and pool variables are kept as states too, so the result is the
whole matrix, with $J_{\mathrm{eff}}$ as its first block.

### 10.1 Simplifications

1. **One set of constants.** Every spread of section 6 is zero except the
   setpoints': all neurons share $\tau_d$, $c$, $K$, the $\tau_{a_k}$,
   $\tau_{rec_1}$, $\tau_{rec_2}$, $\tau_{rel}$ and $n_2$. The setpoint
   $a_{0_i}$ may still differ between neurons; that costs nothing, because it
   enters only through the operating point (18).
2. **Every route alike, depression on every route**, so the shorthand of
   section 3 applies: each neuron has one releasable pool $b_{i1}$ and one
   reserve $b_{i2}$. (The code keeps one copy of a neuron's pools per route.
   With routes alike the copies are equal at all times, and the matrix below is
   the code's Jacobian acting on states with equal copies; section 10.6.)
3. **No facilitation**, as in the paper's network.
4. **$c$ and $K$ are the same for both cell types**, as in the paper's network.
5. **The input and the noise are additive** and do not depend on the state, so
   they do not appear in the Jacobian. They act on it only through the state
   at which it is evaluated.

Cell types are not otherwise simplified away: Dale's law and the E and I
weight statistics are inside $W$, which is left general.

### 10.2 Linearization

Write the state as vectors of length $N$: $x$, the SFA variables
$a_1, \dots, a_K$, and the pools $b_1$ and $b_2$. Following the previous paper,
let

$$
h_i = x_i - a_{0_i} - \frac{c}{K} \sum_{k=1}^{K} a_{ik}, \qquad r_i = \phi(h_i) \tag{18}
$$

be the operating point of neuron $i$ on its nonlinearity and its rate, and
define the diagonal matrices

$$
\Phi' = \mathrm{diag}\big( \phi'(h_i) \big), \qquad R = \mathrm{diag}(r_i), \qquad B_1 = \mathrm{diag}(b_{i1}), \qquad B_2 = \mathrm{diag}(b_{i2}) .
$$

A small perturbation of the state changes the operating point, through it the
rate, and through the rate and the releasable pool the synaptic output
$\theta = R\, b_1$:

$$
\delta h = \delta x - \frac{c}{K} \sum_{k=1}^{K} \delta a_k, \qquad
\delta r = \Phi'\, \delta h, \qquad
\delta \theta = B_1\, \delta r + R\, \delta b_1 \tag{19}
$$

Every equation of the model depends on $x$ and on the $a_k$ only through the
rate, so $\Phi'\, \delta h$ is the one channel by which the dendritic and
adaptation states act. Substituting (19) into the perturbed (S1), (S3), (S5)
and (S6):

$$
\begin{aligned}
\tau_d\, \delta \dot{x} &= -\delta x + W \left( B_1 \Phi'\, \delta h + R\, \delta b_1 \right) \\
\tau_{a_k}\, \delta \dot{a}_k &= -\delta a_k + \Phi'\, \delta h, \qquad k = 1, \dots, K \\
\delta \dot{b}_1 &= -\left( \frac{B_2}{\tau_{rec_1}} + \frac{R}{\tau_{rel}} \right) \delta b_1 + \frac{I - B_1}{\tau_{rec_1}}\, \delta b_2 - \frac{B_1}{\tau_{rel}}\, \Phi'\, \delta h \\
\delta \dot{b}_2 &= \frac{B_2}{n_2 \tau_{rec_1}}\, \delta b_1 - \left( \frac{I}{\tau_{rec_2}} + \frac{I - B_1}{n_2 \tau_{rec_1}} \right) \delta b_2
\end{aligned} \tag{20}
$$

### 10.3 The block form

Order the state as $z = (x;\ a_1;\ \dots;\ a_K;\ b_1;\ b_2)$, with $(K + 3) N$
entries. Let $\mathbf{1}_K$ be the column of $K$ ones,
$T_a = \mathrm{diag}(1/\tau_{a_1}, \dots, 1/\tau_{a_K})$, and $\otimes$ the
Kronecker product, so that $\mathbf{1}_K^{\top} \otimes M$ is the block $M$
repeated $K$ times side by side. With the **effective connectivity**

$$
W_{\mathrm{eff}} = W B_1 \Phi', \qquad J_{\mathrm{eff}} = \frac{1}{\tau_d} \left( -I + W B_1 \Phi' \right) \tag{21}
$$

equations (20) read $\delta \dot{z} = J\, \delta z$ with

$$
J =
\begin{pmatrix}
J_{xx} & J_{xa} & J_{xb_1} & 0 \\
J_{ax} & J_{aa} & 0 & 0 \\
J_{b_1 x} & J_{b_1 a} & J_{b_1 b_1} & J_{b_1 b_2} \\
0 & 0 & J_{b_2 b_1} & J_{b_2 b_2}
\end{pmatrix} \tag{22}
$$

$$
\begin{aligned}
J_{xx} &= J_{\mathrm{eff}} = \tfrac{1}{\tau_d} \left( -I + W B_1 \Phi' \right)
  & J_{xa} &= -\tfrac{c}{K \tau_d}\; \mathbf{1}_K^{\top} \otimes W B_1 \Phi'
  & J_{xb_1} &= \tfrac{1}{\tau_d}\, W R \\
J_{ax} &= T_a \mathbf{1}_K \otimes \Phi'
  & J_{aa} &= -T_a \otimes I - \tfrac{c}{K}\; T_a \mathbf{1}_K \mathbf{1}_K^{\top} \otimes \Phi' \\
J_{b_1 x} &= -\tfrac{1}{\tau_{rel}}\, B_1 \Phi'
  & J_{b_1 a} &= \tfrac{c}{K \tau_{rel}}\; \mathbf{1}_K^{\top} \otimes B_1 \Phi' \\
J_{b_1 b_1} &= -\left( \tfrac{1}{\tau_{rec_1}} B_2 + \tfrac{1}{\tau_{rel}} R \right)
  & J_{b_1 b_2} &= \tfrac{1}{\tau_{rec_1}} (I - B_1) \\
J_{b_2 b_1} &= \tfrac{1}{n_2 \tau_{rec_1}}\, B_2
  & J_{b_2 b_2} &= -\tfrac{1}{\tau_{rec_2}} I - \tfrac{1}{n_2 \tau_{rec_1}} (I - B_1)
\end{aligned} \tag{23}
$$

Every block except those holding $W$ is diagonal or made of diagonal blocks: a
neuron's adaptation and pools are coupled to that neuron alone. The same
matrix is shorter when the channel through the rate is written once:

$$
J = D + U\, \Phi'\, V^{\top} \tag{24}
$$

$$
D =
\begin{pmatrix}
-\tfrac{1}{\tau_d} I & 0 & \tfrac{1}{\tau_d} W R & 0 \\
0 & -T_a \otimes I & 0 & 0 \\
0 & 0 & J_{b_1 b_1} & J_{b_1 b_2} \\
0 & 0 & J_{b_2 b_1} & J_{b_2 b_2}
\end{pmatrix},
\qquad
U =
\begin{pmatrix}
\tfrac{1}{\tau_d} W B_1 \\
T_a \mathbf{1}_K \otimes I \\
-\tfrac{1}{\tau_{rel}} B_1 \\
0
\end{pmatrix},
\qquad
V =
\begin{pmatrix}
I \\
-\tfrac{c}{K}\, \mathbf{1}_K \otimes I \\
0 \\
0
\end{pmatrix}
$$

$V^{\top} \delta z = \delta h$ reads the operating point out of the state,
$\Phi'$ turns it into a change of rate, and $U$ says what a change of rate does
to each variable. $D$ is what is left with the rates held fixed: the leaks, the
pools' own dynamics, and the pools' action on the dendrites.

**What the form says.**

- **The nonlinearity enters only as $\Phi'$**, the slope of each neuron's
  nonlinearity at its own operating point, and it scales **columns**: neuron
  $j$'s slope multiplies what neuron $j$ sends. The setpoint (the bias) appears
  nowhere else.
- **Depression scales the same columns by $B_1$.** As in the previous paper,
  $W_{\mathrm{eff}} = W B_1 \Phi'$: a neuron with an empty releasable pool, or
  one sitting on a flat part of its nonlinearity, is removed from the effective
  connectivity for as long as that lasts. In the chain only $B_1$ scales the
  columns; the reserve acts through $J_{b_1 b_2}$, by slowing the refill.
- **SFA does not change $W_{\mathrm{eff}}$ at a given state.** It appears as
  the same columns with the opposite sign and the weight $c/K$ ($J_{xa}$), as a
  negative feedback of each neuron on itself ($J_{aa}$), and through where it
  has moved the operating point $h$.
- **The pools feed back through $W R$**: a change in a neuron's releasable
  pool is transmitted in proportion to its rate.
- **$J$ depends on the state.** It is a matrix at each moment of a run,
  $J(z(t))$; it is one matrix only at a fixed point.

### 10.4 The conditions as special cases

Removing a variable removes its rows and columns from (22); nothing is
clamped (section 4).

**No adaptation** ($K = 0$, no pools). $h_i = x_i - a_{0_i}$, the state is $x$
alone, and

$$
J = \frac{1}{\tau_d} \left( -I + W \Phi' \right) \tag{25}
$$

**Single timescale** ($K = 1$, single pool). Delete the reserve's row and
column and put $B_2 = I$ in $J_{b_1 b_1}$, which is (9) linearized. The state
is $(x;\ a_1;\ b_1)$ and

$$
J =
\begin{pmatrix}
\tfrac{1}{\tau_d} \left( -I + W B_1 \Phi' \right) & -\tfrac{c}{\tau_d}\, W B_1 \Phi' & \tfrac{1}{\tau_d}\, W R \\
\tfrac{1}{\tau_{a_1}} \Phi' & -\tfrac{1}{\tau_{a_1}} \left( I + c\, \Phi' \right) & 0 \\
-\tfrac{1}{\tau_{rel}}\, B_1 \Phi' & \tfrac{c}{\tau_{rel}}\, B_1 \Phi' & -\left( \tfrac{1}{\tau_{rec_1}} I + \tfrac{1}{\tau_{rel}} R \right)
\end{pmatrix} \tag{26}
$$

**Chained multiple timescale** ($K = 3$, chain). Equations (22) and (23) as
written.

**Product multiple timescale** differs in the synapse only. With
$\theta = R\, b_1 b_2$ and two independent pools, $W_{\mathrm{eff}} = W B_1 B_2 \Phi'$
(both pools scale the columns); the dendrites hear both pools,
$J_{xb_1} = \tfrac{1}{\tau_d} W R B_2$ and $J_{xb_2} = \tfrac{1}{\tau_d} W R B_1$;
each pool is driven by the rate,
$J_{b_m x} = -\tfrac{1}{\tau_{rel_m}} B_m \Phi'$; and the pools do not see each
other, $J_{b_m b_m} = -\left( \tfrac{1}{\tau_{rec_m}} I + \tfrac{1}{\tau_{rel_m}} R \right)$
with $J_{b_1 b_2} = J_{b_2 b_1} = 0$. In (24) that is a different $D$ and $U$
and the same $\Phi' V^{\top}$.

### 10.5 Adaptation as a gain on each neuron's column

Because adaptation and the pools are local to a neuron, they can be eliminated
neuron by neuron, which shows what they do to the effective connectivity at
each timescale. Look for a mode $\delta z\, e^{st}$ of (20). The SFA rows give
$\delta a_k = \Phi' \delta h / (1 + s \tau_{a_k})$, so by (19)

$$
\delta r = \Phi' \left( I + c\, A(s)\, \Phi' \right)^{-1} \delta x, \qquad A(s) = \frac{1}{K} \sum_{k=1}^{K} \frac{1}{1 + s \tau_{a_k}}
$$

and the pool rows give $\delta b_{i1} = -H_i(s)\, b_{i1}\, \delta r_i / \tau_{rel}$,
where $H_i(s)$ is the first diagonal entry of the inverse of neuron $i$'s
two-by-two pool block,

$$
H_i(s) = \frac{s + d_{i2}}{(s + d_{i1})(s + d_{i2}) - \dfrac{(1 - b_{i1})\, b_{i2}}{n_2 \tau_{rec_1}^2}},
\qquad d_{i1} = \frac{b_{i2}}{\tau_{rec_1}} + \frac{r_i}{\tau_{rel}},
\quad d_{i2} = \frac{1}{\tau_{rec_2}} + \frac{1 - b_{i1}}{n_2 \tau_{rec_1}} .
$$

The dendritic row then closes on $\delta x$ alone: $s$ is an eigenvalue of $J$
when

$$
\det \Big[ (1 + s \tau_d)\, I - W\, G(s) \Big] = 0, \qquad G(s) = \mathrm{diag}\big( g_i(s) \big) \tag{27}
$$

$$
g_i(s) = \underbrace{\frac{\phi'(h_i)}{1 + c\, A(s)\, \phi'(h_i)}}_{\text{slope, reduced by SFA}} \;\; \underbrace{b_{i1} \left( 1 - \frac{r_i\, H_i(s)}{\tau_{rel}} \right)}_{\text{depression}} \tag{28}
$$

(for $s$ that is not already an eigenvalue of a neuron's own SFA or pool
block). This is the no-adaptation problem, $\det[(1 + s \tau_d) I - W \Phi'] = 0$,
with each neuron's slope replaced by a gain that depends on the timescale of
the mode:

- **Fast modes** ($|s|$ much larger than every adaptation rate): $A \to 0$,
  $H_i \to 0$, and $g_i \to b_{i1}\, \phi'(h_i)$. Adaptation is frozen and the
  effective connectivity is $W B_1 \Phi'$, the previous paper's.
- **Slow modes** ($s \to 0$): $A \to 1$, so SFA divides the slope by
  $1 + c\, \phi'$; and with the pools at their steady state for the rate $r_i$,
  the depression factor becomes the slope of the sustained-output curve of
  section 5, $d\theta / dr$, which is smaller than $b_{i1}$ and goes to zero as
  the output approaches its ceiling (14).
- **In between**, each SFA timescale switches on in turn as $|s|$ falls below
  $1/\tau_{a_k}$, and the reserve after the releasable pool. Multiple
  timescales make $g_i(s)$ a graded function of timescale where a single one
  makes it a step.
- **The conditions differ only in $g_i(s)$.** No adaptation: $g_i = \phi'(h_i)$.
  Single timescale: $K = 1$ and
  $H_i = 1 / (s + 1/\tau_{rec_1} + r_i / \tau_{rel})$. Product: the depression
  factor of (28) is replaced by
  $b_{i1} b_{i2} \big( 1 - \sum_m r_i / [\tau_{rel_m} (s + 1/\tau_{rec_m} + r_i/\tau_{rel_m})] \big)$.

Equation (27) describes the eigenvalues of $J$ frozen at one state. They are
the Lyapunov exponents only at a fixed point. Along a trajectory the exponents
come from the product of the propagators of $J(z(t))$, and the eigenvalues of
the frozen matrix neither equal nor bound them.

### 10.6 Checked against the code

**Checked on 2026-11-13** by
`scripts/explorations/check_chained_MTS_jacobian_block_form.m`, which assembles
(24) from the formulas above, with the nonlinearity and every block coded
independently of the model class, and compares it with
`SRNNCellTypePairs.compute_Jacobian_fast` at two random states per case:

- One cell type, one route, 12 neurons, where the code's state is a
  permutation of $z$ (the code orders it $a$, then $b$, then $x$, with the SFA
  variables timescale by timescale and the pools pool by pool, as here):
  element by element, for no adaptation (25), single timescale (26), chained
  MTS (22) and product MTS. Largest difference $1.5 \times 10^{-16}$ of the
  largest entry.
- Two cell types and four chain routes, 16 neurons, where the code's state
  holds two copies of each neuron's pools (128 states against 96). With $T$
  the matrix that writes each pool into both copies,
  $J_{\mathrm{code}}\, T = T J$ to $2.4 \times 10^{-16}$: the block form is the
  code's Jacobian on states with equal copies. (The code's matrix has 32 more
  eigenvalues, for differences between the copies; those obey a neuron's pool
  block alone, which is stable, and a run never excites them.)
- The reduction (27): at each of the 96 eigenvalues of that $J$, the smallest
  singular value of $(1 + s \tau_d) I - W G(s)$ is at most $2 \times 10^{-12}$
  of its largest (median $10^{-15}$), against $10^{-5}$ or more a short
  distance away from an eigenvalue.

### 10.7 Why the no-adaptation edge hardly depends on the bias

**The observation** (block C, `docs/notes/Block_C_results_2026-11-13.md`).
Without adaptation $\lambda_1$, measured with the memory task's input on, turns
positive between gain 0.33 and 0.5 at every bias $S_c$ from $-0.6$ to $0.4$,
and between 0.5 and 0.7 at $S_c = 0.6$; the pilot found the network silent at
$S_c = 0.8$. With adaptation the edge lies at higher gain, and further out the
lower the bias. The network is the paper's with the noise off: $N = 500$,
balanced mean weights, a per-neuron setpoint spread of 0.1, and the
**piecewise** sigmoid of section 8 ($S_a = 0.8$: slope 1 for $|h| \le 0.4$,
falling linearly to 0 at $|h| = 0.6$). Gain is $g$, a multiplier on $W$, and
$S_c$ is the center of the setpoints $a_{0_i}$; a lower bias makes a neuron
fire at a lower input. What follows is a reading of (25), not a new
simulation. Three numbers are quoted for the network, read off the built
model: the bulk radius of $W$ at gain 1 is $R_W = 4.80$ (the code's `R`); the
eigenvalues of $W$ at network seed 3001 have largest real part 4.90 and no
outlier beyond the bulk; the row sums of $W$ have standard deviation 4.67.

**The expectation.** "For a given gain I would have expected a certain
eigenvalue abscissa of the Jacobian, and that scaling the bias appropriately
would put that abscissa on the edge." By (25) that is exactly right when every
neuron sits at the same operating point $\bar{h}$. Then $\Phi' = \phi'(\bar{h}) I$,
the eigenvalues of $J$ are those of $W$ scaled by $g\, \phi'(\bar{h})$ and
shifted, the abscissa is $(g R_W\, \phi'(\bar{h}) - 1)/\tau_d$, and for any
$g R_W > 1$ there is a bias that puts it on zero: the one with
$\phi'(\bar{h}) = 1/(g R_W)$. The edge would then be a curve in the plane of
gain and bias that bends strongly, because the bias would have to push every
neuron further onto the flat part of its nonlinearity as the gain rises.

**How the bias actually enters.** Only through $\Phi' = \mathrm{diag}(\phi'(h_i))$
with $h_i = x_i - a_{0_i}$, one slope per neuron, scaling that neuron's column
of $W$. So, as the question says, the bias does not move the eigenvalues of
one random matrix; it changes which matrix $W \Phi'$ is. For a random $W$
whose columns are scaled by slopes that are not correlated with them, the bulk
of the spectrum of $W \Phi'$ fills a disk of radius

$$
\rho \approx g R_W \sqrt{\langle \phi'^2 \rangle}, \qquad \langle \phi'^2 \rangle = \frac{1}{N} \sum_i \phi'(h_i)^2 \tag{29}
$$

so the edge condition of the frozen matrix is $g R_W \sqrt{\langle \phi'^2 \rangle} = 1$:
what matters is the mean squared slope over the **distribution** of operating
points, not the slope at one. (On this network's $W$ at gain 1, keeping a
random 40 % or 17 % of the columns and zeroing the rest gives a largest
eigenvalue modulus of 3.29 or 2.17, against 3.02 and 1.92 from (29).) Three
consequences follow from (25) and (29).

1. **Below $g R_W = 1$ the bias cannot matter.** $\phi' \le 1$, so
   $\rho \le g R_W$: for $g < 1/R_W = 0.21$ the frozen matrix is stable
   wherever the neurons sit. The bias can only lower slopes; it cannot
   destabilize.
2. **At the measured edge only a minority of neurons are on the slope.**
   Gains 0.33 to 0.5 are $g R_W = 1.6$ to $2.4$, so by (29)
   $\langle \phi'^2 \rangle$ is between 0.40 and 0.17 there. For this
   nonlinearity, whose slope is 1 or 0 except in the rounded corners, that is
   about the fraction of neurons in the linear range: a fifth to two fifths.
   If (29) holds, the network at the edge is mostly saturated, at zero or at
   one, at every bias.
3. **The network spreads its own operating points, by an amount set by the
   gain.** $x_i = u_i + g \sum_j w_{ij} r_j$, and with balanced mean weights
   the recurrent sum has a mean near zero and a standard deviation across
   neurons of about $g R_W\, q$, with $q = \sqrt{\langle r^2 \rangle}$ the
   root-mean-square rate. With rates spread over $[0, 1]$, $q$ is 0.6 to 0.7,
   so at $g R_W \approx 1.75$ the operating points are spread by about 1.1.
   That is as wide as the whole range of bias tested ($\pm 0.6$), about twice
   the half-width of the sloped part of the nonlinearity, and ten times the
   setpoint spread of 0.1.

Put together: the bias shifts the center of a distribution of operating points
that the recurrent input has already made wider than the nonlinearity. Sliding
a wide distribution across a narrow window changes the fraction inside the
window only a little. The gain acts differently: it multiplies $W$ directly,
and although the spread it causes lowers $\langle \phi'^2 \rangle$, it does so
only in proportion to $1/g$, so $\rho$ still grows, as $\sqrt{g}$. If the
operating points are taken as Gaussian with mean $m$ and a standard deviation
$s \approx g R_W\, q$ wider than the sloped region, then
$\langle \phi'^2 \rangle \approx \Delta\, e^{-m^2/2s^2} / (\sqrt{2\pi}\, s)$
with $\Delta = \int \phi'^2\, dh$ ($0.93$ for this nonlinearity, $2/3$ for the
logistic one), and (29) puts the edge at

$$
g R_W \approx \frac{\sqrt{2\pi}\; q}{\Delta}\; e^{m^2 / 2 s^2} \tag{30}
$$

For $q = 0.65$ and $m = 0$ this is $g R_W \approx 1.75$, $g \approx 0.36$,
inside the measured bracket; and shifting the center by $m = 0.6$ at
$s \approx 1.1$ multiplies it by $e^{0.15} \approx 1.16$, less than one step of
the gain grid (a factor of 1.4 to 1.5). The strong bias dependence of the
single-operating-point picture is replaced by a weak one, which a grid this
coarse would not show. Equation (30) is an estimate of scale, not a fit: $q$
and $s$ themselves depend on the gain and the bias and were not solved for.

**Why saturation cannot do the bias's work at low bias, and can at high
bias.** The rate lies in $[0, 1]$, so the two ends are not alike. A network
pushed toward **high** rates by a low bias cannot hold all its neurons at one:
with $r_j \approx 1$ the recurrent input has a spread of $g R_W \approx 1.6$ to
$2.4$ and both signs, which returns a large share of neurons to the slope or to
zero. A network pushed toward **low** rates by a high bias can stay there:
silent neurons send nothing, $q \to 0$, the spread collapses, and every slope
is zero. That is the one way the bias can stabilize at $g R_W > 1$, and it is
where the data show an effect: the edge moves out at $S_c = 0.6$ and the
network is silent at $S_c = 0.8$.

**With adaptation.** By (27) and (28) the same $W$ is scaled by
$g_i(s)$ in place of $\phi'(h_i)$, and $g_i$ depends on the rate, which the
bias sets, in a way $\phi'$ does not. The depression factor is at most
$b_{i1}$, which falls with the rate: at a sustained rate near one it is 0.11
for the single pool and 0.054 for the chain, and the product's $b_{i1} b_{i2}$
is 0.012 (section 5); near $r = 0.05$ they are 0.71, 0.63 and 0.51. SFA lowers
the operating point by $c\, r_i$, most for the neurons that fire most, and
divides the slope of slow modes by $1 + c\, \phi'$. A low bias means high
rates, hence small column gains and an edge at higher gain; a high bias means
low rates, full pools, little adaptation, and an edge near the no-adaptation
one. That is the direction the maps show, and the order of the conditions at
the lowest bias (edge at gain 1.4 to 2 for the single timescale, 2 to 2.8 for
the chain, beyond 4 for the product) is the order of those three factors. The
equations give the direction and the order; they do not give the edge's
position, which depends on the rates and pool states the network settles into.

**Where this is only a heuristic.**

- **$\lambda_1$ is not an eigenvalue.** Equations (25) to (30) describe $J$
  frozen at a state. $\lambda_1$ is the growth rate of the product of
  propagators of $J(x(t))$ along a trajectory that is driven by the task's
  input and, past the edge, chaotic. The frozen abscissa equals $\lambda_1$
  only at a fixed point, and $J$ is non-normal, so a time-varying $J$ with a
  negative abscissa at every moment can still have $\lambda_1 > 0$, and the
  reverse.
- **Equation (29) is the circular law used outside its assumptions.** It needs
  the slopes to be independent of the columns they scale and the entries of
  $W$ to be independent with zero mean. Here $W$ is sparse with E and I means,
  and the slope of neuron $j$ depends on its inputs. For $W$ itself the means
  produce no outlier at this seed; for $W \Phi'$ the mean part cancels only
  while the E and I populations have the same slope statistics, which the
  symmetric preset makes likely and nothing guarantees.
- **The Gaussian in (30) is a convenience.** With this nonlinearity the
  rates pile up at zero and one, and the distribution of operating points is
  the solution of a self-consistent problem (the rates set the spread that
  sets the rates), which the estimate does not solve. Dynamic mean-field
  theory solves it for an infinite network with independent zero-mean
  weights and no input; this network is finite, sparse, driven, and obeys
  Dale's law.
- **The input.** The task's input moves the operating points in time. At a
  high bias it may be what keeps neurons on the slope at all, and $\lambda_1$
  with the input on need not equal $\lambda_1$ without it.

**What would test it.** None of these needs a large run.

1. **The edge condition (29) directly.** At the cells on either side of the
   no-adaptation edge, evaluate $g R_W \sqrt{\langle \phi'^2 \rangle}$ along the
   trajectory, averaged over time. The explanation predicts a value near one at
   the edge at every bias, with $\langle \phi'^2 \rangle$ nearly constant across
   bias at a fixed gain. The stored fractions of pinned and silent neurons give
   a first look without any run: one minus their sum should be roughly
   $1/(g R_W)^2$ at the edge, and should change little with bias.
2. **The spread of operating points.** The standard deviation of $h_i$ across
   neurons should grow in proportion to the gain and be about one at the edge.
   If it is much narrower than the sloped region, the explanation is wrong.
3. **A bias that is large compared with the spread.** Extend the bias well
   beyond $\pm g R_W q$ (to $\pm 2$, say). The edge should bend toward higher
   gain on the high-bias side, as it begins to at $S_c = 0.6$, and stay put on
   the low-bias side.
4. **A finer gain grid** between 0.22 and 0.7 would resolve the weak
   dependence that (30) predicts, 10 to 20 %, which the present grid cannot.
5. **The frozen spectrum against $\lambda_1$.** Compare the abscissa of
   $\frac{1}{\tau_d}(-I + g W \Phi')$ at sampled states with $\lambda_1$ in the
   same cells (the eigenvalue-heatmap stage does this at the preset). Where
   they disagree in sign, the first and last caveats above are the place to
   look.
6. **For the adapting conditions,** the same with $G(0)$ and with $B_1 \Phi'$
   in place of $\Phi'$: whether the slow gain or the frozen-adaptation gain
   tracks the edge says which timescale of (28) sets it.
