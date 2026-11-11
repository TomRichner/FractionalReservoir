# Depression and facilitation: model options with a vesicle-pool chain

Discussion note for Tom, 2026-11-11. Not code, not manuscript text.

**This note does not replace the project's equations file.** The one
authoritative statement of the model is
[`docs/EquationsParametersDocs/Equations_stability_paper.md`](../EquationsParametersDocs/Equations_stability_paper.md).
Nothing here is implemented except section 1. Section 2 is the decided chain
(`docs/Plan_A_chained_vessicles.md`; equations (C1) to (C4) of
`../StochasticPlasticDynamicalSystemPaper/adaptation_normalization.md`).
Sections 3 to 6 are options.

**What was checked and what was not.** No simulation was run. Steady states
were derived by hand and, for sections 2 and 4, checked by substituting numbers
back into the equations (residuals at round-off, four rates). Jacobian entries
are hand derivatives, not checked against finite differences. Statements about
time courses come from linearisation only. Sources are listed at the end, each
marked as read in the corpus (`PDF2md`) or recalled from memory.

## Conventions used in every section

- Notation is the manuscript's (`Manuscript5.md`, Methods, "Network model").
- **Per route.** Equations are written per presynaptic neuron $i$, as in the
  manuscript. In the code, every depression and facilitation variable is stored
  per route (presynaptic cell type to postsynaptic cell type), one row per
  presynaptic neuron. So $\theta_j$ in the recurrent sum is the copy on the
  route from the type of $j$ to the type of $i$. When a neuron's routes carry
  identical parameters, its copies are identical trajectories.
- $r_i$ is the firing rate **before** depression. In every option here SFA,
  depression and facilitation are driven by $r_i$, never by $\theta_i$. No
  option departs from this.
- **Notation added for facilitation** (the manuscript has none):

  | Symbol | Meaning | Range, start |
  |---|---|---|
  | $g_{il}$, $l = 1, \dots, L$ | facilitation variable: release probability relative to rest | $[1, G_l]$, starts at 1 |
  | $L$ | number of facilitation components (the equations file calls it $n_f$) | |
  | $\tau_{dec_l}$ | decay time of facilitation back to 1 | |
  | $\tau_{fac_l}$ | build-up time of facilitation per unit rate | |
  | $G_l$ | ceiling of $g_{il}$, $G_l \ge 1$ | |
  | $h_i$ (section 5 only) | residual-calcium trace, a low-pass copy of the rate | $[0, 1]$, starts at 0 |
  | $\tau_h$ (section 5 only) | time constant of $h_i$ | |
  | $\bar g_i$ (section 6 only) | mean of the $g_{il}$ | $[1, G]$ |

  $g$, $\tau_{dec}$, $\tau_{fac}$ and $G$ follow the code's names. $h$ is new;
  the sketch in Plan A calls it $c$, which collides with the SFA strength $c$.
- "1TS" and "MTS" are the single-timescale and multiple-timescale adaptation
  conditions. "States per synapse" means per presynaptic neuron per route.

---

## 1. The model as implemented today

Source: `SRNNCellTypePairs.dynamics_fast`, `get_params`,
`compile_synapse_config`, `initialize_state`.

$$
\tau_d\, dx_i = \Big( -x_i + u_i + \sum_{j=1}^{N} w_{ij}\, \theta_j \Big)\, dt + \sigma_u\, dW_i \tag{1.1}
$$

$$
r_i = \phi\Big( x_i - a_{0_i} - \frac{c}{K} \sum_{k=1}^{K} a_{ik} \Big) \tag{1.2}
$$

$$
\theta_i = r_i \prod_{m=1}^{M} b_{im} \prod_{l=1}^{L} g_{il} \tag{1.3}
$$

$$
\frac{da_{ik}}{dt} = \frac{-a_{ik} + r_i}{\tau_{a,ik}}, \qquad k = 1, \dots, K \tag{1.4}
$$

$$
\frac{db_{im}}{dt} = \frac{1 - b_{im}}{\tau_{rec_m}} - \frac{b_{im}\, r_i}{\tau_{rel_m}}, \qquad m = 1, \dots, M \tag{1.5}
$$

$$
\frac{dg_{il}}{dt} = \frac{1 - g_{il}}{\tau_{dec_l}} + \frac{(G_l - g_{il})\, r_i}{\tau_{fac_l}}, \qquad l = 1, \dots, L \tag{1.6}
$$

**Symbols.** $x_i$ dendritic potential; $u_i$ external input; $w_{ij}$ weight
from $j$ to $i$; $\theta_j$ synaptic output of $j$; $\tau_d$ dendritic time
constant; $\sigma_u$ input-referred noise amplitude; $W_i$ independent Wiener
processes; $r_i \in [0, 1]$ firing rate before depression; $\phi$ sigmoid with
range 0 to 1; $a_{0_i}$ setpoint; $c$ total SFA strength; $K$ number of SFA
variables $a_{ik}$ with time constants $\tau_{a,ik}$; $M$ number of depression
variables $b_{im}$ with recovery times $\tau_{rec_m}$ and release times
$\tau_{rel_m}$; $L$ number of facilitation variables $g_{il}$ with
$\tau_{dec_l}$, $\tau_{fac_l}$, $G_l$.

**What the code does, exactly.**

- Depression components multiply, and facilitation components multiply (1.3).
- Every ODE (1.4) to (1.6) is driven by the pre-depression rate $r_i$.
- $c/K$: the code forms `c_eff = c ./ max(1, n_a)`. Each $a_{ik} \to r_i$, so
  the steady SFA term is $c\, r_i$ for any $K$. With $K = 0$ the term is absent.
- Depression and facilitation are per route. SFA is per neuron.
- Empty products are 1: $M = 0$ gives no depression, $L = 0$ no facilitation.
- Start: $a_{ik} = 0$, $b_{im} = 1$, $g_{il} = 1$, $x_i$ small Gaussian.
- Ranges: $a_{ik} \in [0, 1]$; $b_{im} \in (0, 1]$ (at $b = 0$ the derivative is
  $1/\tau_{rec} > 0$; at $b = 1$ it is $-r/\tau_{rel} \le 0$);
  $g_{il} \in [1, G_l]$ (at $g = 1$ the derivative is
  $(G-1)\, r/\tau_{fac} \ge 0$; at $g = G$ it is $(1-G)/\tau_{dec} \le 0$). The
  code requires $G_l \ge 1$.

**Steady state at constant rate $r$.** Setting each derivative to zero:

$$
a_k = r, \qquad
b_m = \frac{1}{1 + r\, \tau_{rec_m}/\tau_{rel_m}}, \qquad
g_l = \frac{1 + G_l\, \kappa_l\, r}{1 + \kappa_l\, r}, \quad \kappa_l = \frac{\tau_{dec_l}}{\tau_{fac_l}} \tag{1.7}
$$

Check of $g_l$: with (1.7), $1 - g = -(G-1)\kappa r/(1 + \kappa r)$ and
$G - g = (G-1)/(1 + \kappa r)$, so (1.6) gives
$\frac{(G-1)\, r}{1 + \kappa r}\big( -\kappa/\tau_{dec} + 1/\tau_{fac} \big) = 0$.
Limits: $g = 1$ at $r = 0$, $g \to G$ as $r \to \infty$.

**Options in the code, all off in the manuscript's preset.**

- *Geometric-mean combine* (`std.combine = 'geomean'`): the depression factor
  in (1.3) becomes $\big(\prod_m b_{im}\big)^{1/M}$. The ODEs are unchanged.
- *Route weight scale* (`route.scale = s`): multiplies $w_{ij}$ on that route.
- *Zero floor* (`std_zero_floor`): the depression factor $D$ is replaced by
  $(D - p_{min})/(1 - p_{min})$, where
  $p_{min} = \prod_m \tau_{rel_m}/(\tau_{rec_m} + \tau_{rel_m})$ is its
  steady value at $r = 1$ (combined the same way under the geometric mean). It
  acts on depression only, not on $g$.
- *Per-neuron SFA ladders* (`tau_a_spread`): this one **is** on in the preset
  (0.25, log-normal), which is why $\tau_{a,ik}$ carries a neuron index.

**The three conditions in the manuscript's preset**
(`..._dualStd_3cond_mu7revisedMedium` in `srnn_param_preset.m`):

| Condition | $K$ | $\tau_a$ (nominal, s) | $M$ | $\tau_{rec}$ (s) | $\tau_{rel}$ (s) | $L$ |
|---|---:|---|---:|---|---|---:|
| No adaptation | 0 | none | 0 | none | none | 0 |
| 1TS (`sfa1_std1`) | 1 | 0.25 | 1 | 2 | 0.25 | 0 |
| MTS (`sfa3_std2`) | 3 | 0.25, 1.58, 10 | 2 | 2, 4 | 0.25, 0.5 | 0 |

$c = 0.5$ on both cell types. Depression is on all four routes, identical on
each. There is no facilitation.

**Code versus the manuscript's Methods.** I found no numerical difference for
the manuscript's preset. The differences are in form and in what exists:

1. The manuscript gives each neuron $M$ depression variables. The code gives
   each neuron $M$ per **route**, so two copies per neuron in a two-type
   network. With identical routes the copies are equal, so the result is the
   same.
2. The code has facilitation (1.6) and the second product in (1.3). The
   manuscript has neither.
3. The code has the geometric-mean combine, the route scale and the zero
   floor. The manuscript mentions none, and the preset uses none.
4. The code computes $dx = (\dots)/\tau_d\, dt + (\sigma_u/\tau_d)\, dW$ and
   puts the setpoint inside $\phi$ (`S_c`). Both are the manuscript's form
   rewritten; no difference.
5. The equations file writes one ceiling $G$. The code allows one per
   component, $G_l$.

**Relation to the literature.** The product form with facilitation multiplying
only the output is Varela et al. (1997): $A = A_0 F D_1 D_2$, each factor with
its own dynamics, none entering another's equation (corpus). Equation (1.6) is
the Tsodyks-Markram utilisation equation in other variables. Barak and Tsodyks
(2007, corpus) write $du/dt = (U - u)/t_f + U(1 - u)R$. Put $g = u/U$:
$dg/dt = (1 - g)/t_f + U(1/U - g)R$. That is (1.6) with $G = 1/U$,
$\tau_{dec} = t_f$ and build-up rate $U R$. What today's code does **not**
take from that model is the release term $u\, x\, R$ in the resource equation:
in (1.5), depletion does not depend on $g$.

---

## 2. Chained depression, no facilitation (the decided model)

$$
\tau_d\, dx_i = \Big( -x_i + u_i + \sum_{j=1}^{N} w_{ij}\, \theta_j \Big)\, dt + \sigma_u\, dW_i \tag{2.1}
$$

$$
r_i = \phi\Big( x_i - a_{0_i} - \frac{c}{K} \sum_{k=1}^{K} a_{ik} \Big) \tag{2.2}
$$

$$
\theta_i = r_i\, b_{i1} \tag{2.3}
$$

$$
\frac{da_{ik}}{dt} = \frac{-a_{ik} + r_i}{\tau_{a,ik}}, \qquad k = 1, \dots, K \tag{2.4}
$$

$$
\frac{db_{i1}}{dt} = \frac{b_{i2}\,(1 - b_{i1})}{\tau_{rec_1}} - \frac{b_{i1}\, r_i}{\tau_{rel}} \tag{2.5}
$$

$$
\frac{db_{i2}}{dt} = \frac{1 - b_{i2}}{\tau_{rec_2}} - \frac{1}{n_2}\, \frac{b_{i2}\,(1 - b_{i1})}{\tau_{rec_1}} \tag{2.6}
$$

Facilitation: none.

**Symbols.** $x_i$ dendritic potential; $u_i$ external input; $w_{ij}$ weight
from $j$ to $i$; $\theta_j$ synaptic output of $j$; $\tau_d$ dendritic time
constant; $\sigma_u$ noise amplitude; $W_i$ Wiener processes; $r_i \in [0, 1]$
firing rate before depression; $\phi$ sigmoid; $a_{0_i}$ setpoint; $c$ total
SFA strength; $K$ SFA variables $a_{ik}$ with time constants $\tau_{a,ik}$;
$b_{i1}$ filled fraction of the releasable pool; $b_{i2}$ filled fraction of
the reserve; $\tau_{rec_1}$ refill time of the releasable pool from a full
reserve; $\tau_{rec_2}$ refill time of the reserve; $\tau_{rel}$ release time;
$n_2$ capacity of the reserve relative to the releasable pool. Both pools
start at 1.

**1TS and MTS.**

- 1TS: the releasable pool only. No reserve variable exists. (2.5) with
  $b_{i2}$ replaced by the constant 1, which is (1.5) with $M = 1$. One state
  per synapse.
- MTS: (2.5) and (2.6). Two states per synapse.
- The shared equation is (2.5). MTS is 1TS plus one mechanism: the refill term
  is multiplied by $b_{i2}$.

**Ranges.** $b_{i1} \in [0, 1]$: at $b_1 = 0$ the derivative is
$b_2/\tau_{rec_1} \ge 0$, at $b_1 = 1$ it is $-r/\tau_{rel} \le 0$.
$b_{i2} \in [0, 1]$: at $b_2 = 0$ the derivative is $1/\tau_{rec_2} > 0$, at
$b_2 = 1$ it is $-(1 - b_1)/(n_2 \tau_{rec_1}) \le 0$.

**Steady state at constant rate $r$.** Write
$k_1 = r\, \tau_{rec_1}/\tau_{rel}$, $k_2 = r\, \tau_{rec_2}/(n_2 \tau_{rel})$,
$s = 1 + k_1 + k_2$.

- (2.5) $= 0$ gives $b_2 (1 - b_1)/\tau_{rec_1} = b_1 r/\tau_{rel}$, that is
  $b_2 (1 - b_1) = k_1 b_1$.
- (2.6) $= 0$, with that flux substituted, gives
  $(1 - b_2)/\tau_{rec_2} = b_1 r/(n_2 \tau_{rel})$, that is
  $b_2 = 1 - k_2 b_1$.
- Together: $k_2 b_1^2 - s\, b_1 + 1 = 0$. The root in $[0, 1]$ is

$$
b_1 = \frac{2}{s + \sqrt{s^2 - 4 k_2}}, \qquad b_2 = 1 - k_2\, b_1 \tag{2.7}
$$

Checks. The discriminant is non-negative:
$s^2 - 4k_2 \ge (1 + k_2)^2 - 4 k_2 = (1 - k_2)^2$. $b_2 \ge 0$ because
$s + \sqrt{\cdot} \ge 1 + k_2 + |1 - k_2| \ge 2 k_2$. Numerically, (2.7) put
back into (2.5) and (2.6) leaves residuals of order $10^{-17}$.

**Reduction to the single pool.** As $n_2 \to \infty$, $k_2 \to 0$, so
$b_2 \to 1$ and $b_1 \to 2/(2s) = 1/(1 + k_1)$, which is (1.7) with $M = 1$.
Dynamically, the drain term of (2.6) vanishes, $b_2$ stays at 1, and (2.5)
becomes (1.5).

**Low and high rates.** To first order in $r$, $b_1 \approx 1 - k_1$: the
reserve does not appear, so 1TS and MTS coincide at low rate and at the onset
of activity. For large $r$, $b_1 \approx 1/(k_1 + k_2)$, so the sustained
output $r\, b_1$ rises to the ceiling

$$
\theta_{max} = \Big( \frac{\tau_{rec_1}}{\tau_{rel}} + \frac{\tau_{rec_2}}{n_2\, \tau_{rel}} \Big)^{-1} \tag{2.8}
$$

**Chosen parameters.** $n_2 = 2$, $\tau_{rec_1} = 2$ s, $\tau_{rec_2} = 5$ s,
$\tau_{rel} = 0.25$ s. Then $k_1 = 8r$, $k_2 = 10r$, $\theta_{max} = 1/18$.

| $r$ | 1TS $b_1$ | MTS $b_1$ | MTS $b_2$ | 1TS output $r b_1$ | MTS output $r b_1$ |
|---:|---:|---:|---:|---:|---:|
| 0.10 | 0.556 | 0.420 | 0.580 | 0.056 | 0.042 |
| 0.25 | 0.333 | 0.200 | 0.500 | 0.083 | 0.050 |
| 0.50 | 0.200 | 0.106 | 0.472 | 0.100 | 0.053 |
| 1.00 | 0.111 | 0.054 | 0.458 | 0.111 | 0.054 |

By decision, 1TS and MTS are not matched at steady state.

**Linearisation, for comparison with the options.** Only $b_1$ depends on the
rate and only $b_1$ feeds the output.

| | $\partial/\partial b_1$ | $\partial/\partial b_2$ | $\partial/\partial r$ |
|---|---|---|---|
| $\dot b_1$ | $-b_2/\tau_{rec_1} - r/\tau_{rel}$ | $(1 - b_1)/\tau_{rec_1}$ | $-b_1/\tau_{rel}$ |
| $\dot b_2$ | $b_2/(n_2 \tau_{rec_1})$ | $-1/\tau_{rec_2} - (1 - b_1)/(n_2 \tau_{rec_1})$ | $0$ |
| $\theta$ | $r$ | $0$ | $b_1$ |

---

## The two choices behind the options

Every way of adding facilitation to section 2 answers two separate questions.

**Choice 1: does release probability change depletion?**

- *Output only.* $g$ multiplies $\theta$ and nothing else. This is today's
  code and Varela's form.
- *Release-coupled.* $g$ also multiplies the release term of the
  releasable-pool equation. This is the Tsodyks-Markram form.

**Choice 2: what is the second facilitation variable in MTS?**

- *None.* One $g$ in both conditions.
- *Serial.* A calcium trace drives $g$. The second variable sits upstream.
- *Parallel.* Two $g$ with different time constants, combined.

The options below are four points on that grid. The pairings not shown follow
by deleting or adding $g$ in one term.

| | Output only | Release-coupled |
|---|---|---|
| One $g$ | A (section 3) | B (section 4) |
| Serial, $h \to g$ | not shown | C (section 5) |
| Parallel, $g_1, g_2$ | not shown | D (section 6) |

Illustrative facilitation numbers in the tables below use the only
facilitation parameters in the repository, from the `single_neuron_stf` preset:
$\tau_{dec} = 6$ s, $\tau_{fac} = 2.5$ s, $G = 1/0.35 \approx 2.86$, so
$\kappa = 2.4$. They are for illustration, not proposals.

---

## 3. Option A: one facilitation variable beside the chain, output only

The minimal option. Today's facilitation, unchanged, on a chain route.

$$
\tau_d\, dx_i = \Big( -x_i + u_i + \sum_{j=1}^{N} w_{ij}\, \theta_j \Big)\, dt + \sigma_u\, dW_i \tag{3.1}
$$

$$
r_i = \phi\Big( x_i - a_{0_i} - \frac{c}{K} \sum_{k=1}^{K} a_{ik} \Big) \tag{3.2}
$$

$$
\theta_i = r_i\, g_i\, b_{i1} \tag{3.3}
$$

$$
\frac{da_{ik}}{dt} = \frac{-a_{ik} + r_i}{\tau_{a,ik}}, \qquad k = 1, \dots, K \tag{3.4}
$$

$$
\frac{db_{i1}}{dt} = \frac{b_{i2}\,(1 - b_{i1})}{\tau_{rec_1}} - \frac{b_{i1}\, r_i}{\tau_{rel}} \tag{3.5}
$$

$$
\frac{db_{i2}}{dt} = \frac{1 - b_{i2}}{\tau_{rec_2}} - \frac{1}{n_2}\, \frac{b_{i2}\,(1 - b_{i1})}{\tau_{rec_1}} \tag{3.6}
$$

$$
\frac{dg_i}{dt} = \frac{1 - g_i}{\tau_{dec}} + \frac{(G - g_i)\, r_i}{\tau_{fac}} \tag{3.7}
$$

**Symbols.** $x_i$ dendritic potential; $u_i$ external input; $w_{ij}$ weight
from $j$ to $i$; $\theta_j$ synaptic output of $j$; $\tau_d$ dendritic time
constant; $\sigma_u$ noise amplitude; $W_i$ Wiener processes; $r_i \in [0, 1]$
firing rate before depression; $\phi$ sigmoid; $a_{0_i}$ setpoint; $c$ total
SFA strength; $K$ SFA variables $a_{ik}$, time constants $\tau_{a,ik}$;
$b_{i1} \in [0,1]$ releasable pool and $b_{i2} \in [0,1]$ reserve, both
starting at 1; $\tau_{rec_1}$, $\tau_{rec_2}$ their refill times; $\tau_{rel}$
release time; $n_2$ relative reserve capacity; $g_i \in [1, G]$ facilitation,
starting at 1; $\tau_{dec}$ its decay time; $\tau_{fac}$ its build-up time;
$G \ge 1$ its ceiling.

**1TS and MTS.**

- Depression: as in section 2. 1TS has $b_{i1}$ only ((3.5) with
  $b_{i2} \equiv 1$); MTS adds (3.6). Shared equation: (3.5).
- Facilitation: **the same in both conditions**, one variable obeying (3.7).
  There is no 1TS versus MTS distinction for facilitation in this option.
- States per synapse: 2 in 1TS ($b_1$, $g$), 3 in MTS ($b_1$, $b_2$, $g$).

**Coupling.** None. $g$ multiplies the output (3.3) only. Mathematically, the
depression subsystem (3.5), (3.6) and the facilitation equation (3.7) are
independent of each other and are both driven by $r_i$. Physically, a
facilitated synapse transmits more per unit of resource and does not use the
resource faster. That is consistent only if $b$ is not a vesicle count (see
the biophysical reading below).

**Steady state.** $b_1$, $b_2$ from (2.7), $g$ from (1.7), independently:

$$
\theta_{ss}(r) = r\; \frac{1 + G \kappa r}{1 + \kappa r}\; b_1(r), \qquad \kappa = \frac{\tau_{dec}}{\tau_{fac}} \tag{3.8}
$$

Closed form. The sustained ceiling (2.8) is multiplied by up to $G$.

**Limits.** $G = 1$: $g \equiv 1$, section 2 exactly. $n_2 \to \infty$:
section 1 with $M = 1$, $L = 1$ exactly.

**Ranges.** All three invariance arguments of sections 1 and 2 hold unchanged,
because the equations are unchanged.

**Jacobian relative to section 2.**

| New or changed | Entry |
|---|---|
| $g$ depends on rate | $\partial \dot g/\partial r = (G - g)/\tau_{fac}$, giving $(g, x)$ and $(g, a)$ blocks |
| $g$ self-coupling | $\partial \dot g/\partial g = -1/\tau_{dec} - r/\tau_{fac}$ |
| $g$ feeds the output | $\partial \theta/\partial g = r\, b_1$, a new $(x, g)$ block |
| Existing output derivatives | $\partial\theta/\partial b_1 = r g$, $\partial\theta/\partial r = g\, b_1$ (factor $g$ added) |
| $(b, g)$ and $(g, b)$ | zero |

According to Plan A's table of the code, these facilitation blocks already
exist for product routes. I did not read `compute_Jacobian_fast` myself.

**Biophysical reading.** Depression and facilitation as separate
multiplicative factors is a descriptive form (Varela et al. 1997, corpus). It
is coherent with a vesicle chain only if facilitation is not an increase in
vesicle release, for example a postsynaptic effect. It is also coherent if
depression is not depletion: Campagnola et al. (2022, corpus) report that
their mouse cortical data were mostly better fitted by release-independent
depression than by vesicle depletion. But then the chain's two pools lose
their literal reading too. Caveat level: high. This option is a description,
not a mechanism.

**For this project.**

- Strengths: no new equation; the facilitation code path exists; exact
  reduction to sections 1 and 2; lifts Plan A's refusal of `stf` on chain
  routes with no new mathematics.
- Weaknesses: facilitation raises the sustained gain of a route by up to $G$
  with nothing limiting it, which matters in a paper about stability. It says
  nothing about 1TS versus MTS for facilitation. It sits awkwardly beside a
  depression model whose stated virtue is that it is a pool model.

---

## 4. Option B: one facilitation variable, release-coupled

The Tsodyks-Markram coupling placed on the two-pool chain.

$$
\tau_d\, dx_i = \Big( -x_i + u_i + \sum_{j=1}^{N} w_{ij}\, \theta_j \Big)\, dt + \sigma_u\, dW_i \tag{4.1}
$$

$$
r_i = \phi\Big( x_i - a_{0_i} - \frac{c}{K} \sum_{k=1}^{K} a_{ik} \Big) \tag{4.2}
$$

$$
\theta_i = r_i\, g_i\, b_{i1} \tag{4.3}
$$

$$
\frac{da_{ik}}{dt} = \frac{-a_{ik} + r_i}{\tau_{a,ik}}, \qquad k = 1, \dots, K \tag{4.4}
$$

$$
\frac{db_{i1}}{dt} = \frac{b_{i2}\,(1 - b_{i1})}{\tau_{rec_1}} - \frac{g_i\, b_{i1}\, r_i}{\tau_{rel}} \tag{4.5}
$$

$$
\frac{db_{i2}}{dt} = \frac{1 - b_{i2}}{\tau_{rec_2}} - \frac{1}{n_2}\, \frac{b_{i2}\,(1 - b_{i1})}{\tau_{rec_1}} \tag{4.6}
$$

$$
\frac{dg_i}{dt} = \frac{1 - g_i}{\tau_{dec}} + \frac{(G - g_i)\, r_i}{\tau_{fac}} \tag{4.7}
$$

The only difference from option A is the factor $g_i$ in the release term of
(4.5).

**Symbols.** $x_i$ dendritic potential; $u_i$ external input; $w_{ij}$ weight
from $j$ to $i$; $\theta_j$ synaptic output of $j$; $\tau_d$ dendritic time
constant; $\sigma_u$ noise amplitude; $W_i$ Wiener processes; $r_i \in [0, 1]$
firing rate before depression; $\phi$ sigmoid; $a_{0_i}$ setpoint; $c$ total
SFA strength; $K$ SFA variables $a_{ik}$, time constants $\tau_{a,ik}$;
$b_{i1} \in [0,1]$ releasable pool and $b_{i2} \in [0,1]$ reserve, both
starting at 1; $\tau_{rec_1}$, $\tau_{rec_2}$ refill times; $\tau_{rel}$
release time at resting release probability; $n_2$ relative reserve capacity;
$g_i \in [1, G]$ release probability relative to rest, starting at 1;
$\tau_{dec}$, $\tau_{fac}$, $G \ge 1$ its decay time, build-up time and
ceiling.

**1TS and MTS.**

- Depression: 1TS has $b_{i1}$ only ((4.5) with $b_{i2} \equiv 1$); MTS adds
  (4.6). Shared equation: (4.5).
- Facilitation: the same in both conditions, (4.7). No 1TS versus MTS
  distinction for facilitation.
- States per synapse: 2 in 1TS, 3 in MTS.

**Coupling.** The flux leaving the releasable pool is
$g\, b_1\, r/\tau_{rel}$, and the output (4.3) is $\tau_{rel}$ times that
flux. So the output is proportional to what is released, and what is released
is what depletes. Physically, $g/\tau_{rel}$ is the fraction of the releasable
pool used per unit of rate: a release probability. Facilitation raises it, and
the pool then empties faster. Mathematically, the structure becomes a cascade:
$r \to g$, and $(r, g) \to b_1 \leftrightarrow b_2$. $g$ is still independent
of $b$.

**Steady state.** $g$ depends on $r$ alone, so (1.7) still holds. In (4.5),
$g$ and $r$ appear only as the product $g r$. So the depression steady state
is that of section 2 evaluated at the **effective rate** $\rho = g(r)\, r$:

$$
b_1 = \frac{2}{s + \sqrt{s^2 - 4 k_2}}, \quad b_2 = 1 - k_2 b_1, \qquad
k_1 = \frac{\rho\, \tau_{rec_1}}{\tau_{rel}}, \; k_2 = \frac{\rho\, \tau_{rec_2}}{n_2 \tau_{rel}}, \; s = 1 + k_1 + k_2 \tag{4.8}
$$

Derivation: (4.5) $= 0$ gives $b_2(1 - b_1) = k_1 b_1$ with the new $k_1$;
(4.6) is unchanged and gives $b_2 = 1 - k_2 b_1$ with the new $k_2$; the
quadratic is the same. The output is then

$$
\theta_{ss}(r) = r\, g\, b_1 = \rho\, b_1(\rho) = \Theta_{2}\big( g(r)\, r \big) \tag{4.9}
$$

where $\Theta_2(\rho) = \rho\, b_1(\rho)$ is the steady output of section 2.
**At steady state, release-coupled facilitation is a stretch of the rate
axis, $r \to g(r)\, r$, and nothing else.** Closed form.

Consequences:

- The sustained ceiling is unchanged: $\Theta_2$ saturates at (2.8) whatever
  $G$ is. Facilitation moves the output toward the ceiling at lower rates; it
  does not raise the ceiling. (For the chosen depression parameters I checked
  numerically that $\Theta_2$ is increasing on $0 < \rho \le 4$.)
- At low rate, to second order,
  $\theta_{ss} \approx r\,[\,1 + ((G - 1)\kappa - \tau_{rec_1}/\tau_{rel})\, r\,]$.
  The synapse is net facilitating at low rate only if
  $(G - 1)\,\tau_{dec}/\tau_{fac} > \tau_{rec_1}/\tau_{rel}$, which is 8 with
  the chosen depression parameters. The same condition holds for option A.

Check, with the chosen chain and the illustrative facilitation numbers
(residuals of (4.5), (4.6) at these values are of order $10^{-17}$):

| $r$ | $g$ | No facilitation, $r b_1$ (sec. 2) | Option A, $r g b_1$ | Option B: $b_1$ | Option B: $r g b_1$ |
|---:|---:|---:|---:|---:|---:|
| 0.10 | 1.36 | 0.042 | 0.057 | 0.334 | 0.045 |
| 0.25 | 1.70 | 0.050 | 0.085 | 0.123 | 0.052 |
| 0.50 | 2.01 | 0.053 | 0.106 | 0.054 | 0.054 |
| 1.00 | 2.31 | 0.054 | 0.125 | 0.024 | 0.055 |

In option A the sustained output more than doubles at high rate. In option B
it stays under the ceiling $1/18 = 0.056$.

**Limits.** $G = 1$: section 2 exactly. $n_2 \to \infty$: the single pool with
Tsodyks-Markram coupling. Option B does **not** reduce to section 1 with
facilitation on, because section 1 has no $g$ in its depletion term; it
reduces to section 1 only with $G = 1$ and $n_2 \to \infty$.

*Strict Tsodyks-Markram* is the special case in which one number, the resting
utilisation $U$, sets the release time, the build-up time and the ceiling:
$\tau_{rel} = \tau_{fac}$ and $G = 1/U$ (from the change of variables at the
end of section 1). Option B leaves the three free.

**Ranges.** $b_1$: at 0 the derivative is $b_2/\tau_{rec_1} \ge 0$; at 1 it is
$-g r/\tau_{rel} \le 0$. $b_2$ and $g$: their equations are unchanged, so the
arguments of sections 1 and 2 hold. All states stay in range.

**Jacobian relative to section 2.** Everything listed for option A, and in
addition:

| New or changed | Entry |
|---|---|
| New block $(b_1, g)$ | $\partial \dot b_1/\partial g = -b_1\, r/\tau_{rel}$, one entry per neuron |
| $(b_1, b_1)$ | $-b_2/\tau_{rec_1} - g\, r/\tau_{rel}$ (factor $g$) |
| $\partial \dot b_1/\partial r$ | $-g\, b_1/\tau_{rel}$ (factor $g$), in $(b_1, x)$ and $(b_1, a)$ |
| $(g, b)$ | zero |

One new block. `jacobian_times` would need the same.

**Biophysical reading.** This is the standard account: release probability
rises with activity, and a higher release probability depletes the releasable
pool faster. Tsodyks, Pawelzik and Markram (1998, corpus) state that
facilitation and depression are "intricately interconnected since stronger
facilitation leads to higher $U_{SE}^1$ values which in turn leads to stronger
depression". Abbott and Nelson (2000, corpus) say the same of a raised release
probability. That at high rate the sustained response is set by the refill
rate and not by release probability is in Zador and Dobrunz (1997, corpus).
Caveat level: moderate. It is a mean-field rate form of a spiking model, and
the Campagnola result above questions depletion as the main cortical
mechanism in mouse.

**For this project.**

- Strengths: facilitation becomes self-limiting, so it cannot raise the
  sustained gain above the depression ceiling; the pool reading of the chain
  is kept whole; closed-form steady state; exact reduction to section 2 at
  $G = 1$; one extra Jacobian block.
- Weaknesses: it is a **change of behaviour** from today's facilitation, so
  it would have to be opt-in to keep `single_neuron_stf` and the regression
  references bit-identical. It still says nothing about 1TS versus MTS for
  facilitation. With the chosen depression parameters a route is net
  depressing unless $(G-1)\kappa > 8$.

---

## 5. Option C: serial facilitation, calcium then release probability, release-coupled

The literal chain analogue. Plan A's sketch (section 9, item 4), normalised
and with the symbol changed.

$$
\tau_d\, dx_i = \Big( -x_i + u_i + \sum_{j=1}^{N} w_{ij}\, \theta_j \Big)\, dt + \sigma_u\, dW_i \tag{5.1}
$$

$$
r_i = \phi\Big( x_i - a_{0_i} - \frac{c}{K} \sum_{k=1}^{K} a_{ik} \Big) \tag{5.2}
$$

$$
\theta_i = r_i\, g_i\, b_{i1} \tag{5.3}
$$

$$
\frac{da_{ik}}{dt} = \frac{-a_{ik} + r_i}{\tau_{a,ik}}, \qquad k = 1, \dots, K \tag{5.4}
$$

$$
\frac{db_{i1}}{dt} = \frac{b_{i2}\,(1 - b_{i1})}{\tau_{rec_1}} - \frac{g_i\, b_{i1}\, r_i}{\tau_{rel}} \tag{5.5}
$$

$$
\frac{db_{i2}}{dt} = \frac{1 - b_{i2}}{\tau_{rec_2}} - \frac{1}{n_2}\, \frac{b_{i2}\,(1 - b_{i1})}{\tau_{rec_1}} \tag{5.6}
$$

$$
\frac{dg_i}{dt} = \frac{1 - g_i}{\tau_{dec}} + \frac{(G - g_i)\, h_i}{\tau_{fac}} \tag{5.7}
$$

$$
\frac{dh_i}{dt} = \frac{-h_i + r_i}{\tau_h} \tag{5.8}
$$

**Symbols.** $x_i$ dendritic potential; $u_i$ external input; $w_{ij}$ weight
from $j$ to $i$; $\theta_j$ synaptic output of $j$; $\tau_d$ dendritic time
constant; $\sigma_u$ noise amplitude; $W_i$ Wiener processes; $r_i \in [0, 1]$
firing rate before depression; $\phi$ sigmoid; $a_{0_i}$ setpoint; $c$ total
SFA strength; $K$ SFA variables $a_{ik}$, time constants $\tau_{a,ik}$;
$b_{i1} \in [0,1]$ releasable pool and $b_{i2} \in [0,1]$ reserve, both
starting at 1; $\tau_{rec_1}$, $\tau_{rec_2}$ refill times; $\tau_{rel}$
release time at resting release probability; $n_2$ relative reserve capacity;
$g_i \in [1, G]$ release probability relative to rest, starting at 1;
$\tau_{dec}$, $\tau_{fac}$, $G \ge 1$ its decay time, build-up time and
ceiling; $h_i \in [0, 1]$ residual-calcium trace, starting at 0; $\tau_h$ its
time constant.

**1TS and MTS.**

- Depression: 1TS has $b_{i1}$ only ((5.5) with $b_{i2} \equiv 1$); MTS adds
  the reserve (5.6). Shared equation: (5.5).
- Facilitation: 1TS has $g_i$ only, obeying (5.7) with $h_i$ replaced by
  $r_i$, which is today's equation (1.6). MTS adds the calcium trace (5.8).
  Shared equation: (5.7).
- So in both mechanisms **MTS is 1TS plus one variable that modifies one term
  of the shared equation**: $b_{i2}$ multiplies the refill of $b_{i1}$; $h_i$
  replaces the drive of $g_i$. In both, a limit recovers 1TS:
  $n_2 \to \infty$ and $\tau_h \to 0$.
- States per synapse: 2 in 1TS ($b_1$, $g$), 4 in MTS ($b_1$, $b_2$, $g$,
  $h$). $h_i$ depends only on $r_i$, so it is the same on every route of a
  neuron and could be stored once per neuron.

(The other nesting is possible: keep $h$ in 1TS and make $g$ an
instantaneous function of it. I show the one in which 1TS facilitation is
today's equation, as 1TS depression is.)

**Relation to Plan A's sketch.** The sketch has $\dot c = -c/\tau_{ca} + r$.
Put $h = c/\tau_{ca}$: this is (5.8), with the sketch's $\tau_{fac}$ equal to
$\tau_{ca}$ times mine. Same model. The normalised form is used here because
$h \to r$ at steady state, which makes the 1TS limit and the steady state
immediate, and because $c$ is already the SFA strength.

**Coupling.** As option B: $g$ multiplies both the output and the release
term. The causal order is $r \to h \to g \to (b_1, \theta)$. In MTS the rate
no longer enters the $g$ equation directly.

**Steady state.** (5.8) $= 0$ gives $h = r$. Then (5.7) is (4.7), so $g$ is
(1.7), and $b_1$, $b_2$, $\theta_{ss}$ are (4.8), (4.9). **The steady state
is option B's at every constant rate**, and the facilitation steady state is
the same in 1TS and MTS. The calcium stage changes the time course only. This
is the facilitation counterpart of $c/K$ for SFA. It is unlike the depression
chain, whose reserve does change the steady state.

**What the second variable does to the time course.** Linearise (5.7), (5.8)
about a constant rate $r_0$, with $\lambda_g = 1/\tau_{dec} + r_0/\tau_{fac}$
and $A = (G - g_0)/\tau_{fac}$. In the Laplace variable $\sigma$:

$$
\frac{\delta g}{\delta r} = \frac{A}{\sigma + \lambda_g} \quad \text{(1TS)}, \qquad
\frac{\delta g}{\delta r} = \frac{A}{(\sigma + \lambda_g)(1 + \sigma \tau_h)} \quad \text{(MTS)} \tag{5.9}
$$

Same gain at $\sigma = 0$. The MTS response to a brief rate pulse is a
difference of two exponentials: it starts from zero, peaks later, then decays.
So the serial chain adds a **delay and a steeper low-pass**, not a slow tail.
Whichever of $\tau_h$ and $1/\lambda_g$ is longer sets the decay.

**Limits.** $\tau_h \to 0$: option B. $G = 1$: section 2 (with $h$ decoupled).
$n_2 \to \infty$: single pool. As option B, it does not reduce to section 1
with facilitation on.

**Ranges.** $h$: at 0 the derivative is $r/\tau_h \ge 0$; at 1 it is
$(r - 1)/\tau_h \le 0$, since $r \le 1$. $g$: at 1 the derivative is
$(G - 1)h/\tau_{fac} \ge 0$; at $G$ it is $(1 - G)/\tau_{dec} \le 0$. $b_1$,
$b_2$: as option B. All states stay in range.

**Jacobian relative to section 2.**

| New or changed | Entry |
|---|---|
| $h$ depends on rate | $\partial \dot h/\partial r = 1/\tau_h$, giving $(h, x)$ and $(h, a)$ |
| $h$ self-coupling | $-1/\tau_h$ |
| New block $(g, h)$ | $\partial \dot g/\partial h = (G - g)/\tau_{fac}$ |
| $g$ self-coupling | $-1/\tau_{dec} - h/\tau_{fac}$ |
| $g$ depends on rate | **no longer** in MTS (it does in 1TS) |
| New block $(b_1, g)$ | $-b_1\, r/\tau_{rel}$ |
| $(b_1, b_1)$, $\partial \dot b_1/\partial r$ | factor $g$, as option B |
| $g$ feeds the output | $\partial\theta/\partial g = r\, b_1$ |
| $h$ feeds the output | no |

Two new state blocks and three new coupling blocks.

**A remark on $h$.** Equation (5.8) is the SFA equation (5.4). If $\tau_h$
equals one of a neuron's $\tau_{a,ik}$, then $h_i = a_{ik}$ at all times (both
start at 0). Facilitation could then be driven by an existing SFA variable at
no cost in states. I note it; I do not propose it. Somatic and terminal
calcium are different compartments.

**Biophysical reading.** The residual-calcium account of facilitation:
calcium left in the terminal after a spike raises the probability of release
at the next. Tsodyks, Pawelzik and Markram (1998, corpus) motivate their
utilisation variable this way, citing Bertram et al. (1996). Torres and
Kappen (2013, corpus, citing Mongillo et al. 2008) and Hu et al. (2021,
corpus) describe facilitation as mediated by residual calcium. The two-stage
form, a calcium trace followed by a slower site with its own kinetics, is how
I recall the model of Dittman, Kreitzer and Regehr (2000) and later work on a
slow calcium sensor (**from memory**; not in the corpus). Caveat level:
moderate for the picture, high for the equations. $h$ as a linear filter of a
normalised rate is a convenience, and nothing here is fitted.

**For this project.**

- Strengths: the cleanest nesting, and the same construction as the reserve;
  1TS facilitation is today's equation; equal steady state in 1TS and MTS, so
  a difference between conditions is a difference in time course; a calcium
  variable is available if calcium-dependent recovery is ever wanted.
- Weaknesses: the second variable adds a lag, not a second timescale of
  decay, so "multiple-timescale facilitation" would be a weak description of
  it. Four states per synapse in MTS. A band of $N$ eigenvalues near
  $-1/\tau_h$ joins the spectrum; harmless unless $\tau_h$ is the slowest
  time constant in the network. Whether the lag changes stability at all is
  untested.

---

## 6. Option D: parallel facilitation, two timescales averaged, release-coupled

Facilitation and augmentation side by side, normalised as SFA is.

$$
\tau_d\, dx_i = \Big( -x_i + u_i + \sum_{j=1}^{N} w_{ij}\, \theta_j \Big)\, dt + \sigma_u\, dW_i \tag{6.1}
$$

$$
r_i = \phi\Big( x_i - a_{0_i} - \frac{c}{K} \sum_{k=1}^{K} a_{ik} \Big) \tag{6.2}
$$

$$
\theta_i = r_i\, \bar g_i\, b_{i1}, \qquad \bar g_i = \frac{1}{L} \sum_{l=1}^{L} g_{il} \tag{6.3}
$$

$$
\frac{da_{ik}}{dt} = \frac{-a_{ik} + r_i}{\tau_{a,ik}}, \qquad k = 1, \dots, K \tag{6.4}
$$

$$
\frac{db_{i1}}{dt} = \frac{b_{i2}\,(1 - b_{i1})}{\tau_{rec_1}} - \frac{\bar g_i\, b_{i1}\, r_i}{\tau_{rel}} \tag{6.5}
$$

$$
\frac{db_{i2}}{dt} = \frac{1 - b_{i2}}{\tau_{rec_2}} - \frac{1}{n_2}\, \frac{b_{i2}\,(1 - b_{i1})}{\tau_{rec_1}} \tag{6.6}
$$

$$
\frac{dg_{il}}{dt} = \frac{1 - g_{il}}{\tau_{dec_l}} + \frac{(G - g_{il})\, r_i}{\tau_{fac_l}}, \qquad l = 1, \dots, L \tag{6.7}
$$

**Symbols.** $x_i$ dendritic potential; $u_i$ external input; $w_{ij}$ weight
from $j$ to $i$; $\theta_j$ synaptic output of $j$; $\tau_d$ dendritic time
constant; $\sigma_u$ noise amplitude; $W_i$ Wiener processes; $r_i \in [0, 1]$
firing rate before depression; $\phi$ sigmoid; $a_{0_i}$ setpoint; $c$ total
SFA strength; $K$ SFA variables $a_{ik}$, time constants $\tau_{a,ik}$;
$b_{i1} \in [0,1]$ releasable pool and $b_{i2} \in [0,1]$ reserve, both
starting at 1; $\tau_{rec_1}$, $\tau_{rec_2}$ refill times; $\tau_{rel}$
release time at resting release probability; $n_2$ relative reserve capacity;
$L$ facilitation variables $g_{il} \in [1, G]$, starting at 1, with decay
times $\tau_{dec_l}$ and build-up times $\tau_{fac_l}$; $G \ge 1$ their common
ceiling; $\bar g_i \in [1, G]$ their mean, the release probability relative to
rest.

Equation (6.7) is today's (1.6). What changes is how the components combine:
their **mean**, where the code takes their product. The mean is the same as a
sum of enhancements with the budget split,
$\bar g = 1 + \sum_l (g_l - 1)/L$.

**1TS and MTS.**

- Depression: 1TS has $b_{i1}$ only ((6.5) with $b_{i2} \equiv 1$); MTS adds
  (6.6). Shared equation: (6.5).
- Facilitation: 1TS has $L = 1$, so $\bar g = g_{i1}$ and the model is option
  B. MTS has $L = 2$: a fast and a slow component, the slow one with both
  $\tau_{dec}$ and $\tau_{fac}$ longer by the same factor. Shared equation:
  (6.7) for $l = 1$.
- This is nested the way SFA is ($K = 1$ versus $K = 3$), not the way the
  reserve is. MTS is not 1TS plus a mechanism: it splits the same budget over
  two time constants. 1TS is recovered when the two time constants are equal.
- States per synapse: 2 in 1TS ($b_1$, $g_1$), 4 in MTS ($b_1$, $b_2$, $g_1$,
  $g_2$).

**Coupling.** As option B, with $\bar g$ in place of $g$: it multiplies the
output and the release term.

**Steady state.** Each component obeys (1.7):
$g_l = (1 + G \kappa_l r)/(1 + \kappa_l r)$, $\kappa_l = \tau_{dec_l}/\tau_{fac_l}$.
If the components share one ratio $\kappa$, every $g_l$ has the same steady
value, so their mean equals it:

$$
\bar g(r) = \frac{1 + G \kappa r}{1 + \kappa r} \quad \text{for every } L \tag{6.8}
$$

Then $b_1$, $b_2$, $\theta_{ss}$ are (4.8), (4.9) with $\rho = \bar g(r)\, r$.
With equal $\kappa$, **the steady state is option B's at every constant rate
and is the same in 1TS and MTS**, exactly as $c/K$ does for SFA. With unequal
$\kappa_l$ the steady state is still closed form (the mean of $L$ terms) but
differs between conditions.

**What the second variable does to the time course.** Linearising (6.7) about
$r_0$, with $\lambda_l = 1/\tau_{dec_l} + r_0/\tau_{fac_l}$ and
$A_l = (G - g_0)/\tau_{fac_l}$:

$$
\frac{\delta \bar g}{\delta r} = \frac{1}{L} \sum_{l=1}^{L} \frac{A_l}{\sigma + \lambda_l} \tag{6.9}
$$

At $\sigma = 0$ each term is $(G - g_0)\kappa/(1 + \kappa r_0)$, the same for
every $l$, so the gain equals the 1TS gain. The response to a brief pulse is a
**sum** of exponentials: an immediate rise and a decay with a fast and a slow
part. This is the structure the SFA ladder has, and the opposite of option
C's cascade.

**Limits.** $L = 1$: option B. $G = 1$: section 2. Equal time constants: 1TS.
Replacing the mean by the product and deleting $\bar g$ from (6.5) gives
section 1's facilitation.

**Ranges.** Each $g_l \in [1, G]$ by the argument of section 1, and a mean of
numbers in $[1, G]$ is in $[1, G]$. $b_1$, $b_2$: as option B.

**Jacobian relative to section 2.**

| New or changed | Entry |
|---|---|
| Each $g_l$ depends on rate | $(G - g_l)/\tau_{fac_l}$, giving $(g, x)$ and $(g, a)$ |
| $g_l$ self-coupling | $-1/\tau_{dec_l} - r/\tau_{fac_l}$, diagonal; no $g_1$ to $g_2$ coupling |
| Each $g_l$ feeds the output | $\partial\theta/\partial g_l = r\, b_1/L$ |
| New block $(b_1, g_l)$ | $-b_1\, r/(L\, \tau_{rel})$ |
| $(b_1, b_1)$, $\partial \dot b_1/\partial r$ | factor $\bar g$ |

The $g$ blocks are those of the existing product code with $1/L$ in place of
the product-rule factor, plus the one new $(b_1, g)$ block.

**Biophysical reading.** Synaptic enhancement has components on different
timescales: facilitation over tens to hundreds of milliseconds, augmentation
over seconds, post-tetanic potentiation over tens of seconds or more (**from
memory**, Zucker and Regehr 2002; not in the corpus). Tsodyks, Pawelzik and
Markram (1998, corpus) note that a facilitation factor decays "possibly at
several rates", citing Mallart and Martin (1967) and Zengel and Magleby
(1982). Caveat level: moderate for parallel components, **high for the
averaging**. As I recall it, the classical descriptions combine augmentation
with facilitation multiplicatively, not as a mean (from memory, and I am not
certain of it). The mean is this project's normalisation, chosen so that
adding a component does not add strength.

**For this project.**

- Strengths: a genuine second timescale of decay; equal steady state in 1TS
  and MTS; the same logic as $c/K$, which the manuscript already explains;
  the component equation is the code's existing one.
- Weaknesses: not a chain, so facilitation and depression would be built on
  two different principles in one model; MTS is a budget split, not an added
  mechanism; the averaging has no physiological basis I can cite; four states
  per synapse.

---

## 7. Candidates considered and not presented as options

- **Strict Tsodyks-Markram utilisation.** Not a separate option: it is option
  B with $\tau_{rel} = \tau_{fac}$ and $G = 1/U$ (sections 1 and 4).
- **Product of facilitation components** (today's code with $L = 2$). Adding a
  component multiplies the steady enhancement, and the ceiling becomes $G^2$.
  It confounds timescale count with strength, the reason the product was left
  for depression.
- **Output-only versions of options C and D.** Delete $g$ from the release
  term. They inherit option A's weakness (nothing limits the sustained gain)
  and gain nothing else.
- **Calcium-dependent recovery**: the calcium trace also speeds the refill of
  the releasable pool (as I recall Dittman and Regehr 1998; from memory). It
  changes depression, not facilitation, and adds a third coupling. Out of
  scope; option C's $h$ would make it possible later.
- **A slow variable downstream of $g$**, for example a slowly rising resting
  level. It nests like the reserve, but I know of no physiological account of
  it, and it behaves much like option D.

---

## 8. Comparison and recommendation

| | Sec. 1, today | Sec. 2, chain | A: one $g$, output only | B: one $g$, coupled | C: serial $h \to g$, coupled | D: parallel mean, coupled |
|---|---|---|---|---|---|---|
| States per synapse, 1TS | 1 (paper: no fac.) | 1 | 2 | 2 | 2 | 2 |
| States per synapse, MTS | 2 (paper: no fac.) | 2 | 3 | 3 | 4 | 4 |
| Depression: 1TS nested in MTS | yes, as an added factor ($\tau_{rel_2} \to \infty$) | yes, $n_2 \to \infty$ | yes | yes | yes | yes |
| Facilitation: 1TS versus MTS | by $L$, product | none | identical | identical | nested, $\tau_h \to 0$ | budget split, as SFA |
| Facilitation steady state equal in 1TS and MTS | no | n/a | yes (same variable) | yes (same variable) | yes | yes, if equal $\kappa$ |
| Facilitation changes depletion | no | n/a | no | yes | yes | yes |
| Sustained output at high rate | $M = 2$: falls with rate; facilitation multiplies it by up to $G^L$ | rises to $1/18$ | times $G$ | unchanged | unchanged | unchanged |
| Closed-form steady state | yes | yes (quadratic) | yes | yes (rate stretch) | yes (as B) | yes (as B, equal $\kappa$) |
| Second variable gives | n/a | slow recovery, deeper depression | n/a | n/a | lag, steeper filter | slow decay tail |
| Reduces to section 2 | no | is it | $G = 1$ | $G = 1$ | $G = 1$ | $G = 1$ |
| Reduces to section 1 with facilitation on | is it | n/a | yes, $n_2 \to \infty$ | no | no | no |
| New Jacobian blocks beyond section 2 | n/a | n/a | $g$ blocks (exist in code) | plus $(b_1, g)$ | plus $h$ blocks, $(g, h)$, $(b_1, g)$ | plus $(b_1, g_l)$ |
| What the paper would have to explain | product deepens depression | reserve; steady states differ | why facilitation does not deplete | release probability and depletion | calcium stage; that it adds a lag | why a mean; facilitation and augmentation |

### Recommendation

1. **For this paper: section 2, no facilitation.** The paper
   uses none, and Plan A's refusal of `stf` on chain routes can stand.
2. **The decision that matters is the coupling, not the second variable.** If
   facilitation is to sit beside the chain, I recommend release coupling
   (option B's equation (4.5)). Reasons:
   - It is what makes the chain a pool model throughout: the output is
     proportional to the flux leaving the releasable pool.
   - It is self-limiting. Facilitation cannot lift the sustained gain above
     the depression ceiling. For a stability paper that is the safer
     behaviour, and it is the documented one (Tsodyks et al. 1998; Zador and
     Dobrunz 1997).
   - Its steady state is one line: section 2 at the rate $g(r)\, r$.
   - It costs one Jacobian block.
3. **For a two-variable facilitation, option C is the chain analogue and the
   one I would build**, on top of B. It is the only option in which the
   facilitation comparison has the same shape as the depression comparison:
   one shared equation, one added variable, a limit that returns 1TS. It
   keeps today's equation as 1TS.
4. **But be clear about what C buys.** Its second variable delays
   facilitation; it does not add a slow component. If "multiple-timescale
   facilitation" is meant in the sense the paper uses for SFA, a broader
   range of decay times, then option D is the form that does it, at the price
   of an averaging rule with no physiological source. I would not build both.

So: B now if anything, C if a second variable is wanted for symmetry with the
reserve, D only if a slow tail is the point.

### Open questions for Tom

1. Is facilitation to appear in this paper at all, or is this for the code
   and a later paper?
2. Coupling: is release coupling accepted as the form on chain routes? It
   would need to be an opt-in so that today's facilitation (and
   `single_neuron_stf`) stays bit-identical. Should it also be allowed on
   product routes?
3. What should MTS mean for facilitation: a lag (C) or a slow tail (D)? This
   decides between them, and I could not settle it from the mathematics.
4. Equal steady state. Options C and D give facilitation the same steady
   state in 1TS and MTS. The depression chain, by decision, does not. Is that
   asymmetry acceptable, or wanted?
5. Which routes. Tsodyks et al. (1998, corpus) describe facilitation as
   prominent at pyramidal-to-interneuron synapses, and Wu and Zenke (2021,
   corpus) use E-to-I facilitation as a stabiliser. Facilitation on E to I
   only makes an E neuron's two routes differ, and then the reservoir's
   `'synaptic'` readout is undefined (it requires a neuron's routes to be
   identical); memory capacity would need the `'rate'` readout.
6. Parameters. None are proposed here. $\tau_h$, $\tau_{dec}$, $\tau_{fac}$,
   $G$ would need values and a source. With the chosen depression, a route is
   net facilitating at low rate only if $(G - 1)\tau_{dec}/\tau_{fac} > 8$.
7. In option C, is $h$ stored per neuron or per route, and is the identity
   with an SFA variable (section 5, remark) of any interest?
8. The per-neuron spread (`std_spread`) would have to say whether it covers
   the facilitation time constants.

### Not verified

- No option was simulated. Claims about lag, tails and stability rest on the
  linearisations (5.9) and (6.9).
- Jacobian entries were not checked against finite differences.
- That the code's facilitation Jacobian blocks exist as described is taken
  from Plan A's table; I read `dynamics_fast`, `get_params`,
  `compile_synapse_config` and `initialize_state`, not
  `compute_Jacobian_fast` or `jacobian_times`.
- Monotonic rise of the section 2 output was checked numerically for the
  chosen parameters only.
- The manuscript's Supplemental Methods were not read; the comparison in
  section 1 is against Methods, "Network model", only.

## Sources

Read in the corpus (`PDF2md/StabilityPaper`):

- Tsodyks, Pawelzik and Markram (1998), "Neural networks with dynamic
  synapses": the utilisation equation, its calcium motivation, the coupling
  of facilitation to depression, "several rates" of decay, the mean-field
  rate equations.
- Barak and Tsodyks (2007), "Persistent activity in neural networks with
  dynamic synapses": the rate form $du/dt = (U-u)/t_f + U(1-u)R$,
  $dx/dt = (1-x)/t_r - u x R$.
- Varela et al. (1997), "A quantitative description of short-term
  plasticity...": $A = A_0 F D_1 D_2$, additive facilitation decaying to 1.
- Abbott and Nelson (2000), "Synaptic plasticity: taming the beast":
  depression as depletion of the readily releasable pool; higher release
  probability reduces later availability.
- Zador and Dobrunz (1997), "Dynamic synapses in the cortex": at high rate
  the steady response reflects the refill rate.
- Campagnola et al. (2022), "Local connectivity and synaptic dynamics in
  mouse and human neocortex": release-independent depression fitted most
  mouse connections better than vesicle depletion.
- Wu and Zenke (2021), "Nonlinear transient amplification...": E-to-E
  depression and E-to-I facilitation as stabilisers.
- Torres and Kappen (2013); Hu et al. (2021): facilitation described as
  mediated by residual calcium (one sentence each).

From memory, not in the corpus, to be checked before any of it is cited:

- Zucker and Regehr (2002), "Short-term synaptic plasticity": the timescales
  of facilitation, augmentation and post-tetanic potentiation.
- Dittman, Kreitzer and Regehr (2000); Dittman and Regehr (1998): a
  calcium-driven facilitation variable and calcium-dependent recovery.
- That classical descriptions (Magleby and Zengel) combine augmentation and
  facilitation multiplicatively. I am not certain of this.
- That a reserve pool refilling a readily releasable pool is a standard
  picture of the terminal. The corpus supports the readily releasable pool
  (Abbott and Nelson; Zador and Dobrunz) but I found no paper there on the
  reserve.
