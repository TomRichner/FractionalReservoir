# Notes on Multi-Timescale STD and Varela et al. (1997)

Working notes on whether the product form $\left(\prod_{m=1}^{M} b_{jm}\right) r_j$ is a
reasonable model of multi-timescale short-term depression, and what Varela et al. (1997)
actually specifies.

Source PDF and OCR markdown: `Varela_1997.pdf` / `Varela_1997.md` in this directory.

> **Supersedes** the caveats in `../MTS_STD_product_form_assessment.md`. That earlier
> document was written before reading Varela in detail; its "split the depression budget
> across timescales" and weighted-sum recommendations have been **rejected** — see
> [Design decisions](#design-decisions) below. Its core defense of the product form stands.

---

## 1. Our model

$$
\frac{dx_i}{dt} = \frac{-x_i + u_i + \sum_{j} w_{ij}\left(\prod_{m=1}^{M} b_{jm}\right) r_j}{\tau_d},
\qquad
\frac{db_{im}}{dt} = \frac{1-b_{im}}{\tau_{rec_m}} - \frac{b_{im}\, r_i}{\tau_{rel}}
$$

with $\tau_{rec}$ a $1\times M$ vector and $\tau_{rel}$ currently a **shared scalar**
(`SRNNModel2.m:1107` validates `tau_b_E_rel` as scalar; default `0.25`, though
`parameter_table.md` documents $1/2$ — worth reconciling).

## 2. What Varela actually specifies

Varela's model (their Eqs. 1–8) is a product of independent depression factors:

$$A = A_0\, F\, D_1 D_2 D_3 \tag{Varela Eq. 1}$$

Each depression factor has **two** parameters, not one:

| Varela eq. | Rule | Role |
|---|---|---|
| Eq. 7 | $D \to D\,d_m$ at each spike | per-spike depression **strength** |
| Eq. 8 | $\tau_{D_m}\,\dot{D} = 1 - D$ between spikes | **recovery** time constant |

Their line 87 states it directly: *"The different components of depression ($D_1$, $D_2$,
$D_3$) had different constant factors ($d_1$, $d_2$, $d_3$) **and** time constants."*

Equation 8 alone looks like a single time constant per component only because the second
parameter, $d_m$, is written as a per-spike multiplicative jump rather than as a rate.

## 3. Mapping Varela onto our rate model

Varela is an event-based (per-spike) model; ours is a rate model. Take the mean-field limit:
with Poisson spiking at $\nu$ Hz, a spike arrives in $dt$ with probability $\nu\,dt$ and
multiplies $D$ by $d_m$, so $\mathbb{E}[\Delta D] = -(1-d_m) D\,\nu\,dt$. Combining with Eq. 8:

$$
\boxed{\;\frac{dD_m}{dt} = \frac{1-D_m}{\tau_{D_m}} - (1-d_m)\,\nu\,D_m\;}
$$

This is **functionally identical** to our $b_{im}$ equation. The identification is

$$
\tau_{rec_m} \leftrightarrow \tau_{D_m},
\qquad
\frac{1}{\tau_{rel_m}} \leftrightarrow (1-d_m)\, r_{\max}
$$

where $r_{\max}$ (Hz) converts our dimensionless $r \in [0,1]$ into a firing rate.

**Key conclusion:** $\tau_{rel}$ *is* Varela's $d_m$ — it is the release/use strength $U_m$
written as a reciprocal time constant. It is not extra structure that Varela lacks. Our model
and Varela's both carry two parameters per depression component; this is also just the
Tsodyks–Markram $(U, \tau_{rec})$ pair under a different parameterization.

### Fitted parameter values

From the Varela control fits (strengths from their pooled statistics, line 249; $d_2$ from the
Fig. 12/13 legends; time constants from the Fig. 10–13 legends):

| Component | $\tau_{rec}$ | $d_m$ | $1-d_m$ | Character |
|---|---|---|---|---|
| $D_1$ (fast) | 0.44–0.64 s | 0.86–0.88 | ~0.13 | strong per spike, short-lived |
| $D_2$ (slow) | 4.3–9.2 s | ~0.98 | ~0.02 | weak per spike, long-lived |

The per-spike strengths differ by roughly **6×**. Converting with $r_{\max} = 20$ Hz:

$$
\tau_{rec} \approx [0.5,\; 6]\ \text{s}, \qquad \tau_{rel} \approx [0.39,\; 2.5]\ \text{s}
$$

A single shared $\tau_{rel}$ forces $d_1 = d_2$, which is exactly the constraint the data
reject — and because the slow component has the longest $\tau_{rec}$, forcing it to the fast
component's strength is what drives the near-silencing seen with naive shared values:

| $\tau_{rel}$ | $\prod_m b_m^*$ at $r=1$ |
|---|---|
| Varela-derived, per-component $[0.39, 2.5]$ s | 0.13 |
| shared 0.5 s (as documented) | 0.038 |
| shared 0.25 s (current code default) | 0.013 |

using $b_m^* = \dfrac{\tau_{rel_m}}{\tau_{rel_m} + r\,\tau_{rec_m}}$.

## 4. Two points Varela settles in favor of the product form

**The product is empirically validated, not merely convenient.** Their Figure 6 measures
recovery directly after a 15 s, 20 Hz train: *"In each case data were well fit by the product
of two exponentials (solid curves) but were not as well fit by a single exponential"*
(line 179). Product-vs-single was tested explicitly.

Analytically this is consistent with the weak-depression expansion — with
$\delta_m = 1 - b_{im}(0)$,

$$
\prod_{m}\left(1 - \delta_m e^{-t/\tau_{rec_m}}\right) = 1 - \sum_m \delta_m e^{-t/\tau_{rec_m}} + \mathcal{O}(\delta^2)
$$

so the product reduces to the sum-of-exponentials recovery that experiments are fit with, and
deviates only under strong depression.

**Non-monotonic (band-pass) transmission is a real cortical property, not an artifact.** For
$M \ge 2$ the steady-state transmitted signal $r\prod_m b_m^*$ is non-monotonic in $r$,
decaying as $r^{1-M}$ at high rate. Writing $A_m = \tau_{rec_m}/\tau_{rel_m}$, the $M=2$ peak
is at

$$
r^* = \frac{1}{\sqrt{A_1 A_2}}
$$

With the Varela-derived values ($A_1 \approx 1.30$, $A_2 \approx 2.4$) this gives
$r^* \approx 0.57$ of max rate, i.e. ~11 Hz — and Varela's discussion reports that cortical
synapses *"respond poorly to high-frequency (>10–20 Hz) fluctuations in presynaptic firing
rates"* (line 295). Previously logged as a caveat; it should be treated as a documented
feature.

## 5. Mechanistic origin

Varela's pharmacology argues that **both** components are presynaptic and release-dependent:

- Partial postsynaptic blockade (CNQX) left short-term plasticity unchanged (their Fig. 10)
- Blocking AMPA receptor desensitization (cyclothiazide) left it unchanged (Fig. 11)
- Reduced extracellular Ca²⁺ weakened both components (Fig. 12); adenosine abolished $D_2$ (Fig. 13)

This rules out postsynaptic desensitization as the slow component's substrate.

## 6. How multi-timescale STD is modeled generally

Three families, distinguished by what the timescales physically represent:

| Family | Form | Interpretation | Rate-limiting? |
|---|---|---|---|
| **Multiplicative gates in series** (Varela; ours) | $\prod_m b_m$ | one pathway, several gates | no — band-pass for $M\ge2$ |
| **Parallel heterogeneous subpopulations** | $\sum_m w_m b_m$, $\sum w_m = 1$ | a *mixture* of synapse types | yes |
| **Serial vesicle-pool cascade** | linear chain RRP ← recycling ← reserve | depletion-based, resource-conserving | yes |

The pool cascade is the most mechanistic but is a linear chain, not reducible to a product.
The parallel-subpopulation form is appropriate when a projection is a mixture of synapse
types rather than one pathway with multiple gates.

**Interpretation caveat for the product:** for $M = 1$, $b$ is legitimately "fraction of
resources available." For $M > 1$ the product is *not* a resource-conservation model — there
is no single pool, and no individual $b_m$ is a vesicle fraction. The product form should be
labelled phenomenological.

---

## Design decisions

Recorded so these are not relitigated:

1. **Keep the product form.** Directly precedented (Varela), empirically validated against
   single-exponential recovery, and correct in the weak-depression limit.

2. **Reject the normalized weighted-sum** $b_{\text{eff}} = \sum_m w_m b_m$ with
   $\sum w_m = 1$. Not wanted here.

3. **Greater suppressive depth with larger $M$ is desirable, not a bug.** The product gives
   $\prod_m b_m^* = \prod_m (1 + r\,\tau_{rec_m}/\tau_{rel_m})^{-1}$, which decreases
   monotonically with each added timescale by construction. This is the intended behavior —
   multi-timescale STD *should* suppress more deeply than single-timescale STD. There is
   therefore **no depression "budget"** to be conserved across timescales, and no analogue
   here of the SFA convention $c_E = 0.15/3$.

   Note this is independent of decision 4: making $\tau_{rel}$ a vector lets each component
   take its own empirically-grounded strength; it does **not** impose a conserved total.

4. **Make $\tau_{rel}$ a per-timescale vector** — see next steps.

---

## Next step: vectorize `tau_b_E_rel` / `tau_b_I_rel`

**Not being pursued right now.** Recorded for when it is.

Change $\tau_{rel} \to \tau_{rel_m}$:

$$
\frac{db_{im}}{dt} = \frac{1-b_{im}}{\tau_{rec_m}} - \frac{b_{im}\, r_i}{\tau_{rel_m}}
$$

This makes the model exactly Varela's in rate form, so the fitted $(d_m, \tau_{D_m})$ values
drop straight in.

Touch points in `src/SRNNModel2.m`:

- **Defaults** (~line 61) — `tau_b_E_rel = 0.25` and the `tau_b_I_rel` counterpart become
  $1\times M$ vectors.
- **Validation** (~line 1107) — currently `error(... 'tau_b_E_rel must be a scalar ...')`;
  becomes the same `numel == n_b_E` check already applied to `tau_b_E_rec` at line 1102, plus
  a `reshape(..., 1, [])`. Scalar input should broadcast to all $M$ for backward
  compatibility.
- **Dynamics** (~line 1791) — `db_E_dt = (1 - b_E)./tau_b_E_rec - (b_E .* r(E_indices))/tau_b_E_rel;`
  the final `/` becomes `./` with an implicitly-expanded $1\times M$ row vector (`b_E` is
  $n_E \times M$, so this broadcasts correctly with no reshape).
- **`Pmin` bounds** (~lines 1646, 2754) — `prod(tau_b_E_rel ./ (tau_b_E_rec + tau_b_E_rel))`
  is already elementwise and generalizes unchanged.
- **Jacobian** (~line 2653 onward) — the $\partial \dot b / \partial b$ and
  $\partial \dot b / \partial r$ blocks pick up per-$m$ $1/\tau_{rel_m}$; verify against a
  finite-difference check.
- **Parameter export** (~line 869) and any `parfor`-sliced parameter structs.

Also worth doing at the same time: reconcile the `tau_b_E_rel` default (`0.25` in code) with
`parameter_table.md` ($1/2$).

---

## Reference

Varela, J. A., Sen, K., Gibson, J., Fost, J., Abbott, L. F., & Nelson, S. B. (1997).
A quantitative description of short-term plasticity at excitatory synapses in layer 2/3 of rat
primary visual cortex. *Journal of Neuroscience*, 17(20), 7926–7940.
