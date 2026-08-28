# Is the Product Form a Reasonable Model of Multi-Timescale STD?

Assessment of the synaptic output term $\theta_j = r_j \prod_{m=1}^{M} b_{jm}$ used in
[`../Equations_stability_paper.md`](../Equations_stability_paper.md), with

$$
\frac{db_{im}}{dt} = \frac{1-b_{im}}{\tau_{rec_m}} - \frac{b_{im}\, r_i}{\tau_{rel}}, \qquad m = 1, \dots, M
$$

Short answer: **yes, the product form is a well-precedented phenomenological model of
multi-timescale STD** — but as written it has two consequences that should be decided
deliberately rather than inherited.

---

## Why the product form is defensible

It is essentially Varela et al. (1997, *J Neurosci*, L2/3 rat V1), who fit cortical
depression as a product of independent depression factors $D_1 \cdot D_2$ with distinct
recovery constants, each knocked down multiplicatively per spike and recovering
exponentially to 1. The term $\prod_m b_{jm}$ is the rate-model continuum version of
exactly that.

The deeper justification: recovery from a depressed state, with $\delta_m = 1 - b_{im}(0)$, is

$$
\prod_{m=1}^{M} \left(1 - \delta_m e^{-t/\tau_{rec_m}}\right)
= 1 - \sum_{m=1}^{M} \delta_m e^{-t/\tau_{rec_m}} + \mathcal{O}(\delta^2)
$$

so in the weakly-depressed limit the product **reproduces the sum-of-exponentials recovery
that experiments are actually fit with**. Deviations appear only under strong depression,
and there they are multiplicative rather than additive — which is the physically sensible
direction if the timescales represent gates in series (vesicle availability × release-site
availability $\times$ presynaptic $\mathrm{Ca}^{2+}$ channel state).

---

## Caveat 1: with a shared $\tau_{rel}$, depression strength scales with $M$

Each $b_m$ carries its own full $-b_m r/\tau_{rel}$ drive, so near $b = 1$ the initial
depression rate of the product is $M \cdot r/\tau_{rel}$. Adding timescales does not just
broaden the temporal spectrum — it multiplies the depth.

At steady state for constant rate $r$:

$$
b_m^* = \frac{\tau_{rel}}{\tau_{rel} + r\,\tau_{rec_m}},
\qquad
\prod_{m=1}^{M} b_m^* = \prod_{m=1}^{M} \frac{\tau_{rel}}{\tau_{rel} + r\,\tau_{rec_m}}
$$

With $r = 1$, $\tau_{rel} = 0.5$ s, and SFA-style logspaced
$\tau_{rec} = [0.1,\, 1,\, 10]$ s:

| $\tau_{rec_m}$ (s) | $r\,\tau_{rec_m}/\tau_{rel}$ | $b_m^*$ |
|---|---|---|
| 0.1 | 0.2 | 0.833 |
| 1 | 2 | 0.333 |
| 10 | 20 | 0.048 |
| **product** | — | **0.013** |

That is a 75× attenuation, versus 3× for the single-timescale case. (The code default is
`tau_b_E_rel = 0.25`, which makes it ~290×. A now-deleted `parameter_table.md` listed $1/2$.)
The slow timescale alone nearly silences the synapse.

This is the same fairness problem already solved for SFA with $c_E = 0.15/3$. The analogue
here is a **per-timescale release constant** $\tau_{rel_m} = M\,\tau_{rel}$ (or per-timescale
$U_m$ summing to a fixed budget), so that $M = 1$ and $M = 3$ are comparable in total
depression and differ only in spectrum. At present `tau_b_E_rel` is validated as a scalar
(`SRNNModel2.m:1107`), so this would be a small API change.

---

## Caveat 2: $M \ge 2$ destroys the Tsodyks–Markram rate-limiting property

The signature property of single-timescale STD is that transmitted output *saturates*:

$$
r \, b^* = \frac{r\,\tau_{rel}}{\tau_{rel} + r\,\tau_{rec}} \;\xrightarrow[r \to \infty]{}\; \frac{\tau_{rel}}{\tau_{rec}}
\qquad \text{(independent of } r\text{)}
$$

For the product,

$$
r \prod_{m=1}^{M} b_m^* \;=\; r \prod_{m=1}^{M} \frac{\tau_{rel}}{\tau_{rel} + r\,\tau_{rec_m}}
\;\sim\; r^{\,1-M} \qquad (r \to \infty)
$$

so for $M \ge 2$ the steady-state synaptic transmission is **non-monotonic** — it rises,
peaks, then collapses toward zero. That is a qualitatively different network element: a gain
that actively shuts off at high rate rather than clamping.

For a stability / edge-of-chaos question this may well be a feature (it is a very strong
stabilizer, and a natural burst-terminator), but it should be a stated modeling choice rather
than a side effect, because it is the property most likely to drive whatever $\lambda_1$
result comes out.

If multi-exponential recovery is wanted *while preserving* rate-limiting, the alternative is a
convex combination,

$$
b_{\text{eff},i} = \sum_{m=1}^{M} w_m b_{im}, \qquad \sum_{m=1}^{M} w_m = 1
$$

which stays in $[0,1]$, recovers multi-exponentially by construction, keeps
$r \cdot b_{\text{eff}} \to \text{const}$, and makes depression strength independent of $M$
for free.

---

## Caveat 3: interpretation

For $M = 1$, $b$ is legitimately "fraction of resources available." For $M > 1$ the product is
*not* a resource-conservation model — there is no single pool, and no individual $b_m$ is a
vesicle fraction.

Genuine depletion-based multi-exponential recovery would require a serial pool cascade
(readily-releasable ← recycling ← reserve), in which refilling of the fast pool depends on the
slow pool's state. That is more biophysical and is *not* reducible to a product. Worth one
sentence in the equations doc so the product form is clearly flagged as phenomenological.

---

## Recommendation

Keep the product — it matches the standard phenomenology, and the weak-depression limit is
provably the right multi-exponential recovery. But:

1. **Split the depression budget across timescales** via per-$m$ release constants, so that
   $M$-comparisons are fair.
2. **State explicitly** in the equations doc that $M \ge 2$ yields non-monotonic steady-state
   transmission, since that is the dynamically consequential difference from classic STD.

---

## Reference

Varela, J. A., Sen, K., Gibson, J., Fost, J., Abbott, L. F., & Nelson, S. B. (1997).
A quantitative description of short-term plasticity at excitatory synapses in layer 2/3 of rat
primary visual cortex. *Journal of Neuroscience*, 17(20), 7926–7940.
