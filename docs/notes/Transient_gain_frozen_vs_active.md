# Transient gain with adaptation frozen vs active: what the first measurement says

*2026-09-13, R5611351, Claude Code session `90c5825b-...`, repo
`FractionalReservoir`, branch `main`. Sections 1-6 were written after the
fast smoke (HEAD `e0810e4`); section 7 adds the medium and 3-s-horizon runs
the same evening. Data:
`data/transient_gain_smoke` (fast mode: n = 500, one seed, T = 20 s, six
regular states per regime, horizon 1 s, 11 min). Figures:
`figs/transient_gain_smoke`. Code: `run_transient_gain`,
`SRNNCellTypePairs.transient_gain`, `fig_transient_gain`,
`fig_transient_gain_excursions`. Background: `Non_normal_amplification.md`
(the measure and why the dendritic block), `Session_2026-09-12_..._handoff.md`
(the discussion that led here).*

## 1. What was measured

At a sampled state S₀ of the noisy trajectory, propagate the n dendritic
unit vectors through the tangent equation and read the x-in / x-out block
of the propagator every 20 ms:

    G(t) = ‖ P_x Φ(t) P_xᵀ ‖₂        (worst case over dendritic directions)

plus, from the same n × n block, the noise-average gain ‖·‖_F/√n, the gain
along the E/I difference mode, the E/I sum mode and the leading Lyapunov
direction at that state, and the top right singular vector at the peak
(the optimal input pattern). Three propagators Φ per state:

| name | Φ | what it is |
|---|---|---|
| J_xx frozen | e^{J_xx t} | the dendritic block fixed at S₀: the conventional rate-network linearisation, adaptation and depression held constant |
| J frozen | e^{J t} | the full Jacobian fixed at S₀: a, b, g respond linearly, the state does not move |
| active | Φ from J(S(t)) along the stored trajectory | what a tangent perturbation actually does |

Adjacent pairs isolate one effect each: J_xx frozen → J frozen is
adaptation's linear feedback at a fixed operating point; J frozen → active
is nonstationarity (the state moves). G(0) = 1 by construction; the
propagator is verified against `expm` to 2e-6 on small networks.

## 2. Numbers (1-s horizon, medians over six states)

| regime | λ₁ | J_xx frozen | J frozen | active | e^{λ₁·1 s} |
|---|---|---|---|---|---|
| no adaptation | +4.0 | 47 000 | 47 000 | 1 200 | 55 |
| single-timescale (sfa1_std1) | +0.44 | 11 000 | 143 | 30 | 1.6 |
| multiple-timescale (sfa3_std2) | −0.11 | 39 | 6.7 | 8.6 | 0.9 |

The same numbers as ratios, which is how to read them:

| regime | J_xx frozen / J frozen | J frozen / active | active / e^{λ₁ t} |
|---|---|---|---|
| no adaptation | 1 (no a, b states: N = n) | 39 | ≈ 22 |
| single-timescale | 78 | 4.8 | ≈ 20 |
| multiple-timescale | 5.8 | 0.78 | ≈ 10 |

Direction readings (active propagator, at the peak): the optimal input
pattern has a participation ratio of ~60 of 500 neurons, E fraction ~0.5,
and |cos| ≤ 0.15 with the leading Lyapunov direction, in every regime. The
E/I difference mode is the best of the named directions in every regime;
the E/I sum mode is the worst and falls below 1 once adaptation is on; the
noise-average gain is ~1 in the two adapted regimes.

Peak times: the two frozen propagators and the active one in the two
unstable regimes peak AT THE HORIZON (t_peak = 1.0 s), i.e. they are still
growing. Only the multiple-timescale active curve has a genuine peak
(0.18 s) followed by decay.

## 3. How to understand it

### 3.1 The frozen picture is a runaway, not a transient, wherever α(J_xx) > 0

In the no- and single-timescale regimes the dendritic block at a typical
state has spectral abscissa α ≈ +7 and +6 s⁻¹ (`fig_transient_amplification`).
Its propagator e^{J_xx t} therefore grows like e^{αt} times a non-normal
prefactor without bound; the "peak gain" is just the gain at the horizon.
The ω − α margin of the earlier figure is the initial slope of that curve.
So the static rate-network view does not describe a transient at all in
those regimes: it describes exponential divergence at a frozen operating
point that the real system never sits at. Comparisons involving the frozen
propagators must be made at a fixed t, not at "the peak".

### 3.2 Adaptation's linear feedback at the operating point is large

Holding the state fixed but letting a and b respond (J frozen) cuts the
1-s gain by 78× in the single-timescale regime and 5.8× in the
multiple-timescale one. Mechanism: a dendritic perturbation δx raises δr,
which raises δa on τ_a and depletes δb on τ_rec/τ_rel, and the effective
input x − c·Σa and the synaptic output r·Πb both drop. This is the effect
the instantaneous ω cannot see (it is t = 0⁺), and it is the cleanest
statement that adaptation belongs in the conceptual model of the recurrent
dynamics: same weights, same state, two orders of magnitude less
amplification in a second.

Why the reduction is smaller in the multiple-timescale regime (5.8× vs
78×): its J_xx frozen gain is already small (39 vs 11 000) because the
operating point differs -- STD has depressed W_eff and SFA has moved
neurons to a shallower part of φ, the ω = 46 vs 74 story -- so there is
less left for the feedback to remove within 1 s. The 10-s SFA rung has not
acted within this horizon at all.

### 3.3 Nonstationarity is a second, separate reduction, present even without adaptation

Without adaptation the two frozen propagators are identical (a sanity check
of the code: N = n), yet the active gain is 39× smaller than the frozen
one. The linearisation at any single point of a chaotic trajectory has
α ≈ +7, but the trajectory diverges only at λ₁ = +4: the unstable
directions rotate as the state moves, and no direction stays aligned with
the local worst case for long. Freezing the Jacobian overstates
amplification in every regime, adapted or not. In the single-timescale
regime the J-frozen curve's late log-slope is about +3 s⁻¹ against
λ₁ = +0.44: the trajectory is more stable than any of the states along it,
which is what "stabilisation through motion" looks like when adaptation
keeps moving the operating point.

### 3.4 What actually amplifies: 10-20× along a specific pattern

The honest quantity is the active gain divided by e^{λ₁ t}: the transient
amplification over and above asymptotic divergence or decay. It is 22, 20
and 10 in the three regimes. That is the Murphy-Miller / Hennequin
balanced-amplification quantity for this network, and adaptation reduces
it only by about 2×; the four-orders-of-magnitude spread in the frozen
column is mostly the change in λ₁ (+4 → −0.1) and in the operating point,
not a change in non-normality. In the stable regime the whole 8.6× is a
transient: peak at 0.18 s, then decay, the textbook shape. In the other two
regimes the transient rides on exponential divergence and never peaks
within a second.

### 3.5 The amplifying pattern is not the chaotic direction

The optimal input pattern is a ~60-neuron E/I-imbalanced pattern with
|cos| ≤ 0.15 to the leading Lyapunov direction (whose x part we compared;
the full leading vector lives ~80% in the SFA states in the adapted
regimes). Two different mechanisms, then: chaos expands along a slow,
adaptation-dominated direction; transient amplification is a fast
dendritic E/I pattern. The E/I difference mode being the best named
direction everywhere confirms the balanced-amplification prediction; the
E/I sum mode being damped (< 1) in the adapted regimes says a uniform push
engages every neuron's adaptation and is cancelled -- relevant to the
stimulation idea (Next 3), where the prior was that the sum mode is the
direction that engages adaptation without a large transient.

The noise-average gain ~1 in the adapted regimes means the network's own
additive noise is not amplified on average; only structured perturbations
are, 10-30× along the optimal pattern. Adaptation makes the network
selectively sensitive rather than globally quiet.

### 3.6 Active can exceed J frozen (the stable regime: 8.6 vs 6.7)

The plan expected the frozen gains to bound the active one; they do not.
The local Lyapunov rate fluctuates around λ₁ (positive about a third of
the time in this regime, in excursions of a few hundred ms), and a moving
state passes through stretches with more expansion than the sampled point
had; time-varying linear systems can amplify more than any of their frozen
coefficients allow. The same ordering appeared on the 40- and 60-neuron
test networks. With six states from one seed the 8.6 vs 6.7 difference is
inside the state-to-state spread, so it is a caution ("do not call the
frozen gains bounds; report all three"), not yet a result.

## 4. What this means for the paper

* The conventional rate-network linearisation (J_xx) overstates
  amplification by orders of magnitude in the adapted regimes, and by 40×
  even without adaptation. Say so, with the three-curve figure.
* Adaptation's contribution at the operating point is the J_xx-frozen vs
  J-frozen ratio (78× and 5.8× at 1 s). Its contribution to the asymptotic
  rate is λ₁ (+4 → +0.44 → −0.11). Its contribution to the non-normal
  prefactor is modest (22 → 20 → 10). Keep these three apart; the frozen
  picture conflates them.
* The number to quote for "how much does a dendritic perturbation grow
  before adaptation reverses it" is the active peak in the stable regime:
  8.6× at 0.18 s (one seed; get the medium-mode value). In the unstable
  regimes quote the gain at a fixed t together with e^{λ₁ t}.
* The amplifying pattern is an E/I difference pattern over ~60 neurons,
  distinct from the Lyapunov direction; the E/I sum pattern is damped.
* Keep the ω − α figure as the conventional view, captioned as the
  initial slope of the frozen curve.

## 5. Caveats

One seed, six states, 1-s horizon. Tangent-space linearisation along a
noisy trajectory: G is an upper envelope of real amplification (saturation
ignored). The horizon is too short for the 10-s SFA rung. The excursion
sampling (Next 2) found only 0-2 onset / quiet states per regime in 20 s;
the onset-vs-quiet contrast is untested.

## 6. What to do next

1. Run the stage at medium (T = 40 s, 2 seeds, 12 regular + up to 8 onset
   + 8 quiet, horizon 2 s, ~1 h):
   `run_transient_gain('run_mode','medium','out_dir', fullfile('data','topk_med','transient_gain'))`,
   then the two figures with `run_dir` = `data/topk_med`. This gives the
   resolved peak in the stable regime, the seed spread, and the first real
   onset-vs-quiet contrast.
2. Record the spectral abscissa of the FULL J at every sampled state (a few
   `eigs` per state, or the late log-slope of the J-frozen curve) next to
   λ₁: "the trajectory is more stable than any of its states" is a claim
   worth a number.
3. A 3-s horizon on fewer states, so the 10-s rung begins to act and the
   single-timescale curve can be seen to peak or not.
4. Next 3 (stimulation directions): the E/I sum mode's damping supports the
   prior; the intermittent-but-stable networks are listed by
   `find_intermittent_stable` (multiple-timescale regime at
   `mu_IE_relative` 2.06 or `mu_EE_relative` 13.4-15.7; single-timescale at
   n = 190).
## 7. Medium run and 3-s horizon (2026-09-13 evening): the picture holds, with numbers

Two further runs after the smoke, both with the stage now recording α(J_xx),
ω(J_xx) and α(J) at every sampled state (commit `f18d83f`):

* **Medium** (`data/topk_med/transient_gain`, `figs/topk_med`): T = 40 s,
  2 seeds, 12 regular + up to 8 onset + 8 quiet states per seed, horizon
  2 s, 51 min. 24 regular states per regime.
* **3-s horizon** (`data/transient_gain_h3`, `figs/transient_gain_h3`):
  T = 40 s, 1 seed, 4 regular states per regime, horizon 3 s, 21 min.

(The first medium attempt died after 51 min at the assembly step -- the
per-seed sample arrays have different lengths and were concatenated
horizontally; fixed in `43ef72e`, which also saves the raw trials before
assembly and makes the stage test run two seeds.)

### 7.1 The operating point vs the trajectory (medium, medians over 24-40 states)

| regime | λ₁ (trajectory) | α(J_xx) at state | ω(J_xx) at state | α(J) at state |
|---|---|---|---|---|
| no adaptation | +3.35 [2.63, 4.08] | +6.0 [1.9, 10.7] | 80 | +6.0 (= α(J_xx), N = n) |
| single-timescale | +0.48 | +5.9 [2.7, 10.0] | 80 | +1.8 [−0.3, +4.7] |
| multiple-timescale | −0.11 | +0.7 [−0.8, +3.9] | 50 | +0.5 [−0.2, +1.9] (21 of 31 states unconverged) |

Three facts in one table. (i) Adaptation's linear feedback at the operating
point lowers the frozen growth rate from +5.9 to +1.8 s⁻¹ (single) and
from +0.7 to +0.5 (multiple): that is what the J_xx-frozen vs J-frozen gap
is made of. (ii) **The trajectory is more stable than any of its states**:
α(J) at a typical state is +1.8 where λ₁ is +0.48, and +0.5 where λ₁ is
−0.11. Without adaptation the same gap (+6.0 vs +3.35) exists and is pure
nonstationarity; with adaptation the operating point keeps moving away
from its own unstable directions. (iii) ω − α is ~75 s⁻¹ in every regime,
i.e. the instantaneous non-normal margin barely changes with adaptation
even though the finite-time gain changes by orders of magnitude -- the
margin is the wrong number to summarise the effect with.

Caveat on α(J): `eigs(J, 6, 'largestreal')` on the sparse 4000 × 4000
often converges only some of the six; the largest converged real part is
taken, and in the stable regime two thirds of the states returned none.
Treat it as an estimate (the J-frozen curve's late log-slope agrees:
~+1.5 to +2 s⁻¹ in the single-timescale regime).

### 7.2 Gains (medium, horizon 2 s, 24 regular states, median [min, max])

| regime | J_xx frozen | J frozen | active | active t_peak | active / e^{2λ₁} |
|---|---|---|---|---|---|
| no adaptation | 3×10⁶ [5×10², 10¹⁰] | same | 11 600 [650, 1.4×10⁶] | 2.0 (horizon) | 14 |
| single-timescale | 10⁶ [10³, 4×10⁹] | 119 [4, 6×10⁴] | 48 [23, 410] | 1.94 [0.8, 2.0] | 18 |
| multiple-timescale | 22 [5, 2×10⁴] | 4.8 [3.3, 69] | 5.0 [3.6, 12.4] | 0.14 [0.10, 1.4] | 6 |

With the 3-s horizon (one seed, 4 states), the active median curve reads:

| regime | peak | at 1 s | 2 s | 3 s | noise-average at 3 s |
|---|---|---|---|---|---|
| no adaptation | 3×10⁶ at 3 s | 820 | 23 000 | 3×10⁶ | 1.4×10⁵ |
| single-timescale | 155 at 2.7 s | 37 | 55 | 122 | 6 |
| multiple-timescale | 5.1 at 0.14 s | 2.5 | 1.0 | 0.50 | 0.03 |

So: the frozen numbers past ~1 s are astronomically large and meaningless
where α > 0 (the range spans seven decades because e^{αt} with α ∈ [2, 10]
does); the active gain is what a perturbation does. In the **stable
regime a dendritic perturbation is amplified 5× within 140 ms, is back to
its initial size by 2 s, and is at half by 3 s**; the isotropic
(noise-average) response is at 3% by then. The 10-s SFA rung produces no
second rise within 3 s (there is a small shoulder near 1 s in the medium
figure). In the single-timescale regime the active gain never peaks: it
rides on λ₁ = +0.48 with a transient prefactor of ~20-50 that is still
growing at 3 s; the J-frozen curve, which the smoke had at 143 at 1 s, is
at 3.6×10⁶ by 3 s -- the frozen operating point is unstable at +1.8 s⁻¹
and the trajectory is not.

The transient prefactor active / e^{λ₁t} at 2 s is 14 / 18 / 6: the
same ~2× reduction by adaptation as the smoke showed, against four decades
in the frozen column. The direction readings are unchanged: participation
~20 (no adaptation) to ~70 (adapted) neurons, E fraction ~0.5, |cos| with
the leading Lyapunov direction 0.04-0.06 in every regime and every
propagator; E/I difference the best named direction, E/I sum below 1 when
adapted; noise-average ≈ 1 (single) and < 1 (multiple).

### 7.3 Onsets vs quiet states (Next 2): a real but modest effect, not along the Lyapunov direction

Found in 2 × 40 s: 3 / 19 / 4 onsets and 0 / 0 / 3 quiet stretches (≥ 1 s
of negative local rate) for no / single / multiple-timescale adaptation.
A 1-s quiet stretch does not occur in the single-timescale regime (the
local rate is positive about half the time), so the planned onset-vs-quiet
contrast is empty there. Post hoc, the regular states split by the local
rate at the sample (stored) give the contrast instead:

| regime | onset: n, active G_max (2 s) | regular with local rate < 0: n, G_max | rank-sum p | G along v_Lyap at the peak, onset vs contracting |
|---|---|---|---|---|
| single-timescale | 16, 62 [27, 200] | 13, 42 | 0.09 | 4.4 vs 3.1 |
| multiple-timescale | 4, 8.2 [6.2, 11.0] | 24, 5.0 | 0.02 | 0.87 vs 0.72 |
| multiple-timescale, 1-s quiets | -- | 3, 5.2 [4.6, 6.7] | 0.11 (vs onsets) | -- |

Onset states carry ~1.5× more worst-case gain than contracting states,
significant in the stable regime with 4 vs 24 states and marginal in the
single-timescale one. But the optimal direction at an onset has |cos|
0.04-0.05 with the leading Lyapunov direction, the same as anywhere else,
and the gain along the Lyapunov direction itself is small (4.4 vs 62 worst
case; 0.87, i.e. contraction, in the stable regime). So the excursions are
NOT the network recruiting its worst-case non-normal mode. The state at an
onset is somewhat more amplifying in every direction (a higher-gain
operating point: more neurons on the steep part of φ, less depressed
synapses), and the divergence itself happens along a slow,
adaptation-dominated direction that a dendritic perturbation barely
projects onto. The mechanistic account of intermittency is therefore
"operating-point excursions", not "transient amplification along the
unstable direction".

### 7.4 Updated reading for the paper

Everything in §4 stands; two sharpenings. First, the number to quote in
the stable regime is now 5.0 [3.6, 12.4] at 0.14 s over 24 states and 2
seeds, and "gone by 2 s". Second, the operating-point table (§7.1) is the
cleanest single statement of adaptation's role: at a typical state the
recurrent block alone is unstable at +6 s⁻¹, adaptation's feedback at that
state brings it to +1.8 or +0.5, and the trajectory's own rate is +0.48 or
−0.11 -- each step is a mechanism the frozen picture lacks.

### 7.5 What is still open

* `fig_transient_gain_excursions` should draw the post-hoc split
  (onset vs regular-contracting) as well as onset vs 1-s quiets, and the
  quiet threshold should be a `cfg` option (0.5 s would exist in the
  single-timescale regime). Not done; the numbers above came from a
  console analysis of the saved samples.
* α(J) via `eigs` converges poorly on the stable regime's Jacobian; a
  shift-invert or the late log-slope of the J-frozen curve would be more
  reliable.
* Only two seeds; production mode (3 seeds, 3-s horizon, 20 regular states)
  is the version for the manuscript.
* Next 3 (stimulation directions): the E/I sum mode is damped in the
  adapted regimes and the amplifying pattern is a ~70-neuron E/I-difference
  pattern, which supports the prior that a uniform push engages adaptation
  without a large transient. The candidate networks are in §6 of
  `Non_normal_amplification.md` / `find_intermittent_stable`.
