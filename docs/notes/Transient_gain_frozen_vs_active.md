# Transient gain with adaptation frozen vs active: what the first measurement says

*2026-09-13, R5611351, Claude Code session `90c5825b-...`, repo
`FractionalReservoir`, branch `main`, HEAD `e0810e4`. Data:
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
