# Non-normal transient amplification, and how adaptation reduces it

*2026-09-12. Companion to `Lyapunov_estimation_methods.md` §4.4 and to
`fig_transient_amplification`. What was measured on the fast smoke
(`data/topk_smoke_fast`), why the measure is taken on the dendritic block
rather than the full Jacobian, what that implies, and how the reduction by
adaptation could be quantified for the paper.*

## 1. The measure

For a linear system dz/dt = J z, the spectral abscissa α(J) = max Re λ(J)
gives the asymptotic growth rate, and the **numerical abscissa**
ω(J) = max eig((J + Jᵀ)/2) gives the largest possible *instantaneous* growth
rate of the Euclidean norm: d‖z‖/dt ≤ ω ‖z‖, with equality for the top
eigenvector of the symmetric part. For a normal J (Jᵀ commuting with J) the
two coincide; for a non-normal J, ω > α, and the gap is the room for a
perturbation to grow before the eigen-decay wins (Trefethen & Embree,
*Spectra and Pseudospectra*). Every eigenvalue can be in the left half plane
while ω is large and positive: that is transient amplification, and it is
the mechanism Murphy & Miller call *balanced amplification* and Hennequin,
Vogels & Gerstner exploit for movement generation in stable E/I networks.

The instantaneous Jacobian of the SRNN is evaluated along the trajectory
(150 states per condition in the eig-heatmap stage), so ω and α are time
series; `fig_transient_amplification` shows both, and the distribution of
ω − α per condition.

## 2. What was found (fast smoke, one seed, n = 500, medium of 40 states)

| regime | ω(J_xx) (1/s) | α(J_xx) (1/s) | ω − α | λ₁ |
|---|---|---|---|---|
| no adaptation | 82 [79, 86] | 7 | ≈ 82 | +3.94 |
| single-timescale SFA + STD | 73 [69, 77] | 6 | ≈ 73 | +0.63 |
| multiple-timescale SFA + STD | 46 [44, 48] | 0.3 | ≈ 46 | −0.11 |

(median [IQR] over sampled states.) Two things stand out. First, the
instantaneous growth bound is one to two orders of magnitude above the
spectral abscissa in every regime: this is a strongly non-normal network,
as a Dale's-law E/I network must be. Second, adaptation lowers ω from 82 to
46 s⁻¹, a 44% reduction, with the single-timescale regime in between.

Where the size comes from: J_xx = (−I + W diag(θ′))/τ_d with τ_d = 0.1 s,
so ω = (−1 + λ_max(sym(W diag θ′)))/τ_d. A Dale's-law W has a rank-one mean
structure (every E column positive, every I column negative) whose largest
*singular* value grows like √n·μ while its outlier *eigenvalues* stay near
the block means; the symmetric part inherits the singular value, so
λ_max(sym W_eff) ≈ 5–9 and ω ≈ 40–80 s⁻¹. The direction that achieves it
is the E-versus-I difference mode: a perturbation that raises E and lowers
I is amplified through the loop before the balanced sum mode decays. That
is Murphy & Miller's balanced amplification, and it is physical, not an
artefact.

## 3. Why the dendritic block J_xx, not the full Jacobian

The first implementation took ω on the full N × N Jacobian (N = 4000: x,
three SFA variables and eight STD variables per neuron). It read 50–180 s⁻¹
and ranked the regimes plausibly, but the number is not trustworthy, for a
reason that is structural rather than numerical:

**The numerical abscissa is not invariant to a diagonal rescaling of the
state.** If z = D y with D diagonal, J becomes D⁻¹ J D, whose eigenvalues
(hence α) are unchanged, but whose symmetric part is not: ω depends on the
units in which each state variable is measured. The full Jacobian mixes
rows in units of 1/τ_d = 10 s⁻¹ (x), 1/τ_a = 0.1–4 s⁻¹ (SFA) and
1/τ_rec + r/τ_rel ≈ 0.3–5 s⁻¹ (STD), and its off-diagonal couplings are
asymmetric by construction: J(x, a) = −c_eff W diag(φ′)/τ_d against
J(a, x) = diag(φ′)/τ_a, differing by a factor of order W/(c τ_d) · τ_a. The
symmetric part of that pair is dominated by whichever block happens to have
the larger units, so ω on the full J measures how we chose to write the
state vector, not how the dynamics amplify. There is no canonical
normalisation for a state whose components are a potential, a dimensionless
adaptation current and a synaptic resource.

The dendritic block J_xx has one unit throughout, and it is exactly the
object Hennequin et al. and Murphy & Miller analyse: a rate network with
effective connectivity W_eff = W diag(θ′), where θ′ = φ′(x_eff) · Π b · Π g
is the slope of the synaptic output at the operating point. So the
restriction is what makes the number comparable to the literature.

**What the restriction costs.** J_xx freezes the adaptation and depression
states. Its ω therefore captures how adaptation changes the *operating
point* (the slope φ′ where each neuron sits, and the depression factor
Π b multiplying each synapse) but not adaptation's *dynamic* negative
feedback on the timescale of the transient. The 82 → 46 s⁻¹ reduction is
the static part: multiple-timescale STD depresses the effective synaptic
gain and SFA moves neurons to shallower parts of the nonlinearity, and
together they shrink the singular values of W_eff. The dynamic part (the
adaptation current rising during a transient and cutting it short) is
invisible to this measure, and it is the part that operates on the 0.25 s
to 10 s timescales where a transient of a few τ_d would be most affected.
So the figure *understates* adaptation's effect on transient amplification,
in a direction that makes the claim conservative.

## 4. Implications

* The three regimes differ in the *bound* on transient growth by a factor
  of ~2, on top of differing in the asymptotic sign. "Stable but supports
  transient divergence" (the Introduction's claim) has a number attached:
  in the multiple-timescale regime α ≈ 0.3 s⁻¹ and λ₁ ≈ −0.11 s⁻¹, yet a
  well-chosen dendritic perturbation can grow at up to 46 s⁻¹ for a
  fraction of τ_d.
* The bound is achieved by the E/I difference mode. Whether the network's
  own fluctuations actually excite that mode is a separate question, and
  is what the trajectory-side measures answer (the fraction of time the
  local Lyapunov rate is positive; the 0.2 s finite-time exponent; the
  positive-excursion length, all stored per sweep job).
* Because ω scales with the singular values of W_eff, it will track the
  synaptic-gain and μ sweeps almost by construction. The informative
  comparison is *across adaptation regimes at fixed W*, which is how the
  stage samples (same seed for all conditions).

## 5. How to quantify adaptation's reduction, for the paper

In increasing order of effort and completeness:

1. **ω − α on J_xx, per condition, paired by network** (done). Report the
   median over states and seeds and the paired reduction relative to no
   adaptation. Cheap; static; understates the effect.
2. **A scale-free non-normality index** on the same block: the ratio of
   the largest singular value to the spectral radius of W_eff,
   σ_max(W_eff)/ρ(W_eff), or Henrici's departure from normality
   ‖W_eff‖_F² − Σ|λ_i|². These separate "the gain went down" from "the
   matrix became more normal", which ω alone conflates.
3. **Maximal transient gain with adaptation active**: G(t) =
   ‖P_x exp(J t) P_xᵀ‖₂, the largest growth of a *dendritic* perturbation
   into a *dendritic* response through the full system over time t, with
   the adaptation and depression states free to respond. Because the input
   and output are both in x, this is invariant to the state scaling that
   spoils ω on the full J, and it contains the dynamic feedback that J_xx
   omits. Compute by propagating the n-column basis [I_n; 0] through the
   tangent equation with the matrix-free `jacobian_times` for t up to ~1 s
   (about 20 s per sampled state at N = 4000; 20 states per condition is
   affordable) and taking the top singular value of the x-rows. The peak of
   G(t) over t, and the time of the peak, per condition, is the number I
   would put in the paper: "a dendritic perturbation is amplified at most
   G_max-fold, at t_peak, before adaptation reverses it".
4. **Hennequin's evoked energy**, the integral of ‖Δr‖² over the response
   to an optimal initial condition, obtained from a Lyapunov equation for
   a *stable* linearisation. Cleanest for a fixed stable J_xx (n = 500, so
   the Lyapunov solve is milliseconds), but the frozen block is not stable
   in every regime (α > 0 in two of three), and on the full J it inherits
   the scaling problem unless restricted to x-in/x-out as in item 3.

Items 1 and 2 are a few lines in `run_eig_heatmap`'s `sample_eigenvalues`;
item 3 is a small stage of its own and is the one that answers the
question as posed ("how much does adaptation reduce transient
amplification") with adaptation actually acting.

## References

* Trefethen, L. N. & Embree, M. (2005). *Spectra and Pseudospectra*.
  Princeton. (Numerical abscissa; bounds on ‖e^{tA}‖.)
* Murphy, B. K. & Miller, K. D. (2009). Balanced amplification: a new
  mechanism of selective amplification of neural activity patterns.
  *Neuron* 61, 635–648.
* Hennequin, G., Vogels, T. P. & Gerstner, W. (2012). Non-normal
  amplification in random balanced neuronal networks. *Phys. Rev. E* 86,
  011909.
* Hennequin, G., Vogels, T. P. & Gerstner, W. (2014). Optimal control of
  transient dynamics in balanced networks supports generation of complex
  movements. *Neuron* 82, 1394–1406.
* Goldman, M. S. (2009). Memory without feedback in a neural network.
  *Neuron* 61, 621–634. (Feedforward/non-normal structure hidden in
  recurrent networks.)
