# What we learned about transient amplification and adaptation (brief report, 2026-09-13)

*Summary of the day's three runs of `run_transient_gain` on the paper's
network (n = 500, `celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25`):
the fast smoke (1 seed, 1-s horizon), the medium run (2 seeds, 2-s horizon,
`data/topk_med/transient_gain`, `figs/topk_med`) and a 3-s-horizon run
(1 seed, `data/transient_gain_h3`). Full detail and tables:
`Transient_gain_frozen_vs_active.md`; the measure and the reasoning behind
it: `Non_normal_amplification.md`. Commits `478052a` .. `4061bcf`.*

## The measure

The x-in / x-out transient gain G(t) = ‖P_x Φ(t) P_xᵀ‖₂: a unit perturbation
of the dendritic states only, measured on the dendritic states only, with
the propagator Φ taken three ways at each sampled state -- the dendritic
block J_xx frozen (the conventional rate-network linearisation), the full
Jacobian J frozen (adaptation responds, the state does not move), and
active along the trajectory (what a perturbation actually does). The same
n × n block also gives the noise-average gain, the gain along the E/I
difference and sum modes and along the leading Lyapunov direction, and the
optimal input pattern at the peak.

## Six things we learned

1. **The frozen rate-network picture describes a runaway, not a transient.**
   J_xx at a typical state is unstable at about +6 s⁻¹ in the no- and
   single-timescale regimes, so e^{J_xx t} grows without bound and its
   "peak" is the horizon; across states it spans seven decades. The ω − α
   margin of the eig-stage figure is the initial slope of that runaway and
   is ~75 s⁻¹ in every regime, so it cannot summarise adaptation's effect.

2. **The trajectory is more stable than any of its states, in two steps.**
   Single-timescale regime, medians over 40 states: recurrent block alone
   +5.9 s⁻¹, full J with adaptation's feedback +1.8, the trajectory itself
   +0.48. Stable regime: +0.7, +0.5, −0.11. The first step is adaptation's
   linear negative feedback at the operating point; the second is motion
   (the operating point moves away from its own unstable directions), and
   that second step exists even without adaptation (+6.0 vs +3.35). This
   table is the cleanest statement that adaptation belongs in the
   conceptual model of the recurrent dynamics.

3. **What actually amplifies is 10-20× along a specific pattern, and
   adaptation changes that only modestly.** Active gain over e^{λ₁t}: 14 /
   18 / 6 for no / single / multiple-timescale adaptation. The decades of
   difference in the frozen columns are mostly the change in λ₁ (+4 → −0.1)
   and in the operating point, not a change in non-normality. Stable
   regime, the number to quote: a dendritic perturbation is amplified
   5.0× [3.6, 12.4] at 0.14 s (24 states, 2 seeds), is back to its initial
   size by 2 s and at half by 3 s, with the isotropic response at 3%. The
   single-timescale regime has no peak within 3 s: its transient rides on
   λ₁ > 0.

4. **The amplifying pattern is not the chaotic direction.** The optimal
   input is a ~70-neuron E/I-difference pattern with |cos| ≈ 0.05 to the
   leading Lyapunov direction, which lives mostly in the SFA states. Chaos
   expands along a slow, adaptation-dominated direction; transient
   amplification is a fast dendritic E/I pattern. The E/I difference mode
   is the best named direction in every regime (balanced amplification
   holds), the E/I sum mode is damped below 1 once adapted, and the
   noise-average gain is ~1: adaptation makes the network selectively
   sensitive, not globally quiet.

5. **Intermittency is operating-point excursions, not recruited non-normal
   modes.** States at the onset of a local-Lyapunov-rate excursion carry
   ~1.5× more worst-case gain than locally contracting states (8.2 vs 5.0,
   p = 0.02 in the stable regime; 62 vs 42, p = 0.09 single-timescale), but
   the amplifying direction at an onset is no more aligned with the
   Lyapunov direction than anywhere else, and the gain along that direction
   is small (contraction, in the stable regime). The excursion is the state
   passing through a higher-gain operating point; the divergence happens
   along a direction a dendritic perturbation barely touches.

6. **For the stimulation idea.** A uniform (E/I sum) push is damped once
   adaptation is on while the difference pattern amplifies, which supports
   the prior that the sum mode engages adaptation without a large
   transient. `find_intermittent_stable` lists the stable-but-intermittent
   networks to build that analysis on (multiple-timescale regime at
   `mu_IE_relative` 2.06 or `mu_EE_relative` 13.4-15.7; single-timescale at
   n = 190).

## What to change in how we talk about it

Keep the J_xx ω/α figure as the conventional view, captioned as the initial
slope of the frozen curve. Report the three-propagator curves and the
operating-point table (item 2) as the evidence that adaptation is part of
the dynamics, not a correction to them. Quote the active peak in the stable
regime; in the unstable regimes quote gains at a fixed t together with
e^{λ₁t}. Keep the asymptotic rate, the operating-point feedback and the
non-normal prefactor apart; the frozen picture conflates them.

## Caveats and loose ends

Two seeds (production mode has three); tangent-space linearisation, so
every gain is an upper envelope; α(J) via `eigs` converged on only a third
of the stable regime's states (the J-frozen curve's late slope agrees);
no 1-s quiet stretch exists in the single-timescale regime, so the
onset contrast there is against locally contracting regular states, from
a console analysis that the excursion figure does not yet draw. The
horizon should be 3 s in the manuscript run.
