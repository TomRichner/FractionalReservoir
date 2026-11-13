---
title: "Block B pilots: what they show and what they leave to decide"
header-includes:
  - \usepackage{float}
  - \floatplacement{figure}{H}
---

**What this is:** the nine block B pilot reports (0a, 0b, 1a, 1b, 2, 3, 4, 5, 6), edited into one report for planning block B.\
**Date:** 2026-11-12, all pilots run that day.\
**Status:** exploratory throughout. Mostly one seed and one window per case; nothing was tuned to produce an effect.

The figures are linked from `figs/explore_B/`, which is not in git; on another machine this Markdown reads without them (the PDF has them).

Two networks recur. **run13**: 350 neurons in five weakly coupled modules, noise on, brief localized IED-like discharges; as built it has one SFA timescale and one depression pool, so "run13 as built" is single-timescale. **The bursting network** (`bursting_pairs`): 50 neurons, no noise, large spontaneous population bursts; as built it has three SFA timescales on E and one depression pool. "Chained MTS" is three SFA timescales with depression as a releasable pool refilled from a reserve; "product MTS" is the same SFA with depression as the product of two independent pools. Where a network's own values exist they were kept; constants only the MTS forms have came from the mu7chain bundle's defaults.

# Summary

1. **Removing adaptation gives no fair "adaptation ablated" arm on either network.** Removing depression, or replacing adaptation by constants equal to its time averages, saturates or silences the unstimulated network; the one removal that keeps comparable events (SFA removed on run13) tests SFA alone. Clamping adaptation from stimulus onset is the only candidate left (pilot 3).
2. **Clamping works today, exactly, by script, with no change to the model** (pilot 2). With adaptation clamped, a stimulus that quiets the bursting network instead drives it to saturation (mean E rate 0.906 against 0.359 intact); on the paper's network free adaptation absorbs about two thirds of the rate increase. One seed each: feasibility and direction, not size.
3. **Stimulation still suppresses aberrant activity under chained MTS, on both networks.** run13: events fall by 83 % at input 0.20, as under run13 as built, lower in 8 of 8 seeds (pilot 1a). Bursting network: every all-neuron step from 0.1 to 0.4 stops the bursts (pilot 1b).
4. **Suppression is steady elevated firing, not quiet.** The mean rate rises when the events stop (bursting network under chained MTS 0.060 to 0.200; run13 mean E rate +62 %), with adaptation recruited. "Suppressed" has to mean "aberrant events replaced by steady, adapted firing" (pilots 1a, 1b, 5).
5. **Chained and product MTS are not alike on these networks.** At the mu7chain constants the product nearly silences both: run13 keeps 1 event in 100 s, the bursting network none. Chained MTS keeps both networks' aberrant activity (pilots 0a, 0b).
6. **One measure, excess co-activation, works for bursts and IED-like events,** read with the mean rate and the fraction of saturated neurons. It orders both networks as their own detectors do and falls under stimulation. But by the proposed threshold the paper's own chained-MTS network is hyperexcitable at baseline too (pilot 4).
7. **The form of the input matters at 2 Hz only.** Pulses at 10 and 50 Hz act as DC of the same mean; at 2 Hz every pulse evokes one locked population event and more adaptation is recruited (pilot 5).
8. **In the paper's parameter space, the stored $\lambda_1$ and mean rate do not identify aberrant activity.** The high-rate networks are mostly excitatory populations pinned at maximum, and the one clear population burster looks stable and low-rate by its scalars. A parameter-space version needs an event measure computed inside each job (pilot 6).

## Decisions these pilots put to Tom

1. **Which network carries block B:** run13, the bursting network, or both as two kinds of aberrant activity? Both keep their aberrant activity under chained MTS and both are suppressed by tonic input (pilots 0a, 0b, 1a, 1b). Also: is `bursting_pairs` the bursting network you remember (pilot 0b could not confirm it)?
2. **Does product MTS appear in block B?** It has almost no baseline events on either network at the mu7chain constants. Pilot 0b's options: match the product's second pool to the chain, or show chained MTS only.
3. **Is clamping from stimulus onset the "adaptation ablated" arm?** Pilot 3 rules out removal and constants; pilot 2 shows clamping is exact. If yes, the clamped arm needs its own no-stimulus control, because clamping alone changes the dynamics (pilot 2).
4. **What "pre-stimulation value" means in an oscillating network:** the value at the instant of onset (which catches one phase of the burst cycle), the time average over a pre-stimulus window, or several onset phases (pilot 2).
5. **Is "steady elevated firing" an acceptable meaning of suppressed,** or must the stimulated state also be shown more stable, by a Lyapunov or perturbation measure (pilot 1b)?
6. **The measure of aberrant activity and "hyperexcitable at baseline":** adopt excess co-activation with its companion (pilot 4)? If so, accept that the paper's chained-MTS network counts as hyperexcitable at baseline, or add a criterion it does not measure, such as event amplitude.
7. **The form of tonic input:** DC only, or also pulses at a low rate? Pilot 5 calls 2 Hz the clinical stimulation rate (please confirm); at that rate the network follows every pulse, so a burst count alone cannot tell suppression from entrainment.
8. **The parameter-space version:** which state counts (saturated, chaotic, bursting), whether to sample a targeted region rather than the whole grid, and how many networks (pilot 6 gives 145 to about 1600 for 50 hyperexcitable networks, by definition).
9. **The random-generator defect** (see "Defects"): fix the model so it names its generator? Until then a network rebuilt in an ordinary session is not the one a parallel worker simulated (pilots 1a, 6).

# The pilots at a glance

| Pilot | Question | Answer | How solid |
|--|--------|----------|----|
| 0a | Does run13 still make IED-like events with chained MTS? | Yes, fewer and smaller (12 against 30 per 100 s); product MTS nearly silent (1). | 1 seed, 100 s |
| 0b | Does the bursting network still burst with chained MTS? | Yes, same rhythm, smaller bursts (34 against 33 per 100 s); product MTS does not burst. | 1 seed, 100 s, deterministic |
| 1a | Does stimulation suppress run13's events under chained MTS? | Yes, by 83 % at input 0.20, the same fraction as run13 as built. | 8 seeds, 150 s per arm, paired |
| 1b | Does stimulation stop the bursts under chained MTS? | Yes, at every amplitude tried, by steady elevated firing with adaptation recruited. | 1 seed, 60 s on, deterministic |
| 5 | DC against pulses of equal mean? | Same at 10 and 50 Hz; at 2 Hz one locked event per pulse and more recruitment. | 1 seed, 60 s on, deterministic |
| 2 | Can adaptation be clamped at its pre-stimulation value without changing the model? | Yes, exactly, by script. Clamped, the stimulus no longer quiets. | 1 seed per network, 1 onset |
| 3 | Is there a fair ablated baseline by removing adaptation? | No: removal or constants saturate or silence. Clamping is the only candidate. | 1 seed, 100 s, 33 runs |
| 4 | One measure of aberrant activity for both kinds? | Excess co-activation with mean rate and saturated fraction; constants chosen after looking. | 16 traces, 2 networks, 1 seed |
| 6 | How common is hyperexcitability in the paper's parameter space? | 10, 11 and 2 of 32 (chaotic, high-rate, both) with chained MTS; the scalars miss bursting. | 32 networks, 20-s runs |

# Pilot 0a: run13 with chained MTS

**Question.** Does run13 still produce brief, localized discharges when its single-timescale adaptation is replaced by chained MTS, and does the chained form behave like the product form?

**Method.**

- run13 as built (single-timescale) against chained MTS and product MTS on the same network, setpoints, initial state and noise path.
- The MTS conditions keep run13's SFA time constant as the fastest rung, its total SFA strength and its depression constants; the slower rungs (1.58 and 10 s), the reserve (refill 5 s) and the product's second pool come from mu7chain.
- 10 s settling, 100 s observed, seed 42; events by the IED pilot's fixed detector (a module recruited when more than 20 % of its neurons cross threshold).
- Regression: run13 as built, at its own length, reproduces the documented run13 exactly (11 events at the documented times).

**Result.**

- Chained MTS keeps the discharges, rarer and smaller but still local (median 8.0 % of neurons, one module each, none global; median width 0.23 s against 0.20 s).
- Product MTS nearly silences them. The product depresses excitatory transmission more than the chain, the likely reason.
- Mean rates fall with MTS, and within-module correlation falls (0.24, 0.17, 0.11).

| | run13 as built | chained MTS | product MTS |
|---|---|---|---|
| events in 100 s (per min) | 30 (18.0) | 12 (7.2) | 1 (0.6) |
| mean E / I rate | 0.125 / 0.138 | 0.091 / 0.095 | 0.072 / 0.072 |
| mean depression factor, E outputs | 0.607 | 0.585 | 0.557 |

![run13 under its as-built adaptation, chained MTS and product MTS (columns), 100 s. Look at the bottom row (module recruitment, triangles are detected events): many events as built, fewer under the chain, one under the product, and the depression factor (row 5) lowest under the product.](../../figs/explore_B/run13_chain_mts/run13_chain_mts.png){width=100%}

Full report: [figs/explore_B/run13_chain_mts/report.md](../../figs/explore_B/run13_chain_mts/report.md)

# Pilot 0b: the bursting network with chained MTS

**Question.** Does the bursting network still burst when depression is the chained vesicle pool, and does the chain behave like the product?

**Method.**

- The bursting network as built (three SFA timescales on E, one depression pool on E to E and E to I) against chained MTS and product MTS; the network's own constants kept, reserve and second pool from mu7chain.
- No input, no noise, so the bursts are spontaneous and the run deterministic; 10 s settling, 100 s, seed [42 42].
- A burst: population mean E rate above 0.2 for at least 0.1 s.

**Result.**

- As built and under chained MTS the network bursts at the same rhythm (every 3.0 and 2.95 s); under product MTS it does not burst at all, but settles into a small regular oscillation with the synapses held depressed.
- Chained MTS bursts are smaller and shorter: 70 % of neurons active at the peak against 82 %, peak mean E rate 0.39 against 0.59.
- So, as on run13, chain and product are not alike on a first try. Comparing them here would compare a bursting with a non-bursting network.

| | as built | chained MTS | product MTS |
|---|---|---|---|
| bursts in 100 s | 33 | 34 | 0 |
| peak mean E rate | 0.59 | 0.39 | 0.10 |
| mean E / I rate | 0.081 / 0.087 | 0.061 / 0.063 | 0.037 / 0.039 |
| depression factor range | 0.60 to 0.94 | 0.63 to 0.84 | 0.67 to 0.76 |

![The bursting network as built, under chained MTS and under product MTS (columns), 100 s. Rows: every neuron's rate, population mean rates with bursts shaded, mean adaptation term, mean depression factor. Look at the product column: no bursts, and a depression factor that never recovers above 0.76.](../../figs/explore_B/bursting_chain_mts/bursting_chain_mts.png){width=100%}

Full report: [figs/explore_B/bursting_chain_mts/report.md](../../figs/explore_B/bursting_chain_mts/report.md)

# Pilot 1a: stimulation of run13 under chained MTS

**Question.** Does the sustained uniform input that cuts run13's events still do so under chained MTS, does it recruit adaptation the same way, and is the effect as large?

**Method.**

- run13 as built and chained MTS (product MTS left out: about one event per 100 s). The existing experiment's stimulus: a constant uniform input to all 350 neurons at 0.20 (also 0.10 and 0.30), on from t = 0 after 10 s of settling.
- Paired within seed: stimulation off and on, 150 s each, same network, setpoints, initial state and noise path; the existing experiment's 8 seeds.
- Paired differences tested with an exact sign-flip test over all 256 sign patterns.
- Regression: the existing protocol re-run for run13 as built reproduces the documented per-seed table to $5\times10^{-15}$ relative.

**Result.**

- Yes, by the same fraction: at 0.20 the event rate falls by 0.83 under both, lower with stimulation in 8 of 8 seeds in each (p = 0.0078, the smallest possible with 8 seeds).
- Graded with amplitude under both forms (table).
- Adaptation is recruited in every seed, a little more under the chain: E SFA term +61 % against +51 %, depression factor -20 % against -14 %. Mean E rate rises by 62 % and 51 %.
- Remaining events are no larger. Holding the unstimulated threshold fixed, or dropping the first 10 s after onset (stimulus onset often evokes one large event), changes the reductions only slightly (0.80 and 0.87 under the chain).

| Mean over 8 seeds | run13 as built | chained MTS |
|---|---|---|
| events per min, off | 20.3 | 15.3 |
| events per min, on at 0.10 / 0.20 / 0.30 | 9.3 / 3.55 / 1.55 | 5.7 / 2.55 / 1.00 |
| relative reduction at 0.10 / 0.20 / 0.30 | 0.54 / 0.83 / 0.92 | 0.63 / 0.83 / 0.93 |
| E SFA term, off to on (0.20) | 0.062 to 0.094 | 0.043 to 0.069 |
| depression factor, off to on (0.20) | 0.616 to 0.527 | 0.603 to 0.485 |

Seed 42 is not typical: in pilot 0a it gave 30 and 12 events per 100 s (18.0 and 7.2 per min), while the 8-seed means are 20.3 and 15.3 per min. Seed 42's own off arm here, 150 s long, gives 18.8 and 8.8 per min; the noise path depends on the run's length, so a 150-s run does not repeat a 100-s run.

![Eight seeds, off against on at 0.20 (one line per seed) for run13 as built (left) and chained MTS (middle), events against amplitude (right); bottom row, the SFA term, depression factor and mean E rate. Look at every line falling in the top row while every line in the bottom row moves the same way: more SFA, more depression, higher rate.](../../figs/explore_B/p1a_stim_run13/p1a_summary.png){width=100%}

Full report: [figs/explore_B/p1a_stim_run13/report.md](../../figs/explore_B/p1a_stim_run13/report.md)

# Pilot 1b: stimulation of the bursting network under chained MTS

**Question.** Does a tonic stimulus stop the bursts under chained MTS as it does with the network as built, and is adaptation recruited, or is bursting simply replaced by high firing?

**Method.**

- The bursting network as built and under chained MTS; product MTS under the reference stimulus only.
- One run per case: 40 s without input (baseline), then the stimulus from 40 to 100 s; measured at 10 to 40 s and 50 to 100 s. Seed [42 42].
- Stimuli: the manuscript stimulation figure's (DC 0.4 to every neuron plus weak white noise); the same without noise; DC 0.2 and 0.1; DC 0.4 to 7 of the 35 E neurons.
- Beside the burst count: the largest population mean E rate and its peak-to-peak range, because the count misreads a steady rate at threshold (a rule was added: above threshold for more than 5 s is sustained firing, not a burst).

**Result.**

- Yes: under chained MTS every case stops the bursts. With the reference stimulus the largest mean E rate falls from 0.39 to 0.21 and its peak-to-peak range from 0.39 to 0.01. (The reference count reads 0.76 per s, but these are noise-driven crossings of a steady rate sitting at 0.200, the threshold.)
- Suppression is steady elevated firing, not quiet and not saturation: no neuron's mean rate exceeds 0.9.
- Adaptation is recruited and the effect is graded with amplitude (table).
- Partial stimulation separates the conditions: chained MTS goes steady (mean E rate 0.120); the network as built turns its bursts into a small 0.88-Hz oscillation. One choice of 7 neurons only.
- Every condition fires one large transient at onset, product MTS included (which had no bursts); product MTS then settles at 0.170.

| Chained MTS, before / during | mean E rate | largest mean E rate | adaptation term | depression factor |
|---|---|---|---|---|
| DC 0.4 | 0.060 / 0.200 | 0.386 / 0.205 | 0.030 / 0.099 | 0.752 / 0.473 |
| DC 0.2 | 0.060 / 0.133 | 0.386 / 0.136 | 0.030 / 0.066 | 0.752 / 0.560 |
| DC 0.1 | 0.060 / 0.097 | 0.386 / 0.098 | 0.030 / 0.048 | 0.752 / 0.631 |
| as built, DC 0.4 | 0.081 / 0.297 | 0.588 / 0.305 | 0.041 / 0.147 | 0.816 / 0.570 |

![Chained MTS, one column per stimulus (reference, DC 0.4, 0.2, 0.1, 0.4 to 20 % of E); onset at 40 s (red). Look at the bursts giving way to a flat, higher rate in every column, with the adaptation term stepping up and the depression factor stepping down. Grey stripes in the left column are threshold crossings of a steady rate, orange bands sustained firing, not bursts.](../../figs/explore_B/p1b_stim_bursting/p1b_chain_mts.png){width=100%}

Full report: [figs/explore_B/p1b_stim_bursting/report.md](../../figs/explore_B/p1b_stim_bursting/report.md)

# Pilot 5: DC against pulses of equal mean

**Question.** At the same time-averaged input, do DC and a pulse train have the same effect on the bursts, and on the adaptation they recruit, under chained MTS?

**Method.**

- The bursting network under chained MTS, the network as built as a reference; 40 s baseline, input on from 40 to 100 s on every neuron, in phase; seed [42 42], no noise.
- Mean input 0.1 and 0.2; DC, or rectangular pulses at 2, 10 and 50 Hz, 20 % duty, amplitude 5 times the mean. Means verified exact.
- Integrated at 2000 Hz (at the usual 400 Hz a 50-Hz pulse would be 1.6 samples wide).
- Pilot 5 chose 2 Hz as "the rate of clinical chronic subthreshold stimulation" (its statement, to be confirmed).

**Result.**

- Every form stops the network's own bursts. 10 and 50 Hz act as DC (mean rate within 2.5 % of DC, depression and adaptation within 0.1 %).
- At 2 Hz every pulse evokes exactly one population event, locked to the pulse (vector strength 1.00, onset 34 to 62 ms after the pulse starts), and nothing in between: entrainment, not the network's own rhythm. The evoked events are about as large as the spontaneous bursts at mean 0.1 and larger at 0.2 (peak mean E rate 0.62 against 0.39).
- 2 Hz recruits clearly more depression and adaptation than DC of the same mean (table); the network as built behaves the same way.
- A correction to the brief: the dendritic time constant (0.1 s) has its corner near 1.6 Hz, not 10 Hz, so 10 Hz is filtered too, only less than 50 Hz (pulse-locked swing of the mean E rate 0.13 at 10 Hz against 0.024 at 50 Hz, mean 0.2).

| Chained MTS, mean 0.2, during | bursts per s | mean E rate | depression factor | adaptation term |
|---|---|---|---|---|
| before (no input) | 0.33 | 0.060 | 0.752 | 0.030 |
| DC | 0 | 0.133 | 0.560 | 0.066 |
| 2 Hz | 2.00 (locked) | 0.166 | 0.459 | 0.082 |
| 10 Hz | 0 | 0.136 | 0.551 | 0.068 |
| 50 Hz | 0 | 0.133 | 0.560 | 0.066 |

Both mean levels already suppress fully, so the pilot does not test whether the forms differ near the suppression threshold. Pilot 1b's lowest amplitude (0.1) also suppressed fully, so neither pilot locates that threshold on this network.

![Chained MTS at mean 0.1, columns DC, 2, 10, 50 Hz; bottom row, 2 s of input (grey) with the mean E rate (orange). Look at the 2-Hz column: one event per pulse, a lower depression factor and a higher adaptation term than DC; 10 and 50 Hz are indistinguishable from DC.](../../figs/explore_B/p5_dc_vs_pulses/p5_chain_mts_mean0p1.png){width=100%}

Full report: [figs/explore_B/p5_dc_vs_pulses/report.md](../../figs/explore_B/p5_dc_vs_pulses/report.md)

# Pilot 2: clamping adaptation at its pre-stimulation value

**Question.** Can the arm Brian asked for (stimulation on, adaptation held at its value at stimulus onset) be run with the model as it stands, and can it be trusted?

**Method.**

- The model runs normally to the last sample before the input changes; from there it continues with its own equations, integrator and noise path, but with the rate of change of every adaptation variable (SFA, depression pools, reserve) set to zero. The clamped variables still act on rates and synapses; the dendritic state evolves.
- Two networks, one seed each: the paper's network under chained MTS (500 neurons, noise on; a step input, positive to 15 % of E, negative to 15 % of I, from 10 s); the bursting network as built (+0.4 to a random half of E from 20 s).
- Four arms each: intact or clamped, with or without the stimulus.
- Checks: agreement before onset, drift of the clamped variables, continuity at the restart, and whether the continuation with nothing clamped reproduces the model's own run.

**Result.**

- A script of a few lines is enough, for chained and product depression, with or without noise. Every check is exact: the continuation with nothing clamped reproduces the model's run with a largest difference of 0; the clamped variables drift by 0.
- Bursting network: intact, the stimulus stops the bursts and leaves steady firing; clamped, the network runs to saturation. This is the direction block B predicts.
- Paper's network: clamped, the stimulus raises the rate much more; free adaptation absorbs about two thirds of the increase.
- Clamping is not neutral by itself. Clamped without a stimulus, the paper's network keeps its mean rate but its rate fluctuates almost twice as much (SD 0.076 against 0.046); the bursting network, clamped between two bursts when adaptation was near its lowest, falls silent.
- What not to try: making every adaptation time constant infinite clamped nothing on the bursting network's synapses (depression drifted by 0.69, no error), because parameter names differ between synapse models.

| Mean E / I rate after onset | paper's network | bursting network |
|---|---|---|
| intact, stimulus | 0.293 / 0.253 | 0.359 / 0.223 |
| clamped, stimulus | 0.459 / 0.438 | 0.906 / 0.984 |
| intact, no stimulus | 0.202 / 0.211 | 0.082 / 0.088 |
| clamped, no stimulus | 0.210 / 0.231 | 0.007 / 0.016 |

Rates and traces only: Lyapunov exponents of the clamped network cannot be read from the model's estimators as they are.

![Bursting network, onset at 20 s. Left, adaptation intact: the stimulus (solid) stops the bursts as adaptation rises and the synapses depress. Right, clamped from onset: the stimulated network saturates near 1, and the unstimulated clamped network (dashed) falls silent.](../../figs/explore_B/p2_frozen_state/p2_frozen_bursting.png){width=90%}

![Paper's network under chained MTS, onset at 10 s, same layout. Look at the top row: the stimulated rate settles near 0.29 with adaptation intact and climbs to about 0.46 clamped, while the clamped adaptation and depression (rows 2 and 3) are flat.](../../figs/explore_B/p2_frozen_state/p2_frozen_mu7chain.png){width=90%}

Full report: [figs/explore_B/p2_frozen_state/report.md](../../figs/explore_B/p2_frozen_state/report.md)

# Pilot 3: removing adaptation, and a fair ablated baseline

**Question.** Unstimulated, what does each network do with adaptation removed, and does any way of removing it leave comparable aberrant activity?

**Method.**

- Both networks, chained MTS mainly, each as built too; 10 s settling, 100 s, seed 42, no input. Every ablation is a set of model arguments.
- SFA removed; depression removed; both removed; slow parts removed (fastest SFA timescale and first depression pool only).
- Replaced by constants ("mean-matched"): both removed, each depressed route's weights scaled by the intact run's average depression factor, and each neuron given a constant input of minus its time-averaged SFA term. The constant was taken as the plain time average, the rate-weighted average, four values between, and on run13 one per neuron.
- Measures: events or bursts per minute, mean rates, fraction of neurons pinned above 0.9 or silent below 0.01.

**Result.**

- No ablation that removes depression, and no constant, leaves a usable baseline. The only one with comparable aberrant activity is SFA removed on run13 (11 events per 100 s against 12 intact), which tests SFA alone and leaves depression, the larger part there, intact.
- Depression removed, or both removed: both networks saturate (bursting network mean E rate 0.83 to 0.92; run13 1.00, 99 to 100 % of neurons pinned).
- Constants: the plain time average of the depression factor (0.75 bursting, 0.58 run13) is set mostly by quiet periods when synapses have recovered, so it brakes too little and both networks saturate. The rate-weighted average silences the bursting network and leaves run13 firing steadily with no events. Between the two, the bursting network goes from silent (up to 0.71) to saturated (0.75) with no bursts anywhere.
- Slow parts removed with total SFA strength kept increases aberrant activity: the bursting network bursts twice as often (70 against 34 per 100 s), and run13 becomes exactly run13 as built (30 events against 12).
- The reason, as the pilot puts it: the events depend on adaptation changing in time (it ends each burst or event), so any constant removes either the events or the brake.

![Every unstimulated run: aberrant events per minute (left) and mean E rate (right), bursting network above, run13 below, coloured by regime (green comparable, dark red saturation). The dashed line is intact chained MTS. Look at how few bars are green.](../../figs/explore_B/p3_ablation/summary.png){width=100%}

![Bursting network under chained MTS: intact and five ablations (columns). Removing depression, or replacing adaptation by its averages, pins the network; bursts survive only when the slow parts alone are removed.](../../figs/explore_B/p3_ablation/bursting_ablations.png){width=100%}

Full report: [figs/explore_B/p3_ablation/report.md](../../figs/explore_B/p3_ablation/report.md)

# Pilot 4: one measure of aberrant activity

**Question.** Is there one measure, computed from the per-neuron rates alone (no modules, no per-network thresholds), that agrees with each network's detector, falls under stimulation, and does not mistake a saturated or silent network for a healthy one? What should "hyperexcitable at baseline" mean in its terms?

**Method.**

- 16 traces: the six unstimulated ones of pilots 0a and 0b; the four arms of pilot 2 on the bursting network; run13 as built and chained MTS without and with input 0.20 (50 s each); the paper's network without adaptation and under chained MTS.
- Ten candidates: four from the population mean rate, five from co-activation (the fraction of neurons above their own baseline at once), and synchrony.
- Judged on ordering against the detectors, falling under stimulation, readings for saturated and silent arms, stability across window halves and constants.
- The constants of the recommended measure were chosen once, by stated reasoning, after looking at these traces.

**Result.**

- **Recommended: excess co-activation A,** the time mean of the fraction of neurons simultaneously above their own threshold, in excess of the level independent neurons would reach by chance. Read with a companion: the mean rate and the fraction of neurons pinned above 0.9. A alone reads 0 for saturated and silent networks.
- It orders both networks as their detectors do (table) and falls under stimulation: bursting network 47.9 to 0; run13 as built 1.32 to 0.13; run13 chained MTS 0.54 to 0.02 (thousandths).
- A four-way class (silent if mean rate below 0.02; saturated if a quarter of neurons are pinned; aberrant if A is at least 0.25 thousandths; otherwise healthy) classes all 16 traces as described, with two readings to note. The bursting network under product MTS, which the burst detector calls quiet, is aberrant by A: a regular low-amplitude network-wide oscillation (A measures co-activation, not amplitude). And the paper's chained-MTS network reads A = 1.85, above run13 as built: it has brief network-wide events of its own.
- Proposed "hyperexcitable at baseline": unstimulated, at least 50 s after settling, classed aberrant. The threshold sits between run13 product MTS (0.10) and chained MTS (0.54).
- The other candidates each fail: mean-rate measures cannot see run13's events (about 8 % of neurons); synchrony reads high for a stimulated network whose bursts stopped; the burst detector reads a steady rate above 0.2 as one unbroken burst.

| A (thousandths), detector count | as built | chained MTS | product MTS |
|---|---|---|---|
| bursting network | 47.5, 33 | 37.8, 34 | 3.1, 0 |
| run13 | 1.34, 30 | 0.54, 12 | 0.10, 1 |

The ordering holds for every constant tried, but the value does not: on run13, A changes two- to threefold per unit of its constant k, so the threshold belongs to k = 3. Between window halves A changes under 1 % on the bursting network, about 30 % on run13 and fivefold on the paper's network. A cannot see aberrant activity that is the network's majority state.

![The six unstimulated traces (left bursting network, right run13; rows as built, chained MTS, product MTS). Each pair: mean E rate with the detector's events (triangles), and co-activation with the chance level dashed. Look at whether co-activation peaks above the dashed line sit under the triangles.](../../figs/explore_B/p4_aberrant_measure/p4_recommended_unstimulated.png){width=100%}

![The four arms of pilot 2 on the bursting network. Only the unstimulated intact arm has co-activation; the steady, saturated and silent arms all read A = 0 and are told apart by the mean rate and pinned fraction in the titles.](../../figs/explore_B/p4_aberrant_measure/p4_recommended_bursting_arms.png){width=100%}

Full report: [figs/explore_B/p4_aberrant_measure/report.md](../../figs/explore_B/p4_aberrant_measure/report.md)

# Pilot 6: hyperexcitability in the paper's parameter space

**Question.** Of the networks sampled across the paper's connectivity parameter space, how many are hyperexcitable at baseline, by which definition, and how many stay so with adaptation on? That sets whether a parameter-space stimulation experiment is worth building and how large it must be.

**Method.**

- The joint sample of the mu7chain essential run, already on disk: 32 networks, four conditions (none, single-timescale, chained MTS, product MTS), seven axes at three levels; 20-s runs, Lyapunov window 10 to 20 s. Only scalars are stored, no time series.
- Definitions: chaotic ($\lambda_1 > 0$); high rate (mean rate above 0.40, a gap in the distribution); both.
- Six networks rebuilt and simulated for 30 s, after checking each against its stored values. 95 % Wilson intervals; the 1-D sweeps (112 networks per condition) counted as a cross-check.

**Result.**

- With chained MTS, 10 of 32 are chaotic (interval 18 to 49 %), 11 of 32 high-rate (20 to 52 %), 2 of 32 both (2 to 20 %); without adaptation 10, 12 and 2. The answer depends strongly on the definition.
- Across this wide sample adaptation moves almost no network across a threshold. Near the paper's network it does: on the 1-D sweeps, chaotic and high-rate falls from 26 of 112 without adaptation to 4 of 112 with chained MTS.
- The stored scalars do not identify aberrant activity. The two "chaotic and high-rate" networks have every E neuron pinned at maximum, with the chaos in a sparse I population; network 761, "stable and low-rate" under chained MTS, is a regular population burster (a full burst about every 4 s), the closest thing to a paroxysmal discharge in the sample.
- High rate follows the excitatory fraction and the realised E:I weight balance; chaos follows the gain on all weights.
- Sample size for 50 hyperexcitable networks with chained MTS: about 145 (high rate), 160 (chaotic), 800 (both), about 1600 (stable-looking bursters like 761, 1 of 32), with wide ranges.

"Chaotic" is fragile here: stable networks with adaptation sit at $\lambda_1$ near -0.1 per s, and with the threshold at -0.1 the chained-MTS count jumps from 10 to 18. One rebuilt network was still drifting at 20 s.

![The 32 networks under each condition, $\lambda_1$ against mean rate; red fill marks chaotic and high-rate. Look at how little the classification changes across panels, and at the numbered networks that were rebuilt.](../../figs/explore_B/p6_hyperexcitability/p6_lle_vs_rate.png){width=100%}

![Six rebuilt networks, 30 s each, without adaptation (left pair) and with chained MTS (right pair). Look at the bottom right: network 761, stable and low-rate by its stored scalars, bursts every 4 s; in the top two rows the E population is saturated in both conditions.](../../figs/explore_B/p6_hyperexcitability/p6_timeseries_examples.png){width=100%}

Full report: [figs/explore_B/p6_hyperexcitability/report.md](../../figs/explore_B/p6_hyperexcitability/report.md)

# What the pilots mean together for block B

## Which network carries the result

- **Supported:** both networks keep their aberrant activity under chained MTS without retuning (0a, 0b), and tonic input suppresses it under chained MTS as it does with each network as built (1a, 1b). So either can carry block B, and both together give the two kinds of aberrant activity the to-do list asks for. They differ in what "as built" means: run13 as built is single-timescale, so on run13 the chained-MTS comparison is also a 1TS against MTS comparison; the bursting network as built already has three SFA timescales, so there it isolates the depression form.
- **Ruled out, at the mu7chain constants:** product MTS as a third condition on either network. It has 1 event per 100 s on run13 and no bursts on the bursting network, so a stimulation effect cannot be measured in it. Whether that is the product form or the second pool's constants (release 0.5 s) is not known; 0b suspects the constants.
- **Unknown:** the bursting network has been run on one seed and is strictly periodic; nothing yet shows it bursts over other seeds. The paper's own network was not a stimulation substrate in these pilots, though pilot 4 finds it has brief network-wide events of its own.

## The causal experiment and the "adaptation ablated" arm

- **Ruled out** (pilot 3): removal of depression (saturation), removal of both (saturation), adaptation replaced by any constant (saturation or silence, with no comparable events at any constant tried), and removal of the slow parts (more aberrant activity, not less). SFA removed on run13 keeps a comparable baseline but is a partial ablation of the smaller part.
- **Supported** (pilot 2): clamping from stimulus onset, which by construction leaves the unstimulated baseline untouched, is exact and cheap. It is the only fair arm found. Pilot 3 warned that a network clamped between events holds recovered synapses and may run away once driven; pilot 2 saw exactly that on the bursting network.
- So the design becomes intact against clamped, each with and without the stimulus. The clamped no-stimulus cell is needed because clamping alone changes the dynamics (silence on the bursting network; doubled rate fluctuation on the paper's network).
- **First evidence for causality, thin:** on the bursting network the stimulus quiets the network only when adaptation is free (pilot 2). One seed, one onset, and a stimulus (half of E at 0.4) different from pilot 1b's.
- **Unknown:** the result over seeds and onset phases; how "pre-stimulation value" is defined in an oscillating network (decision 4); whether clamped saturation should count as "not suppressed" or as a different failure (pilot 4's class would call it saturated, not aberrant).

## The measure of aberrant activity and "hyperexcitable at baseline"

- **Supported:** excess co-activation with its companion (pilot 4) reads both kinds of aberrant activity, falls under stimulation on both networks, and its class separates the arms of the causal experiment: the stimulated intact arm reads healthy, the clamped stimulated arm saturated, the clamped unstimulated arm silent.
- **Bears on pilot 1b:** suppression in every pilot is steady elevated firing. The burst count misreads such a state when it sits near threshold (1b, 5); A reads it as 0 and the companion reports the rise in rate. Pilot 5 adds a requirement for pulsed input: the measure must say whether events are locked to the pulses. Pilot 4 did not test A on pulsed input.
- **Bears on pilot 6:** the stored $\lambda_1$ and mean rate cannot find aberrant activity, and A is designed to run inside a sweep job (about 0.1 s for 500 neurons over 100 s, nothing stored). Pilot 6's definitions (chaotic, high-rate) and pilot 4's (aberrant by A) are not the same thing and have not been compared on the same networks.
- **Unknown or uncomfortable:** the threshold rests on one IED network with a margin of about a factor of 2; the constants were chosen after looking. By the proposed threshold the paper's own chained-MTS network is hyperexcitable at baseline (decision 6). And A calls product MTS on the bursting network aberrant where the burst detector calls it quiet; the two disagree because A measures co-activation, not amplitude.

## The form of the tonic input

- **Supported:** DC is a fair stand-in for pulses at 10 Hz and above (pilot 5).
- **At 2 Hz the choice is scientific:** the network follows every pulse and recruits more depression and adaptation than DC. If the paper shows a pulse train at a clinical rate, the result to report is "own rhythm replaced by pulse-locked responses", and the measure has to say so.
- **Unknown:** behaviour near the suppression threshold, which neither pilot 1b nor pilot 5 reached on the bursting network (0.1 already suppresses fully); pulses on run13; a seizure-like drive.

## The parameter-space version

- **Supported as feasible** (pilot 6), but not on the stored scalars. It needs an event measure inside each job (pilot 4's A is the candidate), stored beside $\lambda_1$ and mean rate.
- **Sampling:** a uniform sample of the whole grid is wasteful for the rare, relevant classes. Most high-rate networks are excitatory populations locked at maximum, probably not what stimulation is meant to treat; the one burster was at 50 % excitatory and low gain. Pilot 6 suggests sampling the region where the chosen state occurs (for example 50 % excitatory, gain 1 or above, E:I weight balance near 0.5).
- **Unknown:** how many networks are aberrant by A at baseline, and whether the 20-s essential runs are long enough (one network still drifting at 20 s).
- **Must be fixed first:** a sampled network can be rebuilt exactly only with the same random generator the parallel workers used (see Defects).

## Chain and product MTS

The two depression forms were expected to behave alike; on both networks they did not (0a, 0b): the product depresses more and nearly silences the aberrant activity. This matters beyond block B, since every analysis runs both forms. Pilot 6 shows the same split by scale: across the wide joint sample the counts change little between chained and product MTS (chaotic 10 and 8, high-rate 11 and 11 of 32), but on the 1-D sweeps near the paper's network they differ widely (chaotic 43 and 22, high-rate 32 and 8 of 112). Whether this is the product form or the second pool's constants is not resolved.

# Defects and loose ends found on the way (not fixed)

- **The model seeds its random numbers without naming the generator** (pilots 6 and 1a, independently). Parallel workers default to Threefry, an ordinary session to the Mersenne twister, so the same seeds build a different network and noise path. Pilot 6: network 761 rebuilt plainly gave $\lambda_1$ +0.04 against -0.10 stored; with the generator set, all 12 runs matched to 1e-15. Pilot 1a: a pool run gave 13 events for seed 42 where 16 are documented. Consequences: exact rebuilds of sweep networks need the generator set; a serial and a pool sweep simulate different networks; the E:I weight-balance colouring of a sweep figure is computed on a different draw than was simulated.
- **The essential run's 1-D sweeps ran 16 jobs per condition, not 8** (pilot 6): the reps axis got four levels between 1 and 2 (1, 1.33, 1.67, 2), each its own seed, where 4 levels by 2 reps were meant. This is why those sweeps took twice the estimated time.
- **The burst-count rule misreads a steady rate near its threshold** (pilots 1b and 5): a rate at 0.200 becomes noise-driven "bursts", a rate grazing or just above 0.2 becomes one burst spanning the whole input period; a saturated network counts as one unending burst (pilots 3, 4).
- **The header of the manuscript's stimulation figure** (`fig_stim_engages_adaptation`) describes the bursting network the wrong way round: it bursts with no input and goes steady under the DC step (pilots 0b, 1b).
- **The IED code's default input no longer builds** (`ied_build_network`; pilot 0a).
- **The IED stimulation code's "SFA feedback"** (`ied_stim_experiment`) lacks the division by the number of SFA timescales: right for run13 as built, three times too high under MTS (pilot 1a).
- The SRA1 noise path depends on the run's length, so a longer run does not repeat a shorter one's opening (pilot 0a).
- Making adaptation time constants infinite does not clamp product-route synapses, silently (pilot 2).
- `generate_stimulus` pads the settling interval with zeros, so an intrinsic drive is off before t = 0 (pilot 3); `build_from_preset` leaves the class's default stimulus on a preset with no input configuration (pilot 0b); the linear input interpolant spreads a step at t = 0 over the last settling step (pilot 1a).
- Stimulus onset produces one large, widely recruited response that the IED detector counts as an event (pilots 1a, 1b).
- The run13 data folder named in its README (`data/IED_model/20260930/run13`) is not on this machine (pilot 0a); a stale comment in the `bursting_pairs` preset gives the SFA strength as 0.5/3 where the value is 0.5 (pilot 0b).

# Where things are

All figure folders are under `figs/explore_B/` (local, not in git); each holds `report.md`, `report.pdf`, the figures and the data. Scripts are under `scripts/explorations/block_B_pilots/` unless stated.

| Pilot | Figures and report | Scripts | Rerun time |
|-|-----|--------|-----|
| 0a | `run13_chain_mts/` | `scripts/examples/run13_chain_mts_example/` | about 75 s, serial |
| 0b | `bursting_chain_mts/` | `scripts/examples/bursting_chain_mts_example/` | about 19 s, serial |
| 1a | `p1a_stim_run13/` | `p1a_stim_run13/` | about 8.5 min, pool of 4 |
| 1b | `p1b_stim_bursting/` | `p1b_stim_bursting/` | about 3 min, serial |
| 2 | `p2_frozen_state/` | `p2_frozen_state/` | 102 s, serial |
| 3 | `p3_ablation/` | `p3_ablation/` | about 8 min, serial |
| 4 | `p4_aberrant_measure/` | `p4_aberrant_measure/` (the measure: `coactivation_excess.m`) | about 65 s, serial |
| 5 | `p5_dc_vs_pulses/` | `p5_dc_vs_pulses/` | about 12.5 min, serial |
| 6 | `p6_hyperexcitability/` | `p6_hyperexcitability/` | 20 s table; about 5 min rebuild check; 1 min figures |

The condition builders the pilots share are `run13_condition_model.m` and `bursting_condition_model.m`, in the two example folders.
