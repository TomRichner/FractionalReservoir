---
title: "Block C: memory capacity across gain and bias"
subtitle: "A quick test on the chained-MTS network, four conditions, two seeds"
---

# Summary

**The question** (to-do list, block C; Tom): is the memory-capacity result
a property of adaptation, or of one well-tuned network? And is chained
multiple-timescale (MTS) adaptation the most robust of the conditions,
with consistently higher memory capacity than single-timescale (1TS)
adaptation or none?

**What was run.** Memory capacity and the largest Lyapunov exponent
$\lambda_1$ on the block A network, over a grid of 11 values of the overall
connection gain (0.1 to 4.0) by 7 values of the nonlinearity's bias $S_c$
(-0.6 to 0.6), for the four conditions, on two network seeds: 154 networks,
616 network-condition pairs, none failed. This is a quick test (Tom,
2026-11-13); more seeds and a finer grid are the stage's `medium` and
`production` sizes, for later.

**What it shows.**

1. **Adaptation widens the range of gain over which the network has memory,
   about fourfold.** At the preset's bias, memory capacity stays at or above
   half the best any condition reaches (7.5) over 0.55 octaves of gain
   without adaptation, and over 2.0 to 2.4 octaves with any of the three
   kinds of adaptation. Over the whole map, 9 % of the networks reach that
   level without adaptation and 22 to 27 % with adaptation. Without
   adaptation memory capacity is about 0.1 for every gain above 0.5,
   at every bias.
2. **In every condition memory is best just on the stable side of
   $\lambda_1 = 0$.** The maps show a ridge of high memory capacity that
   follows the edge of stability. What adaptation changes is where that edge
   lies: without adaptation it is at gain 0.33 to 0.5 at every bias; with
   adaptation it moves to higher gain, and further the lower the bias, so the
   stable region with good memory is several times larger.
3. **Chained MTS and 1TS are about level.** Their peaks are the same (15.0
   and 15.0); chained MTS is slightly wider along gain at the preset bias
   (2.40 against 2.03 octaves) and reaches half the overall best on a
   slightly larger share of the map (27 % against 23 %); grid point by grid
   point it is higher in both seeds at 40 % of points and lower in both at
   22 %, with a median difference of +0.03. **This test does not support
   "chained MTS is the most robust" as a claim over 1TS.** It supports
   "adaptation, of any of these kinds, makes memory robust to gain and
   bias".
4. **Product MTS** (the manuscript's present form) has a lower peak (12.8)
   but the flattest map (it reaches half its own best on 39 % of the map,
   the most of any condition) and the stable region reaching furthest in
   gain: at bias -0.6 it is still stable at gain 4.0, the end of the grid.
   Its two seeds agree least (r = 0.61).
5. **The preset network sits on the edge for 1TS and chained MTS.** At gain
   1, bias 0.2, $\lambda_1$ is within 0.1 of zero for both, and the two
   seeds differ widely (1TS: memory capacity 8.0 and 1.3, $\lambda_1$ of
   opposite sign; chained MTS: 10.2 and 5.8). Product MTS is safely stable
   there and has the highest memory capacity of the four at that one point
   (13.4 and 12.2). **So a comparison of conditions at the preset network
   alone measures, in part, how far each condition's edge of stability lies
   from that point.** This bears directly on the manuscript's Figure 6,
   which compares at one network.

What this means for the manuscript is left to Tom; a proposed Results
paragraph and caption that state findings 1, 2 and 5, and say plainly that
chained MTS and 1TS are not separated, are in the manuscript repo,
`docs/Block_C_proposed_results.md`.

# 1. Question

From the to-do list, block C: the memory result in the current draft is
measured at one network, so a reader can take it as that network having
been tuned. The claim to establish is that over a range of connectivities
memory capacity has a higher mean and a lower variance with MTS adaptation
than with 1TS or none. Brian's "expanded" experiment asks the same thing
in terms of recurrent gain. Tom added a second axis, the nonlinearity's
bias, and the question whether chained MTS is the most robust.

# 2. Method

**Network.** The `mu7chain` network of block A with noise off, preset

```
celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0_chainStd_4cond_mu7chain
```

500 neurons, half excitatory, in-degree 100, mean weights 7 and -7 on the
four routes (in units of the random-matrix scale), a per-neuron spread of
0.25 in the adaptation and the depression constants. Four conditions, built
on the same seeds in every job, so that at each grid point they share the
connection matrix, the input weights and the input sequence and differ only
in their adaptation:

| Condition | In the code | Adaptation |
|---|---|---|
| No adaptation | `no_adaptation` | none |
| Single timescale (1TS) | `sfa1_std1` | one SFA timescale, one depression pool |
| Chained MTS | `sfa3_stdc2` | three SFA timescales; releasable pool and reserve in a chain |
| Product MTS | `sfa3_std2` | three SFA timescales; two independent pools multiplied |

**The two axes.** *Gain* is `level_of_chaos`, a multiplier on every
connection weight (the preset has 1). *Bias* is $S_c$, the centre of the
per-neuron draw of the nonlinearity's setpoint (the manuscript's $a_0$; the
preset has 0.2, with a per-neuron standard deviation of 0.1, which is
kept). A lower bias makes a neuron fire at a lower input.

**Memory capacity** by the task and scoring of the manuscript's figure,
unchanged (`mc_protocol`): a scalar input of independent values held for
0.3 s each; 10 s washout, 600 s training, 150 s test; a ridge readout from
each neuron's synaptic output reconstructs the input 1, 2, ... holds back,
to 15 s; memory capacity is the sum of the test $R^2$ over delays.
Integrator `sra1` at 200 Hz. The sweep integrates the task in segments and
keeps only the readout samples, which gives the same $R^2$ to the bit as the
existing path (`test_mc_lean_path`) at about a hundredth of the memory.

**$\lambda_1$** on the same network, by a separate 30 s run with the first
30 s of the same input, accumulated over the last 20 s (top-K estimator,
K = 5): the exponent while the task's input is on. The value with no input
is stored beside it and not used here.

**Also stored:** the mean rates of the excitatory and inhibitory neurons and
the fractions of neurons with mean rate above 0.9 (pinned) and below 0.01
(silent).

**Grid and seeds.** Gains 0.1, 0.15, 0.22, 0.33, 0.5, 0.7, 1.0, 1.4, 2.0,
2.8, 4.0; biases -0.6 to 0.6 in steps of 0.2; network seeds 3001 and 3002
(the input sequence's seed follows the network's; the input weights are
the same in every job, as in the manuscript's figure). The grid was set
from a pilot (gains 0.1 to 3, biases -0.2 to 0.8, one seed, two conditions),
which put the no-adaptation peak on its low-bias edge and found bias 0.8
silent; gains 2.8 and 4.0 were added after the first seed's map showed the
adapting conditions still stable at gain 2.0 at low bias (FD-010, decisions
14 and 18).

**Statistics.** Per grid point, the median over the two seeds. "Half the
best of any condition" is 7.49 (the best median cell is 14.98, chained
MTS). A width along gain is the extent, in octaves, of the stretch of gains
at one bias where the median memory capacity is at or above that level; a
width along bias likewise at the best gain. A width that reaches the edge of
the grid is a lower bound, marked "open". A Gaussian fit against log gain
was made too; it fits poorly where the curve is flat and then collapses, so
the measured widths are reported.

**Code and data.** Stage `run_mc_gain_sweep`; figures by
`make_block_C_figures`. Result (local, not in git):

```
data/mc_gain_sweep_main/mc_gain_sweep/mc_gain_sweep_results.mat
```

Figures and tables: `figs/block_C/main/three_conditions/` and
`four_conditions/`. Decisions: FD-010 in `docs/ISSUES.md`.

# 3. Results

## 3.1 The maps

![Memory capacity (top), $\lambda_1$ with the task's input on (middle) and
the spread over the two seeds (bottom), gain by bias, four conditions.
Green dashed: the edge between cells with $\lambda_1 > 0$ and
$\lambda_1 \le 0$. Star: the preset network (gain 1, bias 0.2). Memory
capacity on a square-root colour scale shared by all panels; $\lambda_1$ on
a symmetric-log scale clipped at $\pm 10$ per second.](../../figs/block_C/main/four_conditions/Fig_MC_gain_bias_maps.png)

Read the top row against the middle one. In every condition the darkest
cells (highest memory capacity) lie along the green edge, on its stable
side. To the right of the edge memory capacity falls to about 0.1 (no
adaptation) or to below 1 (with adaptation) within one or two grid steps.
To the left it falls gently: at the lowest gain, 0.1, it is 3 to 5 without
adaptation and 5 to 7 with it.

**Where $\lambda_1$ turns positive**, gains between which the median over
seeds changes sign:

| Bias | None | 1TS | Chained MTS | Product MTS |
|---|---|---|---|---|
| -0.6 | 0.33–0.5 | 1.4–2 | 2–2.8 | stable to 4.0 |
| -0.4 | 0.33–0.5 | 1.4–2 | 2–2.8 | 2.8–4 |
| -0.2 | 0.33–0.5 | 1–1.4 | 1.4–2 | 2–2.8 |
| 0 | 0.33–0.5 | 1–1.4 | 1–1.4 | 2–2.8 |
| 0.2 | 0.33–0.5 | 0.7–1 | 0.7–1 | 1.4–2 |
| 0.4 | 0.33–0.5 | 0.7–1 | 0.7–1 | 0.7–1 |
| 0.6 | 0.5–0.7 | 0.7–1 | 0.7–1 | 0.7–1 |

Without adaptation the edge does not move with bias. With adaptation it
moves to higher gain as the bias falls: over the range of bias by about a
factor of two for 1TS, nearly three for chained MTS, and more than four for
product MTS, whose edge lies furthest out at low bias.
Each column changes sign once. The two seeds disagree on the sign of
$\lambda_1$ at 1 cell without adaptation, 3 with 1TS, 1 with chained MTS and
6 with product MTS, all next to the edge.

**On the stable side, $\lambda_1$ is pinned with MTS.** With chained and
product MTS, $\lambda_1$ is about -0.10 per second over the whole stable
region whatever the gain and bias (the light blue plateau): the slowest
adaptation timescale, not the distance to instability. With 1TS it is
about -0.3 to -1, and without adaptation it runs from -10 to near 0 as the gain
rises. As the to-do list anticipated, for MTS $\lambda_1$ says little about
how close a network is to the edge.

## 3.2 Per condition

Medians over the two seeds. "Half best" is 7.49, half the best median cell
of any condition.

| | None | 1TS | Chained MTS | Product MTS |
|---|---|---|---|---|
| Peak memory capacity | 10.2 | 15.0 | 15.0 | 12.8 |
| ... at (gain, bias) | (0.33, -0.2) | (1, -0.2) | (1.4, -0.2) | (1, 0.2) |
| Median over the map | 0.12 | 5.4 | 6.2 | 5.6 |
| Interquartile range over the map | 4.40 | 7.39 | 7.43 | 2.72 |
| 5th percentile over the map | 0.08 | 0.09 | 0.10 | 0.11 |
| Share of map at or above half its own best | 21 % | 23 % | 27 % | 39 % |
| Share of map at or above half best (7.49) | 9 % | 23 % | 27 % | 22 % |
| Width along gain at bias 0.2 (octaves) | 0.55 | 2.03 | 2.40 | 2.17 |
| Widest width along gain (octaves, at bias) | 0.77 (0) | 2.19 (0.4) | 2.40 (0.2) | 2.17 (0.2) |
| Width along bias at the best gain | 0.76, open | 0.72, open | 0.56, open | 0.44 |

The interquartile ranges are large for no adaptation, 1TS and chained MTS
because each map holds two populations, the stable side and the unstable
side; they measure where the edge falls on the grid, not scatter within
either side. The share of the map above a level, and the width along gain,
are the robustness measures.

![Left: memory capacity of every network (both seeds, every grid point)
per condition. Right: the share of the map at or above each memory
capacity.](../../figs/block_C/main/four_conditions/Fig_MC_robustness.png)

The right panel reads directly: without adaptation fewer than 40 % of the
networks have any memory at all (above about 0.5); with adaptation 65 to
85 % have. Above a memory capacity of about 8 the three adapting conditions
are close, and none separates from the others.

## 3.3 Paired comparisons

At each of the 77 grid points the conditions share the network, so they can
be compared point by point:

| | Higher in both seeds | Lower in both seeds | Median difference |
|---|---|---|---|
| Chained MTS minus none | 74 % | 10 % | +1.10 |
| Chained MTS minus 1TS | 40 % | 22 % | +0.03 |
| Chained MTS minus product MTS | 62 % | 26 % | +0.78 |
| 1TS minus none | 62 % | 14 % | +0.46 |

Chained MTS against 1TS is the one comparison that does not come out
clearly: higher in both seeds at 31 of 77 points, lower in both at 17, and
mixed at the remaining 29, with a median difference near zero.

## 3.4 The preset network

| | None | 1TS | Chained MTS | Product MTS |
|---|---|---|---|---|
| Memory capacity, seed 3001 / 3002 | 0.10 / 0.11 | 7.97 / 1.27 | 10.16 / 5.83 | 13.35 / 12.15 |
| $\lambda_1$ (1/s), seed 3001 / 3002 | 3.32 / 3.60 | -0.02 / +0.38 | +0.06 / +0.07 | -0.10 / -0.10 |

At the preset network no adaptation is chaotic; 1TS and chained MTS sit on
the edge, where a change of seed moves memory capacity by a factor of two
(chained) to six (1TS); product MTS is stable and best. Each of these
networks is one grid cell; on the map the same conditions reach 15 a step or
two to the side. The ranking of conditions at the preset network is
therefore a ranking of how far each one's edge of stability is from that
point, as much as of what each can do.

## 3.5 Seeds agree

| | Pearson r over grid points | Spearman | Largest difference at one point |
|---|---|---|---|
| None | 0.97 | 0.67 | 4.1 at (0.5, 0.6) |
| 1TS | 0.90 | 0.92 | 12.1 at (1, 0) |
| Chained MTS | 0.94 | 0.90 | 8.7 at (2, -0.4) |
| Product MTS | 0.61 | 0.77 | 12.4 at (2.8, 0) |

The large differences are all at cells on the edge, where one seed is
stable and the other not (bottom row of the maps). Away from the edge the
two seeds agree closely, which is why two seeds are enough for this test.

## 3.6 Silent and pinned networks

No cell reaches 90 % of neurons silent or pinned, so the low memory
capacities are not from a dead or saturated network. The most silent cells
(gain 0.1, bias 0.6: 73 to 74 % of neurons silent) still have memory capacity
3 to 5; the most pinned (no adaptation at gain 0.1, bias -0.6: 65 %) has
3.2. The memory capacity of about 0.1 to the right of the edge is chaos,
not silence.

## 3.7 The excitation- and inhibition-dominated networks

_(The two gain-only runs, mean weights 8.4 : 5.6 and 5.6 : 8.4 at bias 0.2,
one seed, are running; this section is filled in when they end.)_

# 4. Caveats

- **Two seeds.** Enough to see the structure, which is consistent across
  the grid; not enough for confidence intervals on the widths or for the
  paired comparison between chained MTS and 1TS, which is close. The
  `medium` size of the stage (16 gains by 9 biases by 5 seeds, about 12 h on
  12 workers) is the run that would settle it.
- **The grid is coarse along gain** (steps of a factor of about 1.4), so a
  width along gain is known to within about half an octave at each end. The
  ridge of high memory is one or two cells wide, so its height on this grid
  depends on where the cells fall relative to the edge.
- **Some widths are open.** The widths along bias reach bias -0.6 for three
  conditions, and product MTS's stable region reaches gain 4.0 at bias -0.6.
  Neither changes the comparison along gain at the preset's bias.
- **Product MTS has a lower peak but a farther edge.** Its two pools
  multiplied give stronger depression at a given rate than the chain does
  (block A, "Strength: settled, the two conditions are not matched at steady
  state"), which plausibly moves its edge out; this test does not separate
  that from the timescale structure.
- **The 1TS condition** uses, since F3 (FD-010), the fastest rung of the
  neuron's MTS adaptation ladder; its depression is the one-pool route of
  block A. Neither is matched in strength to the MTS conditions at steady
  state, as block A decided.
- **The input weights are the same in every job**, as in the manuscript's
  figure; only the network and the input sequence vary with the seed.
- **Noise is off**, as the design has it.

# 5. What is not done

- The two E:I networks (section 3.7) are running.
- More seeds and a finer grid (`medium`, `production`) are for later, Tom's
  to start.
- Whether the manuscript's Figure 6 changes is Tom's decision; the proposal
  is in the manuscript repo.
