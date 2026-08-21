# User Notes

Notes TR asks the agent to write down so they are not forgotten, and so he can
come back to them later if he decides to. **Nothing here is a commitment.** An
entry is a reminder, not a task: it is not a backlog item, it is not approval to
act, and it does not become work until TR says so.

This file is for TR to read. The agent writes an entry only when asked to, and
otherwise leaves it alone. Commit messages remain the primary record of what was
actually done.

Each entry is stamped with the date, branch @ commit, hostname, and the agent
and session that wrote it. Newest first.

---

## Rebuild the single-neuron STF methods figure (deliberately outside the refactorRunAll work)

| | |
|---|---|
| Noted | 2026-08-21 · `refactorRunAll` @ `22a91ee` · R5611351 · Claude Code (Opus 5), session e22d2fab |

TR wants this figure back eventually, but explicitly **parked it outside the
`refactorRunAll` refactor**. Recorded here so the investigation does not have to be
repeated. **Not started.**

**The orphan files.** Three files sit in
`scripts/presentations/Stability_Manuscript/fig_adaptation_methods/panel_A/` with no
script that names them:

```
sfa_std_stf_single_neuron_example_figure_1.png
sfa_std_stf_single_neuron_example_figure_1.svg
sfa_std_stf_single_neuron_example_f_1.fig
```

They are **not** the paper figure. The paper uses the 4-column SFA/STD figure
`sfa_std_single_neuron_example_figure_5.png`, produced by the surviving
`test_single_neuron_adaptation.m` in the same folder.

**What produced them.** `test_single_neuron_stf.m`, formerly at
`scripts/presentations/Stability_Manuscript/fig_adaptation_methods/test_single_neuron_stf.m`
(note: one level up from `panel_A/`).

| Commit | What it did |
|---|---|
| `390d86a` | "Stability_Manuscript: add single-neuron SFA/STD/STF methods figure" — added the script, 201 lines |
| `53f2dfd` | "use hard sigmoid (S_c=0.5) in STF single-neuron figure" |
| `60c2992` | "SRNNModelCellTypes: ragged per-type multi-timescale SFA" — **last version**; recover the file from here |

It was then deleted in the `refactor` cleanup and is **not present on any branch tip**,
local or remote (checked all 13 branches). `git show 60c2992:<path>` is the only way
back to it. Looking on another computer will not help.

**Seven columns, not four:** No adaptation, SFA, STD, **STF**, SFA+STD, STD+STF,
SFA+STD+STF.

**Why it cannot simply be restored — two independent blockers.**

1. *It never used a model class.* It hand-built an `n=1, K=1` params struct with
   `W = 0` and called `SRNNModelCellTypes.dynamics_fast_ct` / `unpack_states_ct`
   directly, with `SRNNModelBase.piecewiseSigmoid`. **`SRNNModelCellTypes` and
   `SRNNModelBase` are both deleted.** Nothing it calls exists.
2. *Its STF equations are superseded.* The facilitation form it documented is not the
   one `SRNNCellTypePairs` implements:

   | | deleted script | current `SRNNCellTypePairs` |
   |---|---|---|
   | facilitation state | `p`, rest `p0 = 0.35` | `g`, ceiling `G` |
   | STD depletion | `db/dt = (1−b)/τ_rec − (p·b·r)/τ_rel` — **coupled to `p`** | `db/dt = (1−b)/τ_rec − (b·r)/τ_rel` — independent of `g` |
   | synaptic gain | `eff = (p/p0)·b` | `g·b`, with `dg/dt = (1−g)/τ_dec + (G−g)·r/τ_fac` |

   So the old figure asserts that facilitation modulates depletion. The current model
   says it does not. Restoring the image would publish a claim the code contradicts.

**What rebuilding actually costs.** A new figure on `SRNNCellTypePairs`' facilitation,
not a restoration. It needs `tau_dec`, `tau_fac`, `G` values that **no preset carries**,
plus a decision on which route facilitates. Note the paper's target preset
`celltype_pairs_Sc0p2_noise0p025_dualStd` has **no STF on any route**, so an STF panel
would document a mechanism the paper does not run — worth deciding whether the paper
wants facilitation in scope at all before drawing it.

**Old values, for reference** (from `60c2992`, exaggerated for legibility):
`step_amp 0.5`, `t_on 5`, `t_off 15`, `T_range [-10 20]`, `fs 400`, `tau_d 0.1`,
piecewise `S_a 1.0` / `S_c 0.5`, `tau_a 3`, `c 1.0`, `tau_b_rec 2`, `tau_b_rel 0.3`,
`p0 0.35`, `tau_f 6`, `kappa 0.4`.

---

## Update CLAUDE.md's equations, and write an SFA version of the parameter table

| | |
|---|---|
| Noted | 2026-08-20 11:05 · `dev` @ `fb3639e` · R5469844 · Claude Code (Opus 5), session 5f04ea68 |

Three related pieces of documentation drift, raised by TR. **Not started.**

1. **The equation block in `CLAUDE.md` ("What this project is") is out of date.**
   It is written for `SRNNModel2`: a scalar `b_i` depressing every outgoing synapse,
   no per-route synapses, no STF, and no noise term (`dx_i` is shown as a plain ODE
   even though the model now supports additive Wiener noise on `x`).

2. **`SRNNCellTypePairs` is the primary class and CLAUDE.md should say so.** It is
   currently presented alongside `SRNNModel2` as a duck-typed sibling, with
   "default to `SRNNModel2` + `ParamSpaceAnalysis2`" as the closing advice. That no
   longer reflects where the work is.

3. **Write an SFA version of `docs/EquationsParametersDocs/cell_type_pair_parameter_table.md`.**
   That document covers STD and STF per route, and the additive Wiener noise on `x`,
   but has **no mention of SFA at all** — no `a_{i,k}` states, no `c`, no `tau_a`.
   Wanted: the same treatment for adaptation, in the per-cell-type notation.

---

## Separate the network seed from the noise seed for reps

| | |
|---|---|
| Noted | 2026-08-20 10:32 · `dev` @ `13d66d6` · R5469844 · Claude Code (Opus 5), session 5f04ea68 |
| Origin | FR-005, raised 2026-08-13 · `dev` @ `ac59b42` · R5456622 |

`noise_seed` derives from `rng_seeds(1)`, so a rep varies the network **and**
the noise realisation together. That is the right default — it is what makes
the noise path shared across the four adaptation conditions at a grid point —
but it means the within-level spread in a reps histogram mixes two sources and
neither can be attributed.

Wanted: a way to hold one fixed while varying the other, so "how much of this
spread is the network and how much is the noise" becomes answerable.
