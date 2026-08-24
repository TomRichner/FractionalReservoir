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

## ~~`SRNNCellTypePairs` cannot build a ONE-cell-type network~~ — FIXED 2026-08-23

| | |
|---|---|
| Noted | 2026-08-22 · `refactorRunAll` @ `bbc013b` · R5611351 · Claude Code (Opus 5), session e22d2fab |
| Found | during `refactorRunAll`, porting the Sompolinsky and single-neuron figures |
| **Fixed** | 2026-08-23 on `refactorRunAll`. `build_network` now calls `rmt.set_types(...)`. Verified bit-identical at `D = 2` across all 14 two-type presets. |

> **This entry is now history, not a to-do.** It is kept because the diagnosis
> below is what a later session would otherwise re-derive, and because the
> *shape* of the bug — a class configured piecemeal where an atomic setter was
> required, failing with an error that blamed the caller — is worth recognising
> again.
>
> Three things changed since it was written:
>
> 1. **The method is `build_network`, not `build_W`.** There is no `build_W`;
>    the name below is wrong throughout.
> 2. **The plotting paths at `C = 1` were the flagged unknown, and they are all
>    clean** — `plot`, `plot_celltypes`, `plot_W`, `plot_W_spectrum`,
>    `plot_weight_histogram`, `plot_eigenvalues`, Benettin and QR.
> 3. **A larger defect surfaced while fixing this.** `fig_adaptation_methods`
>    variant A was building the paper's whole `n = 500` recurrent network for a
>    figure captioned "one unconnected neuron", because `paper_config` handed it
>    `cfg.preset_name`. That is what the new `single_neuron_dualStd` preset
>    fixes. See `singleCellTypeRefactor.md` §3c.
>
> The workarounds described at the end of this entry are **gone**:
> `sompolinsky_pairs` is now one type named `'all'`, and `single_neuron_stf` is
> genuinely `n = 1`. Regression: `scripts/tests/test_pairs_single_celltype.m`.

A latent defect, verified empirically rather than inferred. It is a one-line fix,
but it is in `SRNNCellTypePairs`, so it was left alone rather than folded into a
refactor about something else.

### What happens

```matlab
SRNNCellTypePairs('n_cellTypes',1, 'cell_type_names',{'E'}, 'f',1, ...
                  'mu_tilde_relative',0, 'sigma_tilde_relative',1, ...)
```

constructs **correctly** — the object holds `n_cellTypes = 1`, `numel(f) = 1`,
`mu_tilde_relative` of size `[1 1]`. Then `build()` throws:

```
RMTBlocks:InconsistentTypes
numel(f) = 2 but mu_tilde is [1 1]. Use set_types(f, mu_tilde, sigma_tilde)
to change the number of cell types.
```

The error names the fix but blames the wrong caller, which is what makes it
confusing: nothing in the user's call said 2.

### Why

`SRNNCellTypePairs.build_W`, `src/SRNNCellTypePairs.m:1090-1094`, assigns the
generator **piecemeal**:

```matlab
rmt = RMTBlocks(obj.n);           % constructs at D = 2 (RMTBlocks.m:155, f = [0.5 0.5])
rmt.alpha       = obj.alpha;
rmt.f           = obj.f;          % scalar 1 -> the D=2 setter expands it to [1 0]
rmt.mu_tilde    = obj.mu_tilde;   % 1x1 -> now inconsistent with numel(f) = 2
rmt.sigma_tilde = obj.sigma_tilde;
```

`RMTBlocks.set_types`'s own docstring (`RMTBlocks.m:563`) says exactly why this
fails:

> This is the only way to **CHANGE D**. Setting `f`, `mu_tilde` and
> `sigma_tilde` one at a time works fine when D is unchanged, but changing D
> piecemeal would leave the object transiently inconsistent, so use this instead.

`build_W` is doing the thing that docstring forbids. It works for D = 2 only
because that is what `RMTBlocks` already defaults to, so D never actually changes.

**Confirmed directly**, without editing anything: `r.f = 1` leaves
`numel(r.f) == 2` with `f = [1 0]`, whereas `r.set_types(1, 0, 1)` builds a D = 1
network correctly (200x200 dense, 50.3% negative weights).

### The fix

Replace the three assignments with the atomic setter:

```matlab
rmt.set_types(obj.f, obj.mu_tilde, obj.sigma_tilde);
```

`alpha` and `zrs_mode` stay as they are — they are not part of D.

**But that alone does not make D = 1 usable in the sweep pipeline.** Two
scalar aliases assert two types: `f_E` calls `assert_two_types()`, and `tau_a_E`
additionally requires `cell_type_names{1} == 'E'`. Since `f_E` is the
fraction-excitatory sweep axis and `tau_a_E` is the tau-sensitivity axis, a
one-type preset would still fail there. Fixing `build_W` buys *construction*,
not a swept one-type model.

### Two other constraints that bite alongside it

Both are deliberate and correct, but they compound with the above when trying to
build something small:

- `n >= n_cellTypes` (`SRNNCellTypePairs.m:854`) — so `n = 1` needs D = 1, which
  is the broken path.
- `0 < indegree <= n` (`:857`) — `indegree = 0` is rejected, so "no recurrence"
  cannot be expressed by disconnecting; it must be expressed with zero weights.

### How the refactor worked around it

Neither workaround needs model-class changes, which is why the defect stayed
unfixed:

- **`sompolinsky_pairs`** — the Dale-free random network for figure 1 panel A
  uses **two statistically identical zero-mean types** named `A`/`B`. Same
  network; `R` comes out exactly equal to `level_of_chaos` (0.900 / 1.600 /
  2.500 measured), and the weights are 50.0% negative as required. The traces
  concatenate both types, since together they are the whole network.
- **`single_neuron_stf`** and the single-neuron methods figure — `n = 2`,
  `indegree = 1`, `mu` and `sigma` both zero, giving `W = [0 0; 0 0]` and
  `N_sys_eqs = 5`. The smallest expressible unconnected network; only the E
  neuron is plotted, and with `W == 0` the other cannot influence it.

The cost of the workaround is cosmetic but real: type names are either arbitrary
(`A`/`B`) or slightly dishonest (`E`/`I` for a lone neuron), and `plot_data` is
keyed by type so any "all neurons" trace has to concatenate across types.

---

## Memory capacity is still on `SRNNModel2` — the ESN class needs porting to `SRNNCellTypePairs`

| | |
|---|---|
| Noted | 2026-08-22 · `refactorRunAll` @ `40b578a` · R5611351 · Claude Code (Opus 5), session e22d2fab |
| Raised by | TR, at the end of the `refactorRunAll` refactor |

**The one part of the paper still on a different model class.** Everything else was
moved onto `SRNNCellTypePairs` during the refactor; memory capacity could not be,
and this is why. **Not started.**

### What blocks it

`src/SRNN_ESN_reservoir.m` is declared

```matlab
classdef SRNN_ESN_reservoir < SRNNModel2
```

so it inherits `SRNNModel2`'s vocabulary wholesale — `n_a_E`, `n_b_E`,
`mu_E_tilde_relative`, `tau_b_E_rec`, and that class's state packing. None of
those exist on `SRNNCellTypePairs`. The two classes are **duck-typed siblings,
not a hierarchy** (CLAUDE.md), so there is no base to re-parent to: the port is a
rewrite of the ESN layer against the other class, not a one-line `classdef` edit.

### What that costs the paper right now

- `fig_memory_capacity` and `fig_memory_capacity_example` show a **different
  network from every other figure** — n = 300, f = 0.6, `level_of_chaos` 2.0,
  *logistic* nonlinearity at `S_c = 0.35`, no noise — against the paper's
  n = 500, piecewise, `S_c = 0.2`, `sigma_u_noise = 0.025`, dual-STD network.
- The methods section **has to say so**. Both figures' generated READMEs now
  carry that caveat, and `paper_config.m` names `mc_esn` as a deliberate
  per-figure override rather than hiding it.
- `mc_esn` is the **only `SRNNModel2` preset left** in `srnn_param_preset.m` that
  the paper uses. Killing it is the whole point of this note.

### The four real problems a port has to solve

1. **What does the readout read?** MC uses `readout_signal = 'synaptic'`, i.e.
   `br = b.*r`, so the linear readout can see the STD state. On
   `SRNNCellTypePairs` the synaptic output is **per route** —
   `plot_data.synaptic_output.E.E`, `.E.I`, … — so there is no single `br`
   matrix. A presynaptic neuron has a *different* effective output per
   postsynaptic target. Decide deliberately: read `r` only; read one route; or
   read the per-neuron product across the routes it actually projects on. This
   is a modelling choice, not a mechanical translation, and it changes what MC
   measures.
2. **Conditions.** MC's four regimes are spelled `n_a_E` / `n_b_E`. On Pairs they
   are an `n_a` row plus a whole `synapse_config`. `srnn_adaptation_conditions`
   already returns both spellings under the same four names, so the fix is to
   take conditions from there rather than from the hardcoded `condition_args` in
   `run_memory_capacity`.
3. **`verify_shared_build`.** It currently checks 88 `SRNNModel2` properties and
   asserts exactly `n_a_E` / `n_b_E` / `tau_a_E` differ. It needs a Pairs
   equivalent, or to be made class-agnostic by reading the property list off
   `meta.class.fromName` the way `ParamSpaceAnalysis2.srnn_property_info` does.
   Do not drop the check — it is what makes the paired-trial comparison honest.
4. **`SRNN_ESN_reservoir` does not go through `run()`.** It has its own
   `run_reservoir_esn`, which is precisely why CLAUDE.md warns that the
   protected `integrate()` helper exists "to keep it from silently missing
   integrator work". A Pairs port must route through `integrate()` or it will
   quietly miss the SDE solvers — and the paper's preset is stochastic, so that
   matters immediately.

### Where the code is now

- `scripts/memory_capacity/run_memory_capacity.m` — the paired-trial ensemble, a
  function with its own `mc_run_config` cost table. Already takes a preset.
- `scripts/presentations/Stability_Manuscript/fig_memory_capacity_example/run_memory_capacity_example.m`
  — the single-network compute half of the example figure.
- `src/srnn_param_preset.m`, case `'mc_esn'` — the network, carrying physics only
  (the protocol lives in `mc_run_config`, deliberately).

Because both entry points already take `preset_name` and both resolve their
model class from it, the port is contained: once an ESN exists that speaks
`SRNNCellTypePairs`, pointing `paper_config.cfg.mc_preset` at the main preset is
the only change those two figures need.

---

## Should STF scale STD depletion? `SRNNCellTypePairs` says no; Tsodyks–Markram says yes

| | |
|---|---|
| Noted | 2026-08-21 · `refactorRunAll` @ `22a91ee` · R5611351 · Claude Code (Opus 5), session e22d2fab |
| Raised by | TR, while settling `tau_rel` for the rebuilt STF methods figure |

**To look into after the refactor.** Parked deliberately; the figure ships with
`tau_rel = 0.3` in the meantime.

The two depression equations differ in whether facilitation feeds back into depletion:

```
deleted SRNNModelCellTypes:  db/dt = (1-b)/tau_rec - (p * b * r)/tau_rel
current SRNNCellTypePairs:   db/dt = (1-b)/tau_rec - (    b * r)/tau_rel
```

**The old form is the Tsodyks–Markram one.** In TM, the resource variable depletes in
proportion to the *utilization* — release probability times available resource — so a
facilitated synapse consumes its pool faster. The current class drops that factor, so
`b` and `g` evolve independently and only multiply at the output.

Why it may not matter: with **no STF anywhere** (which is every preset in use, the
paper's `celltype_pairs_Sc0p2_noise0p025_dualStd` included), the missing factor is a
constant and is fully absorbed into `tau_rel`. The two forms are then the same model
with a rescaled constant. **The divergence only appears once facilitation is on.**

Why it may: any future STF work — a facilitating route, an E→I vs I→I comparison, the
STF figure itself — silently gets a synapse whose depression does not respond to its
own facilitation. That is a modelling claim, not an implementation detail, and it is
currently unstated anywhere.

Concrete consequence already in hand: the rebuilt STF methods figure uses
`tau_rel = 0.3` verbatim from the old script, which — with the `p = p0 = 0.35` factor
gone — is about **2.9× stronger depression at rest** than the archived
`sfa_std_stf_single_neuron_example_figure_1.png` shows. Rest could have been matched
with `tau_rel = 0.3/0.35 = 0.857`, but nothing reproduces the *acceleration* as `p`
rises toward 1. TR chose the literal 0.3 as the simpler thing to explain.

Worth deciding: is the decoupling intentional (a deliberate simplification worth
documenting in CLAUDE.md and the equations doc), or an oversight from the port that
should be restored?

---

## ~~Rebuild the single-neuron STF methods figure~~ — SUPERSEDED, now in the refactor

| | |
|---|---|
| Noted | 2026-08-21 · `refactorRunAll` @ `22a91ee` · R5611351 · Claude Code (Opus 5), session e22d2fab |
| **Superseded** | 2026-08-21, same session. TR reversed the parking decision: the STF figure **is** being rebuilt as part of `refactorRunAll`, on `SRNNCellTypePairs` behind a new preset matching the old parameters. See `refactorRunAll.md` §6. |

> **This entry is now history, not a to-do.** It is kept because the archaeology below
> (which commit holds the last version, why it cannot be restored verbatim) is still
> what a later session would otherwise have to re-derive. The *decision* it recorded —
> "leave this outside the refactor" — no longer holds.
>
> One finding since: the facilitation parameters **do** map exactly. Writing the old
> `dp/dt = (p0−p)/tau_f + kappa(1−p)r` with gain `u = p/p0` gives
> `du/dt = (1−u)/tau_f + kappa(1/p0 − u)r`, which is the current class's
> `dg/dt = (1−g)/tau_dec + (G−g)r/tau_fac` with `tau_dec = tau_f = 6`,
> `tau_fac = 1/kappa = 2.5`, `G = 1/p0 = 2.857`. Only the **STD coupling** fails to
> map: the old depletion term carried a factor `p`, the current one does not.

Original entry follows.

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
