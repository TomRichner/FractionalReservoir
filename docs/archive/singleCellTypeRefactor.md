# Single-cell-type support in `SRNNCellTypePairs`

> ## ⚠ ARCHIVED — this refactor is COMPLETE. Do not treat it as a plan.
>
> **Status: implemented 2026-08-23** (`refactorRunAll` @ `0134e9b`). The
> `build_network` fix landed, `srnn_adaptation_conditions` took its cell-type
> count argument, and the three `C = 1` presets are live: `sompolinsky_pairs`,
> `single_neuron_stf`, `single_neuron_dualStd`. Every §6 verification came back
> ANSWERED, `W` was bit-identical at `C = 2`, and the regression is frozen in
> `scripts/tests/test_pairs_single_celltype.m`.
>
> **One conclusion here is a standing design rule, and it is easy to misread as
> a bug** — §4, "What this does NOT unlock". `f_E` and the four
> `mu_*_relative` block aliases deliberately **error** at `C = 1`
> (`SRNNCellTypePairs:NotTwoTypes`), because the fraction excitatory and the
> four blocks are meaningless with one cell type. The boundary is **`C = 1` for
> figures, `C >= 2` for sweeps**, and those errors are what enforce it. Do not
> "fix" them, and do not teach the sweep scripts to select axes conditionally.
> That rule now lives in `CLAUDE.md`, which is the maintained statement of it.
>
> **This file is kept as a record of the reasoning, not as a description of the
> code.** It is archived, not maintained; its `src/` paths predate the
> 2026-09-02 reorganization. Nothing here should be acted on.
>
> Archived 2026-09-02.

What it takes to make `n_cellTypes = 1` work, and what it lets us delete.

Investigated 2026-08-22 and re-probed 2026-08-23 on `refactorRunAll`. **Every
claim below about what does and does not work was verified by running it**, not
inferred from reading.

---

## TL;DR

**The fix is one line.** Everything else in the class already handles `C = 1`
correctly — `RMTBlocks`, `validate()`, `complete_type_defaults`, `n_per_type`,
the `R` dependent, every plotting method, and `type_colors` (which documents its
`n == 1` behaviour explicitly). One statement in `build_network` bypasses the
only supported way to change the number of types, and that is the whole defect.

```matlab
% src/SRNNCellTypePairs.m, in build_network, replacing lines 1092-1094:
rmt.set_types(obj.f, obj.mu_tilde, obj.sigma_tilde);
```

Beyond that, a **second** change to `srnn_adaptation_conditions` makes the four
adaptation regimes usable at `C = 1`, which is what a figure actually wants.

The sweep pipeline is a separate question and **should not** be made
`C = 1`-capable — see §5.

**A third, larger defect surfaced while probing this.** `fig_adaptation_methods`
variant A is not a single neuron at all — see §3c. Fixing it is the main *payoff*
of `C = 1`, not a side errand.

---

## 1. What is actually broken

### The symptom

```matlab
SRNNCellTypePairs('n_cellTypes',1, 'cell_type_names',{'E'}, 'f',1, ...
                  'mu_tilde_relative',0, 'sigma_tilde_relative',1, ...)
```

**constructs fine.** Verified — the object holds:

| | |
|---|---|
| `n_cellTypes` | `1` |
| `numel(f)` | `1` |
| `size(mu_tilde_relative)` | `[1 1]` |
| `size(mu_tilde)` | `[1 1]` |
| `n_per_type` | `200` (all neurons in the one type) |
| `n_a` | `0` (correct 1-element row) |
| `c` | `0.05` (1-element) |
| `tau_a` | `1x1 cell` |
| `R` (dependent) | `1.0000` — **the class's own RMT maths works at C = 1** |

So `validate()`, `complete_type_defaults`, `expand_block` and the dependent
properties are all already correct. Then `build()` throws:

```
RMTBlocks:InconsistentTypes
numel(f) = 2 but mu_tilde is [1 1]. Use set_types(f, mu_tilde, sigma_tilde)
to change the number of cell types.
```

The message names the fix but blames the wrong caller — nothing in the user's
call said 2 — which is why this reads like a bug in the calling code.

### The cause

`src/SRNNCellTypePairs.m:1090-1095`, inside **`build_network`** (there is no
method called `build_W`; an earlier draft of this document called it that):

```matlab
rmt = RMTBlocks(obj.n);           % constructs at D = 2  (RMTBlocks.m:155)
rmt.alpha       = obj.alpha;
rmt.f           = obj.f;          % scalar 1 -> the D=2 setter expands it to [1 0]
rmt.mu_tilde    = obj.mu_tilde;   % 1x1 -> now inconsistent with numel(f) = 2
rmt.sigma_tilde = obj.sigma_tilde;
rmt.zrs_mode    = obj.zrs_mode;
```

`RMTBlocks.set_types`'s own docstring (`RMTBlocks.m:563`) states the rule this
breaks:

> This is the only way to **CHANGE D**. Setting `f`, `mu_tilde` and
> `sigma_tilde` one at a time works fine when D is unchanged, but changing D
> piecemeal would leave the object transiently inconsistent, so use this instead.

`build_network` does exactly what that forbids. It only works at `D = 2` because
that is what `RMTBlocks` already defaults to, so `D` never actually changes.

Confirmed in isolation:

```
r.f = 1              ->  numel(r.f) == 2,  f = [1 0]     (scalar-f rule expands)
r.set_types(1, 0, 1) ->  numel(r.f) == 1,  builds correctly
```

### `RMTBlocks` itself is fine at D = 1

Verified end to end:

| check | zero-mean (`mu=0, sigma=1, alpha=1`) | nonzero-mean (`mu=5.5, sigma=1.5, alpha=0.2`) | `n = 1`, zero weights |
|---|---|---|---|
| `W` | `200x200`, 49.7% negative | builds | `0` (`1x1`) |
| `R` | 1.0000 (at `sigma = 1/sqrt(n)`) | 32.5269 | 0 |
| `lambda_O` | `0` (no outlier, as expected) | `220` | 0 |
| `plot_W`, `plot_spectrum` | OK | — | — |
| `is_column_uniform` | `1` | — | — |

No `RMTBlocks` change is needed. Notably `n = 1` works, which is what makes a
genuinely single-neuron model expressible.

---

## 2. The fix

### 2a. `build_network` — required, one line

```matlab
rmt = RMTBlocks(obj.n);
rmt.alpha    = obj.alpha;
rmt.set_types(obj.f, obj.mu_tilde, obj.sigma_tilde);   % atomic; the only way to change D
rmt.zrs_mode = obj.zrs_mode;
```

`alpha` and `zrs_mode` stay plain assignments — they are not part of `D`, and
keeping `alpha` **before** `set_types` preserves the existing order (`set.alpha`
is the only one of these setters that can touch the RNG, via `update_sparsity`).

This is also a *tightening*: `set_types` validates that `mu_tilde` and
`sigma_tilde` are both `D x D`, so a mismatched pair now fails naming both
matrices instead of complaining about `f`.

**Risk: none, measured.** At `D = 2` the three assignments and `set_types` do the
same thing in the same order (`set_types` assigns the matrices first, then `f`),
and none of the setters consumes RNG. See §6.1 for the numbers.

### 2b. `srnn_adaptation_conditions` — needed to *use* C = 1

`src/srnn_adaptation_conditions.m:83-88` hardcodes two-element rows:

```matlab
sfa_row = [n_a_sfa 0];
struct('name','no_adaptation', 'n_a',[0 0], ...)
```

At `C = 1` the model wants a **1-element** `n_a`. Add a fourth argument
`n_cell_types` (default 2), consistent with the `synapse_config` and `n_a_sfa`
arguments already there, and then:

```matlab
sfa_row  = [n_a_sfa, zeros(1, n_cell_types-1)];
zero_row = zeros(1, n_cell_types);
```

which is `[3 0]` / `[0 0]` unchanged at `C = 2`. The `SRNNModel2` branch should
**error** if handed `n_cell_types ~= 2` rather than ignoring it — that class is
intrinsically two populations, and this repo's convention is that a name may not
lie (cf. `assert_first_type_is_E`).

**`srnn_param_preset` sources `C` from `d.n_cellTypes`, not `numel(d.f)`.** Every
`SRNNCellTypePairs` preset carries `n_cellTypes` explicitly; the `'default'`
preset has no `f` field at all, so `numel(d.f)` — which an earlier draft of this
document proposed — throws `Unrecognized field name "f"`.

```matlab
if isfield(d, 'n_cellTypes'); C = d.n_cellTypes; else; C = 2; end
conditions = srnn_adaptation_conditions(model_class, std_routes, n_a_sfa, C);
```

`scripts/tests/srnn_param_preset_old.m` is a frozen equivalence reference that
calls this with 1 and 2 arguments; the default keeps it working, and it only
enumerates two-type presets, so `test_srnn_param_preset_equivalence` is
unaffected.

### 2c. The scalar aliases — leave alone

`f_E`, `mu_EE_relative`, `mu_EI_relative`, `mu_IE_relative`, `mu_II_relative`
call `assert_two_types()` and correctly error at `C = 1` (verified). They *should*:
"the fraction excitatory" and "the E←I block" are meaningless with one type.

`tau_a_E` is different — it calls `assert_first_type_is_E()`, which permits any
number of types, and **already works at C = 1** (verified: returns `[]` when
`n_a = 0`). No change.

---

## 3. What this lets us delete or fix

### 3a. `sompolinsky_pairs` — the clearest win

Figure 1 panel A is a Dale-free random network: one population, zero-mean
Gaussian weights, `tanh`. It is *conceptually* one cell type. It is currently
**two statistically identical types named `A`/`B`** purely to dodge this defect.

| | now | after |
|---|---|---|
| `n_cellTypes` | 2 | 1 |
| `cell_type_names` | `{'A','B'}` | `{'all'}` |
| `f` | `[0.5 0.5]` | `1` |
| `mu_tilde_relative` | `zeros(2)` | `0` |
| `sigma_tilde_relative` | `ones(2)` | `1` |

`{'all'}` rather than `{'E'}`: the weights take both signs, so an E/I name would
be the same lie the `A`/`B` names were dodging.

Deletions that follow:

- **`all_type_states()`** in `fig_introductory_concepts.m:246` — a helper that
  exists only to concatenate `plot_data.x.A` and `plot_data.x.B` back into the
  single population they always were. ~12 lines including its comment.
- **The `readme_two_types()` section** in the same file — ~10 lines of generated
  README explaining why the types are named `A`/`B`.
- **The preset's 15-line comment block** justifying the two-type construction.
- The docstring paragraph in `fig_introductory_concepts` making the same point.

Net: roughly **40 lines**, all of which are explanation of a workaround rather
than of the science.

**Expected cost that did not materialise.** Changing `D` changes how the
generator partitions its blocks, so `W` was expected to be a different draw, and
TR accepted that on 2026-08-23. **It is bit-identical.** Measured at
`n = 120, indegree = 40, rng_seeds = [42 42]`: `sum(W) = -0.563793163494`,
`‖W‖_F = 13.5057856038`, `nnz = 4856`, before and after. The reason is that the
two types had *identical* zero-mean statistics, so the per-block scaling was
uniform and the underlying draw never depended on the partition at all.
Figure 1 panel A is unchanged, LLEs included (−0.0997 / −0.0011 / +0.0989 at
gammas 0.9 / 1.6 / 2.5, exactly as before). `test_pairs_single_celltype.m`
carries the checksum.

This does **not** generalise: a preset whose blocks differ (every E/I preset)
would genuinely redraw if its type count changed. It holds here only because the
`A`/`B` split was cosmetic in the first place — which is the same fact that made
it a workaround.

### 3b. `single_neuron_stf` — a genuinely single neuron

Currently `n = 2`, `indegree = 1`, `mu` and `sigma` both zero, giving
`W = [0 0; 0 0]` — two neurons where one is wanted, because `n >= n_cellTypes`
forces `n >= 2` at `C = 2`.

After the fix: `C = 1`, `n = 1`, `indegree = 1`, zero weights → `W = 0` (1x1).
Verified: SFA(3) + dual STD + STF on such a model gives `N_sys_eqs = 7` and
correctly shaped `plot_data` throughout.

Deletions:

- **`readme_n2()`** in `fig_adaptation_methods.m` — ~10 lines explaining why
  `n = 2` and why only the E neuron is plotted.
- The "two neurons where one is wanted" paragraph in the preset comment.
- The `first_row(...)` indexing in `row_trace` becomes trivially correct rather
  than "take neuron 1 of 2, and trust that `W == 0` means the other cannot
  influence it".

`{'E'}` is honest here — it is the neuron whose synapse depresses — and it keeps
the `tau_a_E` alias valid.

**Caveat worth checking:** `indegree = 1` with `n = 1` means a self-connection is
possible. With `mu = sigma = 0` the weight is zero either way, so it does not
matter here — but if anyone later wants a *connected* one-neuron model, the
self-coupling semantics need thinking about. Not a blocker for these figures.

### 3c. `fig_adaptation_methods` variant A is not a single neuron — the real find

Variant A ("SFA + STD, on the paper's preset") calls `build_from_preset` with the
`celltype_pairs_Sc0p2_noise0p025_dualStd` preset and **overrides no network
parameters**, so it builds that preset *whole*: `n = 500`, `indegree = 100`,
fully recurrent, `level_of_chaos = 1`. It then plots neuron 1.

The generated README states the contradiction outright — "One unconnected neuron"
in `WHAT IT SHOWS`, `n 500` / `indegree 100` under `PARAMETERS AS RUN` — because
`figure_settings` reads off the **built object** rather than the preset struct.
That is exactly the class of error that helper was written to catch, and it
caught this one; nobody read it.

The figure demonstrates none of the mechanisms:

| column | what it shows | what it should show |
|---|---|---|
| No adaptation | network chaos, `r` saturating between 0 and 1 | a clean step response |
| SFA only | chaotic bursts | a decaying rate with visible `a` |
| STD only | silent, `r ~ 0`, `b` relaxing to 1 | a step with depression |
| SFA + STD | silent | both mechanisms |

**Fix: a new `single_neuron_dualStd` preset**, derived from the dualStd preset
and changing only `C -> 1`, `n -> 1`, `indegree -> 1`, and the connectivity
blocks to zero. Everything the figure is *about* — `c = 0.5/3`, three `tau_a`,
piecewise `S_a = 0.8` / `S_c = 0.20`, `tau_d = 0.1`, and the dual STD route
`tau_rec = [2 4]` / `tau_rel = [0.25 0.5]` — is inherited unchanged.

**Noise stays on** (`sigma_u_noise = 0.025`), per TR. At `n = 1` unconnected,
`x_noise_std = sigma_u/sqrt(2*tau_d) = 0.056` against a 0.5 step — about 11%,
visible but legible, and honest about the model the paper characterises.

This preset is only expressible *because* of the `C = 1` fix, which is why the
two are done together.

### 3d. `fig_stim_engages_adaptation` — NOT a workaround

Its `concat_types()` helper (`:175`) looks similar but is legitimate: that
network genuinely has E and I populations, and "all neurons" genuinely means
concatenating them. **Leave it.** Worth stating so a later cleanup does not
delete it by pattern-match.

### 3e. `UserNotes.md`

The entry *"`SRNNCellTypePairs` cannot build a ONE-cell-type network"* becomes
history and should be struck through with a pointer to the fixing commit, the
way the STF entry already was.

---

## 4. What this does NOT unlock

Being able to *build* `C = 1` is not the same as being able to *sweep* it.

- **`f_E` is a sweep axis.** `run_sensitivity_analysis` and
  `run_param_space_analysis` both sweep `ctx.f_param`, which is `'f_E'` on
  Pairs — and `f_E` correctly errors at `C = 1`. The fraction excitatory is
  meaningless with one type.
- **The four `mu_*_relative` blocks are sweep axes.**
  `run_sensitivity_analysis:51` sweeps all four; at `C = 1` there is one block,
  and it is just "the weight mean".
- `resolve_run_context:120` sets `f_param = 'f_E'` for any Pairs model,
  unconditionally.

So a `C = 1` preset run through `run_all_analyses` would fail on the `f_E`
sweep. **That is correct behaviour and should stay** — see §5.

---

## 5. Recommendation: fix the class, not the pipeline

Three tiers, in decreasing order of value:

| Tier | Change | Verdict |
|---|---|---|
| **1** | `build_network` uses `set_types` | **Do it.** One line, measured no-op at `D = 2`, removes a latent trap and a confusing error. |
| **2** | `srnn_adaptation_conditions` takes `C`; the three presets and two figures follow | **Do it with tier 1.** Without it a `C = 1` model cannot use the four named regimes, and §3c stays broken. |
| **3** | Make the sweep pipeline `C = 1`-capable | **Do not.** |

Tier 3 is not worth it because a one-type network has nothing the sweeps are
*for*. The seven 1-D axes are `n`, `f_E`, `level_of_chaos` and four connectivity
blocks; at `C = 1` two of those are meaningless and four collapse into one. What
remains — `n` and `level_of_chaos` — is a Sompolinsky gain sweep, which is a
figure, not a parameter-space study. Adding conditional axis selection to the
sweep scripts would buy an analysis nobody wants and add a branch to code that
is currently uniform.

The honest boundary is: **`C = 1` is for figures, `C >= 2` is for sweeps**, and
the `assert_two_types` errors are what enforce it. Tiers 1 and 2 do not blur that
line; tier 3 would.

---

## 6. Verification — results

### 6.1 `D = 2` invariance — ANSWERED, bit-identical

Measured on the real preset, before and after the change, at `n = 120`,
`indegree = 40`, `rng_seeds = [42 42]`:

| | `celltype_pairs_Sc0p2_noise0p025_dualStd` |
|---|---|
| `sum(W)` | `-1.14156058489` |
| `‖W‖_F` | `29.7665779012` |
| `nnz(W)` | `4856` |
| `R` | `2.23192327` |
| first three nonzeros | `[0.5125129582  0.4115299747  0.5217920366]` |

Identical in both directions. A frozen table of these checksums for **every**
Pairs preset is carried by `scripts/tests/test_pairs_single_celltype.m`, so the
check is a test rather than a one-off. Two presets are exempt by design:
`sompolinsky_pairs` changes draw when it moves to `C = 1` (§3a), and
`single_neuron_stf` has `W == 0` either way.

### 6.2 `C = 1` builds and runs — ANSWERED

A Sompolinsky-shaped model at `C = 1`, `n = 200`, `indegree = 200`,
`level_of_chaos = 1.6`: `R = 1.600000` exactly, realized spectral radius
`1.6657`, 49.9% negative weights, Benettin `LLE = -0.0272`.

### 6.3 `C = 1` figure parity — ANSWERED, exact

`fig_introductory_concepts` at `C = 1` reproduces the `C = 2` result **exactly**:
LLE = −0.0997 / −0.0011 / +0.0989 at gammas 0.9 / 1.6 / 2.5, the same numbers as
before the change, because `W` is bit-identical (§3a). A distributional check was
planned; an equality check was available instead, and is what the test does.

Note gamma = 1.6 sits essentially at zero (−0.0011) and always did — that is a
property of the figure's chosen gains, not of this change.

### 6.4 Paths not previously exercised at `C = 1` — ANSWERED, all clean

These were unreachable while `build()` failed, and were this document's stated
"only genuine unknown". All were re-probed on 2026-08-23 with the fix applied,
and every one passes:

| path | result |
|---|---|
| `plot()` | OK |
| `plot_celltypes()` | OK (both at `n = 200` and at `n = 1` with SFA+STD+STF) |
| `plot_W`, `plot_W_spectrum`, `plot_weight_histogram` | OK |
| `plot_eigenvalues(times)` with `store_full_state` | OK |
| `compile_synapse_config` with a single `E.E` route | OK — `N_sys_eqs` correct |
| Benettin and QR Lyapunov | both OK |

`type_colors(1)` returns the type-1 warm colour and `draw_order(1) = 1`, which is
why the plotting held up.

Still untested: `ParamSpaceAnalysis2.effective_param` on a `C = 1` result. It is
unreachable from any supported path — §4 and §5 say `C = 1` is never swept — so
this is a boundary, not a gap.

### 6.5 Full regression

The nine-plus tests listed in the implementation plan, plus
`make_all_paper_figures`, since three presets change shape.

---

## 7. Effort

| Item | Size |
|---|---|
| `build_network` one-line fix | minutes |
| `srnn_adaptation_conditions` `C` argument | ~10 lines, 1 live call site |
| New `test_pairs_single_celltype.m` | ~120 lines |
| `sompolinsky_pairs` + `fig_introductory_concepts` cleanup | ~40 lines deleted |
| `single_neuron_stf` + `fig_adaptation_methods` cleanup | ~20 lines deleted |
| New `single_neuron_dualStd` preset + variant A repoint (§3c) | ~40 lines added |
| `UserNotes.md`, `CLAUDE.md` | minutes |

Half a day. The plotting risk that dominated the original estimate is gone —
§6.4 came back clean.

**Urgency changed.** The original note said this was tidiness only. §3c makes it
a **correctness** fix: one manuscript figure currently shows a chaotic
500-neuron network while claiming to show one unconnected neuron.
