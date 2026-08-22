# Single-cell-type support in `SRNNCellTypePairs`

What it would take to make `n_cellTypes = 1` work, and what could then be deleted.

Investigated 2026-08-22 on `refactorRunAll`. **Every claim below about what does
and does not work was verified by running it**, not inferred from reading.

---

## TL;DR

**The fix is one line.** Everything else in the class already handles `C = 1`
correctly — `RMTBlocks`, `validate()`, `complete_type_defaults`, `n_per_type`,
the `R` dependent, and even `type_colors` (which documents its `n == 1`
behaviour explicitly). One statement in `build_W` bypasses the only supported
way to change the number of types, and that is the whole defect.

```matlab
% src/SRNNCellTypePairs.m, in build_W, replacing lines 1092-1094:
rmt.set_types(obj.f, obj.mu_tilde, obj.sigma_tilde);
```

Beyond that, a **second, optional** change to `srnn_adaptation_conditions` makes
the four adaptation regimes usable at `C = 1`. Together they let two figure
presets drop their two-type workarounds and remove roughly 60 lines of
explanation-of-a-workaround from the repo.

The sweep pipeline is a separate question and probably **should not** be made
`C = 1`-capable — see §5.

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

`src/SRNNCellTypePairs.m:1090-1095`, inside `build_W`:

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

`build_W` does exactly what that forbids. It only works at `D = 2` because that
is what `RMTBlocks` already defaults to, so `D` never actually changes.

Confirmed in isolation:

```
r.f = 1              ->  numel(r.f) == 2,  f = [1 0]     (scalar-f rule expands)
r.set_types(1, 0, 1) ->  numel(r.f) == 1,  builds correctly
```

### `RMTBlocks` itself is fine at D = 1

Verified end to end:

| check | zero-mean (`mu=0, sigma=1, alpha=1`) | nonzero-mean (`mu=5.5, sigma=1.5, alpha=0.2`) |
|---|---|---|
| `W` | `200x200`, 50.2% negative | builds |
| `R` | 14.1421 | 32.5269 |
| `lambda_O` | `0` (no outlier, as expected for zero mean) | `220` |
| realized spectral radius | 14.38 | 221.8 |
| `plot_W`, `plot_spectrum` | OK | — |
| `is_column_uniform` | `1` | — |

No `RMTBlocks` change is needed.

---

## 2. The minimal fix

### 2a. `build_W` — required, one line

```matlab
rmt = RMTBlocks(obj.n);
rmt.alpha    = obj.alpha;
rmt.set_types(obj.f, obj.mu_tilde, obj.sigma_tilde);   % atomic; the only way to change D
rmt.zrs_mode = obj.zrs_mode;
```

`alpha` and `zrs_mode` stay as plain assignments — they are not part of `D`.
`set_types` validates that `mu_tilde` and `sigma_tilde` are both `D x D`, so this
is also a *tightening*: a mismatched pair currently reaches `RMTBlocks` and fails
with a message about `f`, where it would now fail naming both matrices.

**Risk: essentially none.** At `D = 2` the three assignments and `set_types` do
the same thing in the same order (`set_types` assigns the matrices first, then
`f`). Every existing preset is `D = 2`, so bit-identical output is expected —
but it must be *checked*, not assumed. See §6.

### 2b. `srnn_adaptation_conditions` — needed to *use* C = 1

`src/srnn_adaptation_conditions.m:83-88` hardcodes two-element rows:

```matlab
sfa_row = [n_a_sfa 0];
struct('name','no_adaptation', 'n_a',[0 0], ...)
```

At `C = 1` the model wants a **1-element** `n_a`. The helper does not currently
know `C`. Two options:

- **(a) Take `C` as an argument.** It already takes `model_class`,
  `synapse_config` and `n_a_sfa`; a fourth `n_cell_types` (default 2) is
  consistent and explicit.
- **(b) Derive `C` from the `synapse_config` route names.** Fragile — the
  `no_adaptation` and `sfa_only` conditions pass an *empty* struct, so there are
  no route names to count.

**(a).** Then:

```matlab
sfa_row = [n_a_sfa, zeros(1, C-1)];
zero_row = zeros(1, C);
```

which is `[3 0]` / `[0 0]` unchanged at `C = 2`.

`srnn_param_preset` would pass `numel(d.f)` — it already resolves the preset
struct before calling, so `C` is in hand.

### 2c. The scalar aliases — leave alone

`f_E`, `mu_EE_relative`, `mu_EI_relative`, `mu_IE_relative`, `mu_II_relative`
call `assert_two_types()` and correctly error at `C = 1` (verified). They *should*:
"the fraction excitatory" and "the E←I block" are meaningless with one type.

`tau_a_E` is different — it calls `assert_first_type_is_E()`, which permits any
number of types, and **already works at C = 1** (verified: returns `[]` when
`n_a = 0`). No change.

---

## 3. What could then be deleted

### 3a. `sompolinsky_pairs` — the clearest win

Figure 1 panel A is a Dale-free random network: one population, zero-mean
Gaussian weights, `tanh`. It is *conceptually* one cell type. It is currently
**two statistically identical types named `A`/`B`** purely to dodge this defect.

| | now | after |
|---|---|---|
| `n_cellTypes` | 2 | 1 |
| `cell_type_names` | `{'A','B'}` | `{'E'}` — or any honest single name |
| `f` | `[0.5 0.5]` | `1` |
| `mu_tilde_relative` | `zeros(2)` | `0` |
| `sigma_tilde_relative` | `ones(2)` | `1` |

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

### 3b. `single_neuron_stf` and the single-neuron methods figure

Currently `n = 2`, `indegree = 1`, `mu` and `sigma` both zero, giving
`W = [0 0; 0 0]` — two neurons where one is wanted, because `n >= n_cellTypes`
forces `n >= 2` at `C = 2`.

After the fix: `C = 1`, `n = 1`, `indegree = 1`, zero weights → `W = 0` (1x1).
A genuinely single neuron.

Deletions:

- **`readme_n2()`** in `fig_adaptation_methods.m` — ~10 lines explaining why
  `n = 2` and why only the E neuron is plotted.
- The "two neurons where one is wanted" paragraph in the preset comment.
- The `first_row(...)` indexing in `row_trace` becomes trivially correct rather
  than "take neuron 1 of 2, and trust that `W == 0` means the other cannot
  influence it".

Net: roughly **20 lines**, plus the figure stops needing a caveat about what it
is actually showing.

**Caveat worth checking:** `indegree = 1` with `n = 1` means a self-connection is
possible. With `mu = sigma = 0` the weight is zero either way, so it does not
matter here — but if anyone later wants a *connected* one-neuron model, the
self-coupling semantics need thinking about. Not a blocker for these figures.

### 3c. `fig_stim_engages_adaptation` — NOT a workaround

Its `concat_types()` helper (`:175`) looks similar but is legitimate: that
network genuinely has E and I populations, and "all neurons" genuinely means
concatenating them. **Leave it.** Worth stating so a later cleanup does not
delete it by pattern-match.

### 3d. `UserNotes.md`

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
| **1** | `build_W` uses `set_types` | **Do it.** One line, no behaviour change at `D = 2`, removes a latent trap and a confusing error. |
| **2** | `srnn_adaptation_conditions` takes `C` | **Do it with tier 1.** Small, and without it a `C = 1` model cannot use the four named regimes — which is what a figure would want. |
| **3** | Make the sweep pipeline `C = 1`-capable | **Do not.** |

Tier 3 is not worth it because a one-type network has nothing the sweeps are
*for*. The seven 1-D axes are `n`, `f_E`, `level_of_chaos` and four connectivity
blocks; at `C = 1` two of those are meaningless and four collapse into one. What
remains — `n` and `level_of_chaos` — is a Sompolinsky gain sweep, which is a
figure, not a parameter-space study. Adding conditional axis selection to the
sweep scripts would buy an analysis nobody wants and add a branch to code that
is currently uniform.

The honest boundary is: **`C = 1` is for figures, `C >= 2` is for sweeps**, and
the `assert_two_types` errors are what enforce it. Fixing tiers 1 and 2 does not
blur that line; fixing tier 3 would.

---

## 6. Verification plan

The fix is one line, so the test burden is mostly *proving nothing changed at
D = 2*.

1. **Bit-identical `W` at `D = 2`.** For each existing Pairs preset, build with a
   fixed `rng_seeds` before and after the change and assert `isequal(W_before,
   W_after)`. This is the load-bearing check: every current preset and every
   committed figure depends on it. A dedicated `scripts/tests/test_pairs_D1.m`
   should carry it rather than leaving it to the figure suite.
2. **`C = 1` builds and runs.** Construct the Sompolinsky network at `C = 1`,
   `build()`, `run()`, and assert `R == level_of_chaos` and ~50% negative
   weights — the same two numbers that currently validate the two-type
   workaround (measured: `R` = 0.900 / 1.600 / 2.500 exactly, 50.0% negative).
3. **`C = 1` figure parity.** `fig_introductory_concepts` at `C = 1` should
   reproduce the `C = 2` LLEs to within seed noise: currently −0.0997 / −0.0011 /
   +0.0989 at gammas 0.9 / 1.6 / 2.5. Note the *networks are not identical* —
   the RNG stream differs when `D` changes — so this is a distributional check,
   not `isequal`.
4. **Paths not yet exercised at `C = 1`.** These were unreachable while `build()`
   failed, so they are **unverified**, not known-good:
   - `plot()` and `plot_celltypes()` — `type_colors(1)` is handled explicitly and
     `draw_order(1) = 1`, so both look safe, but neither has run.
   - `compile_synapse_config` with a single route `E.E`.
   - `plot_eigenvalues`, `plot_W`, `plot_W_spectrum`, `plot_weight_histogram`.
   - `ParamSpaceAnalysis2.srnn_property_info` — unaffected (it reads the class,
     not an instance), but `effective_param` on a `C = 1` result is untested.
5. **Full regression.** The nine tests that currently pass, plus
   `make_all_paper_figures` twice, since two figures change preset shape.

---

## 7. Estimated effort

| Item | Size |
|---|---|
| `build_W` one-line fix | minutes |
| `srnn_adaptation_conditions` `C` argument | ~10 lines, 2 call sites |
| New `test_pairs_D1.m` (D=2 invariance + D=1 build/run) | ~80 lines |
| `sompolinsky_pairs` preset + `fig_introductory_concepts` cleanup | ~40 lines deleted, 1 figure re-run |
| `single_neuron_stf` preset + `fig_adaptation_methods` cleanup | ~20 lines deleted, 2 figure variants re-run |
| Verify §6.4 paths, fix whatever surfaces | unknown — the real risk |
| `UserNotes.md` strike-through, `CLAUDE.md` note | minutes |

Half a day if §6.4 is clean; the plotting paths are the only genuine unknown.

**It is not urgent.** The workarounds are correct, measured, and documented, and
they produce the right networks. This is a tidiness-and-honesty change — it
removes ~60 lines of explaining-around-a-defect and lets two figures say plainly
what they are — not a correctness fix.

---

## Appendix: verification commands used

```matlab
% RMTBlocks at D = 1 (works today)
r = RMTBlocks(200); r.alpha = 1; r.set_types(1, 0, 1);
r.W; r.R; r.lambda_O; r.plot_W(); r.plot_spectrum(); r.is_column_uniform();

% The piecemeal path that build_W takes (fails)
r = RMTBlocks(200); r.alpha = 1;
r.f = 1;                       % -> numel(r.f) == 2, f = [1 0]
r.mu_tilde = 0; r.sigma_tilde = 1;
r.W                            % RMTBlocks:InconsistentTypes

% SRNNCellTypePairs at C = 1: constructs, resolves, fails only in build()
m = SRNNCellTypePairs('n_cellTypes',1, 'cell_type_names',{'E'}, 'f',1, ...
    'mu_tilde_relative',0, 'sigma_tilde_relative',1, 'n',200, 'indegree',200, ...
    'activation','tanh', 'lya_method','none', 'T_range',[0 1]);
m.n_per_type    % 200
m.R             % 1.0000
m.tau_a_E       % []   (works -- does not require two types)
m.f_E           % SRNNCellTypePairs:NotTwoTypes (correct)
m.build()       % RMTBlocks:InconsistentTypes
```
