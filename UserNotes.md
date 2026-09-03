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

## Presets should write out their conditions, not name a regime set

| | |
|---|---|
| Noted | 2026-09-03 · `refactorRunAll` @ `f570cf3` · R5611351 · Claude Code (Opus 5), session e22d2fab |
| Raised by | TR, after tracing where `sfa1_std1` and `sfa3_std2` are actually defined for the 3-condition preset |
| Status | **Recorded, not implemented.** Architecture; no science changes. Separable from the `log_ladder` work. |

### How it works now

Every `srnn_param_preset` case sets five locals — `d`, `model_class`,
`std_routes`, `sfa_timescales`, `regimes` — and then **one** call at the bottom
of the file (`srnn_param_preset.m:657`) expands them:

```matlab
conditions = srnn_adaptation_conditions(model_class, ...
    'synapse_config', std_routes, 'sfa_timescales', sfa_timescales, ...
    'n_cell_types', n_cell_types, 'regimes', regimes);
```

So a preset never names its conditions. It names a regime **set**
(`'standard'` / `'timescales'` / `'single_multi'`), and
`srnn_adaptation_conditions` turns that into the named structs. Answering "what
is `sfa1_std1`?" means reading the preset for the ladder and routes, then the
other file for how they are combined — four hops.

### What is genuinely good about it, and must survive any change

- **The one-timescale regimes are DERIVED.** `first_timescale_routes(sc)` keeps
  each route's first `tau_rec`/`tau_rel` pair, so retuning a preset's depression
  updates its `_oneTS` regimes automatically. Written out per preset, that
  derivation would be duplicated and would drift — which is the exact failure
  this repo keeps paying for.
- **A regime name means one thing everywhere.** `no_adaptation` is identical
  across all three sets, and `sfa3_std2` is asserted byte-identical to
  `sfa_and_std`. Two commits (`6201b2f`, `f570cf3`) went out of their way to
  keep that true.

### What is awkward

- **The physics is split across two files.** The preset owns the ingredients
  (routes, ladder); the other file owns the recipes.
- **Three locals exist only to feed one call.** `std_routes`, `sfa_timescales`
  and `regimes` are not network physics in the `d` sense — they are arguments to
  a factory that happens to live at the bottom of the file.
- **It creates a trap the structure does not enforce.** `tau_a` and
  `synapse_config` must NOT appear in `d`: they are condition-owned, and a
  preset carrying them would be silently overridden by whichever condition is
  applied. That is documented at `srnn_param_preset.m:77-81`, but it is a rule
  you have to know rather than one the shape prevents.

### The alternative TR wants

One helper builds ONE condition; the preset composes them, so it reads
self-contained:

```matlab
conditions = { ...
    adaptation_regime('no_adaptation'), ...
    adaptation_regime('sfa1_std1', 'tau_a', log_ladder(0.25, 10, 1), 'std', std_1), ...
    adaptation_regime('sfa3_std2', 'tau_a', log_ladder(0.25, 10, 3), 'std', std_2) };
```

Physics stays in the preset, boilerplate stays shared. This is the same
independence argument that decided the `log_ladder` change: presets should
restate rather than inherit, because a preset is where a network is described.

### What it costs, and the question to answer first

- **All ten presets get rewritten**, plus the regime-set machinery deleted or
  reduced to condition builders.
- **"`sfa_only` means the same thing everywhere" becomes a CONVENTION rather
  than a mechanism.** Today it is guaranteed by construction. After, it is
  guaranteed by whoever writes the next preset — so the guarantee should move
  into a test that compares regimes of the same name across every preset.
- **The `_oneTS` derivation has to keep working.** If a preset writes out
  `std_1` by hand, it can drift from `std_2`. The builder should derive it
  (`first_timescale_routes`) rather than let the preset restate it — otherwise
  this trades one duplication for a worse one.

That last point is the real design question: how much does the preset spell out,
and how much does the builder still derive? Composing from
`log_ladder(0.25, 10, 1)` and a derived `std_1` is probably the line — the
preset states the ENDPOINTS and the COUNT, the builder handles the shapes.

Related: [[use-c-for-cell-type-count]] is the other naming/structure cleanup
waiting on a deliberate pass.

---

## CLAUDE.md leaves the impression that the paper's STD is on E→E and I→I; it is on all four routes

| | |
|---|---|
| Noted | 2026-09-02 · `refactorRunAll` @ `d6b7fea` · R5611351 · Claude Code (Opus 5), session e22d2fab |
| Raised by | TR, while establishing whether the MC port's `'synaptic'` readout can take a single route |
| Status | **Recorded, not implemented.** Documentation only; the code is correct. |

### What the paper's preset actually does

`srnn_param_preset` builds the routes explicitly for both dualStd presets
(lines 184-188 and 242-246), and they are **identical on all four**:

```matlab
dual_std = struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]);
std_routes.E.E.std = dual_std;
std_routes.E.I.std = dual_std;
std_routes.I.E.std = dual_std;
std_routes.I.I.std = dual_std;
```

Confirmed by resolving the conditions rather than reading the source: the
`std_only` and `sfa_and_std` conditions of
`celltype_pairs_Sc0p2_noise0p025_dualStd_4cond` each carry all four routes at
`tau_rec = [2 4]`, `tau_rel = [0.25 0.5]`. `no_adaptation` and `sfa_only` carry
no `synapse_config` at all.

### Why CLAUDE.md misleads without being wrong

Neither sentence is false on its own, which is what makes this hard to spot:

- **Line 225** — "`synapse_config = []` means *the default E→E and I→I routes*"
  — is **accurate**. Verified: `srnn_adaptation_conditions('SRNNCellTypePairs')`
  with no config returns `E->E:STD  I->I:STD`.
- **Line 224** — "`'celltype_pairs'` is the `SRNNCellTypePairs` preset (E/I with
  STD on E→E and I→I)" — describes a preset that **no longer exists**.
  `srnn_param_preset('celltype_pairs')` now throws
  `srnn_param_preset:RetiredPreset`.

Together they leave a reader believing the Pairs network depresses E→E and I→I.
CLAUDE.md documents the **default** and the **retired** preset, and never says
that every preset the paper actually runs overrides that default. The one thing
a reader most needs — what the paper's network does — is the thing not stated.

### Why it mattered here

It is not cosmetic. The whole question for the memory-capacity port was whether
a neuron has ONE synaptic output or several. Under CLAUDE.md's picture, an E
neuron has STD on E→E but not E→I, so `b·r` is ambiguous and the `'synaptic'`
readout is undefined — the port would need a per-route decision. Under what the
code actually does, all of a neuron's outgoing routes share one `(tau_rec,
tau_rel)`, the `b` trajectories are bit-identical, and reading route 1 is exact.
Opposite conclusions from the same paragraph.

### Suggested fix

Two sentences in the `src/presets/` section: correct line 224 to drop the
retired `'celltype_pairs'`, and add that the dualStd presets supply their own
`std_routes` covering **all four** routes, so the E→E/I→I default applies only
when a preset does not override it. Worth doing alongside the MC port, whose
route-redundancy check is the thing that depends on it.

Related: [[esn-port-to-celltype-pairs]] if that note gets written.

---

## Use `C` for the cell-type count everywhere, and stop calling it `D`

| | |
|---|---|
| Noted | 2026-09-02 · `refactorRunAll` @ `e4a03e5` · R5611351 · Claude Code (Opus 5), session e22d2fab |
| Raised by | TR, after asking whether `C` and `D` were the same quantity |
| Status | **Recorded, not implemented.** Naming only; no behaviour changes. |

### The two names

They are the same number, and equal by construction:

| symbol | layer | property |
|---|---|---|
| `C` | model | `SRNNCellTypePairs.n_cellTypes` |
| `D` | connectivity | `RMTBlocks.n_types`, defined as `numel(f)` |

`SRNNCellTypePairs.build_network` sets one from the other in a single call —
`rmt.set_types(obj.f, obj.mu_tilde, obj.sigma_tilde)` — and `obj.f` is the
`1 x C` fractions row, so `D = C` always. No configuration makes them differ.

Neither is a MATLAB identifier. `C` lives only in comments and CLAUDE.md; `D`
lives only in `RMTBlocks` comments and error messages. So this is a comment and
documentation change, not a rename of any property.

### Why `D` is the one to drop

`D` is **overloaded inside `RMTBlocks` itself**, from two legitimate sources:

- `RMTBlocks.m:26` — "Harris Eq. (6) writes `W = S o (A*D + M)` with `D`
  diagonal." Here `D` is Harris's **diagonal Dale sign matrix**, `n x n`.
- `RMTBlocks.m:5` — "full `D x D` block statistics". Here `D` is the **type
  count**, `C`.

Also at `RMTBlocks.m:241` and `:247` in the Dale sense. One file, one letter,
two meanings, and the Dale one is the paper's notation — so it is the count
that should give way.

### What the change is

1. In `RMTBlocks` comments and error messages, replace the count-`D` with `C`
   (or spell it `n_types`, which is the actual property). Leave the four
   Dale-matrix `D`s alone; they quote Harris and are correct.
2. Same substitution in CLAUDE.md, which currently uses `C` in the
   `SRNNCellTypePairs` section and `D` in the `RMTBlocks` sentence describing
   the same quantity.
3. The still-unwritten `RMTBlocks` doc should be written in `C` from the start,
   with one sentence noting that Harris's `D` is a different object.

### Cost, and the one thing to be careful of

Small — 79 lines in `RMTBlocks.m` contain the letter `D`, but the great
majority are unrelated words, so this is a read-and-judge pass over comments,
not a `sed`. **Do not blanket-substitute**: the Dale-matrix occurrences must
survive, and `set_types`' error messages mention both the dimension and the
paper's form.

Worth doing when the `RMTBlocks` doc is written, so the doc and the comments
agree from day one rather than being reconciled later.

---

## `tau_a` should be the independent property and `n_a` Dependent on it, not both settable

| | |
|---|---|
| Noted | 2026-08-25 · `refactorRunAll` @ `05f6827` · R5611351 · Claude Code (Opus 5), session e22d2fab |
| Raised by | TR, after the `logspace(a,b,1)` finding while planning `n_a = 1` conditions |
| Status | **Recorded, not implemented.** Cleanup of the model classes; no science changes. |

### The asymmetry

The synapse side already does this correctly. The neuron side does not:

| | source of truth | how the count is obtained |
|---|---|---|
| synapse | `tau_rec` vector inside `synapse_config` | **derived**: `n_b = numel(tau_rec)` (`SRNNCellTypePairs.m:1049`) |
| neuron | `n_a` **and** `tau_a`, both public and settable | **duplicated**, and the two must agree |

`validate()` (`:948-953`) enforces `numel(tau_a{q}) == n_a(q)` exactly, and
`complete_type_defaults` (`:819`) papers over the duplication by auto-filling
`tau_a{q}` from `n_a(q)` whenever `tau_a{q}` is empty.

Proposal: `tau_a` becomes the only settable property and

```matlab
n_a = cellfun(@numel, tau_a);        % Dependent, read-only
```

which is exactly what the synapse side already does.

### Why this is the right fix and not just tidiness

The auto-fill is **lossy by construction**: a count cannot say *which*
timescales. That is the actual disease, and the `logspace(log10(0.25),
log10(10), n_a)` bug is only its most visible symptom — at `n_a = 1` MATLAB's
`logspace(a,b,1)` returns `10^b`, so you silently get the **slow 10 s** end.

Special-casing `n_a == 1` to 0.25 would not fix this, it would only change which
arbitrary answer you get: for a single timescale, 0.25, the geometric middle
1.58 and 10 are all defensible, and no count can express the choice. Making
`tau_a` authoritative removes the question by making the caller state the
answer, and it makes a run record self-describing — the saved config then says
*which* timescales, not merely how many.

### Do NOT wrap it in a `neuronal_config` struct

TR raised this as an option, by analogy with `synapse_config`. The analogy does
not hold. `synapse_config` is nested because it is genuinely multi-dimensional:
`C x C` routes, times `{std, stf}`, times several parameters each. The neuron
side is **flat** — one vector per cell type — which the existing `1 x C` cell
`tau_a` already expresses exactly. A wrapper would add a level of nesting with
no gain in expressiveness, and it would break `tau_a_E`, the scalar alias the
tau-sensitivity sweep actually drives.

Revisit only if per-type mechanisms multiply (intrinsic currents, per-type
noise). Today the only per-type quantities are `tau_a` and `c`.

### Cost

- **~8 sites set `n_a` directly** and would have to state `tau_a` instead:
  `srnn_adaptation_conditions.m:105-108`, `fig_adaptation_methods`'s
  `combo_config`, and the `scripts/tests/example_SRNNCellTypePairs*` family.
- **Presets lose the "just say 3" convenience.** Mitigate with a helper —
  `srnn_sfa_timescales(K)` returning the standard ladder — so the default lives
  in one place rather than being inlined at every call site.
- **`SRNNModel2` has the identical shape** (`n_a_E` + `tau_a_E`, auto-filled at
  `:1398`). The classes are duck-typed siblings sharing no implementation, so
  changing only Pairs widens the gap — and memory capacity still runs on
  `SRNNModel2`.
- `figure_settings` and the PSA plotters read `n_a` but never write it, so a
  Dependent read-only property serves them unchanged.

### This is NOT needed to add the `[1,1]` and `[3,1]` conditions

Verified by building: a **condition may carry any model property**, so having
the new conditions state `tau_a` explicitly alongside `n_a` bypasses the
auto-fill entirely, today, with no class change:

```matlab
% note the DOUBLE brace -- struct() with a cell value otherwise makes a struct ARRAY
struct('name','sfa1_std1', 'n_a',[1 0], 'tau_a',{{0.25, []}}, 'synapse_config',sc1)
% builds to:  n_a = [1 0], tau_a{1} = 0.25    (not 10)
```

So this entry is cleanup that can land whenever, and the conditions are not
blocked on it. Sequencing note: it lands better **after** the `c/K` decision (see
the entry above), because if that goes in, `K = numel(tau_a{q})` is then the
same single source of truth.

### Related, found while investigating

`ParamSpaceAnalysis2.condition_field_names()` returns a hardcoded
`{'n_a_E','n_a_I','n_b_E','n_b_I'}` — **`SRNNModel2` names only**. So on
`SRNNCellTypePairs` neither `n_a` nor `synapse_config` is in the ban-list, and
`add_grid_parameter('n_a', ...)` would not hit the guard CLAUDE.md describes as
protecting the conditions from being silently overridden by a grid axis.
Probably inert in practice — a row-valued axis would fail elsewhere — but the
protection is not actually there for the class the paper sweeps. Small, separate
fix, and it becomes moot for `n_a` if `n_a` turns Dependent (a Dependent
property with no setter cannot be a grid axis at all).

---

## SFA scaling should be `c/K`, not `c`, so the adaptation budget is invariant to the number of timescales

| | |
|---|---|
| Noted | 2026-08-25 · `refactorRunAll` @ `10e5fbc` · R5611351 · Claude Code (Opus 5), session e22d2fab |
| Raised by | TR, while planning two new adaptation conditions with `n_a = 1` |
| Status | **Recorded, not implemented.** A change to the dynamics of both model classes; would invalidate every existing preset value and archived run. |

### The proposal

The rate currently subtracts an **unnormalized sum** over SFA timescales:

```
r_i = phi( x_i - c * sum_k a_{i,k} )          % today
r_i = phi( x_i - (c/K) * sum_k a_{i,k} )      % proposed, K = n_a
```

With the division, `c` becomes **the total adaptation budget** and stays
meaningful as `K` changes. Without it, `c` is a per-timescale scaling and the
total scales linearly with `K`.

### Why it matters, exactly

Every adaptation state relaxes to the rate: `da_k/dt = (-a_k + r)/tau_k`, so at
steady state `a_k -> r` for **every** timescale, whatever its `tau`. Therefore

```
sum_k a_k  ->  K * r
```

and the total steady-state adaptation is `c*K*r` today versus `(c/K)*K*r = c*r`
under the proposal — **exactly** K-invariant, not approximately. The tau values
set how fast each component gets there, not where it ends up, so they do not
enter this at all.

### What today's presets do instead

They hand-divide. `celltype_pairs_Sc0p2_noise0p025_dualStd` carries
`c = [0.5/3, 0]` and the sweeps run `n_a = 3`, giving a total budget of 0.5.
That is the proposed normalization applied manually, at one hardcoded `K`. It is
correct only as long as `n_a` is 3 — and `n_a` is a **condition** field, so
nothing stops a condition changing it.

Measured on the dualStd preset (`c = 0.1667`, budget `= c * n_a`):

| `n_a` | auto-filled `tau_a` | total budget today | under `c/K` with `c = 0.5` |
|---|---|---|---|
| 3 | `[0.25  1.581  10]` | **0.500** | 0.500 |
| 1 | `10` | **0.167** | 0.500 |

So a proposed `[1,1]` condition would silently have **one third** the adaptation
of `[3,1]`, confounding "how many timescales" with "how much adaptation". That
is what this note exists to prevent.

### Cost of implementing it

- **Every preset's `c` must lose its hand-division**: `0.5/3 -> 0.5`,
  `0.15/3 -> 0.15`, and so on. Miss one and that network gets 3x weaker
  adaptation with no error.
- **Both model classes**, which share no implementation. In
  `SRNNCellTypePairs`: `:1824` (dynamics), `:1920` (Jacobian's `x_eff`), `:1942`
  (the `a -> x` block) and `:1986` (the route block). In `SRNNModel2`: `:2373`
  and `:3421`, plus the plotting readouts around `:2054-2069`. The Jacobian
  blocks carry `c` too, so a partial change would leave the analytic Jacobian
  disagreeing with the dynamics — which `test_SRNNCellTypePairs`'s
  finite-difference check would catch, and nothing else would.
- **Archived runs become incomparable.** `resolved_defaults` records `c`, not the
  budget, so an old and a new run with the same recorded `c` would mean
  different networks. `same_config` would not notice.
- `n_a = 0` must not divide by zero — the branch is already guarded by
  `n_a > 0`, but the new code has to keep it that way.

### Related

Same investigation surfaced a separate defect: `complete_type_defaults`
(`SRNNCellTypePairs.m:819`, and `SRNNModel2.m:1398`) auto-fills
`tau_a = logspace(log10(0.25), log10(10), n_a)`, and MATLAB's
`logspace(a,b,1)` returns `10^b` — so `n_a = 1` yields the **slow** 10 s end
rather than the fast 0.25 s one. Verified by building. TR wants the fast end.
Both issues have to be settled before the `[1,1]` and `[3,1]` conditions mean
anything, and neither is implemented yet.

The **STD side needs no analogous normalization**, and the reason is not that
its depth is K-invariant — it very much is not. Two things differ:

1. **Nothing is auto-filled or hidden.** There is no `n_b` count to set: `n_b` is
   `numel(tau_rec)` read off the `synapse_config` a condition carries
   (`SRNNCellTypePairs.m:1049`), so the timescales are always written out
   explicitly. Choosing one instead of two is an explicit edit to a route, never
   a silent consequence of a count — which is the whole failure mode above.
2. **Depression enters as a PRODUCT, not a sum**, so there is no budget being
   split to normalize back up. Adding a timescale multiplies rather than
   subdividing, and that is deliberate: with both dualStd pairs sharing
   `tau_rec/tau_rel = 8`, the steady state is `1/(1+8r)^2` against the parent's
   `1/(1+8r)` — 0.086 versus 0.29 at `r = 0.3`. Deeper depression is the reason
   the preset exists.

So changing `n_b` DOES change depression depth substantially. That is intended
and visible, where changing `n_a` changes adaptation strength unintentionally
and invisibly. Going from `tau_rec = [2 4]` to a single `tau_rec = 2` is exactly
the parent preset `celltype_pairs_Sc0p2_noise0p025`, which makes `[.,1]` and
`[.,2]` a clean pairing.

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
>    fixes. See `docs/archive/singleCellTypeRefactor.md` §3c.
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

## ~~Memory capacity is still on `SRNNModel2`~~ — DONE 2026-09-02

| | |
|---|---|
| Noted | 2026-08-22 · `refactorRunAll` @ `40b578a` · R5611351 · Claude Code (Opus 5), session e22d2fab |
| Raised by | TR, at the end of the `refactorRunAll` refactor |
| Re-confirmed | 2026-08-25 · `refactorRunAll` @ `05f6827`. TR asked for this issue again; it was already here, so it was re-checked rather than duplicated. Still open and still accurate: all four blockers below stand, `mc_esn` is still the only `SRNNModel2` preset the paper uses, and both entry points still take a preset. The single-cell-type work (`0134e9b`) did not touch it. |
| **Implemented** | 2026-09-02 · `refactorRunAll` @ `7966fc7`. `SRNN_ESN_reservoir < SRNNCellTypePairs`, `mc_esn` deleted, replaced by `mc_pairs_dualStd`. The paper pipeline is now single-class. |

> ### ✅ Done. Kept for the reasoning, not as a description of the code.
>
> The blockers below were **smaller than they look**, which is the useful part
> to carry forward. Nothing about the MC protocol needed `SRNNModel2`: the ESN
> overrides `build_stimulus` (protected on both classes), writes properties that
> are `protected` SetAccess on both, and calls `integrate` / `build_noise` /
> `get_params` / `dynamics_fast`, all present on both. The genuine coupling was
> **eight lines**.
>
> Two things the note below did not anticipate:
>
> - **The `'synaptic'` readout was the real problem**, not the inheritance.
>   Depression is per *route* on Pairs, so a neuron can transmit a different
>   signal to each target and `b·r` may not exist. Resolved by requiring, per
>   presynaptic type, that a type's outgoing routes carry identical STD **and**
>   STF — then reading route 1, which is exact rather than approximate. See
>   `SRNN_ESN_reservoir.assert_route_redundancy` and
>   `scripts/tests/test_esn_route_redundancy.m`.
> - **The conditions inverted in our favour.** `tau_a` is the settable property
>   on Pairs and `n_a` is Dependent — the reverse of `SRNNModel2` — so the
>   `rmfield_if(preset, {'tau_a_E'})` dance *disappeared* from both entry points
>   rather than needing translation.
>
> `plot_esn_timeseries` was worse than estimated: ~20 E/I-coupled sites, not 5.
>
> **Still open, and deliberately so:** `mc_pairs_dualStd` runs with
> `sigma_u_noise = 0`. The mechanism to turn it on exists (`mc_run_config` now
> selects the integrator from the preset, and takes an explicit override), but
> the protocol has not been re-validated with noise. And the preset's physics is
> a starting point carried over from `mc_esn` — TR expects to tune it.

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
| **Superseded** | 2026-08-21, same session. TR reversed the parking decision: the STF figure **is** being rebuilt as part of `refactorRunAll`, on `SRNNCellTypePairs` behind a new preset matching the old parameters. See `docs/archive/refactorRunAll.md` §6. |

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
