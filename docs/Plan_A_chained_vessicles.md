# Plan A: chained vesicle-pool depression

**Status: plan only. Nothing here has been implemented.** Written 2026-11-10
by the lead agent at the end of the session that produced the decisions, the
test audit and the regression references, so that the work can be picked up
after the conversation is compacted. **Revised 2026-11-11 with Tom's answers
to the questions in section 9.** Still under discussion: facilitation (a
note with the options is at
`docs/notes/Depression_facilitation_model_options_2026-11-11.md`), and how a
MATLAB session's parallel-pool size is set and honoured (item 9). Neither
blocks reading the plan; the second must be settled before workers run.

This is the implementation plan for **block A** of
`../StochasticPlasticDynamicalSystemPaper/docs/final_draft_to_do.md`
("Make 1TS versus MTS apples to apples"). It covers issues FD-004, FD-005,
FD-006, FD-007, FD-015, and the lead-in to FD-009 and FD-014, in
`docs/ISSUES.md`.

(The file name keeps the spelling it was asked for under; "vesicle" is the
word used in the text.)

## 0. How to use this plan

**If you are the lead agent picking this up:**

1. Read, in order: `CLAUDE.md` ("Current work"), block A of the to-do list,
   `docs/ISSUES.md`, the top two entries of `docs/WORKLOG.md`, this file,
   `../StochasticPlasticDynamicalSystemPaper/docs/MATLAB_worker_agents.md`,
   and `../StochasticPlasticDynamicalSystemPaper/adaptation_normalization.md`
   (the equations, C1 to C4).
2. Read the code named in section 3 before briefing anyone. Line numbers
   there are as of commit `fde6e88` and will drift; the function names will
   not.
3. Put the questions in section 9 to Tom. **Do not start until he says go.**
4. Then work through section 7 in order. You orchestrate, review and merge;
   workers build.

**What the plan assumes is true** (all verified on 2026-11-10; re-check the
first two before starting):

- The working tree is clean and `run_all_tests('quick')` passes 53 of 53.
- The three regression references exist and pass on the current code
  (`scripts/tests/fixtures/regression_{rapid,stages,fast}_digest.mat`).
- **The model class has not been touched.** The references describe the code
  before this work, which is what makes them worth having.

## 1. What is being built, in one page

Depression gets a second model. Today a route's depression factor is the
product of independent components. The new model is two pools in series.

For presynaptic neuron $i$ on a route:

$$
\frac{db_{i1}}{dt} = \frac{b_{i2}\,(1 - b_{i1})}{\tau_{rec_1}} - \frac{b_{i1}\, r_i}{\tau_{rel}} \tag{C1}
$$

$$
\frac{db_{i2}}{dt} = \frac{1 - b_{i2}}{\tau_{rec_2}} - \frac{1}{n_2}\,\frac{b_{i2}\,(1 - b_{i1})}{\tau_{rec_1}} \tag{C2}
$$

$$
\theta_i = r_i\, b_{i1} \tag{C3}
$$

$b_{i1}$ is the releasable pool, $b_{i2}$ the reserve, $n_2$ the reserve's
capacity relative to the releasable pool. Both start at 1. Both stay in
$[0, 1]$: at $b_1 = 0$ its derivative is $b_2/\tau_{rec_1} \ge 0$, at
$b_1 = 1$ it is $-r/\tau_{rel} \le 0$, and likewise for $b_2$.

**Decisions already made** (to-do list, block A; do not reopen without Tom):

| Decision | Value |
|---|---|
| Both models supported | The chain is an opt-in per route. The product stays the default and **bit-identical**. The geometric-mean option stays. |
| 1TS | A route with one depression component: the current single-pool equation, one state per presynaptic neuron. **No reserve variable exists**, so no unused state and no zero eigenvalue. |
| MTS | The same releasable-pool equation plus a dynamic reserve: two states per route, the same count as today's two-timescale product. |
| Reserve size $n_2$ | 2 |
| $\tau_{rec_1}$, $\tau_{rel}$ | 2 s, 0.25 s (unchanged, so (C1) is identical in 1TS and MTS) |
| $\tau_{rec_2}$ | **5 s** (Tom, 2026-11-11; the 10 s first proposed is dropped) |
| Routes | all four carry depression, as now |
| Strength | 1TS and MTS are **not** matched at steady state, by decision. They coincide at the onset of activity (reserve full). No frozen-reserve control (FD-008 dropped). |
| Sustained output | Rises to a ceiling under the chain, where the product's falls with rate. Accepted. |
| Per-neuron spread | Log-normal, $\sigma = 0.25$, on the depression time constants and on $n_2$. Independent draws per presynaptic neuron, shared by all its outgoing routes. **Allowed on all three kinds of route: chain, one-component, and multi-component product** (Tom, 2026-11-11). 1TS is jittered with the same $\tau_{rec_1}$ and $\tau_{rel}$ draws as MTS. On a chain route the reserve refill is kept slower than the releasable-pool refill, with a minimum separation. **Default is no spread, so existing presets run bit for bit as today.** |
| `std_zero_floor` | Refused on chain routes. |
| MTS chain condition | Named **`sfa3_stdc2`** (Tom, 2026-11-11). |
| Conditions analysed | **Four**: no adaptation, 1TS, chained MTS, product MTS. Tom and Brian want the chain and the product side by side during blocks A to D; the final paper will probably show three. |
| Figure sets | **Two from every run**: a three-condition set (no adaptation, 1TS, chained MTS) and a four-condition set with product MTS as a fourth column or colour. In the three-condition styling the chained MTS takes over the colour product MTS has now; product MTS gets a new colour (green). |
| Reruns | **Agents never run `medium` or larger.** Only Tom starts those and only Tom says when they are done. The plan ends with one `fast` run. |
| Parallel pools | Automatic pool creation is **off** machine-wide (Tom, 2026-11-11). Up to 12 workers in total across the lead and all workers. Tom keeps one MATLAB session with a one-worker pool on the `Just1` profile to hold a toolbox licence seat; it does not count against the 12. |

Steady state at constant rate, with $k_1 = r\tau_{rec_1}/\tau_{rel}$,
$k_2 = r\tau_{rec_2}/(n_2\tau_{rel})$, $s = 1 + k_1 + k_2$:

$$
b_1 = \frac{2}{s + \sqrt{s^2 - 4k_2}}, \qquad b_2 = 1 - k_2 b_1 \tag{C4}
$$

(This form avoids cancellation at small $k_2$ and reduces to $1/(1+k_1)$ at
$k_2 = 0$.) With the chosen parameters ($\tau_{rec_2} = 5$ s), $k_1 = 8r$ and
$k_2 = 10r$: $b_1$ is 0.420, 0.200, 0.106, 0.054 at $r$ = 0.1, 0.25, 0.5,
1.0 (1TS: 0.556, 0.333, 0.200, 0.111; current product: 0.309, 0.111, 0.040,
0.012), and the sustained output $r\,b_1$ rises to a ceiling of $1/18$.

## 2. The one rule that shapes the implementation: branch, do not generalise

The acceptance test for "the product stays the default and bit-identical" is
the ultrarapid regression, and it compares **bytes**. Reordering a floating
point expression changes the last bit and fails it.

So at every site, **leave the existing product code textually as it is and
add a branch for chain routes**:

```matlab
if params.std_chain(pre, post)
    % new chain code
else
    % the existing lines, untouched
end
```

Do not rewrite the product path into a generic helper that both models call,
however tidy that looks. A shared helper is fine only where the existing code
already goes through one (`combine_std_pair`, `apply_std_zero_floor_pair`).

If a refactor of the product path is ever wanted, it is a separate commit,
justified on its own, and it replaces the references deliberately.

## 3. The code as it stands

All in `src/model/SRNNCellTypePairs.m` unless stated. Depression state is
stored **per route** `(pre, post)` as an `npre x nb` matrix, one row per
presynaptic neuron, one column per component, at `layout.b{pre, post}`,
initialised to ones.

| Function (line at `fde6e88`) | What it does with depression today | What the chain needs |
|---|---|---|
| `compile_synapse_config` (1349) | Validates `std.tau_rec` / `tau_rel` / `combine`; fills `pair.n_b`, `tau_b_rec`, `tau_b_rel`, `std_geomean`, `scale` | New `std.model` and `std.n2`; fill `pair.std_chain`, `pair.std_n2` |
| `get_params` (854; STD fields near 906) | Copies the pair fields into `params` | Copy the new fields; later, the per-neuron matrices |
| `validate` (1173) | Cross-field checks | Refuse `std_zero_floor` with any chain route; validate the spread |
| `make_state_layout` (2107), `initialize_state` (2137) | `nb` columns per route, ones | **No change.** A chain route has `n_b = 2`: column 1 is $b_1$, column 2 is $b_2$ |
| `dynamics_fast` (2259) | `db = (1-b)./tau_rec - b.*r./tau_rel`, vectorised over columns; `depression = combine_std_pair(prod(b,2), ...)` | Chain branch: (C1), (C2); `depression = b(:,1)` |
| `combine_std_pair` (2364) | Takes `P = prod(b)`; product or geometric mean; returns `dD/dP` | Unchanged; not used on chain routes |
| `apply_std_zero_floor_pair` (2345), `std_zero_floor_pmin` (2354) | Product-specific floor | Unchanged; never reached for a chain route (refused in `validate`) |
| `compute_Jacobian_fast` (2390) | `(b,b)` block diagonal; `(b,a)` and `(b,x)` from `-b/tau_rel`; `(x,b)` from the product rule | Chain branch: a 2x2 coupling per neuron; only $b_1$ depends on rate; only $b_1$ feeds $x$ |
| `jacobian_times` (2582) | The same blocks applied matrix-free. **Must equal `compute_Jacobian_fast * Y` to round-off** (its header states the contract; `test_jacobian_times` enforces it) | Mirror the chain branch exactly |
| `unpack_and_compute_states` (3120) | `synaptic_output = facilitation .* depression .* rate`, depression from `prod(b)` | Chain branch: `b(:,1)` |
| `routes_identical` (2804), `route_fields` (2831) | Route signature for comparing routes | Add model, `n2`, and the spread |
| `finite_difference_jacobian` (2852) | Differences `dynamics_fast` | **No change**; it is the chain's independent check |
| `transient_gain` (2873) | Uses `jacobian_times` and the dense `J_xx` | **No change**; inherits |
| `lya_summary` (748) | Leading-vector fraction on `layout.b` counted as "STD" | **No change**; both pools count as depression |
| `plot_celltypes` (1829) | Draws `prod(b)` collapsed across routes | Display: draw the route's factor |
| `build_tau_a_matrix` (1551), `tau_a_spread` / `tau_a_seed` / `tau_a_clipped` (properties near 92, 97, 295) | **The pattern to copy for the spread**: own seeded stream, RNG state saved and restored, minimum log-separation floor, clip count recorded | A sibling for depression |

Outside the model class:

| File | Why it matters |
|---|---|
| `src/model/SRNN_ESN_reservoir.m`: `assert_route_redundancy` (943), `route_signature` (996), `plot_esn_timeseries` (488, draws `prod(bq, 2)` near 674) | The synaptic readout needs a neuron's routes identical; the signature must see the model, `n2` and the spread |
| `src/figures/fig_STD_steady_state.m`, `fig_STD_strength_matching.m`, `fig_example_timeseries.m`, `fig_main1_concepts.m`, `src/figures/helpers/write_manuscript_tables.m`, `src/plotting/plot_SRNN_combined_tseries.m`, `src/analysis/run_eig_heatmap.m`, `src/analysis/run_tau_sensitivity_analysis.m` | Each reads `tau_rec` / `tau_rel` / `n_b_pairs` or recomputes `prod(b)` itself |
| `src/presets/srnn_condition_titles.m`, `full_adaptation_condition.m`, `src/figures/helpers/manuscript_style.m`, `mc_condition_colors.m`, `style_main4_grouped.m`, `style_main5_grouped.m`, `src/util/write_run_parameters_md.m`, and about a dozen `fig_*.m` | Keyed on, or mentioning, the condition name `sfa3_std2` |
| `docs/EquationsParametersDocs/Equations_stability_paper.md` | **The one human-written statement of the equations.** The chain goes here, once. Do not add another copy anywhere. |
| `src/model/SRNNModel2.m`, `src/model/jacobian/compute_J_eff.m` | **Leave alone.** Not used by the paper. |

## 4. Design

### 4.1 Route configuration

```matlab
% 1TS: unchanged, the current single pool
synapse_config.E.E.std = struct('tau_rec', 2, 'tau_rel', 0.25);

% MTS: the chain
synapse_config.E.E.std = struct('model', 'chain', ...
    'tau_rec', [2 5], ...    % [releasable-pool refill, reserve refill]
    'tau_rel', 0.25, ...     % ONE release time: only the releasable pool releases
    'n2', 2);                % reserve capacity / releasable-pool capacity
```

- `std.model`: `'product'` (default, and what an absent field means) or
  `'chain'`.
- Chain routes require exactly two `tau_rec` values, a scalar `tau_rel`, and a
  positive finite `n2`. No default for `n2`: a chain route states it.
- `std.combine` is an error on a chain route.
- `route.scale` stays allowed on a chain route (it multiplies the weights and
  is folded into `params.W` in `get_params`, independent of the model).
- Facilitation (`route.stf`) on a chain route: **refused in this pass**,
  with a clear error, unless Tom decides otherwise (section 9, item 4, has
  the sketch of a two-variable facilitation for later). The paper uses none.
- `params.std_chain` (`C x C` logical) and `params.std_n2` (`C x C`, NaN where
  unused). `n_b_pairs` is 2 for a chain route.

### 4.2 The chain's linearisation

With $u = 1/\tau_{rec_1}$, $v = 1/\tau_{rec_2}$, $w = 1/\tau_{rel}$,
$q = 1/n_2$, per neuron:

| | $\partial/\partial b_1$ | $\partial/\partial b_2$ | $\partial/\partial r$ |
|---|---|---|---|
| $\dot b_1$ | $-b_2 u - r w$ | $(1 - b_1)\,u$ | $-b_1 w$ |
| $\dot b_2$ | $q\, b_2 u$ | $-v - q\,(1 - b_1)\,u$ | $0$ |

and the depression factor $D = b_1$, so $\partial D/\partial b_1 = 1$,
$\partial D/\partial b_2 = 0$.

How these enter the blocks, in the notation of `compute_Jacobian_fast`
($\phi'$ the activation derivative, `c_eff` the SFA strength per timescale):

- `(b, b)`: the 2x2 above for each neuron. Rows and columns are
  timescale-major, like `layout.b`: entry for neuron $i$, components $(m, k)$
  sits at row `(m-1)*npre + i`, column `(k-1)*npre + i`.
- `(b, x)`: $(\partial \dot b_m/\partial r)\,\phi'_i$ on the presynaptic
  neuron's own $x$. Nonzero for $b_1$ only.
- `(b, a)`: $-(\partial \dot b_m/\partial r)\, c_{\mathrm{eff}}\,\phi'_i$ for
  each SFA timescale of that neuron. Nonzero for $b_1$ only.
- `(x, b)`: `W_block * diag(facilitation .* pre_rate .* dD_db_m) / tau_d`.
  Nonzero for $b_1$ only.
- `(x, x)` and `(x, a)`: unchanged in form; they use the depression factor,
  which is now $b_1$.

`jacobian_times` applies the same four blocks to the slices of `Y`. Write
the chain branch of both functions from this table, not one from the other.

### 4.3 Per-neuron spread (FD-015)

Follow the SFA spread (`tau_a_spread`, `build_tau_a_matrix`) closely.

- New properties, proposed names: `std_spread` (`1 x C`, by presynaptic cell
  type, default zeros), `std_seed` (default derived from `rng_seeds(1)` with
  its own offset, well away from the streams that build `W`, the setpoints
  and the SFA ladders), and read-only `std_draws` / `std_clipped` for what
  was realised.
- Draw `z` from that stream as `n x Q` standard normals, **always all Q
  columns for every neuron**, whatever its routes use. The columns are in a
  fixed order that does not depend on the network:

  | Column | 1 | 2 | 3 | 4 | 5 | 6, 7, ... |
  |---|---|---|---|---|---|---|
  | Quantity | $\tau_{rec_1}$ | $\tau_{rel_1}$ | $\tau_{rec_2}$ | $\tau_{rel_2}$ | $n_2$ | $\tau_{rec_3}$, $\tau_{rel_3}$, ... |

  A value is `nominal * exp(sigma * z)`. A one-component route uses columns
  1 and 2; a chain route 1, 2, 3 and 5; a two-component product route 1 to
  4. Because `randn(n, Q)` fills column by column, the first columns are the
  same numbers whatever Q is. So on one network seed the 1TS condition, the
  chained MTS and the product MTS all have the **same** $\tau_{rec_1}$ and
  $\tau_{rel_1}$ per neuron, and the two MTS conditions the same
  $\tau_{rec_2}$. That is what keeps (C1) identical between 1TS and MTS and
  makes the four conditions comparable neuron by neuron.
- One draw per presynaptic neuron, used on every route that neuron has.
- Ordering, **chain routes only**: where a neuron's
  $\log(\tau_{rec_2}/\tau_{rec_1})$ falls below a quarter of the nominal
  log-span, push the two apart symmetrically about their log-midpoint, and
  count it in `std_clipped` (the SFA rule). With medians 2 s and 5 s and
  $\sigma = 0.25$ the nominal log-span is 0.92, so this will trigger for a
  small share of neurons; report the count. A product route needs no
  ordering: the product is symmetric in its components.
- Save and restore the RNG state around the draw.
- `get_params` then carries per-neuron values for the routes of a type whose
  spread is nonzero. **A zero spread must leave `params` exactly as it is
  today**: do not convert the shared rows to matrices unless the ultrarapid
  regression still passes byte for byte. If it does not, keep the row path
  and branch on whether the type has a spread.
- With the spread on, the product path's three sites (`dynamics_fast`,
  `compute_Jacobian_fast`, `jacobian_times`) must accept `npre x nb`
  matrices of time constants where they take `1 x nb` rows today. In
  `dynamics_fast` the existing division already broadcasts. In the two
  Jacobian functions the rows are expanded with `kron(..., ones(npre, 1))`
  and `(1 ./ tau)'`; those need a per-neuron form, behind a branch.

### 4.4 Presets, conditions and the bundle (FD-007)

- A new preset family on the manuscript's baseline network
  (`..._3cond_mu7revisedMedium` is the model to copy; write the main preset
  out in full, chain its twins within the bundle: the noise-free twin for
  memory capacity, the step-input twin for local Lyapunov, the single-neuron
  preset). **Four conditions**, in this order:

  | Name | SFA | Depression | Role |
  |---|---|---|---|
  | `no_adaptation` | none | none | as now |
  | `sfa1_std1` | one timescale | single pool `(2, 0.25)` | 1TS |
  | `sfa3_stdc2` | three timescales | chain: `tau_rec [2 5]`, `tau_rel 0.25`, `n2 2` | chained MTS, the paper's new MTS |
  | `sfa3_std2` | three timescales | product: `tau_rec [2 4]`, `tau_rel [0.25 0.5]` | product MTS, kept for the side-by-side |

  All four routes carry the depression in each adapting condition.
- **The spread is set in the preset, not in a condition** (Tom,
  2026-11-11): `std_spread = [0.25 0.25]`, exactly as `tau_a_spread` is set
  today, with the condition names unchanged. All three adapting conditions
  are then jittered with shared draws, and the preset's name says so
  (`..._stdSpread0p25_...`).
- Every map keyed on condition names must learn `sfa3_stdc2`: titles
  (`srnn_condition_titles`), colours and short labels (`manuscript_style`,
  `mc_condition_colors`), the grouped figure styles, and
  `full_adaptation_condition`, which decides the most-adapted regime by two
  exact regular expressions (read its header on why optional groups cannot
  be used). With four conditions it must return the chained MTS, not the
  product MTS.
- **Colours and titles.** `sfa3_stdc2` takes the colour and the title
  ("Multiple-Timescale Adaptation") that `sfa3_std2` has now. `sfa3_std2`
  gets a new colour (green) and a title that says what it is (for example
  "Multiple-Timescale Adaptation (product)"). Runs made before this change
  have only `sfa3_std2`; regenerating their figures will show it in the new
  colour, which is correct but different from the manuscript's current
  figures. Say so in the commit.
- **Two figure sets from one run.** `make_all_paper_figures` (or the config)
  gains a way to name the conditions to draw, and the bundle's run script
  calls it twice: once for `no_adaptation`, `sfa1_std1`, `sfa3_stdc2` into
  `figs/<bundle>_3cond/`, once for all four into `figs/<bundle>_4cond/`.
  Each figure function must lay out correctly for three or four conditions.
  The step 2 inventory lists which figures assume three.
- A config and run script pair in `scripts/paper/`, copied from
  `mu7revisedMedium_config.m` as `regression_fast_config.m` was.
- Adding presets fails `test_preset_golden` by design. Procedure: confirm the
  only failure is "valid-preset list unchanged", delete the fixture, run
  `make_preset_golden`, force-add it, and say in the commit why. Also add the
  names to `make_preset_golden.m` and `test_preset_conditions.m`, which list
  presets explicitly.

### 4.5 Display and reservoir code (FD-006)

- `SRNN_ESN_reservoir.route_signature` and `SRNNCellTypePairs.route_fields`
  gain the model, `n2` and the spread, so the synaptic readout's redundancy
  check sees them. With one draw per presynaptic neuron shared across its
  routes, a preset with the chain on all four routes still passes it.
- Wherever display code computes `prod(b)` itself, it should take the
  route's factor from one public static on the model class (proposed:
  `SRNNCellTypePairs.std_route_factor(b, params, pre, post)`, which returns
  `b(:,1)` for a chain route and goes through `combine_std_pair` otherwise).
  That static is added in step 1 so the display work can use it.
- `fig_STD_steady_state` and `write_manuscript_tables` state the steady
  state and the equations; both need the chain's (C4 and C1 to C3).

## 5. Tests

Written with the code they cover, by the worker who writes that code.

**`test_std_chain_model.m`** (FD-005), on a small network with a chain route,
a one-component route and an undepressed route:

1. `dynamics_fast` equals a per-neuron loop evaluation of (C1) and (C2),
   written independently of the vectorised code.
2. Steady state under a constant rate equals (C4), for several rates.
3. Stored synaptic output equals rate times $b_1$ exactly.
4. Analytic Jacobian against `finite_difference_jacobian`, at random interior
   states and at states from a trajectory.
5. `jacobian_times` against the assembled Jacobian (relative error at
   round-off), for `K` = 1, several, and `Y = eye(N)`.
6. A very large reserve (`n2 = 1e12`) reproduces the single pool: $b_1$'s
   trajectory against a one-component route, to a stated tolerance; and
   $b_2$ stays at 1.
7. $b_1$, $b_2$ stay in $[0, 1]$ over a run with large rate swings.
8. A chain route with `std_zero_floor` on, with `combine`, with one or three
   `tau_rec` values, with a vector `tau_rel`, or without `n2`, each raises
   its own error.
9. All three Lyapunov estimators run on a chain network and agree as well as
   they do on a product network of the same size.

**`test_std_spread.m`** (FD-015):

1. Zero spread: `params` and a run are bit-identical to a model that never
   mentions the spread.
2. Same seed, same draws; another seed, other draws; `W`, setpoints, SFA
   ladders, the stimulus and `S0` are untouched by turning the spread on.
3. 1TS and MTS built on one network seed have identical $\tau_{rec_1}$ and
   $\tau_{rel}$ per neuron.
4. A neuron's routes carry identical values (so `routes_identical` and the
   reservoir's redundancy check still pass).
5. All drawn values positive; log-mean and log-spread match the nominal and
   $\sigma$ within sampling error at large `n`; the ordering floor holds and
   `std_clipped` counts what it moved.
6. Analytic Jacobian against finite differences, and `jacobian_times`
   against the assembled Jacobian, with the spread on, **on each of the
   three kinds of route**: chain, one-component, two-component product.
7. On one network seed, the product MTS and the chained MTS have the same
   $\tau_{rec_1}$, $\tau_{rel_1}$ and $\tau_{rec_2}$ per neuron.

Both go in the quick suite (each well under 30 s). Existing tests that must
keep passing untouched: `test_std_product_model` (frozen trajectory
checksums of the product model), `test_jacobian_times`, `test_std_geomean`,
`test_SRNNCellTypePairs`, `test_esn_route_redundancy`,
`test_pairs_reproducibility`.

## 6. Gates

Run by the lead on merged `main`, in the lead's MATLAB, after each merge.

| Gate | Command | When | Must show |
|---|---|---|---|
| Ultrarapid | `run_all_tests('test_regression_rapid')` | after every merge that touches `src/model/` or `src/analysis/ParamSpaceAnalysis2.m` | PASS, exact. **Any difference is a defect in the merge**, not a reason to regenerate. |
| Quick suite | `run_all_tests('quick')` | after every merge | all pass |
| Stages | `run_all_tests('test_regression_stages')` | after the display and stage work merges | PASS |
| Full | `regression_check()` launched through a timer, outcome read from `data/regression_fast/<tag>.result.txt` (about 50 min; the MCP client drops a call silent for 30 min) | once, when steps 1 to 4 are all merged | PASS |

**The existing references are never regenerated in this work.** The chain
adds new behaviour on new routes; it must not move a single existing number.
Coverage of the chain itself goes into a **new** reference (step 5).

After a classdef edit the lead's MATLAB needs `clear classes` (it wipes the
base workspace; say so). `ParamSpaceAnalysis2.srnn_property_info` caches the
model's property list and reports a new property as "not a property" until
then.

## 7. Work plan: who builds what, and in what order

### 7.1 Roles

- **Lead (Fable 5.2, this session's successor):** orchestrates. Writes the
  briefs, reviews every diff, runs the gates, merges, keeps `ISSUES.md` and
  `WORKLOG.md`, talks to Tom. Does not build, so its context stays on the
  architecture.
- **Fable 5.2 workers, medium effort, at most 2 at a time**
  (`matlab-fable-1`, `-2`, `-3`): anything that re-implements a class, a
  numerical method, or equations in code. Here: the model class, the
  reservoir class, their tests.
- **Opus 5.5 workers, medium effort, at most 3 at a time**
  (`matlab-opus-1`, `-2`, `-3`): the straightforward work. Here: presets,
  configs and run scripts, condition-name maps, figures and tables.

Rules from `MATLAB_worker_agents.md` that bind this plan:

- **One instance of each named agent at a time.** For N parallel jobs use N
  different agents.
- **Every brief states `pool: N`.** Pool workers summed over the lead's
  MATLAB and every worker's stay at or below 12. Tom's own one-worker
  `Just1` session is outside that count. The lead closes its own pool
  before handing workers out.
- **A worker first verifies its session is its own** (its MATLAB process id
  differs from the lead's, which the brief states; its workspace is empty;
  no pool is up), then opens its pool. A worker running plain unit tests
  can take `pool: 0`; one running sweeps or stage functions needs a pool.
- **Automatic pool creation is off** (Tom, 2026-11-11), so a `parfor`
  reached with no pool runs serially in the client instead of opening the
  default pool. The `pool: 0` gap described in `MATLAB_worker_agents.md` is
  closed; a pool exists only where someone opened one.
- **Name the pool profile.** See section 9, item 9: the default profile is
  currently `Just1` (one worker), so `parpool(parallel.defaultProfile, N)`
  fails for N above 1. Until that is settled, open pools with
  `parpool('Processes', N)`.
- A worker's MATLAB is cold, private, and gone when the worker finishes.
  What the lead needs must be on disk or in the report.
- The worker agents load only in a session started after
  `.claude/agents/` existed. **If they are not listed, restart the session**
  before step 1.

### 7.2 Isolation and merging

Workers that run at the same time each get their own git worktree and
branch (the Agent tool's `isolation: "worktree"`), so one worker's
half-edited class cannot break another's tests. Each worker runs
`setup_paths` from its own worktree root; `data/` and `figs/` are per
worktree and gitignored; the tracked fixtures are present in every worktree.

Workers commit on their branch (stage by name, no AI attribution, messages
that record reasoning and wrong turns). The lead merges into `main`, then
runs the gates. A worker never merges, never regenerates a regression
reference, and never touches a file outside its brief.

Files are assigned so that no two concurrent workers edit the same file:

| File or area | Owner |
|---|---|
| `src/model/SRNNCellTypePairs.m` | one Fable worker at a time, never two |
| `src/model/SRNN_ESN_reservoir.m` | the second Fable worker |
| `src/presets/*`, `scripts/paper/*`, the preset fixture | one Opus worker |
| `src/figures/*`, `src/plotting/*`, `src/util/write_run_parameters_md.m` | Opus workers, split by file |
| `scripts/tests/test_std_chain_model.m`, `test_std_spread.m` | the Fable worker who writes the code under test |
| `docs/*`, `CLAUDE.md`, merges | the lead |

### 7.3 Sequence

The model class is one file and its three depression sites must agree with
each other, so **steps 1 and 3 cannot be split or run in parallel**. Work
that does not touch that file runs beside them.

**Step 0. Pre-flight (lead).** Clean tree; worker agents listed; the lead's
pool closed; quick suite and ultrarapid tier pass; Tom's answers to
section 9 written into the to-do list.

**Step 1. The chain in the model class: FD-004 and FD-005.**
One Fable worker (`matlab-fable-1`), alone, `pool: 4` (the quick suite
reaches `parfor`).

- Edits: `src/model/SRNNCellTypePairs.m` and the new
  `scripts/tests/test_std_chain_model.m`. Nothing else.
- Builds: the route configuration and validation; the chain branch in
  `dynamics_fast`, `compute_Jacobian_fast`, `jacobian_times` and
  `unpack_and_compute_states`; `route_fields`; the public
  `std_route_factor` static; the refusals.
- Suggested commit order, each leaving the suite green: (a) configuration,
  validation and params, with no dynamics yet (a chain route errors as "not
  implemented"); (b) `dynamics_fast` and the tests of the equations and
  steady state; (c) `compute_Jacobian_fast` and the finite-difference test;
  (d) `jacobian_times` and the equality test; (e) unpacking, the factor
  static, the remaining tests.
- Acceptance: `test_std_chain_model` passes; `run_all_tests('quick')` passes
  on its branch; **`test_regression_rapid` passes exactly**.
- Lead: review the diff against section 2 (are the product lines
  untouched?) and section 4.2, merge, run the ultrarapid and quick gates.

**Step 2. In parallel with step 1, on files step 1 does not touch.**

- `matlab-opus-1`, `pool: 0`, read-only: an **inventory**, written to
  `docs/notes/`, of every place outside the model class that (i) reads
  `tau_rec` / `tau_rel` / `n_b_pairs`, (ii) computes `prod(b)`, or (iii)
  keys on or names `sfa3_std2`, or (iv) **assumes exactly three
  conditions** (a fixed number of columns, a hard-coded list of names, a
  three-colour palette indexed by position), each classified as code that
  must change, a name map that must gain an entry, or a comment. Section 3
  lists the files found by search on 2026-11-10; this turns that into
  line-level work for step 4, including what each figure needs in order to
  build for three or for four conditions.
- Nothing else can usefully start before step 1 merges: presets need the
  chain to validate, display code needs the factor static.

**Step 3. After step 1 is merged: three jobs at once.**

- `matlab-fable-1` (`pool: 2`): **the per-neuron spread, FD-015**, in
  `src/model/SRNNCellTypePairs.m` and `scripts/tests/test_std_spread.m`.
  Acceptance as in step 1, plus `test_std_spread`.
- `matlab-fable-2` (`pool: 2`): **the reservoir class**,
  `src/model/SRNN_ESN_reservoir.m` (`route_signature`,
  `assert_route_redundancy`, the time-series plot's depression trace) and
  the chain cases in `scripts/tests/test_esn_route_redundancy.m`. It works
  against the property and field names fixed in section 4, so it does not
  wait for FD-015; the lead reconciles the two at merge.
- `matlab-opus-2` (`pool: 2`): **presets and bundle, FD-007 without the
  spread line**: the preset family, the condition-name maps
  (`srnn_condition_titles`, `full_adaptation_condition`, `manuscript_style`,
  `mc_condition_colors`, the grouped styles), the config and run script,
  the preset lists in the two tests, and the fixture regeneration.
  Acceptance: preset tests pass; a chain model builds and runs from the
  preset; the bundle's config resolves.

Merge order: fable-1 (FD-015), then fable-2, then opus-2. After all three,
the lead adds `std_spread = [0.25 0.25]` to the new preset (one line),
regenerates the preset fixture again if needed, and runs the gates.

**Step 4. After step 3 is merged: display code, FD-006.** Up to three Opus
workers, files split by the step 2 inventory, for example:

- `matlab-opus-1`: `fig_STD_steady_state`, `fig_STD_strength_matching`,
  `write_manuscript_tables` (equations and steady state for the chain).
- `matlab-opus-2`: `fig_example_timeseries`, `plot_SRNN_combined_tseries`,
  `fig_main1_concepts`, and `plot_celltypes` if the lead assigns that one
  method of the model class to an Opus worker; otherwise it goes to a Fable
  worker as a small follow-up.
- `matlab-opus-3`: the figures that name `sfa3_std2`, so they also work on
  the chain condition, and `write_run_parameters_md`.

Across all three: the option that names which conditions a figure draws
(carried by the config through `make_all_paper_figures`), and a layout that
is right for three conditions and for four. The lead fixes that option's
name and where it is read before briefing, so the three workers agree.

Each: `pool: 2`. Acceptance: **every figure builds in both sets**,
three-condition and four-condition, on the new preset, and still builds on
an existing three-condition product preset; `test_regression_stages`
passes. Lead: stages gate, then the full gate.

**Step 5. A reference for the chain itself (lead, with one Fable worker).**
Once the chain passes its unit tests and the gates: a **new** ultrarapid
reference covering a chain network (trajectories, Jacobians, the three
estimators, a sweep, memory capacity, with and without the spread), saved
as its own fixture beside the existing three, with its own check and test.
The existing references are not edited. Two runs first, to confirm it
reproduces.

**Step 6. The four-condition bundle at `fast`, and stop.** In the lead's
desktop MATLAB with the whole machine (pool opened by name, 12 workers):
run the new bundle at `fast` on all four conditions (the manuscript-size
`fast` pipeline took about 67 min on 2026-09-15 with three), then build
both figure sets, three-condition and four-condition. Check every stage
and every figure in both sets, report to Tom, and stop. **No agent runs
`medium` or larger**: those are Tom's to start and Tom's to call done
(FD-009). If the joint sample is to be seeded (section 9, item 5), that
goes in before this run so later runs sample the same points.

**Step 7. The comparison figure, FD-014** (one Opus worker): the two
conditions' depression under step and sinusoidal rates and their
steady-state curves, from
`scripts/explorations/adaptation_normalization/demo_chain_1TS_vs_MTS.m`,
with the chosen parameters. It must show what Brian accepted: the first-cycle
match and the steady-state difference.

**Not in this plan:** the manuscript text (Methods Eqs. 5 and 6, the
Supplemental Methods, the conditions and timescale tables, the paragraph on
steady-state depth). Those are tracked in the to-do list and are done with
Tom once the runs exist. The single equations file
`docs/EquationsParametersDocs/Equations_stability_paper.md` **is** in scope:
the lead updates it when step 1 merges.

### 7.4 What every brief must contain

Workers start with no memory of this conversation. Each brief gives:

1. `pool: N` on the first line.
2. The goal in two sentences, and which issue it serves.
3. What to read first: `CLAUDE.md`, this plan (the sections that apply),
   `adaptation_normalization.md` for the equations.
4. The files it may edit, and that nothing else may be touched.
5. The rule of section 2, where it applies, in full.
6. The acceptance commands, and that a difference in
   `test_regression_rapid` is a defect to find, never a reason to regenerate.
7. Commit rules: on its own branch, stage by name, no AI attribution, the
   message carries the reasoning and what turned out wrong.
8. What to report: what was built, what was run and what it showed, what
   failed with the error text, anything it was unsure of, defects it noticed
   elsewhere (reported, not fixed), and that its pool is closed.

## 8. Risks and traps

- **Bit-identity.** The reason for section 2. If the ultrarapid gate fails
  after a merge, bisect the merge; do not adjust a tolerance.
- **`jacobian_times` drift.** It has no independent check except equality
  with the assembled Jacobian, which in turn rests on finite differences.
  Both tests must cover chain routes, with and without the spread.
- **Stale class cache.** After a classdef edit, `clear classes`; a new
  property otherwise reports as nonexistent from `ParamSpaceAnalysis2`.
- **The condition-name trap.** `full_adaptation_condition` uses two exact
  regular expressions because MATLAB returns no token for a trailing
  optional group. A new name needs its own pattern.
- **Conditions must set the same fields.** `validate_preset_conditions`
  requires every condition of a preset to set the same field set. The
  spread therefore belongs at preset level, not in a condition.
- **The readout.** `SRNN_ESN_reservoir`'s `'synaptic'` readout requires a
  presynaptic neuron's routes to be identical. Per-route draws would break
  it; per-neuron draws shared across routes keep it.
- **Worktrees and the shared machine.** Each worker has its own MATLAB but
  they share the cores and the RAM. `run_memory_capacity` sizes itself from
  free RAM as if alone; do not run two memory-heavy jobs side by side.
- **The default pool profile is `Just1`.** Checked 2026-11-11 from the
  lead's MATLAB: `parallel.defaultProfile` is `Just1`, a one-worker profile.
  Code that opens the default profile now misbehaves: `restart_parpool`
  (called by `run_all_analyses` and `run_all_paper_analyses`) opens a
  one-worker pool, so sweeps run on one worker; and
  `run_lyapunov_spectrum`, `run_transient_gain`,
  `run_numerics_verification` and `wait_for_parpool` call
  `parpool(parallel.defaultProfile, n)` with n up to 12, which errors. The
  six worker-agent prompts use the same call. Until section 9, item 9 is
  settled: open a `Processes` pool by name before running a stage, and do
  not run the full bundle.
- **No pool means serial.** With automatic creation off, forgetting to
  open a pool no longer oversubscribes the machine; it runs slowly on one
  core instead. A job that seems stuck may simply have no pool.
- **The 30-minute limit.** A MATLAB call silent for 30 minutes is dropped by
  the client though MATLAB carries on. Long runs are launched through a
  timer and followed by reading a file; multi-hour runs belong in the lead's
  desktop MATLAB.
- **The joint sample is not reproducible.** `ParamSpaceAnalysis2` seeds it
  with `rng('shuffle')`. No regression tier covers it, and reruns will
  sample different connectivity points each time until that is decided.
- **`build()` reseeds the global random stream** and does not restore it.
  The spread's draw must use its own stream and restore, as the setpoint and
  ladder draws do.

## 9. Tom's answers (2026-11-11), and what is still open

**Settled:**

1. **Where the spread is allowed.** On all three kinds of route: chain,
   one-component, and multi-component product.
2. **The MTS chain condition's name.** `sfa3_stdc2`.
3. **Where the jitter is set.** In the **preset**, not in a condition, as
   the SFA spread is (`tau_a_spread` in the preset's struct). Condition
   names are unchanged: `sfa1_std1` and `sfa3_std2` keep their names with
   the spread on, and no `stdj` names are introduced. The reasons: the
   validator requires every condition of a preset to set the same fields;
   precedence is preset < condition < grid and the spread is a property of
   the network; the regression references use presets that never set it;
   and `same_config` refuses to pool runs whose `std_spread` differs.
5. **A fixed seed for the joint parameter-space sample.** Yes, as an
   option whose default is today's behaviour (FD-016). `ParamSpaceAnalysis2`
   gains a seed for its shuffle; left unset it calls `rng('shuffle')` as
   now, so existing bundles and the three regression references are
   untouched. New bundles set it, so the chain bundle and any later rerun
   sample the same 128 connectivity points and can be compared across
   runs. It is a setting, so the sample can be changed to check that a
   result does not depend on it. A small job for a Fable worker (it edits
   a class), done before the four-condition `fast` run of step 6; it can
   run beside step 3 because it touches only `ParamSpaceAnalysis2.m`,
   `run_param_space_analysis.m` and the new bundle's config.
6. **Reruns.** Agents never run `medium` or larger; only Tom starts them and
   says when they are done. The plan includes one `fast` run after block A is
   built. Analyses run on **four** conditions, and two figure sets are
   produced, three-condition and four-condition (section 4.4).
7. **Automatic pool creation.** Off, machine-wide; Tom has done it.
8. **Reserve refill.** 5 s.

**Still under discussion (the lead's recommendation in each):**

4. **Facilitation.** Recommendation: refuse it on chain routes in this pass.
   If a two-variable facilitation is wanted later, the form that matches the
   chain's logic is residual calcium driving a slower release-probability
   variable, with release probability also scaling how fast the releasable
   pool is used:

   - residual calcium: $\dot c = -c/\tau_{ca} + r$
   - facilitation: $\dot g = (1 - g)/\tau_{dec} + (G - g)\,c/\tau_{fac}$
   - releasable pool: the release term of (C1) becomes $g\,b_1 r/\tau_{rel}$
   - output: $\theta = r\, g\, b_1$

   Today's model drives $g$ by the rate directly and does not let $g$ change
   how fast resources are used, so a facilitated synapse transmits more
   without depleting faster. Not built, not planned for the draft.
**New, found 2026-11-11:**

9. **The default pool profile is `Just1`.** It has one worker, so pool code
   that uses the default profile breaks (section 8). Two ways to settle it:
   set the default profile back to `Processes` and have the licence-holding
   session open `Just1` by name; or leave `Just1` as the default and change
   the five call sites and the six worker-agent files to name `Processes`.
   The first changes nothing in the repo.

## 10. Done means

Code side of block A is done when:

- [ ] `test_std_chain_model` and `test_std_spread` pass and are in the quick
      suite (FD-005, FD-015).
- [ ] The three existing regression references pass unchanged on merged
      `main` (ultrarapid exactly; stages and full with their recorded ignore
      patterns).
- [ ] A chain preset builds, runs, sweeps, and goes through memory capacity
      with the synaptic readout (FD-004, FD-007).
- [ ] The four-condition bundle runs end to end at `fast`, and both figure
      sets build, three-condition and four-condition (FD-006, FD-007).
- [ ] The chain has its own regression reference (step 5).
- [ ] `Equations_stability_paper.md` states the chain, once.
- [ ] `ISSUES.md` and `WORKLOG.md` are current, and the to-do list's block A
      checkboxes are ticked for what is finished.

Then: the comparison figure (FD-014) and the manuscript text, with Tom. The
`medium` reruns (FD-009) are Tom's to start.
