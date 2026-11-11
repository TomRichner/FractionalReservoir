# Inventory: chained depression in display code, and the four-condition figure sets

Plan A, step 2 (serves FD-006 and FD-007). Read-only inventory taken on the
worktree at `c0300fc`; every `file:line` is as the line stands at that commit.
Nothing was changed to produce it.

**What it is for.** The multiple-timescale depression condition `sfa3_std2`
(a product of two depression variables) is joined by `sfa3_stdc2` (the
chained vesicle pools of Plan A, section 1). Every analysis runs four
conditions, in this order: `no_adaptation`, `sfa1_std1`, `sfa3_stdc2`,
`sfa3_std2`. Every figure is built twice from one run: a three-condition set
(`no_adaptation`, `sfa1_std1`, `sfa3_stdc2`, where the chain takes the colour,
title and position `sfa3_std2` has today) and a four-condition set (adding
`sfa3_std2` as a fourth column, in green). This note turns that into
line-level work for step 4.

**Terms used below.**

- *chain route*: `std = struct('model','chain','tau_rec',[2 5],'tau_rel',0.25,'n2',2)`.
  Two `tau_rec`, **one** `tau_rel`, plus `n2`. Its depression factor is
  `b(:,1)` (the releasable pool), not `prod(b,2)`. Its steady state is Plan A
  (C4): with `k1 = r*tau_rec(1)/tau_rel`, `k2 = r*tau_rec(2)/(n2*tau_rel)`,
  `s = 1 + k1 + k2`, `b1 = 2/(s + sqrt(s^2 - 4*k2))`, `b2 = 1 - k2*b1`.
- *the factor static*: `SRNNCellTypePairs.std_route_factor(b, params, pre, post)`
  proposed in Plan A 4.5 and added in step 1. Display code should take the
  factor from it rather than compute `prod(b)` itself.
- *the option*: `cfg.conditions` (a cell of names; empty means all the run
  has), passed by `make_all_paper_figures` to every figure as the name-value
  pair `'conditions'`, resolved by the new helper
  `src/figures/helpers/select_conditions.m` (the lead's decision).
- **A trap common to sections 1 and 2.** A chain route's `tau_rel` is a scalar
  and its `tau_rec` a 2-vector, so every product formula below
  (`1./(1 + r*tau_rec./tau_rel)` then `prod`) still **runs without error** on a
  chain route and returns the product of two independent single pools,
  `1/(1+8r) * 1/(1+20r)`, which is not the chain. None of these sites will fail
  loudly; each must branch on `route.model`. Sites that read the *condition
  struct* see the scalar `tau_rel`; sites that read *compiled params*
  (`params.tau_b_rel{pre,post}`) see whatever step 1 compiles (the product path
  today expands a scalar to `1 x nb`; `test_SRNNCellTypePairs.m:14`).

---

## 1. Reads of depression parameters outside the model class

Searched: `tau_rec`, `tau_rel`, `tau_b_rec`, `tau_b_rel`, `n_b_pairs`,
`std_geomean`, `.std` / `.std.` in `src/` (excluding
`src/model/SRNNCellTypePairs.m`), `scripts/paper/` (excluding `archive/`) and
`scripts/tests/`.

### 1a. Code that must change for a chain route

| # | Site | What it does | What a chain route needs |
|---|---|---|---|
| 1.1 | `src/figures/fig_STD_steady_state.m:96-100` | Reads `sc.(pre).(post).std` of the most-adapted condition; `tau_rec = route.tau_rec(:)'`, `tau_rel = route.tau_rel(:)'`, `rho = tau_rel./tau_rec`, `K = numel(rho)`. | Branch on `route.model`. Chain: `tau_rec(1)`, `tau_rec(2)`, scalar `tau_rel`, `n2`; there is no `rho` per component. Which condition is read comes from `full_adaptation_condition` at `:87` (see 3b). |
| 1.2 | `src/figures/fig_STD_steady_state.m:157-159` | Panel title prints `tau_rec`, `tau_rel`, `rho`. | Print `tau_rec1`, `tau_rec2`, `tau_rel`, `n2`. |
| 1.3 | `src/figures/fig_STD_strength_matching.m:216-226` (`std_route`) | Reads `route.std.tau_rec/tau_rel` of one named condition, `v.rho = tau_rel./tau_rec`. Table at `:174`, `:181` prints them. | Carry `model` and `n2`; `theta_ss` (2.4) must branch. Figure compares the strength-matching presets, which are product-only by purpose; it is not in the manuscript bundle's registry. |
| 1.4 | `src/figures/fig_main1_concepts.m:61-63` | `std = cond.synapse_config.E.E.std` of the most-adapted condition; `b_each = 1./(1+r_ref.*std.tau_rec./std.tau_rel)`. | Chain: `b1` from (C4) at `r_ref = 0.25` (0.200 with the Plan A values). Notes text `:106-109` prints `tau_rec`, `tau_rel` and the product: rewrite for the chain. |
| 1.5 | `src/figures/helpers/write_manuscript_tables.m:153-182` (`route_facts`) | Reads `model.synapse_config` (the raw property): `r.tau_rec = e.std.tau_rec(:)'`, `r.tau_rel = e.std.tau_rel(:)'`. | Also record `e.std.model` and `e.std.n2`. |
| 1.6 | `write_manuscript_tables.m:141-148` | `M.n_b_max = max(numel(tau_rec))`: the "M depression timescales" of the equations. | For a chain `numel(tau_rec) = 2` but the depression factor has one variable; `M` stops meaning "number of factors in the product". Needs a model flag next to it. |
| 1.7 | `write_manuscript_tables.m:226-229` | Emits the product `db_im/dt` equation whenever `n_b_max > 0`. | Chain: emit (C1)-(C3) instead, copied from `docs/EquationsParametersDocs/Equations_stability_paper.md` once step 1 adds them there (no new hand-written copy). |
| 1.8 | `write_manuscript_tables.m:275-306` | When `n_b_max > 1`, prints "the depression timescales enter as a **product** ... SQUARE the depression" and computes `ratios = tau_rec./tau_rel` (`:282`). | Wrong prose for a chain (sustained output rises to a ceiling, Plan A section 1). For a chain the ratios are `[8 20]`, unequal, so the numeric sentence is skipped but the product sentence at `:286-290` is still printed. Branch. |
| 1.9 | `write_manuscript_tables.m:390-398` | Parameter rows: `$M$` (`:392`), `$\tau_{rec}$` and `$\tau_{rel}$` per route (`:396-397`). | Add an `n2` row and a depression-model row for chain routes. |
| 1.10 | `write_manuscript_tables.m:474-502` (`cond_routes`) | Conditions table lists only *which routes* carry `std`. | `sfa3_stdc2` and `sfa3_std2` would print identical rows (same routes, same `tau_a`). Add a column naming the depression model (and the timescales). |

### 1b. Reads that need no change

| # | Site | Why no change |
|---|---|---|
| 1.11 | `src/figures/fig_adaptation_methods.m:236-240` (`combo_config`) | Copies the whole `route.std` struct into the column's config; a chain passes through intact. (The factor it then draws is 2.7.) |
| 1.12 | `src/analysis/run_transient_gain.m:228` | `B.std` is `lya_summary`'s leading-vector STD fraction; both pools count as depression (Plan A section 3: `lya_summary` unchanged). |
| 1.13 | `src/figures/helpers/write_manuscript_tables.m:136-138` | `SRNNModel2` branch (`tau_b_E_rec`); not used by the paper. |
| 1.14 | `src/model/SRNNModel2.m:2287, 2317, 3422, 3495-3512` | Comments in the class Plan A says to leave alone. |
| 1.15 | `src/presets/srnn_param_preset.m` (127 hits, e.g. `:706-713`) | All are preset definitions (writes), none reads. The new family is step 3 (`matlab-opus-2`). |
| 1.16 | Comments only: `src/analysis/run_eig_heatmap.m:289`, `src/analysis/run_tau_sensitivity_analysis.m:143`, `src/figures/fig_STD_steady_state.m:9-10, 62-64, 117-118` | `fig_STD_steady_state`'s comments become product-specific and should be updated with 1.1. |
| 1.17 | False positives: `src/analysis/run_paper_illustrations.m:85`, `src/analysis/run_stim_neurons.m:89` | `.std` there is an *output* field (the saved depression trace), not a parameter; it is a factor site (2.10, 2.12). |
| 1.18 | `scripts/paper/` (non-archive): one hit, `fig_to_do_in_future.md` | A note, not code. |

### 1c. Reserved for other workers

| # | Site | Owner |
|---|---|---|
| 1.19 | `src/model/SRNN_ESN_reservoir.m:547` (`has_std = any(params.n_b_pairs(:) > 0)`) | Step 3, `matlab-fable-2`. Chain has `n_b_pairs = 2`, so this stays true. |
| 1.20 | `src/model/SRNN_ESN_reservoir.m:996-1019` (`route_signature`: `n_b`, `tau_b_rec`, `tau_b_rel`, `std_geomean`) and the description that prints them, `:1020-1038` (`:1025-1026`) | Step 3, `matlab-fable-2`. **Without `model` in the signature a chain route and a product route with `tau_rec [2 5]` and a scalar `tau_rel 0.25` (which the product path broadcasts) would compare equal**, so the redundancy check could pass a mixed preset. |

### 1d. Tests

| # | Site | Status |
|---|---|---|
| 1.21 | `scripts/tests/test_SRNNCellTypePairs.m:12-17` (`n_b_pairs`, `params.tau_b_rel{1,1} == [0.2 0.2]`) | Product config; keep. |
| 1.22 | `scripts/tests/test_jacobian_times.m:48-49` (prints `n_b_pairs`) | Product config; keep. A chain case belongs to step 1's own tests. |
| 1.23 | `scripts/tests/test_psa_pairs_end_to_end.m:148-149` (`n_b_pairs == last digit of the name`) | Holds for `sfa3_stdc2` too (`n_b_pairs = 2`); runs on a product preset; keep. |
| 1.24 | `scripts/tests/test_mu7revised_bundle.m:24` (`isscalar(...std.tau_rec)` for 1TS STD) | Product bundle; keep. |
| 1.25 | `scripts/tests/test_pairs_single_celltype.m:152-153` | Product presets; keep. |
| 1.26 | `scripts/tests/test_route_scale.m:105-112` | Matching presets; keep. |
| 1.27 | `scripts/tests/test_preset_conditions.m:67-71, 104-107` (1TS takes the first `tau_rec`), `:378` (`name_vs_counts`), `:386-393` (`std_count`) | Step 3 (`matlab-opus-2`) adds the new family here. Note that `name_vs_counts` uses the patterns at `:353-354`, which **skip** `sfa3_stdc2`, so "names must not lie" would not cover the chain unless extended (see section 7). |
| 1.28 | `scripts/tests/test_std_product_model.m`, `scripts/tests/test_std_geomean.m` (53 hits between them) | Tests of the product model by purpose; keep. |
| 1.29 | `scripts/tests/test_esn_route_redundancy.m` (13 hits) | Step 3, `matlab-fable-2`. |
| 1.30 | Construction only (writes a product `std`, reads nothing): `regression/regression_rapid_run.m:150-154`, `test_jacobian_times.m:189-194`, `test_route_scale.m:87, 127-132`, `test_sample_hold_mc.m:49-52`, `test_SRNNCellTypePairs.m:142-147, 195-216`, `test_SRNNCellTypePairs_S_c_heterogeneity.m:134-136`, `test_SRNNCellTypePairs_tau_a_spread.m:234, 307-312`, `test_transient_gain.m:186-191`, `test_psa_model_class.m:37`, `test_pairs_single_celltype.m:197`, `test_n_a_dependent.m:9` (comment) | No change. |

---

## 2. Sites that compute the depression factor themselves

For a chain route the factor is `b(:,1)`; through the factor static where a
`b` array and `params` are in hand, through (C4) where the figure is analytic.

### 2a. Must change

| # | Site | What is drawn or computed | What must change |
|---|---|---|---|
| 2.1 | `src/figures/fig_STD_steady_state.m:103-113` | Steady-state `b_each = 1./(1 + r/rho)`, `b_prod = prod(b_each)`, `theta = b_prod.*r`; the peak search at `:109-113`. | Chain: `b1(r)` from (C4); `theta = r*b1` rises monotonically to the ceiling `1/(tau_rec1/tau_rel + tau_rec2/(n2*tau_rel))` (1/18 with the Plan A values), so `has_peak` is false and the "R =" annotation (`:182-191`) does not apply. Also draw `b2`. The dashed "single `b_k`" reference (`b_each(1,:)`) is exactly the 1TS pool, which is a meaningful comparison for the chain. |
| 2.2 | `fig_STD_steady_state.m:121-140` | Step response integrated in closed form per segment (`:127-129`: linear ODE per component), `b_prod_t = prod(b_t)` (`:135`), `ss_b`, `ss_theta` (`:138-139`), `recovery_frac` from `max(tau_rec)` (`:140`). | The chain is bilinear (`b2*(1-b1)`), so there is no closed form: integrate (C1)-(C2) numerically at `cfg.fs` (an explicit step at 1 kHz is ample; time constants are 0.25 s and up). Plateaus from (C4). Recovery is set by `tau_rec2` and the reserve, not `max(tau_rec)` alone. |
| 2.3 | `fig_STD_steady_state.m:168-171, 196-201, 218-222, 237-243` | Axis and legend text `\prod_k b_k`. | `b_1` (releasable pool) for a chain. |
| 2.4 | `src/figures/fig_STD_strength_matching.m:228-231` (`theta_ss`) and the label at `:119` | `scale .* r .* prod(1./(1 + r/rho))`. | Chain: `scale .* r .* b1(r)` from (C4). |
| 2.5 | `src/figures/fig_main1_concepts.m:63-67` | `b_min = prod(b_each)`; `b_levels = linspace(1, b_min, 5)` drive the STD sigmoids and discs; notes at `:107, :111-113` say "frozen product B". | Chain: `b_min = b1(r_ref)` from (C4). Rewrite the notes (the "product must not be identified with a single depression variable" sentence at `:113` is exactly wrong for the chain, whose factor IS one variable). |
| 2.6 | `src/figures/fig_example_timeseries.m:179-197` (`Pb = prod(B, 2)` at `:186`) | Row (e), depression of the shown neurons, route 1. Row label `:94` (`'depression  \Pi b'`), header `:22`. | Factor static (route 1 of each presynaptic type, as now). Label: `b_1` on a chain column, `\Pi b` on a product column. |
| 2.7 | `src/figures/fig_adaptation_methods.m:262-267` | Row `b`: `squeeze(prod(p.b.E.E(1,:,:), 2))`. | Factor static. In this bundle the entry runs `single_neuron_stf` (product, one component), so nothing changes numerically today; it matters only if a chain single-neuron preset is ever drawn. |
| 2.8 | `src/figures/helpers/write_manuscript_tables.m:203-205` | `theta` factor `\prod_{m=1}^{M} b_{im}` in the equation block. | Chain: `b_{i1}`. |
| 2.9 | `write_manuscript_tables.m:291-305` | `1/(1 + rr)^M` depth numbers. | Product-only; skip for a chain (see 1.8). |
| 2.10 | `src/analysis/run_paper_illustrations.m:81-85` (`sb = prod(B(ix,:,:), 2)` at `:84`) | Saves `r.std` for `representative_dynamics` and `single_neuron`. Consumers: `fig_main2_dynamics.m:47-50` (row "Synaptic Depression"), `fig_single_neuron_revised.m:16-17` (label `'STD product'`). | Factor static. Analysis stage: the saved `std` changes meaning per condition, so the label in `fig_single_neuron_revised.m:17` must stop saying "product". |
| 2.11 | `src/analysis/run_stim_evoked.m:194-197` (`dep = prod(B, 2)` at `:196`); header `:26`, `:28` | Population-mean `std` readout. | Factor static. (Not in the manuscript bundle; it is the stimulation-experiment stage.) |
| 2.12 | `src/analysis/run_stim_neurons.m:88-89`; header `:21` | Per-neuron `std` trace. | Factor static. (Not in the bundle.) |

### 2b. No change (and why)

| # | Site | Why |
|---|---|---|
| 2.13 | Consumers of the model's own `synaptic_output`: `src/analysis/ParamSpaceAnalysis2.m:2468-2476`, `fig_example_timeseries.m:151-154`, `fig_adaptation_methods.m:253-256`, `run_paper_illustrations.m:79-80`, `run_stim_evoked.m:197`, `run_stim_neurons.m:84-86` | The model computes it (`unpack_and_compute_states`); correct once step 1 lands. |
| 2.14 | `src/plotting/plot_SRNN_combined_tseries.m:152-170` (`prod(b.I, 2)`, `prod(b.E, 2)` at `:160`, `:168`) | Reads `p.n_b_E`, `plot_data.b.E`: `SRNNModel2`-shaped. Its only caller is `scripts/examples/Single_vs_dual_adaptation_example.m:253`. Unreachable for a chain; leave. |
| 2.15 | `src/plotting/plot_SRNN_paired_pulse_tseries.m:30` | `SRNNModel2`-shaped; leave. |
| 2.16 | `src/model/jacobian/compute_J_eff.m:139, 142` | `SRNNModel2`; Plan A says leave alone. |
| 2.17 | Text only: `src/figures/fig_sine_paired_stim.m:21, 24, 176, 291`, `src/figures/fig_sine_stim_poster.m:20`, `src/analysis/stim_evoked_measures.m:196-197` | They say "depression product" about the `std` readout of 2.11. Update the words with 2.11 if those figures are regenerated on a chain preset. |
| 2.18 | Tests of the product by purpose: `test_pairs_single_celltype.m:81-85`, `test_SRNNCellTypePairs.m:44-46`, `test_std_geomean.m:50, 61`, `test_std_product_model.m:91, 116-118, 139, 291`, `test_route_scale.m:111` | Keep. |

### 2c. Reserved

| # | Site | Owner |
|---|---|---|
| 2.19 | `src/model/SRNN_ESN_reservoir.m:662-689` (plot 6 of `plot_esn_timeseries`, which starts at `:488`; `pr = reshape(prod(bq, 2), ...)` at `:674`) | Step 3, `matlab-fable-2`. |
| 2.20 | `src/model/SRNNCellTypePairs.m:1791` (`plot`: `collapse_routes(p.b, 'prod')`), `:1849` (`plot_celltypes`: same), `:3220-3237` (`collapse_routes`, `prod(value, 2)` at `:3229`) | Model class; the lead assigns (Plan A step 4 suggests Opus or a Fable follow-up). Feeds `fig_stim_engages_adaptation` (`model.plot()` at `:91`) but that figure runs `bursting_pairs` (product). |

---

## 3. Sites keyed on or naming `sfa3_std2`

Classes: **(a)** a name-to-value map that must gain `sfa3_stdc2`; **(b)** a
hard-coded condition list or logic; **(c)** a test expectation; **(d)** a
comment only.

### (a) Maps that must gain `sfa3_stdc2`

| # | Site | Now | Needs |
|---|---|---|---|
| 3.1 | `src/presets/srnn_condition_titles.m:65` | `'sfa3_std2' -> 'Multiple-Timescale Adaptation'` | `sfa3_stdc2 -> 'Multiple-Timescale Adaptation'`; `sfa3_std2` retitled (for example "Multiple-Timescale Adaptation (product)"). **Must land before the new preset can validate**: `validate_preset_conditions` requires every condition name to be a key here. Step 3. |
| 3.2 | `src/figures/helpers/manuscript_style.m:114` (`condition_color`) | `sfa3_std2 -> std_hue` (Okabe-Ito blue) | `sfa3_stdc2 -> std_hue`; `sfa3_std2 -> green`. Many figures index this map **unguarded** (the header at `:76-79` says so), so a missing key throws `MATLAB:Containers:Map:NoKey`. Step 3. See section 7 on green. |
| 3.3 | `manuscript_style.m:143` (`condition_short`) | `sfa3_std2 -> 'SFA+STD'` | A short label for each MTS (used as tick labels by `fig_transient_amplification.m:89-91`, legends by `fig_memory_capacity_example.m:70-75`, `fig_STD_strength_matching.m:151, 157`). Step 3. |
| 3.4 | `src/figures/fig_sine_paired_stim.m:364-365` (local `cond_label` map) | three names -> three titles | Falls back to `srnn_condition_titles` for other names (`:366-370`), so it degrades rather than errors; add the key if the figure is regenerated on the chain preset. |

`mc_condition_colors.m` and `mc_display_names.m` read the two maps above and
need **no** change (their `sfa3_std2` mentions, `mc_condition_colors.m:8, 10`,
are comments).

### (b) Hard-coded lists and logic

| # | Site | Now | Needs |
|---|---|---|---|
| 3.5 | `src/presets/full_adaptation_condition.m:75, 80` (`parse_condition_name`) | Two exact patterns `^sfa(\d+)_std(\d+)_stf(\d+)$`, `^sfa(\d+)_std(\d+)$`; `sfa3_stdc2` matches neither and is ignored. With four conditions it returns `sfa3_std2` (the product). | Accept `stdc` (keep two exact patterns per the header at `:68-73`, or add a third). Then `sfa3_stdc2` and `sfa3_std2` tie on X+Y = 5 and on X = 3, and the strict `>` at `:50` returns whichever is listed first: add an explicit rule preferring the chain. Step 3. **Callers**, all of which then follow the chain: `src/figures/helpers/build_from_preset.m:57` (empty condition name; used by `fig_SFA_steady_state.m:169`, `fig_stim_engages_adaptation.m:80`, `write_manuscript_tables.m:74`), `src/figures/fig_main1_concepts.m:59`, `src/figures/fig_adaptation_methods.m:95`, `src/figures/fig_STD_steady_state.m:87`, `src/figures/fig_memory_capacity_example.m:80-82`, `src/figures/helpers/write_manuscript_tables.m:78`, `src/analysis/run_tau_sensitivity_analysis.m:84`, `src/analysis/check_sensitivity_sim.m:99`, `src/util/write_run_parameters_md.m:351`. |
| 3.6 | `src/figures/fig_main2_dynamics.m:52` | `names = {'no_adaptation','sfa1_std1','sfa3_std2'}` colours the titles positionally | Names from the saved results (`D.results(c).name`), filtered by the option. |
| 3.7 | `src/figures/helpers/style_main4_grouped.m:24` and `:33` | Same three-name list, twice (title colouring; legend labels) | Derive from the legend's own entries or take the selected names. |
| 3.8 | `src/figures/helpers/style_main5_grouped.m:25` | Same list; with `labels` at `:26` | As 3.7. |
| 3.9 | `src/figures/fig_STD_strength_matching.m:61-63, 69, 76, 82, 107, 126` | Reads the `sfa3_std2` route of three matching presets; colours from `sfa3_std2`; occupied rates for three names | Product-only by purpose (the matching study). Not in the manuscript bundle. If kept, take the MTS name from the preset. |
| 3.10 | `src/figures/fig_sine_paired_stim.m:63`, `src/figures/fig_sine_stim_poster.m:48` (`cfg.show_condition = 'sfa3_std2'`), `fig_sine_stim_poster.m:213` (series list) | Defaults | Stimulation figures, not in the bundle. Default to `full_adaptation_condition` when regenerated on a chain preset. |
| 3.11 | `src/analysis/run_stim_neurons.m:32` (`opts.condition = 'sfa3_std2'`; usage text `:5`) | Default | As 3.10. |
| 3.12 | `src/presets/srnn_param_preset.m` (23 occurrences) | Product condition definitions in existing presets | Untouched; the new family is step 3. |

### (c) Test expectations

| # | Site | Effect of the change |
|---|---|---|
| 3.13 | `scripts/tests/test_preset_conditions.m:128-129` (`t('sfa3_std2') == 'Multiple-Timescale Adaptation'`) | **Fails** when `sfa3_std2` is retitled. Step 3 must update it (and `:130-131`), and add `sfa3_stdc2`. |
| 3.14 | `test_preset_conditions.m:59, 81, 94, 115-117, 217` | Order and identity checks on the existing 3/4/7-condition product presets; unaffected. New-family checks go beside them (step 3). |
| 3.15 | Product-preset tests naming `sfa3_std2` as a model to run: `test_c_over_K.m:39, 63`; `test_jacobian_times.m:126-127`; `test_lyapunov_spectrum_stage.m:70`; `test_lyapunov_topk.m:11, 59, 163, 165`; `test_lya_K_auto.m:6, 38`; `test_mu7revised_bundle.m:13`; `test_n_a_dependent.m:82`; `test_pairs_reproducibility.m:15, 18, 27, 92-105`; `test_pairs_single_celltype.m:173`; `test_psa_pairs_end_to_end.m:47, 111, 127-128, 135`; `test_route_scale.m:99`; `test_SRNNCellTypePairs_tau_a_spread.m:10, 32, 70, 268, 270`; `test_sweep_metrics.m:67, 84`; `test_tau_a_EI_alias.m:42`; `test_transient_gain.m:20, 122` | Unaffected: they run existing product presets, which keep their names. |

### (d) Comments only

`src/presets/full_adaptation_condition.m:13, 18`; `src/presets/srnn_condition_titles.m:31, 42, 60`;
`src/util/write_run_parameters_md.m:221`; `src/figures/fig_adaptation_methods.m:93`;
`src/figures/fig_memory_capacity_example.m:79`; `src/figures/fig_SFA_steady_state.m:168`;
`src/figures/fig_sine_paired_stim.m:12`; `src/figures/fig_sine_stim_poster.m:11`;
`src/figures/fig_STD_steady_state.m:83`; `src/figures/fig_stim_engages_adaptation.m:78`;
`src/analysis/run_tau_sensitivity_analysis.m:22, 77`; `src/figures/helpers/build_from_preset.m:7`;
`src/figures/helpers/manuscript_style.m:93, 100`; `src/figures/helpers/mc_condition_colors.m:8, 10`;
`scripts/tests/regression/regression_digest.m:23`;
`scripts/paper/mu7revisedMedium_config.m:21`, `regression_fast_config.m:48`,
`stdScaled_fast_config.m:45`, `stdSingleMatched_fast_config.m:45`,
`stdUsage_fast_config.m:45`, `stdUsage_med_config.m:45`. (22 sites.)

---

## 4. Sites that assume exactly three conditions (or exactly one)

| # | Site | The assumption | Needs |
|---|---|---|---|
| 4.1 | `src/figures/fig_main2_dynamics.m:18` | `assert(numel(D.results)==3, ...)` | Select by name from the saved results; assert the selection, not the stage's count. |
| 4.2 | `fig_main2_dynamics.m:51` | Positional `titles` for three columns | `srnn_condition_titles` by name (and note the stage froze `r.title`, section 7). |
| 4.3 | `fig_main2_dynamics.m:53, 55` | `for c=1:3`; column x positions `.105+(c-1)*.30`, width `.245` | Column geometry from N. |
| 4.4 | `fig_main2_dynamics.m:42` | Figure `1450 x 1200` sized for three | Width from N. |
| 4.5 | `fig_main2_dynamics.m:93-95` | E/I legend placed on `ax_all(2,3)` (the last column) | `ax_all(2,end)`. |
| 4.6 | `fig_main2_dynamics.m:101` | Column dividers at `[.3775 .6775]` | From N. |
| 4.7 | `src/figures/fig_main3_stability.m:10, 13, 73` | Positional `titles`; `for c=1:3` twice | By name; loop over the selection. |
| 4.8 | `fig_main3_stability.m:14-15` | `D.results(c)` (local_lyapunov stage) and `old(c)` (the c-th tile of `fig_eig_heatmap`) both positional | Select both by name: the stage carries four conditions, and the heat-map tiles come from a different stage (`eig_heatmap`) whose order happens to match. |
| 4.9 | `fig_main3_stability.m:15, 35, 80` | Axis x positions `.105+(c-1)*.30`, `.705` for column 3 | From N. |
| 4.10 | `fig_main3_stability.m:32-36` | Colour bar attached to `c==3` | Attach to the last column. |
| 4.11 | `fig_main3_stability.m:54` | `if c==1, ylim([-.5 4]) else ylim([-.5 .5])`: assumes column 1 is `no_adaptation` | Key on `r.name`. |
| 4.12 | `fig_main3_stability.m:57` | Label offset `if c<=2` | Key on name or value. |
| 4.13 | `fig_main3_stability.m:6, 71` | Figure size; dividers `[.3775 .6775]` | From N. |
| 4.14 | `src/figures/fig_grouped_main.m:54` | Overwrites the margin panel's `XTickLabel` with three strings | With four ticks MATLAB **cycles** the three labels silently (section 7). Use the selected names. |
| 4.15 | `fig_grouped_main.m:57, 63` | `tiledlayout(pan,1,3)` for the transient-gain row, then loops over `numel(D.results)` | Four results overflow a 1 x 3 layout (error). Select by name; layout 1 x N. |
| 4.16 | `fig_grouped_main.m:93` | Overwrites the memory-capacity horizon panel's `XTickLabel` with `{'No adaptation','1TS','MTS'}` | As 4.14. |
| 4.17 | `fig_grouped_main.m:122, 125` | Tag regex `'dynamics_r1_c[123]'`, `'stability_r1_c[123]'` | `c\d+`. |
| 4.18 | `fig_grouped_main.m:25` | Canvas `1300 x 950` | Group 2 and 3 copy their source's size; groups 4-6 are resized by the styles. Check after 4.19-4.27. |
| 4.19 | `src/figures/helpers/style_main4_grouped.m:10` | `assert(numel(axs)==9)`: six sensitivity panels plus **three** rate panels from `fig_lle_vs_rate` | `6 + N`. |
| 4.20 | `style_main4_grouped.m:24, 33-37` | Three names, `colored = cell(1,3)`, `for j=1:3` | From the selection. |
| 4.21 | `style_main4_grouped.m:46-50` | `legend(..., [handles; zero_ref], [colored {'Zero Growth'}])`: handle count N + 1 against 3 + 1 labels | Errors at N = 4. |
| 4.22 | `style_main4_grouped.m:6, 54-65` | Canvas `910 x 807.5`; legend space "for all four keys" (3 conditions + zero) | Size and legend space from N. |
| 4.23 | `src/figures/helpers/style_main5_grouped.m:6` | `assert(numel(axs)==4)`: one margin plus **three** gain axes | `1 + N`. |
| 4.24 | `style_main5_grouped.m:15-19` | Four panel slots at `.055+(k-1)*.242` | `1 + N` slots. |
| 4.25 | `style_main5_grouped.m:25-37` | Three names and labels; `gains(k).Title.Color` positional; `margin.XTick=1:3`; `gains(1..3)` titles written literally | From the selection. |
| 4.26 | `style_main5_grouped.m:41-46` | Legend centred under `ordered(3)`, the middle of four panels | Centre from N. |
| 4.27 | `src/figures/helpers/style_main6_grouped.m:6, 8, 20` | `assert(numel(axs)==5)` and exactly **two** "slowest" (tau) axes, both titled `'Multiple-Timescale Adaptation'` | Assumes the tau sweep carries **one** condition: `fig_sfa_EOC_allStd`'s block figure has one subplot per tau condition (4.30), so two tau conditions give three tau axes. Titles by name. |
| 4.28 | `style_main6_grouped.m:52` | `bottom(3).XTickLabel = {'No Adaptation','1TS','MTS'}` | As 4.14. (`:47-48`, `for k=1:3`, is the three MC panels, not conditions.) |
| 4.29 | `src/figures/fig_sfa_EOC_allStd.m:109-115, 120, 161-164, 238` | **One condition**: `cond_use = cond_names{1}` with a warning if more; the density image, band, median and table are for that one | Both MTS conditions: overlay medians and bands (one colour each) on one density panel, or one panel per condition. The block figure at `:196-221` already lays out one subplot per condition. |
| 4.30 | `src/analysis/run_tau_sensitivity_analysis.m:7-8, 12-31` (header), `:83-87` (picks the one `full_name`), `:124` (`set_conditions(condition)`) | Sweeps **one** condition, the most-adapted | Sweep both MTS conditions in one PSA (same network seed per grid point, so paired). Safe, because both carry the same three-rung SFA ladder and the grid axis overrides `tau_a` identically for each (the header's limitation 2 at `:26-31`); assert the ladders are equal. |
| 4.31 | `run_tau_sensitivity_analysis.m:98-105` | Axis choice and `n_elements` from `condition{1}` | Assert identical across the swept conditions. |
| 4.32 | `run_tau_sensitivity_analysis.m:207, 223` (`write_tau_levels`) | `cond = psa.conditions{1}`; `tau_levels.md` names one condition | Loop, or assert equal and name all. |
| 4.33 | `src/util/write_run_parameters_md.m:726-735` | Prints the **first** condition's `synapse_config` and says it is "identical for every condition that has depressing routes" | Already false (1TS has one `tau_rec`, MTS two); false again for chain vs product. Print per condition. |
| 4.34 | `src/figures/fig_single_neuron_revised.m:15, 18-19` | `tiledlayout(6,3)`, `for c=1:3`, `D.results(c)` positional | The `single_neuron` stage runs `cfg.single_neuron_preset`'s conditions (no adaptation / 1TS SFA / 1TS STD), so three is that preset's count. Make it N-general; whether the chain bundle's single-neuron preset gains a chain column is a preset decision (section 7). |
| 4.35 | `src/figures/fig_eig_heatmap.m:119-122` | `floor(sqrt(N))` rows: 3 -> 1 x 3, **4 -> 2 x 2** | Not a fourth column. Decide whether this figure stays square-ish; `fig_main3_stability` copies tiles positionally either way (4.8). |
| 4.36 | `src/figures/fig_STD_strength_matching.m:64-69, 76, 110-113` | Four variants and three conditions, literal | Product matching study; not in the bundle (3.9). |
| 4.37 | `src/figures/fig_sine_paired_stim.m:364-365`, `src/figures/fig_sine_stim_poster.m:213-214` | Three-name label map; two-series list | Stimulation figures, not in the bundle. |
| 4.38 | Tests: `test_lyapunov_spectrum_stage.m:41`, `test_transient_gain_stage.m:45`, `test_numerics_verification_stage.m:51, 55, 70`, `test_preset_conditions.m:94`, `test_psa_pairs_end_to_end.m:47` | `numel(R) == 3` / literal three-name lists | Correct for the three-condition presets they run. A four-condition stage test would be new. |
| 4.39 | Comment: `fig_grouped_main.m:102` ("three-condition learning experiment") | | Group 7 has no condition data. |

False positives checked and excluded: `fig_main1_concepts.m:18-20, 31, 34, 96`
(three Sompolinsky gains and three row labels), `fig_sine_stim_neurons.m:113`
(three neuron classes), `fig_sine_stim_poster.m:157` (three rows),
`style_main6_grouped.m:47-48` (three MC panels).

Already N-general, for reference: `fig_lle_vs_rate`, `fig_local_vs_finite_lle`,
`fig_sensitivity_medians`, `fig_eig_heatmap_imbalance`, `fig_lyapunov_spectrum`,
`fig_transient_amplification`, `fig_transient_gain`,
`fig_transient_gain_excursions`, `fig_local_lyapunov`,
`fig_numerics_verification`, `fig_memory_capacity`,
`fig_memory_capacity_example`, `fig_dc_lle`, `fig_example_timeseries`, and the
PSA plotters (`ParamSpaceAnalysis2.plot` `:958-970`, `plot_sensitivity`
`:1087-1132`, `load_and_make_unit_histograms.m:84-212`).

---

## 5. The manuscript bundle's registry, entry by entry

Registry: `scripts/paper/mu7revisedMedium_config.m:92-160` (28 entries) plus
`scripts/paper/grouped_figure_registry.m:9-15` (8 entries).

### 5.0 What applies to every entry

- **Every registered function must accept `'conditions'`.**
  `make_all_paper_figures.m:149-150` passes one fixed list of name-value pairs
  to every entry, and every figure's `arguments` block is strict.
  `call_figure` (`:257-276`) retries only when the error mentions `run_dir`,
  so an undeclared `conditions` fails the entry. Either every function declares
  `cfg.conditions` (the same uniform contract as `run_dir`; the analytic ones
  accept and ignore it), or the master passes it only to entries that declare
  it. The rows below assume the first.
- **`make_all_paper_figures`** needs: `cfg.conditions` defaulted to `{}` when
  absent (next to the `verbose` default at `:62`), added to the args at
  `:149-150`, and recorded by `write_figure_manifest.m:51-52` (a "Conditions"
  row), so a `_3cond` and a `_4cond` folder say which they are.
- **`select_conditions`** is applied wherever a figure forms its condition
  list. Suggested contract: `[names, idx] = select_conditions(available, requested)`;
  empty `requested` returns `available` unchanged; output in **requested**
  order (that is what puts the chain in product MTS's slot); a requested name
  the run does not have is an **error** (no silent fallback, as
  `resolve_data_file`). The sweep-sheet figures get their columns from
  `ParamSpaceAnalysis2` plotters that iterate `psa.conditions`; that property is
  public, so the figure can set `psa.conditions = psa.conditions(idx)` on the
  loaded object (never saved) and `ParamSpaceAnalysis2` itself needs no edit.
  Results are keyed by name, so a subset is safe.
- **Titles saved by the stages.** Several stages save a display title at
  compute time and the figures draw it (section 7, R3). Where a row says "by
  name", it means titles and colours from `srnn_condition_titles` /
  `manuscript_style` keyed on the saved *name*.
- **Two figure roots.** `grouped_figure_registry(cfg.preset_name, cfg.fig_root)`
  is evaluated once in the config (`mu7revisedMedium_config.m:160`), so its
  `source_fig_root` is baked in. Group 2 (its archive fallback) and group 8
  (`fig_main8_psd.m:21`) read `<source_fig_root>/<entry>/...` written earlier in
  the same pass. Two passes into `figs/<bundle>_3cond` and `figs/<bundle>_4cond`
  need the registry rebuilt per pass.

### 5.1 The 28 individual entries

"Source" is where the figure gets its condition list today. "N-general" is
whether the layout already follows the count. "Read the option at" is where
`select_conditions` would go.

| Entry | Function | Source | N-general | 3 named / 4 | Read the option at |
|---|---|---|---|---|---|
| `fig_introductory_concepts` | `fig_introductory_concepts` | Literal: `sompolinsky_pairs`, `'no_adaptation'` (`:66`) | n/a (three gains) | Nothing; accept and ignore. | arguments block only |
| `fig_FI_curve` | `fig_FI_curve` | None (analytic) | n/a | Accept and ignore. | arguments block only |
| `fig_single_neuron_revised` | `fig_single_neuron_revised` | Run data, `single_neuron` stage, positional (`:18-19`) | **No** (4.34) | Name-keyed loop over the selection; label `:17` must stop saying "STD product" (2.10). Which columns exist depends on the single-neuron preset (section 7). | `:18` (over `{D.results.name}`) |
| `fig_adaptation_methods_stf` | `fig_adaptation_methods` | Preset `single_neuron_stf`, `full_adaptation_condition` (`:95`) | Yes (columns are mechanism combos) | No change for this bundle (product preset). Factor static at `:266` (2.7) for safety. Accept and ignore. | arguments block only |
| `fig_SFA_steady_state` | `fig_SFA_steady_state` | Preset via `build_from_preset(preset)` (`:169`) | n/a | None beyond 3.5 (both MTS share `tau_a`). Accept and ignore. | arguments block only |
| `fig_STD_steady_state` | `fig_STD_steady_state` | Preset, most-adapted condition (`:85-87`) | n/a (one route) | **Chain physics** (1.1, 1.2, 2.1-2.3). For the 4-condition set, draw the product as a second curve set (green), or a second figure; for the 3-condition set the chain only. Use the option to decide which MTS routes to draw. | `:85-87` |
| `fig_stim_engages_adaptation` | `fig_stim_engages_adaptation` | Preset `bursting_pairs`, most-adapted (`:80`) | n/a | None (product preset; `model.plot()` is reserved, 2.20). Accept and ignore. | arguments block only |
| `fig_sensitivity_analysis_allStd` | `fig_sensitivity_analysis_allStd` | Run data: `replot_sensitivity` -> `ParamSpaceAnalysis2.plot_sensitivity` iterates `psa.conditions` | Yes (columns found from the axes, `:327-341`) | Subset `psa.conditions` inside `replot_sensitivity` (`:87`, after `from_dir`), so `replot_sensitivity` gains an argument. Titles come from `srnn_condition_titles` (3.1). Figure width scales with N in `plot_sensitivity` (`:1128`). | `:159` (pass to `replot_sensitivity`) |
| `fig_sensitivity_medians` | `fig_sensitivity_medians` | **Preset** conditions (`:114-119`), unguarded map lookups | Yes (overlaid curves, legend `:344-346`) | Filter `cond_names` at `:115`; needs 3.1-3.2 (NoKey otherwise). | `:115` |
| `fig_param_space_allStd` | `fig_param_space_allStd` | Run data: `replot_param_space_analysis` -> `ParamSpaceAnalysis2.plot` iterates `psa.conditions` | Yes (`nCols = numel(src_ax)`, `:70`) | Subset in `replot_param_space_analysis.m:42` (after `from_dir`), which gains an argument. | `:46` (pass to the replot) |
| `fig_EI_weights_param_space` | `fig_EI_weights_param_space` | Run data via `load_and_make_unit_histograms` (iterates `psa.conditions`, `:84`); the title map from the preset (`:124-128`) | Yes (`ei_metric_sheet` lays out `nGrid` columns) | `load_and_make_unit_histograms` gains a `Conditions` option (it is in `src/plotting/`, shared code: flag to the lead); filter `cnames` at `:125` likewise. | `:102` (pass) and `:125` |
| `fig_lle_vs_rate` | `fig_lle_vs_rate` | Run data, all PSA conditions (`:88-109`) | Yes | Filter `cn` at `:88` (or `cond_names` at `:109`, which is then in **fieldnames order** = load order; use the selection's order instead). | `:88`, `:109` |
| `fig_local_vs_finite_lle` | `fig_local_vs_finite_lle` | Run data (`load_near_default`, `load_joint`) | Yes | Filter `cn` at `:223`, `:248` (or `near.cond_names` at `:71` and reorder). | `:71` |
| `fig_sfa_EOC_allStd` | `fig_sfa_EOC_allStd` | Run data, tau PSA, **first condition only** (`:110`) | Main panel **no** (4.29); block figure yes | Draw every selected tau condition that the run has (two MTS in the 4-set, the chain alone in the 3-set); colour by name; the table names each. Depends on 4.30. | `:109-110` |
| `fig_memory_capacity` | `fig_memory_capacity` | Run data `R.conditions` (`:71`) | Yes | Index `cond_keys`, `H_trials`, `MC_trials`, `R2_*`, `sm.*` rows by the selection; keep only `stats` pairs whose two names are both selected (pairs from `run_memory_capacity.m:241-244`, labelled `"a vs b"`). | `:71` |
| `fig_memory_capacity_example` | `fig_memory_capacity_example` | Run data `S.condition_names` (`:54`) | Yes | Index `results`, `R2` by the selection; `recon_name` (`:80-82`) must be the chain once 3.5 lands; legend uses `condition_short` (3.3). | `:54` |
| `fig_eig_heatmap` | `fig_eig_heatmap` | Run data `D.condition_titles` (`:53`), positional | Yes (but 4 -> 2 x 2, 4.35) | Select by `D.cond_names` (saved by `run_eig_heatmap.m:207`); retitle by name. | `:52-56` |
| `fig_eig_heatmap_imbalance` | `fig_eig_heatmap_imbalance` | Run data `D.cond_names`, `D.condition_titles` (`:61-64`) | Yes (rows x N columns, `:107-108`) | Select columns by name in `ex(e).evals_by_cond`, `lle_by_cond` etc. | `:61-64` |
| `fig_dc_lle` | `fig_dc_lle` | Run data `r.conditions` (`:56`) | Yes (guarded colours, `:69-75`) | Filter `conditions` at `:56` and every per-condition array with it. | `:56` |
| `fig_lyapunov_spectrum` | `fig_lyapunov_spectrum` | Run data `D.results` (`:42`) | Yes | `R = R(idx)` by `{R.name}`. | `:42` |
| `fig_transient_amplification` | `fig_transient_amplification` | Run data `D.cond_names` (`:55-57`) | Yes | Select by name; ticks use `condition_short` (3.3). Group 5 then overwrites them (4.14). | `:55` |
| `fig_transient_gain` | `fig_transient_gain` | Run data `D.results` (`:48`) | Yes | `R = R(idx)`. | `:48` |
| `fig_transient_gain_excursions` | `fig_transient_gain_excursions` | Run data `D.results` (`:43`) | Yes | `R = R(idx)`. | `:43` |
| `fig_local_lyapunov` | `fig_local_lyapunov` | Run data `D.results` (`:37`) | Yes | `R = R(idx)`; colour by name (`:47`). | `:37` |
| `fig_numerics_solver` | `fig_numerics_verification` (`'variant','solver'`) | Run data `R` (`:67`) | Yes | `R = R(idx)` before `n_cond`; the `verdict` table also by name. | `:67` |
| `fig_numerics_lya_method` | `fig_numerics_verification` (`'variant','lya_method'`) | as above | Yes | as above | `:67` |
| `fig_numerics_jacobian` | `fig_numerics_verification` (`'variant','jacobian'`) | as above (`D.results`, `:89`) | Yes | as above | `:67`, `:89` |
| `doc_tables` | `fig_doc_tables` -> `write_manuscript_tables` | Preset (`write_manuscript_tables.m:70-78`) | n/a | **Chain equations and parameters** (1.5-1.10, 2.8-2.9). The parameter table documents ONE condition (the most-adapted, so the chain after 3.5); with both MTS in the run it should state both routes' parameters. The conditions table should list the selected conditions and name each one's depression model. | `fig_doc_tables.m:26-40` (pass through) and `write_manuscript_tables.m:70`, `:82`, `:85` |

### 5.2 The eight grouped entries

All eight are `fig_mainN_grouped(varargin) -> fig_grouped_main(N, varargin{:})`
(`fig_main1_grouped.m:3` and siblings), with `preset_name`,
`source_fig_root`, `pytorch_file`, `human_psd_file` from the registry. **The
wrappers need no edit**; `fig_grouped_main`'s `arguments` block (`:6-18`)
gains `cfg.conditions` and forwards it to every inner call.

| Entry | Path through the helpers | Source | N-general | 3 named / 4 |
|---|---|---|---|---|
| `fig_main1_grouped` | `fig_grouped_main:28-34` -> `fig_main1_concepts(preset, intro_root)` | Archived Sompolinsky FIGs (no conditions) plus the preset's most-adapted STD route (`fig_main1_concepts.m:58-61`) | n/a | Chain STD concept panel (1.4, 2.5). The 4-condition version: decide whether it shows the product as well (section 7). Pass the option into `fig_main1_concepts`, whose positional signature `(preset_name, intro_root)` gains a third argument. |
| `fig_main2_grouped` | `fig_grouped_main:35-38` -> `fig_main2_dynamics(run_dir, source_fig_root)` (numeric path, `representative_dynamics` stage from `run_paper_illustrations`) | Run data, positional | **No** (4.1-4.6, 3.6) | Select by name, N columns, geometry from N; `fig_grouped_main:121-123` tag regex (4.17). The saved `std` row is the factor (2.10). |
| `fig_main3_grouped` | `fig_grouped_main:39-42` -> `fig_main3_stability(run_dir)` -> `fig_eig_heatmap(...)` (`fig_main3_stability.m:5`) plus the `local_lyapunov` stage | Run data, positional, two stages | **No** (4.7-4.13) | Pass the option to `fig_eig_heatmap` and select `D.results` by name; place copied tiles by name, not by `Layout.Tile`; `:124-126` tag regex (4.17). |
| `fig_main4_grouped` | `fig_grouped_main:43-49` -> `fig_sensitivity_medians` (`figs(1)`, 6 axes) and `fig_lle_vs_rate` (`figs(1)`, N axes) -> `style_main4_grouped` | Preset (medians) and run data (rate) | **No** in the style (4.19-4.22, 3.7) | Pass the option to both inner figures; the style takes names from the legend or the selection; `6 + N` axes. |
| `fig_main5_grouped` | `fig_grouped_main:50-80` -> `fig_transient_amplification` (`figs(2)`, the margin) plus an inline gain row from `transient_gain_data.mat` -> `style_main5_grouped` | Run data | **No** (4.14, 4.15, 4.23-4.26, 3.8) | Pass the option to `fig_transient_amplification`; select `D.results` by name at `:63`; layout `1 x N` at `:57`; drop the tick overwrite at `:54` or make it name-driven. |
| `fig_main6_grouped` | `fig_grouped_main:81-98` -> `fig_sfa_EOC_allStd` (both figures) and `fig_memory_capacity` (`figs(1)`, 3 axes) -> `style_main6_grouped` | Run data | **No** (4.16, 4.27-4.29) | Depends on the tau sweep carrying both MTS (4.30). The 3-set shows the chain's tau sweep; the 4-set both. Style: count tau axes from the data, titles by name, ticks by name. |
| `fig_main7_grouped` | `fig_grouped_main:99-102` (external PNG or a placeholder) | None | n/a | Accept and ignore. |
| `fig_main8_grouped` | `fig_grouped_main:103-107` -> `fig_main8_psd(source_fig_root)` | `bursting_psd.fig` written by `fig_stim_engages_adaptation` in the same pass | n/a | None, provided `source_fig_root` is this pass's root (5.0). |

---

## 6. Proposed split into three work packages (step 4)

No file appears in two packages. Each package also writes the `'conditions'`
argument into every file it owns, and runs its figures in both sets (the
step-4 acceptance in Plan A 7.3).

### Reserved for other workers (not in any package)

- `src/model/SRNNCellTypePairs.m`: `plot` (`:1791`), `plot_celltypes`
  (`:1849`), `collapse_routes` (`:3220-3237`), and the factor static itself
  (step 1).
- `src/model/SRNN_ESN_reservoir.m`: `route_signature` (`:996-1019`) and the
  description after it (`:1020-1038`), `assert_route_redundancy` (`:943-994`),
  `plot_esn_timeseries` (depression panel `:662-689`), `:547`
  (step 3, `matlab-fable-2`), with `scripts/tests/test_esn_route_redundancy.m`.
- Step 3, `matlab-opus-2` (presets and bundle): `srnn_condition_titles.m`,
  `full_adaptation_condition.m` (3.5, including the tie rule),
  `manuscript_style.m` (3.2, 3.3), `srnn_param_preset.m`, the new bundle's
  config and run script (including calling `make_all_paper_figures` twice and
  rebuilding the grouped registry per `fig_root`, 5.0),
  `test_preset_conditions.m` (3.13), `make_preset_golden.m` and the fixture.
  Plan A 7.3 also lists "the grouped styles" for step 3; this note recommends
  leaving them to step 4 (WP3), because their only use of the name is the
  three-name lists step 4 must rewrite anyway (section 7, R11).

### WP1: selection plumbing, the sweep sheets and the tau sweep (size L)

`scripts/paper/make_all_paper_figures.m`; `src/figures/helpers/select_conditions.m`
(new); `src/figures/helpers/write_figure_manifest.m`;
`scripts/tests/test_select_conditions.m` (new);
`src/figures/replot/replot_sensitivity.m`;
`src/figures/replot/replot_param_space_analysis.m`;
`src/plotting/param_space_plots/load_and_make_unit_histograms.m` (shared
plotting code; the lead to approve);
`src/figures/fig_sensitivity_analysis_allStd.m`; `fig_sensitivity_medians.m`;
`fig_param_space_allStd.m`; `fig_EI_weights_param_space.m`;
`fig_lle_vs_rate.m`; `fig_local_vs_finite_lle.m`; `fig_sfa_EOC_allStd.m`;
`src/analysis/run_tau_sensitivity_analysis.m` (4.30-4.32);
`src/util/write_run_parameters_md.m` (4.33, and its full-condition pick at
`:341-363` once 3.5 lands).

About 17 files. Small edits in most; the substance is `select_conditions`, the
two replot helpers, the two-condition tau sweep and its figure. **Starts
first or publishes `select_conditions` early**: WP2 and WP3 call it.

### WP2: stage-data figures and the depression-factor sites in analysis and examples (size M)

`src/figures/fig_eig_heatmap.m`; `fig_eig_heatmap_imbalance.m`;
`fig_transient_amplification.m`; `fig_transient_gain.m`;
`fig_transient_gain_excursions.m`; `fig_lyapunov_spectrum.m`;
`fig_local_lyapunov.m`; `fig_numerics_verification.m`;
`fig_memory_capacity.m`; `fig_memory_capacity_example.m`; `fig_dc_lle.m`;
`fig_single_neuron_revised.m`; `fig_example_timeseries.m` (2.6);
`fig_adaptation_methods.m` (2.7); `src/analysis/run_paper_illustrations.m`
(2.10); `src/analysis/run_stim_evoked.m` (2.11);
`src/analysis/run_stim_neurons.m` (2.12, 3.11); and the accept-and-ignore
entries `fig_introductory_concepts.m`, `fig_FI_curve.m`,
`fig_stim_engages_adaptation.m`, `fig_SFA_steady_state.m`.

About 21 files, nearly all a few lines each: `R = R(idx)` by name, titles by
name, factor static. Optional, if the lead wants the stimulation figures
regenerated on the chain: `fig_sine_paired_stim.m`, `fig_sine_stim_poster.m`,
`fig_sine_stim_neurons.m` (3.4, 3.10, 4.37); otherwise they stay untouched.

### WP3: grouped composites and the analytic depression figures (size L)

`src/figures/fig_grouped_main.m`; `src/figures/helpers/style_main4_grouped.m`;
`style_main5_grouped.m`; `style_main6_grouped.m`;
`src/figures/fig_main2_dynamics.m`; `fig_main3_stability.m`;
`fig_main1_concepts.m`; `fig_STD_steady_state.m`;
`src/figures/helpers/write_manuscript_tables.m`; `fig_doc_tables.m`;
`fig_STD_strength_matching.m` (product study; optional, only to stay correct
if ever pointed at a chain).

About 11 files, but the largest share of real layout and physics: N-column
geometry for groups 2-6, the (C4) steady state and a numerical step response,
and the chain equations in the tables. Depends on WP1's `select_conditions`
and on the inner figures' axes counts from WP1 and WP2 (`fig_lle_vs_rate`,
`fig_transient_amplification`, `fig_sfa_EOC_allStd`, `fig_memory_capacity`),
so WP3 should key on names and counts read from the figures, not on numbers.
`fig_main1_grouped.m` .. `fig_main8_grouped.m` and `fig_main8_psd.m` need no
edit and are in no package.

---

## 7. Raised, not fixed

- **R1. Silent wrong numbers on a chain route.** `fig_STD_steady_state.m:99-105`,
  `fig_main1_concepts.m:63-64`, `fig_STD_strength_matching.m:224, 229` and
  `write_manuscript_tables.m:282` all accept a chain route without error and
  compute the product of two independent pools (scalar `tau_rel` broadcasts).
  Each should branch on `route.model`, and error on an unknown model, rather
  than assume the product.
- **R2. `write_run_parameters_md.m:726-735`** prints only the first
  condition's `synapse_config` and states it is identical for every condition
  with depressing routes. That is already false for `sfa1_std1` (one
  `tau_rec`) against `sfa3_std2` (two). Existing defect.
- **R3. Stage data freeze display titles at compute time**, and the figures
  draw the saved titles: `run_eig_heatmap.m:110, 204`,
  `run_paper_illustrations.m:60`, `run_local_lyapunov.m:88`,
  `run_lyapunov_spectrum.m:104`, `run_numerics_verification.m:240`,
  `run_transient_gain.m:149`, `run_stim_evoked.m:120`. This contradicts the
  keyed-by-name design of `srnn_condition_titles` (retitle without
  recomputing). After `sfa3_std2` is retitled, every figure regenerated from an
  existing run keeps "Multiple-Timescale Adaptation" on the **product**.
  Figures should title from the saved *name*.
- **R4. Green for product MTS collides with green already in use.**
  `manuscript_style.m:56` makes the lambda_1 = 0 line green and `:63-64` says
  green is avoided in the condition palette for that reason;
  `style_main6_grouped.m:38` gives the STD state block bluish green
  (`[0 158 115]/255`); `fig_STD_steady_state.m:147` uses green for the
  undepressed reference. A decision for Tom on the exact green.
- **R5. Tick-label overwrites cycle silently.** `fig_grouped_main.m:54, 93` and
  `style_main6_grouped.m:52` assign three labels; with four ticks MATLAB
  repeats them (a wrong label, no error).
- **R6. The grouped registry bakes in `fig_root`** (`mu7revisedMedium_config.m:160`,
  `grouped_figure_registry.m:13`): two passes into different roots need it
  rebuilt per pass, or group 8 reads the other pass's (or a missing)
  `bursting_psd.fig`.
- **R7. `call_figure`** (`make_all_paper_figures.m:257-276`) tolerates only an
  undeclared `run_dir`; an undeclared `conditions` fails the entry. Hence 5.0.
- **R8. `full_adaptation_condition` tie** (`:50`): `sfa3_stdc2` and
  `sfa3_std2` tie on both keys once `stdc` parses; today's strict `>` would
  return whichever the preset lists first. Make the preference explicit.
- **R9. Memory-capacity statistics.** With four conditions
  `run_memory_capacity.m:241` tests six pairs, with no multiple-comparison
  correction; the three-condition figure must show only the three pairs among
  its conditions, and the four-condition figure's table should say how many
  tests were run.
- **R10. `fig_eig_heatmap.m:119-120`** puts four conditions in 2 x 2, which is
  not "a fourth column"; `fig_main3_stability` copies its tiles positionally.
- **R11. Step 3 and the grouped styles.** Plan A 7.3 step 3 gives "the grouped
  styles" to the preset worker for their name maps; the only name use in them
  is the hard-coded three-name lists (3.7, 3.8), which step 4 rewrites. Doing
  both in step 4 touches each file once.
- **R12. `fig_sfa_EOC_allStd.m:110`** draws `cond_names{1}` with only a
  warning when the tau sweep has more conditions: with two MTS it would show
  whichever the preset lists first.
- **R13. `replot_sensitivity.m:54, 86, 108, 111`** print with bare `fprintf`
  (bypasses `vprintf`/`vlog`), and `:107` calls `close all` inside a figure
  pass. Existing.
- **R14. `test_preset_conditions.m:353-354`** (and the same patterns in
  `full_adaptation_condition`) do not parse `sfa3_stdc2`, so the "names must
  not lie" check skips chain conditions unless extended. `std_count` (`:389`)
  counts `numel(tau_rec)`, which is 2 for a chain and so agrees with the name.
- **R15. Order of the step-3 edits.** `validate_preset_conditions` refuses a
  condition name with no title, so `srnn_condition_titles` must learn
  `sfa3_stdc2` in the same commit as the preset, or before it.
- **R16. Decisions the step-4 briefs need from the lead or Tom:**
  whether the chain bundle's single-neuron preset gains a chain column
  (`fig_single_neuron_revised`, 4.34); what the four-condition version of the
  group-1 concept figure shows (one STD panel or two); what
  `fig_STD_steady_state` shows in each set; short labels for the two MTS
  (3.3); and whether `fig_eig_heatmap` should be `1 x N`.

---

## 8. Not in scope (counts only)

Hits per folder, `.m` files, recursive. "params" is the section 1 pattern;
"prod" is `prod(`; "three" is `{'no_adaptation'` or `1:3`.

| Folder | Files | params | prod | `sfa3_std2` | three |
|---|---|---|---|---|---|
| `scripts/paper/archive/` | 36 | 3 | 0 | 14 | 0 |
| `scripts/examples/` | 23 | 107 | 7 | 2 | 0 |
| `scripts/explorations/` | 19 | 86 | 15 | 0 | 3 |
| `scripts/stim_figures/` | 2 | 0 | 0 | 0 | 0 |

(`scripts/tools/` holds one PowerShell script and no hits.)

---

## Counts

| Section | Entries |
|---|---|
| 1. Reads of depression parameters | 30 (1.1-1.30) |
| 2. Sites computing the factor | 20 (2.1-2.20) |
| 3. Sites naming `sfa3_std2` | 15 numbered (3.1-3.15; a: 4, b: 8, c: 3) plus 22 comment-only sites |
| 4. Three-condition (or one-condition) assumptions | 39 (4.1-4.39) |
| 5. Registry entries | 28 individual + 8 grouped |
