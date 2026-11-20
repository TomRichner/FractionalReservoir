# Plan: updating Figure 3 (switching input, expansion, eigenvalue density)

**Not to be executed until Tom says go** (Tom, 2026-11-20). This file collects
what has been decided about the next version of the paper's Figure 3
(`fig:eig_heatmap`; composite `fig_main3_grouped`) and what is still open, so
that one run can apply all of it. Decisions are also recorded in
`../StochasticPlasticDynamicalSystemPaper/docs/final_draft_to_do.md`; the code
history is `ISSUES.md` FD-033 and FD-034.

## Decided

### 1. Timing of the switching run (Tom, 2026-11-20)

| Interval | What happens |
|---|---|
| −15 to −10 s | The network runs with no stimulus; nothing is tracked. |
| −10 to 0 s | Alignment: the tracked directions are carried and re-orthogonalized, still with no stimulus. Their rates do not count. |
| 0 to 40 s | Accumulation: the stimulus is on (eight 5-s steps) and everything the figure reports comes from this window. |

Today the stage aligns over −14 to 1 s and accumulates over 1 to 40 s, so the
first second of stimulus is outside λ₁ and inside the alignment, while panel B
is drawn from 0 s. The reason for the 1-s start is not recorded (it came in
with commit `609afd0`, 2026-09-15). The new timing matches what the
single-run illustration figures already use (alignment over −10 to 0 s).

In settings, for `mu7chain_config.m` (lines 242 to 244 today):

```matlab
cfg.local_lyapunov_accumulation_start_s = 0;    % was 1
cfg.local_lyapunov_warmup_s             = 10;   % was 15
cfg.local_lyapunov_simulation_start_s   = -15;  % unchanged
```

### 2. One trajectory and 100 directions throughout (Tom, 2026-11-20)

Asked for and built as an opt-in (FD-034, commits `b90990f`, `3eb72b6`):
panel C's eigenvalue density and its Ritz values are taken at states of the
switching run, from the same top-100 pass whose local rates panel B sums; K =
100 in B, in the Ritz values and in the exponent ticks. Turned on with
`with_switch_spectra(cfg)`. A first version at `essential`, with the old
timing, is in `figs/fd034_mu7chain_essential_switchK100/`.

Whether this version replaces the paper's Figure 3 is not yet decided (see
Open questions).

### 3. Drawing (Tom, 2026-11-19)

Ritz values as a density, raw norm, an inset per condition, no explanatory
text on the panels, one title line at most; numbers go in the `_table.md`
sidecar and the caption.

## What executing this involves

1. **Change the three settings** above in `mu7chain_config.m`. Find and update
   whatever pins them (config tests, golden fixtures).
2. **Check what else reads the stage.** `fig_local_lyapunov` (marked
   `in_paper`) draws from the same `local_lyapunov` output, so it changes too.
3. **Make the 0-to-40 window consistent in the figure code**: panel B's axis
   already runs 0 to 40 s; the sampled states of panel C and the "shared
   pass" section of the sidecar must use the new window.
4. **Run at `essential`** with `with_switch_spectra`, and repeat the FD-034
   checks under the new timing: panel B against the diagonal of H at the
   sampled states, and the first half of the window against the second (the
   alignment is now 10 s, not 15, and the slowest adaptation timescale is
   10 s).
5. **Gates.** Option-off output against HEAD on the same machine; the
   covering tests; the exact regression tiers on the workstation (still owed
   for FD-032, FD-033 and FD-034).
6. **The paper's version** comes from the `fastmedium` run on the workstation
   (manuscript `docs/FEATURE_REQUESTS.md` FR-003), not from `essential`. If
   the FD-034 version is adopted, `with_switch_spectra` is called at the end
   of `mu7chain_config`, `figures/figure_sources.json` in the manuscript repo
   is repointed and `set_figs.py --sync` run.
7. **Caption and text** (manuscript repo, only after the figure is adopted):
   "accumulated over 1–40 s" becomes 0–40 s; "50 tracked QR directions
   ($q_{50}$)" becomes 100; "In separate baseline simulations" goes, since C
   is then the switching run; the quoted growth rates come from the run that
   made the figure; the colored points are described as instantaneous rates
   of the tracked directions, not as eigenvalues.

Workers: one MATLAB worker, pool allocated by the lead. On the MacBook Air
memory is tight (16 GB): pool 0 or small; only the dense eigendecompositions
use the pool, and they cost about the same serially.

## Open questions (Tom to decide)

- **Does the FD-034 version replace Figure 3**, or stay a variant?
- **K.** 100 throughout (Tom, 2026-11-20) against the earlier record of
  production K = 150 for this figure (figure review, 2026-11-17); the bundle
  itself sets 50 today.
- **Where the 40 states of panel C sit.** Today at 1, 2, …, 40 s, so every
  fifth state falls on a switch instant. Lead's proposal: 0.5, 1.5, …, 39.5 s,
  five per step and none on a switch.
- **The regression bundle.** `regression_fast_config.m` carries the old timing
  too. Changing it means new reference files, which is Tom's call. Lead's
  proposal: leave it as it is.
- **Sampling rate.** `run_local_lyapunov` sets 400 Hz itself; the uniform
  rerun is planned at 500 Hz. Fold that into the same change?
- **The repeated top-100 pass.** The stage makes the pass twice because the
  class cannot return H from its own pass. Removing the repeat needs a small
  edit in `SRNNCellTypePairs.m`, which carries another session's uncommitted
  changes.
- **The no-adaptation panel at K = 100.** With 100 of 500 directions most of
  its Ritz values pile at the neurons' decay rate (about −10 1/s), so its
  inset is one spike and its cloud is faint. Accept, or give that panel a
  different inset range?
- **Small things.** The λ₁ label printed on each panel of the stand-alone C
  figure; red for the no-adaptation condition elsewhere (manuscript FR-004).
