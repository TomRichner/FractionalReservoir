# Manuscript figure style guide

Consensus from Tom's grouped-figure feedback, 2026-09-15. Give this document to future figure agents together with the specific requested change. These are working presentation conventions; newer explicit feedback takes precedence. Figure-specific settings below are not universal defaults.

## Shared conventions

| Element | Convention |
|---|---|
| Fonts | **14 pt** for ticks, axis labels, legends, annotations, and panel letters. Set sizes explicitly, including MATLAB label/title multipliers. |
| Condition headings | Full names: **No adaptation**, **Single-timescale adaptation**, **Multiple-timescale adaptation**. Color each heading with its established condition color. STA/MTA were conversational shorthand, not requested title text. |
| Heading size | Slightly larger than body text when needed: Figures 2 and 3 use 16 pt. Allow wrapping without reducing readability. |
| Panel letters | Parenthesized capitals: **(A), (B), …**, in **regular weight**, not bold. Read across rows, then down. A letter may identify a related block of plots; see the mappings below. |
| Axes | Line width **1.0 pt**. Prefer open axes without top/right borders unless a box is explicitly requested. |
| Tick density | Use a few meaningful ticks. Do not fill small panels with redundant labels. Preserve the scientific range; figure-specific ticks are listed below. |
| Legends | One legend per shared encoding; avoid repeating it in each condition. Usually vertical and inside an upper-right corner. Shared legends go close to the panels they describe, not across the top as a supertitle. |
| Colorbars | Readable 14 pt labels/ticks; no bounding box. |
| Canvas | Compact enough that 14 pt text remains readable at the intended manuscript display size. Reduce excessive canvas dimensions instead of shrinking fonts to fit. Inspect exported images, not only MATLAB windows. |
| Spacing | Align related panels and tighten unnecessary gaps. Leave enough room for labels and legends without clipping. |
| Dividers | Where condition columns need separation, use faint gray vertical lines, width **2.5 pt**. These are layout dividers, not data or reference lines. |
| Mathematical labels | Use LaTeX symbols and correct subscripts, with descriptive words where requested. Link to the authoritative equations rather than inventing new notation. |

Use the existing condition-style/color helpers rather than defining a different condition palette in every figure. E/I neuron colors are a separate encoding: **E reddish, I bluish**, with enough lightness and hue variation to distinguish neurons within each population. Keep each neuron's color consistent across its displayed variables. Color alone should not carry all distinctions in multi-direction gain plots: retain different dashes and markers too.

## Figure-specific decisions

### Figure 1 — concepts

- Three row labels: **(A), (B), (C)**. Bring the eigenvalue and example-trace rows closer together. Row B has no time scale bar or time text: the example time is arbitrary up to rescaling.
- Bottom row, left to right: SFA-shifted sigmoids; SFA-shifted eigenvalue discs; STD-rescaled sigmoids; STD-shrunken discs.
- SFA curves and corresponding discs share a **black → orange** ramp. More orange means a larger leftward disc shift.
- STD curves and corresponding discs share **black → dark blue → teal**. Black is the tallest sigmoid/largest disc; teal is the most flattened sigmoid/smallest disc. The intermediate levels must be visibly distinct dark and medium blues, not near-black shades.
- Both sigmoid x labels read **Dendritic potential**.
- Maximum illustrated STD retains a **nonzero** sigmoid height. Derive and document the chosen attenuation from the model parameters and representative rate; do not treat the illustrative frozen-factor F-I curves as a solved steady-state response.
- Conceptual discs use custom horizontal/vertical axes with **Re/Im** labels, like the upper examples. Hide normal axis boxes/ticks. No disc titles, units, extra axis labels, or legends.
- All fonts 14 pt; example-trace y axes and sigmoid x/y axes width 1.0.

### Figure 2 — representative dynamics

- No supertitle. Six row letters **(A)–(F)**; three full-name, condition-colored column headings. Faint gray 2.5 pt column dividers.
- Descriptive LaTeX y labels: **Dendritic potential, x_i**; **Spike rate, r_i**; **Synaptic output, theta_i**; **SFA** plus the plotted feedback expression; **STD** plus the plotted depression expression; **lambda_1**, without units in this row. Render symbols as LaTeX, not literal underscore text.
- Use the correct expressions from [the model equations](EquationsParametersDocs/Equations_stability_paper.md); label the quantity actually plotted.
- Numeric plots show **half the previously displayed neurons**, selected deterministically with the same identities across rows. This is a reduction of the example subset, not half the network.
- Per-neuron traces: width **0.5 pt**, more distinguishable reddish/bluish E/I palettes. One E/I legend in the upper-right panel only.
- Finite leading-exponent curves: width **1.0 pt**. Row F zero line: **green, dashed, width 0.5 pt**.
- Dendritic-potential limits **[-6, 6]**, ticks **[-5, 0, 5]**. SFA limits **[0, 0.6]**, ticks **[0, 0.5]**. Exponent limits **[-5, 5]**, ticks **[-5, 0, 5]**.
- Hide every x axis. One scale bar in the **lower-left exponent panel**, from **t=5 to 15 s**, at **y=-4.8**; text below: **10 seconds**, 14 pt.
- Future numeric data should start finite-exponent accumulation near zero with prior alignment; the prepared configuration accumulates from t=0 and labels its first completed segment at 0.05 s. This does not authorize extending archived estimates backward.

### Figure 3 — spectra and local expansion

- Four row labels **(A)–(D)**. Full-name condition-colored titles above columns on **single lines**, **16 pt**; other text 14 pt. Gray column dividers.
- Show more than eight actual input neurons; the revised plot uses **24 fixed traces per condition**.
- Rows 3 and 4 share the display interval **[0, 60] s**. Requested future estimation uses **K=100** and accumulation **[1, 60] s**.
- Accumulated leading-exponent limits: no adaptation **[-0.5, 4]**; both adapting conditions **[-0.5, 0.5]**.
- Remove lambda_1 annotations from eigenvalue-density panels. Place the final accumulated value near **t=55 s** on the accumulated-exponent panel, just above or below its final value.
- Rows 2–4: no top/right box. Rows 2–3: hide lower x-axis spine/ticks. Row 3 retains a horizontal **LLE=0** line spanning 0–60 s. Row 4 x ticks point **outward/down**.
- Keep local-rate, accumulated-rate, and entropy interpretations distinct. Current saved data have **K=30** and accumulated estimates over **30–60 s**; wider display axes do not change those facts.

### Figure 4 — sensitivity and rate

- **(A)** labels the sensitivity block spanning rows 1–2; **(B)** labels row 3.
- Fonts 14 pt. Rows 1–2 y ticks: **[-1, 0, 1]**.
- Shared sensitivity legend: **upper right, vertically stacked**, not a horizontal supertitle.
- Row 3 gray summary line: **RGB [0.5, 0.5, 0.5]**, width **2.5 pt**.
- Ratio colorbar has no box and uses readable 14 pt text. Condition titles use condition colors.

### Figure 5 — non-normal margin and directional gain

- **1×4** layout: margin first, then the three condition-specific active-gain panels.
- **(A)** on subplot 1; **(B)** on subplot 2 identifies the gain block spanning subplots 2–4.
- Fonts 14 pt, axes width 1.0. Color margin condition tick labels and gain-panel titles by condition.
- Margin ticks: **[0, 25, 50, 75]**. Current limits **[0, 100]**: the upper bound was an implementation assumption because the request omitted it.
- Preserve all five gain curves and their distinguishing styles: worst case, noise RMS, excite E / inhibit I, excite E / excite I, and Lyapunov direction. Keep shared logarithmic gain limits across conditions.
- Center the shared gain legend **under subplot 3**, immediately below that panel's x label; do not center it between subplots 2 and 3.
- Plot the actual saved time horizon. A longer horizon requires new data, not stretching the existing one-second curves.

### Figure 6 — slow adaptation and memory

- Compact canvas, all fonts 14 pt, axes width 1.0. Current reduced size is **1100×780 px**, down from 1500×1050.
- Regular-weight letters: **(A)** upper-left timescale/LLE; **(B)** upper-right vector fractions; **(C)** lower-left capacity/reconstruction pair; **(D)** lower-right horizon.
- This four-letter mapping preserves **five scientific plots**; it is the implementation interpretation of the requested A–D labeling.
- Fewer y ticks on B, C, D. Remove **“(E and I)”** from the top two x labels.

### Figure 7 — learning

Apply shared typography and condition naming when revising this figure. No additional figure-specific style was agreed in this feedback round. The archived experiment remains provisional evidence.

### Figure 8 — spectral power

- Labels **(A), (B)**; no titles. Panel B is **only an empty bordered box** until its human source is supplied.
- Panel A y label: **Power spectral density of dendritic potential, x**.
- Axes width 1.0; fonts 14 pt. No-stim/stim legend **inside the upper-right of A**.
- Canvas width and height reduced to **60% of the earlier dimensions**, retaining font size. Current target is **780×456 px**. Set the target deterministically; repeated styling must not shrink it again.

## Instructions for implementation agents

1. Read the latest user request and relevant helper before editing. Apply the requested scope; do not silently restyle unrelated figures.
2. Prefer reusable, figure-specific helpers so both saved-data replotting and future runs inherit the same style. Coordinate shared-file edits and MATLAB access with other agents.
3. Run MATLAB through the **MATLAB MCP server only**. Formatting work does not authorize new simulations, model fitting, bootstrap recomputation, or a full analysis run. Inspect source figure functions before calling them; a function named as a plot may still compute statistics.
4. Prefer styling an existing native FIG when available. Preserve scientific XData/YData, pooling, windows, and limits unless explicitly requested otherwise.
5. Distinguish **implemented styling**, **available source data**, and **prepared future analysis**. In particular, Figure 2's current mu7 traces and Figure 8's selected model PSD survive as raster artwork. Individual overlapping neurons cannot be faithfully removed, recolored, or thinned from that flattened image. Apply such changes to future native data and record the current limitation.
6. Do not substitute another run's scientific trajectories because they are easier to edit. Any verified equivalent conceptual archive must have explicit provenance. Keep source hashes, calibration, and limitations in provenance notes rather than adding implementation prose to the figure.
7. Inspect the exported PNG at a realistic display size. Check text/legend clipping, overlapping labels, line visibility, color diversity, and alignment. Verify native font/axis properties and data preservation where applicable. Avoid tests that merely restate constant assignments; focused saved-figure/data checks are useful.
8. Styling helpers should be safe to apply repeatedly: no duplicate letters/legends, repeated canvas scaling, or accumulating offsets.
9. Refresh the combined report from existing outputs, preserving all eight entries and their provenance. Do not rerun analyses merely to refresh the report.
10. Commit completed, scoped source changes at logical steps, as authorized for this figure revision. Coordinate staging of shared files; generated MATLAB figures/reports remain gitignored. Report the commit and distinguish remaining data-dependent work.

## Entry points

- [Grouped dispatcher](../src/figures/fig_grouped_main.m)
- [Current-data / future-run instructions](../scripts/paper/README_mu7revised.md)
- Figure-specific implementations: `src/figures/fig_main1_concepts.m`, `fig_main2_dynamics.m`, `fig_main3_stability.m`, `fig_main8_psd.m`, and `src/figures/helpers/style_main*_grouped.m`.
- E/I palette helpers: [excitatory_colormap.m](../src/plotting/excitatory_colormap.m), [inhibitory_colormap.m](../src/plotting/inhibitory_colormap.m). Per-figure contrast adjustments should preserve E/I identity.

This guide records presentation choices. For the paper's selected figures and scientific narrative, consult the manuscript repository's `figure_plan_todo_status.md` and `Manuscript5_Outline.md`.
