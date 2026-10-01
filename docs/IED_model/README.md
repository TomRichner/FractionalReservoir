# Exploratory IED-like network model

Sequential single-network pilot requested September 30, 2026. The aim is brief,
spontaneous, partially recruited events, rather than a globally synchronized
burst train. This is a model-morphology exploration, not clinical IED validation.

## Selected working base

On October 1, Tom asked to select one and stop the search. **Run13** is the
working base, frozen explicitly in `ied_base_config.m`. It uses 350 neurons
(245 E/105 I), five weakly coupled 70-neuron modules, one SFA timescale per
neuron, and one depression component on E outputs only. E and I SFA medians
are 0.40 and 0.75 s, with log-normal SDs 0.70 and 0.50; total strength is 0.5
on both types. Setpoints have mean0.45/SD0.15 on both types. Wiener amplitude
is0.025; external input is zero.

Within-module connection probability is0.30; between-module probability is
0.003. Absolute Gaussian weight means are E0.60/I−0.40, SD0.08, and bridges
are scaled by0.5. The realized mean indegree is21.37. The frozen constructor
also retains the old random-network placeholder parameters; the replacement
matrix above is the connectivity actually simulated.

Run13 produced11 candidates in50 s, median recruitment8.6%, median dominant-
group waveform FWHM195 ms, and negligible rate saturation (0.011% of neuron-time
samples). Within-module x correlation was0.241, versus0.005 between modules.
Ten of11 detected events involved a single group. The detector recruitment
width was108 ms; this is a different measure from waveform FWHM.
It is the most useful compromise among the tested networks,
not a finished sparse/nonrhythmic IED model: smaller background humps remain.
Run15 retained localized events on one additional realization of these
physical parameters. The later embedded-focus design either spread activity
too broadly (run18) or produced no events with weak bridges (run19).

The existing result is in `data/IED_model/20260930/run13/run.mat`. To replay
the frozen base later without overwriting it:

```matlab
setup_paths();
cfg = ied_base_config();
result = run_ied_experiment(cfg,'base_recheck');
```

Use a new tag for each replay. No new simulation was run merely to freeze it.

The unchanged dynamics are defined in
[the authoritative equations](../EquationsParametersDocs/Equations_stability_paper.md)
and implemented by `SRNNCellTypePairs`. Random runs call that class directly.
`IEDExplorationNetwork` inherits the same equations and only exposes a checked
replacement for a built recurrent weight matrix in modular experiments.
The embedded-focus trials also use a checked per-neuron setpoint replacement;
both replacements rebuild the cached parameters and activation handles.

From the FractionalReservoir root, through the MATLAB MCP:

```matlab
setup_paths();
result = run_ied_experiment(1);
```

Each ID is frozen in `scripts/explorations/IED_model/ied_run_config.m`. Completed
IDs refuse to overwrite their saved trajectory. Start with mu7-derived 1TS SFA
and 1TS STD, 350 neurons, no external input, SRA1 integration at 400 Hz,
input-referred Wiener amplitude 0.025, a 10-s settling/alignment interval and a
50-s observation. Connectivity normalization remains pinned to the original
500-neuron/100-indegree reference for the initial random network.

`tau_a_spread` is a log-normal spread, not an arithmetic time-constant SD.
Realized mean/SD are recorded. Setpoint heterogeneity is normal. Modular
experiments specify absolute Gaussian weight means/SD, within/between edge
probabilities, and bridge scaling; they do not covertly renormalize the result.
Group labels in random networks are diagnostic partitions and are not evidence
of anatomical modules. They contain both E and I neurons.
Ten-module trials use nominal35-neuron groups; per-type rounding gives slightly
unequal realized sizes. Embedded trials have group1 as a245-neuron sparse
background and groups2–4 as three35-neuron dense foci.

## Candidate detector and metrics

Detector v1 smooths dendritic states over 25 ms, uses each neuron's median plus
the larger of 3 robust SD or 0.12 a.u., and detects peaks of group recruitment
above 20% with prominence 10%. Peak widths must be 25 ms–1 s, with 150-ms peak
separation. Group detections within 120 ms are consolidated. Recruitment is the
fraction of all neurons crossing threshold within ±100 ms of the event. Events
recruiting more than 60% are flagged as global. Longer group excursions are
counted separately. The threshold rule is fixed for the pilot; relative
thresholds depend on each run's observed variability and can miss persistently
active or highly irregular regimes. Visual review is essential.

Reports include event count/rate, width, recruitment, global fraction,
inter-event variability, mean rate, saturation, and mean pairwise x correlation.
Top-8 QR growth is tracked in runs01–05; top-25 from run06 onward following
Tom's request, without automatic K expansion. Full 400-Hz state trajectories
and tangent parameters are retained from run06 for later reanalysis. The sum of positive
local tracked rates is a truncated expansion diagnostic in bits/s, not entropy
or information transmission. Event windows (−100 to +200 ms) are compared with
remaining segments descriptively. Alignment and finite observation limit
interpretation of these exploratory exponents; sums with different K are not
directly comparable. Coincident growth is not proof
that it causes an event.

## Outputs

- [Chronological progress](progress_2026_09_30.md): all configurations, metrics,
  plots, observations, and reasons for the next change.
- `data/IED_model/20260930/runNN/run.mat`: saved numerical results, ignored by git.
- `figs/IED_model/20260930/runNN/`: overview, event detail and native model plot,
  ignored by git. Links in the chronological log work on this machine.
- Scripts and Markdown are committed; generated figures/data stay local.
- `summary_2026_09_30.csv` and `events_run13.csv` preserve compact numerical
  comparisons. `summarize_ied_pilot(1:19,13)` rebuilds the saved-data summaries
  and plots without rerunning any network.

One exploratory realization per configuration establishes a candidate regime,
not robustness across networks. No stimulation test is included in this pilot.
