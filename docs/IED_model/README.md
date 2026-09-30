# Exploratory IED-like network model

Sequential single-network pilot requested September 30, 2026. The aim is brief,
spontaneous, partially recruited events, rather than a globally synchronized
burst train. This is a model-morphology exploration, not clinical IED validation.

The unchanged dynamics are defined in
[the authoritative equations](../EquationsParametersDocs/Equations_stability_paper.md)
and implemented by `SRNNCellTypePairs`. Random runs call that class directly.
`IEDExplorationNetwork` inherits the same equations and only exposes a checked
replacement for a built recurrent weight matrix in modular experiments.

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

One exploratory realization per configuration establishes a candidate regime,
not robustness across networks. No stimulation test is included in this pilot.
