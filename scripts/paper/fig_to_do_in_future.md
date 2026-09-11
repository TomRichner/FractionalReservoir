# These are user notes, not a to do list for agents.

- add more PSD figures comparing no adaptation, 1 TS, 3 TS

- fix the sensitivity plots (don't force the hists into range.  Make the darkest color dark gray, not black, fix the xlabel overlap.  Go wiht -50% to +50%

- median sensitivity plots need different colors

- memory capacity needs to be redone on the final reference preset with 15 completed paired trials, verified source data, paired tests, bootstrap intervals, readable cumulative-capacity scaling, and preserved seed-level results and summary statistics

- make the jacobian density plot use f(x) = log(1+log(1+x)) (doing now)

- make separate folders per run+config

- config.m files need to be independent

- time series example should show 1 synaptic output per neuron if all have same STD routes

- time series should show depression for for all if all have same STD routes

- time series should show a subset of neurons

- time series should show local lya

- local lya density histogram over time

- modify the main time-series example to show local LLE and the accumulating finite-time LLE together

- compare the distribution of local LLE across time and trials with the distribution of finite-time LLE across trials; use this to test the proposed combination of transient expansion and long-interval stability in each adaptation condition

- Eliminate SRNNModel2.m

- Eigenvalue density plot for different imbalanced networks or networks with E only SFA

- add a stimulation-versus-no-stimulation IED analysis of the human CSCS recordings, with prespecified detection, artifact masking, channel aggregation, and participant-level summaries

- Verify precision of SRA1 vs ODE45 in noise-free (and maybe with-noise?) simulations and the similarity of LLE.  make a supplemental figure based on a 3 example time series plots with the base preset in the three conditions of the preset.  

- verify the benettin reshooting vs QR methods agree for SRNNCellTypePairs.m (it was done on the SRNNModel2.m).  this can also be done with the SRA1 vs ODE45 check.  

- only tau_a_E was swept for the tau sensitivity analysis.  it should be both E and I for runs in which both E and I have SFA.  We need to fix this.  Can this be done considering how the param space analysis class works? 

- could we go to a single cell type model (but with a dale's law weight matrix)?  this would reduce the std routes and possibly clean up the number of eignevalues in the jacobian.  results should be the same.  might make the tau_a sweep work without modification.  Is single cell type supported?  We ran into that during a refactor previously.

- quiet down the printing to the command line so it stops filling the agent's context through the MATLAB MCP server.  Now priority 0 in the ranked list below: a three-level `verbose` setting (`verbose` / `minimal` / `near-none`), `minimal` by default.

## Additional items from the Manuscript5 planning audit

- rebuild the introductory Sompolinsky figure as one MATLAB figure with eigenspectra above their corresponding stable, near-boundary, and chaotic state trajectories.  This is a small conceptual figure for the Introduction, not a new result, and the Sompolinsky network has no E/I cell types

- consider adding stochastic input to the introductory Sompolinsky examples so the trajectories illustrate the Goldilocks argument in a noisy system.  Preserve a clear progression from overly stable through near-boundary to chaotic dynamics, and decide whether gain or noise is held fixed across the three examples

- rebuild the adaptation-and-network Methods figure as one MATLAB figure: SFA/STD input-output transformations together with the representative multiple-timescale network time series.  Show a subset of neurons, consolidate identical route-level synaptic outputs, and add local and accumulating finite-time LLE traces

- directly analyze the relationship between mean firing rate and LLE across the joint parameter-space results.  Make an LLE-versus-mean-rate scatter or density plot separated by adaptation condition and preferably colored by E:I weight balance.  Report an appropriate nonlinear or rank association per condition.  The intended claim is that firing rate and stability interact but are not directly correlated or interchangeable measures

- normalize effective STD strength between the one-timescale and two-timescale conditions before the final comparison, or add a prespecified strength-matched control.  Plot the steady-state synaptic-output curves over the firing-rate range occupied by the network and record exactly what quantity is matched

- fix the longest-timescale SFA figure as well as the underlying E/I sweep: show the complete LLE distribution without clipping positive values, display the median and uncertainty across trials, keep the zero boundary visible, and save the resolved E and I time constants in the run metadata

- restore an automated central-finite-difference validation of the full analytic Jacobian, including active dendritic, SFA, and STD states.  Compare full matrices and eigenvalues at representative states in all three adaptation conditions, away from activation-function breakpoints, and save tolerances and errors

- make the imbalanced-network Jacobian-occupancy examples a defined comparison: use the three manuscript adaptation conditions, share seeds across conditions, use common real/imaginary limits and density normalization, and annotate each example with its matched finite-time LLE and mean firing rate.  Retain multiple seeds in the supplement if practical

- add uncertainty bands to the main median sensitivity plots and fix the E:I-weight parameter-space title/label overlap.  Consider combining the LLE, firing-rate, and joint E:I parameter-space summaries into one main figure while keeping full distributions in the supplement

- redesign or explicitly validate the bursting-network experiment.  Test the three manuscript adaptation conditions across held-out connectivity and noise seeds; do not select evaluation seeds based on favorable stimulation responses

- repeat the DC-input/LLE analysis in the bursting network and incorporate it into the stimulation figure.  Across a prespecified DC grid, quantify long-interval and local LLE, discharge rate/amplitude/duration/inter-discharge interval, mean rate, mean dendritic state, PSD, synchrony, and SFA/STD engagement

- decide whether the final bursting model remains the separately hand-tuned preset or becomes a minimal hyperexcitable derivative of the reference preset.  If redesigned, change one prespecified axis at a time, freeze the selected configuration on development seeds, and evaluate it on held-out seeds

- consider a 2-Hz pulsed-stimulation condition in addition to DC, with prespecified pulse width, polarity, amplitude, spatial projection, and charge/mean-input matching.  Keep DC as a mechanistic control if it helps isolate the operating-point shift

- if pulsed stimulation is used, separate neural input from a synthetic recording artifact and test recovery of latent PSD, discharge rate, and synchrony after the same artifact-removal logic used for the human data

- align model stimulation outcomes with the human analysis where possible: discharge-like event rate versus IED rate, absolute or log-ratio PSD with stimulation harmonics marked or masked, and a clearly defined synchrony measure.  Treat these as analogous rather than identical unless a shared observation model is implemented

- save a compact provenance artifact for every final figure: source commit, resolved preset/configuration, run mode, seeds, completed and failed trials, input data/run directory, summary statistics, and output paths

- the human stimulation-versus-no-stimulation IED analysis is a cross-repository reminder; its implementation belongs in cscs_stim_stability.  Include prespecified detection or blinded review, stimulation-artifact masking, channel/event aggregation, and participant-level summaries

- cross-repository reminder for cscs_stim_stability: strengthen the human PSD figure with participant-level curves, equal participant weighting, awake and sleep analyses kept separate, physical-unit verification, ICA/artifact-removal sensitivity checks, uncertainty across blocks, and stimulation harmonics marked or masked.  Copy the final figure into the manuscript repository rather than linking it by an absolute path

- cross-repository reminder for train-srnn: complete the MATLAB-aligned three-condition PyTorch rerun, replace the existing two-condition skip-connected figure, and consider a staged experiment that first trains an ESN-style linear readout with recurrent dynamics frozen and then continues with end-to-end recurrent training.  Preserve paired seeds, validation trajectories, AULC statistics, active parameter counts, final held-out results, and run provenance

## Ranked priorities for completing the paper's MATLAB results

This is a handoff for work on the Results of `StochasticPlasticDynamicalSystemPaper/Manuscript5.md`.  Assume the other computer has this repository and the complete `StochasticPlasticDynamicalSystemPaper` repository available side by side.  Before changing an analysis, read the live Results outline in `Manuscript5.md`, the status distinctions in `reports/reality_manuscript_suggested.md`, the model recommendations in `reports/suggest_updates_main_matlab_model.md`, and the relevant stimulation reports in `reports/suggested_updates_bursting_model.md` and `reports/suggested_updates_cscs_reanalysis.md`.  The paper's primary comparison is deliberately simple: no adaptation versus one-timescale adaptation versus multiple-timescale adaptation, with SFA and STD present together in both adapting conditions.  Brian asked that the paper not return to a large SFA-only/STD-only factorial comparison.  Do not broaden the current Results into a mechanistic separation of SFA and STD; that question is reserved for `StochasticPlasticDynamicalSystemPaper/reports/Future_investigations.md`.

The authoritative equations for the MATLAB model are `docs/EquationsParametersDocs/Equations_stability_paper_v2.md`: raw firing rate drives both SFA and STD, recurrent output is firing rate multiplied by the product of all active depression variables, SFA includes the fixed offset, and facilitation is not used.  The present paper preset is `celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25`.  It has 500 neurons, equal E and I populations, SFA on both cell types, two-timescale STD on all four connection routes in the multiple-timescale condition, heterogeneous neuronal setpoints, and additive dendritic noise.  Do not silently substitute an older E-only-SFA preset, `SRNNModel2`, or the separate hand-tuned bursting configuration for the reference model.

This order is dependency-aware.  Do not spend substantial compute polishing or rerunning downstream analyses until the model-defining decision in priority 1 is frozen.  For every final run, save a compact provenance artifact with the source commit, resolved preset and overrides, run mode, seeds, completed and failed trials, input run, summary statistics, and output paths.  A PNG without the data and configuration that produced it is not a completed result.

### 0. Give the whole code base a three-level `verbose` setting, defaulting to the AI-friendly level

**Why this comes before everything else.** Every analysis and figure in this list will be run, watched and debugged through the MATLAB MCP server, and everything MATLAB prints to the command window comes back into the agent's context verbatim.  The model classes print on every build and run ("W created: spectral radius…", "SRNNCellTypePairs built successfully", "Integration complete in…", "Largest Lyapunov Exponent…"), `ParamSpaceAnalysis2` prints per job, and a sweep or a 25-seed ensemble multiplies that by hundreds.  The result is context filled with chatter, more frequent compaction, and an agent that cannot see the few lines that matter.  `run_numerics_verification` had to wrap every `build()` and `run()` in `evalc` to be usable at all; that is a workaround at one call site, not a fix.  Doing this first makes every later priority cheaper to run and to supervise.

**What is wanted.** One setting, set in the config (`paper_config` and every `*_config.m`), threaded to every model class, analysis driver, stage and figure, and respected by all of them.  Three levels:

* `verbose` — everything that prints today.  For a human at the prompt, and for debugging one run.
* `minimal` — **the default.** The optimal level for an agent supervising a run over MCP: one line per stage or major step (what started, what finished, how long, where it wrote), one line per failure with identifier and message, the final summary table, and nothing per model, per job, per seed or per grid point.  Progress on long loops as a single line every N jobs or every few minutes, not per iteration.
* `near-none` — an even sparser level for long development cycles, to postpone compaction as far as possible: only errors and the final one-line outcome of each entry point (run directory, figure count, minutes).  No progress lines at all.

**Work required.** Add a `verbose` property to `SRNNModel2`, `SRNNCellTypePairs` (inherited by `SRNN_ESN_reservoir` and `SRNNNumericsProbe`), and `ParamSpaceAnalysis2`; route every `fprintf`/`disp` in the classes, the integrators, the stages under `src/analysis/`, the two entry points and the figure helpers through one small helper that checks the level, so the decision is made in one place.  Carry the level in `cfg` and in the `ctx` struct `resolve_run_context` builds, so a stage never reads it from anywhere else.  Workers must respect it too (the level travels with the model or `ctx` into `parfor`).  Remove the `evalc` wrappers in `run_numerics_verification` and elsewhere once the classes are quiet by default.  Keep warnings and errors untouched at every level: quiet means fewer lines, never hidden failures.

**Definition of done.** With the default level, a full `run_all_paper_analyses` + `make_all_paper_figures` at `'fast'` returns a command-window transcript short enough to read in one screen per stage; `verbose` reproduces today's output; `near-none` prints only the final lines; `test_run_modes` and the model tests still pass; and CLAUDE.md states the three levels and that `minimal` is the default.

### 1. Freeze the one-timescale versus multiple-timescale comparison, especially STD normalization

**Why this is first and how it affects the manuscript.** The main paper attributes differences among the three conditions to the temporal structure of adaptation.  SFA is already normalized: total SFA coupling is divided across its active timescales, so adding SFA timescales changes temporal structure without increasing its steady-state strength.  STD is not normalized in the same way.  The synaptic output is $r\prod_m b_m$, so adding a second depression variable changes both the number of recovery timescales and the magnitude and shape of steady-state depression.  With the current ratios, the two-timescale attenuation is approximately the square of the one-timescale attenuation.  If left unresolved, the paper can support a claim about the combined one-timescale and multiple-timescale architectures, but not a clean claim that timescale count alone caused the difference.

**Work required.** Read `StochasticPlasticDynamicalSystemPaper/reports/suggest_updates_main_matlab_model.md` before changing the model.  Decide explicitly what quantity should be matched: the steady-state synaptic-output curve over the empirically occupied firing-rate range, the output at one prespecified reference rate, or a new normalized multiplicative parameterization.  The least disruptive option is to retain the biologically motivated product model as primary and add a strength-matched control; changing the model equation itself would require updating the v2 equations, Methods, presets, tests, and every downstream analysis.  Plot the one- and two-timescale steady-state synaptic-output curves together and show the network's occupied rate distribution so the match is visible rather than asserted.

**Definition of done.** The choice and rationale are written down; the preset names make normalized and unnormalized conditions unambiguous; relevant model tests pass; the resolved parameters are saved; and the manuscript wording accurately distinguishes a timescale comparison from a combined timescale-plus-strength comparison.  Only then should the final sensitivity, timescale, memory, eigenvalue, and stimulation analyses be launched.

### 2. Finish the numerical-validity gate before the final production runs

**Why this matters to the manuscript.** `Manuscript5.md` states that 400-Hz SRA1 is sufficiently precise, agrees with `ode45` in noise-free simulations, remains consistent under time-step refinement, and produces a Benettin LLE agreeing with a full-spectrum QR calculation on smaller networks.  The Supplemental Methods also state that the full analytic Jacobian was verified against central finite differences.  These are validation claims, not optional implementation details: every stability result depends on the integrator, shared-noise path, perturbation reshooting, and Jacobian being correct.

**What already exists.** The intended ensemble runner is `scripts/paper/numerics_verification_trials_run.m`, configured by `numerics_verification_trials_config.m` and implemented by `src/analysis/run_numerics_verification.m`.  It runs four paired checks across the three manuscript conditions: noise-free SRA1 reshooting against tight-tolerance `ode45`; noisy SRA1 strong convergence on the same coarsened Brownian path; full-size Benettin LLE with `ode45` versus 400-Hz SRA1; and Benettin versus QR on a reduced network with the same preset physics.  The medium run uses five reshooting networks, 25 LLE/QR networks, and a 13-worker pool.  Run `wait_for_parpool(13)` first because the floating Parallel Computing Toolbox license may be unavailable.  The worker-side RNG is intentionally pinned to `twister` so a seed describes the same network on the client and workers.

**Remaining work.** Before examining the final results, record acceptance criteria for integration error, convergence behavior, paired LLE bias, and Benettin-versus-QR disagreement.  Do not choose thresholds after seeing the plots.  Add or separately preserve the missing full analytic-Jacobian check on representative states from all three manuscript conditions, away from the piecewise activation breakpoints.  `scripts/tests/test_SRNNCellTypePairs.m` contains a central-finite-difference check, but a remembered or transient console result is insufficient; save the configuration, maximum matrix error, eigenvalue comparison, and tolerance.  After the ensemble run, confirm that all three conditions contain exactly five reshoot trials and 25 LLE/QR trials, all required summaries are finite, and all three figure variants succeeded.

**Definition of done.** The `.mat` output, figures, provenance, acceptance criteria, and a short interpretation report are preserved together.  The report says whether each validation claim passed and distinguishes finite-time seed scatter from systematic numerical bias.  Any failed criterion is resolved in code or reflected honestly in the Methods before expensive production analyses proceed.

### 3. Establish the main stability result on the frozen reference model

**Manuscript connection.** This supplies the first Results subsection, currently titled “Multiple-timescale adaptation constrains dynamics near the stability boundary.”  The important quantity is the largest Lyapunov exponent measured along a changing nonlinear trajectory.  Instantaneous Jacobian eigenvalues are useful local descriptions of effective connectivity but cannot establish asymptotic stability when the Jacobian changes with state.  Keep that distinction explicit throughout the analysis and captions.

**Scientific question.** Does multiple-timescale adaptation produce a reproducible dynamical regime that is slightly stable over an extended interval while still allowing transient epochs of local expansion?  The current figures show a representative finite-time LLE and Jacobian occupancy, but they do not yet establish the desired stronger claim about transient expansion versus long-interval contraction across trials.

**Work required.** On the frozen three-condition preset, run paired network realizations with the same structural weights and compatible noise paths across conditions.  Modify the main time-series analysis to show the local LLE and the accumulating finite-time LLE together.  Across trials, retain the entire local-LLE time series rather than only its mean.  Compare the distribution of local LLE across time and trials with the distribution of the final finite-time LLE across trials.  Show how often and for how long local LLE is positive, while reporting whether the extended estimate is negative, near zero, or positive.  Use language such as “transient expansion” only for local positive intervals and reserve “stable” or “chaotic” for an explicitly stated finite-time window.

**Definition of done.** There is a representative trace that explains the quantities, an across-trial summary for all three conditions, prespecified analysis windows, and saved seed-level data.  The evidence must support either the narrow current claim—that longer timescales move the median LLE toward zero—or the stronger claim about transient expansion plus long-interval stability.  If it does not support the stronger claim, keep that idea in the Discussion rather than forcing the interpretation.

### 4. Complete the connectivity-robustness result and directly separate firing rate from stability

**Manuscript connection.** This supplies the subsection “Multiple-timescale adaptation preserves stability across changes in connectivity.”  The section has two linked messages.  First, the multiple-timescale condition should be less sensitive to changes in network structure and synaptic weights.  Second, mean firing rate and dynamical stability are related through the nonlinear operating point but are not directly correlated or interchangeable: a quiet or saturated network can be stable, and networks with similar rates can have different LLEs.  This is also why merely increasing inhibition is an incomplete description of network dynamics.

**Work required.** Rerun the seven one-dimensional sweeps—network size, E-cell fraction, global gain, and the four directed block means—on the final frozen preset with 15 paired networks per level.  Rerun the 64-point joint parameter-space sample with all three conditions paired on the same structural matrix.  Preserve full distributions, not just medians.  Add uncertainty bands to the compact median plots.  Construct a direct LLE-versus-mean-rate scatter, hexbin, or density view for each adaptation condition, preferably colored by realized E:I weight balance rather than only E-cell count.  Quantify the relationship with a prespecified rank or nonlinear association and inspect saturation/floor subgroups; do not summarize the point with Pearson correlation alone if the relationship is visibly nonlinear or multimodal.

**Interpretation constraints.** “Less sensitive” should refer to a defined measure—flatter median response, narrower distribution, a larger fraction remaining in a target LLE/rate region, or another stated statistic—not visual impression alone.  Do not call rate and stability orthogonal in a strict mathematical sense.  The intended conclusion is that they are distinct, interacting outcomes and that neither substitutes for the other.

**Definition of done.** The main figure communicates the robustness comparison and the rate-versus-stability distinction without relying on clipped values.  Full per-level distributions and overflow values remain available in the supplement.  Every plotted point can be traced to a resolved configuration, condition, and seed.

### 5. Correct and rerun the longest-adaptation-timescale analysis

**Manuscript connection.** This is the direct evidence for “Long adaptation timescales tune stable networks toward the edge of chaos.”  The paper is not using this analysis to separate the causal roles of SFA and STD.  It asks a narrower question within the combined multiple-timescale architecture: whether lengthening the slowest SFA component changes proximity to the stability boundary while the rest of the condition is held fixed.

**Known problems.** The earlier analysis changed the longest SFA time constant in only one cell type even though the final reference model places the same SFA ladder on E and I neurons.  The present plotting script also clips a substantial positive part of the LLE distribution, making the result appear more uniformly stable than it is.  The current manuscript Methods already describe an E-and-I sweep, so the final data must actually match that description.

**Work required.** Make the parameter-space machinery vary the slowest time constant for both E and I, keep the fastest value fixed, and recompute the intermediate value by the declared logarithmic spacing.  Save both resolved time-constant vectors at every grid point.  Use enough independent networks and a long enough LLE window to distinguish a shift in the median from finite-time scatter.  Redesign the plot so positive observations and overflow bins are visible, the zero boundary is prominent, and median plus uncertainty is shown over the full range.

**Definition of done.** A manifest proves that both cell types were swept; the full distribution is visible; the median trend and uncertainty are reported; and the manuscript makes only the claim the distribution supports.  Do not reuse the earlier one-cell-type figure under the new Methods description.

### 6. Redo and audit reservoir memory capacity on the final reference preset

**Manuscript connection.** This is the fixed-recurrent-weight half of “Multiple-timescale adaptation improves recurrent computation.”  It is a reservoir-computing analysis: recurrent and input weights are fixed, only a linear readout is trained, and delayed-input reconstruction measures how long the untrained dynamics retain usable information.  It is distinct from the PyTorch HalfCheetah experiment, where recurrent and dynamical parameters are trained end to end; the two analyses share a Results heading because both concern computation, not because they are the same task.

**Scientific question and claim boundary.** The defensible question is whether multiple-timescale adaptation extends fading memory in networks used as drawn, without tuning each recurrent matrix to a target spectral radius or LLE.  The current experiment does not measure the amount of tuning required, so avoid claims such as “requires little tuning” unless a separate tuning comparison is performed.

**Work required.** Run the final reference preset for 15 paired trials under all three conditions, sharing recurrent and input weights appropriately across conditions.  Verify the training/test split, readout signal, delay grid, horizon threshold, bootstrap procedure, and exact paired sign-flip tests.  Preserve the result `.mat`, compact seed-level CSV, configuration, and statistical summary.  Fix the current figure's excessive cumulative-capacity scale and make individual paired horizons readable.  Check that the plotted file comes from this final run rather than an older fallback or a different network size/preset.

**Definition of done.** All 15 paired trials completed; the source run and preset are unambiguous; uncertainty and paired statistics match the saved data; the main figure communicates per-delay performance, cumulative memory, and horizon; and the text describes an untuned fixed-reservoir comparison without conflating it with the PyTorch result.

### 7. Strengthen the Jacobian-occupancy evidence without treating it as the stability test

**Manuscript connection.** Jacobian occupancy visually connects adaptation to time-varying effective connectivity.  The reference heatmap is a representative example supporting the first stability subsection.  It complements—but never replaces—the trajectory-based LLE.  Even many instantaneous spectra cannot by themselves determine the stability of the nonlinear time-varying system.

**Work required.** Regenerate the reference example from the final frozen preset, using the same three conditions and one shared structural matrix.  Preserve the full active Jacobian for each condition; disabled SFA/STD states must be absent rather than retained as artificial zero modes.  Use common real/imaginary limits, binning, smoothing, and the same double-log density transform across panels.  Add a small, prespecified set of structurally inhibition-dominant, reference, and excitation-dominant examples, preferably by changing one interpretable block-mean coordinate or realized E:I weight-balance coordinate.  Use shared seeds across conditions and annotate every example with its matched finite-time LLE and mean rate.  More seeds can appear in the supplement, but do not cherry-pick visually attractive spectra.

**Definition of done.** The main figure clearly labels the heatmap as state occupancy from a representative realization, explains what the zero real-part line means locally, and points readers to the LLE for the stability conclusion.  The imbalance examples have recorded selection rules and provenance, and their visual interpretation is consistent with the quantitative sensitivity results.

### 8. Turn the bursting/stimulation demonstration into a replicated result

**Manuscript connection.** This supports the final Results subsection, “Sustained stimulation recruits adaptation and suppresses discharge-like dynamics.”  The present model result is descriptive: one hand-tuned 50-neuron realization produces recurrent discharge-like events without input and tonic activity with DC input.  It is not the standard reference operating point, does not compare all three manuscript conditions, and does not reproduce the 2-Hz clinical waveform.  Read `StochasticPlasticDynamicalSystemPaper/reports/suggested_updates_bursting_model.md` before redesigning it.

**Model decision.** Either retain the current bursting preset and label it explicitly as a separately configured illustrative hyperexcitable model, or construct a minimal derivative of the final reference preset by changing one prespecified axis at a time.  Do not imply that the reference model spontaneously discharges if it does not.  If redesigning, use development seeds to locate a robust discharge regime, freeze the parameter rule, and reserve independent connectivity/noise seeds for evaluation.  Avoid selecting evaluation seeds because stimulation works well on them.

**Work required.** Evaluate no adaptation, one-timescale adaptation, and multiple-timescale adaptation under matched unstimulated and stimulated conditions.  Quantify discharge-like event rate, amplitude, duration, and inter-event interval; mean firing rate and dendritic state; PSD; synchrony; and SFA/STD state engagement.  Define the discharge detector before comparing conditions and show sensitivity to its principal thresholds.  The current paper does not require SFA-only and STD-only conditions, and adding them should not derail the simplified three-condition story.

**Definition of done.** The stimulation effect is replicated across held-out seeds, the prevalence of spontaneous discharges is reported, all parameter differences from the reference preset are documented, and the result is described as model evidence rather than proof of the mechanism in patients.

### 9. Repeat DC-input/LLE in the bursting network and build the main stimulation figure

**Why this is separate from priority 8.** The existing DC/LLE analysis was performed on the reference network and shows how tonic input changes stability there, while the visually compelling discharge suppression comes from the separate bursting network.  Placing those unrelated analyses side by side would leave a mechanistic gap.  Repeating DC/LLE within the bursting experiment ties the stability metric directly to the observed disappearance of discharge-like events.

**Work required.** Across a prespecified DC-amplitude grid, measure local and accumulating/long-interval LLE together with the discharge, rate, PSD, synchrony, and adaptation-state outcomes defined in priority 8.  Use the same held-out seeds and analysis windows for every metric.  Include the three manuscript adaptation conditions so the data show whether input suppression depends on the adaptive architecture rather than merely forcing the network to another saturated fixed point.  Inspect both increasing and, if relevant, decreasing/recovery phases so rebound or hysteresis is not missed.

**Figure design.** The main composite should include one readable before/during-input time series, an across-seed curve or paired summary of discharge burden and LLE versus DC level, and compact PSD/synchrony summaries.  Mean firing rate must be shown because reduced low-frequency PSD can accompany increased tonic firing.  Additional seed traces and detector checks belong in the supplement.

**Definition of done.** The same model realizations support the statements about stimulation, LLE, discharge suppression, PSD, and synchrony; the DC level and settled windows are prespecified; and the figure distinguishes loss of bursting from simple silencing or saturation.

### 10. Rebuild the two introductory and methodological composite figures in MATLAB

**Figure 1: known Goldilocks concept in the Introduction.** Brian requested a small introductory figure showing that a Sompolinsky network can be overly stable, near the stability/chaos boundary, or chaotic as recurrent gain increases.  It is established background, not a new Result.  The conceptual Sompolinsky model has no E/I cell types and should not be described with Dale-law or E:I language.  Combine `fig_introductory_concepts/eigenspectra/panelA_eigenspectrum.png` with the corresponding `statetraces/panelA_bottom_traces.png`, aligning the three gain conditions column by column.  Consider adding modest stochastic input so the stable regime does not look trivially dead and the Goldilocks idea is visible in a stochastic system, but keep gain as the clear organizing variable and do not let changing noise create the apparent transition.

**Figure 2: model mechanisms in the Methods.** Combine the SFA/STD input-output illustration with the representative adapted-network time series.  This figure explains the model rather than establishing a Results claim.  Show only a readable subset of neurons, display one synaptic-output or depression trace per neuron when all routes are identical, label E and I only in this cell-type model, and add local plus accumulating finite-time LLE if the traces remain interpretable.  The current layout mockups are in `StochasticPlasticDynamicalSystemPaper/combined_fig_placeholders/`; MATLAB should ultimately generate the final composites directly.

**Definition of done.** Both figures are generated by the paper pipeline, use the manuscript palette and terminology, have captions that distinguish conceptual background from analyzed results, and replace the temporary manuscript placeholder paths without changing the intended section placement: Figure 1 in the Introduction and Figure 2 in Methods.

### 11. Treat pulse trains and synthetic stimulation artifacts as an extension unless required for the final claim

**Context and scope.** The current model uses a constant spatially uniform input because it cleanly tests whether sustained drive can engage adaptation without introducing a pulse artifact.  The human comparison uses 2-Hz chronic subthreshold cortical stimulation.  DC is therefore a mechanistic abstraction, not a reproduction of clinical stimulation.  The manuscript can retain that distinction, especially if stimulation is de-emphasized in the title.

**Optional extension.** If time permits after priorities 1–10, add a prespecified 2-Hz pulse condition with stated width, polarity, amplitude, spatial projection, and charge or mean-input matching to DC.  Model the neural input separately from the measurement artifact.  If comparing preprocessing with the human pipeline, create synthetic observed channels from known latent network states, add a stimulation-locked artifact with controlled amplitude and spatial mixing, apply the same artifact-removal logic, and quantify recovery of latent PSD, discharge rate, and synchrony.

**Guardrails and definition of done.** Do not interpret a removed 2-Hz line as a biological response, and do not confuse artifact removal with neural suppression.  This extension is complete only if latent truth, contaminated observation, and reconstructed observation are all retained and compared.  It should not delay the core stability, robustness, timescale, reservoir, and replicated-DC analyses if those are still incomplete.

### 12. Perform final figure assembly and manuscript-consistency checks only after analyses are frozen

**Why this is last.** Polishing figures before the model and data are frozen creates attractive but obsolete artifacts and makes it easy for the manuscript to mix runs.  Several current PNGs have clipped ranges, overlapping titles, inconsistent colors, placeholder captions, or insufficient provenance.  Some manuscript statements are prospective and must be reconciled with the actual completed analyses.

**Work required.** Use one stable condition palette and the exact names “No Adaptation,” “Single-Timescale Adaptation,” and “Multiple-Timescale Adaptation.”  Repair clipped values, overflow handling, label overlap, and illegible full-network traces.  Put concise summary panels in the main text and complete distributions, individual seeds, numerical validation, and detector sensitivity in the supplement.  Confirm every number in `Manuscript5.md` against the saved summary rather than transcribing from a plot.  Ensure each figure uses the final reference preset or is explicitly labeled as the conceptual Sompolinsky model or the separate hyperexcitable model.  Update the manuscript figure index and `Manuscript5_TOC.md` whenever figures or headings move.

**Cross-repository boundaries.** The final PyTorch three-condition run belongs in `train-srnn`; human PSD, participant weighting, artifact checks, synchrony, and IED detection belong in `cscs_stim_stability`.  This MATLAB repository may carry reminders and model-side outcomes, but it must not manufacture patient results or silently substitute a different learned-model equation.  The final CSCS figure should be copied into the manuscript repository instead of referenced through an absolute external path.

**Definition of done.** A clean figure build produces all intended files and a manifest; captions, Methods, Results, supplemental material, and saved statistics agree; temporary placeholders and absolute paths are gone; and another researcher can identify the code commit and data source for every panel without relying on this conversation.
