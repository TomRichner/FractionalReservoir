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
