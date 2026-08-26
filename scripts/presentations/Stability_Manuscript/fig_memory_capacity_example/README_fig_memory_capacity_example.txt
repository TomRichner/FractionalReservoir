Stability_Manuscript figure: example memory capacity
====================================================

Generated: 26-Aug-2026 14:56:17
By:        fig_memory_capacity_example.m

WHAT IT SHOWS
  (a) cumulative memory capacity against delay, all four adaptation
  conditions. (b) per-delay R^2, all four conditions. Below, input
  reconstruction for the SFA+STD condition at a few delays -- target against
  trained readout -- each panel titled with the delay in seconds and its R^2,
  all sharing y-limits.

HOW IT WAS MADE
  Two-step, no re-simulation at plot time. run_memory_capacity_example runs
  the protocol for the four conditions on ONE network and saves the per-delay
  R^2 plus the reconstruction traces to mc_example_data.mat; this function
  renders that file. The full mc_results also holds the complete state time
  series (~0.6 GB per condition), which the figure does not use and which is
  stripped before saving.

SOURCE
  data_file  C:\Users\m218089\Desktop\github_repos\FractionalReservoir\scripts\presentations\Stability_Manuscript\fig_memory_capacity_example\mc_example_data.mat
  preset     mc_esn
  run_mode   fast

PARAMETERS AS RUN
  preset_name     mc_esn
  run_mode        fast
  fs              200
  T_hold          0.3
  input_type      sample_hold
  readout_signal  synaptic
  T_train_sec     60
  T_test_sec      30
  T_wash_sec      10
  d_max_sec       5

FIGURES PRODUCED (in this folder)
  Fig_MC_Example.png
  Fig_MC_Example.svg
  Fig_MC_Example.fig

READING THIS FIGURE
  CONDITION COLOURS come from manuscript_style and are keyed BY NAME, so this
  figure and the ensemble memory-capacity figure cannot drift apart on them.

  MODEL CLASS. SRNN_ESN_reservoir subclasses SRNNModel2, so the
  memory-capacity network is NOT the SRNNCellTypePairs network the rest of the
  paper shows. The ESN readout has not been ported to the Pairs class; the
  methods section must say so.

