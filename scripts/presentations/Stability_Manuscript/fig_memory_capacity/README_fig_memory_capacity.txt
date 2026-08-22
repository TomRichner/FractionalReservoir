Stability_Manuscript figure: memory capacity
============================================

Generated: 22-Aug-2026 11:53:30
By:        fig_memory_capacity.m

WHAT IT SHOWS
  A 1 x 3 strip: (a) cumulative memory capacity against delay, (b) per-delay
  R^2, (c) memory horizon paired across trials. All four adaptation
  conditions, with bootstrap 95% confidence bands.

HOW IT WAS MADE
  Presentation replot -- no simulation is re-run. A saved run_memory_capacity
  ensemble is loaded and plot_memory_capacity_combined assembles the strip
  from the same summary statistics the working figures use. Fig1 and Fig2 are
  shown for reference but not written to disk.

SOURCE
  mat_file  C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\memory_capacity\MC_sample_hold_20260822_103708_trials4_results.mat
  run_tag   MC_sample_hold_20260822_103708_trials4
  preset    mc_esn
  run_mode  fast

PARAMETERS AS RUN
  n               300
  f               0.6
  fs              200
  level_of_chaos  2
  tau_d           0.1
  activation      logistic
  S_c             0.35
  input_type      sample_hold
  T_hold          0.3
  readout_signal  synaptic
  n_trials        4
  T_train_sec     60
  T_test_sec      30
  d_max_sec       5
  ode_solver      rk4
  std_zero_floor  false

FIGURES PRODUCED (in this folder)
  Fig_Memory_Capacity.png
  Fig_Memory_Capacity.svg
  Fig_Memory_Capacity.fig

READING THIS FIGURE
  MODEL CLASS. These figures use SRNN_ESN_reservoir, which subclasses
  SRNNModel2, so the memory-capacity network is NOT the SRNNCellTypePairs
  network every other figure in the paper shows. That is a structural
  constraint -- the ESN readout has not been ported to the Pairs class -- and
  the methods section must say so.

