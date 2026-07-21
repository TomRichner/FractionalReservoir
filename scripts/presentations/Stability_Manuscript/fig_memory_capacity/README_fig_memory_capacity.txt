Stability_Manuscript figure: Memory Capacity
============================================

Generated: 21-Jul-2026 16:16:47
By script: Fig_memory_capacity.m

HOW IT WAS MADE
  This is a presentation replot -- no simulation is re-run. The script
  loads a saved looped_memory_capacity.m results file and calls
  replot_memory_capacity -> plot_memory_capacity to regenerate the
  figures into this folder. See git_provenance.txt for the exact commit.

SOURCE DATASET
  C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\memory_capacity\paper_ready\MC_sample_hold_20260706_194403_trials50_results.mat
  run_tag: MC_sample_hold_20260706_194403_trials50

KEY SETTINGS (from the saved run)
  input_type    : sample_hold
  n (neurons)   : 300
  fs (Hz)       : 200
  n_trials      : 50
  level_of_chaos: 2.5
  T_hold (s)    : 0.3
  conditions    : Baseline, SFA, STD, SFA+STD

FIGURES PRODUCED (in this folder)
  Fig_Memory_Capacity_figure_3.png
  Fig_Memory_Capacity_figure_3.svg
  Fig_Memory_Capacity_f_3.fig

  Only the final combined figure is saved. Fig1 (paired total-MC +
  memory-horizon distributions) and Fig2 (per-delay R^2(d) and cumulative
  MC) are shown for reference during the replot but not written to disk.
  Fig3 = 1x3 combined panel: (a) cumulative MC, (b) per-delay R^2,
         (c) horizon paired trials  [paper-ready]
