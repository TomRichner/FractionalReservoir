Stability_Manuscript figure: Example memory capacity
====================================================

Generated: 08-Jul-2026 11:56:03
By script: Fig_memory_capacity_example.m

HOW IT WAS MADE
  Two-step, no re-simulation at plot time: compute_memory_capacity_example.m
  runs the memory-capacity protocol for the 4 adaptation conditions and saves
  the per-condition mc_results structs to mc_example_data.mat (gitignored).
  This script loads that file and renders the figure, so the look can be
  iterated without re-running the sim. See git_provenance.txt for the commit.

MODEL SETTINGS
  Match looped_memory_capacity.m (c_E = 0.5/3, sample_hold input, n = 300,
  f = 0.6, level_of_chaos = 2.5). See compute_memory_capacity_example.m.

FIGURE PRODUCED (in this folder)
  Fig_MC_Example.png / .svg / .fig
    (a) Cumulative Memory Capacity vs delay (0-10 s), all 4 conditions.
    (b) Per-delay R^2 vs delay (0-10 s), all 4 conditions (legend).
    Below: SFA+STD input reconstruction (target vs readout) for hold-delays
    [2 4 8 16] (delay indices), each titled with the delay in seconds and R^2;
    all reconstruction panels share y-limits [-0.6, 0.6].
