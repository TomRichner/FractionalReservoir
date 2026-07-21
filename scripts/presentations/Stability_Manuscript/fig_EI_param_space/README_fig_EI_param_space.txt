Stability_Manuscript figure: Parameter-space distributions (E:I colored)
=======================================================================

Generated: 21-Jul-2026 15:58:34
By script: Fig_EI_param_space.m

HOW IT WAS MADE
  Presentation replot -- no simulation is re-run. load_and_make_unit_histograms
  builds per-metric (LLE + mean_rate) 1x4 histograms where each bar is a stack
  of per-network patches colored by f (blue_gray_red_colormap; blue = low f /
  inhibition-dominated, red = high f / excitation-dominated), LLERange [-1.5,1.5],
  probability-normalized. Those axes are copied into a single 2x5 figure:
    row 1 = LLE ("Growth Rate", green dashed zero line)
    row 2 = mean firing rate
    columns 1-4 = No Adaptation, SFA, STD, SFA+STD
    column 5   = f-gradient colorbar (upper-right cell); lower-right cell empty
  Cleanups: condition titles raised into column-header position above the top
  row (not bold, enlarged to match the sensitivity figure); extra row spacing;
  y-ticks reduced (row 1: 0, 0.5; row 2: 0, 0.25); vertical gray column
  dividers; panel letters (a)..(h) on the 8 data panels only; fonts matched to
  the MC/sensitivity figures; y-axes linked within each row. The embedded
  colorbar encodes f as an E:I ratio
  (ticks 1:3, 1:2, 2:3, 1:1, 3:2, 2:1, 3:1). See git_provenance.txt.

SOURCE RUN
  C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\param_space\run_all_jul_06_26_22_00
  param_space subfolder used:
    param_space_test_refactor_nLevs_4_jul_07_26_10_44

FIGURES PRODUCED (in this folder)
  Fig_EI_ParamSpace.png / .svg / .fig   (2x5: f-colored distributions + colorbar)
