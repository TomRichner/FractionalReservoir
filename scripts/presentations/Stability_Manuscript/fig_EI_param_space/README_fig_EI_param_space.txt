Stability_Manuscript figure: Parameter-space distributions (E:I colored)
=======================================================================

Generated: 07-Jul-2026 16:15:45
By script: Fig_EI_param_space.m

HOW IT WAS MADE
  Presentation replot -- no simulation is re-run. load_and_make_unit_histograms
  builds per-metric (LLE + mean_rate) 1x4 histograms where each bar is a stack
  of per-network patches colored by f (blue_gray_red_colormap; blue = low f /
  inhibition-dominated, red = high f / excitation-dominated), LLERange [-1.5,1.5],
  probability-normalized. Those axes are copied into a single 2x4 figure:
    row 1 = LLE ("Growth Rate", green dashed zero line)
    row 2 = mean firing rate
    columns = No Adaptation, SFA, STD, SFA+STD
  Cleanups: condition titles only on the top row (not bold); vertical gray
  column dividers; panel letters (a)..(h); fonts matched to the MC/sensitivity
  figures; y-axes linked within each row. A separate colorbar figure encodes f
  as an E:I ratio
  (ticks 1:3, 1:2, 2:3, 1:1, 3:2, 2:1, 3:1). See git_provenance.txt.

SOURCE RUN
  C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\param_space\run_all_jul_06_26_22_00
  param_space subfolder used:
    param_space_test_refactor_nLevs_4_jul_07_26_10_44

FIGURES PRODUCED (in this folder)
  Fig_EI_ParamSpace.png / .svg / .fig   (2x4 f-colored distributions)
  Fig_EI_Colorbar.png / .svg / .fig     (E:I-ratio gradient bar)
