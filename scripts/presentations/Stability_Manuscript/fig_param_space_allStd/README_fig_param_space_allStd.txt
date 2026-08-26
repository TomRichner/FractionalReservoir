Stability_Manuscript figure: parameter-space distributions
==========================================================

Generated: 26-Aug-2026 01:16:57
By:        fig_param_space_allStd.m

WHAT IT SHOWS
  A 2 x N sheet. Row 1 is the distribution of the largest Lyapunov exponent
  over the whole parameter grid, with a green dashed line at zero; row 2 is
  the distribution of mean firing rate. One column per adaptation condition:
  No Adaptation, SFA, STD, SFA + STD.

HOW IT WAS MADE
  Presentation replot -- no simulation is re-run. replot_param_space_analysis
  regenerates psa.plot for LLE and mean_rate into a temporary replot folder;
  those axes are copied into a single combined figure, restyled, and the prep
  folder is deleted. Cleanups: LLE xlabel to lambda_1; condition titles only
  on the top row and not bold; vertical grey column dividers; y-axes linked
  within each row with sparse probability ticks.

SOURCE
  run_dir  C:\Users\m218089\Desktop\github_repos\FractionalReservoir\data\param_space\run_all_aug_25_26_22_08
  preset   celltype_pairs_Sc0p2_noise0p025_dualStd_7cond

FIGURES PRODUCED (in this folder)
  Fig_ParamSpace_allStd.png
  Fig_ParamSpace_allStd.svg
  Fig_ParamSpace_allStd.fig

READING THIS FIGURE
  THE LLE HISTOGRAM RANGE IS FIXED AT [-1.5, 1.5] inside
  ParamSpaceAnalysis2.plot and cannot be set from here. On the aug_14 preset
  the param-space LLEs spanned roughly -10 to +4.8, so the outermost bins
  accumulated a large share of the distribution. Check that against the run
  actually plotted before quoting the shape of the tails.

