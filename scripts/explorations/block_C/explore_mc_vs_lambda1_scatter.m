% EXPLORE_MC_VS_LAMBDA1_SCATTER Memory capacity against lambda_1, every run of the main gain sweep.
%
% Tom asked how correlated memory capacity is with lambda_1, and whether
% the relation is simply linear. The block C maps show the two side by
% side on the gain-by-bias grid; this puts them on one pair of axes. One
% panel per condition, every run a point (11 gains x 7 biases x 2 seeds),
% x = lambda_1 with the task input on (lle_input), y = total memory
% capacity, colored by log10 gain so the path along gain can be followed.
%
% x is linear and clipped at +-5 1/s, as in Tom's linear map: a clipped
% run is drawn at the axis edge as an open marker and counted in the
% panel. The statistics in each panel use the UNCLIPPED values of every
% run with both values finite: Pearson r and Spearman rho over all runs,
% and Pearson r on the stable side (lambda_1 < 0) and on the unstable side
% (lambda_1 > 0) alone, each with its n. If the relation is peaked near
% zero rather than monotone, the two one-sided r have opposite signs and
% the overall r understates it.
%
% NOW A THIN WRAPPER (2026-11-13): Tom wanted the figure in the paper, so it
% was promoted to the pipeline figure fig_mc_vs_lambda1, which
% make_block_C_figures and the bundle mu7chain_gainbias_config draw. This
% script calls it with the exploration's settings and output folder; the
% one implementation is there.
%
% Output: figs/explore_C/mc_vs_lambda1/Fig_MC_vs_lambda1_scatter.{png,svg,fig}
%
% See also: fig_mc_vs_lambda1, fig_mc_gain_bias_maps, read_mc_gain_sweep,
%           explore_mc_maps_linear_scales

root = fileparts(which('setup_paths'));
out = fig_mc_vs_lambda1( ...
    'data_file', fullfile(root, 'data', 'mc_gain_sweep_main', 'mc_gain_sweep', 'mc_gain_sweep_results.mat'), ...
    'conditions', {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2', 'sfa3_std2'}, ...
    'lle_clip', 5, 'color_by', 'gain', ...
    'out_dir', fullfile(root, 'figs', 'explore_C', 'mc_vs_lambda1'), ...
    'fig_tag', 'Fig_MC_vs_lambda1_scatter', 'visible', false);
close(out.figs);
disp(out.stats);
fprintf('%s\n', out.files{:});
