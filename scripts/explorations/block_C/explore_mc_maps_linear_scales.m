% EXPLORE_MC_MAPS_LINEAR_SCALES Block C's gain-by-bias maps with linear color scales.
%
% Tom wanted to compare linear color scales with the report's square-root
% scale for memory capacity and symmetric-log scale for lambda_1. This
% draws fig_mc_gain_bias_maps on the main gain sweep with both scales
% linear and lambda_1 clipped at +-5 1/s; everything else (the spread row,
% the lambda_1 = 0 edge, the star, ticks, titles, layout) is the figure as
% it is. Once with the paper's three conditions, once with all four.
%
% Output: figs/explore_C/linear_scales/{three,four}_conditions/
%         Fig_MC_gain_bias_maps_linear.{png,svg,fig}
%
% See also: fig_mc_gain_bias_maps, read_mc_gain_sweep, symlog_lle

root = fileparts(which('setup_paths'));
D = read_mc_gain_sweep(fullfile(root, 'data', 'mc_gain_sweep_main', 'mc_gain_sweep', ...
    'mc_gain_sweep_results.mat'));
out_root = fullfile(root, 'figs', 'explore_C', 'linear_scales');

sets = {'three_conditions', {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'}
        'four_conditions',  {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2', 'sfa3_std2'}};
for k = 1:size(sets, 1)
    o = fig_mc_gain_bias_maps('data', D, 'conditions', sets{k, 2}, ...
        'mc_scale', 'linear', 'lle_scale', 'linear', 'lle_clip', 5, ...
        'fig_tag', 'Fig_MC_gain_bias_maps_linear', 'visible', false, ...
        'out_dir', fullfile(out_root, sets{k, 1}));
    fprintf('%s\n', o.files{:});
    close(o.figs);
end
