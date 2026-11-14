function out = fig_mc_vs_lambda1(cfg)
% FIG_MC_VS_LAMBDA1 Memory capacity against lambda_1, every run of the gain sweep, one panel per condition (block C).
%
%   out = FIG_MC_VS_LAMBDA1('run_dir', d)
%   out = FIG_MC_VS_LAMBDA1('data_file', f)
%   out = FIG_MC_VS_LAMBDA1('data', D, 'save', false)   % D from read_mc_gain_sweep, or a table
%
% Tom's question: how correlated is memory capacity with lambda_1, and is
% the relation simply linear? The maps (fig_mc_gain_bias_maps) show the two
% side by side on the sweep's grid; this puts them on one pair of axes.
%
% It is fig_mc_sweep_scatter with x = lambda_1 (cfg.lle_field, default
% lle_input: with the memory task's input on), y = total memory capacity,
% colored by gain on the grid's own scale (cfg.color_by = 'gain', the default; the
% name of the result's second axis, 'bias' or 'self_feedback', colors by it
% instead; any quantity of mc_sweep_quantity is accepted), lambda_1 linear
% and clipped at
% +-cfg.lle_clip (default 5 1/s). cfg.size_by names a signed quantity drawn
% as the marker's area, filled above zero and a ring below (the
% self-coupling of a gain x self-coupling grid, the bias); empty, the
% default, draws every marker one size. The drawing, the statistics (Pearson r,
% Spearman rho, and Pearson r on each side of lambda_1 = 0) and out are
% those of fig_mc_sweep_scatter; see it.
%
% out: figs, files, source, and stats (a table: condition, n, pearson_r,
% spearman_rho, n_stable, r_stable, n_unstable, r_unstable, n_clipped,
% n_failed).
%
% See also: fig_mc_sweep_scatter, mc_sweep_quantity, fig_mc_gain_bias_maps,
%           read_mc_gain_sweep, explore_mc_vs_lambda1_scatter

arguments
    cfg.verbose     (1,:) char    = 'minimal'
    cfg.run_dir     (1,:) char    = ''
    cfg.data_file   (1,:) char    = ''
    cfg.data                      = []
    cfg.grid        (1,:) char    = ''   % with run_dir: the grid's tag in the run directory
    cfg.out_dir     (1,:) char    = ''
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.preset_name (1,:) char    = ''   % unused; accepted for a uniform call
    cfg.conditions                = {}
    cfg.lle_clip    (1,1) double  {mustBePositive} = 5   % |lambda_1| clipped here (1/s)
    cfg.lle_field   (1,:) char    {mustBeMember(cfg.lle_field, {'lle_input', 'lle_noinput'})} = 'lle_input'
    cfg.color_by    (1,:) char    = 'gain'   % a quantity of mc_sweep_quantity, which checks it
    cfg.size_by     (1,:) char    = ''       % a signed quantity drawn as marker area; empty = one size
    cfg.fig_tag     (1,:) char    = 'Fig_MC_vs_lambda1'
end

out = fig_mc_sweep_scatter('verbose', cfg.verbose, 'run_dir', cfg.run_dir, ...
    'data_file', cfg.data_file, 'data', cfg.data, 'grid', cfg.grid, ...
    'out_dir', cfg.out_dir, 'save', cfg.save, ...
    'visible', cfg.visible, 'conditions', cfg.conditions, ...
    'x', cfg.lle_field, 'y', 'mc_total', 'color_by', cfg.color_by, 'size_by', cfg.size_by, ...
    'clip', struct(cfg.lle_field, cfg.lle_clip), 'fig_tag', cfg.fig_tag);
end
