function cfg = mu7chain_gainbias_config(run_mode, figure_set)
% MU7CHAIN_GAINBIAS_CONFIG Block C's bundle: memory capacity over gain against the bias and against the self-coupling, on the mu7chain network (FD-010, FD-024).
%
%   cfg = mu7chain_gainbias_config()                       % 'essential', figure set '4cond'
%   cfg = mu7chain_gainbias_config('essential', '3cond')   % the paper's three conditions
%   cfg = mu7chain_gainbias_config('fast', '4cond')        % with product MTS added
%
% WHAT IT RUNS. One stage, the opt-in 'mc_gain_sweep' of
% run_all_paper_analyses (paper_opt_in_stages), on the bundle's noise-free
% memory-capacity twin, all four conditions on the same seeds. Its runs
% (cfg.mc_gain_sweep.grids; mc_gain_sweep_stage_plan), each into its own
% folder of the run directory:
%   mc_gain_sweep_bias          gain x bias (S_c): gains log-spaced 0.1 to
%                               4, biases linear -0.6 to 0.6
%   mc_gain_sweep_selfcoupling  gain x self-coupling (self_coupling, the
%                               diagonal of W in units of R_W; s < 0
%                               negative, s > 0 positive feedback): gains
%                               LINEARLY spaced 0.1 to 4, s linear -4 to 4
%                               (an odd count of levels, so s = 0 is one)
%   mc_gain_sweep_bias_EI3to2,  the two E:I variants of the twin (mean
%   mc_gain_sweep_bias_EI2to3   weights 3 : 2 and 2 : 3), each along gain
%                               only, on the bias grid's gains, at its own
%                               bias; the gain curves' extra rows
% The folders are 'mc_gain_sweep' followed by the grid's tag, then a
% variant's (mc_gain_sweep_folder), so a listing of the run directory shows
% each grid with its variants beside it. Nothing else runs: the sweeps
% stage only creates the run directory and its manifest, which record the
% other stages as skipped. For block A's stages and these in one run, see
% mu7chain_AC_config.
%
% HOW BIG AND HOW LONG. The sizes are mc_gain_sweep_grid's table, the same
% for both grids: levels per axis x networks (seeds 3001, 3002, ...).
%   essential    3 x  3 x  1     24 jobs    about 20 min on 12 workers
%   fast         7 x  7 x  3    336 jobs    about 4 h 25
%   medium(2)   13 x 13 x 10   3640 jobs    about 48 h
%   production  17 x 17 x 10   6120 jobs    about 81 h
% (2 n^2 + 2 n) x networks jobs at n levels, 9.5 worker-minutes a job, all
% of the stage's jobs in one parfor. Agents run 'essential' and 'fast';
% 'medium' and larger are Tom's to start.
%
% DERIVED FROM mu7chain_config, so nothing about the network is restated:
% the preset, the memory-capacity twin (cfg.mc_preset), its integrator
% (cfg.mc_ode_solver), the conditions drawn per set and the sets drawn per
% mode (cfg.figure_sets) are mu7chain_config's. What differs:
%   cfg.stages        {'mc_gain_sweep'}; cfg.sweeps and eig examples empty
%   cfg.mc_gain_sweep the stage's settings: the two grids, the variants of
%                     the bias grid (mc_gain_sweep_stage_plan lists every
%                     field; a grid's gain_levels or second_levels, or
%                     seeds and conditions for the stage, depart from the
%                     table)
%   cfg.run_dir       data/mu7chain_gainbias_<run_mode>
%   cfg.fig_root      figs/mu7chain_gainbias_<run_mode>_<figure_set>
%   cfg.figures       the block C figures only (below)
%
% THE FIGURES (Tom, 2026-11-13), all drawn from the saved tables, each in
% its own folder <fig_root>/<entry>/:
%   fig_mc_maps_bias, fig_mc_maps_selfcoupling
%       the maps of memory capacity and lambda_1 on each grid, LINEAR color
%       scales (memory capacity 0 to its maximum, lambda_1 clipped at +-5
%       1/s): the paper's maps. The square-root and symmetric-log scales of
%       fig_mc_gain_bias_maps' defaults are for reports and not drawn here.
%   fig_mc_gain_curves
%       memory capacity and lambda_1 against gain at the preset's bias (0.2;
%       the nearest level of the bias grid), a row for the network and one
%       per E:I variant; two columns (without memory capacity against
%       lambda_1)
%   fig_mc_robustness_bias, fig_mc_robustness_selfcoupling
%       the distribution of memory capacity over each map
%   fig_mc_vs_lambda1
%       memory capacity against lambda_1, every run of the SELF-COUPLING
%       grid, colored by gain, marker area by |s| (filled s >= 0, rings
%       s < 0, clipped runs open squares)
%   fig_lambda1_vs_bias, fig_rate_vs_lambda1, fig_slope_share_vs_gain
%       explorations on the bias grid: lambda_1 against bias; mean rate
%       against lambda_1 with marker area by bias; the share of neurons on
%       the slope (neither pinned nor silent) against gain, colored by bias
%   mc_statistics_bias, mc_statistics_selfcoupling
%       the statistics tables of each grid (mc_gain_bias_statistics)
% Each entry has its own file tag (Fig_MC_maps_bias, ...), none a prefix of
% another: save_figure_stable deletes <tag>* before it saves, so two
% figures of one function must not share the function's default tag.
%
% See also: mu7chain_gainbias_run, mu7chain_AC_config, mu7chain_config,
%           run_mc_gain_sweep, mc_gain_sweep_stage_plan, mc_gain_sweep_grid,
%           run_all_paper_analyses

arguments
    run_mode (1,:) char = 'essential'
    figure_set (1,:) char {mustBeMember(figure_set, {'3cond', '4cond'})} = '4cond'
end
if ~ismember(run_mode, run_mode_names())
    error('mu7chain_gainbias_config:BadRunMode', 'run_mode ''%s'' is not one of: %s.', ...
        run_mode, strjoin(run_mode_names(), ', '));
end

cfg = mu7chain_config(run_mode, figure_set);

%% Which stage runs: the opt-in block C stage alone
cfg.stages = {'mc_gain_sweep'};
cfg.sweeps = {};
cfg.eig_heatmap_examples = [];

%% The stage's grids
% A grid's tag names its folder and is how a figure asks for it. The
% variant networks belong to the bias grid: a variant's tag names its folder
% and the figure's row, and its preset is the memory-capacity twin's variant
% of that name.
BIAS = 'bias';
SELF = 'selfcoupling';
variant_tags = {'EI3to2', 'EI2to3'};
variants = struct('tag', variant_tags, 'preset', strcat(cfg.mc_preset, '_', variant_tags));
cfg.mc_gain_sweep = struct('grids', struct( ...
    'tag',          {BIAS,        SELF}, ...
    'second_axis',  {'S_c',       'self_coupling'}, ...
    'gain_range',   {[0.1 4],     [0.1 4]}, ...
    'gain_spacing', {'log',       'linear'}, ...
    'second_range', {[-0.6 0.6],  [-4 4]}, ...
    'variants',     {variants,    []}));

%% Where things land -- fixed, and apart from mu7chain_config's
cfg.run_dir  = ['data/mu7chain_gainbias_' run_mode];
cfg.fig_root = ['figs/mu7chain_gainbias_' run_mode '_' figure_set];

%% The figure registry: block C only
% Each entry gets run_dir, out_dir, save, visible, verbose and the set's
% conditions from make_all_paper_figures, then its own arguments. Every
% entry names its grid; the figure reads that grid's result and no other.
LINEAR_SCALES = {'mc_scale', 'linear', 'lle_scale', 'linear', 'lle_clip', 5};
% The scatters' statistics boxes are placed where their panels have fewest
% points (explore_bias_rate_lambda1, where these three were first drawn).
cfg.figures = { ...
    grid_entry('fig_mc_maps_bias',               @fig_mc_gain_bias_maps, true,  BIAS, ...
        [{'fig_tag', 'Fig_MC_maps_bias'}, LINEAR_SCALES]), ...
    grid_entry('fig_mc_maps_selfcoupling',       @fig_mc_gain_bias_maps, true,  SELF, ...
        [{'fig_tag', 'Fig_MC_maps_selfcoupling'}, LINEAR_SCALES]), ...
    grid_entry('fig_mc_gain_curves',             @fig_mc_gain_curves,    true,  BIAS, ...
        {'fig_tag', 'Fig_MC_gain_curves', 'variants', variant_tags, 'mc_vs_lambda1', false}), ...
    grid_entry('fig_mc_robustness_bias',         @fig_mc_robustness,     true,  BIAS, ...
        {'fig_tag', 'Fig_MC_robustness_bias'}), ...
    grid_entry('fig_mc_robustness_selfcoupling', @fig_mc_robustness,     true,  SELF, ...
        {'fig_tag', 'Fig_MC_robustness_selfcoupling'}), ...
    grid_entry('fig_mc_vs_lambda1',              @fig_mc_vs_lambda1,     true,  SELF, ...
        {'fig_tag', 'Fig_MC_vs_lambda1', 'size_by', 'self_coupling'}), ...
    grid_entry('fig_lambda1_vs_bias',            @fig_mc_sweep_scatter,  false, BIAS, ...
        {'fig_tag', 'Fig_lambda1_vs_bias', 'x', 'bias', 'y', 'lle_input', 'color_by', 'gain', ...
         'stats_position', [0.56 0.27]}), ...
    grid_entry('fig_rate_vs_lambda1',            @fig_mc_sweep_scatter,  false, BIAS, ...
        {'fig_tag', 'Fig_rate_vs_lambda1', 'x', 'lle_input', 'y', 'rate_mean', 'color_by', 'gain', ...
         'size_by', 'bias', 'stats_position', [0.56 0.99]}), ...
    grid_entry('fig_slope_share_vs_gain',        @fig_mc_sweep_scatter,  false, BIAS, ...
        {'fig_tag', 'Fig_slope_share_vs_gain', 'x', 'gain', 'y', 'frac_on_slope', 'color_by', 'bias', ...
         'stats_position', [0.50 0.22]}), ...
    grid_entry('mc_statistics_bias',             @statistics_entry,      false, BIAS, ...
        {'tag', 'MC_statistics_bias'}), ...
    grid_entry('mc_statistics_selfcoupling',     @statistics_entry,      false, SELF, ...
        {'tag', 'MC_statistics_selfcoupling'})};
end

%% ------------------------------------------------------------------------
function entry = grid_entry(name, fn, in_paper, grid, args)
% A registry entry that reads one grid of the run directory.
entry = struct('name', name, 'fn', fn, 'in_paper', in_paper, 'args', {[{'grid', grid}, args]});
end

function out = statistics_entry(opts)
% mc_gain_bias_statistics behind the figure contract: one grid's result
% from the run directory, the set's conditions, the tables written to
% out_dir. out.files lists them; there is no figure.
arguments
    opts.run_dir     (1,:) char    = ''
    opts.data_file   (1,:) char    = ''
    opts.grid        (1,:) char    = ''
    opts.out_dir     (1,:) char    = ''
    opts.save        (1,1) logical = true
    opts.visible     (1,1) logical = false     % accepted for a uniform call
    opts.verbose     (1,:) char    = 'minimal'
    opts.conditions                = {}
    opts.preset_name (1,:) char    = ''        % accepted for a uniform call
    opts.reference   (1,:) char    = 'sfa3_stdc2'
    opts.tag         (1,:) char    = 'MC_gain_bias_statistics'
end
D = resolve_mc_gain_sweep('data_file', opts.data_file, 'run_dir', opts.run_dir, 'grid', opts.grid);
out_dir = '';
if opts.save; out_dir = opts.out_dir; end
[~, files] = mc_gain_bias_statistics('data', D, 'conditions', opts.conditions, ...
    'reference', opts.reference, 'out_dir', out_dir, 'tag', opts.tag, 'verbose', opts.verbose);
out = struct('figs', gobjects(0), 'files', {files}, 'source', D.source);
end
