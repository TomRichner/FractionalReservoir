function cfg = mu7chain_gainbias_config(run_mode, figure_set)
% MU7CHAIN_GAINBIAS_CONFIG Block C's bundle: memory capacity over gain against the bias, the self-coupling and two kinds of noise, on the mu7chain network (FD-010, FD-024, FD-027).
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
%   mc_gain_sweep_releasenoise  gain (log-spaced 0.1 to 4) x the release
%                               noise sigma_gamma (release_sigma_all; 0,
%                               then log-spaced to 0.02) at the leak
%                               r_s = 0.002, Wiener off (FD-027)
%   mc_gain_sweep_selfcoupling_releasenoise   self-coupling (linear -4 to
%                               4) x sigma_gamma, at the preset's gain, the
%                               same leak, Wiener off
%   mc_gain_sweep_wienernoise   gain x the Wiener noise sigma_u
%                               (sigma_u_noise; 0, then log-spaced to
%                               0.025), release off
% The three noise grids read the memory through the tau_d-filtered release
% (readout_signal 'release_filtered'); the first four keep the synaptic
% readout. Their fixed arguments (model_args per grid) and the calibration
% of the noise axes are stated where the grids are written, below.
% The folders are 'mc_gain_sweep' followed by the grid's tag, then a
% variant's (mc_gain_sweep_folder), so a listing of the run directory shows
% each grid with its variants beside it. Nothing else runs: the sweeps
% stage only creates the run directory and its manifest, which record the
% other stages as skipped. For block A's stages and these in one run, see
% mu7chain_AC_config.
%
% HOW BIG AND HOW LONG. The sizes are mc_gain_sweep_grid's table, the same
% for every grid: levels per axis x networks (seeds 3001, 3002, ...).
%   mode         levels    jobs   of which release noise   on 12 workers
%   essential    3 x  3 x  1    51      18                about 45 min
%   fast         7 x  7 x  3   777     294                about 11 h 40
%   medium(2)   13 x 13 x 10  8710    3380                about 131 h
%   production  17 x 17 x 10 14790    5780                about 223 h
% (5 n^2 + 2 n) x networks jobs at n levels, all of the stage's jobs in one
% parfor. A job of four conditions is about 9.5 worker-minutes, and about
% 12.9 with the release noise on (its condition-runs measured 157 s against
% 116 s, FD-027 N2; the calibration runs of 2026-11-14 took 170 to 200 s
% with release noise against 124 to 130 s without, the Wiener noise costing
% nothing). The stage's own estimate (plan.worker_minutes) counts 9.5 for
% every job and so falls short by about a quarter. Without the noise grids
% (the first four alone) the stage was 24 jobs, about 20 min, at
% 'essential' and 336, about 4 h 25, at 'fast'. Agents run 'essential' and
% 'fast'; 'medium' and larger are Tom's to start.
%
% DERIVED FROM mu7chain_config, so nothing about the network is restated:
% the preset, the memory-capacity twin (cfg.mc_preset), its integrator
% (cfg.mc_ode_solver), the conditions drawn per set and the sets drawn per
% mode (cfg.figure_sets) are mu7chain_config's. What differs:
%   cfg.stages        {'mc_gain_sweep'}; cfg.sweeps and eig examples empty
%   cfg.mc_gain_sweep the stage's settings: the five grids, the variants of
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
%   fig_mc_maps_<noise grid>, fig_mc_robustness_<noise grid>,
%   mc_statistics_<noise grid>, for releasenoise,
%   selfcoupling_releasenoise and wienernoise (FD-027)
%       the same three for each noise grid, linear color scales; the noise
%       axis drawn as its log ladder, the zero one step below the lowest
%       positive level and the axis marked broken there (mc_axis_spacing).
%       Tags name the noise, then the first axis (Fig_MC_maps_releasenoise_gain,
%       ..._releasenoise_selfcoupling, ..._wienernoise_gain). Widths along
%       the noise axis in the statistics are in sigma's own units. Not
%       marked in_paper (Tom decides after the first maps).
%   fig_noise_levels_wiener, fig_noise_levels_gamma
%       fig_noise_levels_timeseries: the bundle's noisy network, every
%       condition of the set, at zero, the middle and the top of each noise
%       axis (its 'essential' levels): sigma_u with the release off, and
%       sigma_gamma at the leak with sigma_u_noise = 0. They simulate (about
%       a minute each) and read no run data.
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
GAMMA = 'releasenoise';
SELF_GAMMA = 'selfcoupling_releasenoise';
WIENER = 'wienernoise';
variant_tags = {'EI3to2', 'EI2to3'};
variants = struct('tag', variant_tags, 'preset', strcat(cfg.mc_preset, '_', variant_tags));
% The noise grids (FD-027, step P3). Each noise axis is a zero level, then
% log-spaced levels ('zero_log', mc_gain_sweep_grid: the mode's n-level log
% ladder over the range with its lowest rung replaced by 0), and all three
% read the memory through the tau_d-filtered release ('release_filtered'),
% the Wiener grid included, so the three are comparable with each other.
% The ranges were calibrated on 2026-11-14 with eight single lean memory-
% capacity runs: chained MTS (sfa3_stdc2), gain 1, seed 3001, the
% filtered readout, the full mc_protocol (lambda_1 over 6 s only), on the
% noise-free twin with the noise in question switched on:
%   no noise (the lead's measurement)            MC 11.2
%   leak alone, r_s = 0.002                      MC 12.50   (at r_s = 0.01: 16.1)
%   r_s = 0.002, sigma_gamma = 1e-4 / 5e-4 / 2e-3 / 1e-2   MC 12.24 / 10.34 / 5.95 / 1.56
%   Wiener alone, sigma_u = 5e-4 / 2.5e-3 / 1e-2          MC 10.32 / 7.44 / 2.76
%   (sigma_u = 0.025, the lead's measurement)              MC 0.8
% So the leak is r_s = 0.002 (it still moves memory capacity by +12 %; the
% 0.01 of the plan moved it by +44 %), and the ladders are
%   sigma_gamma  1e-4 to 2e-2:   memory capacity within 10 % of the leak's
%                below about 2.4e-4, the lowest positive level at 'fast';
%                about 1 near the top (1.56 at 1e-2)
%   sigma_u      2.5e-4 to 0.025: within 10 % of no noise below about 6e-4
%                (5.4e-4 is the lowest positive level at 'fast'); 0.8 at the
%                top, the paper's noise level
% Levels by mode (rounded to 10 decimals by the grid table):
%   essential  sigma_gamma 0, 1.41e-3, 0.02          sigma_u 0, 2.5e-3, 0.025
%   fast       sigma_gamma 0, 2.42e-4 ... 0.02 (x2.42) sigma_u 0, 5.39e-4 ... 0.025 (x2.15)
%   medium     sigma_gamma 0, 1.56e-4 ... 0.02 (x1.56) sigma_u 0, 3.67e-4 ... 0.025 (x1.47)
LEAK = 0.002;
GAMMA_LADDER = [1e-4 2e-2];
WIENER_LADDER = [2.5e-4 0.025];
release_fixed = {'sigma_u_noise', 0, 'release_rate_all', LEAK};   % the gamma grids: Wiener off, the leak
wiener_fixed  = {'release_rate_all', 0, 'release_sigma_all', 0};  % the Wiener grid: release off
FILTERED = 'release_filtered';
cfg.mc_gain_sweep = struct('grids', [ ...
    grid_spec(BIAS, 'second_axis', 'S_c', 'gain_range', [0.1 4], 'gain_spacing', 'log', ...
        'second_range', [-0.6 0.6], 'variants', variants), ...
    grid_spec(SELF, 'second_axis', 'self_coupling', 'gain_range', [0.1 4], 'gain_spacing', 'linear', ...
        'second_range', [-4 4]), ...
    grid_spec(GAMMA, 'second_axis', 'release_sigma_all', 'gain_range', [0.1 4], 'gain_spacing', 'log', ...
        'second_range', GAMMA_LADDER, 'second_spacing', 'zero_log', ...
        'model_args', release_fixed, 'readout_signal', FILTERED), ...
    grid_spec(SELF_GAMMA, 'first_axis', 'self_coupling', 'first_range', [-4 4], 'first_spacing', 'linear', ...
        'second_axis', 'release_sigma_all', 'second_range', GAMMA_LADDER, 'second_spacing', 'zero_log', ...
        'model_args', release_fixed, 'readout_signal', FILTERED), ...
    grid_spec(WIENER, 'second_axis', 'sigma_u_noise', 'gain_range', [0.1 4], 'gain_spacing', 'log', ...
        'second_range', WIENER_LADDER, 'second_spacing', 'zero_log', ...
        'model_args', wiener_fixed, 'readout_signal', FILTERED)]);
% The time-series figures draw zero, the middle and the top of each noise
% axis: its 'essential' levels, from the grid table.
three = @(axis, ladder) getfield(mc_gain_sweep_grid('essential', 'second_axis', axis, ...
    'second_range', ladder, 'second_spacing', 'zero_log', 'gain_levels', 1), 'second_levels');
gamma_three  = three('release_sigma_all', GAMMA_LADDER);
wiener_three = three('sigma_u_noise', WIENER_LADDER);

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
        {'tag', 'MC_statistics_selfcoupling'}), ...
    ...% The noise grids: maps, robustness and statistics of each. Their
    ...% tags name the noise first, then the first axis, so none is a
    ...% prefix of another or of the tags above.
    grid_entry('fig_mc_maps_releasenoise',               @fig_mc_gain_bias_maps, false, GAMMA, ...
        [{'fig_tag', 'Fig_MC_maps_releasenoise_gain'}, LINEAR_SCALES]), ...
    grid_entry('fig_mc_maps_selfcoupling_releasenoise',  @fig_mc_gain_bias_maps, false, SELF_GAMMA, ...
        [{'fig_tag', 'Fig_MC_maps_releasenoise_selfcoupling'}, LINEAR_SCALES]), ...
    grid_entry('fig_mc_maps_wienernoise',                @fig_mc_gain_bias_maps, false, WIENER, ...
        [{'fig_tag', 'Fig_MC_maps_wienernoise_gain'}, LINEAR_SCALES]), ...
    grid_entry('fig_mc_robustness_releasenoise',         @fig_mc_robustness,     false, GAMMA, ...
        {'fig_tag', 'Fig_MC_robustness_releasenoise_gain'}), ...
    grid_entry('fig_mc_robustness_selfcoupling_releasenoise', @fig_mc_robustness, false, SELF_GAMMA, ...
        {'fig_tag', 'Fig_MC_robustness_releasenoise_selfcoupling'}), ...
    grid_entry('fig_mc_robustness_wienernoise',          @fig_mc_robustness,     false, WIENER, ...
        {'fig_tag', 'Fig_MC_robustness_wienernoise_gain'}), ...
    grid_entry('mc_statistics_releasenoise',             @statistics_entry,      false, GAMMA, ...
        {'tag', 'MC_statistics_releasenoise_gain'}), ...
    grid_entry('mc_statistics_selfcoupling_releasenoise', @statistics_entry,     false, SELF_GAMMA, ...
        {'tag', 'MC_statistics_releasenoise_selfcoupling'}), ...
    grid_entry('mc_statistics_wienernoise',              @statistics_entry,      false, WIENER, ...
        {'tag', 'MC_statistics_wienernoise_gain'}), ...
    ...% The two time-series figures: they simulate, and read no run data.
    struct('name', 'fig_noise_levels_wiener', 'fn', @fig_noise_levels_timeseries, 'in_paper', false, ...
        'args', {{'preset_name', cfg.preset_name, 'level_name', 'sigma_u_noise', ...
        'levels', wiener_three, 'level_label', '\sigma_u', 'model_args', wiener_fixed, ...
        'tag', 'Fig_noise_levels_wiener'}}), ...
    struct('name', 'fig_noise_levels_gamma', 'fn', @fig_noise_levels_timeseries, 'in_paper', false, ...
        'args', {{'preset_name', cfg.preset_name, 'level_name', 'release_sigma_all', ...
        'levels', gamma_three, 'level_label', '\sigma_\gamma', 'model_args', release_fixed, ...
        'tag', 'Fig_noise_levels_gamma'}})};
end

%% ------------------------------------------------------------------------
function g = grid_spec(tag, varargin)
% One grid of the stage, every field of the grid struct present (empty
% where the grid leaves the default), so the grids form one struct array.
g = struct('tag', tag, 'first_axis', [], 'first_range', [], 'first_spacing', [], ...
    'second_axis', [], 'gain_range', [], 'gain_spacing', [], 'second_range', [], ...
    'second_spacing', [], 'model_args', [], 'readout_signal', [], 'variants', []);
for k = 1:2:numel(varargin)
    g.(varargin{k}) = varargin{k + 1};
end
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
