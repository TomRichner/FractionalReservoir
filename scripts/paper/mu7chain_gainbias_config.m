function cfg = mu7chain_gainbias_config(run_mode, figure_set, run_grids)
% MU7CHAIN_GAINBIAS_CONFIG Block C's bundle: memory capacity over gain against the bias, the self-coupling and two kinds of noise, on the mu7chain network (FD-010, FD-024, FD-027).
%
%   cfg = mu7chain_gainbias_config()                       % 'essential', figure set '4cond'
%   cfg = mu7chain_gainbias_config('essential', '3cond')   % the paper's three conditions
%   cfg = mu7chain_gainbias_config('fast', '4cond')        % with product MTS added
%   cfg = mu7chain_gainbias_config('fast', '3cond', {'selfcoupling', 'bias'})
%   cfg = mu7chain_gainbias_config('fast', '3cond', 'all') % every grid defined, as before FD-029
%
% run_grids, optional: the tags of the grids that run, and whose figures
% are registered; absent or empty = the bundle's choice, RUN_GRIDS below
% (the self-coupling grid alone). 'all' = every grid defined (for example
% to draw the figures of a run directory made before FD-029).
%
% WHAT IT RUNS. One stage, the opt-in 'mc_gain_sweep' of
% run_all_paper_analyses (paper_opt_in_stages), on the bundle's noise-free
% memory-capacity twin, all four conditions on the same seeds. Its grids
% (cfg.mc_gain_sweep.grids; mc_gain_sweep_stage_plan) are listed below,
% each with its folder in the run directory. SINCE FD-029 (Tom,
% 2026-11-16) THREE RUN, in every mode: the gain x self-coupling grid and
% the two 1-D noise sweeps at the preset's gain. cfg.mc_gain_sweep.run_grids
% names the grids that run (RUN_GRIDS below), the others stay defined and
% selectable, and the figure registry keeps only the entries whose grids
% run. One tag added to RUN_GRIDS brings a grid back, with its figures (the
% bias grid with its two variants).
%   mc_gain_sweep_selfcoupling  gain x self-coupling (self_coupling, the
%                               diagonal of W in units of R_W; s < 0
%                               negative, s > 0 positive feedback): gains
%                               LOG-spaced 0.05 to 5, s linear -4 to 6
%                               (FD-029; until then gains linear 0.1 to 4
%                               and s -4 to 4). s = 0 is a level at 11 and
%                               21 levels (fastmedium, medium), not at 7
%                               (essential, fast) nor 17 (production).
%   mc_gain_sweep_noise1d_wiener   the Wiener noise sigma_u alone (0, then
%                               log-spaced 2.5e-4 to 0.025), release off,
%                               at the preset's gain only
%   mc_gain_sweep_noise1d_release  the release noise sigma_gamma alone (0,
%                               then log-spaced 1e-4 to 1e-2) at the leak
%                               r_s = 0.002, sigma_u = 0, at the preset's
%                               gain only
% Defined, not run (add to RUN_GRIDS to run):
%   mc_gain_sweep_bias          gain x bias (S_c): gains log-spaced 0.1 to
%                               4, biases linear -0.6 to 0.6
%   mc_gain_sweep_bias_EI3to2,  the two E:I variants of the twin (mean
%   mc_gain_sweep_bias_EI2to3   weights 3 : 2 and 2 : 3), each along gain
%                               only, on the bias grid's gains, at its own
%                               bias; the gain curves' extra rows
%   mc_gain_sweep_releasenoise  gain (log-spaced 0.1 to 4) x the release
%                               noise sigma_gamma (release_sigma_all; 0,
%                               then log-spaced to 0.01) at the leak
%                               r_s = 0.002, Wiener off (FD-027)
%   mc_gain_sweep_selfcoupling_releasenoise   self-coupling (linear -4 to
%                               4) x sigma_gamma, at the preset's gain, the
%                               same leak, Wiener off
%   mc_gain_sweep_wienernoise   gain x the Wiener noise sigma_u
%                               (sigma_u_noise; 0, then log-spaced to
%                               0.025), release off
% The noise grids, 2-D and 1-D, read the memory through the tau_d-filtered
% release (readout_signal 'release_filtered'); the others keep the synaptic
% readout. Their fixed arguments (model_args per grid) and the calibration
% of the noise axes are stated where the grids are written, below.
% The folders are 'mc_gain_sweep' followed by the grid's tag, then a
% variant's (mc_gain_sweep_folder), so a listing of the run directory shows
% each grid with its variants beside it. Nothing else runs: the sweeps
% stage only creates the run directory and its manifest, which record the
% other stages as skipped. For block A's stages and these in one run, see
% mu7chain_AC_config.
%
% HOW BIG AND HOW LONG. The sizes are mc_gain_sweep_grid's table: levels
% per axis x networks (seeds 3001, 3002, ...); the 1-D noise sweeps take
% the table's count of noise levels except 5 at 'essential'. Since FD-029:
%   mode        self-coupling   each 1-D sweep   jobs (release)   on 12 workers
%   essential    7 x  7 x 2        5 x 2          118  (10)       about 1 h 44
%   fast         7 x  7 x 3        7 x 3          189  (21)       about 2 h 48
%   fastmedium  11 x 11 x 5       11 x 5          715  (55)       about 10 h 27
%   medium(2)   21 x 21 x 5       21 x 5         2415 (105)       about 34 h 52
%   production  17 x 17 x 10      17 x 10        3230 (170)       about 46 h 47
% With all five 2-D grids and the two variants, as until FD-029, it was
% (5 n^2 + 2 n) x networks jobs: 51 at essential (58.9 min measured), 777
% at fast (about 12 h 30), and at the sizes of then 846 at fastmedium
% (9 x 9 x 2, about 13 h 40; the run of 2026-11-15 took 19 h 28 with block
% A) and 8710 at medium (13 x 13 x 10).
% All of the stage's jobs run in one
% parfor. A job of four conditions is about 10.24 worker-minutes (the fast
% run of 2026-11-13: 336 jobs without release noise in 286.7 min on 12
% workers), and 1.35 times that with the release noise on (its
% condition-runs measured 157 s against 116 s, FD-027 N2; the calibration
% runs of 2026-11-14 took 170 to 200 s with release noise against 124 to
% 130 s without, the Wiener noise costing nothing). The stage's own
% estimate (plan.worker_minutes) counts both (mc_gain_sweep_grid). Agents
% run 'essential' and 'fast'; 'fastmedium' and larger are Tom's to start.
%
% DERIVED FROM mu7chain_config, so nothing about the network is restated:
% the preset, the memory-capacity twin (cfg.mc_preset), its integrator
% (cfg.mc_ode_solver), the conditions drawn per set and the sets drawn per
% mode (cfg.figure_sets) are mu7chain_config's. What differs:
%   cfg.stages        {'mc_gain_sweep'}; cfg.sweeps and eig examples empty
%   cfg.mc_gain_sweep the stage's settings: the seven grids, which of them
%                     run (run_grids), the variants of
%                     the bias grid (mc_gain_sweep_stage_plan lists every
%                     field; a grid's gain_levels or second_levels, or
%                     seeds and conditions for the stage, depart from the
%                     table)
%   cfg.run_dir       data/mu7chain_gainbias_<run_mode>
%   cfg.fig_root      figs/mu7chain_gainbias_<run_mode>_<figure_set>
%   cfg.figures       the block C figures only (below)
%
% THE FIGURES (Tom, 2026-11-13), all drawn from the saved tables, each in
% its own folder <fig_root>/<entry>/. Only the entries whose grids run are
% registered (since FD-029: fig_mc_maps_selfcoupling,
% fig_mc_robustness_selfcoupling, fig_mc_vs_lambda1,
% mc_statistics_selfcoupling, fig_mc_noise_sweeps, and the two time
% series, which read no grid); the others below come back with their grid.
%   fig_mc_maps_bias, fig_mc_maps_selfcoupling
%       the maps of memory capacity and lambda_1 on each grid: memory
%       capacity on a LINEAR scale (0 to its maximum), lambda_1 on the
%       symmetric-log scale (linear within +-0.05 1/s, clipped at +-5 1/s;
%       Tom, 2026-11-16, FD-029; linear until then). The paper's maps.
%       Since FD-029 in Brian's style (MAP_STYLE below): ticks on the bottom
%       and left only, no lambda_1 = 0 line, no saturated or silent marks,
%       no spread row, lambda_1 on CMasher's guppy (dark center) and memory
%       capacity from white to its center color.
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
%   fig_mc_noise_sweeps
%       the two 1-D noise sweeps at the preset's gain (FD-029): memory
%       capacity (and lambda_1) against sigma_u and against sigma_gamma,
%       one line per condition, with the spread over networks
%       (fig_mc_noise_sweeps), and its table
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
    run_grids = {}
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
NOISE_WIENER  = 'noise1d_wiener';    % the 1-D noise sweeps at the preset's gain (FD-029)
NOISE_RELEASE = 'noise1d_release';
% WHICH GRIDS RUN (FD-029, Tom, 2026-11-16: the paper shows the gain x
% self-coupling grid, and the effect of each noise at the preset, so the
% others are not computed). Every grid below stays defined; this list says
% which run, and the figure registry keeps only the entries whose grids all
% run (and those that read no grid). To bring a grid back, add its tag
% here, e.g. {SELF, NOISE_WIENER, NOISE_RELEASE, BIAS}.
RUN_GRIDS = {SELF, NOISE_WIENER, NOISE_RELEASE};
if ~isempty(run_grids); RUN_GRIDS = run_grids; end
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
%   sigma_gamma  1e-4 to 1e-2:   memory capacity within 10 % of the leak's
%                below about 2.4e-4 (2.15e-4 is the lowest positive level at
%                'fast'); 1.56 at the top, the floor the brief asked for
%   sigma_u      2.5e-4 to 0.025: within 10 % of no noise below about 6e-4
%                (5.4e-4 is the lowest positive level at 'fast'); 0.8 at the
%                top, the paper's noise level
% WHY THE sigma_gamma TOP IS 1e-2 AND NOT HIGHER (lead, after review): at the
% leak r_s the mean size of a release event is sigma_gamma^2 / r_s. At 2e-2
% (first proposed, to bring memory capacity nearer 1) that is 0.2 against
% tau_rel = 0.25: one event empties about a whole pool (an equivalent pool of
% about one vesicle, events about once per 100 s per neuron), which is not
% spontaneous release. At 1e-2 it is 0.05, a pool of about five vesicles,
% and memory capacity is already at its floor.
% Levels by mode (rounded to 10 decimals by the grid table):
%   essential, 1-D sweeps (5 levels)
%              sigma_gamma 0, 3.16e-4 ... 1e-2 (x3.16) sigma_u 0, 7.91e-4 ... 0.025 (x3.16)
%   essential (2-D grids) and fast, 7 levels
%              sigma_gamma 0, 2.15e-4 ... 1e-2 (x2.15) sigma_u 0, 5.39e-4 ... 0.025 (x2.15)
%   fastmedium sigma_gamma 0, 1.58e-4 ... 1e-2 (x1.58) sigma_u 0, 3.96e-4 ... 0.025 (x1.58)
%   medium     sigma_gamma 0, 1.26e-4 ... 1e-2 (x1.26) sigma_u 0, 3.15e-4 ... 0.025 (x1.26)
%   (at the sizes of FD-029; the 2-D noise grids do not run unless named in
%   RUN_GRIDS, the 1-D sweeps do)
LEAK = 0.002;
GAMMA_LADDER = [1e-4 1e-2];
WIENER_LADDER = [2.5e-4 0.025];
release_fixed = {'sigma_u_noise', 0, 'release_rate_all', LEAK};   % the gamma grids: Wiener off, the leak
wiener_fixed  = {'release_rate_all', 0, 'release_sigma_all', 0};  % the Wiener grid: release off
FILTERED = 'release_filtered';
% THE 1-D NOISE SWEEPS (FD-029, Tom, 2026-11-16: "a figure at just the
% preset but with each type of noise varied ... two 1D sweeps, which would
% be much faster and still tell important info about the effect of noise on
% MC"), in place of the 2-D noise grids, which stay defined but do not run.
% Each is a grid with ONE level on its first axis, the preset's gain, and
% the noise on its second, with the 2-D grids' ladders, fixed arguments and
% readout: sigma_u with the release off, and sigma_gamma at the leak r_s
% with sigma_u = 0. Noise levels (with the zero) x networks: the table's
% n x networks (7 x 3, 11 x 5, 21 x 5) except 5 x 2 at 'essential'.
preset_gain = mc_gain_sweep_preset_point(cfg.mc_preset);
noise_1d_count = [];
if strcmp(run_mode, 'essential'); noise_1d_count = 5; end
cfg.mc_gain_sweep = struct('grids', [ ...
    grid_spec(BIAS, 'second_axis', 'S_c', 'gain_range', [0.1 4], 'gain_spacing', 'log', ...
        'second_range', [-0.6 0.6], 'variants', variants), ...
    grid_spec(SELF, 'second_axis', 'self_coupling', 'gain_range', [0.05 5], 'gain_spacing', 'log', ...
        'second_range', [-4 6]), ...
    grid_spec(GAMMA, 'second_axis', 'release_sigma_all', 'gain_range', [0.1 4], 'gain_spacing', 'log', ...
        'second_range', GAMMA_LADDER, 'second_spacing', 'zero_log', ...
        'model_args', release_fixed, 'readout_signal', FILTERED), ...
    grid_spec(SELF_GAMMA, 'first_axis', 'self_coupling', 'first_range', [-4 4], 'first_spacing', 'linear', ...
        'second_axis', 'release_sigma_all', 'second_range', GAMMA_LADDER, 'second_spacing', 'zero_log', ...
        'model_args', release_fixed, 'readout_signal', FILTERED), ...
    grid_spec(WIENER, 'second_axis', 'sigma_u_noise', 'gain_range', [0.1 4], 'gain_spacing', 'log', ...
        'second_range', WIENER_LADDER, 'second_spacing', 'zero_log', ...
        'model_args', wiener_fixed, 'readout_signal', FILTERED), ...
    grid_spec(NOISE_WIENER, 'first_levels', preset_gain, 'second_axis', 'sigma_u_noise', ...
        'second_range', WIENER_LADDER, 'second_spacing', 'zero_log', 'second_count', noise_1d_count, ...
        'model_args', wiener_fixed, 'readout_signal', FILTERED), ...
    grid_spec(NOISE_RELEASE, 'first_levels', preset_gain, 'second_axis', 'release_sigma_all', ...
        'second_range', GAMMA_LADDER, 'second_spacing', 'zero_log', 'second_count', noise_1d_count, ...
        'model_args', release_fixed, 'readout_signal', FILTERED)]);
if isequal(RUN_GRIDS, 'all'); RUN_GRIDS = {cfg.mc_gain_sweep.grids.tag}; end
cfg.mc_gain_sweep.run_grids = cellstr(RUN_GRIDS);
% The time-series figures draw zero, the middle and the top of each noise
% axis: a three-level zero_log ladder over its range, from the grid table.
three = @(axis, ladder) getfield(mc_gain_sweep_grid('essential', 'second_axis', axis, ...
    'second_range', ladder, 'second_spacing', 'zero_log', 'second_count', 3, 'gain_levels', 1), ...
    'second_levels');
gamma_three  = three('release_sigma_all', GAMMA_LADDER);
wiener_three = three('sigma_u_noise', WIENER_LADDER);

%% Where things land -- fixed, and apart from mu7chain_config's
cfg.run_dir  = ['data/mu7chain_gainbias_' run_mode];
cfg.fig_root = ['figs/mu7chain_gainbias_' run_mode '_' figure_set];

%% The figure registry: block C only
% Each entry gets run_dir, out_dir, save, visible, verbose and the set's
% conditions from make_all_paper_figures, then its own arguments. Every
% entry names its grid; the figure reads that grid's result and no other.
% The maps' color scales: memory capacity linear from 0 to its maximum;
% lambda_1 on the symmetric-log scale (symlog_lle: linear within +-0.05 1/s,
% logarithmic beyond, clipped at +-5 1/s).
MAP_SCALES = {'mc_scale', 'linear', 'lle_scale', 'symlog', 'lle_lin', 0.05, 'lle_clip', 5};
% The maps' style (FD-029, from Brian Lundstrom's feedback): no ticks on the
% top and right, no lambda_1 = 0 line, no saturated or silent markers, no
% spread row; lambda_1 on a diverging map with a dark center (lambda_1 = 0
% the darkest band) and memory capacity from white to that center color.
% THE COLORMAP AND THE lambda_1 SCALE ARE TOM'S CHOICE among the FD-029
% candidates (2026-11-16: "guppy with symlog looks best to me"): CMasher's
% guppy (dark_center_colormaps), on the symmetric-log scale. Symlog because
% on the linear +-5 scale the stable cells of the adapting conditions, near
% lambda_1 = -0.1, are as dark as zero itself, so their sign cannot be seen;
% symlog spreads the band around zero so they read blue.
PAPER_LLE_COLORMAP = 'guppy';
MAP_STYLE = [MAP_SCALES, {'ticks_top_right', false, 'show_zero_line', false, 'show_flags', false, ...
    'show_iqr', false, 'lle_colormap', PAPER_LLE_COLORMAP, 'mc_colormap', 'lle_center'}];
% The scatters' statistics boxes are placed where their panels have fewest
% points (explore_bias_rate_lambda1, where these three were first drawn).
cfg.figures = { ...
    grid_entry('fig_mc_maps_bias',               @fig_mc_gain_bias_maps, true,  BIAS, ...
        [{'fig_tag', 'Fig_MC_maps_bias'}, MAP_STYLE]), ...
    grid_entry('fig_mc_maps_selfcoupling',       @fig_mc_gain_bias_maps, true,  SELF, ...
        [{'fig_tag', 'Fig_MC_maps_selfcoupling'}, MAP_STYLE]), ...
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
        [{'fig_tag', 'Fig_MC_maps_releasenoise_gain'}, MAP_STYLE]), ...
    grid_entry('fig_mc_maps_selfcoupling_releasenoise',  @fig_mc_gain_bias_maps, false, SELF_GAMMA, ...
        [{'fig_tag', 'Fig_MC_maps_releasenoise_selfcoupling'}, MAP_STYLE]), ...
    grid_entry('fig_mc_maps_wienernoise',                @fig_mc_gain_bias_maps, false, WIENER, ...
        [{'fig_tag', 'Fig_MC_maps_wienernoise_gain'}, MAP_STYLE]), ...
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
    ...% The two 1-D noise sweeps at the preset's gain, one figure (FD-029).
    struct('name', 'fig_mc_noise_sweeps', 'fn', @fig_mc_noise_sweeps, 'in_paper', true, ...
        'args', {{'grids', {NOISE_WIENER, NOISE_RELEASE}, 'fig_tag', 'Fig_MC_noise_sweeps'}}), ...
    ...% The two time-series figures: they simulate, and read no run data.
    struct('name', 'fig_noise_levels_wiener', 'fn', @fig_noise_levels_timeseries, 'in_paper', false, ...
        'args', {{'preset_name', cfg.preset_name, 'level_name', 'sigma_u_noise', ...
        'levels', wiener_three, 'level_label', '\sigma_u', 'model_args', wiener_fixed, ...
        'tag', 'Fig_noise_levels_wiener'}}), ...
    struct('name', 'fig_noise_levels_gamma', 'fn', @fig_noise_levels_timeseries, 'in_paper', false, ...
        'args', {{'preset_name', cfg.preset_name, 'level_name', 'release_sigma_all', ...
        'levels', gamma_three, 'level_label', '\sigma_\gamma', 'model_args', release_fixed, ...
        'tag', 'Fig_noise_levels_gamma'}})};
cfg.figures = entries_of_grids(cfg.figures, RUN_GRIDS);
end

%% ------------------------------------------------------------------------
function figures = entries_of_grids(figures, run_grids)
% The registry entries whose grids all run, and those that read no grid
% (they simulate, or need no block C data): an entry of a grid that is not
% computed would only fail on a new run directory.
keep = cellfun(@(f) all(ismember(entry_grids(f.args), run_grids)), figures);
figures = figures(keep);
end

function g = entry_grids(args)
% The grids an entry reads, by its 'grid' (one) or 'grids' (several)
% argument; empty when it names none.
g = {};
names = args(1:2:end);
for key = {'grid', 'grids'}
    k = find(strcmp(names, key{1}), 1);
    if ~isempty(k); g = [g, cellstr(args{2 * k})]; end %#ok<AGROW>
end
end

%% ------------------------------------------------------------------------
function g = grid_spec(tag, varargin)
% One grid of the stage, every field of the grid struct present (empty
% where the grid leaves the default), so the grids form one struct array.
g = struct('tag', tag, 'first_axis', [], 'first_range', [], 'first_spacing', [], 'first_levels', [], ...
    'second_axis', [], 'gain_range', [], 'gain_spacing', [], 'second_range', [], ...
    'second_spacing', [], 'second_count', [], 'model_args', [], 'readout_signal', [], 'variants', []);
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
