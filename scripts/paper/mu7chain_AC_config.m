function cfg = mu7chain_AC_config(run_mode, figure_set)
% MU7CHAIN_AC_CONFIG Blocks A and C as one bundle: mu7chain_config's stages and figures, then block C's memory-capacity grids and theirs.
%
%   cfg = mu7chain_AC_config()                       % 'essential', its set '4cond'
%   cfg = mu7chain_AC_config('fast')                 % three conditions, set '3cond' (FD-030)
%
% CONDITIONS (FD-030): mu7chain_config's, by mode: four at 'essential', three
% from 'fast' upward (product MTS not simulated), with the matching figure
% set. Every stage of both blocks and the noise switch simulates them
% (cfg.run_conditions, carried by run_all_paper_analyses).
%
% ONE RUN DIRECTORY AND ONE FIGURE ROOT PER SET hold both blocks, so one
% command per run mode (mu7chain_AC_run) produces block A's figure set and
% block C's (Tom, 2026-11-13). Block B is not part of it.
%
% IT IS A COMPOSITION, and states nothing of its own about either block:
%   mu7chain_config(run_mode, figure_set)            block A: the network,
%       its stages (all of them, or 'essential''s four), sweeps, presets,
%       conditions per set and figure registry
%   mu7chain_gainbias_config(run_mode, figure_set)   block C: the opt-in
%       stage 'mc_gain_sweep' with its seven grids defined (bias,
%       self-coupling, the three 2-D noise grids of FD-027 and the two 1-D
%       noise sweeps of FD-029; the bias grid with two variant networks;
%       cfg.mc_gain_sweep), of which since FD-029 THE GAIN x SELF-COUPLING
%       GRID AND THE TWO 1-D NOISE SWEEPS RUN (its run_grids), and the
%       figure entries of those plus the two noise-level time-series
%       figures
%   the opt-in stage 'noise_switch'                  copies of the network
%       that share the noise for 10 s, then run on noise seeds of their own
%       (run_noise_switch_rollouts; Tom, 2026-11-14: "it should be added to
%       our set of analyses and figures"), on cfg.preset_name, the bundle's
%       conditions, sized by the run mode (noise_switch_stage_plan's table:
%       4 copies x 5 s at 'essential', 8 x 10 s at 'fast', 16 x 20 s at
%       'medium', 32 x 20 s at 'production'), into <run_dir>/noise_switch/
% Block C's stage and then 'noise_switch' are added to block A's stages
% (they run after them, as every opt-in stage does), block C's settings are
% carried over, and block C's figure entries and then the noise-switch
% entry are appended to block A's. Neither of the two configs changes: each
% still runs alone, into its own directories, and neither runs the
% noise-switch stage.
%
% WHY ONE BUNDLE and not the two run one after the other. One call of
% run_all_paper_analyses pre-flights everything (the run mode, the stage
% selection, block C's grids and variant presets, the pool allocation, the
% empty run directory) before any compute; a stage that fails does not cost
% the others; and one run directory with one manifest is what the figures
% of both blocks are drawn from. A queue of the two bundles would check the
% second only after the first had run.
%
% What differs from the parts:
%   cfg.stages    mu7chain_config's selection (every stage where it names
%                 none) followed by 'mc_gain_sweep' and 'noise_switch'
%   cfg.run_dir   data/mu7chain_AC_<run_mode>
%   cfg.fig_root  figs/mu7chain_AC_<run_mode>_<figure_set>; the grouped
%                 main figures read their source figures from this root
%   cfg.figures   block A's entries (main Figures 1 to 6, S1 to S8 and the
%                 appendix), then block C's (all appendix since Plan 2; the
%                 old Fig 7, fig_main_mc_grouped, among them), then the
%                 paper's Figure 7 and S9 (fig_main7_grouped,
%                 fig_supp9_grouped; from 'fast' upward), then
%                 fig_noise_switch_rollouts: one time-series figure per
%                 condition of the set (Fig_noise_switch_timeseries_<c>_copies)
%                 and the summary (Fig_noise_switch_divergence, with its
%                 table), read from <run_dir>/noise_switch/
% A bundle may resize the noise-switch stage with cfg.noise_switch
% (n_copies, t_shared, t_rollout, conditions, model_args;
% noise_switch_stage_plan). This one does not set it.
%
% HOW LONG, on 12 pool workers: block A's stages (about 35 minutes at
% 'essential'; mu7chain_config) and then block C's, since FD-029 the
% self-coupling grid and the two 1-D noise sweeps at 10.24 worker-minutes a
% job of four conditions, 1.35 times that with the release noise on
% (measured on the fast run of 2026-11-13; mu7chain_gainbias_config has the
% table), and three quarters of that for a job of three conditions (FD-030,
% an estimate until measured): 118 jobs of four conditions at 'essential',
% about 1 h 44; of three, 189 at 'fast', about 2 h 06; 715 at 'fastmedium',
% about 7 h 50; 2415 at 'medium', about 26 h 09. Then the
% noise switch (0.5 min at 'fast' on 4 workers, measured 2026-11-14; its
% figure entry about 45 s at 'fast'). The figures add the two noise-level
% time series, about a minute each. So an 'essential' run of this bundle is
% about 2 h 30. (Until FD-029 block C was 51 jobs, about 45 to 59 min, at
% 'essential' and 777, about 12 h 30, at 'fast'.)
%
% 'fastmedium' (Tom, 2026-11-14) WAS SIZED FOR ABOUT 19 HOURS on 12
% workers, with block C at 9 x 9 x 2 and five grids; the run of 2026-11-15
% took 19 h 28. Since FD-029 block C is 715 jobs, about 10 h 27 with four
% conditions, so the bundle was about 15 h 10 (block A 4 h 15, block C
% 10 h 27, figures 0 h 30). SINCE FD-030 every stage simulates three
% conditions and one figure set is drawn: about 11 h 15, block A about
% 3 h 10 (three quarters of 4 h 15, on the assumption that each stage's
% cost is proportional to its conditions; not measured), block C about
% 7 h 50 (the stage plan's estimate), figures about 0 h 15 (one pass).
% The table below is the
% sizing of 2026-11-14. Each stage's estimate is the fast run's measured time
% (2026-11-13) scaled by the ratio of work; the sizes are in each stage's
% table, which gives the reasoning:
%   sweeps              84 min   1-D 11 levels x 7 reps (28), tau 11 x 7 on
%                                medium's 50-s window (16), joint 4 levels,
%                                192 points (40)
%   memory_capacity     15       medium's 15 trials (two rounds)
%   mc_example           5       medium's protocol
%   eig_heatmap         18       T 60 s, 80 samples, examples 60
%   numerics            18       L and C on medium's windows, 5 / 6 trials
%   lyapunov_spectrum    8       T 40 s, 1 seed, K 50
%   transient_gain      76       fast's size
%   local_lyapunov       8       K 50, the bundle's own
%   dc_lle              19       medium's 15 seeds (two rounds)
%   paper_illustrations  3
%   mc_gain_sweep      819       9 x 9 levels x 2 networks: 846 jobs, 324
%                                with release noise (13 h 40); since FD-029
%                                627 (715 jobs, 55 with release noise)
%   noise_switch         3       medium's 16 copies x 20 s
%   two figure passes   30       about 10 min each at 'fast'
% about 18 h 25 in all: block A 4 h 15, block C 13 h 40, figures 0 h 30
% (every row with four conditions and two figure sets, as of then).
% Agents run 'essential' and 'fast'; 'fastmedium' and larger are Tom's to
% start.
%
% See also: mu7chain_AC_run, mu7chain_config, mu7chain_gainbias_config,
%           run_all_paper_analyses, paper_opt_in_stages

arguments
    run_mode (1,:) char = 'essential'
    figure_set (1,:) char = ''   % '' = the mode's own set; checked by mu7chain_config
end

A = mu7chain_config(run_mode, figure_set);            % also checks the run mode and the set
figure_set = A.figure_set;                            % resolved: the mode's own when none was named
C = mu7chain_gainbias_config(run_mode, figure_set);

cfg = A;

%% Stages: block A's, then block C's opt-in stage, then the noise switch
% An empty selection means every stage of the first list, and never an
% opt-in stage, so to add one the first list is written out.
stages_A = A.stages;
if isempty(stages_A); stages_A = paper_stage_names(); end
cfg.stages = [stages_A, C.stages, {'noise_switch'}];
cfg.mc_gain_sweep = C.mc_gain_sweep;
% Block A's memory-capacity stages take the lean path here: the same R^2
% and total memory capacity to the bit (test_mc_lean_path), at tens of MB
% per trial instead of about 17 GB. At 'fast' and above the bundle runs
% every block A stage, memory capacity among them, beside block C's grids
% on the same 12 workers, and the full-trajectory path would not fit in
% memory. mu7chain_config on its own is unchanged (lead, FD-025).
cfg.mc_lean = true;
% The noise-switch rollouts store by whole-sample decimation of the
% integration rate, so the storage rates must divide cfg.fs = 500 (Tom,
% 2026-11-17): every other sample (250 Hz), and every fourth for the
% example series (125 Hz). The old 200 and 100 Hz do not divide 500.
cfg.noise_switch = struct('store_fs', 250, 'example_fs', 125);

%% Where things land -- apart from either part's own directories
cfg.run_dir  = ['data/mu7chain_AC_' run_mode];
cfg.fig_root = ['figs/mu7chain_AC_' run_mode '_' figure_set];

%% The figure registry: block A's entries, then block C's
% The grouped main figures read their source figures from the figure root
% (grouped_figure_registry bakes it into each entry), so the entries built
% for mu7chain_config's root are pointed at this bundle's.
figures_A = A.figures;
for k = 1:numel(figures_A)
    args = figures_A{k}.args;
    at = find(strcmp(args(1:2:end), 'source_fig_root'));
    args(2 * at) = {cfg.fig_root};
    figures_A{k}.args = args;
end
% The noise-switch entry: one call draws a time-series figure per condition
% of the set and the summary, from <run_dir>/noise_switch/ (an error when
% the run does not hold it). Its tags, Fig_noise_switch_timeseries_<c>_copies
% and Fig_noise_switch_divergence, are prefixes of no other tag, and it has
% its own folder. In the appendix (the figure plan of 2026-11-21).
noise_switch = paper_figure_entry('fig_noise_switch_rollouts', @fig_noise_switch_rollouts, 'appendix');
% THE PAPER'S FIGURE 7 AND S9 (Plan 2, integration, 2026-11-21): the two
% memory composites read block A (the tau sweep, memory_capacity,
% mc_example) and block C (the gain x self-coupling grid and, for Figure 7,
% the Wiener noise sweep) from THIS run directory (their mc_run_dir is left
% empty, so it is run_dir), which is why they are entries of this bundle
% and of neither part. Registered only where the run holds all of it: not
% at 'essential', which skips block A's tau sweep and memory stages, and
% not when block C's run_grids leave out either grid. run_mode and
% figure_set name block C's registry entry whose map style they copy.
mc_args = {'preset_name', cfg.preset_name, 'run_mode', run_mode, 'figure_set', figure_set};
block_A_memory = all(ismember({'memory_capacity', 'mc_example'}, cfg.stages)) && ...
    (isempty(cfg.sweeps) || ismember('tau', cfg.sweeps));
grids = cfg.mc_gain_sweep.run_grids;
memory = {};
if block_A_memory && all(ismember({'selfcoupling', 'noise1d_wiener'}, grids))
    memory = {paper_figure_entry('fig_main7_grouped', @fig_main7_memory, 'main', [mc_args, {'fig_tag', 'Fig_Main7_Grouped'}]), ...
              paper_figure_entry('fig_supp9_grouped', @fig_supp9_grouped, 'supplement', [mc_args, {'fig_tag', 'Fig_Supp9_Grouped'}])};
end
cfg.figures = [figures_A, C.figures, memory, {noise_switch}];
names = cellfun(@(f) f.name, cfg.figures, 'UniformOutput', false);
[~, first] = unique(names, 'stable');
shared = unique(names(setdiff(1:numel(names), first)));
assert(isempty(shared), 'mu7chain_AC_config:DuplicateEntry', ...
    'Two parts of the bundle have a figure entry named %s; each entry needs its own folder.', ...
    strjoin(shared, ', '));
end
