function cfg = mu7chain_AC_config(run_mode, figure_set)
% MU7CHAIN_AC_CONFIG Blocks A and C as one bundle: mu7chain_config's stages and figures, then block C's memory-capacity grids and theirs.
%
%   cfg = mu7chain_AC_config()                       % 'essential', figure set '4cond'
%   cfg = mu7chain_AC_config('fast', '3cond')
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
%       stage 'mc_gain_sweep' with its two grids and two variant networks
%       (cfg.mc_gain_sweep), and its figure registry
%   the opt-in stage 'noise_switch'                  copies of the network
%       that share the noise for 10 s, then run on noise seeds of their own
%       (run_noise_switch_rollouts; Tom, 2026-11-14: "it should be added to
%       our set of analyses and figures"), on cfg.preset_name, every
%       condition, sized by the run mode (noise_switch_stage_plan's table:
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
%   cfg.figures   block A's entries, then block C's, then
%                 fig_noise_switch_rollouts: one time-series figure per
%                 condition of the set (Fig_noise_switch_timeseries_<c>_copies)
%                 and the summary (Fig_noise_switch_divergence, with its
%                 table), read from <run_dir>/noise_switch/
% A bundle may resize the noise-switch stage with cfg.noise_switch
% (n_copies, t_shared, t_rollout, conditions, model_args;
% noise_switch_stage_plan). This one does not set it.
%
% HOW LONG, on 12 pool workers: block A's stages (about 35 minutes at
% 'essential'; mu7chain_config) and then block C's (about 20 minutes at
% 'essential', about 4 h 25 at 'fast'; mu7chain_gainbias_config). Agents
% run 'essential' and 'fast'; 'medium' and larger are Tom's to start.
%
% See also: mu7chain_AC_run, mu7chain_config, mu7chain_gainbias_config,
%           run_all_paper_analyses, paper_opt_in_stages

arguments
    run_mode (1,:) char = 'essential'
    figure_set (1,:) char {mustBeMember(figure_set, {'3cond', '4cond'})} = '4cond'
end

A = mu7chain_config(run_mode, figure_set);            % also checks the run mode
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
% its own folder. Not in the manuscript (a possible supplemental figure).
noise_switch = struct('name', 'fig_noise_switch_rollouts', 'fn', @fig_noise_switch_rollouts, ...
    'in_paper', false, 'args', {{}});
cfg.figures = [figures_A, C.figures, {noise_switch}];
names = cellfun(@(f) f.name, cfg.figures, 'UniformOutput', false);
[~, first] = unique(names, 'stable');
shared = unique(names(setdiff(1:numel(names), first)));
assert(isempty(shared), 'mu7chain_AC_config:DuplicateEntry', ...
    'Two parts of the bundle have a figure entry named %s; each entry needs its own folder.', ...
    strjoin(shared, ', '));
end
