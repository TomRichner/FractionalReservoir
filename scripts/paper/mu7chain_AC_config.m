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
% Block C's stage is added to block A's stages (it runs after them, as every
% opt-in stage does), its settings are carried over, and its figure entries
% are appended to block A's. Neither of those two configs changes: each
% still runs alone, into its own directories.
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
%                 none) followed by 'mc_gain_sweep'
%   cfg.run_dir   data/mu7chain_AC_<run_mode>
%   cfg.fig_root  figs/mu7chain_AC_<run_mode>_<figure_set>; the grouped
%                 main figures read their source figures from this root
%   cfg.figures   block A's entries, then block C's
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

%% Stages: block A's, then block C's opt-in stage
% An empty selection means every stage of the first list, and never an
% opt-in stage, so to add one the first list is written out.
stages_A = A.stages;
if isempty(stages_A); stages_A = paper_stage_names(); end
cfg.stages = [stages_A, C.stages];
cfg.mc_gain_sweep = C.mc_gain_sweep;

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
cfg.figures = [figures_A, C.figures];
names_A = cellfun(@(f) f.name, figures_A, 'UniformOutput', false);
names_C = cellfun(@(f) f.name, C.figures, 'UniformOutput', false);
shared = intersect(names_A, names_C);
assert(isempty(shared), 'mu7chain_AC_config:DuplicateEntry', ...
    'Blocks A and C both have a figure entry named %s; each entry needs its own folder.', ...
    strjoin(shared, ', '));
end
