function [run_dir, fig_results] = mu7chain_AC_run(run_mode)
% MU7CHAIN_AC_RUN Blocks A and C in one command: one run, then its figure sets.
%
%   mu7chain_AC_run                 % 'essential': about an hour on 12 workers
%   mu7chain_AC_run('fast')
%   [run_dir, fig_results] = mu7chain_AC_run(run_mode)
%
%   setup_paths is called on the first line. Every setting comes from
%   mu7chain_AC_config, which composes mu7chain_config (block A) and
%   mu7chain_gainbias_config (block C). The session needs a pool allocation
%   first (pool_set(N)); both entry points refuse without.
%
%   run_dir   data/mu7chain_AC_<run_mode>: block A's stage folders, and
%             block C's mc_gain_sweep_bias/, mc_gain_sweep_selfcoupling/,
%             mc_gain_sweep_bias_EI3to2/ and mc_gain_sweep_bias_EI2to3/
%   figures   figs/mu7chain_AC_<run_mode>_3cond   (no adaptation, 1TS, chained MTS: the paper's)
%             figs/mu7chain_AC_<run_mode>_4cond   (the same plus product MTS: for review)
%             'essential' draws the four-condition set alone (cfg.figure_sets)
%
% EVERYTHING IS CHECKED BEFORE ANYTHING RUNS: run_all_paper_analyses
% pre-flights the run mode, the stage selection, block C's grids and
% variant networks, the pool allocation and the run directory, which must
% be absent or empty. The analyses then run ONCE, over all four conditions,
% and the figures are drawn once per set from that one run directory.
%
% Agents run 'essential' and 'fast'. 'medium' and larger are Tom's to start.
% RERUNNING REQUIRES DELETING data/mu7chain_AC_<run_mode> FIRST.
% The analysis transcript is written to <run_dir>/command_window.log and each
% figure pass's to <fig_root>/command_window.log, by vlog.
%
% See also: mu7chain_AC_config, mu7chain_run, mu7chain_gainbias_run
arguments
    run_mode (1,:) char = 'essential'
end
setup_paths();

cfg = mu7chain_AC_config(run_mode);
run_dir = run_all_paper_analyses(cfg);

sets = cfg.figure_sets;
fig_results = cell(1, numel(sets));
for k = 1:numel(sets)
    % Rebuilt per set: the registry's grouped entries name their own fig_root.
    fig_results{k} = make_all_paper_figures(mu7chain_AC_config(run_mode, sets{k}));
end

fprintf('\n========================================================\n');
fprintf('mu7chain BLOCKS A AND C RUN COMPLETE\n');
fprintf('  preset  : %s (%s)\n', cfg.preset_name, cfg.run_mode);
fprintf('  run_dir : %s\n', run_dir);
for k = 1:numel(sets)
    r = fig_results{k};
    fprintf('  figures %s : %d of %d succeeded  (%s)\n', sets{k}, sum([r.ok]), numel(r), ...
        mu7chain_AC_config(run_mode, sets{k}).fig_root);
end
fprintf('========================================================\n');
if nargout == 0; clear run_dir fig_results; end
end
