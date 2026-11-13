function [run_dir, fig_results] = mu7chain_gainbias_run(run_mode)
% MU7CHAIN_GAINBIAS_RUN Block C's bundle: the gain-by-bias memory-capacity sweep, then its figure sets (FD-010).
%
%   mu7chain_gainbias_run                % 'essential': about 30 minutes on 12 workers
%   mu7chain_gainbias_run('fast')        % about 3 h 20; a run for Tom to start
%   [run_dir, fig_results] = mu7chain_gainbias_run(run_mode)
%
%   setup_paths is called on the first line. Every setting comes from
%   mu7chain_gainbias_config. The session needs a pool allocation first
%   (pool_set(N)); both entry points refuse without.
%
%   run_dir   data/mu7chain_gainbias_<run_mode>, holding mc_gain_sweep/ and
%             the variant networks' mc_gain_sweep_EI3to2/ and mc_gain_sweep_EI2to3/
%   figures   figs/mu7chain_gainbias_<run_mode>_3cond   (the paper's three)
%             figs/mu7chain_gainbias_<run_mode>_4cond   (with product MTS)
%
% The stage runs ONCE, over all four conditions; the figures are drawn once
% per set from that one run directory (cfg.figure_sets, as in mu7chain_run).
%
% Agents run 'essential' only; 'fast' and larger are Tom's to start
% (mu7chain_gainbias_config gives the times).
% RERUNNING REQUIRES DELETING data/mu7chain_gainbias_<run_mode> FIRST (the
% entry point refuses a non-empty run directory). An interrupted stage can be
% resumed by calling run_mc_gain_sweep on the same output_dir directly.
arguments
    run_mode (1,:) char = 'essential'
end
setup_paths();

cfg = mu7chain_gainbias_config(run_mode);
run_dir = run_all_paper_analyses(cfg);

sets = cfg.figure_sets;
fig_results = cell(1, numel(sets));
for k = 1:numel(sets)
    fig_results{k} = make_all_paper_figures(mu7chain_gainbias_config(run_mode, sets{k}));
end

fprintf('\n========================================================\n');
fprintf('mu7chain GAIN x BIAS RUN COMPLETE\n');
fprintf('  preset  : %s (%s)\n', cfg.mc_preset, cfg.run_mode);
fprintf('  run_dir : %s\n', run_dir);
for k = 1:numel(sets)
    r = fig_results{k};
    fprintf('  figures %s : %d of %d succeeded  (%s)\n', sets{k}, sum([r.ok]), numel(r), ...
        mu7chain_gainbias_config(run_mode, sets{k}).fig_root);
end
fprintf('========================================================\n');
if nargout == 0; clear run_dir fig_results; end
end
