function [run_dir, fig_results] = mu7chain_gainbias_run(run_mode)
% MU7CHAIN_GAINBIAS_RUN Block C's bundle alone: the two memory-capacity grids, then their figure sets (FD-010, FD-024).
%
%   mu7chain_gainbias_run                % 'essential': about 20 minutes on 12 workers
%   mu7chain_gainbias_run('fast')        % about 4 h 25 on 12 workers
%   [run_dir, fig_results] = mu7chain_gainbias_run(run_mode)
%
%   setup_paths is called on the first line. Every setting comes from
%   mu7chain_gainbias_config. The session needs a pool allocation first
%   (pool_set(N)); both entry points refuse without.
%
%   run_dir   data/mu7chain_gainbias_<run_mode>, holding mc_gain_sweep_bias/,
%             mc_gain_sweep_selfcoupling/ and the variant networks'
%             mc_gain_sweep_bias_EI3to2/ and mc_gain_sweep_bias_EI2to3/
%   figures   figs/mu7chain_gainbias_<run_mode>_3cond   (the paper's three)
%             figs/mu7chain_gainbias_<run_mode>_4cond   (with product MTS)
%
% The stage runs ONCE, over all four conditions; the figures are drawn once
% per set from that one run directory (cfg.figure_sets, as in mu7chain_run).
% For block A's stages and figures in the same run, use mu7chain_AC_run.
%
% Agents run 'essential' and 'fast'; 'medium' and larger are Tom's to start
% (mu7chain_gainbias_config gives the sizes and times).
% RERUNNING REQUIRES DELETING data/mu7chain_gainbias_<run_mode> FIRST (the
% entry point refuses a non-empty run directory). An interrupted stage is
% resumed by run_mc_gain_sweep_stage(cfg, run_dir, mc_gain_sweep_stage_plan(cfg))
% with cfg = mu7chain_gainbias_config(run_mode): finished jobs are skipped.
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
fprintf('mu7chain GAIN x BIAS AND GAIN x SELF-COUPLING RUN COMPLETE\n');
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
