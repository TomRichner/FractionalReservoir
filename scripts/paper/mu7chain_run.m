function [run_dir, fig_results] = mu7chain_run(run_mode)
% MU7CHAIN_RUN The chained-depression bundle: one run, then its figure sets (2026-11-12, FD-007).
%
%   mu7chain_run                 % 'fast': open this file and press Run, or
%                                % launch it through the MATLAB MCP server
%   mu7chain_run('essential')    % a first look: figures 2-4, brief sweeps,
%                                % the four-condition set only (FD-023)
%   [run_dir, fig_results] = mu7chain_run(run_mode)
%
%   setup_paths is called on the first line. Every setting comes from
%   mu7chain_config, which states all of them itself. The session needs a
%   pool allocation first (pool_set(N)); both entry points refuse without.
%
%   run_dir   data/mu7chain_<run_mode>
%   figures   figs/mu7chain_<run_mode>_3cond   (no adaptation, 1TS, chained MTS: the paper's)
%             figs/mu7chain_<run_mode>_4cond   (the same plus product MTS: for review)
%
% The analyses run ONCE, over all four conditions; the figures are drawn
% once per set from that one run directory. Which sets is the config's
% cfg.figure_sets: both, except 'essential', which draws the four-condition
% set alone.
%
% WAS A SCRIPT with run_mode = 'fast' written in it (until FD-023). It is a
% function with that as its default so that another mode is an argument and
% not an edit to this file; run with no argument it does what the script
% did. The one difference: the script left cfg, run_dir and fig_results in
% the base workspace, and the function returns run_dir and fig_results.
%
% Agents run 'essential' and 'fast' only. 'medium' and larger are Tom's to start.
% RERUNNING REQUIRES DELETING data/mu7chain_<run_mode> FIRST.
% The analysis transcript is written to <run_dir>/command_window.log and each
% figure pass's to <fig_root>/command_window.log, by vlog.
arguments
    run_mode (1,:) char = 'fast'
end
setup_paths();

cfg = mu7chain_config(run_mode);
run_dir = run_all_paper_analyses(cfg);

sets = cfg.figure_sets;
fig_results = cell(1, numel(sets));
for k = 1:numel(sets)
    % Rebuilt per set: the grouped registry bakes in its own fig_root.
    fig_results{k} = make_all_paper_figures(mu7chain_config(run_mode, sets{k}));
end

fprintf('\n========================================================\n');
fprintf('mu7chain RUN COMPLETE\n');
fprintf('  preset  : %s (%s)\n', cfg.preset_name, cfg.run_mode);
fprintf('  run_dir : %s\n', run_dir);
for k = 1:numel(sets)
    r = fig_results{k};
    fprintf('  figures %s : %d of %d succeeded  (%s)\n', sets{k}, sum([r.ok]), numel(r), ...
        mu7chain_config(run_mode, sets{k}).fig_root);
end
fprintf('========================================================\n');
if nargout == 0; clear run_dir fig_results; end
end
