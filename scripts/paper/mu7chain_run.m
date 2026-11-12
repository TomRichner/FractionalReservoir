% MU7CHAIN_RUN The chained-depression bundle: one run, two figure sets (2026-11-12, FD-007).
%
%   Open this file, set run_mode, and press Run, or launch it through the
%   MATLAB MCP server. setup_paths is called on the first line. Every setting
%   comes from mu7chain_config, which states all of them itself.
%
%   run_dir   data/mu7chain_<run_mode>
%   figures   figs/mu7chain_<run_mode>_3cond   (no adaptation, 1TS, chained MTS: the paper's)
%             figs/mu7chain_<run_mode>_4cond   (the same plus product MTS: for review)
%
% The analyses run ONCE, over all four conditions; the figures are drawn
% twice from that one run directory, once per set.
%
% Agents run 'fast' only. 'medium' and larger are Tom's to start.
% RERUNNING REQUIRES DELETING data/mu7chain_<run_mode> FIRST.
% The analysis transcript is written to <run_dir>/command_window.log and each
% figure pass's to <fig_root>/command_window.log, by vlog.
setup_paths();
run_mode = 'fast';

cfg = mu7chain_config(run_mode);
run_dir = run_all_paper_analyses(cfg);

sets = {'3cond', '4cond'};
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
