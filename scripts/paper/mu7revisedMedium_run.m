% MU7REVISEDMEDIUM_RUN The mu7revised network at MEDIUM: overnight, unattended (TR, 2026-09-16).
%
%   Open this file and press Run, or launch it through the MATLAB MCP server.
%   setup_paths is called on the first line. Every setting comes from
%   mu7revisedMedium_config(), which states all of them itself.
%
%   run_dir   data/mu7revisedMedium      fig_root  figs/mu7revisedMedium
%
% RERUNNING REQUIRES DELETING data/mu7revisedMedium FIRST.
% The transcript is written to <run_dir>/command_window.log by vlog.
setup_paths();
cfg = mu7revisedMedium_config();
run_dir = run_all_paper_analyses(cfg);
results = make_all_paper_figures(cfg);
fprintf('\n========================================================\n');
fprintf('mu7revised MEDIUM RUN COMPLETE\n');
fprintf('  preset  : %s (%s)\n', cfg.preset_name, cfg.run_mode);
fprintf('  run_dir : %s\n', run_dir);
fprintf('  figures : %d of %d succeeded\n', sum([results.ok]), numel(results));
fprintf('========================================================\n');
