% Fig_memory_capacity.m
% Stability_Manuscript presentation figure: memory capacity.
% Replots a saved looped_memory_capacity run (from data/memory_capacity/
% paper_ready) into this presentation folder using the existing plotting
% pipeline. No simulation is re-run -- edit mat_file to point at a different run.

close all
clear

this_dir     = fileparts(mfilename('fullpath'));
% scripts/presentations/Stability_Manuscript/fig_memory_capacity -> root is 4 levels up
project_root = fileparts(fileparts(fileparts(fileparts(this_dir))));

% Resolve replot_memory_capacity / plot_memory_capacity, capture_git_provenance,
% and the combined-figure helper (plot_memory_capacity_combined, this folder).
setup_paths();

% Source run to replot.
mat_file = fullfile(project_root, 'data', 'memory_capacity', 'paper_ready', ...
    'MC_sample_hold_20260722_154245_trials30_results.mat');

% Write the regenerated manuscript figures next to this script.
out_dir = this_dir;

% Fig1 (total-MC + horizon distributions) and Fig2 (per-delay + cumulative) are
% shown for reference but NOT saved -- only the final combined Fig3 is written.
% Pass '' so replot_memory_capacity displays without saving (plot_memory_capacity
% guards its save on ~isempty(out_dir)).
replot_memory_capacity(mat_file, '');

% Load the saved run once (used for the combined figure and the README below).
S       = load(mat_file, 'results_all');
run_tag = S.results_all.run_tag;
cfg     = S.results_all.settings;

% Fig3: paper-ready 1x3 strip assembled from pieces of Fig1/Fig2 --
% (a) cumulative memory, (b) per-delay memory, (c) horizon paired trials.
fig3 = plot_memory_capacity_combined(S.results_all, out_dir);

% Log the git state alongside the figures so this presentation output can be
% traced back to an exact commit (+ working-tree diff if dirty).
capture_git_provenance(out_dir, project_root);

%% -------------------- Human-readable description --------------------
% Write a plain-text record of how this figure was produced: the source
% dataset, the pipeline used, and the names of the figure files generated.
% Fig3 is saved by save_some_figs_to_folder_2, which appends the figure number:
% <name>_figure_<N>.png/.svg and <name>_f_<N>.fig
fig3_tag = 'Fig_Memory_Capacity';
fig3_num = num2str(fig3.Number);
fig_files = { ...
    [fig3_tag '_figure_' fig3_num '.png']; ...
    [fig3_tag '_figure_' fig3_num '.svg']; ...
    [fig3_tag '_f_' fig3_num '.fig'] };

desc_path = fullfile(out_dir, 'README_fig_memory_capacity.txt');
fid = fopen(desc_path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Stability_Manuscript figure: Memory Capacity\n');
fprintf(fid, '============================================\n\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now')));
fprintf(fid, 'By script: %s.m\n\n', mfilename);

fprintf(fid, 'HOW IT WAS MADE\n');
fprintf(fid, '  This is a presentation replot -- no simulation is re-run. The script\n');
fprintf(fid, '  loads a saved looped_memory_capacity.m results file and calls\n');
fprintf(fid, '  replot_memory_capacity -> plot_memory_capacity to regenerate the\n');
fprintf(fid, '  figures into this folder. See git_provenance.txt for the exact commit.\n\n');

fprintf(fid, 'SOURCE DATASET\n');
fprintf(fid, '  %s\n', mat_file);
fprintf(fid, '  run_tag: %s\n\n', run_tag);

fprintf(fid, 'KEY SETTINGS (from the saved run)\n');
fprintf(fid, '  input_type    : %s\n', cfg.input_type);
fprintf(fid, '  n (neurons)   : %d\n', cfg.n);
fprintf(fid, '  fs (Hz)       : %d\n', cfg.fs);
fprintf(fid, '  n_trials      : %d\n', cfg.n_trials);
fprintf(fid, '  level_of_chaos: %g\n', cfg.level_of_chaos);
fprintf(fid, '  T_hold (s)    : %g\n', cfg.T_hold);
fprintf(fid, '  conditions    : %s\n\n', strjoin(S.results_all.conditions, ', '));

fprintf(fid, 'FIGURES PRODUCED (in this folder)\n');
for k = 1:numel(fig_files)
    fprintf(fid, '  %s\n', fig_files{k});
end
fprintf(fid, '\n  Only the final combined figure is saved. Fig1 (paired total-MC +\n');
fprintf(fid, '  memory-horizon distributions) and Fig2 (per-delay R^2(d) and cumulative\n');
fprintf(fid, '  MC) are shown for reference during the replot but not written to disk.\n');
fprintf(fid, '  Fig3 = 1x3 combined panel: (a) cumulative MC, (b) per-delay R^2,\n');
fprintf(fid, '         (c) horizon paired trials  [paper-ready]\n');

clear cleanup;  % flush + close
fprintf('Description written: %s\n', desc_path);
