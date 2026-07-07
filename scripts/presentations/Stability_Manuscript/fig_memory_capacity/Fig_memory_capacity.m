% Fig_memory_capacity.m
% Stability_Manuscript presentation figure: memory capacity.
% Replots a saved looped_memory_capacity run (from data/memory_capacity/
% paper_ready) into this presentation folder using the existing plotting
% pipeline. No simulation is re-run -- edit mat_file to point at a different run.

this_dir     = fileparts(mfilename('fullpath'));
% scripts/presentations/Stability_Manuscript/fig_memory_capacity -> root is 4 levels up
project_root = fileparts(fileparts(fileparts(fileparts(this_dir))));

% Make replot_memory_capacity / plot_memory_capacity resolvable regardless of
% cwd (replot_memory_capacity in turn adds src/ + its own folder), and src/ so
% capture_git_provenance is on the path.
addpath(fullfile(project_root, 'scripts', 'memory_capacity'));
addpath(fullfile(project_root, 'src'));

% Source run to replot.
mat_file = fullfile(project_root, 'data', 'memory_capacity', 'paper_ready', ...
    'MC_sample_hold_20260706_194403_trials50_results.mat');

% Write the regenerated manuscript figures next to this script.
out_dir = this_dir;

replot_memory_capacity(mat_file, out_dir);

% Log the git state alongside the figures so this presentation output can be
% traced back to an exact commit (+ working-tree diff if dirty).
capture_git_provenance(out_dir, project_root);

%% -------------------- Human-readable description --------------------
% Write a plain-text record of how this figure was produced: the source
% dataset, the pipeline used, and the names of the figure files generated.
S       = load(mat_file, 'results_all');
run_tag = S.results_all.run_tag;
cfg     = S.results_all.settings;

fig_files = { ...
    [run_tag '_Fig1_MC_Distributions.png']; ...
    [run_tag '_Fig1_MC_Distributions.pdf']; ...
    [run_tag '_Fig2_R2_Curves.png']; ...
    [run_tag '_Fig2_R2_Curves.pdf'] };

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
fprintf(fid, '\n  Fig1 = paired total-MC + memory-horizon distributions\n');
fprintf(fid, '  Fig2 = per-delay R^2(d) and cumulative MC (mean +/- 95%% CI)\n');

clear cleanup;  % flush + close
fprintf('Description written: %s\n', desc_path);
