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
