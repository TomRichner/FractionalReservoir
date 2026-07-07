% Fig_sensitivity_analysis.m
% Stability_Manuscript presentation figure: combined LLE sensitivity.
% Regenerates the stacked 1D-sensitivity LLE figure from a saved run_all_<dt>
% run and writes the final figure into this presentation folder. No simulation
% is re-run -- edit data_root to point at a different run_all output.

this_dir     = fileparts(mfilename('fullpath'));
% .../Stability_Manuscript/fig_sensitivity_analysis -> project root is 4 up
project_root = fileparts(fileparts(fileparts(fileparts(this_dir))));

% Resolve the replot pipeline (scripts/run_all_analyses/replot) + src helpers
% (save_some_figs_to_folder_2, capture_git_provenance).
addpath(genpath(fullfile(project_root, 'scripts')));
addpath(genpath(fullfile(project_root, 'src')));

% Source run (a run_all_<dt> folder with 1D_sensitivity_* subdirs).
data_root = fullfile(project_root, 'data', 'param_space', 'run_all_jul_06_26_22_00');
out_dir   = this_dir;   % write the final figure next to this script

% 1) Regenerate per-param LLE figs, 2) assemble the combined figure.
replot_dir = replot_sensitivity(data_root);
assemble_sensitivity_figure(replot_dir, 'LLE');

% Re-open the saved combined figure and save it into the presentation folder
% (png + svg + fig), matching the memory-capacity figure's formats/naming.
combined_fig_path = fullfile(replot_dir, 'figures', 'sensitivity_LLE_combined.fig');
cf = openfig(combined_fig_path, 'visible');
save_some_figs_to_folder_2(out_dir, 'Fig_Sensitivity_LLE', cf.Number, {'png', 'svg', 'fig'});

% Log the git state alongside the figure so this presentation output can be
% traced back to an exact commit (+ working-tree diff if dirty).
capture_git_provenance(out_dir, project_root);

%% -------------------- Human-readable description --------------------
% Write a plain-text record of how this figure was produced: the source run,
% the swept-parameter subfolders used, the pipeline, and the output filenames.
% save_some_figs_to_folder_2 appends the figure number:
%   <name>_figure_<N>.png/.svg and <name>_f_<N>.fig
fig_tag = 'Fig_Sensitivity_LLE';
fig_num = num2str(cf.Number);
fig_files = { ...
    [fig_tag '_figure_' fig_num '.png']; ...
    [fig_tag '_figure_' fig_num '.svg']; ...
    [fig_tag '_f_' fig_num '.fig'] };

sens_dirs = dir(fullfile(data_root, '1D_sensitivity_*'));
sens_dirs = sens_dirs([sens_dirs.isdir]);

desc_path = fullfile(out_dir, 'README_fig_sensitivity_analysis.txt');
fid = fopen(desc_path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Stability_Manuscript figure: LLE Sensitivity (combined)\n');
fprintf(fid, '=======================================================\n\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now')));
fprintf(fid, 'By script: %s.m\n\n', mfilename);

fprintf(fid, 'HOW IT WAS MADE\n');
fprintf(fid, '  Presentation replot -- no simulation is re-run. The script reloads the\n');
fprintf(fid, '  saved 1D-sensitivity PSA objects from a run_all_<dt> run and calls\n');
fprintf(fid, '  replot_sensitivity -> assemble_sensitivity_figure to rebuild the stacked\n');
fprintf(fid, '  LLE figure (rows = swept params, cols = adaptation conditions), then\n');
fprintf(fid, '  saves it here. See git_provenance.txt for the exact commit.\n\n');

fprintf(fid, 'SOURCE RUN\n');
fprintf(fid, '  %s\n', data_root);
fprintf(fid, '  1D_sensitivity subfolders used:\n');
for k = 1:numel(sens_dirs)
    fprintf(fid, '    %s\n', sens_dirs(k).name);
end
fprintf(fid, '\n');

fprintf(fid, 'FIGURES PRODUCED (in this folder)\n');
for k = 1:numel(fig_files)
    fprintf(fid, '  %s\n', fig_files{k});
end
fprintf(fid, '\n  Combined LLE sensitivity: one row per swept parameter\n');
fprintf(fid, '  (f, level_of_chaos, n), one column per adaptation condition.\n');

clear cleanup;  % flush + close
fprintf('Description written: %s\n', desc_path);
