close all
clc

% Fig_param_space.m
% Stability_Manuscript presentation figure: parameter-space distributions.
% Recreates the param-space histogram figures (LLE + mean_rate distributions)
% from a saved run_all_<dt> run and writes the final figures into this
% presentation folder. No simulation is re-run -- edit data_root to point at a
% different run_all output.

this_dir     = fileparts(mfilename('fullpath'));
% .../Stability_Manuscript/fig_param_space -> project root is 4 up
project_root = fileparts(fileparts(fileparts(fileparts(this_dir))));

% Resolve the replot pipeline (scripts/run_all_analyses/replot) + src helpers
% (save_some_figs_to_folder_2, capture_git_provenance).
addpath(genpath(fullfile(project_root, 'scripts')));
addpath(genpath(fullfile(project_root, 'src')));

% Source run (a run_all_<dt> folder with a param_space_* subdir).
data_root = fullfile(project_root, 'data', 'param_space', 'run_all_jul_06_26_22_00');
out_dir   = this_dir;   % write the final figures next to this script

% Start from a clean slate: replot_param_space_analysis saves ALL open figures,
% so any stray figure lingering in the session would pollute the save.
close all force;

% 1) Regenerate the param-space histogram figures (LLE + mean_rate) into a
%    replot_param_space_<dt>/figures/ folder under data_root.
replot_dir = replot_param_space_analysis(data_root);

% 2) Re-open each saved figure, map it to a stable output name by its Name
%    (robust vs the run-dependent figure number), and save png/svg/fig into the
%    presentation folder.
name_map = containers.Map( ...
    {'LLE Distribution', 'mean_rate Distribution'}, ...
    {'Fig_ParamSpace_LLE', 'Fig_ParamSpace_MeanRate'});

fig_listing = dir(fullfile(replot_dir, 'figures', '*.fig'));
saved_bases = {};
for k = 1:numel(fig_listing)
    f = openfig(fullfile(fig_listing(k).folder, fig_listing(k).name), 'visible');
    nm = get(f, 'Name');
    if ~isKey(name_map, nm)
        close(f);
        continue;
    end
    base = name_map(nm);

    % --- cleanup / post-processing goes here (follow-up styling pass) --------

    % Save with a STABLE name: save_some_figs_to_folder_2 suffixes filenames
    % with the figure number, so save then rename to fixed names, after clearing
    % any prior <base> outputs so nothing stale lingers.
    old = dir(fullfile(out_dir, [base '*']));
    for a = 1:numel(old)
        fp = fullfile(old(a).folder, old(a).name);
        if ~old(a).isdir && (endsWith(fp, '.png') || endsWith(fp, '.svg') || endsWith(fp, '.fig'))
            delete(fp);
        end
    end
    save_some_figs_to_folder_2(out_dir, base, f.Number, {'png', 'svg', 'fig'});
    num = num2str(f.Number);
    movefile(fullfile(out_dir, [base '_figure_' num '.png']), fullfile(out_dir, [base '.png']));
    movefile(fullfile(out_dir, [base '_figure_' num '.svg']), fullfile(out_dir, [base '.svg']));
    movefile(fullfile(out_dir, [base '_f_' num '.fig']),      fullfile(out_dir, [base '.fig']));
    saved_bases{end+1} = base; %#ok<SAGROW>
end

% 3) Prep figures exist only to build the finals -- remove the whole replot
%    folder so no extra figs are left behind in the data dir.
if isfolder(replot_dir)
    rmdir(replot_dir, 's');
end

% 4) Log the git state alongside the figures so this output can be traced back
%    to an exact commit (+ working-tree diff if dirty).
capture_git_provenance(out_dir, project_root);

%% -------------------- Human-readable description --------------------
% Write a plain-text record of how these figures were produced.
ps_dirs = dir(fullfile(data_root, 'param_space_*'));
ps_dirs = ps_dirs([ps_dirs.isdir]);

desc_path = fullfile(out_dir, 'README_fig_param_space.txt');
fid = fopen(desc_path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Stability_Manuscript figure: Parameter-space distributions\n');
fprintf(fid, '==========================================================\n\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now')));
fprintf(fid, 'By script: %s.m\n\n', mfilename);

fprintf(fid, 'HOW IT WAS MADE\n');
fprintf(fid, '  Presentation replot -- no simulation is re-run. The script reloads the\n');
fprintf(fid, '  saved param-space PSA object from a run_all_<dt> run and calls\n');
fprintf(fid, '  replot_param_space_analysis (psa.plot for LLE + mean_rate) to rebuild the\n');
fprintf(fid, '  distribution histograms (1 row, one column per adaptation condition),\n');
fprintf(fid, '  then saves them here. See git_provenance.txt for the exact commit.\n\n');

fprintf(fid, 'SOURCE RUN\n');
fprintf(fid, '  %s\n', data_root);
fprintf(fid, '  param_space subfolder(s) used:\n');
for k = 1:numel(ps_dirs)
    fprintf(fid, '    %s\n', ps_dirs(k).name);
end
fprintf(fid, '\n');

fprintf(fid, 'FIGURES PRODUCED (in this folder)\n');
for k = 1:numel(saved_bases)
    fprintf(fid, '  %s.png / .svg / .fig\n', saved_bases{k});
end
fprintf(fid, '\n  Fig_ParamSpace_LLE      = LLE (growth-rate) distribution histograms\n');
fprintf(fid, '  Fig_ParamSpace_MeanRate = mean firing-rate distribution histograms\n');
fprintf(fid, '  (one column per condition: No Adaptation, SFA, STD, SFA+STD).\n');

clear cleanup;  % flush + close
fprintf('Description written: %s\n', desc_path);
