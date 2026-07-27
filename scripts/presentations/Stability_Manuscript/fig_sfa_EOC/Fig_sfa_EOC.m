close all
clc

% Fig_sfa_EOC.m
% Stability_Manuscript presentation figure: SFA edge-of-chaos.
% Replots the tau_a-sensitivity LLE panel -- how the largest Lyapunov exponent
% (growth rate) approaches 0 (edge of chaos) as the maximum SFA adaptation
% timescale tau_a grows. No simulation is re-run: it reloads the saved
% tau_a_E_max PSA object and restyles the single-condition (SFA+STD) panel,
% matching the presentation treatment in fig_sensitivity_analysis.
%
% Edit data_root / tau_dir to point at a different run_all output.

this_dir     = fileparts(mfilename('fullpath'));
% .../Stability_Manuscript/fig_sfa_EOC -> project root is 4 up
project_root = fileparts(fileparts(fileparts(fileparts(this_dir))));

% Resolve save_some_figs_to_folder_2, capture_git_provenance, ParamSpaceAnalysis2.
addpath(genpath(fullfile(project_root, 'scripts')));
addpath(genpath(fullfile(project_root, 'src')));

% Source run + the specific tau_a_E_max sensitivity subfolder (SFA+STD only).
data_root = fullfile(project_root, 'data', 'param_space', 'run_all_jul_06_26_22_00');
tau_dir   = fullfile(data_root, ...
    'tau_sensitivity_tau_timescales_tau_a_E_max_nLevs_15_jul_07_26_08_49');
out_dir   = this_dir;   % write the final figure next to this script

% LLE (growth-rate) y-axis range -- matches the original tau_a run.
lle_range = [-0.3, 0.1];

% Start clean: plot_sensitivity operates on the current session's figures.
close all force;

% --- Reload the tau_a PSA and regenerate the single LLE panel --------------
psa_file = fullfile(tau_dir, 'psa_object.mat');
if ~isfile(psa_file)
    error('Fig_sfa_EOC:NoPSA', 'Missing PSA object:\n  %s', psa_file);
end
S   = load(psa_file);
psa = S.psa_tau_a;

psa.plot_sensitivity('metric', 'LLE', 'hist_range', lle_range);
cf = gcf;

% --- Presentation restyle (same knobs as Fig_sensitivity_analysis.m) --------
tick_fs   = 14;    % tick numbers
label_fs  = 15.4;  % axis labels
clim_frac = 0.4;   % darken imagesc: cap CLim at total_reps*clim_frac
% Colormap ramps white (0 counts) -> 90% black (max), not pure black, so the
% blue median line stays visible over the darkest cells.
dark_cmap    = repmat(linspace(1, 0.1, 256)', 1, 3);
median_alpha = 0.35;   % blue median line transparency
median_lw    = 3;      % blue median line width
zeroline_lw  = 2;      % green dashed zero line width

ax = findobj(cf, 'Type', 'axes');
set(ax, 'FontSize', tick_fs);
box(ax, 'off');
colormap(ax, dark_cmap);

% Darken the histogram density (shared black->white scale from plot_sensitivity).
cl = get(ax, 'CLim');
set(ax, 'CLim', [0, cl(2) * clim_frac]);

% Blue median line: more transparent + thinner. (imagesc is Type 'image', the
% zero line is 'constantline', so 'line' is the median.)
ml = findobj(ax, 'Type', 'line');
for m = 1:numel(ml)
    mc = get(ml(m), 'Color');
    if numel(mc) < 4; mc(4) = 1; end
    mc(4) = median_alpha;
    set(ml(m), 'Color', mc, 'LineWidth', median_lw);
end
% Green dashed zero line: thinner.
set(findobj(ax, 'Type', 'constantline'), 'LineWidth', zeroline_lw);

% Labels: lambda_1 -> "Growth Rate"; x -> "max \tau_a (s)"; drop the title.
ylabel(ax, 'Growth Rate', 'Interpreter', 'none', 'FontSize', label_fs);
xlabel(ax, 'max $\tau_a$ (s)', 'Interpreter', 'latex', 'FontSize', label_fs);
title(ax, '');

% Fewer y-ticks + a tighter view around the data.
ylim(ax, [-0.25, 0.05]);
yticks(ax, -0.2:0.1:0.1);

% --- Save ONLY the cleaned panel, with a STABLE name -----------------------
% save_some_figs_to_folder_2 suffixes filenames with the (run-dependent) figure
% number; save, then rename to fixed names so re-runs overwrite cleanly. First
% clear any prior Fig_SFA_EOC outputs (stable or numbered) so nothing stale
% lingers.
fig_tag = 'Fig_SFA_EOC';
old = dir(fullfile(out_dir, [fig_tag '*']));
for a = 1:numel(old)
    fp = fullfile(old(a).folder, old(a).name);
    if ~old(a).isdir && (endsWith(fp, '.png') || endsWith(fp, '.svg') || endsWith(fp, '.fig'))
        delete(fp);
    end
end
save_some_figs_to_folder_2(out_dir, fig_tag, cf.Number, {'png', 'svg', 'fig'});
num = num2str(cf.Number);
movefile(fullfile(out_dir, [fig_tag '_figure_' num '.png']), fullfile(out_dir, [fig_tag '.png']));
movefile(fullfile(out_dir, [fig_tag '_figure_' num '.svg']), fullfile(out_dir, [fig_tag '.svg']));
movefile(fullfile(out_dir, [fig_tag '_f_' num '.fig']),      fullfile(out_dir, [fig_tag '.fig']));

% Log the git state alongside the figure so this presentation output can be
% traced back to an exact commit (+ working-tree diff if dirty).
capture_git_provenance(out_dir, project_root);

%% -------------------- Human-readable description --------------------
desc_path = fullfile(out_dir, 'README_fig_sfa_EOC.txt');
fid = fopen(desc_path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Stability_Manuscript figure: SFA Edge of Chaos (tau_a sensitivity)\n');
fprintf(fid, '=================================================================\n\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now')));
fprintf(fid, 'By script: %s.m\n\n', mfilename);

fprintf(fid, 'HOW IT WAS MADE\n');
fprintf(fid, '  Presentation replot -- no simulation is re-run. The script reloads the\n');
fprintf(fid, '  saved tau_a_E_max PSA object (psa_object.mat) and calls\n');
fprintf(fid, '  ParamSpaceAnalysis2.plot_sensitivity(''metric'',''LLE'') to rebuild the\n');
fprintf(fid, '  single-condition (SFA+STD) LLE-vs-tau_a panel, then restyles and saves\n');
fprintf(fid, '  it here. See git_provenance.txt for the exact commit.\n\n');

fprintf(fid, 'SOURCE RUN\n');
fprintf(fid, '  %s\n', data_root);
fprintf(fid, '  tau sensitivity subfolder used:\n');
fprintf(fid, '    %s\n\n', 'tau_sensitivity_tau_timescales_tau_a_E_max_nLevs_15_jul_07_26_08_49');

fprintf(fid, 'FIGURES PRODUCED (in this folder)\n');
fprintf(fid, '  Fig_SFA_EOC.png\n');
fprintf(fid, '  Fig_SFA_EOC.svg\n');
fprintf(fid, '  Fig_SFA_EOC.fig\n\n');

fprintf(fid, '  Single panel: largest Lyapunov exponent (growth rate) vs the maximum\n');
fprintf(fid, '  SFA adaptation timescale tau_a (the last, largest of 3 log-spaced\n');
fprintf(fid, '  tau_a_E elements, swept 5..60 s). As tau_a grows the growth rate rises\n');
fprintf(fid, '  toward 0 (the green dashed edge-of-chaos line) but stays negative.\n');
fprintf(fid, '  x-axis relabelled tau_a_E(last) -> "max \\tau_a (s)"; ylabel\n');
fprintf(fid, '  lambda_1 -> "Growth Rate"; condition title ("SFA + STD") removed.\n');
fprintf(fid, '  imagesc CLim capped at total_reps*%.2g; colormap white -> 90%% black so\n', clim_frac);
fprintf(fid, '  the blue median line stays visible. Blue median: alpha %.2g, width %d;\n', median_alpha, median_lw);
fprintf(fid, '  green zero line width %d; axis box removed. LLE range [%.2g, %.2g].\n', ...
    zeroline_lw, lle_range(1), lle_range(2));

clear cleanup;  % flush + close
fprintf('Description written: %s\n', desc_path);
