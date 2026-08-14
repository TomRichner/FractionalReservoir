close all
clc

% Fig_sfa_EOC_allStd.m
% Stability_Manuscript presentation figure: SFA edge-of-chaos, for the
% SRNNCellTypePairs "allStd" runs. Replots the tau_a-sensitivity LLE panel -- how
% the largest Lyapunov exponent (lambda_1) approaches 0 (edge of chaos) as the
% maximum SFA adaptation timescale tau_a grows. No simulation is re-run: it
% reloads the saved tau_a_E_max PSA object and restyles the single-condition
% (SFA+STD) panel, matching the presentation treatment in
% fig_sensitivity_analysis_allStd.
%
% Sibling of ../fig_sfa_EOC/Fig_sfa_EOC.m, which is left untouched. Two things
% differ:
%
%   1. THE LOADER. The original did `S = load(psa_object.mat); psa = S.psa_tau_a`.
%      That variable name no longer exists: run() now writes the object itself,
%      always as `psa`, and ParamSpaceAnalysis2.from_dir is the supported way in
%      (it selects the saved object by CLASS, so it reads old psa_tau_a files
%      too). The original would fail outright on this run.
%   2. The tau subfolder is located by GLOB rather than by its timestamped name,
%      so re-pointing data_root at the medium run needs no other edit.
%
% The x-axis follows the sweep automatically -- tau_a_E(end) now runs over
% [1, 30] s rather than the original [5, 60] s.
%
% See also: Fig_sensitivity_analysis_allStd, Fig_param_space_allStd

this_dir = fileparts(mfilename('fullpath'));
% Depth-independent project root (CLAUDE.md convention).
setup_paths();
project_root = fileparts(which('setup_paths'));

% Source run. Swap this one line to regenerate against the medium run.
data_root = fullfile(project_root, 'data', 'param_space', 'run_all_aug_14_26_12_04');
out_dir   = this_dir;   % write the final figure next to this script

% --- Locate the tau sensitivity subfolder ----------------------------------
tau_listing = dir(fullfile(data_root, 'tau_sensitivity_*'));
tau_listing = tau_listing([tau_listing.isdir]);
if isempty(tau_listing)
    error('Fig_sfa_EOC_allStd:NoTauDir', ...
        'No tau_sensitivity_* subfolder found in:\n  %s', data_root);
end
if numel(tau_listing) > 1
    warning('Fig_sfa_EOC_allStd:MultipleTauDirs', ...
        'Found %d tau_sensitivity_* subfolders; using the newest.', numel(tau_listing));
    [~, newest] = max([tau_listing.datenum]);
    tau_listing = tau_listing(newest);
end
tau_dir = fullfile(tau_listing.folder, tau_listing.name);

% --- Presentation constants ------------------------------------------------
% Kept identical to the original figure by decision, so the two read the same.
% Worth knowing when reading the output: on this preset the tau LLEs run
% -0.23 .. +0.26 with a median of +0.002, so roughly half the distribution is
% POSITIVE. hist_range tops out at 0.1 and the view is cropped at 0.05, which
% puts the +inf overflow band off-screen -- the panel will look entirely
% sub-zero when it is not. Widen lle_range / y_view together to change that.
lle_range = [-0.3, 0.1];    % histogram range passed to plot_sensitivity
y_view    = [-0.25, 0.05];  % visible y-range (crops the overflow bands)
y_ticks   = -0.2:0.1:0.1;
fig_position = [457 637 264 252];

tick_fs   = 14;    % tick numbers
label_fs  = 15.4;  % axis labels
clim_frac = 0.4;   % darken imagesc: cap CLim at total_reps*clim_frac
% Colormap ramps white (0 counts) -> 90% black (max), not pure black, so the
% blue median line stays visible over the darkest cells.
dark_cmap    = repmat(linspace(1, 0.1, 256)', 1, 3);
median_alpha = 0.35;   % blue median line transparency
median_lw    = 3;      % blue median line width
zeroline_lw  = 2;      % green dashed zero line width

% Start clean: plot_sensitivity operates on the current session's figures.
close all force;

% --- Reload the tau PSA and regenerate the single LLE panel ----------------
psa = ParamSpaceAnalysis2.from_dir(tau_dir);

psa.plot_sensitivity('metric', 'LLE', 'hist_range', lle_range);
cf = gcf;

% --- Presentation restyle (same knobs as Fig_sensitivity_analysis_allStd) ---
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

% Labels: keep lambda_1; x -> "max \tau_a (s)"; drop the title. latex here (not
% tex) so the symbol matches the latex xlabel below.
ylabel(ax, '$\lambda_1$', 'Interpreter', 'latex', 'FontSize', label_fs);
xlabel(ax, 'max $\tau_a$ (s)', 'Interpreter', 'latex', 'FontSize', label_fs);
title(ax, '');

% Fewer y-ticks + a tighter view around the data.
ylim(ax, y_view);
yticks(ax, y_ticks);

set(cf, 'Position', fig_position);

% --- Save ONLY the cleaned panel, with a STABLE name -----------------------
fig_tag = 'Fig_SFA_EOC_allStd';
save_figure_stable(out_dir, fig_tag, cf);

% Log the git state alongside the figure so this presentation output can be
% traced back to an exact commit (+ working-tree diff if dirty).
capture_git_provenance(out_dir, project_root);

%% -------------------- Human-readable description --------------------
desc_path = fullfile(out_dir, 'README_fig_sfa_EOC_allStd.txt');
fid = fopen(desc_path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Stability_Manuscript figure: SFA Edge of Chaos (tau_a sensitivity, allStd)\n');
fprintf(fid, '=========================================================================\n\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now')));
fprintf(fid, 'By script: %s.m\n\n', mfilename);

fprintf(fid, 'HOW IT WAS MADE\n');
fprintf(fid, '  Presentation replot -- no simulation is re-run. The script reloads the\n');
fprintf(fid, '  saved tau_a_E_max PSA object via ParamSpaceAnalysis2.from_dir and calls\n');
fprintf(fid, '  plot_sensitivity(''metric'',''LLE'') to rebuild the single-condition\n');
fprintf(fid, '  (SFA+STD) LLE-vs-tau_a panel, then restyles and saves it here. See\n');
fprintf(fid, '  git_provenance.txt for the exact commit.\n\n');

fprintf(fid, 'SOURCE RUN\n');
fprintf(fid, '  %s\n', data_root);
fprintf(fid, '  tau sensitivity subfolder used:\n');
fprintf(fid, '    %s\n\n', tau_listing.name);

fprintf(fid, 'FIGURES PRODUCED (in this folder)\n');
fprintf(fid, '  %s.png\n', fig_tag);
fprintf(fid, '  %s.svg\n', fig_tag);
fprintf(fid, '  %s.fig\n\n', fig_tag);

fprintf(fid, '  Single panel: largest Lyapunov exponent (lambda_1) vs the maximum\n');
fprintf(fid, '  SFA adaptation timescale tau_a (the last, largest of 3 log-spaced\n');
fprintf(fid, '  tau_a_E elements). x-axis relabelled tau_a_E(last) -> "max \\tau_a (s)";\n');
fprintf(fid, '  ylabel kept as \\lambda_1 (latex, matching the xlabel); condition title\n');
fprintf(fid, '  ("SFA + STD") removed.\n');
fprintf(fid, '  imagesc CLim capped at total_reps*%.2g; colormap white -> 90%% black so\n', clim_frac);
fprintf(fid, '  the blue median line stays visible. Blue median: alpha %.2g, width %d;\n', median_alpha, median_lw);
fprintf(fid, '  green zero line width %d; axis box removed. LLE histogram range\n', zeroline_lw);
fprintf(fid, '  [%.2g, %.2g], view cropped to [%.2g, %.2g].\n\n', ...
    lle_range(1), lle_range(2), y_view(1), y_view(2));

fprintf(fid, 'READING THIS PANEL -- IMPORTANT\n');
fprintf(fid, '  The histogram range and y-view are kept identical to the original\n');
fprintf(fid, '  figure by choice. On this preset that crops real data: the tau LLEs\n');
fprintf(fid, '  span about -0.23 to +0.26 with a MEDIAN of +0.002, so roughly half the\n');
fprintf(fid, '  distribution is positive, while the view stops at %.2g and the +inf\n', y_view(2));
fprintf(fid, '  overflow band sits off-screen. The panel therefore looks entirely\n');
fprintf(fid, '  sub-zero when it is not. To show the full spread, widen lle_range and\n');
fprintf(fid, '  y_view (and y_ticks) at the top of the script.\n');

clear cleanup;  % flush + close
fprintf('Description written: %s\n', desc_path);


%% ==================== Local helper ====================
function save_figure_stable(out_dir, fig_tag, fh)
% SAVE_FIGURE_STABLE Save one figure as <fig_tag>.{png,svg,fig} in out_dir.
%
% save_some_figs_to_folder_2 suffixes filenames with the (run-dependent) figure
% number; save, then rename to fixed names so re-runs overwrite cleanly. First
% clear any prior outputs for THIS tag (stable or numbered) so nothing stale
% lingers. Saving by fh.Number means other open figures are untouched.
%
% The renames are GUARDED because save_some_figs_to_folder_2 catches export
% failures, warns, and carries on (see its comment: roughly one run in four, a
% live figure reaches a state where every rasterizing path throws). An
% unguarded movefile would then abort this script.
    old = dir(fullfile(out_dir, [fig_tag '*']));
    for a = 1:numel(old)
        fp = fullfile(old(a).folder, old(a).name);
        if ~old(a).isdir && (endsWith(fp, '.png') || endsWith(fp, '.svg') || endsWith(fp, '.fig'))
            delete(fp);
        end
    end

    save_some_figs_to_folder_2(out_dir, fig_tag, fh.Number, {'png', 'svg', 'fig'});

    num = num2str(fh.Number);
    renames = { ...
        [fig_tag '_figure_' num '.png'], [fig_tag '.png']; ...
        [fig_tag '_figure_' num '.svg'], [fig_tag '.svg']; ...
        [fig_tag '_f_' num '.fig'],      [fig_tag '.fig']};
    for r = 1:size(renames, 1)
        src = fullfile(out_dir, renames{r, 1});
        if isfile(src)
            movefile(src, fullfile(out_dir, renames{r, 2}));
        else
            warning('Fig_sfa_EOC_allStd:ExportMissing', ...
                'Export did not produce %s; skipping that format for ''%s''.', ...
                renames{r, 1}, fig_tag);
        end
    end
end
