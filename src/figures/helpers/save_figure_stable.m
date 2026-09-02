function save_figure_stable(out_dir, fig_tag, fh)
% SAVE_FIGURE_STABLE Save one figure as <fig_tag>.{png,svg,fig} in out_dir.
%
%   SAVE_FIGURE_STABLE(out_dir, fig_tag, fh)
%
% save_some_figs_to_folder_2 suffixes filenames with the (run-dependent) figure
% number; save, then rename to fixed names so re-runs overwrite cleanly. First
% clear any prior outputs for THIS tag (stable or numbered) so nothing stale
% lingers. Saving by fh.Number means other open figures are untouched.
%
% The renames are GUARDED because save_some_figs_to_folder_2 catches export
% failures, warns, and carries on (see its comment: roughly one run in four, a
% live figure reaches a state where every rasterizing path throws). An
% unguarded movefile would then abort the calling script and lose the figures
% that did export.
%
% Extracted verbatim from the local subfunction of
% Fig_sensitivity_analysis_allStd.m so fig_sensitivity_medians can share it.
%
% See also: save_some_figs_to_folder_2
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
            warning('save_figure_stable:ExportMissing', ...
                'Export did not produce %s; skipping that format for ''%s''.', ...
                renames{r, 1}, fig_tag);
        end
    end
end
