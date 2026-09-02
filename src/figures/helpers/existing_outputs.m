function files = existing_outputs(out_dir, fig_tag)
% EXISTING_OUTPUTS The figure files actually on disk for one tag.
%
%   files = EXISTING_OUTPUTS(out_dir, 'Fig_Sensitivity_LLE_core')
%
% Returns a cellstr of file NAMES (not paths) matching <fig_tag>.{png,svg,fig}.
%
% Called AFTER save_figure_stable so the caller can report what was really
% written, and so the manifest lists the files that exist rather than the
% ones that were intended. That distinction matters here: roughly one run in
% four, save_some_figs_to_folder_2 hits a figure state where a rasterizing path
% throws; it warns and carries on, so a run can legitimately be missing one
% format. A README claiming all three would be wrong.
%
% The master script uses a non-empty return as its per-figure success check.

arguments
    out_dir (1,:) char
    fig_tag (1,:) char
end

files = {};
for ext = {'png', 'svg', 'fig'}
    name = sprintf('%s.%s', fig_tag, ext{1});
    if isfile(fullfile(out_dir, name))
        files{end+1} = name; %#ok<AGROW>
    end
end
end
