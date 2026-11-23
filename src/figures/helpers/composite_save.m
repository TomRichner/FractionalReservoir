function out = composite_save(fig, cfg, out_dir, fig_tag, heading, caption, parts, sources)
% COMPOSITE_SAVE Save a composite figure and its _table.md: a caption paragraph, then the sidecar of every figure it is built from.
%
%   out = COMPOSITE_SAVE(fig, cfg, out_dir, fig_tag, heading, caption, parts, sources)
%
%   cfg       the composite's options (.save)
%   heading   the table file's title line
%   caption   what the panels are: the paragraph a caption is written from
%   parts     n x 2 cell: {'Panel A', table text of the source figure}
%             (composite_part reads it), written as sections in this order
%   sources   cellstr, the data files read (out.source)
%
% out.figs, out.files (images written and <fig_tag>_table.md), out.source.
% Nothing is written with cfg.save false.
%
% See also: composite_figure, composite_part, save_figure_stable

out = struct('figs', fig, 'files', {{}}, 'source', {sources});
if ~cfg.save; return; end
sz = getappdata(fig, 'composite_size');
if ~isempty(sz) && any(abs(fig.Position(3:4) - sz) > 1e-6); fig.Position(3:4) = sz; drawnow; end
true_minus_ticks(fig);   % the true minus on every numeric tick label, last before export
save_figure_stable(out_dir, fig_tag, fig, canvas_export_args());   % the whole canvas, untrimmed
out.files = existing_outputs(out_dir, fig_tag);
name = [fig_tag '_table.md'];
fid = fopen(fullfile(out_dir, name), 'w');
if fid < 0; return; end
closer = onCleanup(@() fclose(fid));
fprintf(fid, '# %s\n\n%s\n\n', heading, caption);
fprintf(fid, 'Data read:\n\n');
for k = 1:numel(sources); fprintf(fid, '- `%s`\n', sources{k}); end
fprintf(fid, '\n');
for k = 1:size(parts, 1)
    txt = regexprep(parts{k, 2}, '^# ', '', 'once');             % the source's own title line becomes part of the section's
    fprintf(fid, '## %s: %s\n\n', parts{k, 1}, strtrim(txt));
    fprintf(fid, '\n');
end
out.files{end + 1} = name;
end
