function out = fresp_fig_save(fig, cfg, out_dir, fig_tag, data_file, heading, caption, head, rows)
% FRESP_FIG_SAVE Save one frequency-response figure and its _table.md sidecar; return the figure's `out` struct.
%
%   out = FRESP_FIG_SAVE(fig, cfg, out_dir, fig_tag, data_file, heading, caption, head, rows)
%
%   fig        the figure
%   cfg        the figure function's options (.save, .verbose)
%   out_dir, fig_tag   where, and the file stem (save_figure_stable)
%   data_file  the .mat the figure read (out.source, and named in the table)
%   heading    the table file's title line
%   caption    what the panels show: everything that is not drawn on them
%   head, rows the table: 1 x m cellstr and an n x m cell of char or numbers
%              (numbers are written with three significant digits)
%
% out.figs, out.files (the images written and <fig_tag>_table.md),
% out.source. Nothing is written with cfg.save false.
%
% See also: save_figure_stable, existing_outputs, fresp_fig_load

out = struct('figs', fig, 'files', {{}}, 'source', data_file);
lines = cell(1, size(rows, 1));
for r = 1:size(rows, 1)
    cells = cell(1, size(rows, 2));
    for q = 1:size(rows, 2)
        v = rows{r, q};
        if ischar(v)
            cells{q} = v;
        elseif isscalar(v) && isreal(v)
            cells{q} = sprintf('%.3g', v);
        else
            cells{q} = mat2str(v, 3);
        end
    end
    lines{r} = sprintf('| %s |', strjoin(cells, ' | '));
end
if cfg.save
    save_figure_stable(out_dir, fig_tag, fig);
    out.files = existing_outputs(out_dir, fig_tag);
    table_name = [fig_tag '_table.md'];
    fid = fopen(fullfile(out_dir, table_name), 'w');
    if fid > 0
        closer = onCleanup(@() fclose(fid));
        fprintf(fid, '# %s\n\n%s\n\nData: `%s` (stage `frequency_response`).\n\n', heading, caption, data_file);
        fprintf(fid, '| %s |\n', strjoin(head, ' | '));
        fprintf(fid, '|%s\n', repmat('---|', 1, numel(head)));
        fprintf(fid, '%s\n', lines{:});
        out.files{end + 1} = table_name;
    end
end
vprintf(cfg.verbose, 'verbose', '%s\n', lines{:});
end
