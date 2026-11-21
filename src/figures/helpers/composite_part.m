function [fig, table_text, out] = composite_part(fn, cfg, varargin)
% COMPOSITE_PART Call one figure function for a composite: its figure, invisible and unsaved, and the text of its _table.md.
%
%   [fig, table_text, out] = COMPOSITE_PART(@fig_outlier_dose, cfg)
%   [fig, table_text, out] = COMPOSITE_PART(@fig_eigdens_radial, cfg, 'network', 'EEx1p5', 'panels', {'cumulative'})
%
%   fn    a figure function with the standard signature and the option
%         `images` (false with save: the sidecar alone is written)
%   cfg   the composite's options: run_dir, conditions, save, verbose are
%         passed on; further name-value pairs follow
%
% With cfg.save the source function is asked for its sidecar only, into a
% temporary folder that is removed again, and the text is returned for
% composite_save; without, table_text is ''. The figure is the source's
% own, to be handed to composite_figure (which closes it).
%
% See also: composite_figure, composite_save

args = [{'run_dir', cfg.run_dir, 'conditions', cfg.conditions, 'verbose', cfg.verbose, 'visible', false}, varargin];
table_text = '';
if cfg.save
    tmp = tempname; mkdir(tmp);
    cleaner = onCleanup(@() rmdir(tmp, 's'));
    out = fn(args{:}, 'save', true, 'images', false, 'out_dir', tmp);
    d = dir(fullfile(tmp, '*_table.md'));
    if ~isempty(d); table_text = fileread(fullfile(tmp, d(1).name)); end
else
    out = fn(args{:}, 'save', false);
end
fig = out.figs(1);
end
