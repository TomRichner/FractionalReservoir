function e = paper_figure_entry(name, fn, group, args)
% PAPER_FIGURE_ENTRY One row of a paper bundle's figure registry, with its group (Plan 2, package A).
%
%   e = PAPER_FIGURE_ENTRY(name, fn, group)
%   e = PAPER_FIGURE_ENTRY(name, fn, group, args)
%
%   name   the entry's name: its folder under the figure root, unique in the
%          bundle (char)
%   fn     the figure function (handle), called as fn(standard args, args{:})
%   group  'main' | 'supplement' | 'appendix' (paper_figure_groups)
%   args   the entry's own name-value arguments (cell, default {})
%
%   e      struct with fields, in this order: name, fn, group, in_paper, args.
%          in_paper is DERIVED from the group (true for main and supplement)
%          and kept so that existing callers that read it (the manifest, the
%          bundle tests, ied_stim_figures) read what they did.
%
% Every entry of every paper bundle (mu7chain_config, mu7chain_AC_config,
% mu7chain_gainbias_config, ied_stim_config, and the helpers that add to
% them: grouped_figure_registry, with_spectral_supplement,
% with_switch_spectra, with_ritz_overlay) is built here, so that the field
% order is one (struct arrays of entries can be concatenated) and an unknown
% group is an error where the row is written. The regression bundle and the
% older bundles keep their plain structs without a group.
%
% Errors :BadGroup for a group paper_figure_groups does not name.
%
% See also: paper_figure_groups, make_all_paper_figures
arguments
    name (1,:) char
    fn (1,1) function_handle
    group (1,:) char
    args cell = {}
end
G = paper_figure_groups();
k = find(strcmp({G.name}, group), 1);
if isempty(k)
    error('paper_figure_entry:BadGroup', 'Entry %s: group ''%s'' is not one of %s.', ...
        name, group, strjoin({G.name}, ', '));
end
e = struct('name', name, 'fn', fn, 'group', group, 'in_paper', G(k).in_paper, 'args', {args});
end
