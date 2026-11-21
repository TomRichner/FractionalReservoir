function G = paper_figure_groups()
% PAPER_FIGURE_GROUPS The three groups of the paper's figures and their limits, in one place (Plan 2, package A).
%
%   G = PAPER_FIGURE_GROUPS()
%
%   G  1 x 3 struct array, in the order the groups are reported:
%        name      'main' | 'supplement' | 'appendix'
%        title     the heading used in figure_groups.md
%        limit     the most figures the paper may have in the group, counted
%                  over the bundles that make the paper (Inf: no limit)
%        in_paper  true for the groups the manuscript references (main and
%                  supplement); a registry entry's in_paper is derived from
%                  this, never set by hand
%
% THE RULES (Tom, 2026-11-21; the manuscript repository's
% docs/Advancing_the_narrative/figures_main_supp_appen.md, "The rules"):
% main, eight figures and no more; supplement, about ten; appendix, the rest.
% Every main and supplementary figure is referenced in the manuscript text;
% an appendix figure never is, but is made in the same style, so that it can
% be promoted or demoted without rework. "About ten" is held here as a hard
% limit of ten: raise it here, with Tom's word, if the plan grows.
%
% Who reads it: paper_figure_entry (validates a row's group and derives
% in_paper), make_all_paper_figures (figure_groups.md and the summary line's
% counts by group), and test_paper_figure_groups (the guard on the bundles).
% Nothing else names the groups or their limits.
%
% See also: paper_figure_entry, make_all_paper_figures, grouped_figure_registry,
%           test_paper_figure_groups
G = struct( ...
    'name',     {'main', 'supplement', 'appendix'}, ...
    'title',    {'Main figures', 'Supplementary figures', 'Appendix figures'}, ...
    'limit',    {8, 10, Inf}, ...
    'in_paper', {true, true, false});
end
