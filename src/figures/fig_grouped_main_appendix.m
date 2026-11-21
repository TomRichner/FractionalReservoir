function out = fig_grouped_main_appendix(varargin)
% FIG_GROUPED_MAIN_APPENDIX The old composites of fig_grouped_main that left the main group, as appendix entries.
%
%   out = FIG_GROUPED_MAIN_APPENDIX('group', k, <fig_grouped_main's name-value options>)
%
% k is 1, 5, 6, 7 or 8: the old Figure 1 (concept panels), the old
% Figure 5 (margin and directional gains), the old Figure 6 (timescales and
% memory), and the two placeholders (learning, PSD sheet). Everything else
% is passed to fig_grouped_main(k, ...) unchanged. grouped_figure_registry
% registers each in the group 'appendix' with a fig_tag of its own.
%
% WHY THIS WRAPPER (Plan 2, integration, 2026-11-21). Until then each group
% had a wrapper fig_main<k>_grouped.m, and the registry entry of that name
% called it. Since the figure plan, the entries fig_main1_grouped and
% fig_main4_grouped to fig_main7_grouped draw new composites
% (fig_main1_concept_networks, fig_main4_connectivity, fig_main5_mechanism,
% fig_main6_amplification, fig_main7_memory), so a function
% fig_main5_grouped drawing the OLD Figure 5 would be a second, misleading
% path. The wrappers of groups 1 and 4 to 8 were removed; groups 2 and 3
% keep theirs (they are still the paper's Figures 2 and 3), group 4 is
% reached through fig_main4_connectivity (its panels A and B), and the old
% groups 1, 5, 6, 7 and 8 through this function.
%
% See also: fig_grouped_main, grouped_figure_registry
k = find(strcmp(varargin(1:2:end), 'group'), 1);
if isempty(k)
    error('fig_grouped_main_appendix:NoGroup', 'Name the old composite with ''group'', k (1, 5, 6, 7 or 8).');
end
group = varargin{2 * k};
varargin(2 * k - 1:2 * k) = [];
if ~(isnumeric(group) && isscalar(group) && ismember(group, [1 5 6 7 8]))
    error('fig_grouped_main_appendix:BadGroup', 'group must be 1, 5, 6, 7 or 8 (the old composites).');
end
out = fig_grouped_main(group, varargin{:});
end
