function names = mc_display_names(keys)
%MC_DISPLAY_NAMES Condition keys -> display strings, for the MC figures.
%
%   names = MC_DISPLAY_NAMES(keys)
%
% keys is the cell array a memory-capacity run saved in results_all.conditions.
% Since 2026-09-02 those are the project's snake_case KEYS (no_adaptation,
% sfa_only, std_only, sfa_and_std) rather than the MC figures' own display
% strings ('Baseline', 'SFA', 'STD', 'SFA+STD'), so an MC run directory names
% its conditions the same way every other run does and the display text lives
% in one place for all 19 figures.
%
% UNKNOWN KEYS PASS THROUGH UNCHANGED, deliberately. A saved run directory can
% name anything -- including the pre-2026-09-02 display strings, which are not
% keys in the map -- and a figure regenerated from an older MC run should still
% label its axes rather than erroring or printing blanks. This is the same
% isKey() guard srnn_condition_titles' own header asks callers to apply.
%
% It is a shared helper rather than a local function in each plotter because
% plot_memory_capacity and plot_memory_capacity_combined both need it, and a
% duplicated four-line map is exactly what srnn_condition_titles was created to
% eliminate.
%
% See also: srnn_condition_titles, plot_memory_capacity,
%           plot_memory_capacity_combined

titles = srnn_condition_titles();
names = cell(size(keys));
for k = 1:numel(keys)
    key = keys{k};
    if ischar(key) && isKey(titles, key)
        names{k} = titles(key);
    else
        names{k} = key;
    end
end
end
