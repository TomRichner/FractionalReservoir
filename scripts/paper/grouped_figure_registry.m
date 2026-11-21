function F = grouped_figure_registry(preset_name, source_fig_root, pytorch_file, human_psd_file)
% GROUPED_FIGURE_REGISTRY Append these entries after their source figures.
%
% The paper's grouped main figures 1 to 6 (fig_main1_grouped ..
% fig_main6_grouped), each in the group 'main' (Plan 2, package A,
% 2026-11-21; paper_figure_groups, paper_figure_entry). Groups 7 and 8 are
% NOT registered (FD-031 item 3, 2026-11-17): fig_main7_grouped and
% fig_main8_grouped are placeholders (a learning figure and a PSD sheet) that
% are not the paper's Figs 7 and 8, so a run's manifest would list figures
% the manuscript does not show. The paper's Fig 7 is fig_main_mc_grouped
% (mu7chain_gainbias_config) and its Fig 8 the combined stimulation figure
% (the IED stimulation bundle). The two functions stay in src/figures and
% can be drawn by hand.
% pytorch_file and human_psd_file were read by groups 7 and 8 only; they are
% still accepted and passed so callers need not change.
arguments
    preset_name char
    source_fig_root char
    pytorch_file char = ''
    human_psd_file char = ''
end
GROUPS = 1:6;
F = cell(1, numel(GROUPS));
for k = 1:numel(GROUPS)
    name = sprintf('fig_main%d_grouped', GROUPS(k));
    F{k} = paper_figure_entry(name, str2func(name), 'main', ...
        {'preset_name', preset_name, 'source_fig_root', source_fig_root, ...
        'pytorch_file', pytorch_file, 'human_psd_file', human_psd_file});
end
end
