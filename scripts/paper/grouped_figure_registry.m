function F = grouped_figure_registry(preset_name, source_fig_root, pytorch_file, human_psd_file)
% GROUPED_FIGURE_REGISTRY The composite figures of block A, in their groups: main Figures 1 to 4, supplementary S1 to S5, and the old composites in the appendix.
%
%   F = GROUPED_FIGURE_REGISTRY(preset_name, source_fig_root)
%   F = GROUPED_FIGURE_REGISTRY(preset_name, source_fig_root, pytorch_file, human_psd_file)
%
% mu7chain_config appends these after its single figures (some composites
% copy a single figure's saved output, so they come after their sources).
% The groups are the figure plan's (the manuscript repository's
% docs/Advancing_the_narrative/figures_main_supp_appen.md; Plan 2,
% integration, 2026-11-21). Where the other composites are registered:
%   with_spectral_supplement   main Figures 5 and 6, S6, S7, S8 (they read
%                              the spectral stages)
%   with_switch_spectra        rewrites fig_main3_grouped into the paper's
%                              switching Figure 3 and adds its appendix
%                              versions (fig_main3_grouped_ritz, _baseline)
%   mu7chain_AC_config         main Figure 7 and S9 (they read blocks A and C)
%   ied_stim_config            main Figure 8 and S10 (the stimulation bundle)
%
% THE ENTRIES, in order:
%   fig_main1_grouped   main   fig_main1_concept_networks: four example
%                              networks (stage concept_examples) above the
%                              sigmoid and circle diagrams
%   fig_main2_grouped   main   fig_grouped_main group 2 (fig_main2_grouped)
%   fig_main3_grouped   main   fig_grouped_main group 3 (fig_main3_grouped);
%                              with_switch_spectra makes it the switching
%                              composite with row D and no Ritz overlay
%   fig_main4_grouped   main   fig_main4_connectivity: the sweeps, the joint
%                              sample, and panel C (leading direction against
%                              lambda_1 with the uncoupled reference)
%   fig_supp1_grouped .. fig_supp5_grouped   supplement   S1 to S5
%                              (fig_supp3_grouped is
%                              fig_supp3_connectivity_gallery)
%   the OLD COMPOSITES that left the main group, appendix, each drawn by
%   fig_grouped_main_appendix (fig_grouped_main's groups 1, 5, 6, 7, 8) under
%   a tag of its own, so no file of theirs carries a paper figure's name:
%     fig_main1_concepts_old      the old Figure 1 (concept panels)
%     fig_main5_transient_old     the old Figure 5 (margin and directional gains)
%     fig_main6_memory_old        the old Figure 6 (timescales and memory)
%     fig_main7_placeholder       the learning placeholder (FD-031 item 3)
%     fig_main8_psd_placeholder   the PSD sheet placeholder (FD-031 item 3)
%
% pytorch_file and human_psd_file are read by the two placeholders only.
%
% See also: mu7chain_config, paper_figure_entry, fig_grouped_main,
%           fig_grouped_main_appendix, with_spectral_supplement, with_switch_spectra
arguments
    preset_name char
    source_fig_root char
    pytorch_file char = ''
    human_psd_file char = ''
end
P = {'preset_name', preset_name};
old = {'preset_name', preset_name, 'source_fig_root', source_fig_root};
F = { ...
    paper_figure_entry('fig_main1_grouped', @fig_main1_concept_networks, 'main', P), ...
    paper_figure_entry('fig_main2_grouped', @fig_main2_grouped, 'main', old), ...
    paper_figure_entry('fig_main3_grouped', @fig_main3_grouped, 'main', old), ...
    paper_figure_entry('fig_main4_grouped', @fig_main4_connectivity, 'main', [P, {'fig_tag', 'Fig_Main4_Grouped'}]), ...
    paper_figure_entry('fig_supp1_grouped', @fig_supp1_grouped, 'supplement', P), ...
    paper_figure_entry('fig_supp2_grouped', @fig_supp2_grouped, 'supplement', P), ...
    paper_figure_entry('fig_supp3_grouped', @fig_supp3_connectivity_gallery, 'supplement', P), ...
    paper_figure_entry('fig_supp4_grouped', @fig_supp4_grouped, 'supplement', P), ...
    paper_figure_entry('fig_supp5_grouped', @fig_supp5_grouped, 'supplement', P), ...
    legacy('fig_main1_concepts_old',    1, 'Fig_Main1_Concepts_old', old), ...
    legacy('fig_main5_transient_old',   5, 'Fig_Main5_Transient_old', old), ...
    legacy('fig_main6_memory_old',      6, 'Fig_Main6_Memory_old', old), ...
    legacy('fig_main7_placeholder',     7, 'Fig_Main7_Placeholder', [old, {'pytorch_file', pytorch_file}]), ...
    legacy('fig_main8_psd_placeholder', 8, 'Fig_Main8_PSD_Placeholder', [old, {'human_psd_file', human_psd_file}])};
end

function e = legacy(name, group, tag, args)
e = paper_figure_entry(name, @fig_grouped_main_appendix, 'appendix', [{'group', group}, args, {'fig_tag', tag}]);
end
