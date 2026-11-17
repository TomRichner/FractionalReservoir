% MAKE_FIG_STIM_COMBINED_V2 The paper's stimulation figure in the mirror
% layout (Tom, 2026-11-17): clinical A traces, B rate, C PSD; model D traces,
% E rate, F PSD. The output contains clinical recordings: never git add -f it.
%
%   stage = 1  existing model figure (stim_20261005, events/min panel)
%   stage = 2  block run stim_blocks_20261117 (rate normalised as panel B),
%              model figure written by fig_ied_stim_blocks
if ~exist('stage','var'), stage = 2; end
setup_paths();
clinical_dir = 'C:/Users/m218089/Desktop/gitlab_repos/cscs_stim_stability/results/manuscript_stim/tica_awake/panels';
fr = 'C:/Users/m218089/Desktop/github_repos/FractionalReservoir';
clinical_note = 'cscs_stim_stability commit f0a5bab (k = 3 SD detection, 12 awake participants; Fig_Awake_Main_Normalized)';
if stage == 1
    model_fig = [fr '/figs/IED_model/stim_20261005/fig_ied_stim.fig'];
    out_dir = [fr '/figs/stim_combined_v2_stage1'];
    model_note = 'IED model run stim_20261005, amplitude 0.20, 8 realizations (fig_ied_stim_stats.md beside it); stage 1, existing panels only';
else
    model_fig = [fr '/figs/IED_model/stim_blocks_20261117/fig_ied_stim_blocks.fig'];
    out_dir = [fr '/figs/stim_combined_v2'];
    model_note = 'IED model block run stim_blocks_20261117, amplitude 0.20, 12 networks, 4 + 4 blocks each (fig_ied_stim_blocks_stats.md beside it)';
end
out = fig_stim_combined(clinical_dir, model_fig, 'out_dir', out_dir, 'layout', 'mirror', ...
    'clinical_note', clinical_note, 'model_note', model_note);
fprintf('Wrote %s\n', strjoin(fullfile(out_dir, out.files), ', '));
