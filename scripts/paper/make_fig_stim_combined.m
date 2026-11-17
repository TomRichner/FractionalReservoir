% MAKE_FIG_STIM_COMBINED The paper's stimulation figure (clinical over model).
% Machine-specific defaults live here, not in src/figures/fig_stim_combined.m.
% The clinical panels come from the cscs_stim_stability repo (read only; the
% approved figure is commit f0a5bab, Tom 2026-11-17). The output contains
% clinical recordings: never git add -f it.
%
% To swap in a newer model version, point model_fig at the .fig that
% fig_ied_stim writes for that run.
setup_paths();
root = fileparts(which('setup_paths'));
clinical_dir = 'C:/Users/m218089/Desktop/gitlab_repos/cscs_stim_stability/results/manuscript_stim/tica_awake/panels';
model_fig = 'C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/IED_model/stim_20261005/fig_ied_stim.fig';
out_dir = 'C:/Users/m218089/Desktop/github_repos/FractionalReservoir/figs/stim_combined_20261117';
out = fig_stim_combined(clinical_dir, model_fig, 'out_dir', out_dir, ...
    'clinical_note', 'cscs_stim_stability commit f0a5bab (k = 3 SD detection, 12 awake participants; Fig_Awake_Main_Normalized)', ...
    'model_note', 'IED model run stim_20261005, amplitude 0.20, 8 realizations (fig_ied_stim_stats.md beside it)');
fprintf('Wrote %s\n', strjoin(fullfile(out_dir, out.files), ', '));
