function replot_dirs = replot_all_analyses(opts)
% REPLOT_ALL_ANALYSES Regenerate the raw sweep figures for a finished run.
%
%   replot_dirs = REPLOT_ALL_ANALYSES()
%   replot_dirs = REPLOT_ALL_ANALYSES('run_dir', d)
%
% Reruns the three replot helpers -- tau sensitivity, 1-D sensitivity (plus the
% stacked sheets), and parameter space -- against a prior run_all_<dt>
% directory, without re-running any simulation. Each writes into its own
% replot_<...>_<dt>/ subfolder, so originals are never overwritten.
%
% THIS IS THE DIAGNOSTIC VIEW, not the manuscript. It regenerates the plain PSA
% plots, which are what to look at when checking a sweep. The publication
% figures are make_all_paper_figures, which restyles these into the manuscript's
% sheets and writes them into the figure folders.
%
% WAS A SCRIPT with a hardcoded data_root pointing at run_all_mar_02_26_17_12 --
% a directory that no longer exists, so it could not run at all. It now resolves
% the run the same way every figure does.
%
% See also: make_all_paper_figures, resolve_run_dir, replot_sensitivity,
%           replot_tau_sensitivity, replot_param_space_analysis

arguments
    opts.run_dir        (1,:) char   = ''
    opts.preset_name    (1,:) char   = 'celltype_pairs_Sc0p2_noise0p025_dualStd_4cond'
    opts.lle_hist_range (1,2) double = [-1, 1]
    opts.assemble       (1,1) logical = true
end

setup_paths();

run_dir = resolve_run_dir('run_dir', opts.run_dir, 'preset_name', opts.preset_name);

fprintf('=== Replotting all analyses ===\n');
fprintf('  run_dir        : %s\n', run_dir);
fprintf('  LLE hist_range : [%g, %g]\n\n', opts.lle_hist_range(1), opts.lle_hist_range(2));
t0 = tic;

replot_dirs = struct();

fprintf('[1/3] tau sensitivity...\n');
replot_dirs.tau = replot_tau_sensitivity(run_dir, opts.lle_hist_range);

fprintf('[2/3] 1-D sensitivity...\n');
replot_dirs.sensitivity = replot_sensitivity(run_dir, opts.lle_hist_range);
if opts.assemble
    % Two sheets, split by what the parameters mean. Skipped with a warning on
    % a model class that has no mu blocks, so this costs nothing on SRNNModel2.
    [~, model_class] = srnn_param_preset(opts.preset_name);
    if strcmp(model_class, 'SRNNCellTypePairs'); f_param = 'f_E'; else; f_param = 'f'; end
    assemble_sensitivity_figure(replot_dirs.sensitivity, 'LLE', ...
        {f_param, 'n', 'level_of_chaos'}, 'network');
    assemble_sensitivity_figure(replot_dirs.sensitivity, 'LLE', ...
        {'mu_EE_relative', 'mu_EI_relative', 'mu_IE_relative', 'mu_II_relative'}, ...
        'mu_blocks');
end

fprintf('[3/3] parameter space...\n');
replot_dirs.param_space = replot_param_space_analysis(run_dir);

fprintf('\n=== Replot complete in %.2f min ===\n', toc(t0)/60);
end
