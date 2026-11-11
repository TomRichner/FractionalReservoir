function cfg = regression_fast_config(tag)
% REGRESSION_FAST_CONFIG The regression bundle: the manuscript's pipeline on a small network, in minutes (2026-11-10, FD-002).
%
%   cfg = regression_fast_config()        % run_dir data/regression_fast/reference
%   cfg = regression_fast_config('abc')   % run_dir data/regression_fast/abc
%
% WHAT IT IS FOR. A change to the model class can alter results it was not
% meant to touch. This bundle runs every analysis stage and every figure the
% manuscript's bundle (mu7revisedMedium_config) runs, on a network one fifth
% the size at run_mode 'fast', so the whole pipeline can be rerun before and
% after a change and the outputs compared (regression_digest,
% test_regression_bundle). It is a test fixture: nothing it produces goes in
% the paper.
%
% HOW IT RELATES TO THE MANUSCRIPT BUNDLE. This file is mu7revisedMedium_config
% with four things changed and nothing else:
%   * the four presets are the regression ones (regression_pairs_* and
%     single_neuron_regression: the manuscript network at n = 100,
%     indegree = 20; see their case blocks in srnn_param_preset);
%   * run_mode is 'fast';
%   * run_dir and fig_root are under data/regression_fast/ and
%     figs/regression_fast/, in a subfolder named by TAG;
%   * the figure registry is the same list.
% Every other setting below is the manuscript bundle's, so the same code
% paths run. If a stage or figure is added to the manuscript bundle, add it
% here too.
%
% TAG names the run. run_all_paper_analyses refuses a non-empty run
% directory, so each run needs its own: 'reference' (the default) is the run
% the saved digest was made from; a comparison run uses another tag.
%
% See also: regression_fast_run, mu7revisedMedium_config, srnn_param_preset

arguments
    tag (1,:) char = 'reference'
end
if isempty(regexp(tag, '^[A-Za-z0-9_\-]+$', 'once'))
    error('regression_fast_config:BadTag', ...
        'tag must be letters, digits, underscore or hyphen (got ''%s'').', tag);
end

cfg = struct();

%% The experiment
% The manuscript's 3-condition network (mu_tilde_relative 7 on every route,
% SFA on both cell types, per-neuron setpoints and SFA ladders, no external
% input, Wiener noise) at n = 100, indegree = 20. Regimes: no_adaptation /
% sfa1_std1 / sfa3_std2, depression as the product of independent components.
% See the preset, which is written out in full; its three companions chain
% to it.
%
% 'fast': 5 levels x 5 reps per 1-D sweep, 7 tau levels, 128 joint samples,
% 5 memory-capacity trials, 10-s sensitivity runs and 20-s joint runs.
cfg.preset_name = 'regression_pairs_noStim_noise0p025_dualStd_3cond';
run_mode = 'fast';   % FIXED: the regression bundle is the fast run
cfg.run_mode    = run_mode;

%% Where things land -- both fixed, so this cannot touch any other run
cfg.run_dir  = ['data/regression_fast/' tag];   % analyses write, figures read
cfg.fig_root = ['figs/regression_fast/' tag];   % overwritten in place

% Figures are built and saved but do not pop up: a new figure window raises
% itself and takes keyboard focus, and a run draws dozens.
cfg.visible_figures = false;
cfg.verbose = 'minimal';   % 'verbose' | 'minimal' | 'near-none' -- how much the run prints (verbose_level)

%% Presets per stage and figure
% Memory capacity on THIS network with the Wiener process OFF (TR 2026-09-14):
% the preset below chains to cfg.preset_name and sets sigma_u_noise = 0 (the
% within-bundle chaining rule is in its case block). mc_run_config would pick
% rk4 at sigma = 0, so the integrator is named explicitly: sra1, the same
% solver as every other stage.
cfg.mc_preset            = 'regression_pairs_noStim_noise0_dualStd_3cond';
cfg.mc_ode_solver        = 'sra1';
% The local-Lyapunov stage (top-K local rates and local KS entropy under a
% stimulus staircase, TR 2026-09-14) runs on the steps5s twin: a new random
% step every 5 s, never returning to zero. Same bundle, chained preset.
cfg.local_lyapunov_preset = 'regression_pairs_steps5s_noise0p025_dualStd_3cond';
cfg.bursting_preset      = 'bursting_pairs';
cfg.sompolinsky_preset   = 'sompolinsky_pairs';
cfg.stf_preset           = 'single_neuron_stf';
% One unconnected reference E neuron; exactly one SFA or one STD timescale.
% Noise and heterogeneity are disabled to isolate the mechanisms.
cfg.single_neuron_preset = 'single_neuron_regression';
cfg.illustrations = true;
cfg.illustration_step_amp = 0.25;   % TR, 2026-09-16: half the mu7revised amplitude
cfg.local_lyapunov_K = 50; % explicit also in fast mode
cfg.local_lyapunov_T = 40;
cfg.local_lyapunov_step_duration_s = 5; %8 steps over40s
cfg.local_lyapunov_accumulation_start_s = 1;
cfg.local_lyapunov_warmup_s = 15;
cfg.local_lyapunov_simulation_start_s = -15; % full alignment [-14,1]
cfg.illustration_lya_start_s = 0; % first finite estimate after the [0,0.05] s segment
cfg.illustration_lya_warmup_s = 10; % alignment over [-10,0] s; simulation starts -15
cfg.illustration_display_window = [0 30];
% The step occupies the MIDDLE THIRD of the display window, [10 20] s, and
% 25 neurons of each type are saved (50 drawn); the mu7revised bundle's step
% started at 15 s and stayed on, with 4 per type saved (TR, 2026-09-15).
cfg.illustration_step_window = [10 20];
cfg.illustration_neurons_per_type = 25;
cfg.transient_gain_horizon_s = 5;
cfg.transient_gain_duration_s = 40; % fast: sample starts [25,35] s, 5-s horizon
if strcmp(run_mode,'production'), cfg.transient_gain_duration_s=60; end

% The gains shared by the two halves of figure 1 panel A; they must match.
cfg.panelA_gammas = [0.9, 1.6, 2.5];

%% The figure registry
% Order is the order make_all_paper_figures runs them. in_paper marks the ones in
% the manuscript; the rest are kept working and regenerated. Each entry is
% {name, handle, in_paper, extra-args}; the extra args are appended to the
% standard fields (run_dir, out_dir, save, visible).
F = {};
F = add(F, 'fig_introductory_concepts',       @fig_introductory_concepts,       true, ...
        {'preset_name', cfg.sompolinsky_preset, 'gammas', cfg.panelA_gammas});
% fig_energy_landscape and fig_EI_param_space are not registered (TR,
% 2026-09-14): the manuscript uses fig_EI_weights_param_space.
% Revised main 2 reads the saved representative_dynamics stage below.
F = add(F, 'fig_FI_curve',                    @fig_FI_curve,                    true, {});
F = add(F, 'fig_single_neuron_revised', @fig_single_neuron_revised, false, {});
F = add(F, 'fig_adaptation_methods_stf',      @fig_adaptation_methods,          false, ...
        {'variant', 'sfa_std_stf', 'preset_name', cfg.stf_preset});
F = add(F, 'fig_SFA_steady_state',            @fig_SFA_steady_state,            false, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_STD_steady_state',            @fig_STD_steady_state,            false, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_stim_engages_adaptation',     @fig_stim_engages_adaptation,     true, ...
        {'preset_name', cfg.bursting_preset});
F = add(F, 'fig_sensitivity_analysis_allStd', @fig_sensitivity_analysis_allStd, true, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_sensitivity_medians',         @fig_sensitivity_medians,         false, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_param_space_allStd',          @fig_param_space_allStd,          true, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_EI_weights_param_space',      @fig_EI_weights_param_space,      true, ...
        {'preset_name', cfg.preset_name});
% Rate vs stability (Codex sec. 4) and local vs finite-time LLE across trials
% (Codex sec. 3), 2026-09-14: both read the sweeps of the run directory.
F = add(F, 'fig_lle_vs_rate',                @fig_lle_vs_rate,                true, {});
F = add(F, 'fig_local_vs_finite_lle',        @fig_local_vs_finite_lle,        true, {});
F = add(F, 'fig_sfa_EOC_allStd',              @fig_sfa_EOC_allStd,              true, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_memory_capacity',             @fig_memory_capacity,             true, {});
F = add(F, 'fig_memory_capacity_example',     @fig_memory_capacity_example,     true, {});
% Double-log colour scale, log10(1 + log10(1 + density)), TR 2026-09-10: the
% single log leaves the dense core saturated and the sparse outer cloud
% invisible on this network. The figure's default stays 'log'.
F = add(F, 'fig_eig_heatmap',                 @fig_eig_heatmap,                 false, ...
        {'density_scale', 'loglog'});
% E:I imbalance examples (2026-09-14): rows = mu_EE x 0.5 / 1 / 1.5 on one seed,
% columns = regimes, each panel annotated with its matched lambda_1, mean rate
% and realised weight balance. From the same eig_heatmap stage.
F = add(F, 'fig_eig_heatmap_imbalance',       @fig_eig_heatmap_imbalance,       true, ...
        {'density_scale', 'loglog'});
F = add(F, 'fig_dc_lle',                      @fig_dc_lle,                      false, {});
% Top-K Lyapunov measures, TR 2026-09-12: the spectrum stage's figure and the
% numerical-vs-spectral abscissa (transient amplification) from the eig stage.
F = add(F, 'fig_lyapunov_spectrum',          @fig_lyapunov_spectrum,          false, {});
F = add(F, 'fig_transient_amplification',    @fig_transient_amplification,    false, {});
% Transient gain, adaptation frozen vs active, and the onset-vs-quiet contrast
% (TR 2026-09-13): from the transient_gain stage.
F = add(F, 'fig_transient_gain',             @fig_transient_gain,             false, {});
F = add(F, 'fig_transient_gain_excursions',  @fig_transient_gain_excursions,  false, {});
% Local Lyapunov exponents and local KS entropy under the stimulus staircase
% (TR 2026-09-14): from the local_lyapunov stage on cfg.local_lyapunov_preset.
F = add(F, 'fig_local_lyapunov',             @fig_local_lyapunov,             true, {});
% Numerical-method verification, TR 2026-09-10: SRA1 reshot against a 1e-10
% ode45 reference and against itself on a shared Brownian path, and Benettin
% vs QR on a reduced network. Two entries so each variant has its own folder.
F = add(F, 'fig_numerics_solver',             @fig_numerics_verification,       false, ...
        {'variant', 'solver'});
F = add(F, 'fig_numerics_lya_method',         @fig_numerics_verification,       false, ...
        {'variant', 'lya_method'});
% Check J (2026-09-14): analytic Jacobian vs central finite differences, with
% the pre-registered acceptance threshold.
F = add(F, 'fig_numerics_jacobian',          @fig_numerics_verification,       false, ...
        {'variant', 'jacobian'});
% The generated equation and conditions tables, an ordinary entry so its
% failures count in the headline number.
F = add(F, 'doc_tables',                      @fig_doc_tables,                  true, ...
        {'preset_name', cfg.preset_name});
F = [F, grouped_figure_registry(cfg.preset_name,cfg.fig_root)];
cfg.figures = F;
end

%% ------------------------------------------------------------------------
function F = add(F, name, fn, in_paper, args)
F{end+1} = struct('name', name, 'fn', fn, 'in_paper', in_paper, 'args', {args});
end
