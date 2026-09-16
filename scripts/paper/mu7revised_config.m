function cfg = mu7revised_config(run_mode)
% MU7REVISED_CONFIG Independent mu=7 bundle for the next full analysis run.
% Defaults to FAST. Use mu7revised_config('medium') for manuscript estimates.
% Every network/stage setting is explicit here; no parent config is called.
% Relative data/figure paths keep the bundle portable between computers.
% New analyses: midpoint input-step dynamics, separate 1TS single-neuron
% mechanisms, and a genuine 5-s transient horizon (40-s gain trajectories).
% Current archived mu7 data are never overwritten or reused for this bundle.
arguments
    run_mode (1,:) char {mustBeMember(run_mode,{'fast','medium','medium2','production'})} = 'fast'
end

cfg = struct();

%% The experiment
% The 3-condition sfaEI network at a BASELINE mu_tilde_relative OF 7 on every
% route (TR, 2026-09-15; the mu5 bundle of 2026-09-14 ran at 5, the mu8p25
% family was 1.5 x 5.5): SFA on BOTH cell types (c = [0.5 0.5], the same
% ladder on I as on E), a per-neuron setpoint S_c_i = 0.2 + 0.1*randn,
% per-neuron SFA ladders (tau_a_spread 0.25) and NO external input -- the
% Wiener process alone drives it. Regimes: no_adaptation / sfa1_std1 /
% sfa3_std2 (STD not strength-matched; the matched bundles are stdScaled_*,
% stdUsage_*, stdSingleMatched_*). The four mu sweeps and the All Weights
% (level_of_chaos) sweep span -75% .. +75% of the preset, in the 1-D sweeps
% and the joint grid alike, and the joint sample is 128 points in every mode.
% See the preset, which is written out in full; its two twins below chain to it.
%
% 'fast': 5 levels x 5 reps per 1-D sweep (TR, 2026-09-15), 7 tau levels, 128
% joint samples, 5 MC trials (the sign-flip floor is p = 0.0625), 10-s
% sensitivity runs and 20-s joint runs. About three and a half hours for the
% whole pipeline on the n = 500 network; switch run_mode to 'medium' for
% manuscript numbers (15 MC trials, exact sign-flip over 2^15).
cfg.preset_name = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0p025_dualStd_3cond_mu7revised';
cfg.run_mode    = run_mode;

%% Where things land -- both fixed, so this cannot touch any other run
cfg.run_dir  = ['data/mu7revised_' run_mode];   % analyses write, figures read
cfg.fig_root = ['figs/mu7revised_' run_mode];   % overwritten in place

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
cfg.mc_preset            = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0_dualStd_3cond_mu7revised';
cfg.mc_ode_solver        = 'sra1';
% The local-Lyapunov stage (top-K local rates and local KS entropy under a
% stimulus staircase, TR 2026-09-14) runs on the steps5s twin: a new random
% step every 5 s, never returning to zero. Same bundle, chained preset.
cfg.local_lyapunov_preset = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_steps5s_noise0p025_dualStd_3cond_mu7revised';
cfg.bursting_preset      = 'bursting_pairs';
cfg.sompolinsky_preset   = 'sompolinsky_pairs';
cfg.stf_preset           = 'single_neuron_stf';
% One unconnected reference E neuron; exactly one SFA or one STD timescale.
% Noise and heterogeneity are disabled to isolate the mechanisms.
cfg.single_neuron_preset = 'single_neuron_mu7revised';
cfg.illustrations = true;
cfg.illustration_step_amp = 0.5;
cfg.local_lyapunov_K = 100; % explicit also in fast mode
cfg.local_lyapunov_accumulation_start_s = 1;
cfg.local_lyapunov_warmup_s = 15;
cfg.local_lyapunov_simulation_start_s = -15; % full alignment [-14,1]
cfg.illustration_lya_start_s = 0; % first finite estimate after the [0,0.05] s segment
cfg.illustration_lya_warmup_s = 10; % alignment over [-10,0] s; simulation starts -15
cfg.illustration_display_window = [0 30]; % step exactly midway, at 15 s
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
