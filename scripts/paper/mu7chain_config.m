function cfg = mu7chain_config(run_mode, figure_set)
% MU7CHAIN_CONFIG The manuscript network with chained vesicle-pool depression: four conditions, two figure sets (2026-11-12, FD-007).
%
%   cfg = mu7chain_config()                   % run_mode 'fast', figure set '4cond'
%   cfg = mu7chain_config('medium')           % Tom starts medium and larger himself
%   cfg = mu7chain_config('fast', '3cond')    % the paper's figure pass
%   cfg = mu7chain_config('fast', '4cond')    % the review pass, product MTS added
%
% WHAT CHANGED FROM mu7revisedMedium_config. The network is the manuscript's
% baseline, value for value (see the preset's case block), and every stage
% setting below is mu7revisedMedium_config's. Four things differ:
%   * the presets are the mu7chain family, whose FOUR conditions are
%       no_adaptation / sfa1_std1 (1TS) / sfa3_stdc2 (chained MTS, the
%       paper's MTS) / sfa3_std2 (product MTS, kept for the side-by-side);
%     every analysis runs all four;
%   * the run mode is an ARGUMENT (default 'fast'). Agents run 'fast' only;
%     'medium' and larger are started by Tom with this same file;
%   * cfg.param_space_seed is fixed, so reruns of the bundle sample the same
%     joint parameter-space points (FD-016);
%   * there are TWO FIGURE SETS from one run directory, chosen by the second
%     argument:
%       '3cond'  cfg.conditions = {no_adaptation, sfa1_std1, sfa3_stdc2}
%                fig_root figs/mu7chain_<run_mode>_3cond      (the paper's)
%       '4cond'  cfg.conditions = all four, product MTS last
%                fig_root figs/mu7chain_<run_mode>_4cond      (Tom's review)
%     The analyses ignore the set: run_dir is data/mu7chain_<run_mode> for
%     both, and run_all_paper_analyses is called once. make_all_paper_figures
%     is called once per set (mu7chain_run).
%
% cfg.conditions is the option's fixed name (Plan A 4.4): a cell of condition
% names to draw, in order; make_all_paper_figures passes it to every entry as
% 'conditions'. The single-neuron figure is the exception: its run holds the
% single-neuron preset's conditions, so its entry names its own list per set
% (cfg.single_neuron_conditions, below), which overrides the bundle's.
%
% WHY THE REGISTRY IS BUILT PER SET. grouped_figure_registry bakes the figure
% root into its entries (the grouped figures read their source figures from
% it), so a registry built for one fig_root and reused for the other would
% have the grouped figures read the wrong pass's outputs, or none
% (inventory R6). Each call here builds the registry for its own fig_root.
%
% THE DEPRESSION SPREAD IS NOT IN YET. The preset's name says stdSpread0p25;
% the property arrives with FD-015 and the lead adds it to the preset then.
% A run made before that has no depression spread.
%
% See also: mu7chain_run, mu7revisedMedium_config, srnn_param_preset,
%           run_all_paper_analyses, make_all_paper_figures

arguments
    run_mode (1,:) char = 'fast'
    figure_set (1,:) char {mustBeMember(figure_set, {'3cond', '4cond'})} = '4cond'
end
if ~ismember(run_mode, run_mode_names())
    error('mu7chain_config:BadRunMode', 'run_mode ''%s'' is not one of: %s.', ...
        run_mode, strjoin(run_mode_names(), ', '));
end

cfg = struct();

%% The experiment
% The manuscript's sfaEI network at mu_tilde_relative 7 on every route, SFA on
% both cell types with per-neuron ladders (tau_a_spread 0.25), per-neuron
% setpoints S_c_i = 0.2 + 0.1 randn, no external input (the Wiener process
% alone drives it), and four depression regimes on all four routes. The
% preset is written out in full; its three companions chain to it.
cfg.preset_name = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain';
cfg.run_mode    = run_mode;

% A fixed sample seed for the joint parameter-space grid (FD-016): two runs
% of this bundle at the same mode simulate the same 128 points. Any fixed
% nonnegative integer would do; this one is used nowhere else.
cfg.param_space_seed = 7202;

%% Which conditions are drawn, and where
% Order is the order the figures lay them out. In the 3-condition set the
% chain takes the product's place, colour and title; the 4-condition set adds
% the product as a fourth column, in green (manuscript_style).
switch figure_set
    case '3cond'
        cfg.conditions = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
    case '4cond'
        cfg.conditions = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2', 'sfa3_std2'};
end
cfg.figure_set = figure_set;

%% Where things land -- fixed, so this cannot touch any other run
cfg.run_dir  = ['data/mu7chain_' run_mode];                       % analyses write, figures read
cfg.fig_root = ['figs/mu7chain_' run_mode '_' figure_set];        % overwritten in place

% Figures are built and saved but do not pop up: a new figure window raises
% itself and takes keyboard focus, and a run draws dozens.
cfg.visible_figures = false;
cfg.verbose = 'minimal';   % 'verbose' | 'minimal' | 'near-none' -- how much the run prints (verbose_level)

%% Presets per stage and figure
% Memory capacity on THIS network with the Wiener process OFF: the preset
% below chains to cfg.preset_name and sets sigma_u_noise = 0. mc_run_config
% would pick rk4 at sigma = 0, so the integrator is named explicitly: sra1,
% the same solver as every other stage.
cfg.mc_preset            = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0_chainStd_4cond_mu7chain';
cfg.mc_ode_solver        = 'sra1';
% The local-Lyapunov stage runs on the steps5s twin: a new random step every
% 5 s, never returning to zero. Same bundle, chained preset.
cfg.local_lyapunov_preset = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_steps5s_noise0p025_chainStd_4cond_mu7chain';
cfg.bursting_preset      = 'bursting_pairs';
cfg.sompolinsky_preset   = 'sompolinsky_pairs';
cfg.stf_preset           = 'single_neuron_stf';
% One unconnected reference E neuron; one mechanism per condition: one SFA
% timescale, one depression pool, or the chained pool alone (four conditions).
% Noise and heterogeneity are disabled to isolate the mechanisms.
cfg.single_neuron_preset = 'single_neuron_mu7chain';
cfg.illustrations = true;
cfg.illustration_step_amp = 0.25;   % TR, 2026-09-16: half the mu7revised amplitude
cfg.local_lyapunov_K = 50; % explicit also in fast mode
cfg.local_lyapunov_T = 40;
cfg.local_lyapunov_step_duration_s = 5; % 8 steps over 40 s
cfg.local_lyapunov_accumulation_start_s = 1;
cfg.local_lyapunov_warmup_s = 15;
cfg.local_lyapunov_simulation_start_s = -15; % full alignment [-14,1]
cfg.illustration_lya_start_s = 0; % first finite estimate after the [0,0.05] s segment
cfg.illustration_lya_warmup_s = 10; % alignment over [-10,0] s; simulation starts -15
cfg.illustration_display_window = [0 30];
% The step occupies the MIDDLE THIRD of the display window, [10 20] s, and
% 25 neurons of each type are saved (50 drawn).
cfg.illustration_step_window = [10 20];
cfg.illustration_neurons_per_type = 25;
cfg.transient_gain_horizon_s = 5;
cfg.transient_gain_duration_s = 40; % sample starts [25,35] s, 5-s horizon
if strcmp(run_mode,'production'), cfg.transient_gain_duration_s=60; end

% The gains shared by the two halves of figure 1 panel A; they must match.
cfg.panelA_gammas = [0.9, 1.6, 2.5];

%% The figure registry
% mu7revisedMedium_config's registry, entry for entry. Order is the order
% make_all_paper_figures runs them. in_paper marks the ones in the manuscript;
% the rest are kept working and regenerated. Each entry is
% {name, handle, in_paper, extra-args}; the extra args are appended to the
% standard fields (run_dir, out_dir, save, visible).
F = {};
F = add(F, 'fig_introductory_concepts',       @fig_introductory_concepts,       true, ...
        {'preset_name', cfg.sompolinsky_preset, 'gammas', cfg.panelA_gammas});
F = add(F, 'fig_FI_curve',                    @fig_FI_curve,                    true, {});
% The single-neuron figure has its OWN condition list in each set (the
% lead's decision, 2026-11-12): its run holds the single-neuron preset's
% conditions, not the network's. The paper's version keeps its three columns
% (no adaptation, one SFA timescale, one depression pool); the review version
% adds chained depression alone. An entry's own args follow the bundle-wide
% 'conditions' in the call, so these override cfg.conditions.
switch figure_set
    case '3cond'
        sn_conditions = {'no_adaptation', 'sfa1_std0', 'sfa0_std1'};
    case '4cond'
        sn_conditions = {'no_adaptation', 'sfa1_std0', 'sfa0_std1', 'sfa0_stdc2'};
end
cfg.single_neuron_conditions = sn_conditions;
F = add(F, 'fig_single_neuron_revised', @fig_single_neuron_revised, false, ...
        {'conditions', sn_conditions});
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
F = add(F, 'fig_lle_vs_rate',                @fig_lle_vs_rate,                true, {});
F = add(F, 'fig_local_vs_finite_lle',        @fig_local_vs_finite_lle,        true, {});
F = add(F, 'fig_sfa_EOC_allStd',              @fig_sfa_EOC_allStd,              true, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_memory_capacity',             @fig_memory_capacity,             true, {});
F = add(F, 'fig_memory_capacity_example',     @fig_memory_capacity_example,     true, {});
F = add(F, 'fig_eig_heatmap',                 @fig_eig_heatmap,                 false, ...
        {'density_scale', 'loglog'});
F = add(F, 'fig_eig_heatmap_imbalance',       @fig_eig_heatmap_imbalance,       true, ...
        {'density_scale', 'loglog'});
F = add(F, 'fig_dc_lle',                      @fig_dc_lle,                      false, {});
F = add(F, 'fig_lyapunov_spectrum',          @fig_lyapunov_spectrum,          false, {});
F = add(F, 'fig_transient_amplification',    @fig_transient_amplification,    false, {});
F = add(F, 'fig_transient_gain',             @fig_transient_gain,             false, {});
F = add(F, 'fig_transient_gain_excursions',  @fig_transient_gain_excursions,  false, {});
F = add(F, 'fig_local_lyapunov',             @fig_local_lyapunov,             true, {});
F = add(F, 'fig_numerics_solver',             @fig_numerics_verification,       false, ...
        {'variant', 'solver'});
F = add(F, 'fig_numerics_lya_method',         @fig_numerics_verification,       false, ...
        {'variant', 'lya_method'});
F = add(F, 'fig_numerics_jacobian',          @fig_numerics_verification,       false, ...
        {'variant', 'jacobian'});
F = add(F, 'doc_tables',                      @fig_doc_tables,                  true, ...
        {'preset_name', cfg.preset_name});
% Built for THIS set's fig_root (see "why the registry is built per set").
F = [F, grouped_figure_registry(cfg.preset_name, cfg.fig_root)];
cfg.figures = F;
end

%% ------------------------------------------------------------------------
function F = add(F, name, fn, in_paper, args)
F{end+1} = struct('name', name, 'fn', fn, 'in_paper', in_paper, 'args', {args});
end
