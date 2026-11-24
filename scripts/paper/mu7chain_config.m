function cfg = mu7chain_config(run_mode, figure_set)
% MU7CHAIN_CONFIG The manuscript network with chained vesicle-pool depression: four conditions at 'essential', three from 'fast' upward (2026-11-12, FD-007; FD-030).
%
%   cfg = mu7chain_config()                   % run_mode 'fast', its figure set '3cond'
%   cfg = mu7chain_config('medium')           % Tom starts medium and larger himself
%   cfg = mu7chain_config('fast', '3cond')    % the paper's figure pass
%   cfg = mu7chain_config('fast', '4cond')    % product MTS added: only for a run
%                                             % made before FD-030, which simulated it
%   cfg = mu7chain_config('essential')        % a first look: figures 2-4 and brief
%                                             % sweeps; four conditions, set '4cond'
%
% CONDITIONS BY MODE (FD-030, Tom, 2026-11-16). Product MTS is not in the
% paper, so from 'fast' upward the analyses simulate three conditions
% (cfg.run_conditions) and the three-condition set is the one drawn
% (cfg.figure_sets); 'essential' keeps all four and the four-condition set.
% Stated once, under "Which conditions are SIMULATED" below. A figure_set
% left empty is the mode's own.
%
% 'essential' (Tom, 2026-11-12, FD-023): A LOOK AT A NETWORK BEFORE COMMITTING
% TO LONG RUNS, sized for about 35 minutes on 12 pool workers, simulation and
% figures together (an estimate from the preliminary fast run's stage
% timings; not measured when this was written). It is the one mode of this
% bundle that does not run every stage or draw every figure:
%   runs     sweeps (the seven 1-D sensitivity sweeps and the joint
%            parameter-space sample), eig_heatmap (the reference example
%            only), local_lyapunov, paper_illustrations, and the opt-in
%            lyapunov_spectrum_n80 that every mode runs (FD-032)
%   skips    the tau sweep, memory_capacity, mc_example, numerics,
%            lyapunov_spectrum, transient_gain, dc_lle (the data of main
%            figures 5 and 6 and of the supplementary checks)
%   sizes    1-D sweeps 4 levels x 2 reps (fast: 5 x 5); joint sample 32
%            points (fast: 128); every window, integrator and stage size
%            otherwise is fast's (analysis_run_config, run_mode_names)
%   figures  ONE set, four conditions, in figs/mu7chain_essential_4cond: the
%            grouped main figures 2, 3 and 4, the sensitivity sheets and
%            medians, the joint-sample histograms, the joint sample coloured
%            by realised E:I weight balance, and the stand-alone versions of
%            what the grouped figures compose (eigenvalue heat map, lambda_1
%            against rate). mu7chain_run makes that one pass; the 3cond set
%            can still be asked for by name and is not made by default.
%   reading  TWO REPS PER LEVEL: the interquartile bands of figure 4 are the
%            spread of two runs and mean almost nothing. The run shows the
%            regime (where lambda_1 and the rates sit, which way each
%            parameter pushes), not trends, and nothing from it goes in the
%            paper.
% The selection is three config fields read by run_all_paper_analyses
% (paper_stage_plan): cfg.stages, cfg.sweeps, cfg.eig_heatmap_examples. In
% every other mode the last two are empty, which means everything, and
% cfg.stages names every stage plus the opt-in lyapunov_spectrum_n80
% (FD-032), so nothing is skipped.
%
% THE SPECTRAL SUPPLEMENT (FD-040 to FD-046; Tom, 2026-11-20: "integrate
% these figures into the full run"), every mode: the opt-in stages
% outlier_spectra, frequency_response, eigenvalue_density and
% outlier_structures and their figure entries (all 'appendix'), added at
% the end by with_spectral_supplement. Sizes by run mode are the stages'
% own (outlier_spectra_stage_plan: 2 seeds and 20 s windows at 'essential'
% and 'fast', 10 seeds and 200 s for reference and E-to-E x 1.5 above;
% outlier_structures_stage_plan: fifteen networks, 1 seed at 'essential'
% and 'fast', 5 above).
%
% FIGURE 3 IS THE SWITCHING VERSION (Tom, 2026-11-20: "we will use the
% switching version"; docs/Plan_Figure3_update.md, step 6), every mode:
% with_switch_spectra(cfg, cfg.local_lyapunov_K, 40, 1, 'as_default', true)
% at the end. So the local-Lyapunov stage tracks K = 150 directions
% (cfg.local_lyapunov_K, in every mode since 2026-11-21; 100 from
% 2026-11-20; see the note where it is set) and saves the Jacobian's
% spectra and H at 40 states of its own switching run, at 1, 2, ..., 40 s;
% the registered fig_main3_grouped is the four-row composite (input,
% expansion, heat map with the Ritz overlay, the real-axis histogram),
% saved as fig_main3_grouped/Fig_Main3_Grouped.png (the name FD-031 gave
% the paper's Figure 3), and the baseline composite (FD-031's three-row
% paper layout on the separate eig_heatmap run) stays as
% fig_main3_grouped_baseline (Fig_Main3_Grouped_baseline.png), in the
% appendix. K also changes fig_local_lyapunov, which reads the same stage.
%
% THE RANKED LYAPUNOV SPECTRUM (FD-032), every mode: the stage
% lyapunov_spectrum_n80 (top 240 exponents per condition on the 80-neuron
% preset cfg.lyapunov_spectrum_n80_preset, into
% <run_dir>/lyapunov_spectrum_n80/) and the supplementary figure entry
% fig_lyapunov_ranked_spectrum, in both figure sets and at 'essential'.
%
% WHAT CHANGED FROM mu7revisedMedium_config. The network is the manuscript's
% baseline, value for value (see the preset's case block), and every stage
% setting below is mu7revisedMedium_config's. Four things differ:
%   * the presets are the mu7chain family, whose FOUR conditions are
%       no_adaptation / sfa1_std1 (1TS) / sfa3_stdc2 (chained MTS, the
%       paper's MTS) / sfa3_std2 (product MTS, kept for the side-by-side);
%     every analysis ran all four until FD-030; since then all four at
%     'essential' and the first three from 'fast' upward;
%   * the run mode is an ARGUMENT (default 'fast'). Agents run 'fast' only;
%     'medium' and larger are started by Tom with this same file;
%   * cfg.param_space_seed is fixed, so reruns of the bundle sample the same
%     joint parameter-space points (FD-016);
%   * there are TWO FIGURE SETS, chosen by the second argument (since FD-030
%     a mode draws one of them, cfg.figure_sets; a three-condition run can
%     draw only '3cond'):
%       '3cond'  cfg.conditions = {no_adaptation, sfa1_std1, sfa3_stdc2}
%                fig_root figs/mu7chain_<run_mode>_3cond      (the paper's)
%       '4cond'  cfg.conditions = all four, product MTS last
%                fig_root figs/mu7chain_<run_mode>_4cond      (Tom's review)
%     The analyses ignore the set (they simulate cfg.run_conditions): run_dir
%     is data/mu7chain_<run_mode> for both, and run_all_paper_analyses is
%     called once. make_all_paper_figures is called once per set of the mode
%     (mu7chain_run).
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
% THE DEPRESSION SPREAD (stdSpread0p25 in the preset's name) is in the
% preset since 2026-11-12: std_spread = [0.25 0.25] (FD-015). A run made
% before that (data/mu7chain_fast_prelim and its copies) has none.
%
% See also: mu7chain_run, mu7revisedMedium_config, srnn_param_preset,
%           run_all_paper_analyses, make_all_paper_figures

arguments
    run_mode (1,:) char = 'fast'
    figure_set (1,:) char {mustBeMemberOrEmpty(figure_set)} = ''
end
if ~ismember(run_mode, run_mode_names())
    error('mu7chain_config:BadRunMode', 'run_mode ''%s'' is not one of: %s.', ...
        run_mode, strjoin(run_mode_names(), ', '));
end

%% Which conditions are SIMULATED, and which figure sets are drawn, by mode
% THE ONE PLACE this bundle family states them (FD-030; Tom, 2026-11-16:
% product MTS "will not go into the paper. We will only keep 4 conditions
% for the essential run and for test."). mu7chain_gainbias_config and
% mu7chain_AC_config take both from here.
%   essential               all four (run_conditions empty = every condition
%                           of the preset, so the run is exactly what it was
%                           before FD-030): no_adaptation, sfa1_std1,
%                           sfa3_stdc2, sfa3_std2; the four-condition set
%   fast and every larger   three: no_adaptation, sfa1_std1, sfa3_stdc2;
%                           product MTS (sfa3_std2) is not simulated; the
%                           three-condition figure set only
% cfg.run_conditions is read by run_all_paper_analyses, which passes it to
% every stage (bundle_run_conditions, preset_with_conditions) and records it
% in run_manifest (conditions_run); the preset and its four conditions are
% unchanged. It is listed in the preset's order. A three-condition run saves
% about a quarter of every stage's compute. The 4cond set can still be asked
% for by name at any mode (to draw a run made before FD-030); on a
% three-condition run make_all_paper_figures refuses it before drawing.
% (medium2 and production follow medium: the lead's decision under FD-030.)
THREE = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
if strcmp(run_mode, 'essential')
    run_conditions = {};
    figure_sets = {'4cond'};
else
    run_conditions = THREE;
    figure_sets = {'3cond'};
end
% With no figure set named, the mode's first (the one it draws).
if isempty(figure_set); figure_set = figure_sets{1}; end

cfg = struct();

%% The experiment
% The manuscript's sfaEI network at mu_tilde_relative 7 on every route, SFA on
% both cell types with per-neuron ladders (tau_a_spread 0.25), per-neuron
% setpoints S_c_i = 0.2 + 0.1 randn, no external input (the Wiener process
% alone drives it), and four depression regimes on all four routes. The
% preset is written out in full; its three companions chain to it.
cfg.preset_name = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain';
cfg.run_mode    = run_mode;

% THE INTEGRATION RATE, every run mode (Tom, 2026-11-17, FD-031 item 14): 500
% Hz. A 400/500/800-Hz test showed second-order step error; 500 Hz halves it
% at about 1.3x the cost. Stated once here; mu7chain_gainbias_config and
% mu7chain_AC_config inherit it. Runs made before this (400 Hz at production
% and medium, 200 Hz in several stages at fast and fastmedium) are left as
% they are.
% run_all_paper_analyses hands it as 'fs' to every stage: the sweeps
% (replacing the analysis_run_config mode table's rate), eig_heatmap,
% local_lyapunov, lyapunov_spectrum, transient_gain, dc_lle, the paper
% illustrations, memory capacity and mc_example, block C's grids
% (protocol.fs) and the noise-switch rollouts. The opt-in stages built since
% (lyapunov_spectrum_n80, outlier_spectra, outlier_structures) read it and
% cfg.ode_solver in their plans; frequency_response and eigenvalue_density
% simulate nothing. One stage keeps its own rates by design: numerics (an fs
% ladder and a solver comparison).
cfg.fs = 500;

% THE INTEGRATOR, every stage and every run mode (Tom, 2026-11-17): SRA1,
% Roessler's strong-order-1.5 scheme for additive noise. Without noise it is
% Ralston's second-order Runge-Kutta, so the noise-free stages (memory
% capacity, block C's grids, the noise-off Lyapunov spectrum, the
% single-neuron illustration) use the same scheme and the Methods state one.
% Handed beside cfg.fs to the same stages (stage_solver); numerics excepted.
cfg.ode_solver = 'sra1';

% A fixed sample seed for the joint parameter-space grid (FD-016): two runs
% of this bundle at the same mode simulate the same 128 points. Any fixed
% nonnegative integer would do; this one is used nowhere else.
cfg.param_space_seed = 7202;

%% Which stages run (FD-023)
% Empty means every stage, every sweep and the eig-heatmap stage's default
% examples: every mode but 'essential'. See the header for what 'essential'
% is for. Its stages are the ones main figures 2, 3 and 4 read:
%   figure 2   representative_dynamics, saved by paper_illustrations
%   figure 3   local_lyapunov, and the reference example of eig_heatmap
%   figure 4   the 1-D sensitivity sweeps and the joint sample
%
% PLUS, IN EVERY MODE, the opt-in stage 'lyapunov_spectrum_n80' (FD-032; Tom,
% 2026-11-19: "it should get made when we do a run"): the top 240 Lyapunov
% exponents per condition on the 80-neuron version of this network
% (cfg.lyapunov_spectrum_n80_preset), for the supplementary ranked spectrum
% and histogram. An opt-in stage runs only when cfg.stages names it, so
% from 'fast' upward cfg.stages is every stage of the first list written
% out (paper_stage_names), then this one; nothing is skipped, as before.
% Its settings (seeds by mode; K, window and rate fixed) are in
% lyapunov_spectrum_n80_stage_plan.
% PLUS, IN EVERY MODE (Plan 2, integration, 2026-11-21), the opt-in stages
% 'concept_examples' (the four example networks of Figure 1 and the gallery
% S3; concept_examples_stage_plan) and 'uncoupled_reference' (the growth
% rate of the nearly uncoupled network, the dotted line of Figure 4C;
% uncoupled_reference_stage_plan). Both take minutes at any size.
is_essential = strcmp(run_mode, 'essential');
NEW_STAGES = {'lyapunov_spectrum_n80', 'concept_examples', 'uncoupled_reference'};
cfg.stages = [paper_stage_names(), NEW_STAGES];
cfg.sweeps = {};
cfg.eig_heatmap_examples = [];
if is_essential
    cfg.stages = [{'sweeps', 'eig_heatmap', 'local_lyapunov', 'paper_illustrations'}, NEW_STAGES];
    cfg.sweeps = {'sensitivity', 'param_space'};
    % The reference network alone: fig_eig_heatmap (and figure 3 through it)
    % reads the reference; the two E:I-imbalanced examples feed only
    % fig_eig_heatmap_imbalance, which this mode does not draw.
    cfg.eig_heatmap_examples = struct('label', 'reference', 'overrides', struct());
end

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
% The sets mu7chain_run draws for this mode, one make_all_paper_figures pass
% each, and the conditions the analyses simulate: stated by mode at the top
% of this function. 'essential' makes the four-condition set alone (Tom,
% 2026-11-12: "the 3 condition figures don't need to be generated"); every
% larger mode the three-condition set alone (Tom, 2026-11-16, FD-030; until
% then both sets).
cfg.figure_sets = figure_sets;
cfg.run_conditions = run_conditions;

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
% would pick rk4 at sigma = 0, so the integrator is named explicitly: the
% bundle's cfg.ode_solver (sra1), the same solver as every other stage.
% The lean path (no state trajectory held): at 500 Hz the full path's
% 760-s state series is 2.5 times its 200-Hz size (about 16 GB a trial
% then), and the lean path gives the same numbers (test_mc_lean_path).
cfg.mc_preset            = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0_chainStd_4cond_mu7chain';
cfg.mc_ode_solver        = cfg.ode_solver;
cfg.mc_lean              = true;
% The local-Lyapunov stage runs on the steps5s twin: a new random step every
% 5 s, never returning to zero. Same bundle, chained preset.
cfg.local_lyapunov_preset = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_steps5s_noise0p025_chainStd_4cond_mu7chain';
% The ranked Lyapunov spectrum (stage lyapunov_spectrum_n80, FD-032) runs on
% the 80-neuron twin: this network with only its size changed. Same bundle.
cfg.lyapunov_spectrum_n80_preset = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain_n80';
cfg.bursting_preset      = 'bursting_pairs';
cfg.sompolinsky_preset   = 'sompolinsky_pairs';
cfg.stf_preset           = 'single_neuron_stf';
% One unconnected reference E neuron; one mechanism per condition: one SFA
% timescale, one depression pool, or the chained pool alone (four conditions).
% Noise and heterogeneity are disabled to isolate the mechanisms.
cfg.single_neuron_preset = 'single_neuron_mu7chain';
cfg.illustrations = true;
cfg.illustration_step_amp = 0.25;   % TR, 2026-09-16: half the mu7revised amplitude
% Top-K of the local-Lyapunov stage (Figure 3's switching data: the
% expansion row, the Ritz values and the exponent ticks): K = 150 IN EVERY
% RUN MODE (Tom, 2026-11-21: "go with 150"; Plan 2, package A). Explicit in
% every mode, overriding run_local_lyapunov's per-mode default, and passed
% to with_switch_spectra at the end, so the stage, the switching spectra
% and the registered Figure 3 all use it.
% History, recorded because the decisions were made on two machines over
% four days: on 2026-11-17 (FD-031 item 8, the PC) K was 50 at 'essential',
% 'fast' and 'fastmedium' and 150 at 'medium', 'medium2' and 'production'
% (about 11-12 min for three conditions on an idle machine); on 2026-11-20
% (the MacBook) K = 100 in every mode replaced both; on 2026-11-21 K = 150
% in every mode replaced that. What FD-031 found still holds: the
% multiple-timescale expansion does not converge in K, so Figure 3 calls it
% a lower bound, and fig_main3_stability writes the K the run used into the
% figure's provenance. The cost grows about linearly in K over 100 (the
% top-K pass is matrix-free, N x K per step), so the stage is about 1.5
% times its K = 100 time; not measured when this was written.
cfg.local_lyapunov_K = 150;
cfg.local_lyapunov_T = 40;
cfg.local_lyapunov_step_duration_s = 5; % 8 steps over 40 s
% Timing of the switching run (Tom, 2026-11-20; docs/Plan_Figure3_update.md):
% [-15,-10] s the network runs with no stimulus and nothing is tracked,
% [-10,0] s the tracked directions align (no stimulus), [0,40] s accumulation
% with the stimulus on. Before: accumulation from 1 s, alignment [-14,1] s.
cfg.local_lyapunov_accumulation_start_s = 0;
cfg.local_lyapunov_warmup_s = 10;
cfg.local_lyapunov_simulation_start_s = -15; % alignment [-10,0]
cfg.illustration_lya_start_s = 0; % first finite estimate after the [0,0.05] s segment
cfg.illustration_lya_warmup_s = 10; % alignment over [-10,0] s; simulation starts -15
% 20 s shown, in every run mode (Tom, 2026-11-22, second round on Figure 2:
% 30 s made the figure too tight left to right). Was [0 30] with the step
% over [10 20] s. The single-neuron illustration of the same stage (S1,
% panel A) shares this window.
cfg.illustration_display_window = [0 20];
% The step occupies the MIDDLE THIRD of the display window, [20/3 40/3] s
% (about 6.7 to 13.3 s), and 25 neurons of each type are saved.
cfg.illustration_step_window = [20/3 40/3];
cfg.illustration_neurons_per_type = 25;
cfg.transient_gain_horizon_s = 5;
cfg.transient_gain_duration_s = 40; % sample starts [25,35] s, 5-s horizon
if strcmp(run_mode,'production'), cfg.transient_gain_duration_s=60; end

% The gains shared by the two halves of figure 1 panel A; they must match.
cfg.panelA_gammas = [0.9, 1.6, 2.5];

%% The figure registry
% mu7revisedMedium_config's registry, entry for entry. Order is the order
% make_all_paper_figures runs them. Each entry is
% add(F, name, handle, group, extra-args), built by paper_figure_entry; the
% extra args are appended to the standard fields (run_dir, out_dir, save,
% visible).
%
% THE GROUPS (Plan 2, package A, 2026-11-21; paper_figure_groups and the
% manuscript repository's figures_main_supp_appen.md): every entry is
% 'main', 'supplement' or 'appendix', and in_paper is derived from it (main
% and supplement). The main figures here are the grouped composites
% (grouped_figure_registry: main Figures 1 to 4, supplementary S1 to S5
% and the leading-direction figure,
% the old composites in the appendix; with_spectral_supplement adds main 5
% and 6 and S6, S7 (and S8 in the appendix since Plan 3: its panels are
% S6 D and E), mu7chain_AC_config main 7 and S9). Every single
% figure below is 'appendix': it is in the plan's appendix, or a panel of a
% supplementary composite. Until 2026-11-21 several of them were in_paper
% true, as panels of the old supplement.
F = {};
% The two entries that simulate while drawing take the bundle's rate and
% integrator too (Tom, 2026-11-17: every paper simulation at SRA1, 500 Hz).
sim_args = {'fs', cfg.fs, 'ode_solver', cfg.ode_solver};
F = add(F, 'fig_introductory_concepts',       @fig_introductory_concepts,       'appendix', ...
        [{'preset_name', cfg.sompolinsky_preset, 'gammas', cfg.panelA_gammas}, sim_args]);
F = add(F, 'fig_FI_curve',                    @fig_FI_curve,                    'appendix', {});
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
F = add(F, 'fig_single_neuron_revised', @fig_single_neuron_revised, 'appendix', ...
        {'conditions', sn_conditions});
F = add(F, 'fig_adaptation_methods_stf',      @fig_adaptation_methods,          'appendix', ...
        {'variant', 'sfa_std_stf', 'preset_name', cfg.stf_preset});
F = add(F, 'fig_SFA_steady_state',            @fig_SFA_steady_state,            'appendix', ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_STD_steady_state',            @fig_STD_steady_state,            'appendix', ...
        {'preset_name', cfg.preset_name});
% What 1TS and the chained MTS share and where they differ (FD-014): the
% first-cycle match and the steady-state difference under a sinusoidal rate,
% rate steps, and the steady state against rate. Not in mu7revisedMedium's
% registry, which has no chain. Built from the preset's route models, no
% run directory; the 4cond set adds the product in green.
F = add(F, 'fig_std_1TS_vs_MTS',              @fig_std_1TS_vs_MTS,              'appendix', ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_stim_engages_adaptation',     @fig_stim_engages_adaptation,     'appendix', ...
        [{'preset_name', cfg.bursting_preset}, sim_args]);
F = add(F, 'fig_sensitivity_analysis_allStd', @fig_sensitivity_analysis_allStd, 'appendix', ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_sensitivity_medians',         @fig_sensitivity_medians,         'appendix', ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_param_space_allStd',          @fig_param_space_allStd,          'appendix', ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_EI_weights_param_space',      @fig_EI_weights_param_space,      'appendix', ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_lle_vs_rate',                @fig_lle_vs_rate,                'appendix', {});
F = add(F, 'fig_local_vs_finite_lle',        @fig_local_vs_finite_lle,        'appendix', {});
F = add(F, 'fig_sfa_EOC_allStd',              @fig_sfa_EOC_allStd,              'appendix', ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_memory_capacity',             @fig_memory_capacity,             'appendix', {});
F = add(F, 'fig_memory_capacity_example',     @fig_memory_capacity_example,     'appendix', {});
F = add(F, 'fig_eig_heatmap',                 @fig_eig_heatmap,                 'appendix', ...
        {'density_scale', 'loglog'});
F = add(F, 'fig_eig_heatmap_imbalance',       @fig_eig_heatmap_imbalance,       'appendix', ...
        {'density_scale', 'loglog'});
% The same eigenvalues projected on the real axis, against the projection of
% a uniform disk (Tom, 2026-11-19; additional figure, not in the paper). From
% the eig_heatmap stage; without the stage's saved gains it draws the
% unit-gain disk.
F = add(F, 'fig_eig_real_histogram',          @fig_eig_real_histogram,          'appendix', {});
F = add(F, 'fig_dc_lle',                      @fig_dc_lle,                      'appendix', {});
F = add(F, 'fig_lyapunov_spectrum',          @fig_lyapunov_spectrum,          'appendix', {});
% The ranked top-240 spectrum and its histogram on the 80-neuron network
% (FD-032, supplementary), from the stage lyapunov_spectrum_n80, which this
% bundle runs in every mode.
F = add(F, 'fig_lyapunov_ranked_spectrum',   @fig_lyapunov_ranked_spectrum,   'appendix', {});
F = add(F, 'fig_transient_amplification',    @fig_transient_amplification,    'appendix', {});
F = add(F, 'fig_transient_gain',             @fig_transient_gain,             'appendix', {});
F = add(F, 'fig_transient_gain_excursions',  @fig_transient_gain_excursions,  'appendix', {});
F = add(F, 'fig_local_lyapunov',             @fig_local_lyapunov,             'appendix', {});
F = add(F, 'fig_numerics_solver',             @fig_numerics_verification,       'appendix', ...
        {'variant', 'solver'});
F = add(F, 'fig_numerics_lya_method',         @fig_numerics_verification,       'appendix', ...
        {'variant', 'lya_method'});
F = add(F, 'fig_numerics_jacobian',          @fig_numerics_verification,       'appendix', ...
        {'variant', 'jacobian'});
F = add(F, 'doc_tables',                      @fig_doc_tables,                  'appendix', ...
        {'preset_name', cfg.preset_name});
% Built for THIS set's fig_root (see "why the registry is built per set").
F = [F, grouped_figure_registry(cfg.preset_name, cfg.fig_root)];

% 'essential' keeps the entries whose data its stages make, and no other: a
% figure whose stage did not run errors by design (resolve_data_file), so it
% must not be listed. Kept by NAME from the registry above, so each entry is
% the one every other mode draws, in the same order. What each reads:
%   fig_sensitivity_analysis_allStd   the 1-D sweeps (the sensitivity sheets)
%   fig_sensitivity_medians           the 1-D sweeps
%   fig_param_space_allStd            the joint sample (histograms)
%   fig_EI_weights_param_space        the joint sample, coloured by realised
%                                     E:I weight balance
%   fig_lle_vs_rate                   the joint sample and the 1-D sweeps
%   fig_eig_heatmap                   eig_heatmap, the reference example
%   fig_eig_real_histogram            eig_heatmap, the reference example
%   fig_main2_grouped                 representative_dynamics (paper_illustrations)
%   fig_main3_grouped                 local_lyapunov and eig_heatmap
%   fig_main4_grouped                 the 1-D sweeps and the joint sample
%   fig_lyapunov_ranked_spectrum      lyapunov_spectrum_n80 (FD-032)
%   fig_main1_grouped, fig_supp3_grouped   concept_examples (Plan 2)
%   fig_supp1_grouped                 single_neuron (paper_illustrations) and
%                                     the preset's adaptation models
%   fig_supp4_grouped                 the 1-D sweeps
%   fig_supp5_grouped                 the joint sample and the 1-D sweeps
% (main Figure 5, S6, S7 and, in the appendix, S8 are added by
% with_spectral_supplement in
% every mode; Figures 6 and 7, S2, S9 and fig_supp_leading_grouped need
% stages this mode skips: the leading-direction figure reads the tau sweep.)
if is_essential
    keep = {'fig_sensitivity_analysis_allStd', 'fig_sensitivity_medians', ...
            'fig_param_space_allStd', 'fig_EI_weights_param_space', 'fig_lle_vs_rate', ...
            'fig_eig_heatmap', 'fig_eig_real_histogram', 'fig_lyapunov_ranked_spectrum', ...
            'fig_main1_grouped', 'fig_main2_grouped', 'fig_main3_grouped', 'fig_main4_grouped', ...
            'fig_supp1_grouped', 'fig_supp3_grouped', 'fig_supp4_grouped', 'fig_supp5_grouped'};
    names = cellfun(@(f) f.name, F, 'UniformOutput', false);
    missing = setdiff(keep, names);
    assert(isempty(missing), 'mu7chain_config:EssentialRegistry', ...
        'The essential set names entries the registry lacks: %s.', strjoin(missing, ', '));
    F = F(ismember(names, keep));
end
cfg.figures = F;

%% The spectral supplement and the switching Figure 3: this bundle's defaults
% (see the header). Last, so that both see the registry above, the
% essential selection included.
cfg = with_spectral_supplement(cfg);
cfg = with_switch_spectra(cfg, cfg.local_lyapunov_K, 40, 1, 'as_default', true);
end

%% ------------------------------------------------------------------------
function F = add(F, name, fn, group, args)
F{end+1} = paper_figure_entry(name, fn, group, args);
end

function mustBeMemberOrEmpty(figure_set)
% '' (the mode's own set) or one of the two sets.
if ~isempty(figure_set)
    mustBeMember(figure_set, {'3cond', '4cond'});
end
end
