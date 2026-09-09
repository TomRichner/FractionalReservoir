function cfg = single_multi_TS_independent_config()
% SINGLE_MULTI_TS_INDEPENDENT_CONFIG Single vs multi timescale, self-contained.
%
%   run_dir = run_all_paper_analyses(single_multi_TS_independent_config());
%   results = make_all_paper_figures(single_multi_TS_independent_config());
%
% EVERY SETTING IS STATED HERE. This does not call paper_config and does not
% inherit from any other *_config; read this file and you know the whole run.
% Earlier configs were chained modifications of paper_config, and the cost of
% that surfaced on 2026-09-09: two experiments that changed the network had
% silently run memory capacity on a DIFFERENT network (cfg.mc_preset =
% 'mc_pairs_dualStd', n = 300, mu 3/-4), because that line lived in the parent.
% Going forward, every *_config.m is independent (TR).
%
% MEMORY CAPACITY RUNS ON THE SAME PRESET AS THE SWEEPS. That is the one
% substantive difference from the paper_config lineage, where MC has its own
% smaller network. Consequences worth knowing:
%   - n = 500 rather than 300, so the readout has 500 features. MC is well
%     posed only while N_train = T_train/T_hold exceeds that; at 'medium'
%     that is 300/0.3 = 1000 hold-samples, fine, but at 'fast' it is 200 and
%     the readout is UNDERDETERMINED. Do not read MC numbers from a 'fast' run
%     of this config. T_train is mc_run_config's per-mode table, not a cfg knob.
%   - MC runs the preset's own conditions, so its four-regime sheet becomes the
%     three-regime one every other figure shows.
%   - The 'synaptic' readout requires all of a presynaptic type's routes to
%     carry identical STD, which this preset's all-four-routes config does.
%
% The four figures that are DELIBERATELY different networks keep their own
% presets, as before: two single-neuron mechanism cartoons, the Sompolinsky
% reproduction and the bursting network make points the 500-neuron recurrent
% network cannot make. They are named explicitly below, not inherited.
%
% RERUNNING: delete data/single_multi_TS_independent first.
% run_all_paper_analyses refuses a run directory that is not absent or empty.
%
% See also: single_multi_TS_independent_run, run_all_paper_analyses,
%           make_all_paper_figures, srnn_param_preset

cfg = struct();

%% The experiment
% The paper's 3-condition network with mean connectivity 50% stronger on all
% four routes (mu_tilde_relative 8.25). Regimes: no_adaptation / sfa1_std1 /
% sfa3_std2. To run the baseline network instead, name
% 'celltype_pairs_Sc0p2_noise0p025_dualStd_3cond' here and change the two
% roots below so the runs do not collide.
cfg.preset_name = 'celltype_pairs_Sc0p2_noise0p025_dualStd_3cond_mu8p25';
cfg.run_mode    = 'medium';

%% Where things land -- both fixed, so this cannot touch any other run
cfg.run_dir  = 'data/single_multi_TS_independent';   % analyses write, figures read
cfg.fig_root = 'figs/single_multi_TS_independent';   % overwritten in place

% Figures are built and saved but do not pop up: a new figure window raises
% itself and takes keyboard focus, and a run draws dozens.
cfg.visible_figures = false;

%% Presets per stage and figure
cfg.mc_preset            = cfg.preset_name;         % memory capacity on THIS network
cfg.bursting_preset      = 'bursting_pairs';
cfg.sompolinsky_preset   = 'sompolinsky_pairs';
cfg.stf_preset           = 'single_neuron_stf';
% One unconnected neuron carrying the paper's c, SFA ladder and dual STD. It is
% a separate preset because handing it cfg.preset_name once built the whole
% 500-neuron network for a figure captioned "one unconnected neuron".
cfg.single_neuron_preset = 'single_neuron_dualStd';

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
F = add(F, 'fig_energy_landscape',            @fig_energy_landscape,            false, ...
        {'gammas', cfg.panelA_gammas});
F = add(F, 'fig_example_timeseries',          @fig_example_timeseries,          true, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_FI_curve',                    @fig_FI_curve,                    true, {});
F = add(F, 'fig_adaptation_methods',          @fig_adaptation_methods,          true, ...
        {'variant', 'sfa_std', 'preset_name', cfg.single_neuron_preset});
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
F = add(F, 'fig_EI_param_space',              @fig_EI_param_space,              true, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_EI_weights_param_space',      @fig_EI_weights_param_space,      true, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_sfa_EOC_allStd',              @fig_sfa_EOC_allStd,              true, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_memory_capacity',             @fig_memory_capacity,             true, {});
F = add(F, 'fig_memory_capacity_example',     @fig_memory_capacity_example,     true, {});
F = add(F, 'fig_eig_heatmap',                 @fig_eig_heatmap,                 false, {});
F = add(F, 'fig_dc_lle',                      @fig_dc_lle,                      false, {});
% The generated equation and conditions tables, an ordinary entry so its
% failures count in the headline number.
F = add(F, 'doc_tables',                      @fig_doc_tables,                  true, ...
        {'preset_name', cfg.preset_name});
cfg.figures = F;
end

%% ------------------------------------------------------------------------
function F = add(F, name, fn, in_paper, args)
F{end+1} = struct('name', name, 'fn', fn, 'in_paper', in_paper, 'args', {args});
end
