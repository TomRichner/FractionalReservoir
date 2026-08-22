function cfg = paper_config(opts)
% PAPER_CONFIG The single place that names what the paper is built from.
%
%   cfg = PAPER_CONFIG()
%   cfg = PAPER_CONFIG('run_mode', 'fast')
%   cfg = PAPER_CONFIG('preset_name', 'celltype_pairs_uniform_std_n500')
%
% THIS IS THE ONE FILE TO EDIT. Change the preset here and both master scripts
% follow: run_all_paper_analyses sweeps that network, make_all_paper_figures
% resolves the matching run and rebuilds every figure from it.
%
% TWO ORTHOGONAL KNOBS, as everywhere else in this project:
%   preset_name  WHICH NETWORK   (src/srnn_param_preset.m)
%   run_mode     HOW MUCH COMPUTE (scripts/run_all_analyses/analysis_run_config.m)
%
% FIGURE PRESET OVERRIDES. Four figures are DELIBERATELY a different network
% from the paper's operating point, and each names its own preset below. They
% are not oversights: a single-neuron mechanism cartoon, a Sompolinsky
% reproduction, a small bursting network and the memory-capacity reservoir are
% all making points that the 500-neuron production network cannot make. Every
% other figure inherits cfg.preset_name.
%
% See also: run_all_paper_analyses, make_all_paper_figures, srnn_param_preset

arguments
    opts.preset_name (1,:) char = 'celltype_pairs_Sc0p2_noise0p025_dualStd'
    opts.run_mode    (1,:) char = 'production'
    opts.run_dir     (1,:) char = ''      % '' -> resolve by preset at figure time
end

cfg = struct();
cfg.preset_name = opts.preset_name;
cfg.run_mode    = opts.run_mode;
cfg.run_dir     = opts.run_dir;

% Presets for the figures that are deliberately different networks.
cfg.mc_preset          = 'mc_esn';              % SRNN_ESN_reservoir; see below
cfg.bursting_preset    = 'bursting_pairs';
cfg.sompolinsky_preset = 'sompolinsky_pairs';
cfg.stf_preset         = 'single_neuron_stf';

% The gains shared by the two halves of figure 1 panel A. They MUST match, so
% they live here rather than as a literal in each figure.
cfg.panelA_gammas = [0.9, 1.6, 2.5];

%% The figure registry
% Order is the order make_all_paper_figures runs them, and is roughly the order
% they appear in the manuscript. `in_paper` marks the ones currently used; the
% rest are kept working and regenerated, per TR, but are not in the manuscript.
%
% Every entry is  {name, handle, extra-args}  where the extra args are appended
% to the standard cfg fields (run_dir, out_dir, save, visible).
F = {};
F = add(F, 'fig_introductory_concepts',      @fig_introductory_concepts,      true, ...
        {'preset_name', cfg.sompolinsky_preset, 'gammas', cfg.panelA_gammas});
F = add(F, 'fig_energy_landscape',           @fig_energy_landscape,           false, ...
        {'gammas', cfg.panelA_gammas});
F = add(F, 'fig_example_timeseries',         @fig_example_timeseries,         true, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_FI_curve',                   @fig_FI_curve,                   true, {});
F = add(F, 'fig_adaptation_methods',         @fig_adaptation_methods,         true, ...
        {'variant', 'sfa_std', 'preset_name', cfg.preset_name});
F = add(F, 'fig_adaptation_methods_stf',     @fig_adaptation_methods,         false, ...
        {'variant', 'sfa_std_stf', 'preset_name', cfg.stf_preset});
F = add(F, 'fig_SFA_steady_state',           @fig_SFA_steady_state,           false, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_STD_steady_state',           @fig_STD_steady_state,           false, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_stim_engages_adaptation',    @fig_stim_engages_adaptation,    true, ...
        {'preset_name', cfg.bursting_preset});
F = add(F, 'fig_sensitivity_analysis_allStd', @fig_sensitivity_analysis_allStd, true, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_sensitivity_medians',        @fig_sensitivity_medians,        false, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_param_space_allStd',         @fig_param_space_allStd,         true, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_EI_param_space',             @fig_EI_param_space,             true, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_sfa_EOC_allStd',             @fig_sfa_EOC_allStd,             true, ...
        {'preset_name', cfg.preset_name});
F = add(F, 'fig_memory_capacity',            @fig_memory_capacity,            true, {});
F = add(F, 'fig_memory_capacity_example',    @fig_memory_capacity_example,    true, {});
F = add(F, 'fig_eig_heatmap',                @fig_eig_heatmap,                false, {});
cfg.figures = F;
end

%% ------------------------------------------------------------------------
function F = add(F, name, fn, in_paper, args)
F{end+1} = struct('name', name, 'fn', fn, 'in_paper', in_paper, 'args', {args});
end
