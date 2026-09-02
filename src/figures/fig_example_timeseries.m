function out = fig_example_timeseries(cfg)
% FIG_EXAMPLE_TIMESERIES Example SRNN time series at the paper's operating point.
%
%   out = FIG_EXAMPLE_TIMESERIES()
%   out = FIG_EXAMPLE_TIMESERIES('preset_name', p, 'out_dir', d)
%
% One realization of the model the rest of the paper analyses: the stimulus, the
% dendritic state x, the firing rate r, the synaptic output, and the adaptation
% and depression variables, under SFA + STD.
%
% PORTED. This was SRNNModel2 at class defaults with the LOGISTIC nonlinearity,
% n = 300, level_of_chaos = 1 and no noise -- while every sweep figure in the
% paper is SRNNCellTypePairs, piecewise, n = 500, S_c = 0.2, sigma_u_noise =
% 0.025. It was the "here is what the model does" figure showing a different
% model from the one being characterised. It now takes the same preset as the
% sweeps, so the example and the statistics describe one network.
%
% See also: paper_config, srnn_param_preset, SRNNCellTypePairs

arguments
    cfg.preset_name (1,:) char    = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    cfg.out_dir     (1,:) char    = ''
    cfg.condition   (1,:) char    = 'sfa_and_std'
    cfg.rng_seeds   (1,2) double  = [1 2]
    cfg.T_range     (1,2) double  = [0 20]
    cfg.save        (1,1) logical = true
    cfg.visible     (1,1) logical = true
    cfg.run_dir     (1,:) char    = ''      % unused; accepted for a uniform call
end

setup_paths();
out_dir = default_out_dir(cfg.out_dir, mfilename('fullpath'));

%% Build the model from the preset, under one adaptation condition
% build_from_preset merges preset + condition, selects the integrator the
% preset's noise requires, and builds. See its header for the three traps it
% exists to close (integrator selection, flat argument expansion, and the fact
% that adaptation counts live in the CONDITION, never in a preset).
model = build_from_preset(cfg.preset_name, cfg.condition, ...
    'rng_seeds', cfg.rng_seeds, 'T_range', cfg.T_range, ...
    'fs', 400, 'lya_method', 'none');

model.run();

[fig_handle, ~] = model.plot();
if ~cfg.visible; set(fig_handle, 'Visible', 'off'); end

%% Save
out = struct('figs', fig_handle, 'files', {{}}, 'source', 'simulated inline');
if cfg.save
    save_figure_stable(out_dir, 'fig_example_timeseries', fig_handle);
    out.files = existing_outputs(out_dir, 'fig_example_timeseries');

end
end




