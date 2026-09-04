function cfg = sing_multi_TS_50percentStronger_std3_config(opts)
% SING_MULTI_TS_50PERCENTSTRONGER_STD3_CONFIG Mu x1.5, three STD timescales.
%
%   run_dir = run_all_paper_analyses(sing_multi_TS_50percentStronger_std3_config());
%   results = make_all_paper_figures(sing_multi_TS_50percentStronger_std3_config());
%
% sing_multi_TS_50percentStronger_config with the preset swapped for
% celltype_pairs_Sc0p2_noise0p025_tripleStd_3cond_mu8p25: the same mu x1.5
% network, but depression on THREE timescales (tau_rel [0.25 0.5 1], tau_rec
% 4*tau_rel) rather than two, so the full regime is sfa3_std3. Note the
% depression ratio rho = tau_rel/tau_rec is 0.25 here against 0.125 in the
% two-timescale presets -- see the preset for what that changes.
%
% Both roots are fixed and named for this experiment. Delete data/... before
% rerunning: run_all_paper_analyses refuses a run directory that is not absent
% or empty.
%
% See also: sing_multi_TS_50percentStronger_config,
%           sing_multi_TS_50percentStronger_std3_run, paper_config,
%           srnn_param_preset

arguments
    opts.preset_name (1,:) char = 'celltype_pairs_Sc0p2_noise0p025_tripleStd_3cond_mu8p25'
    opts.run_mode    (1,:) char = 'medium'
    opts.run_dir     (1,:) char = 'data/sing_multi_TS_50percentStronger_std3'
    opts.fig_root    (1,:) char = 'figs/sing_multi_TS_50percentStronger_std3'
    opts.visible_figures (1,1) logical = false
end

cfg = paper_config( ...
    'preset_name', opts.preset_name, ...
    'run_mode',    opts.run_mode, ...
    'run_dir',     opts.run_dir, ...
    'fig_root',    opts.fig_root, ...
    'visible_figures', opts.visible_figures);
end
