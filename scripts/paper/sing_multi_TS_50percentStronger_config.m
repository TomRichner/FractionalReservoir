function cfg = sing_multi_TS_50percentStronger_config(opts)
% SING_MULTI_TS_50PERCENTSTRONGER_CONFIG Single vs multi timescale, mu x1.5.
%
%   run_dir = run_all_paper_analyses(sing_multi_TS_50percentStronger_config());
%   results = make_all_paper_figures(sing_multi_TS_50percentStronger_config());
%
% single_multi_TS_config with ONE change: the preset is
% celltype_pairs_Sc0p2_noise0p025_dualStd_3cond_mu8p25, which is the paper's
% 3-condition network with the mean connectivity 50% stronger on all four
% routes (mu_tilde_relative 5.5 -> 8.25) and nothing else different. Same
% regimes, same run mode, so a run here lines up against a single_multi_TS run
% condition for condition.
%
% Both roots are fixed and named for this experiment, so it cannot touch the
% paper's output or the single_multi_TS runs. Delete data/... before rerunning:
% run_all_paper_analyses refuses a run directory that is not absent or empty.
%
% See also: single_multi_TS_config, sing_multi_TS_50percentStronger_run,
%           paper_config, srnn_param_preset

arguments
    opts.preset_name (1,:) char = 'celltype_pairs_Sc0p2_noise0p025_dualStd_3cond_mu8p25'
    opts.run_mode    (1,:) char = 'medium'
    opts.run_dir     (1,:) char = 'data/sing_multi_TS_50percentStronger'
    opts.fig_root    (1,:) char = 'figs/sing_multi_TS_50percentStronger'
    opts.visible_figures (1,1) logical = false
end

cfg = paper_config( ...
    'preset_name', opts.preset_name, ...
    'run_mode',    opts.run_mode, ...
    'run_dir',     opts.run_dir, ...
    'fig_root',    opts.fig_root, ...
    'visible_figures', opts.visible_figures);
end
