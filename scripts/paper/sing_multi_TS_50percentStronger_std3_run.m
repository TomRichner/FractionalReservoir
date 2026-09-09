% SING_MULTI_TS_50PERCENTSTRONGER_STD3_RUN Mu x1.5 with three STD timescales.
%
%   Open this file and press Run. setup_paths is called on the first line.
%
% Every setting comes from sing_multi_TS_50percentStronger_std3_config().
%
%   preset    celltype_pairs_Sc0p2_noise0p025_tripleStd_3cond_mu8p25
%   run_mode  medium
%   run_dir   data/sing_multi_TS_50percentStronger_std3
%   fig_root  figs/sing_multi_TS_50percentStronger_std3
%
% Same network as sing_multi_TS_50percentStronger_run, with depression on three
% timescales (sfa3_std3) rather than two. That is the ONLY difference: rho is
% 0.125 in both, so their sfa1_std1 controls are identical.
%
% RERUNNING REQUIRES DELETING data/sing_multi_TS_50percentStronger_std3 FIRST.
%
% See also: sing_multi_TS_50percentStronger_std3_config,
%           sing_multi_TS_50percentStronger_run, run_all_paper_analyses,
%           make_all_paper_figures

setup_paths();

cfg = sing_multi_TS_50percentStronger_std3_config();

run_dir = run_all_paper_analyses(cfg);
results = make_all_paper_figures(cfg);

fprintf('\n========================================================\n');
fprintf('SINGLE vs MULTI TIMESCALE (mu x1.5, STD x3) RUN COMPLETE\n');
fprintf('  preset  : %s (%s)\n', cfg.preset_name, cfg.run_mode);
fprintf('  run_dir : %s\n', run_dir);
fprintf('  figures : %d of %d succeeded\n', sum([results.ok]), numel(results));
fprintf('  to rerun: delete %s first\n', cfg.run_dir);
fprintf('========================================================\n');
