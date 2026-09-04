% SING_MULTI_TS_50PERCENTSTRONGER_RUN Single vs multi timescale at mu x1.5.
%
%   Open this file and press Run. setup_paths is called on the first line.
%
% Every setting comes from sing_multi_TS_50percentStronger_config(); edit that.
%
%   preset    celltype_pairs_Sc0p2_noise0p025_dualStd_3cond_mu8p25
%   run_mode  medium
%   run_dir   data/sing_multi_TS_50percentStronger
%   fig_root  figs/sing_multi_TS_50percentStronger
%
% The network is the paper's 3-condition one with mu_tilde_relative 50%
% stronger on every route (5.5 -> 8.25); regimes and everything else are
% unchanged, so results compare directly against single_multi_TS_run.
%
% RERUNNING REQUIRES DELETING data/sing_multi_TS_50percentStronger FIRST.
% run_all_paper_analyses refuses a run directory that is not absent or empty.
%
% run_all_paper_analyses wraps each stage and make_all_paper_figures reports
% per figure, counting success as files on disk -- read the summary tables at
% the end rather than assuming "finished" means "worked".
%
% See also: sing_multi_TS_50percentStronger_config, single_multi_TS_run,
%           run_all_paper_analyses, make_all_paper_figures

setup_paths();

cfg = sing_multi_TS_50percentStronger_config();

run_dir = run_all_paper_analyses(cfg);
results = make_all_paper_figures(cfg);

fprintf('\n========================================================\n');
fprintf('SINGLE vs MULTI TIMESCALE (mu x1.5) RUN COMPLETE\n');
fprintf('  preset  : %s (%s)\n', cfg.preset_name, cfg.run_mode);
fprintf('  run_dir : %s\n', run_dir);
fprintf('  figures : %d of %d succeeded\n', sum([results.ok]), numel(results));
fprintf('  to rerun: delete %s first\n', cfg.run_dir);
fprintf('========================================================\n');
