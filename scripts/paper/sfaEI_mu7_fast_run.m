% SFAEI_MU7_FAST_RUN The mu = 7 network (SFA on E and I, per-neuron S_c, SFA spread 0.25, no external input), fast, every analysis, MC noise-free on sra1, local-Lyapunov staircase stage.
%
%   Open this file and press Run. setup_paths is called on the first line.
%
% Every setting comes from sfaEI_mu7_fast_config(), which states ALL of them
% itself -- it does not call paper_config. Read that one file and you know the
% whole run.
%
%   preset    celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0p025_dualStd_3cond_mu7
%             (mu_tilde_relative 7, SFA on E AND I, per-neuron S_c with sigma
%             0.1, tau_a_spread 0.25, Wiener process only -- no steps)
%   run_mode  fast
%   run_dir   data/sfaEI_mu7_fast
%   fig_root  figs/sfaEI_mu7_fast
%
% Memory capacity runs on THE SAME NETWORK with the noise off (the chained
% ..._noise0_... twin) on sra1; the local-Lyapunov stage on the ..._steps5s_...
% twin (a new random step every 5 s).
%
% RERUNNING REQUIRES DELETING data/sfaEI_mu7_fast FIRST.
% run_all_paper_analyses refuses a run directory that is not absent or empty.
%
% See also: sfaEI_mu7_fast_config, run_all_paper_analyses, make_all_paper_figures

setup_paths();

cfg = sfaEI_mu7_fast_config();

run_dir = run_all_paper_analyses(cfg);
results = make_all_paper_figures(cfg);

fprintf('\n========================================================\n');
fprintf('SINGLE vs MULTI TIMESCALE (sfaEI mu 7, fast) RUN COMPLETE\n');
fprintf('  preset  : %s (%s)\n', cfg.preset_name, cfg.run_mode);
fprintf('  MC on   : %s\n', cfg.mc_preset);
fprintf('  run_dir : %s\n', run_dir);
fprintf('  figures : %d of %d succeeded\n', sum([results.ok]), numel(results));
fprintf('  to rerun: delete %s first\n', cfg.run_dir);
fprintf('========================================================\n');
