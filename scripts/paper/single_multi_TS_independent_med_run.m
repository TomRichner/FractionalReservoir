% SINGLE_MULTI_TS_INDEPENDENT_MED_RUN Single vs multi timescale, self-contained.
%
%   Open this file and press Run. setup_paths is called on the first line.
%
% Every setting comes from single_multi_TS_independent_med_config(), which states
% ALL of them itself -- it does not call paper_config. Read that one file and
% you know the whole run.
%
%   preset    celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25
%             (mu x1.5, SFA on E AND I, per-neuron S_c with sigma 0.1)
%   run_mode  medium
%   run_dir   data/single_multi_TS_independent_med
%   fig_root  figs/single_multi_TS_independent_med
%
% Memory capacity runs on THE SAME PRESET as the sweeps here, not on the
% separate mc_pairs_dualStd network the paper_config lineage uses. At 'medium'
% that is 15 paired trials with an exact sign-flip test -- see the config.
%
% The same experiment as single_multi_TS_independent_run at 'fast'; both roots
% are distinct, so the two never collide.
%
% RERUNNING REQUIRES DELETING data/single_multi_TS_independent_med FIRST.
% run_all_paper_analyses refuses a run directory that is not absent or empty.
%
% See also: single_multi_TS_independent_med_config, run_all_paper_analyses,
%           make_all_paper_figures

setup_paths();

cfg = single_multi_TS_independent_med_config();

run_dir = run_all_paper_analyses(cfg);
results = make_all_paper_figures(cfg);

fprintf('\n========================================================\n');
fprintf('SINGLE vs MULTI TIMESCALE (independent config, medium) RUN COMPLETE\n');
fprintf('  preset  : %s (%s)\n', cfg.preset_name, cfg.run_mode);
fprintf('  MC on   : %s\n', cfg.mc_preset);
fprintf('  run_dir : %s\n', run_dir);
fprintf('  figures : %d of %d succeeded\n', sum([results.ok]), numel(results));
fprintf('  to rerun: delete %s first\n', cfg.run_dir);
fprintf('========================================================\n');
