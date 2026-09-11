% NUMERICS_VERIFICATION_TRIALS_RUN Run the numerics stage on many seeds, in parallel, and its figures.
%
%   Open this file and press Run. setup_paths is called on the first line.
%
% Every setting comes from numerics_verification_trials_config(). The stage
% is called directly (run_all_paper_analyses would also run the sweeps), with
% the config's two trial counts and worker count, then make_all_paper_figures gets the three-entry
% registry.
%
%   preset    celltype_pairs_sfaEI_Sc0p2sig0p1_noise0p025_dualStd_3cond_mu8p25
%   run_mode  medium, 5 reshoot / 25 LLE seeds, 12 workers  (roughly 40 min)
%   run_dir   data/numerics_verification_trials
%   fig_root  figs/numerics_verification_trials
%
% Can be rerun in place.
%
% See also: numerics_verification_trials_config, run_numerics_verification,
%           fig_numerics_verification

setup_paths();

cfg = numerics_verification_trials_config();
project_root = fileparts(which('setup_paths'));
run_dir = fullfile(project_root, cfg.run_dir);

mat_file = run_numerics_verification( ...
    'preset_name', cfg.preset_name, ...
    'run_mode',    cfg.run_mode, ...
    'n_trials_reshoot', cfg.n_trials_reshoot, ...
    'n_trials_lle',     cfg.n_trials_lle, ...
    'n_workers',        cfg.n_workers, ...
    'out_dir',     fullfile(run_dir, 'numerics_verification'));

results = make_all_paper_figures(cfg);

fprintf('\n========================================================\n');
fprintf('NUMERICS VERIFICATION TRIALS RUN COMPLETE\n');
fprintf('  preset  : %s (%s, %d reshoot / %d LLE seeds)\n', cfg.preset_name, cfg.run_mode, ...
    cfg.n_trials_reshoot, cfg.n_trials_lle);
fprintf('  data    : %s\n', mat_file);
fprintf('  figures : %d of %d succeeded\n', sum([results.ok]), numel(results));
fprintf('========================================================\n');
