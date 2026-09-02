% FAST4_RUN Run the whole pipeline fast, on the 4-condition network.
%
%   Open this file and press Run. Nothing else is needed -- not even a prior
%   setup_paths, which the first line does.
%
% The smoke-test twin of reproduce_paper_run: same two entry points, same
% figure registry, but every setting comes from fast4_config() instead of
% paper_config(). **fast4_config.m is the file to edit** if you want a
% different preset or run mode here.
%
%   preset    celltype_pairs_Sc0p2_noise0p025_dualStd_4cond
%   run_mode  fast
%   run_dir   data/fast_4
%   fig_root  figs/fast_4
%
% WHAT IT IS FOR. Exercising the pipeline end to end when you care whether the
% code runs, not what the numbers say -- after a refactor, after touching a
% model class, before committing to hours of 'production'. The figures it makes
% are structurally correct and far too coarse to read.
%
% IT CANNOT TOUCH THE PAPER'S OUTPUT. Both roots are named in fast4_config, so
% nothing here writes into figs/paper or into a dated run directory. That is the
% point of a separate config rather than paper_config('run_mode', 'fast'): no
% forgotten argument turns a smoke test into something that clobbers the real
% set.
%
% RERUNNING REQUIRES DELETING data/fast_4 FIRST.
%
%   rmdir(fullfile(fileparts(which('setup_paths')), 'data', 'fast_4'), 's')
%
% run_all_paper_analyses REFUSES a run directory that is not absent or empty,
% with TargetNotEmpty. That is deliberate, not an inconvenience: a reused run
% directory accumulates -- each sweep appends its own timestamped subfolder
% while the top-level manifest is replaced, leaving a manifest that describes
% one run in a directory holding several, and figures that match a sweep folder
% by parameter name then find two. figs/fast_4 needs no such care; it is
% overwritten in place.
%
% A HANDOFF, NOT A LOOKUP -- the one behavioural difference from
% reproduce_paper_run worth knowing. There, run_dir is empty, so the figures
% resolve the newest run matching the preset and could in principle follow a
% different run that finished in between. Here run_dir is named, so both halves
% address the same directory by name and there is nothing to race.
%
% COST: roughly a quarter hour, most of it the sweeps. Four conditions rather
% than seven, so each sweep does a little over half the work of the paper's.
%
% ERROR ISOLATION. run_all_paper_analyses wraps each of its four stages, so one
% failure is reported and the queue continues. make_all_paper_figures reports
% per figure and counts success as files actually on disk. Read the two summary
% tables at the end; do not assume that "it finished" means "everything worked".
%
% See also: fast4_config, reproduce_paper_run, run_all_paper_analyses,
%           make_all_paper_figures

setup_paths();

cfg = fast4_config();

run_dir = run_all_paper_analyses(cfg);
results = make_all_paper_figures(cfg);

fprintf('\n========================================================\n');
fprintf('FAST4 RUN COMPLETE\n');
fprintf('  preset  : %s (%s)\n', cfg.preset_name, cfg.run_mode);
fprintf('  run_dir : %s\n', run_dir);
fprintf('  figures : %d of %d succeeded\n', sum([results.ok]), numel(results));
fprintf('  to rerun: delete %s first\n', cfg.run_dir);
fprintf('========================================================\n');
