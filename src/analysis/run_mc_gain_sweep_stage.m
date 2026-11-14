function detail = run_mc_gain_sweep_stage(cfg, run_dir, plan)
% RUN_MC_GAIN_SWEEP_STAGE A bundle's mc_gain_sweep stage: every grid and variant network of its plan, as one pool of jobs.
%
%   detail = RUN_MC_GAIN_SWEEP_STAGE(cfg, run_dir, plan)
%
%   cfg      the bundle's config (its verbosity; the plan already holds
%            what the config asked for)
%   run_dir  the bundle's run directory
%   plan     from mc_gain_sweep_stage_plan(cfg)
%   detail   the result folders, '; '-separated
%
% The run half of the opt-in stage 'mc_gain_sweep' (paper_opt_in_stages).
% The plan's runs (each grid, and each variant network of a grid) go to
% run_mc_gain_sweep_set together, each into its own folder of the run
% directory (plan.runs.folder). It prepares every run before any job
% computes, so a run that cannot start fails the stage in seconds, and runs
% the jobs of all of them in one parfor. A job that fails is recorded in its
% run's table and costs no other job; a stage that is interrupted keeps its
% finished job files, and the same call resumes it.
%
% See also: mc_gain_sweep_stage_plan, run_mc_gain_sweep_set, run_all_paper_analyses

arguments
    cfg     (1,1) struct
    run_dir (1,:) char
    plan    (1,1) struct
end

runs = plan.runs;
arg_lists = arrayfun(@(r) [r.args, {'output_dir', fullfile(run_dir, r.folder)}], runs, ...
    'UniformOutput', false);
vprintf(cfg.verbose, 'minimal', ...
    '  mc_gain_sweep: %d runs (%s), %d jobs, about %.0f worker-minutes (divide by the pool''s workers)\n', ...
    numel(runs), strjoin({runs.label}, ', '), plan.n_jobs, plan.worker_minutes);
out_dirs = run_mc_gain_sweep_set(arg_lists);
detail = strjoin(out_dirs, '; ');
end
