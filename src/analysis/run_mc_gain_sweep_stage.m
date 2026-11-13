function detail = run_mc_gain_sweep_stage(cfg, run_dir, plan)
% RUN_MC_GAIN_SWEEP_STAGE A bundle's mc_gain_sweep stage: the main sweep, then each variant network's.
%
%   detail = RUN_MC_GAIN_SWEEP_STAGE(cfg, run_dir, plan)
%
%   cfg      the bundle's config (unused beyond its verbosity; the plan
%            already holds what the config asked for)
%   run_dir  the bundle's run directory
%   plan     from mc_gain_sweep_stage_plan(cfg)
%   detail   the result folders, '; '-separated
%
% The run half of the opt-in stage 'mc_gain_sweep' (paper_opt_in_stages).
% Every run of the plan is tried in turn, so one that fails does not cost
% the others their hours; the stage then throws with every failure listed,
% and run_all_paper_analyses records it as failed.
%
% See also: mc_gain_sweep_stage_plan, run_mc_gain_sweep, run_all_paper_analyses

arguments
    cfg     (1,1) struct %#ok<INUSA>  the uniform signature of an opt-in stage's run
    run_dir (1,:) char
    plan    (1,1) struct
end

runs = plan.runs;
out_dirs = {};
failures = {};
for k = 1:numel(runs)
    try
        out_dirs{end+1} = run_mc_gain_sweep(runs(k).args{:}, ...
            'output_dir', fullfile(run_dir, runs(k).folder)); %#ok<AGROW>
    catch err
        failures{end+1} = sprintf('%s: %s: %s', runs(k).label, err.identifier, err.message); %#ok<AGROW>
        vfail('  mc_gain_sweep run %s FAILED: %s\n', runs(k).label, err.message);
    end
end
if ~isempty(failures)
    error('run_mc_gain_sweep_stage:RunsFailed', '%d of %d run(s) failed:\n%s', ...
        numel(failures), numel(runs), strjoin(failures, '\n'));
end
detail = strjoin(out_dirs, '; ');
end
