function T = mc_gain_sweep_job_rows(job)
% MC_GAIN_SWEEP_JOB_ROWS One job of the gain-by-bias sweep as rows of its results table.
%
%   T = MC_GAIN_SWEEP_JOB_ROWS(job)
%
%   job   the struct run_mc_gain_sweep saves in a jobs/job_*.mat file: one
%         (gain, bias, seed) with a value per condition
%   T     one row per condition, with the columns of the stage's
%         consolidated table (run_mc_gain_sweep's header lists them)
%
% THE ONE MAPPING from a job file to table rows. The stage consolidates its
% finished jobs with it, and read_mc_gain_sweep_jobs reads a run still in
% progress with it, so the two cannot disagree on a column. A field added
% to the job struct is added here and nowhere else.
%
% See also: run_mc_gain_sweep, read_mc_gain_sweep_jobs

arguments
    job (1,1) struct
end

n = numel(job.conditions);
per_job = @(v) repmat(v, n, 1);
T = table(per_job(job.gain), per_job(job.bias), per_job(job.seed), per_job(job.seed_stim), ...
    job.conditions(:), job.mc_total(:), job.horizon(:), job.R2_delay, ...
    job.lle_input(:), job.lle_noinput(:), job.rate_E(:), job.rate_I(:), ...
    job.frac_pinned(:), job.frac_silent(:), job.seconds(:), job.failed(:), job.message(:), ...
    'VariableNames', {'gain', 'bias', 'seed', 'seed_stim', 'condition', 'mc_total', 'horizon', ...
    'R2_delay', 'lle_input', 'lle_noinput', 'rate_E', 'rate_I', 'frac_pinned', 'frac_silent', ...
    'seconds', 'failed', 'message'});
end
