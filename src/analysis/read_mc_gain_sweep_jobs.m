function [D, info] = read_mc_gain_sweep_jobs(folder)
% READ_MC_GAIN_SWEEP_JOBS Read a gain-by-bias run that is still in progress, from its job files (block C).
%
%   D = READ_MC_GAIN_SWEEP_JOBS(folder)
%   [D, info] = READ_MC_GAIN_SWEEP_JOBS(folder)
%
% folder is a run's jobs\ folder, its mc_gain_sweep\ folder, or the run
% directory above it. Every finished job_*.mat there is read (read-only:
% nothing in the folder is written, moved or deleted; a job still being
% saved is a *.tmp file and is not matched) and consolidated into the table
% run_mc_gain_sweep writes at the end, one row per (gain, bias, seed,
% condition), then passed through read_mc_gain_sweep, the one reader. So
% D has exactly the shape every block C figure and mc_gain_bias_statistics
% take, and the figures can be drawn from the main run before it ends.
%
% D.gains, D.biases and D.seeds are the levels FOUND so far, not the run's
% planned grid; a partial run's map has gaps (grey cells) where jobs have
% not finished. settings carries what the job files hold: conditions,
% preset_name, delay_s (from the protocol in the job's signature, computed
% as run_mc_gain_sweep computes it), and preset_gain / preset_bias from the
% signature's model_args when they set level_of_chaos / S_c (else
% read_mc_gain_sweep's defaults: gain 1, the preset's S_c), plus partial =
% true and n_job_files.
%
% info: n_files, n_read, n_unreadable (a file that fails to load, listed in
% info.unreadable and skipped with a warning), n_signatures (should be 1;
% more than 1 is an error, as run_mc_gain_sweep refuses to mix them).
%
% WHY A COPY: run_mc_gain_sweep consolidates with its local function
% job_rows, which is not callable from outside, and that file is not to be
% edited while the main run uses it. job_rows below is the same code; if
% the stage's job struct gains a field, add it here too
% (test_mc_gain_bias_figures checks this reader against the stage's own
% consolidated table on a synthetic job folder).
%
% See also: read_mc_gain_sweep, run_mc_gain_sweep, make_block_C_figures

folder = char(folder);
job_dir = folder;
if isfolder(fullfile(folder, 'mc_gain_sweep', 'jobs'))
    job_dir = fullfile(folder, 'mc_gain_sweep', 'jobs');
elseif isfolder(fullfile(folder, 'jobs'))
    job_dir = fullfile(folder, 'jobs');
end
if ~isfolder(job_dir)
    error('read_mc_gain_sweep_jobs:NoFolder', 'No such folder: %s', job_dir);
end
L = dir(fullfile(job_dir, 'job_*.mat'));
info = struct('job_dir', job_dir, 'n_files', numel(L), 'n_read', 0, ...
    'n_unreadable', 0, 'unreadable', {{}}, 'n_signatures', 0);
if isempty(L)
    error('read_mc_gain_sweep_jobs:NoJobs', 'No job_*.mat in %s (yet).', job_dir);
end

rows = cell(numel(L), 1);
sig = {};
cond_names = {};
for k = 1:numel(L)
    f = fullfile(L(k).folder, L(k).name);
    try
        S = load(f, 'job');
        job = S.job;
    catch err
        info.n_unreadable = info.n_unreadable + 1;
        info.unreadable{end+1} = f;
        warning('read_mc_gain_sweep_jobs:Unreadable', 'Skipped %s: %s', f, err.message);
        continue
    end
    rows{k} = job_rows(job);
    if isempty(sig) || ~any(cellfun(@(s) isequaln(s, job.signature), sig))
        sig{end+1} = job.signature; %#ok<AGROW>
    end
    if isempty(cond_names); cond_names = job.conditions; end
end
info.n_read = info.n_files - info.n_unreadable;
info.n_signatures = numel(sig);
if info.n_read == 0
    error('read_mc_gain_sweep_jobs:NoJobs', 'No readable job file in %s.', job_dir);
end
if numel(sig) > 1
    error('read_mc_gain_sweep_jobs:MixedRuns', ...
        '%s holds job files of %d different runs (preset, protocol, integrator or conditions differ).', ...
        job_dir, numel(sig));
end
T = vertcat(rows{:});

sg = sig{1};
settings = struct();
settings.conditions  = cond_names;
settings.preset_name = sg.preset_name;
settings.partial     = true;
settings.n_job_files = info.n_read;
if isfield(sg, 'protocol')
    P = sg.protocol;
    try
        hold_len  = max(1, round(P.T_hold * P.fs));
        d_max_eff = max(1, floor(round(P.d_max_sec * P.fs) / hold_len));
        settings.delay_s = (1:d_max_eff) * (hold_len / P.fs);
    catch
    end
end
if isfield(sg, 'model_args') && ~isempty(sg.model_args)
    ma = sg.model_args;
    for i = 1:2:numel(ma) - 1
        if strcmp(ma{i}, 'level_of_chaos'); settings.preset_gain = ma{i + 1}; end
        if strcmp(ma{i}, 'S_c');            settings.preset_bias = ma{i + 1}; end
    end
end
D = read_mc_gain_sweep(struct('T', T, 'settings', settings));
D.source = sprintf('%s (partial: %d job files)', job_dir, info.n_read);
end

%% ------------------------------------------------------------------------
function T = job_rows(job)
% One table row per condition of a job. THE SAME CODE as run_mc_gain_sweep's
% local job_rows (see the header).
n = numel(job.conditions);
T = table(repmat(job.gain, n, 1), repmat(job.bias, n, 1), repmat(job.seed, n, 1), ...
    repmat(job.seed_stim, n, 1), job.conditions(:), job.mc_total(:), job.horizon(:), ...
    job.R2_delay, job.lle_input(:), job.lle_noinput(:), job.rate_E(:), job.rate_I(:), ...
    job.frac_pinned(:), job.frac_silent(:), job.seconds(:), job.failed(:), job.message(:), ...
    'VariableNames', {'gain', 'bias', 'seed', 'seed_stim', 'condition', 'mc_total', 'horizon', ...
    'R2_delay', 'lle_input', 'lle_noinput', 'rate_E', 'rate_I', 'frac_pinned', 'frac_silent', ...
    'seconds', 'failed', 'message'});
end
