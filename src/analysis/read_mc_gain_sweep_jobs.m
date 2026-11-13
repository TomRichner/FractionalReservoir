function [D, info] = read_mc_gain_sweep_jobs(folder)
% READ_MC_GAIN_SWEEP_JOBS Read a gain-by-bias run that is still in progress, from its job files (block C).
%
%   D = READ_MC_GAIN_SWEEP_JOBS(folder)
%   [D, info] = READ_MC_GAIN_SWEEP_JOBS(folder)
%
% folder is a run's jobs folder, its mc_gain_sweep folder, or the run
% directory above it. Every finished job_*.mat there is read and
% consolidated into the table run_mc_gain_sweep writes at its end (the same
% mapping, mc_gain_sweep_job_rows), then passed through read_mc_gain_sweep.
% So D has the shape every block C figure and mc_gain_bias_statistics take,
% and they can be drawn before the run ends.
%
% READ-ONLY: nothing in the folder is written, moved or deleted. A job
% still being saved is a *.tmp file and is not matched.
%
% D.gains, D.biases and D.seeds are the levels FOUND, not the run's planned
% grid; a partial run's map has gaps (grey cells) where jobs have not
% finished. D.settings carries what the job files hold: conditions,
% preset_name, delay_s (mc_delay_axis of the signature's protocol),
% preset_gain and preset_bias (mc_gain_sweep_preset_point of the
% signature's preset and model_args), partial = true and n_job_files.
%
% info: job_dir, n_files, n_read, n_unreadable (a file that fails to load
% is listed in info.unreadable and skipped with a warning), n_signatures
% (1; job files of different runs are an error here, as run_mc_gain_sweep
% refuses to mix them).
%
% See also: read_mc_gain_sweep, run_mc_gain_sweep, mc_gain_sweep_job_rows,
%           make_block_C_figures

arguments
    folder (1,:) char
end

job_dir = folder;
if isfolder(fullfile(folder, 'mc_gain_sweep', 'jobs'))
    job_dir = fullfile(folder, 'mc_gain_sweep', 'jobs');
elseif isfolder(fullfile(folder, 'jobs'))
    job_dir = fullfile(folder, 'jobs');
end
if ~isfolder(job_dir)
    error('read_mc_gain_sweep_jobs:NoFolder', 'No such folder: %s', job_dir);
end
listing = dir(fullfile(job_dir, 'job_*.mat'));
info = struct('job_dir', job_dir, 'n_files', numel(listing), 'n_read', 0, ...
    'n_unreadable', 0, 'unreadable', {{}}, 'n_signatures', 0);
if isempty(listing)
    error('read_mc_gain_sweep_jobs:NoJobs', 'No job_*.mat in %s (yet).', job_dir);
end

rows = cell(numel(listing), 1);
signatures = {};
cond_names = {};
for k = 1:numel(listing)
    f = fullfile(listing(k).folder, listing(k).name);
    try
        loaded = load(f, 'job');
        job = loaded.job;
    catch err
        info.unreadable{end+1} = f;
        warning('read_mc_gain_sweep_jobs:Unreadable', 'Skipped %s: %s', f, err.message);
        continue
    end
    rows{k} = mc_gain_sweep_job_rows(job);
    if ~any(cellfun(@(s) isequaln(s, job.signature), signatures))
        signatures{end+1} = job.signature; %#ok<AGROW>
    end
    if isempty(cond_names); cond_names = job.conditions; end
end
info.n_unreadable = numel(info.unreadable);
info.n_read = info.n_files - info.n_unreadable;
info.n_signatures = numel(signatures);
if info.n_read == 0
    error('read_mc_gain_sweep_jobs:NoJobs', 'No readable job file in %s.', job_dir);
end
if info.n_signatures > 1
    error('read_mc_gain_sweep_jobs:MixedRuns', ...
        '%s holds job files of %d different runs (preset, protocol, integrator or conditions differ).', ...
        job_dir, info.n_signatures);
end
T = vertcat(rows{:});

signature = signatures{1};
settings = struct();
settings.conditions  = cond_names;
settings.preset_name = signature.preset_name;
settings.delay_s     = mc_delay_axis(signature.protocol);
[settings.preset_gain, settings.preset_bias] = ...
    mc_gain_sweep_preset_point(signature.preset_name, signature.model_args);
settings.partial     = true;
settings.n_job_files = info.n_read;

D = read_mc_gain_sweep(struct('T', T, 'settings', settings));
D.source = sprintf('%s (partial: %d job files)', job_dir, info.n_read);
end
