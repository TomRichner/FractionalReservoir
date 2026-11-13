%% test_mc_gain_sweep_job_files.m
% A job file run_mc_gain_sweep wrote EARLIER is still read and resumed by
% the code as it stands (FD-010). The stage's runs take hours and are
% resumed across sessions and code changes, so a change to the job struct,
% the job file name or the signature must not orphan the files on disk.
%
% Uses a real run under data/ (this checkout's, or the main checkout's when
% this is a git worktree) and SKIPS, with a line saying so, when there is
% none. The run's folder is only read; the resume check works on a copy.
%
% Checks, on the first job file of the run:
%   1. mc_gain_sweep_job_rows maps the job to one row per condition, every
%      column the job's own values.
%   2. read_mc_gain_sweep_jobs, on the run's folder, gives for that job the
%      rows read_mc_gain_sweep gives for its consolidated rows; and, when
%      the run has written its consolidated table and it holds the job, the
%      rows of that table.
%   3. RESUME: with the file copied into an empty output directory,
%      run_mc_gain_sweep asked for that one job finds it by name, accepts
%      its signature, runs nothing and consolidates the same rows.
%
% Asserts; prints one PASS or SKIP line. A few seconds, serial.
% Assumes setup_paths has run.

job_dir = find_job_dir();
if isempty(job_dir)
    fprintf('SKIP: test_mc_gain_sweep_job_files -- no gain-sweep run with job files under data/.\n');
    return
end
listing = dir(fullfile(job_dir, 'job_*.mat'));
job_file = fullfile(listing(1).folder, listing(1).name);
loaded = load(job_file, 'job');
job = loaded.job;
n_cond = numel(job.conditions);

%% 1. The job as table rows
rows = mc_gain_sweep_job_rows(job);
assert(height(rows) == n_cond, 'one row per condition');
assert(all(rows.gain == job.gain | (isnan(rows.gain) & isnan(job.gain))) && all(rows.bias == job.bias) ...
    && all(rows.seed == job.seed) && all(rows.seed_stim == job.seed_stim), 'the job''s own gain, bias and seeds');
assert(isequal(rows.condition, job.conditions(:)) && isequal(rows.message, job.message(:)));
assert(isequaln(rows.R2_delay, job.R2_delay));
per_condition = {'mc_total', 'horizon', 'lle_input', 'lle_noinput', 'rate_E', 'rate_I', ...
    'frac_pinned', 'frac_silent', 'seconds', 'failed'};
for f = per_condition
    assert(isequaln(rows.(f{1}), job.(f{1})(:)), 'column %s is not the job''s', f{1});
end
assert(isequal(rows.Properties.VariableNames, [{'gain', 'bias', 'seed', 'seed_stim', 'condition', ...
    'mc_total', 'horizon', 'R2_delay'}, per_condition(3:end), {'message'}]), 'the columns, in order');

%% 2. The reader of a run in progress agrees with the one reader
is_job = @(T) T.gain == job.gain & T.bias == job.bias & T.seed == job.seed;
key = {'gain', 'bias', 'seed', 'condition'};
before = dir(job_dir);
D_run = read_mc_gain_sweep_jobs(job_dir);
after = dir(job_dir);
kept = ismember({before.name}, {after.name});     % a run in progress may add files
[~, where] = ismember({before(kept).name}, {after.name});
assert(isequal([before(kept).datenum], [after(where).datenum]), 'the reader must not touch the job folder');
got = sortrows(D_run.T(is_job(D_run.T), :), key);
D_job = read_mc_gain_sweep(rows);
want = sortrows(D_job.T, key);
assert(height(got) == n_cond && isequaln(got, want), ...
    'read_mc_gain_sweep_jobs and read_mc_gain_sweep disagree on %s', listing(1).name);
results_file = fullfile(fileparts(job_dir), 'mc_gain_sweep_results.mat');
checked_consolidated = false;
if isfile(results_file)
    D_all = read_mc_gain_sweep(results_file);
    from_file = sortrows(D_all.T(is_job(D_all.T), :), key);
    if ~isempty(from_file)
        assert(isequaln(from_file, want), 'the consolidated table holds other rows for %s', listing(1).name);
        checked_consolidated = true;
    end
end

%% 3. Resume from a copy
tmp = tempname;
mkdir(fullfile(tmp, 'mc_gain_sweep', 'jobs'));
cleanup = onCleanup(@() rmdir(tmp, 's'));
copyfile(job_file, fullfile(tmp, 'mc_gain_sweep', 'jobs', listing(1).name));
sig = job.signature;
[~, f] = run_mc_gain_sweep('preset_name', sig.preset_name, 'protocol', sig.protocol, ...
    'ode_solver', sig.ode_solver, 'model_args', sig.model_args, 'conditions', sig.conditions, ...
    'lya_K', sig.lya(1), 'lya_T_sec', sig.lya(2), 'lya_acc_sec', sig.lya(3), ...
    'gain_levels', job.gain, 'bias_levels', job.bias, 'seeds', job.seed, ...
    'output_dir', tmp, 'retry_failed', false, 'use_parallel', false, 'verbose', 'near-none');
resumed = load(f);
assert(resumed.settings.n_jobs == 1 && resumed.settings.n_run == 0, ...
    'the stage did not take %s for its job: it ran %d', listing(1).name, resumed.settings.n_run);
assert(isequaln(resumed.T, rows), 'the resumed run consolidated other rows');

consolidated_txt = 'not yet written, or without this job';
if checked_consolidated; consolidated_txt = 'agrees'; end
fprintf(['PASS: test_mc_gain_sweep_job_files -- %s of %s: rows, both readers, resume ' ...
    '(consolidated table: %s).\n'], listing(1).name, job_dir, consolidated_txt);

%% ------------------------------------------------------------------------
function job_dir = find_job_dir()
% The jobs folder of a gain-sweep run under data/, the main block C run
% first. In a git worktree (<repo>/.claude/worktrees/<name>) data/ is the
% main checkout's: it is gitignored, so a worktree has none of its own.
job_dir = '';
root = fileparts(which('setup_paths'));
roots = {root};
in_worktree = regexp(root, '^(.*)[\\/]\.claude[\\/]worktrees[\\/][^\\/]+$', 'tokens', 'once');
if ~isempty(in_worktree); roots{end+1} = in_worktree{1}; end
for r = 1:numel(roots)
    runs = dir(fullfile(roots{r}, 'data', 'mc_gain_sweep_*'));
    names = sort({runs([runs.isdir]).name});
    names = [names(strcmp(names, 'mc_gain_sweep_main')), names(~strcmp(names, 'mc_gain_sweep_main'))];
    for k = 1:numel(names)
        candidate = fullfile(roots{r}, 'data', names{k}, 'mc_gain_sweep', 'jobs');
        if ~isempty(dir(fullfile(candidate, 'job_*.mat')))
            job_dir = candidate;
            return
        end
    end
end
end
