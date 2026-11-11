% TEST_POOL_MODULE The pool module (src/util/pool_*.m, issue FD-017).
%
% All parallel-pool handling goes through five functions and a session
% allocation held in an environment variable of the MATLAB process. This
% checks the behaviour the rest of the code and the multi-agent working
% arrangement rely on:
%
%   1. pool_set(N) opens exactly N workers on the module's profile, records
%      the allocation and the label, and writes the registry entry
%   2. pool_status lists this session with the allocated count and the count
%      of workers the operating system reports
%   3. pool_ensure reopens a pool that went down, at the allocated size
%   4. pool_restart gives a NEW pool of the same size
%   5. pool_set to another size replaces the pool
%   6. pool_set(0) closes it; pool_ensure is then empty and silent
%   7. with the allocation UNSET, pool_ensure returns empty and warns once
%      per call site, opens nothing, and both entry points refuse with their
%      own error identifier before creating anything
%   8. pool_release removes the registry entry and leaves the allocation 0
%   9. a registry entry of a process that no longer exists is pruned
%
% IT LEAVES THE SESSION AS IT FOUND IT: the allocation, label and priority
% option are put back, and a pool of the size that was up is reopened (none
% if none was). It does restart the session's pool, so do not run it beside
% something that is using that pool. It never holds more than 2 workers of
% its own.
%
% About 50 s, nearly all of it pools starting and stopping, so it is on the
% slow list of run_all_tests.
%
% See also: pool_set, pool_ensure, pool_restart, pool_release, pool_status

c = pool_config();
my_pid = c.pid;
entry_file = fullfile(c.registry_dir, sprintf('session_%d.json', my_pid));
saved = struct('workers', getenv(c.env_workers), 'label', getenv(c.env_label), ...
    'low', getenv(c.env_low_priority), 'had_entry', isfile(entry_file), 'n_up', 0);
p0 = gcp('nocreate');
if ~isempty(p0); saved.n_up = p0.NumWorkers; end
clear p0
quiet = 'near-none';
all_passed = true;
fprintf('test_pool_module: profile ''%s'', budget %d, registry %s\n', c.profile, c.budget, c.registry_dir);

try
    %% 1. pool_set(N) opens N
    fprintf('\n-- pool_set --\n');
    pool = pool_set(2, 'label', 'test_pool_module', 'verbose', quiet);
    up = gcp('nocreate');
    all_passed = check('pool_set(2) opens 2 workers', ~isempty(up) && up.NumWorkers == 2 && pool == up) && all_passed;
    all_passed = check('...on the module''s profile, not the default', strcmp(up.Cluster.Profile, c.profile)) && all_passed;
    c2 = pool_config();
    all_passed = check('...records the allocation and the label', ...
        c2.is_set && c2.allocation == 2 && strcmp(c2.label, 'test_pool_module')) && all_passed;
    E = pool_registry('read');
    mine = E([E.pid] == my_pid);
    all_passed = check('...writes this session''s registry entry', ...
        isscalar(mine) && mine.workers == 2 && strcmp(mine.label, 'test_pool_module')) && all_passed;

    %% 2. pool_status
    fprintf('\n-- pool_status --\n');
    T = pool_status(quiet);
    row = T(T.this, :);
    all_passed = check('lists this session once, allocated 2', height(row) == 1 && row.allocated == 2) && all_passed;
    if isnan(row.running)
        fprintf('  note  running workers cannot be counted from the operating system here; registry only\n');
    else
        all_passed = check('counts 2 running workers from the operating system', row.running == 2) && all_passed;
    end

    %% 3. pool_ensure reopens at the allocation
    fprintf('\n-- pool_ensure --\n');
    all_passed = check('returns the pool that is up', pool_ensure() == up) && all_passed;
    delete(gcp('nocreate'));
    pool = pool_ensure('verbose', quiet);
    all_passed = check('reopens a pool that went down, at the allocated size', ...
        ~isempty(pool) && pool.NumWorkers == 2) && all_passed;

    %% 4. pool_restart
    fprintf('\n-- pool_restart --\n');
    before = gcp('nocreate');
    pool = pool_restart(quiet);
    all_passed = check('gives a fresh pool', ~isvalid(before) && ~isempty(pool)) && all_passed;
    all_passed = check('...of the same size', ~isempty(pool) && pool.NumWorkers == 2) && all_passed;

    %% 5. another size replaces the pool
    fprintf('\n-- pool_set to another size --\n');
    pool = pool_set(1, 'verbose', quiet);
    up = gcp('nocreate');
    all_passed = check('pool_set(1) replaces the 2-worker pool with 1 worker', ...
        ~isempty(up) && up.NumWorkers == 1 && pool == up && pool_config().allocation == 1) && all_passed;
    T = pool_status(quiet);
    row = T(T.this, :);
    all_passed = check('pool_status follows: allocated 1', row.allocated == 1 && (isnan(row.running) || row.running == 1)) && all_passed;

    %% 6. pool_set(0)
    fprintf('\n-- pool_set(0) --\n');
    pool_set(0, 'verbose', quiet);
    c0 = pool_config();
    all_passed = check('closes the pool and records allocation 0', ...
        isempty(gcp('nocreate')) && c0.is_set && c0.allocation == 0) && all_passed;
    lastwarn('');
    pool = pool_ensure();
    [~, wid] = lastwarn;
    all_passed = check('pool_ensure is then empty, without a warning', isempty(pool) && isempty(wid)) && all_passed;
    T = pool_status(quiet);
    row = T(T.this, :);
    all_passed = check('pool_status: allocated 0, none running', ...
        height(row) == 1 && row.allocated == 0 && (isnan(row.running) || row.running == 0)) && all_passed;

    %% 7. no allocation set
    fprintf('\n-- allocation unset --\n');
    setenv(c.env_workers, '');
    clear pool_ensure                   % forget which call sites were warned
    wids = cell(1, 2);
    pools = cell(1, 2);
    warn_state = warning('off', 'backtrace');
    for k = 1:2                         % ONE call site, reached twice
        lastwarn('');
        pools{k} = pool_ensure();
        [~, wids{k}] = lastwarn;
    end
    warning(warn_state);
    all_passed = check('pool_ensure returns empty', isempty(pools{1}) && isempty(pools{2})) && all_passed;
    all_passed = check('...warns pool_ensure:NoAllocation', strcmp(wids{1}, 'pool_ensure:NoAllocation')) && all_passed;
    all_passed = check('...once per call site', isempty(wids{2})) && all_passed;
    all_passed = check('...and opens nothing', isempty(gcp('nocreate'))) && all_passed;

    root = fileparts(which('setup_paths'));
    probe_run = fullfile(root, 'data', '__pool_probe__');
    probe_fig = fullfile(root, 'figs', '__pool_probe__');
    all_passed = check('run_all_paper_analyses refuses (NoPoolAllocation)', ...
        throws_id(@() run_all_paper_analyses(paper_config('run_mode', 'fast', 'run_dir', 'data/__pool_probe__')), ...
        'run_all_paper_analyses:NoPoolAllocation')) && all_passed;
    all_passed = check('...and creates no run directory', ~isfolder(probe_run)) && all_passed;
    all_passed = check('make_all_paper_figures refuses (NoPoolAllocation)', ...
        throws_id(@() make_all_paper_figures(paper_config('fig_root', 'figs/__pool_probe__')), ...
        'make_all_paper_figures:NoPoolAllocation')) && all_passed;
    all_passed = check('...and creates no figure directory', ~isfolder(probe_fig)) && all_passed;
    if isfolder(probe_run); rmdir(probe_run, 's'); end
    if isfolder(probe_fig); rmdir(probe_fig, 's'); end

    %% 8. pool_release
    fprintf('\n-- pool_release --\n');
    pool_set(0, 'verbose', quiet);
    had = isfile(entry_file);
    pool_release(quiet);
    cr = pool_config();
    all_passed = check('removes this session''s registry entry', had && ~isfile(entry_file)) && all_passed;
    all_passed = check('...and leaves the allocation at 0', cr.is_set && cr.allocation == 0) && all_passed;
    T = pool_status(quiet);
    all_passed = check('pool_status no longer lists this session', ~any(T.this)) && all_passed;

    %% 9. a stale entry is pruned
    fprintf('\n-- stale registry entry --\n');
    P = pool_registry('processes');
    if ispc && ~P.ok
        fprintf('  note  the operating system could not be asked for its processes; pruning not checked\n');
    else
        dead_pid = 4194300;             % above any process id in use; checked below
        if ispc
            while any(P.pid == dead_pid); dead_pid = dead_pid - 4; end
        else
            while system(sprintf('kill -0 %d 2>/dev/null', dead_pid)) == 0; dead_pid = dead_pid - 1; end
        end
        stale_file = fullfile(c.registry_dir, sprintf('session_%d.json', dead_pid));
        if ~isfolder(c.registry_dir); mkdir(c.registry_dir); end
        fid = fopen(stale_file, 'w');
        fprintf(fid, '%s\n', jsonencode(struct('pid', dead_pid, 'label', 'gone', 'workers', 7, ...
            'started', '2026-01-01 00:00:00', 'host', 'nowhere')));
        fclose(fid);
        E = pool_registry('read');
        seen = any([E.pid] == dead_pid);
        T = pool_status(quiet);
        all_passed = check('an entry of a process that no longer exists is read, then pruned', ...
            seen && ~isfile(stale_file) && ~any(T.pid == dead_pid)) && all_passed;
        if isfile(stale_file); delete(stale_file); end
    end
catch test_error
    restore_session(c, saved);
    rethrow(test_error);
end

%% Leave the session as it was found
restore_session(c, saved);
now_up = gcp('nocreate');
n_now = 0;
if ~isempty(now_up); n_now = now_up.NumWorkers; end
all_passed = check('the session''s allocation and pool are back as they were found', ...
    strcmp(getenv(c.env_workers), saved.workers) && strcmp(getenv(c.env_label), saved.label) && ...
    strcmp(getenv(c.env_low_priority), saved.low) && n_now == saved.n_up && ...
    isfile(entry_file) == saved.had_entry) && all_passed;

if all_passed
    fprintf('\ntest_pool_module: ALL TESTS PASSED\n');
else
    fprintf(2, '\ntest_pool_module: FAILURES ABOVE\n');
end

%% ------------------------------------------------------------------------
function restore_session(c, saved)
% Put back the three environment variables, the pool size and the registry
% entry. pool_set is the only way to open a pool, and it records what it
% opens, so the variables are restored AFTER it.
up = gcp('nocreate');
n_up = 0;
if ~isempty(up); n_up = up.NumWorkers; end
if saved.n_up > 0
    if n_up ~= saved.n_up
        pool_set(saved.n_up, 'verbose', 'near-none');
    end
elseif n_up > 0
    delete(up);
end
setenv(c.env_workers, saved.workers);
setenv(c.env_label, saved.label);
setenv(c.env_low_priority, saved.low);
if saved.had_entry
    pool_registry('write');
else
    pool_registry('remove');
end
clear pool_ensure
end

function ok = check(label, ok)
if ok; tag = 'PASS'; else; tag = 'FAIL'; end
fprintf('  %s  %s\n', tag, label);
end

function tf = throws_id(fn, want)
try
    fn();
    tf = false;
catch ME
    tf = strcmp(ME.identifier, want);
    if ~tf; fprintf(2, '     got %s: %s\n', ME.identifier, ME.message); end
end
end
