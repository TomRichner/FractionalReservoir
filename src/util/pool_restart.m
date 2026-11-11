function pool = pool_restart(verbose)
% POOL_RESTART Close this session's pool and open a fresh one of the same size.
%
%   pool = pool_restart()
%   pool = pool_restart(verbose)     % 'verbose' | 'minimal' | 'near-none'
%
% Replaces restart_parpool (deleted with FD-017). That function reopened the
% DEFAULT profile, so the size it came back with was whatever the machine's
% default happened to be: 10 or 14 workers when it was written, and ONE once
% the default became 'Just1', which silently put every sweep of
% run_all_analyses on a single worker.
%
% The size now comes from the session:
%   a pool is up        closed, and the same number of workers reopened. If
%                       that number is not the recorded allocation (the pool
%                       was opened by hand), it BECOMES the allocation: the
%                       session is holding that many, and the record should
%                       say so.
%   none is up          opened at the allocation (pool_ensure): nothing with
%                       allocation 0, and the serial warning with none set.
% An empty result means the caller's parfor loops run serially.
%
% WHY RESTART AT ALL (kept from restart_parpool). run_all_analyses runs
% three multi-hour analyses back to back against one long-lived pool. Worker
% processes that have been churning ~1 GB arrays for hours are a plausible
% source of the failures seen in the aug_13 overnight run, where the tau
% stage returned "Out of memory" on 88 of 195 runs and a worker crash dump
% appeared at the moment that stage began -- with system memory nowhere near
% exhausted (page file peaked at 7.8 GB against a 128 GB commit limit).
% That diagnosis is NOT established. A fresh pool per stage is cheap
% insurance either way: it costs ~10 s per stage and bounds whatever state
% the workers accumulate.
%
% See also: pool_set, pool_ensure, pool_release, pool_status, run_all_analyses

arguments
    verbose = 'minimal'
end

pool = [];
if isempty(which('parpool')); return; end
c = pool_config();
old = gcp('nocreate');
if ~isempty(old)
    n = old.NumWorkers;
    vprintf(verbose, 'verbose', '[pool_restart] shutting down existing pool (%d workers)...\n', n);
    delete(old);
    if ~c.is_set || c.allocation ~= n
        vprintf(verbose, 'minimal', ...
            '[pool_restart] the pool that was up had %d workers; recording that as this session''s allocation\n', n);
        setenv(c.env_workers, sprintf('%d', n));
    end
end

pool = pool_ensure('verbose', verbose);
if isempty(pool); return; end
pool_registry('write');

ram = '';
if usejava('jvm')
    try
        [~, sys] = memory;
        ram = sprintf(', client RAM available %.1f GB', sys.PhysicalMemory.Available / 2^30);
    catch
        % memory() is Windows-only; not worth failing a run over.
    end
end
vprintf(verbose, 'verbose', '[pool_restart] fresh pool: %d workers, %d thread(s), profile ''%s''%s\n', ...
    pool.NumWorkers, pool.NumThreads, pool.Cluster.Profile, ram);
end
