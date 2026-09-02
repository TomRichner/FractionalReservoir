function pool = restart_parpool()
%RESTART_PARPOOL Tear down any existing parallel pool and start a fresh one.
%
%   pool = restart_parpool()
%
%   Uses the DEFAULT cluster profile (parallel.defaultProfile), so it honours
%   whatever the machine is configured for rather than hard-coding a worker
%   count -- this repo is used on more than one machine, and their profiles
%   differ (10 workers here, 14 elsewhere).
%
%   Why restart at all. run_all_analyses runs three multi-hour analyses back to
%   back against one long-lived pool. Worker processes that have been churning
%   ~1 GB arrays for hours are a plausible source of the failures seen in the
%   aug_13 overnight run, where the tau stage returned "Out of memory" on 88 of
%   195 runs and a worker crash dump appeared at the moment that stage began --
%   with system memory nowhere near exhausted (page file peaked at 7.8 GB
%   against a 128 GB commit limit).
%
%   That diagnosis is NOT established. A fresh pool per stage is cheap
%   insurance either way: it costs ~30 s per stage and bounds whatever state
%   the workers accumulate.
%
%   See also: run_all_analyses, parpool, parallel.defaultProfile

old = gcp('nocreate');
if ~isempty(old)
    fprintf('[restart_parpool] shutting down existing pool (%d workers)...\n', ...
        old.NumWorkers);
    delete(old);
end

profile_name = parallel.defaultProfile;
pool = parpool(profile_name);
fprintf('[restart_parpool] fresh pool: %d workers, %d thread(s), profile ''%s''\n', ...
    pool.NumWorkers, pool.NumThreads, profile_name);

if usejava('jvm')
    try
        [~, sys] = memory;
        fprintf('[restart_parpool] client RAM available: %.1f GB\n', ...
            sys.PhysicalMemory.Available / 2^30);
    catch
        % memory() is Windows-only; not worth failing a run over.
    end
end
end
