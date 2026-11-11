function pool_release(verbose)
% POOL_RELEASE Close this session's pool and give its allocation back.
%
%   pool_release()
%   pool_release(verbose)     % 'verbose' | 'minimal' | 'near-none'
%
% What a session does before it finishes, or when whoever divides the
% machine wants its workers for something else:
%   * closes the pool that is up, whatever its size;
%   * removes this session's entry from the machine registry, so
%     pool_status in other sessions stops listing it;
%   * leaves the allocation at 0, an explicit "holds no pool": library code
%     reached afterwards runs serially WITHOUT the no-allocation warning,
%     and nothing reopens a pool until pool_set(N) is called again.
%
% See also: pool_set, pool_ensure, pool_restart, pool_status

arguments
    verbose = 'minimal'
end
c = pool_config();
n = 0;
if ~isempty(which('parpool'))
    old = gcp('nocreate');
    if ~isempty(old)
        n = old.NumWorkers;
        delete(old);
    end
end
setenv(c.env_workers, '0');
pool_registry('remove');
vprintf(verbose, 'minimal', '[pool] released (%d workers closed); this session''s allocation is now 0\n', n);
end
