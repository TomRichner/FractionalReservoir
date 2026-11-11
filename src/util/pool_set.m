function pool = pool_set(n, opts)
% POOL_SET Record this MATLAB session's pool allocation and open exactly that many workers.
%
%   pool_set(4)                                  % hold 4 workers
%   pool_set(0)                                  % hold none: close any pool
%   pool_set(6, 'label', 'matlab-fable-1')       % name the session for pool_status
%   pool_set(12, 'low_priority', true)           % workers below normal OS priority
%
% THE FIRST THING A SESSION DOES after setup_paths, once, with the number
% whoever divides the machine gave it. From then on every pool in the
% session comes from this allocation: pool_ensure() opens it where a parfor
% needs it, pool_restart() restarts it at the same size, pool_release()
% gives it back, pool_status() shows every session's share.
%
% What it does:
%   * records N (and the label and priority option, when given) in
%     environment variables of this MATLAB process (pool_config), so the
%     allocation survives `clear all` and `clear classes`. A launcher may
%     set the same variables before MATLAB starts.
%   * a pool of another size that is up is CLOSED AND REPLACED by one of N;
%     a pool of N that is up is kept. N = 0 closes whatever is up.
%   * opens on pool_config().profile, never the default profile. A failure
%     to open (no free licence seat, N above what the profile allows) is an
%     ERROR here, with the allocation left recorded so pool_ensure retries.
%   * writes this session's entry in the machine registry (pool_registry).
%
% pool_set(0) IS an allocation: it says "this session holds no pool", which
% silences pool_ensure's warning and satisfies the entry points' check. A
% session with nothing set is a different state: unknown, so it warns.
%
% The budget (pool_config().budget workers over all sessions) is not
% enforced here or anywhere: see pool_status.
%
% 'low_priority' runs the session's pool workers at BelowNormal priority so
% the desktop and remote sessions stay responsive under a full pool. It is
% remembered and re-applied each time the pool is reopened; false restores
% Normal. Windows only; a no-op elsewhere. It changes scheduling, never a
% result.
%
% See also: pool_ensure, pool_restart, pool_release, pool_status, pool_config

arguments
    n                 (1,1) double {mustBeInteger, mustBeNonnegative}
    opts.label        (1,:) char
    opts.low_priority (1,1) logical
    opts.verbose                    = 'minimal'   % 'verbose' | 'minimal' | 'near-none'
end
if ~isempty(which('getCurrentTask')) && ~isempty(getCurrentTask())
    error('pool_set:OnWorker', 'pool_set was called on a pool worker; an allocation belongs to a client session.');
end

c = pool_config();
setenv(c.env_workers, sprintf('%d', n));
if isfield(opts, 'label'); setenv(c.env_label, opts.label); end
if isfield(opts, 'low_priority'); setenv(c.env_low_priority, sprintf('%d', opts.low_priority)); end
c = pool_config();

old = gcp('nocreate');
kept = ~isempty(old) && old.NumWorkers == n;
if ~isempty(old) && ~kept
    vprintf(opts.verbose, 'verbose', '[pool] closing the %d-worker pool that is up\n', old.NumWorkers);
    delete(old);
end

if n == 0
    pool = gcp('nocreate');
else
    pool = pool_ensure('strict', true, 'verbose', opts.verbose);
    if kept && isfield(opts, 'low_priority')
        pool_registry('priority', opts.low_priority);
    end
end
pool_registry('write');

if n == 0
    vprintf(opts.verbose, 'minimal', '[pool] this session holds no pool (allocation 0)\n');
else
    vprintf(opts.verbose, 'minimal', '[pool] this session holds %d workers on profile ''%s''%s\n', ...
        pool.NumWorkers, c.profile, label_txt(c));
end
end

function s = label_txt(c)
s = '';
if ~isempty(c.label); s = sprintf(', label ''%s''', c.label); end
if c.low_priority; s = [s ', lowered priority']; end
end
