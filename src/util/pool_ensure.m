function pool = pool_ensure(opts)
% POOL_ENSURE This session's parallel pool, opened at its allocated size if it is down.
%
%   pool = pool_ensure()
%   pool = pool_ensure('strict', true)      % an opening failure is an error
%
% The call library code makes before a parfor. What it returns:
%
%   a pool is up            that pool, as it is, whatever its size
%   allocation N >= 1       a pool of N workers on pool_config().profile,
%                           opened now
%   allocation 0            empty, quietly: the session was told to hold no
%                           pool, and the parfor runs serially in the client
%   NO allocation set       empty, with a LOUD WARNING (pool_ensure:NoAllocation)
%                           that the code will run serially, and how to fix it
%   on a pool worker        empty, quietly (a worker cannot hold a pool)
%
% So an empty result always means "run serially", and callers need no other
% branch: with automatic pool creation off (the machine's setting), a parfor
% reached with no pool runs in the client.
%
% SERIAL IS THE SAME COMPUTATION BUT NOT THE SAME BYTES. Pool workers run
% one computational thread each and the client runs many, so linear algebra
% done in the client rounds differently: on 2026-11-11 the stages regression
% run with allocation 0 differed from its reference in 466 leaves at 1e-13
% to 1e-15 relative, while pools of 4 and of 14 workers gave identical
% bytes. The SIZE of a pool moves nothing; having one or not does, in the
% last bits. The byte-exact stages and full regression tiers therefore need
% a pool up (the ultrarapid tier passes either way).
%
% IT NEVER OPENS THE DEFAULT PROFILE AND NEVER INVENTS A SIZE. The size is
% the session's allocation (pool_set), nothing else: a stage's own idea of
% how many workers it would like no longer opens anything. That is the point
% of the module: several MATLAB sessions share one machine, and only whoever
% divides the machine between them knows what each may take.
%
% The NoAllocation warning is raised at most ONCE PER CALL SITE (the file
% and line that called pool_ensure) per MATLAB session, so a loop over
% stages does not repeat it; `clear pool_ensure` resets that memory.
%
% If the pool cannot be opened (no free licence seat, say) the default is a
% warning (pool_ensure:OpenFailed) and an empty result, so an unattended run
% carries on serially rather than dying; pool_set and wait_for_parpool pass
% 'strict' to get the error instead.
%
% The pool is opened with the profile's own settings, idle timeout included:
% a pool left idle shuts down, and the next pool_ensure reopens it.
%
% See also: pool_set, pool_restart, pool_release, pool_status, pool_config

arguments
    opts.strict  (1,1) logical = false
    opts.verbose               = 'minimal'   % 'verbose' | 'minimal' | 'near-none'
end
persistent warned_sites
if isempty(warned_sites); warned_sites = {}; end

pool = [];
if isempty(which('parpool')) || ~isempty(getCurrentTask())
    return          % no toolbox, or running on a pool worker
end
pool = gcp('nocreate');
if ~isempty(pool); return; end

c = pool_config();
if ~c.is_set
    st = dbstack(1);
    if isempty(st)
        site = 'the command line';
    else
        site = sprintf('%s (line %d)', st(1).name, st(1).line);
    end
    if ~ismember(site, warned_sites)
        warned_sites{end + 1} = site;
        warning('pool_ensure:NoAllocation', ...
            ['NO POOL ALLOCATION IS SET for this MATLAB session, so the parallel loops reached ' ...
             'from %s run SERIALLY in the client: far more slowly, and equal to a pool run ' ...
             'only to round-off (the client is multithreaded, pool workers are not). ' ...
             'No pool was opened. %s'], site, c.how_to_set);
    end
    return
end
if c.allocation == 0; return; end

try
    pool = parpool(c.profile, c.allocation);
catch ME
    if opts.strict; rethrow(ME); end
    warning('pool_ensure:OpenFailed', ...
        ['Could not open the %d-worker pool allocated to this session on profile ''%s'' (%s). ' ...
         'Running SERIALLY in the client.'], c.allocation, c.profile, strtok(ME.message, newline));
    pool = [];
    return
end
pool_registry('write');
n_lowered = 0;
if c.low_priority; n_lowered = pool_registry('priority', true); end
vprintf(opts.verbose, 'verbose', '[pool] opened %d workers on profile ''%s''%s\n', ...
    pool.NumWorkers, c.profile, priority_txt(n_lowered));
end

function s = priority_txt(n_lowered)
if n_lowered > 0; s = sprintf(', %d at lowered priority', n_lowered); else; s = ''; end
end
