function T = pool_status(verbose)
% POOL_STATUS Who holds how many pool workers across this machine.
%
%   pool_status()
%   T = pool_status()
%   T = pool_status('near-none')     % the table only, nothing printed
%
% One row per MATLAB session that has a registry entry (it called pool_set)
% or that has pool workers running:
%
%   pid         the session's process id (feature('getpid') in that MATLAB)
%   label       what the session called itself in pool_set; '(unregistered)'
%               for a session with workers and no entry (a pool opened
%               without this module: Tom's one-worker licence-holding
%               session is the standing example)
%   allocated   the workers the session recorded with pool_set; NaN for an
%               unregistered session
%   running     the pool workers ACTUALLY RUNNING under that session,
%               counted from the operating system (pool_registry: MATLAB
%               processes started with -parallelserver whose parent is the
%               session), NOT taken from the registry. A session allocated 4
%               whose pool is down shows 4 and 0. NaN where the operating
%               system cannot be asked (not Windows, or the query failed);
%               the table is then the registry alone.
%   started     when that MATLAB started ('' when unregistered)
%   this        true for the session making the call
%
% PRUNING. An entry whose process no longer exists is deleted from the
% registry here and left out: on Windows when no MATLAB client has that
% process id, elsewhere when the process id is not alive.
%
% THE BUDGET IS REPORTED, NOT ENFORCED. The printed total compares the sum
% of the allocations with pool_config().budget, and a warning
% (pool_status:OverBudget) is raised when it is exceeded. Nothing refuses to
% open a pool on that account: exceeding it slows the machine, it does not
% break a run, and the aim is only to leave a couple of cores free.
% Unregistered sessions are outside the budget and are not counted in it.
%
% See also: pool_set, pool_ensure, pool_restart, pool_release, pool_registry

arguments
    verbose = 'minimal'
end
c = pool_config();
my_pid = c.pid;
E = pool_registry('read');
P = pool_registry('processes');

%% Prune entries whose session is gone
keep = true(1, numel(E));
for k = 1:numel(E)
    if P.ok
        alive = any(P.pid == E(k).pid & ~P.is_worker);
    elseif ispc
        alive = true;                      % cannot tell: keep
    else
        alive = system(sprintf('kill -0 %d 2>/dev/null', E(k).pid)) == 0;
    end
    if ~alive && E(k).pid ~= my_pid
        pool_registry('remove', E(k).pid);
        keep(k) = false;
    end
end
E = E(keep);

%% One row per registered session, then per unregistered session with workers
pid = [E.pid]';
label = string({E.label}');
allocated = [E.workers]';
started = string({E.started}');
running = nan(numel(E), 1);
if P.ok
    for k = 1:numel(E)
        running(k) = nnz(P.is_worker & P.parent == E(k).pid);
    end
    others = setdiff(unique(P.parent(P.is_worker)), pid);
    for k = 1:numel(others)
        pid(end + 1, 1) = others(k); %#ok<AGROW>
        label(end + 1, 1) = "(unregistered)"; %#ok<AGROW>
        allocated(end + 1, 1) = NaN; %#ok<AGROW>
        started(end + 1, 1) = ""; %#ok<AGROW>
        running(end + 1, 1) = nnz(P.is_worker & P.parent == others(k)); %#ok<AGROW>
    end
end
this = pid == my_pid;
T = table(pid, label, allocated, running, started, this);

%% Report
total_alloc = sum(allocated, 'omitnan');
vprintf(verbose, 'minimal', '[pool_status] %d MATLAB session(s) on %s; allocated %d of a budget of %d pool workers\n', ...
    height(T), host_txt(), total_alloc, c.budget);
if ~P.ok
    vprintf(verbose, 'minimal', '  (running workers cannot be counted from the operating system here: registry only)\n');
end
vprintf(verbose, 'minimal', '  %8s  %-22s %9s %8s  %-19s\n', 'pid', 'label', 'allocated', 'running', 'started');
for k = 1:height(T)
    mark = '';
    if this(k); mark = '  <- this session'; end
    vprintf(verbose, 'minimal', '  %8d  %-22s %9s %8s  %-19s%s\n', pid(k), label(k), ...
        num_txt(allocated(k)), num_txt(running(k)), started(k), mark);
end
if P.ok
    vprintf(verbose, 'minimal', '  running in total: %d (unregistered sessions included)\n', sum(running));
end
if total_alloc > c.budget
    warning('pool_status:OverBudget', ...
        ['The sessions on this machine have allocated %d pool workers between them, above the ' ...
         'budget of %d. Nothing is stopped on that account, but the machine will be slow to ' ...
         'respond; release or shrink a pool (pool_release, pool_set).'], total_alloc, c.budget);
end
end

function s = num_txt(x)
if isnan(x); s = '-'; else; s = sprintf('%d', x); end
end

function s = host_txt()
s = getenv('COMPUTERNAME');
if isempty(s); s = 'this machine'; end
end
