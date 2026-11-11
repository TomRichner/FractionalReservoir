function out = pool_registry(action, varargin)
% POOL_REGISTRY The machine-level record of who holds a pool, and what the operating system says.
%
%   pool_registry('write')            write or refresh THIS session's entry
%   pool_registry('remove')           delete this session's entry
%   pool_registry('remove', pid)      delete another session's entry (pruning)
%   E = pool_registry('read')         every entry, as a struct array
%   P = pool_registry('processes')    the MATLAB processes the OS reports
%   n = pool_registry('priority', lowered)   set this session's workers' priority
%
% A helper of the pool module; call pool_set / pool_status, not this.
%
% THE REGISTRY is one small JSON file per MATLAB session,
% <pool_config().registry_dir>/session_<pid>.json, holding the process id,
% the session's label, the workers it has ALLOCATED, when that MATLAB
% started, and the host. It exists so that several MATLAB sessions on one
% machine (a lead agent's and its workers') can see each other. It is for
% VISIBILITY ONLY: nothing reads it to decide whether a pool may open.
%
% It is not trusted for what is actually running. A session that died
% leaves its file behind, and a pool can be opened or closed without this
% module. So 'processes' asks the operating system: on Windows a pool
% worker is a MATLAB.exe with -parallelserver on its command line whose
% parent process id is the client's feature('getpid') (checked 2026-11-11,
% FD-017). P.ok is false where that question cannot be answered (any other
% platform, or a failed query), and callers then fall back to the registry
% alone:
%   P.ok          the listing below is real
%   P.pid         every MATLAB process
%   P.parent      its parent's process id
%   P.is_worker   true for a pool worker
%
% 'priority' lowers (true) or restores (false) the operating-system priority
% of this session's pool workers, BelowNormal against Normal, so that a
% machine running a full pool stays responsive at the desktop. It returns
% the number of workers changed. Windows only: elsewhere it does nothing and
% returns 0. It does not pin workers to cores.
%
% See also: pool_config, pool_set, pool_status

c = pool_config();
my_pid = c.pid;
switch action
    case 'write'
        if ~isfolder(c.registry_dir); mkdir(c.registry_dir); end
        entry = struct('pid', my_pid, 'label', c.label, 'workers', 0, ...
            'started', session_start(), 'host', host_name());
        if c.is_set; entry.workers = c.allocation; end
        % Written whole to a temporary name and then moved, so a reader in
        % another session never sees half a file.
        target = entry_file(c, my_pid);
        tmp = [target '.tmp'];
        fid = fopen(tmp, 'w');
        if fid < 0
            warning('pool_registry:WriteFailed', ...
                'Could not write the pool registry entry %s; pool_status in other sessions will not list this one.', target);
            out = '';
            return
        end
        fprintf(fid, '%s\n', jsonencode(entry));
        fclose(fid);
        movefile(tmp, target, 'f');
        out = target;

    case 'remove'
        pid = my_pid;
        if ~isempty(varargin); pid = varargin{1}; end
        target = entry_file(c, pid);
        if isfile(target); delete(target); end
        out = target;

    case 'read'
        out = struct('pid', {}, 'label', {}, 'workers', {}, 'started', {}, 'host', {});
        listing = dir(fullfile(c.registry_dir, 'session_*.json'));
        for k = 1:numel(listing)
            try
                e = jsondecode(fileread(fullfile(listing(k).folder, listing(k).name)));
                out(end + 1) = struct('pid', double(e.pid), 'label', char(e.label), ...
                    'workers', double(e.workers), 'started', char(e.started), ...
                    'host', char(e.host)); %#ok<AGROW>
            catch
                % Unreadable or half-written by a session that is exiting: not
                % worth failing a status call over.
            end
        end

    case 'processes'
        out = os_processes();

    case 'priority'
        out = set_priority(my_pid, varargin{1});

    otherwise
        error('pool_registry:BadAction', 'Unknown action ''%s''.', action);
end
end

%% ------------------------------------------------------------------------
function f = entry_file(c, pid)
f = fullfile(c.registry_dir, sprintf('session_%d.json', pid));
end

function s = host_name()
s = getenv('COMPUTERNAME');
if isempty(s); s = getenv('HOSTNAME'); end
end

function s = session_start()
% When this MATLAB started, from the JVM's own start time (the JVM starts
% with MATLAB). Without a JVM, the time of this call.
t = datetime('now');
if usejava('jvm')
    try
        ms = java.lang.management.ManagementFactory.getRuntimeMXBean().getStartTime();
        t = datetime(double(ms) / 1000, 'ConvertFrom', 'posixtime', 'TimeZone', 'local');
    catch
    end
end
s = char(t, 'yyyy-MM-dd HH:mm:ss');
end

function P = os_processes()
P = struct('ok', false, 'pid', [], 'parent', [], 'is_worker', false(0, 1));
if ~ispc; return; end
cmd = ['powershell -NoProfile -Command "Get-CimInstance Win32_Process ' ...
    '-Filter \"Name=''MATLAB.exe''\" | ForEach-Object { ''{0}|{1}|{2}'' -f ' ...
    '$_.ProcessId,$_.ParentProcessId,($_.CommandLine -match ''parallelserver'') }"'];
[status, txt] = system(cmd);
if status ~= 0; return; end
tok = regexp(txt, '(\d+)\|(\d+)\|(True|False)', 'tokens');
if isempty(tok); return; end      % this MATLAB at least must be listed
tok = vertcat(tok{:});
P.pid       = str2double(tok(:, 1));
P.parent    = str2double(tok(:, 2));
P.is_worker = strcmp(tok(:, 3), 'True');
P.ok        = true;
end

function n = set_priority(my_pid, lowered)
n = 0;
if ~ispc; return; end
P = os_processes();
if ~P.ok; return; end
ids = P.pid(P.is_worker & P.parent == my_pid);
if isempty(ids); return; end
if lowered; cls = 'BelowNormal'; else; cls = 'Normal'; end
cmd = sprintf(['powershell -NoProfile -Command "Get-Process -Id %s | ' ...
    'ForEach-Object { $_.PriorityClass = ''%s'' }"'], ...
    strjoin(compose('%d', ids(:)'), ','), cls);
[status, txt] = system(cmd);
if status ~= 0
    warning('pool_registry:PriorityFailed', ...
        'Could not set the pool workers'' priority to %s: %s', cls, strtrim(txt));
    return
end
n = numel(ids);
end
