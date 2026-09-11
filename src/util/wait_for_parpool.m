function pool = wait_for_parpool(n_workers, opts)
% WAIT_FOR_PARPOOL Poll the network licence until a Parallel Computing Toolbox
% seat is free, then start a pool, which holds the seat for this session.
%
%   pool = WAIT_FOR_PARPOOL()               % min(12, cores) workers
%   pool = WAIT_FOR_PARPOOL(8)
%   pool = WAIT_FOR_PARPOOL(12, 'poll_s', 20, 'max_wait_s', 4*3600)
%
% The toolbox is a 15-seat floating licence here and on 2026-09-11 every seat
% was taken. A pool needs ONE seat (for the client, not one per worker), and
% once the pool is up the seat stays checked out until the pool is deleted or
% MATLAB exits -- so run this first, then launch whatever needs the pool.
%
% Each poll asks the licence server through scripts/tools/pct_licenses.ps1
% (lmutil lmstat; the script's exit code is the number of free seats). When
% it reports a free seat, or cannot answer, parpool is attempted; if someone
% else took the seat first the attempt fails and polling resumes. An existing
% pool is returned as is, whatever its size.
%
% HOST GATE. The script shells out to lmutil at a path on one workstation and
% this is meant to sit and wait for hours on that machine; on any other host
% it errors rather than guessing (TR, 2026-09-11). Override with 'host'.
%
% See also: restart_parpool, parpool, scripts/tools/pct_licenses.ps1

arguments
    n_workers       (1,1) double = 0            % 0 -> min(12, cores)
    opts.poll_s     (1,1) double = 20           % seconds between licence checks
    opts.max_wait_s (1,1) double = Inf          % give up (error) after this long
    opts.host       (1,:) char   = 'R5456622'   % the workstation this is written for
    opts.verbose    (1,1) logical = true
end

host = getenv('COMPUTERNAME');
if ~strcmpi(host, opts.host)
    error('wait_for_parpool:WrongHost', ...
        'wait_for_parpool is written for %s; this is %s. Pass ''host'' to override.', ...
        opts.host, host);
end
if n_workers <= 0
    n_workers = min(12, feature('numcores'));
end

pool = gcp('nocreate');
if ~isempty(pool)
    if opts.verbose
        fprintf('[wait_for_parpool] a pool is already up (%d workers); using it.\n', pool.NumWorkers);
    end
    return
end

script = fullfile(fileparts(which('setup_paths')), 'scripts', 'tools', 'pct_licenses.ps1');
t0 = tic;
n_polls = 0;
while true
    free = pct_free_seats(script);
    if isnan(free) || free > 0
        try
            pool = parpool(parallel.defaultProfile, n_workers);
            if opts.verbose
                fprintf('[wait_for_parpool] pool up: %d workers after %.1f min (%d polls).\n', ...
                    pool.NumWorkers, toc(t0) / 60, n_polls);
            end
            return
        catch ME
            last = strtok(ME.message, newline);
        end
    else
        last = 'licence server reports 0 free seats';
    end
    if toc(t0) > opts.max_wait_s
        error('wait_for_parpool:Timeout', ...
            'No pool after %.0f min (%s).', toc(t0) / 60, last);
    end
    n_polls = n_polls + 1;
    if opts.verbose && (n_polls == 1 || mod(n_polls, 15) == 0)
        fprintf('[wait_for_parpool] %s; polling every %g s, %.0f min so far (%s)\n', ...
            last, opts.poll_s, toc(t0) / 60, datestr(now, 'HH:MM')); %#ok<TNOW1,DATST>
    end
    pause(opts.poll_s);
end
end

function free = pct_free_seats(script)
% Free seats per the licence server, or NaN if the question cannot be answered.
free = NaN;
if ~ispc || ~isfile(script); return; end
[status, ~] = system(sprintf( ...
    'powershell -NoProfile -ExecutionPolicy Bypass -File "%s" -Quiet', script));
if status >= 0 && status < 200      % exit code = free seats; 254/255 are the script's own errors
    free = status;
end
end
