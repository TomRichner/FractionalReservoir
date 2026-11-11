function c = pool_config()
% POOL_CONFIG The pool module's constants and this session's allocation, read fresh.
%
%   c = pool_config()
%
% THE ONE PLACE the profile name and the worker budget are written. Every
% other pool_* function, and the few callers that need to know whether an
% allocation is set (run_all_paper_analyses, make_all_paper_figures,
% wait_for_parpool), read them from here.
%
%   c.profile       the cluster profile every pool is opened on. NOT the
%                   default profile: on Tom's workstation that is 'Just1',
%                   one worker, kept on purpose so that legacy code cannot
%                   start a large pool (FD-017).
%   c.budget        pool workers summed over every MATLAB on the machine
%                   that the working arrangement aims to stay within. It is
%                   REPORTED by pool_status, never enforced.
%   c.is_set        true when this session has an allocation, 0 included
%   c.allocation    the allocated worker count; [] when none is set
%   c.pid           this MATLAB's process id
%   c.label         the session's label ('' when none), an agent's name say
%   c.low_priority  true when the session's workers are to run at lowered
%                   operating-system priority
%   c.registry_dir  where each session writes its one small registry file
%   c.env_*         the names of the environment variables that hold the above
%   c.how_to_set    the sentence error and warning messages end with
%
% WHY ENVIRONMENT VARIABLES. The allocation has to survive `clear all` and
% `clear classes` (needed after every classdef edit), which a persistent or
% a global does not, and a launcher has to be able to set it before MATLAB
% starts. An environment variable of the MATLAB process does both. It is
% CLIENT-SIDE ONLY: pool workers copy the client's environment once, when
% the pool starts, and never see a later setenv, which is why pool_set
% replaces a pool whose size no longer matches and why no pool_* function
% does anything on a worker. An empty value counts as unset.
%
% See also: pool_set, pool_ensure, pool_restart, pool_release, pool_status

c = struct();
c.profile = 'Processes';
c.budget  = 12;

c.env_workers      = 'SRNN_POOL_WORKERS';
c.env_label        = 'SRNN_POOL_LABEL';
c.env_low_priority = 'SRNN_POOL_LOW_PRIORITY';

if ispc && ~isempty(getenv('LOCALAPPDATA'))
    c.registry_dir = fullfile(getenv('LOCALAPPDATA'), 'matlab_pool_registry');
elseif ~ispc && ~isempty(getenv('HOME'))
    c.registry_dir = fullfile(getenv('HOME'), '.matlab_pool_registry');
else
    c.registry_dir = fullfile(tempdir, 'matlab_pool_registry');
end

raw = strtrim(getenv(c.env_workers));
n = str2double(raw);
c.is_set = ~isempty(raw) && isfinite(n) && n >= 0 && n == round(n);
if c.is_set
    c.allocation = n;
else
    c.allocation = [];
    if ~isempty(raw)
        warning('pool_config:BadAllocation', ...
            ['The environment variable %s holds ''%s'', which is not a whole number of ' ...
             'workers; this session is treated as having no pool allocation.'], ...
            c.env_workers, raw);
    end
end
% The process id other sessions and the operating system know this MATLAB
% by. feature('getpid') and not matlabProcessID, which the Code Analyzer
% suggests: the latter is recent, and this repo is run on older releases too.
c.pid          = feature('getpid'); %#ok<FEATGPID>
c.label        = getenv(c.env_label);
c.low_priority = strcmp(strtrim(getenv(c.env_low_priority)), '1');

c.how_to_set = ['Run pool_set(N) first to hold N pool workers in this MATLAB session, ' ...
    'or pool_set(0) to state that it holds none.'];
end
