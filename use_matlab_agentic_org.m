function org_root = use_matlab_agentic_org()
% USE_MATLAB_AGENTIC_ORG Put the shared MatlabAgenticOrg pool module on the MATLAB path.
%
%   use_matlab_agentic_org()
%   org_root = use_matlab_agentic_org()
%
% Call it from this repository's path setup (setup_paths.m, startup.m, or
% whatever a session runs first). It adds <MatlabAgenticOrg>/matlab to the
% path, which holds pool_set, pool_ensure, pool_restart, pool_release and
% pool_status: the one copy of the parallel-pool module that every MATLAB
% repository on this machine uses, so that MATLAB sessions started from
% different repositories share one registry and one worker budget.
%
% WHERE THE CLONE IS. There is one clone of MatlabAgenticOrg per machine,
% and it is found, in this order, from
%   1. the environment variable MATLAB_AGENTIC_ORG, when set;
%   2. the one-line file .matlab_agentic_org in the user's home folder
%      (%USERPROFILE% on Windows, $HOME elsewhere), which MatlabAgenticOrg's
%      setup_repo.py writes.
% Neither lives in this repository, so a git worktree of it (a plain
% checkout of tracked files) finds the same clone as the main checkout.
% When the clone cannot be found this ERRORS, with the fix, rather than
% leaving a session that would fail later at its first pool_set.
%
% THIS FILE IS A COPY. Its source is bootstrap/use_matlab_agentic_org.m in
% MatlabAgenticOrg, and setup_repo.py puts it here. Track it in this
% repository (worktrees need it). Do not edit it here: change the source and
% rerun setup_repo.py.

fix = ['Clone git@github.com:TomRichner/MatlabAgenticOrg.git once on this machine and run, ' ...
    'from a shell,  uv run <the clone>/setup_repo.py <this repository>  ' ...
    '(it records where the clone is). See the README of MatlabAgenticOrg.'];

org_root = strtrim(getenv('MATLAB_AGENTIC_ORG'));
source = 'the environment variable MATLAB_AGENTIC_ORG';
if isempty(org_root)
    if ispc; home = getenv('USERPROFILE'); else; home = getenv('HOME'); end
    pointer = fullfile(home, '.matlab_agentic_org');
    source = pointer;
    if ~isfile(pointer)
        error('use_matlab_agentic_org:NotSetUp', ...
            ['This repository takes its parallel-pool module (pool_set, pool_status, ...) from a ' ...
             'shared clone of MatlabAgenticOrg, and this machine has no record of where that ' ...
             'clone is: %s does not exist and MATLAB_AGENTIC_ORG is not set. %s'], pointer, fix);
    end
    org_root = strtrim(fileread(pointer));
end

module_dir = fullfile(org_root, 'matlab');
if ~isfile(fullfile(module_dir, 'pool_set.m'))
    error('use_matlab_agentic_org:CloneMissing', ...
        ['%s says the MatlabAgenticOrg clone is at "%s", but there is no matlab/pool_set.m ' ...
         'there. %s'], source, org_root, fix);
end
addpath(module_dir);

% Two copies of the module on the path shadow each other silently, and
% which one wins depends on path order.
copies = which('pool_set', '-all');
if numel(copies) > 1
    warning('use_matlab_agentic_org:Shadowed', ...
        ['More than one pool_set is on the MATLAB path, so one copy of the pool module shadows ' ...
         'the other:%s\nRemove the copy that is not in %s.'], ...
        sprintf('\n    %s', copies{:}), module_dir);
end
end
