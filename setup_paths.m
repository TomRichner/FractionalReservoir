function projectRoot = setup_paths()
%SETUP_PATHS Add this repository's src/ and scripts/ trees to the MATLAB path.
%
% Call once per MATLAB session: from the repository root in a fresh session
% (this file lives there, so it resolves without any addpath), or from anywhere
% once the root is on the path. Idempotent.
%
% Optionally returns the project root, which is also recoverable elsewhere as
%   project_root = fileparts(which('setup_paths'));
%
% Also puts the shared parallel-pool module on the path (pool_set,
% pool_ensure, pool_restart, pool_release, pool_status). Since 2026-11-19
% that module lives in the MatlabAgenticOrg repository, one clone per
% machine, used by every MATLAB repository on it; use_matlab_agentic_org
% (beside this file) finds the clone and errors with the fix when this
% machine has not been set up.
%
% Deliberately does NOT genpath the project root: data/, figs/ and docs/ must
% stay off the path (data/param_space/run_all_*/ holds copies of the launcher
% scripts, which would otherwise shadow the originals in scripts/).

projectRoot = fileparts(mfilename('fullpath'));
srcPath     = fullfile(projectRoot, 'src');
scriptsPath = fullfile(projectRoot, 'scripts');

if ~isfolder(srcPath)
    error('setup_paths:MissingSrc', ...
        'Could not find src directory at %s', srcPath);
end

addpath(projectRoot);               % keeps setup_paths itself resolvable
addpath(genpath(srcPath));
if isfolder(scriptsPath)
    addpath(genpath(scriptsPath));
end
use_matlab_agentic_org();
end
